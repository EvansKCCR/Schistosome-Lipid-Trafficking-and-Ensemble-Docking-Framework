#!/usr/bin/env python3
"""
Multi-stage ensemble docking and substrate-specificity screening pipeline.

Modules
-------
1. DiffDock-L / SMINA output clustering and representative-pose selection.
2. Ligand and receptor preparation script generation for Meeko/AutoDock Vina.
3. AutoDock Vina refinement log parsing and ranking.
4. Post-refinement protein-ligand contact analysis, including CD36RP-style automated Vina mode splitting and complex construction.
5. CEH-like enzyme catalytic-geometry engine plus optional flexdock_v3 flexible-output parsing, clustering, representative selection, and complex construction.

The pipeline is intentionally modular: each step can be run independently, or the
`workflow` command can orchestrate all available steps in one project directory.

Typical workflow
----------------
python multistage_ensemble_docking_pipeline.py workflow \
  --diffdock-csv DiffDock_L_results.csv \
  --sdf-dir diffdock_sdf \
  --receptor-pdb receptor.pdb \
  --ligand-pdbqt-dir ligand_pdbqt \
  --vina-log-dir vina_logs \
  --complex-dir vina_complexes \
  --engine ceh \
  --outdir docking_screen
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# General utilities
# ---------------------------------------------------------------------------

NUMERIC_NA = {"", "NA", "NaN", "nan", "None", None}


def mkdir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def safe_stem(text: str) -> str:
    stem = Path(str(text)).stem
    return re.sub(r"[^A-Za-z0-9._-]+", "_", stem).strip("_") or "model"


def to_float(value, default=np.nan) -> float:
    if value in NUMERIC_NA:
        return default
    try:
        x = float(value)
        return x if math.isfinite(x) else default
    except Exception:
        return default


def robust_scale(series: pd.Series, direction: str = "higher_better") -> pd.Series:
    x = pd.to_numeric(series, errors="coerce")
    if x.notna().sum() == 0:
        return pd.Series(0.5, index=series.index)
    med = x.median()
    x = x.fillna(med)
    p5 = np.nanpercentile(x, 5)
    p95 = np.nanpercentile(x, 95)
    if not np.isfinite(p5) or not np.isfinite(p95) or p95 == p5:
        return pd.Series(0.5, index=series.index)
    s = ((x.clip(p5, p95) - p5) / (p95 - p5)).clip(0, 1)
    if direction == "lower_better":
        s = 1 - s
    return s


def weighted_sum(df: pd.DataFrame, weights: Dict[str, float]) -> pd.Series:
    cols = [c for c in weights if c in df.columns and pd.to_numeric(df[c], errors="coerce").notna().any()]
    if not cols:
        return pd.Series(0.5, index=df.index)
    w = np.array([weights[c] for c in cols], dtype=float)
    w = w / w.sum()
    vals = df[cols].apply(pd.to_numeric, errors="coerce").fillna(0.5)
    return (vals * w).sum(axis=1)


def write_excel(path: Path, sheets: Dict[str, pd.DataFrame]):
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        for name, df in sheets.items():
            if df is None or df.empty:
                continue
            safe = re.sub(r"[\\/*?:\[\]]", "_", str(name))[:31]
            df.to_excel(writer, sheet_name=safe, index=False)


# ---------------------------------------------------------------------------
# DiffDock-L / SMINA clustering and top-pose selection
# ---------------------------------------------------------------------------

DIFFDOCK_COLUMN_ALIASES = {
    "receptor": "Receptor",
    "rec": "Receptor",
    "protein": "Receptor",
    "ligand": "Ligand",
    "lig": "Ligand",
    "prediction": "Prediction",
    "pose": "Prediction",
    "sdf": "Prediction",
    "file": "Prediction",
    "filename": "Prediction",
    "diffdock confidence": "DiffDock_Confidence",
    "diffdock_confidence": "DiffDock_Confidence",
    "confidence": "DiffDock_Confidence",
    "smina minimized affinity": "SMINA_Minimized_Affinity",
    "smina_minimized_affinity": "SMINA_Minimized_Affinity",
    "minimized affinity": "SMINA_Minimized_Affinity",
    "smina minimized rmsd": "SMINA_Minimized_RMSD",
    "smina_minimized_rmsd": "SMINA_Minimized_RMSD",
    "minimized rmsd": "SMINA_Minimized_RMSD",
    "smina affinity": "SMINA_Affinity",
    "smina_affinity": "SMINA_Affinity",
}


def normalize_diffdock_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = [str(c).strip() for c in df.columns]
    rename = {}
    for c in df.columns:
        key = c.lower().strip()
        if key in DIFFDOCK_COLUMN_ALIASES:
            rename[c] = DIFFDOCK_COLUMN_ALIASES[key]
    df = df.rename(columns=rename)
    required = ["Receptor", "Ligand", "Prediction"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"DiffDock/SMINA table is missing required column(s): {missing}")
    for c in ["DiffDock_Confidence", "SMINA_Minimized_Affinity", "SMINA_Minimized_RMSD", "SMINA_Affinity"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    if "SMINA_Minimized_Affinity" not in df.columns:
        df["SMINA_Minimized_Affinity"] = np.nan
    if "SMINA_Minimized_RMSD" not in df.columns:
        df["SMINA_Minimized_RMSD"] = np.nan
    if "DiffDock_Confidence" not in df.columns:
        df["DiffDock_Confidence"] = np.nan
    return df


def rank_pose_table(df: pd.DataFrame, confidence_mode: str = "auto") -> pd.DataFrame:
    out = df.copy()
    out["affinity_scaled"] = robust_scale(out["SMINA_Minimized_Affinity"], "lower_better")
    out["rmsd_scaled"] = robust_scale(out["SMINA_Minimized_RMSD"], "lower_better")
    if confidence_mode == "higher_better":
        conf_scaled = robust_scale(out["DiffDock_Confidence"], "higher_better")
    elif confidence_mode == "lower_better":
        conf_scaled = robust_scale(out["DiffDock_Confidence"], "lower_better")
    else:
        # DiffDock-style confidence exports vary by workflow; infer from correlation
        # with SMINA affinity where possible. If higher confidence tracks more negative
        # affinity, treat confidence as higher_better.
        mask = out["DiffDock_Confidence"].notna() & out["SMINA_Minimized_Affinity"].notna()
        if mask.sum() >= 5:
            corr = np.corrcoef(out.loc[mask, "DiffDock_Confidence"], out.loc[mask, "SMINA_Minimized_Affinity"])[0, 1]
            mode = "lower_better" if np.isfinite(corr) and corr > 0 else "higher_better"
        else:
            mode = "higher_better"
        conf_scaled = robust_scale(out["DiffDock_Confidence"], mode)
        out["DiffDock_Confidence_direction_inferred"] = mode
    out["global_pose_score"] = weighted_sum(out, {
        "affinity_scaled": 0.55,
        "rmsd_scaled": 0.25,
        "conf_scaled": 0.20,
    }) if False else (0.55*out["affinity_scaled"] + 0.25*out["rmsd_scaled"] + 0.20*conf_scaled)
    out["confidence_scaled"] = conf_scaled
    return out


def find_pose_file(sdf_dir: Path, prediction: str) -> Optional[Path]:
    pred = str(prediction)
    candidates = [sdf_dir / pred, sdf_dir / Path(pred).name]
    if not pred.lower().endswith(".sdf"):
        candidates += [sdf_dir / f"{pred}.sdf", sdf_dir / f"{Path(pred).stem}.sdf"]
    for c in candidates:
        if c.exists():
            return c
    # recursive loose fallback
    target = Path(pred).name
    for c in sdf_dir.rglob("*.sdf"):
        if c.name == target or c.stem == Path(target).stem:
            return c
    return None


def _rdkit_heavy_rmsd_matrix(paths: List[Path], largest_fragment: bool = True) -> Optional[np.ndarray]:
    try:
        from rdkit import Chem
        from rdkit.Chem import rdMolAlign
        from rdkit.ML.Cluster import Butina
    except Exception:
        return None

    def load_mol(p: Path):
        suppl = Chem.SDMolSupplier(str(p), removeHs=False, sanitize=False)
        mol = next((m for m in suppl if m is not None), None)
        if mol is None:
            return None
        try:
            Chem.SanitizeMol(mol)
            Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
        except Exception:
            pass
        if largest_fragment:
            frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
            if frags:
                mol = max(frags, key=lambda x: x.GetNumHeavyAtoms())
        return mol

    mols = [load_mol(p) for p in paths]
    if any(m is None for m in mols):
        return None
    n = len(mols)
    mat = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i+1, n):
            try:
                val = float(rdMolAlign.GetBestRMS(mols[i], mols[j]))
            except Exception:
                val = 1e6
            mat[i, j] = mat[j, i] = val
    return mat


def _cluster_from_distance_matrix(dist: np.ndarray, cutoff: float) -> List[int]:
    n = dist.shape[0]
    if n == 0:
        return []
    if n == 1:
        return [1]
    try:
        from scipy.cluster.hierarchy import linkage, fcluster
        from scipy.spatial.distance import squareform
        labels = fcluster(linkage(squareform(dist), method="average"), cutoff, criterion="distance")
        return [int(x) for x in labels]
    except Exception:
        # greedy fallback
        labels = [0] * n
        next_id = 1
        for i in range(n):
            if labels[i]:
                continue
            labels[i] = next_id
            for j in range(i+1, n):
                if dist[i, j] <= cutoff:
                    labels[j] = next_id
            next_id += 1
        return labels


def _score_space_clusters(df: pd.DataFrame, min_aff_bin=0.25, min_rmsd_bin=0.5) -> List[int]:
    aff = pd.to_numeric(df["SMINA_Minimized_Affinity"], errors="coerce")
    rmsd = pd.to_numeric(df["SMINA_Minimized_RMSD"], errors="coerce")
    aff = aff.fillna(aff.median() if aff.notna().any() else 0)
    rmsd = rmsd.fillna(rmsd.median() if rmsd.notna().any() else 0)
    def bins(x, step):
        if len(x) == 0:
            return []
        origin = float(np.nanmin(x))
        spread = float(np.nanstd(x)) / 2.0
        step = max(step, spread) if np.isfinite(spread) else step
        return np.floor((x - origin) / step).astype(int).tolist()
    keys = list(zip(bins(aff.to_numpy(float), min_aff_bin), bins(rmsd.to_numpy(float), min_rmsd_bin)))
    ordered_keys = {k: i+1 for i, k in enumerate(sorted(set(keys)))}
    return [ordered_keys[k] for k in keys]


def cluster_diffdock(
    csv_path: Path,
    sdf_dir: Path,
    outdir: Path,
    rmsd_cutoff: float = 2.0,
    top_k_clusters: int = 3,
    top_n_poses: int = 10,
    min_aff_cutoff: Optional[float] = None,
    max_rmsd_cutoff: Optional[float] = None,
    confidence_mode: str = "auto",
    drop_zero_scores: bool = False,
    largest_fragment: bool = True,
) -> Dict[str, pd.DataFrame]:
    mkdir(outdir)
    selected_dir = mkdir(outdir / "selected_poses")
    df = pd.read_csv(csv_path)
    df = normalize_diffdock_columns(df)
    df = rank_pose_table(df, confidence_mode=confidence_mode)

    if drop_zero_scores and "SMINA_Affinity" in df.columns:
        df = df[(pd.to_numeric(df["SMINA_Affinity"], errors="coerce") != 0) | df["SMINA_Affinity"].isna()].copy()

    memberships = []
    clusters = []
    representatives = []
    topN = []
    dropped = []

    for (receptor, ligand), g in df.groupby(["Receptor", "Ligand"], dropna=False):
        g = g.copy().reset_index(drop=True)
        pose_paths = []
        keep_rows = []
        for idx, row in g.iterrows():
            pose_path = find_pose_file(sdf_dir, row["Prediction"])
            if pose_path is None:
                dropped.append({"Receptor": receptor, "Ligand": ligand, "Prediction": row["Prediction"], "Reason": "SDF_not_found"})
                continue
            pose_paths.append(pose_path)
            keep_rows.append(idx)
        if not keep_rows:
            continue
        gg = g.loc[keep_rows].copy().reset_index(drop=True)
        dist = _rdkit_heavy_rmsd_matrix(pose_paths, largest_fragment=largest_fragment)
        if dist is None:
            labels = _score_space_clusters(gg)
            clustering_method = "score_space_fallback"
        else:
            labels = _cluster_from_distance_matrix(dist, rmsd_cutoff)
            clustering_method = "ligand_heavy_atom_rmsd"
        gg["ClusterID"] = labels
        gg["PosePath"] = [str(p) for p in pose_paths]
        gg["ClusteringMethod"] = clustering_method
        memberships.append(gg)

        # cluster summaries
        for cluster_id, cg in gg.groupby("ClusterID"):
            best_idx = cg.sort_values(["global_pose_score", "SMINA_Minimized_Affinity", "SMINA_Minimized_RMSD"], ascending=[False, True, True]).index[0]
            best = gg.loc[best_idx]
            summary = {
                "Receptor": receptor,
                "Ligand": ligand,
                "ClusterID": int(cluster_id),
                "ClusterSize": int(cg.shape[0]),
                "BestPrediction": best["Prediction"],
                "BestPosePath": best["PosePath"],
                "BestGlobalPoseScore": float(best["global_pose_score"]),
                "Best_Minimized_Affinity": float(best["SMINA_Minimized_Affinity"]) if pd.notna(best["SMINA_Minimized_Affinity"]) else np.nan,
                "Median_Minimized_Affinity": float(cg["SMINA_Minimized_Affinity"].median()) if cg["SMINA_Minimized_Affinity"].notna().any() else np.nan,
                "Median_Minimized_RMSD": float(cg["SMINA_Minimized_RMSD"].median()) if cg["SMINA_Minimized_RMSD"].notna().any() else np.nan,
                "ClusteringMethod": clustering_method,
            }
            clusters.append(summary)
        cdf = pd.DataFrame([r for r in clusters if r["Receptor"] == receptor and r["Ligand"] == ligand])
        if not cdf.empty:
            cdf = cdf.sort_values(["ClusterSize", "BestGlobalPoseScore", "Best_Minimized_Affinity"], ascending=[False, False, True])
            reps = cdf.head(top_k_clusters).copy()
            reps["RepresentativeRank"] = np.arange(1, len(reps)+1)
            representatives.append(reps)

        # top-N docking-ready pose pool
        pool = gg.copy()
        if min_aff_cutoff is not None:
            pool = pool[pool["SMINA_Minimized_Affinity"] <= min_aff_cutoff]
        if max_rmsd_cutoff is not None:
            pool = pool[pool["SMINA_Minimized_RMSD"] <= max_rmsd_cutoff]
        if pool.empty:
            pool = gg.copy()
        pool = pool.sort_values(["global_pose_score", "SMINA_Minimized_Affinity", "SMINA_Minimized_RMSD"], ascending=[False, True, True]).head(top_n_poses).copy()
        pool.insert(0, "TopPoseRank", np.arange(1, len(pool)+1))
        topN.append(pool)
        # copy selected SDFs
        for _, row in pool.iterrows():
            src = Path(row["PosePath"])
            dst = selected_dir / f"{safe_stem(str(receptor))}__{safe_stem(str(ligand))}__rank{int(row['TopPoseRank']):02d}__{src.name}"
            if src.exists() and not dst.exists():
                shutil.copy2(src, dst)

    out = {
        "membership": pd.concat(memberships, ignore_index=True) if memberships else pd.DataFrame(),
        "clusters": pd.DataFrame(clusters),
        "representatives": pd.concat(representatives, ignore_index=True) if representatives else pd.DataFrame(),
        "top_poses": pd.concat(topN, ignore_index=True) if topN else pd.DataFrame(),
        "dropped": pd.DataFrame(dropped),
    }
    for name, table in out.items():
        table.to_csv(outdir / f"diffdock_{name}.tsv", sep="\t", index=False)
    write_excel(outdir / "diffdock_pose_selection.xlsx", out)
    plot_diffdock(out, outdir / "plots")
    return out


def plot_diffdock(tables: Dict[str, pd.DataFrame], plot_dir: Path):
    mkdir(plot_dir)
    top = tables.get("top_poses", pd.DataFrame())
    if not top.empty and "SMINA_Minimized_Affinity" in top.columns:
        fig, ax = plt.subplots(figsize=(8, 5), dpi=200)
        vals = pd.to_numeric(top["SMINA_Minimized_Affinity"], errors="coerce").dropna()
        if not vals.empty:
            ax.hist(vals, bins=20)
            ax.set_xlabel("SMINA minimized affinity (kcal/mol)")
            ax.set_ylabel("Pose count")
            ax.set_title("Selected DiffDock-L/SMINA pose affinity distribution")
            fig.tight_layout()
            fig.savefig(plot_dir / "selected_pose_affinity_distribution.png")
        plt.close(fig)
    cl = tables.get("clusters", pd.DataFrame())
    if not cl.empty:
        data = cl.sort_values("ClusterSize", ascending=False).head(30)
        labels = data.apply(lambda r: f"{r['Receptor']}|{r['Ligand']}|C{int(r['ClusterID'])}", axis=1)
        fig, ax = plt.subplots(figsize=(9, max(4, 0.25*len(data))), dpi=200)
        ax.barh(range(len(data)), data["ClusterSize"])
        ax.set_yticks(range(len(data)))
        ax.set_yticklabels(labels, fontsize=7)
        ax.invert_yaxis()
        ax.set_xlabel("Cluster size")
        ax.set_title("Major pose clusters")
        fig.tight_layout()
        fig.savefig(plot_dir / "pose_cluster_sizes.png")
        plt.close(fig)


# ---------------------------------------------------------------------------
# Receptor/ligand/Vina script generation
# ---------------------------------------------------------------------------


def generate_ligand_preparation_script(outdir: Path, selected_sdf_dir: str = "selected_poses") -> Path:
    mkdir(outdir)
    script = outdir / "prepare_ligands_meeko.sh"
    script.write_text(f"""#!/bin/bash
set -euo pipefail
mkdir -p ligand_sdf_H ligand_pdbqt logs
for sdf in {selected_sdf_dir}/*.sdf; do
    [ -e "$sdf" ] || continue
    base=$(basename "$sdf" .sdf)
    echo "Processing $base"
    obabel "$sdf" -O "ligand_sdf_H/${{base}}_H.sdf" -h > "logs/${{base}}_obabel.log" 2>&1
    mk_prepare_ligand.py -i "ligand_sdf_H/${{base}}_H.sdf" -o "ligand_pdbqt/${{base}}.pdbqt" > "logs/${{base}}_meeko.log" 2>&1
    test -s "ligand_pdbqt/${{base}}.pdbqt" && echo "OK $base" || echo "FAILED $base"
done
""")
    script.chmod(0o755)
    return script


def generate_receptor_preparation_script(
    outdir: Path,
    receptor_pdb: str,
    box_center: Optional[Tuple[float, float, float]] = None,
    box_size: Optional[Tuple[float, float, float]] = None,
    flex_residues: str = "",
    receptor_name: str = "receptor",
    box_enveloping: Optional[str] = None,
    padding: float = 4.0,
    mode: str = "auto",
) -> Path:
    """Generate a receptor/grid preparation script.

    `mode=general-rigid` follows the project-specific prepare_receptor.txt
    pattern: a rigid receptor is prepared with a docking box enveloping a
    selected ligand pose, e.g.

      mk_prepare_receptor.py -i pdb_fixed/RP1_receptor.pdb \
        -p receptor_grid/RP1_receptor_grid.pdbqt \
        --box_enveloping ligand_pdbqt/rank1_rp1_ce18_1.pdbqt \
        --padding 4 -v receptor_grid/RP1_box.txt

    `mode=explicit-box` preserves the earlier center/size-grid behaviour and
    supports optional flexible residues for CEH-like workflows.
    """
    mkdir(outdir)
    script = outdir / "prepare_receptor_grid.sh"
    flex = f" -f {flex_residues} -a" if flex_residues else ""
    mode_norm = (mode or "auto").lower()
    use_enveloping = bool(box_enveloping) and mode_norm in {"auto", "general", "general-rigid", "enveloping"}

    if use_enveloping or mode_norm in {"general", "general-rigid", "enveloping"}:
        if not box_enveloping:
            box_enveloping = "ligand_pdbqt/rank1_pose.pdbqt"
        script.write_text(f"""#!/bin/bash
set -euo pipefail
mkdir -p receptor_grid
mk_prepare_receptor.py -i "{receptor_pdb}" -p "receptor_grid/{receptor_name}_receptor_grid.pdbqt" --box_enveloping "{box_enveloping}" --padding {padding} -v "receptor_grid/{receptor_name}_box.txt"{flex}
""")
    else:
        center = box_center or (0.0, 0.0, 0.0)
        size = box_size or (24.0, 24.0, 24.0)
        script.write_text(f"""#!/bin/bash
set -euo pipefail
mkdir -p receptor_grid
mk_prepare_receptor.py -i "{receptor_pdb}" -o "receptor_grid/{receptor_name}" -p -v --box_size {size[0]} {size[1]} {size[2]} --box_center {center[0]} {center[1]} {center[2]}{flex}
""")
    script.chmod(0o755)
    return script


def generate_vina_refinement_script(
    outdir: Path,
    receptor_rigid: str = "receptor_grid/receptor_rigid.pdbqt",
    receptor_flex: Optional[str] = None,
    box_file: str = "receptor_grid/receptor_box.txt",
    ligand_dir: str = "ligand_pdbqt",
    out_subdir: str = "vina_results",
    log_subdir: str = "vina_logs",
    exhaustiveness: int = 64,
    num_modes: int = 50,
    energy_range: float = 3.0,
    seed: int = 42,
) -> Path:
    mkdir(outdir)
    flex_line = f"        --flex {receptor_flex} \\\n" if receptor_flex else ""
    script = outdir / "run_vina_refinement.sh"
    script.write_text(f"""#!/bin/bash
set -euo pipefail
RECEPTOR=\"{receptor_rigid}\"
BOX=\"{box_file}\"
OUTDIR=\"{out_subdir}\"
LOGDIR=\"{log_subdir}\"
mkdir -p "$OUTDIR" "$LOGDIR"
CENTER_X=$(grep center_x "$BOX" | awk '{{print $3}}')
CENTER_Y=$(grep center_y "$BOX" | awk '{{print $3}}')
CENTER_Z=$(grep center_z "$BOX" | awk '{{print $3}}')
SIZE_X=$(grep size_x "$BOX" | awk '{{print $3}}')
SIZE_Y=$(grep size_y "$BOX" | awk '{{print $3}}')
SIZE_Z=$(grep size_z "$BOX" | awk '{{print $3}}')
for lig in {ligand_dir}/*.pdbqt; do
    [ -e "$lig" ] || continue
    base=$(basename "$lig" .pdbqt)
    echo "Running Vina refinement for $base"
    mkdir -p "$OUTDIR/$base"
    vina \\
        --receptor "$RECEPTOR" \\
{flex_line}        --ligand "$lig" \\
        --center_x "$CENTER_X" --center_y "$CENTER_Y" --center_z "$CENTER_Z" \\
        --size_x "$SIZE_X" --size_y "$SIZE_Y" --size_z "$SIZE_Z" \\
        --exhaustiveness {exhaustiveness} \\
        --num_modes {num_modes} \\
        --energy_range {energy_range} \\
        --seed {seed} \\
        --out "$OUTDIR/$base/${{base}}_vina_out.pdbqt" \\
        > "$LOGDIR/${{base}}.log" 2>&1
done
""")
    script.chmod(0o755)
    return script


# ---------------------------------------------------------------------------
# Vina parsing and summary
# ---------------------------------------------------------------------------

MODE_LINE_RE = re.compile(r"^\s*(\d+)\s+(-?\d+(?:\.\d+)?)\s+(\d+(?:\.\d+)?)\s+(\d+(?:\.\d+)?)\s*$")


def parse_vina_log(log_path: Path) -> List[Dict[str, object]]:
    rows = []
    ligand = log_path.stem
    with log_path.open(errors="ignore") as handle:
        for line in handle:
            m = MODE_LINE_RE.match(line)
            if not m:
                continue
            rows.append({
                "Ligand": ligand,
                "LogFile": str(log_path),
                "Mode": int(m.group(1)),
                "Affinity": float(m.group(2)),
                "RMSD_lb": float(m.group(3)),
                "RMSD_ub": float(m.group(4)),
            })
    return rows


def trimmed_mean(values: Sequence[float], fraction: float = 0.10) -> float:
    vals = sorted([float(v) for v in values if np.isfinite(float(v))])
    if not vals:
        return np.nan
    k = int(math.floor(len(vals) * fraction))
    if len(vals) - 2*k <= 0:
        return float(np.mean(vals))
    return float(np.mean(vals[k:len(vals)-k]))


def parse_vina_results(log_dir: Path, outdir: Path, trim_fraction: float = 0.10) -> Dict[str, pd.DataFrame]:
    mkdir(outdir)
    logs = sorted(log_dir.glob("*.log"))
    long_rows = []
    for p in logs:
        long_rows.extend(parse_vina_log(p))
    long = pd.DataFrame(long_rows)
    if long.empty:
        long.to_csv(outdir / "vina_modes_long.tsv", sep="\t", index=False)
        return {"vina_modes_long": long, "vina_summary": pd.DataFrame()}
    summaries = []
    for ligand, g in long.groupby("Ligand"):
        g = g.sort_values("Mode")
        top = g.iloc[0]
        summaries.append({
            "Ligand": ligand,
            "Top1_Affinity": top["Affinity"],
            "Top1_RMSD_lb": top["RMSD_lb"],
            "Top1_RMSD_ub": top["RMSD_ub"],
            "Best_Affinity": g["Affinity"].min(),
            "Mean_Affinity": g["Affinity"].mean(),
            "TrimmedMean_Affinity": trimmed_mean(g["Affinity"].values, trim_fraction),
            "Mean_RMSD_lb": g["RMSD_lb"].mean(),
            "TrimmedMean_RMSD_lb": trimmed_mean(g["RMSD_lb"].values, trim_fraction),
            "NumModes": int(g.shape[0]),
        })
    summary = pd.DataFrame(summaries).sort_values(["Top1_Affinity", "Best_Affinity"], ascending=[True, True])
    long.to_csv(outdir / "vina_modes_long.tsv", sep="\t", index=False)
    summary.to_csv(outdir / "vina_summary.tsv", sep="\t", index=False)
    write_excel(outdir / "vina_refinement_summary.xlsx", {"modes_long": long, "summary": summary})
    plot_vina(summary, outdir / "plots")
    return {"vina_modes_long": long, "vina_summary": summary}


def plot_vina(summary: pd.DataFrame, plot_dir: Path):
    mkdir(plot_dir)
    if summary.empty:
        return
    data = summary.sort_values("Top1_Affinity").head(30)
    fig, ax = plt.subplots(figsize=(9, max(4, 0.28*len(data))), dpi=200)
    ax.barh(range(len(data)), data["Top1_Affinity"])
    ax.set_yticks(range(len(data)))
    ax.set_yticklabels(data["Ligand"], fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("Top Vina affinity (kcal/mol)")
    ax.set_title("Top refined Vina affinities")
    fig.tight_layout()
    fig.savefig(plot_dir / "top_vina_affinities.png")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Lightweight PDB parser for contact/CEH analysis
# ---------------------------------------------------------------------------

@dataclass
class AtomRec:
    record: str
    name: str
    resname: str
    chain: str
    resid: int
    x: float
    y: float
    z: float
    element: str
    line: str

    @property
    def coord(self) -> np.ndarray:
        return np.array([self.x, self.y, self.z], dtype=float)

    @property
    def residue_label(self) -> str:
        c = f"-{self.chain}" if self.chain.strip() else ""
        return f"{self.resname}{c}-{self.resid}"


def infer_element(atom_name: str, element_field: str = "") -> str:
    raw = (element_field or "").strip().upper()
    if raw:
        return raw[:2] if raw.startswith(("CL", "BR")) else raw[:1]
    letters = re.sub(r"[^A-Za-z]", "", atom_name).upper()
    if letters.startswith("CL"):
        return "CL"
    if letters.startswith("BR"):
        return "BR"
    return letters[:1] if letters else ""


def parse_pdb_atoms(path: Path) -> List[AtomRec]:
    atoms = []
    with path.open(errors="ignore") as handle:
        for line in handle:
            rec = line[:6].strip()
            if rec not in {"ATOM", "HETATM"}:
                continue
            try:
                atom = AtomRec(
                    record=rec,
                    name=line[12:16].strip(),
                    resname=line[17:20].strip(),
                    chain=line[21:22].strip(),
                    resid=int(line[22:26]),
                    x=float(line[30:38]),
                    y=float(line[38:46]),
                    z=float(line[46:54]),
                    element=infer_element(line[12:16], line[76:78] if len(line) >= 78 else ""),
                    line=line.rstrip("\n"),
                )
            except Exception:
                continue
            atoms.append(atom)
    return atoms


def is_water(atom: AtomRec) -> bool:
    return atom.resname.upper() in {"HOH", "WAT", "SOL", "TIP", "TIP3", "TIP3P"}


def split_protein_ligand(atoms: List[AtomRec], ligand_resname: Optional[str] = None) -> Tuple[List[AtomRec], List[AtomRec]]:
    protein = [a for a in atoms if a.record == "ATOM" and not is_water(a)]
    if ligand_resname:
        ligand = [a for a in atoms if a.resname.upper() == ligand_resname.upper()]
        return protein, ligand
    groups = defaultdict(list)
    for a in atoms:
        if a.record == "HETATM" and not is_water(a):
            groups[(a.chain, a.resid, a.resname)].append(a)
    if not groups:
        nonprotein = [a for a in atoms if a.record != "ATOM" and not is_water(a)]
        return protein, nonprotein
    key = max(groups, key=lambda k: sum(1 for a in groups[k] if a.element != "H"))
    return protein, groups[key]


def distance(a: AtomRec, b: AtomRec) -> float:
    return float(np.linalg.norm(a.coord - b.coord))


def angle_degrees(a: np.ndarray, b: np.ndarray) -> float:
    na = np.linalg.norm(a)
    nb = np.linalg.norm(b)
    if na == 0 or nb == 0:
        return np.nan
    return float(np.degrees(np.arccos(np.clip(np.dot(a, b) / (na * nb), -1, 1))))


def min_dist_atoms(a: Sequence[AtomRec], b: Sequence[AtomRec]) -> float:
    if not a or not b:
        return np.nan
    A = np.array([x.coord for x in a])
    B = np.array([x.coord for x in b])
    return float(np.linalg.norm(A[:, None, :] - B[None, :, :], axis=2).min())


def find_atoms(atoms: Sequence[AtomRec], resid: int, names: Sequence[str]) -> List[AtomRec]:
    names_u = {n.upper() for n in names}
    return [a for a in atoms if a.resid == resid and a.name.upper() in names_u]


def find_carbonyl_pairs(ligand: Sequence[AtomRec]) -> List[Tuple[AtomRec, AtomRec, float]]:
    carbons = [a for a in ligand if a.element == "C"]
    oxygens = [a for a in ligand if a.element == "O"]
    pairs = []
    for c in carbons:
        for o in oxygens:
            d = distance(c, o)
            if 1.15 <= d <= 1.35:
                pairs.append((c, o, d))
    return pairs


def triangular_score(value: float, ideal: float, tolerance: float) -> float:
    if not np.isfinite(value):
        return 0.0
    return max(0.0, 1.0 - abs(value - ideal) / tolerance)


def ceh_competence_score(ser_c: float, bd_angle: float, oxy_contacts: int, ser_his: float, acid_his: float) -> float:
    dist_s = triangular_score(ser_c, 3.0, 2.5)
    angle_s = triangular_score(bd_angle, 105.0, 75.0) if np.isfinite(bd_angle) and 90 <= bd_angle <= 120 else 0.0
    oxy_s = min(max(float(oxy_contacts), 0.0) / 2.0, 1.0)
    serhis_s = triangular_score(ser_his, 3.0, 2.0)
    acidhis_s = triangular_score(acid_his, 2.8, 1.5)
    triad_s = 0.6*serhis_s + 0.4*acidhis_s
    return 0.30*dist_s + 0.20*angle_s + 0.30*oxy_s + 0.20*triad_s


def generic_contact_analysis(pdb_path: Path, ligand_resname: Optional[str], cutoff: float = 4.5) -> Tuple[pd.DataFrame, pd.DataFrame]:
    atoms = parse_pdb_atoms(pdb_path)
    protein, ligand = split_protein_ligand(atoms, ligand_resname)
    rows = []
    if not protein or not ligand:
        return pd.DataFrame(), pd.DataFrame()
    res_groups = defaultdict(list)
    for a in protein:
        res_groups[(a.chain, a.resid, a.resname)].append(a)
    for key, ratoms in res_groups.items():
        md = min_dist_atoms(ratoms, ligand)
        if np.isfinite(md) and md <= cutoff:
            chain, resid, resname = key
            rows.append({"Pose": pdb_path.name, "Residue": f"{resname}-{chain}-{resid}" if chain else f"{resname}-{resid}", "MinDistance": md})
    df = pd.DataFrame(rows)
    if df.empty:
        return df, pd.DataFrame()
    freq = df.groupby("Residue").agg(ContactCount=("Pose", "count"), MinDistance=("MinDistance", "min"), MeanDistance=("MinDistance", "mean")).reset_index()
    freq["PoseFraction"] = 1.0
    return df, freq


def analyze_ceh_pose(
    pdb_path: Path,
    ligand_resname: Optional[str] = None,
    ser_resid: int = 207,
    his_resid: int = 508,
    acid_resid: int = 388,
    oxyanion_resids: Sequence[int] = (127, 128, 129),
    contact_cutoff: float = 4.5,
) -> Optional[Dict[str, object]]:
    atoms = parse_pdb_atoms(pdb_path)
    protein, ligand = split_protein_ligand(atoms, ligand_resname)
    if not protein or not ligand:
        return None
    ser_og = find_atoms(protein, ser_resid, ["OG"])
    his_ne2 = find_atoms(protein, his_resid, ["NE2"])
    his_nd1 = find_atoms(protein, his_resid, ["ND1"])
    acid_oe = find_atoms(protein, acid_resid, ["OE1", "OE2", "OD1", "OD2"])
    if not ser_og:
        return None
    ser_his = min_dist_atoms(ser_og, his_ne2)
    acid_his = min_dist_atoms(acid_oe, his_nd1)
    pairs = find_carbonyl_pairs(ligand)
    if not pairs:
        return None
    oxy_n = []
    for r in oxyanion_resids:
        oxy_n.extend(find_atoms(protein, r, ["N"]))
    candidates = []
    for c, o, co_len in pairs:
        ser_dist = min(distance(s, c) for s in ser_og)
        best_ser = min(ser_og, key=lambda s: distance(s, c))
        bd = angle_degrees(best_ser.coord - c.coord, o.coord - c.coord)
        oxy_dists = [distance(o, n) for n in oxy_n]
        oxy_contacts = sum(1 for d in oxy_dists if 2.5 <= d <= 3.5)
        score = ceh_competence_score(ser_dist, bd, oxy_contacts, ser_his, acid_his)
        candidates.append((score, c, o, co_len, ser_dist, bd, oxy_contacts, min(oxy_dists) if oxy_dists else np.nan))
    best = max(candidates, key=lambda x: x[0])
    score, c, o, co_len, ser_dist, bd, oxy_contacts, oxy_min = best
    nac = bool(ser_dist <= 4.0 and np.isfinite(bd) and 90 <= bd <= 120 and oxy_contacts >= 2 and ser_his <= 3.5 and acid_his <= 3.5)
    if nac:
        state = "catalytically-competent"
    elif ser_dist <= 4.0 and ser_his <= 3.5:
        state = "near-reactive"
    elif ser_dist <= contact_cutoff:
        state = "binding"
    else:
        state = "distant"
    contacts, freq = generic_contact_analysis(pdb_path, ligand_resname, contact_cutoff)
    contact_list = [] if contacts.empty else sorted(contacts["Residue"].unique().tolist())
    return {
        "Pose": pdb_path.name,
        "Path": str(pdb_path),
        "LigandResname": ligand[0].resname if ligand else "",
        "CarbonylCandidates": len(candidates),
        "CarbonylCName": c.name,
        "CarbonylOName": o.name,
        "CarbonylBondLength": co_len,
        "SerCarbonylDistance": ser_dist,
        "BD_Angle": bd,
        "SerHisDistance": ser_his,
        "AcidHisDistance": acid_his,
        "OxyanionContacts": oxy_contacts,
        "OxyanionMinDistance": oxy_min,
        "NAC": nac,
        "State": state,
        "CompetenceScore": score,
        "ContactCount": len(contact_list),
        "Contacts": ";".join(contact_list),
    }


def analyze_complexes(
    complex_dir: Path,
    outdir: Path,
    engine: str = "general",
    ligand_resname: Optional[str] = None,
    recursive: bool = False,
    pattern: str = "*.pdb",
    contact_cutoff: float = 4.5,
    ser_resid: int = 207,
    his_resid: int = 508,
    acid_resid: int = 388,
) -> Dict[str, pd.DataFrame]:
    mkdir(outdir)
    iterator = complex_dir.rglob(pattern) if recursive else complex_dir.glob(pattern)
    files = sorted(p for p in iterator if p.is_file())
    contact_rows = []
    ceh_rows = []
    skipped = []
    for p in files:
        contacts, freq = generic_contact_analysis(p, ligand_resname, contact_cutoff)
        if not contacts.empty:
            contact_rows.append(contacts)
        if engine.lower() in {"ceh", "hydrolase", "ceh-like", "alpha-beta-hydrolase"}:
            row = analyze_ceh_pose(p, ligand_resname, ser_resid, his_resid, acid_resid, contact_cutoff=contact_cutoff)
            if row is None:
                skipped.append({"Pose": p.name, "Reason": "CEH geometry not parsed or ligand/carbonyl not found"})
            else:
                ceh_rows.append(row)
    contacts_all = pd.concat(contact_rows, ignore_index=True) if contact_rows else pd.DataFrame()
    if not contacts_all.empty:
        freq = contacts_all.groupby("Residue").agg(Frequency=("Pose", "nunique"), MinDistance=("MinDistance", "min"), MeanDistance=("MinDistance", "mean")).reset_index()
        denom = max(1, len(files))
        freq["PoseFraction"] = freq["Frequency"] / denom
        freq = freq.sort_values(["PoseFraction", "MinDistance"], ascending=[False, True])
    else:
        freq = pd.DataFrame()
    ceh = pd.DataFrame(ceh_rows)
    if not ceh.empty:
        ceh = ceh.sort_values(["CompetenceScore", "NAC", "SerCarbonylDistance", "BD_Angle"], ascending=[False, False, True, False]).reset_index(drop=True)
        ceh.insert(0, "Rank", np.arange(1, len(ceh)+1))
    skipped_df = pd.DataFrame(skipped)
    contacts_all.to_csv(outdir / "pose_contacts_long.tsv", sep="\t", index=False)
    freq.to_csv(outdir / "contact_frequency.tsv", sep="\t", index=False)
    ceh.to_csv(outdir / "ceh_catalytic_geometry.tsv", sep="\t", index=False)
    skipped_df.to_csv(outdir / "skipped_complexes.tsv", sep="\t", index=False)
    write_excel(outdir / "post_refinement_analysis.xlsx", {"contacts_long": contacts_all, "contact_frequency": freq, "ceh_catalytic_geometry": ceh, "skipped": skipped_df})
    plot_postrefinement(freq, ceh, outdir / "plots")
    return {"contacts_long": contacts_all, "contact_frequency": freq, "ceh_catalytic_geometry": ceh, "skipped": skipped_df}


def plot_postrefinement(freq: pd.DataFrame, ceh: pd.DataFrame, plot_dir: Path):
    mkdir(plot_dir)
    if freq is not None and not freq.empty:
        data = freq.sort_values("PoseFraction", ascending=False).head(30).sort_values("PoseFraction")
        fig, ax = plt.subplots(figsize=(8, max(4, 0.27*len(data))), dpi=200)
        ax.barh(range(len(data)), data["PoseFraction"])
        ax.set_yticks(range(len(data)))
        ax.set_yticklabels(data["Residue"], fontsize=7)
        ax.set_xlabel("Pose fraction")
        ax.set_title("Protein-ligand contact hotspots")
        fig.tight_layout()
        fig.savefig(plot_dir / "contact_hotspots.png")
        plt.close(fig)
    if ceh is not None and not ceh.empty:
        top = ceh.head(30)
        fig, ax = plt.subplots(figsize=(8, 5), dpi=200)
        ax.scatter(top["SerCarbonylDistance"], top["BD_Angle"], s=35)
        ax.axvspan(0, 4.0, alpha=0.08)
        ax.axhspan(90, 120, alpha=0.08)
        ax.set_xlabel("Ser OG to ligand carbonyl C distance (Å)")
        ax.set_ylabel("Bürgi-Dunitz-like angle (°)")
        ax.set_title("CEH catalytic geometry of refined poses")
        fig.tight_layout()
        fig.savefig(plot_dir / "ceh_catalytic_geometry_scatter.png")
        plt.close(fig)
        fig, ax = plt.subplots(figsize=(9, max(4, 0.26*len(top))), dpi=200)
        ax.barh(range(len(top)), top["CompetenceScore"])
        ax.set_yticks(range(len(top)))
        ax.set_yticklabels(top["Pose"], fontsize=7)
        ax.invert_yaxis()
        ax.set_xlabel("Catalytic competence score")
        ax.set_title("Top CEH-like refined poses")
        fig.tight_layout()
        fig.savefig(plot_dir / "ceh_top_competence_scores.png")
        plt.close(fig)


# ---------------------------------------------------------------------------
# Final substrate specificity integration
# ---------------------------------------------------------------------------


def final_screening_report(
    outdir: Path,
    diffdock_tables: Optional[Dict[str, pd.DataFrame]] = None,
    vina_tables: Optional[Dict[str, pd.DataFrame]] = None,
    post_tables: Optional[Dict[str, pd.DataFrame]] = None,
    engine: str = "general",
) -> pd.DataFrame:
    mkdir(outdir)
    score = pd.DataFrame()
    if vina_tables and not vina_tables.get("vina_summary", pd.DataFrame()).empty:
        score = vina_tables["vina_summary"].copy()
        score["LigandKey"] = score["Ligand"].astype(str).map(safe_stem)
        score["vina_affinity_scaled"] = robust_scale(score["Top1_Affinity"], "lower_better")
        score["vina_rmsd_scaled"] = robust_scale(score["Top1_RMSD_lb"], "lower_better")
    if diffdock_tables and not diffdock_tables.get("top_poses", pd.DataFrame()).empty:
        dd = diffdock_tables["top_poses"].copy()
        dd["LigandKey"] = dd["Ligand"].astype(str).map(safe_stem)
        dd_best = dd.sort_values("global_pose_score", ascending=False).groupby("LigandKey").first().reset_index()
        cols = ["LigandKey", "Receptor", "Ligand", "Prediction", "ClusterID", "global_pose_score", "SMINA_Minimized_Affinity", "SMINA_Minimized_RMSD"]
        dd_best = dd_best[[c for c in cols if c in dd_best.columns]]
        if score.empty:
            score = dd_best.copy()
        else:
            score = score.merge(dd_best, on="LigandKey", how="outer", suffixes=("_vina", "_diffdock"))
    if post_tables:
        if engine.lower() in {"ceh", "hydrolase", "ceh-like", "alpha-beta-hydrolase"} and not post_tables.get("ceh_catalytic_geometry", pd.DataFrame()).empty:
            ceh = post_tables["ceh_catalytic_geometry"].copy()
            ceh["LigandKey"] = ceh["Pose"].astype(str).map(lambda x: safe_stem(x).split("__")[0])
            # Use best catalytic pose per ligand-like prefix; if names do not encode ligand, still include all pose-level rows.
            ceh_best = ceh.sort_values("CompetenceScore", ascending=False).groupby("LigandKey").first().reset_index()
            keep = ["LigandKey", "Pose", "State", "NAC", "CompetenceScore", "SerCarbonylDistance", "BD_Angle", "SerHisDistance", "AcidHisDistance", "OxyanionContacts"]
            ceh_best = ceh_best[[c for c in keep if c in ceh_best.columns]]
            if score.empty:
                score = ceh_best.copy()
            else:
                score = score.merge(ceh_best, on="LigandKey", how="outer")
        elif not post_tables.get("contact_frequency", pd.DataFrame()).empty:
            pass
    if score.empty:
        score.to_csv(outdir / "substrate_specificity_screen.tsv", sep="\t", index=False)
        return score
    components = pd.DataFrame(index=score.index)
    if "vina_affinity_scaled" in score.columns:
        components["vina_affinity"] = score["vina_affinity_scaled"]
    if "vina_rmsd_scaled" in score.columns:
        components["vina_rmsd"] = score["vina_rmsd_scaled"]
    if "global_pose_score" in score.columns:
        components["diffdock_pose"] = pd.to_numeric(score["global_pose_score"], errors="coerce")
    if "CompetenceScore" in score.columns:
        components["ceh_catalytic"] = pd.to_numeric(score["CompetenceScore"], errors="coerce")
    if engine.lower() in {"ceh", "hydrolase", "ceh-like", "alpha-beta-hydrolase"}:
        weights = {"ceh_catalytic": 0.45, "vina_affinity": 0.30, "diffdock_pose": 0.15, "vina_rmsd": 0.10}
    else:
        weights = {"vina_affinity": 0.45, "diffdock_pose": 0.30, "vina_rmsd": 0.15, "ceh_catalytic": 0.10}
    for c in weights:
        if c not in components:
            components[c] = np.nan
    components = components.apply(pd.to_numeric, errors="coerce")
    for c in components:
        components[c] = components[c].fillna(components[c].median() if components[c].notna().any() else 0.5)
    score["SubstrateSpecificityScore"] = weighted_sum(components, weights)
    score = score.sort_values("SubstrateSpecificityScore", ascending=False)
    score.to_csv(outdir / "substrate_specificity_screen.tsv", sep="\t", index=False)
    write_excel(outdir / "substrate_specificity_screen.xlsx", {"screen": score})
    if not score.empty:
        top = score.head(25)
        label_col = "Ligand" if "Ligand" in top.columns else "LigandKey"
        fig, ax = plt.subplots(figsize=(9, max(4, 0.3*len(top))), dpi=200)
        ax.barh(range(len(top)), top["SubstrateSpecificityScore"])
        ax.set_yticks(range(len(top)))
        ax.set_yticklabels(top[label_col].astype(str), fontsize=7)
        ax.invert_yaxis()
        ax.set_xlabel("Integrated specificity score")
        ax.set_title("Top substrate candidates")
        fig.tight_layout()
        mkdir(outdir / "plots")
        fig.savefig(outdir / "plots" / "top_substrate_specificity_scores.png")
        plt.close(fig)
    return score



# ---------------------------------------------------------------------------
# General rigid-receptor Vina pose splitting, complex construction, contacts
# ---------------------------------------------------------------------------

def _format_template(template: str, protein: str, ligand: str) -> str:
    return template.format(protein=protein, ligand=ligand, Protein=protein, Ligand=ligand)


def _atom_lines_for_complex(path: Path, ligand: bool = False) -> List[str]:
    lines = []
    for line in path.read_text(errors="ignore").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            if ligand and line.startswith("ATOM"):
                line = "HETATM" + line[6:]
            lines.append(line.rstrip())
    return lines


def construct_receptor_ligand_complex(receptor_pdb: Path, ligand_pose_pdb: Path, out_pdb: Path) -> Path:
    """Create a simple receptor-ligand complex PDB from receptor + ligand pose."""
    mkdir(out_pdb.parent)
    receptor_lines = _atom_lines_for_complex(receptor_pdb, ligand=False)
    ligand_lines = _atom_lines_for_complex(ligand_pose_pdb, ligand=True)
    with out_pdb.open("w") as out:
        for line in receptor_lines:
            out.write(line + "\n")
        if receptor_lines:
            out.write("TER\n")
        for line in ligand_lines:
            out.write(line + "\n")
        out.write("END\n")
    return out_pdb


def split_vina_pdbqt_modes(vina_pdbqt: Path, outdir: Path, keep_existing: bool = True) -> List[Path]:
    """Split a multi-MODEL Vina PDBQT into individual PDB pose files.

    OpenBabel is used when available, matching the batch_contact_cd36rp.py
    workflow. If OpenBabel is unavailable, pre-existing *_pose*.pdb files in
    the output directory are returned so a previous split can still be analyzed.
    """
    mkdir(outdir)
    if keep_existing:
        existing = sorted(outdir.glob(f"{safe_stem(vina_pdbqt.name)}_pose*.pdb"))
        if existing:
            return existing
    obabel = shutil.which("obabel")
    if obabel is None:
        return sorted(outdir.glob("*_pose*.pdb"))
    before = set(outdir.glob("*.pdb"))
    tmp_base = outdir / f"{safe_stem(vina_pdbqt.name)}_poses.pdb"
    subprocess.run([obabel, str(vina_pdbqt), "-O", str(tmp_base), "-m"], check=True)
    after = sorted(set(outdir.glob("*.pdb")) - before)
    if not after:
        after = sorted(outdir.glob(f"{tmp_base.stem}*.pdb"))
    renamed = []
    for i, path in enumerate(sorted(after), start=1):
        new_path = outdir / f"{safe_stem(vina_pdbqt.name)}_pose{i}.pdb"
        if path != new_path:
            if new_path.exists():
                new_path.unlink()
            path.rename(new_path)
        if new_path.stat().st_size > 100:
            renamed.append(new_path)
    return renamed


def batch_contact_general(
    project_root: Path,
    outdir: Path,
    proteins: Sequence[str],
    ligands: Sequence[str],
    receptor_template: str = "{protein}/receptor.pdb",
    workdir_template: str = "{protein}/{ligand}",
    vina_glob: str = "*.pdbqt",
    contact_cutoff: float = 4.5,
    freq_threshold: float = 0.5,
    ligand_resname: Optional[str] = None,
) -> Dict[str, pd.DataFrame]:
    """General/CD36RP-style automated Vina post-processing.

    For every protein-ligand pair, this module splits Vina PDBQT modes,
    builds receptor-ligand complex PDB files, and computes residue contact
    frequencies. It is a parameterized version of batch_contact_cd36rp.py.
    """
    mkdir(outdir)
    complexes_dir = mkdir(outdir / "complexes")
    all_contacts = []
    stats_rows = []
    skipped = []

    for protein in proteins:
        receptor = project_root / _format_template(receptor_template, protein, "")
        if not receptor.is_file():
            skipped.append({"Protein": protein, "Ligand": "*", "Reason": f"missing receptor: {receptor}"})
            continue
        for ligand in ligands:
            workdir = project_root / _format_template(workdir_template, protein, ligand)
            if not workdir.is_dir():
                skipped.append({"Protein": protein, "Ligand": ligand, "Reason": f"missing workdir: {workdir}"})
                continue
            vina_files = sorted(workdir.glob(vina_glob))
            if not vina_files:
                skipped.append({"Protein": protein, "Ligand": ligand, "Reason": f"no Vina PDBQT matched {vina_glob}"})
                continue
            pair_pose_dir = mkdir(outdir / "poses" / protein / ligand)
            pair_complex_dir = mkdir(complexes_dir / protein / ligand)
            pair_contacts = []
            for vina_file in vina_files:
                poses = split_vina_pdbqt_modes(vina_file, pair_pose_dir)
                if not poses:
                    skipped.append({"Protein": protein, "Ligand": ligand, "Reason": f"no split poses for {vina_file}"})
                    continue
                for pose in poses:
                    complex_path = pair_complex_dir / f"{safe_stem(pose.name)}_complex.pdb"
                    construct_receptor_ligand_complex(receptor, pose, complex_path)
                    contacts, _ = generic_contact_analysis(complex_path, ligand_resname, contact_cutoff)
                    if contacts.empty:
                        continue
                    contacts["Protein"] = protein
                    contacts["Ligand"] = ligand
                    contacts["Complex"] = complex_path.name
                    pair_contacts.append(contacts)
            if pair_contacts:
                pair_df = pd.concat(pair_contacts, ignore_index=True)
                pair_df.to_csv(outdir / f"{protein}_{ligand}_contacts_all.tsv", sep="\t", index=False)
                n_poses = max(1, pair_df["Pose"].nunique())
                freq = (pair_df.groupby("Residue")
                        .agg(ContactCount=("Pose", "nunique"), MinDistance=("MinDistance", "min"), MeanDistance=("MinDistance", "mean"))
                        .reset_index())
                freq["Frequency"] = freq["ContactCount"] / float(n_poses)
                freq_f = freq[freq["Frequency"] > freq_threshold].sort_values("Frequency", ascending=False)
                freq.to_csv(outdir / f"{protein}_{ligand}_contacts_frequency.tsv", sep="\t", index=False)
                freq_f.to_csv(outdir / f"{protein}_{ligand}_contacts_filtered.tsv", sep="\t", index=False)
                stats_rows.append({
                    "Protein": protein,
                    "Ligand": ligand,
                    "PoseCount": n_poses,
                    "MeanFrequency": float(freq["Frequency"].mean()) if not freq.empty else np.nan,
                    "HighFreqResidues": int((freq["Frequency"] >= 0.8).sum()) if not freq.empty else 0,
                    "Entropy": float(-np.sum((freq["Frequency"]/freq["Frequency"].sum()) * np.log((freq["Frequency"]/freq["Frequency"].sum()) + 1e-9))) if not freq.empty and freq["Frequency"].sum() > 0 else np.nan,
                })
                all_contacts.append(pair_df)

    contacts_all = pd.concat(all_contacts, ignore_index=True) if all_contacts else pd.DataFrame()
    stats = pd.DataFrame(stats_rows)
    skipped_df = pd.DataFrame(skipped)
    contacts_all.to_csv(outdir / "ALL_contact_summary.tsv", sep="\t", index=False)
    stats.to_csv(outdir / "binding_statistics.tsv", sep="\t", index=False)
    skipped_df.to_csv(outdir / "skipped_pairs.tsv", sep="\t", index=False)
    write_excel(outdir / "general_vina_contact_analysis.xlsx", {"contacts": contacts_all, "statistics": stats, "skipped": skipped_df})
    plot_general_contact_summary(contacts_all, stats, outdir / "plots", freq_threshold)
    return {"contacts": contacts_all, "statistics": stats, "skipped": skipped_df}


def plot_general_contact_summary(contacts_all: pd.DataFrame, stats: pd.DataFrame, plot_dir: Path, freq_threshold: float):
    mkdir(plot_dir)
    if contacts_all is not None and not contacts_all.empty:
        for protein in sorted(contacts_all["Protein"].dropna().unique()):
            sub = contacts_all[contacts_all["Protein"] == protein]
            if sub.empty:
                continue
            pose_counts = sub.groupby("Ligand")["Pose"].nunique().to_dict()
            freq_rows = []
            for (ligand, residue), g in sub.groupby(["Ligand", "Residue"]):
                denom = max(1, int(pose_counts.get(ligand, 1)))
                freq_rows.append({"Ligand": ligand, "Residue": residue, "Frequency": g["Pose"].nunique() / denom})
            freq = pd.DataFrame(freq_rows)
            freq_f = freq[freq["Frequency"] > freq_threshold]
            if not freq_f.empty:
                pivot = freq_f.pivot_table(index="Residue", columns="Ligand", values="Frequency", fill_value=0)
                fig, ax = plt.subplots(figsize=(max(6, 1.2*pivot.shape[1]), max(5, 0.25*pivot.shape[0])), dpi=200)
                im = ax.imshow(pivot.values, aspect="auto", interpolation="nearest")
                ax.set_xticks(range(pivot.shape[1])); ax.set_xticklabels(pivot.columns, rotation=45, ha="right")
                ax.set_yticks(range(pivot.shape[0])); ax.set_yticklabels(pivot.index, fontsize=7)
                ax.set_title(f"{protein} contact frequency heatmap (> {freq_threshold})")
                fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Pose fraction")
                fig.tight_layout(); fig.savefig(plot_dir / f"{protein}_contact_heatmap_filtered.png"); plt.close(fig)
    if stats is not None and not stats.empty:
        fig, ax = plt.subplots(figsize=(8, max(4, 0.3*len(stats))), dpi=200)
        data = stats.sort_values("MeanFrequency", ascending=False)
        labels = data.apply(lambda r: f"{r['Protein']}|{r['Ligand']}", axis=1)
        ax.barh(range(len(data)), data["MeanFrequency"])
        ax.set_yticks(range(len(data))); ax.set_yticklabels(labels, fontsize=7)
        ax.invert_yaxis(); ax.set_xlabel("Mean contact frequency"); ax.set_title("General docking contact summary")
        fig.tight_layout(); fig.savefig(plot_dir / "binding_statistics_mean_frequency.png"); plt.close(fig)


# ---------------------------------------------------------------------------
# CEH flexible Vina ensemble wrapper for flexdock_v3
# ---------------------------------------------------------------------------

def generate_flexdock_v3_script(
    outdir: Path,
    flex_out: str,
    receptor: str,
    ligand_sdf: str,
    results_dir: str = "flexdock_v3_results",
    expected_ligand_atoms: int = 47,
    his_cutoff: float = 0.75,
    ligand_cutoff: float = 2.0,
) -> Path:
    """Write a shell script that runs the packaged flexdock_v3 CEH workflow."""
    mkdir(outdir)
    script = outdir / "run_flexdock_v3_ceh.sh"
    script.write_text(f"""#!/bin/bash
set -euo pipefail
# CEH-like flexible Vina ensemble analysis: parse modes, cluster HIS508,
# cluster ligand poses, select representatives, and build competent complexes.
PYTHONPATH="$(pwd):$PYTHONPATH" python -m flexdock_v3.run_pipeline \\
  --flex-out "{flex_out}" \\
  --receptor "{receptor}" \\
  --ligand-sdf "{ligand_sdf}" \\
  --results-dir "{results_dir}" \\
  --expected-ligand-atoms {expected_ligand_atoms} \\
  --his-cutoff {his_cutoff} \\
  --ligand-cutoff {ligand_cutoff}
""")
    script.chmod(0o755)
    return script


def run_flexdock_v3_wrapper(
    flex_out: Path,
    receptor: Path,
    ligand_sdf: Path,
    results_dir: Path,
    expected_ligand_atoms: int = 47,
    his_cutoff: float = 0.75,
    ligand_cutoff: float = 2.0,
) -> None:
    """Run flexdock_v3 if it is importable from the current Python path."""
    from flexdock_v3.run_pipeline import build_parser, run_pipeline  # type: ignore
    args = build_parser().parse_args([
        "--flex-out", str(flex_out),
        "--receptor", str(receptor),
        "--ligand-sdf", str(ligand_sdf),
        "--results-dir", str(results_dir),
        "--expected-ligand-atoms", str(expected_ligand_atoms),
        "--his-cutoff", str(his_cutoff),
        "--ligand-cutoff", str(ligand_cutoff),
    ])
    run_pipeline(args)

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def add_common_workflow_args(p):
    p.add_argument("--engine", default="general", choices=["general", "ceh", "hydrolase"], help="Post-refinement analysis engine.")
    p.add_argument("--outdir", default="ensemble_docking_results", type=Path)
    p.add_argument("--ligand-resname", default=None)
    p.add_argument("--contact-cutoff", type=float, default=4.5)
    p.add_argument("--ser-resid", type=int, default=207)
    p.add_argument("--his-resid", type=int, default=508)
    p.add_argument("--acid-resid", type=int, default=388)


def build_parser():
    ap = argparse.ArgumentParser(description="Multi-stage ensemble docking and substrate-specificity screening pipeline")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("cluster-diffdock", help="Cluster/select DiffDock-L + SMINA poses")
    p.add_argument("--csv", required=True, type=Path)
    p.add_argument("--sdf-dir", required=True, type=Path)
    p.add_argument("--outdir", default="diffdock_selection", type=Path)
    p.add_argument("--rmsd-cutoff", type=float, default=2.0)
    p.add_argument("--top-k-clusters", type=int, default=3)
    p.add_argument("--top-n-poses", type=int, default=10)
    p.add_argument("--min-aff-cutoff", type=float, default=None)
    p.add_argument("--max-rmsd-cutoff", type=float, default=None)
    p.add_argument("--confidence-mode", choices=["auto", "higher_better", "lower_better"], default="auto")
    p.add_argument("--drop-zero-scores", action="store_true")
    p.add_argument("--no-largest-fragment", action="store_true")

    p = sub.add_parser("make-scripts", help="Generate ligand/receptor/Vina shell scripts")
    p.add_argument("--outdir", default="docking_scripts", type=Path)
    p.add_argument("--receptor-pdb", required=True)
    p.add_argument("--box-center", nargs=3, type=float, default=None)
    p.add_argument("--box-size", nargs=3, type=float, default=None)
    p.add_argument("--box-enveloping", default=None, help="Ligand PDBQT used to define an enveloping rigid-receptor box.")
    p.add_argument("--padding", type=float, default=4.0, help="Padding for --box-enveloping mode.")
    p.add_argument("--receptor-mode", choices=["auto", "general-rigid", "explicit-box"], default="auto")
    p.add_argument("--flex-residues", default="")
    p.add_argument("--receptor-name", default="receptor")
    p.add_argument("--selected-sdf-dir", default="selected_poses")
    p.add_argument("--ligand-dir", default="ligand_pdbqt")
    p.add_argument("--receptor-rigid", default="receptor_grid/receptor_rigid.pdbqt")
    p.add_argument("--receptor-flex", default=None)
    p.add_argument("--box-file", default="receptor_grid/receptor_box.txt")
    p.add_argument("--exhaustiveness", type=int, default=64)
    p.add_argument("--num-modes", type=int, default=50)

    p = sub.add_parser("parse-vina", help="Parse AutoDock Vina refinement logs")
    p.add_argument("--log-dir", required=True, type=Path)
    p.add_argument("--outdir", default="vina_parsed", type=Path)
    p.add_argument("--trim-fraction", type=float, default=0.10)

    p = sub.add_parser("analyze-complexes", help="Analyze refined receptor-ligand complex PDBs")
    p.add_argument("--complex-dir", required=True, type=Path)
    p.add_argument("--recursive", action="store_true")
    p.add_argument("--pattern", default="*.pdb")
    add_common_workflow_args(p)

    p = sub.add_parser("batch-contact-general", help="Split rigid Vina modes, build complexes, and run CD36RP/general contact analysis")
    p.add_argument("--project-root", default=".", type=Path)
    p.add_argument("--proteins", nargs="+", required=True)
    p.add_argument("--ligands", nargs="+", required=True)
    p.add_argument("--receptor-template", default="{protein}/receptor.pdb")
    p.add_argument("--workdir-template", default="{protein}/{ligand}")
    p.add_argument("--vina-glob", default="*.pdbqt")
    p.add_argument("--freq-threshold", type=float, default=0.5)
    p.add_argument("--outdir", default="general_contact_analysis", type=Path)
    p.add_argument("--ligand-resname", default=None)
    p.add_argument("--contact-cutoff", type=float, default=4.5)

    p = sub.add_parser("make-flexdock-ceh-script", help="Generate a flexdock_v3 CEH flexible-output analysis script")
    p.add_argument("--outdir", default="ceh_flexdock_scripts", type=Path)
    p.add_argument("--flex-out", required=True)
    p.add_argument("--receptor", required=True)
    p.add_argument("--ligand-sdf", required=True)
    p.add_argument("--results-dir", default="flexdock_v3_results")
    p.add_argument("--expected-ligand-atoms", type=int, default=47)
    p.add_argument("--his-cutoff", type=float, default=0.75)
    p.add_argument("--ligand-cutoff", type=float, default=2.0)

    p = sub.add_parser("run-flexdock-ceh", help="Run flexdock_v3 CEH ensemble parsing/clustering/complex-building if the module is available")
    p.add_argument("--flex-out", required=True, type=Path)
    p.add_argument("--receptor", required=True, type=Path)
    p.add_argument("--ligand-sdf", required=True, type=Path)
    p.add_argument("--results-dir", default="flexdock_v3_results", type=Path)
    p.add_argument("--expected-ligand-atoms", type=int, default=47)
    p.add_argument("--his-cutoff", type=float, default=0.75)
    p.add_argument("--ligand-cutoff", type=float, default=2.0)

    p = sub.add_parser("workflow", help="Run all available stages and generate scripts for external tools")
    p.add_argument("--diffdock-csv", type=Path, default=None)
    p.add_argument("--sdf-dir", type=Path, default=None)
    p.add_argument("--receptor-pdb", default=None)
    p.add_argument("--box-center", nargs=3, type=float, default=None)
    p.add_argument("--box-size", nargs=3, type=float, default=None)
    p.add_argument("--box-enveloping", default=None, help="Ligand PDBQT used to define an enveloping rigid-receptor box.")
    p.add_argument("--padding", type=float, default=4.0, help="Padding for --box-enveloping mode.")
    p.add_argument("--receptor-mode", choices=["auto", "general-rigid", "explicit-box"], default="auto")
    p.add_argument("--flex-residues", default="")
    p.add_argument("--vina-log-dir", type=Path, default=None)
    p.add_argument("--complex-dir", type=Path, default=None)
    p.add_argument("--recursive-complexes", action="store_true")
    p.add_argument("--top-k-clusters", type=int, default=3)
    p.add_argument("--top-n-poses", type=int, default=10)
    p.add_argument("--rmsd-cutoff", type=float, default=2.0)
    add_common_workflow_args(p)
    return ap


def command_main(args):
    if args.cmd == "cluster-diffdock":
        cluster_diffdock(args.csv, args.sdf_dir, args.outdir, args.rmsd_cutoff, args.top_k_clusters, args.top_n_poses, args.min_aff_cutoff, args.max_rmsd_cutoff, args.confidence_mode, args.drop_zero_scores, not args.no_largest_fragment)
        print(f"Wrote DiffDock-L pose selection outputs to {args.outdir}")
    elif args.cmd == "make-scripts":
        lp = generate_ligand_preparation_script(args.outdir, args.selected_sdf_dir)
        rp = generate_receptor_preparation_script(args.outdir, args.receptor_pdb, tuple(args.box_center) if args.box_center else None, tuple(args.box_size) if args.box_size else None, args.flex_residues, args.receptor_name, args.box_enveloping, args.padding, args.receptor_mode)
        vp = generate_vina_refinement_script(args.outdir, args.receptor_rigid, args.receptor_flex, args.box_file, args.ligand_dir, exhaustiveness=args.exhaustiveness, num_modes=args.num_modes)
        print(f"Wrote scripts:\n  {lp}\n  {rp}\n  {vp}")
    elif args.cmd == "parse-vina":
        parse_vina_results(args.log_dir, args.outdir, args.trim_fraction)
        print(f"Wrote Vina summaries to {args.outdir}")
    elif args.cmd == "analyze-complexes":
        analyze_complexes(args.complex_dir, args.outdir, args.engine, args.ligand_resname, args.recursive, args.pattern, args.contact_cutoff, args.ser_resid, args.his_resid, args.acid_resid)
        print(f"Wrote post-refinement analysis to {args.outdir}")
    elif args.cmd == "batch-contact-general":
        batch_contact_general(args.project_root, args.outdir, args.proteins, args.ligands, args.receptor_template, args.workdir_template, args.vina_glob, args.contact_cutoff, args.freq_threshold, args.ligand_resname)
        print(f"Wrote general/CD36RP contact analysis to {args.outdir}")
    elif args.cmd == "make-flexdock-ceh-script":
        sp = generate_flexdock_v3_script(args.outdir, args.flex_out, args.receptor, args.ligand_sdf, args.results_dir, args.expected_ligand_atoms, args.his_cutoff, args.ligand_cutoff)
        print(f"Wrote CEH flexdock_v3 script: {sp}")
    elif args.cmd == "run-flexdock-ceh":
        run_flexdock_v3_wrapper(args.flex_out, args.receptor, args.ligand_sdf, args.results_dir, args.expected_ligand_atoms, args.his_cutoff, args.ligand_cutoff)
        print(f"Wrote flexdock_v3 CEH outputs to {args.results_dir}")
    elif args.cmd == "workflow":
        outdir = mkdir(args.outdir)
        diffdock_tables = None
        vina_tables = None
        post_tables = None
        if args.diffdock_csv and args.sdf_dir:
            diffdock_tables = cluster_diffdock(args.diffdock_csv, args.sdf_dir, outdir / "01_diffdock_pose_selection", args.rmsd_cutoff, args.top_k_clusters, args.top_n_poses)
        scripts_dir = mkdir(outdir / "02_scripts")
        generate_ligand_preparation_script(scripts_dir, "../01_diffdock_pose_selection/selected_poses")
        if args.receptor_pdb:
            generate_receptor_preparation_script(scripts_dir, args.receptor_pdb, tuple(args.box_center) if args.box_center else None, tuple(args.box_size) if args.box_size else None, args.flex_residues, "receptor", args.box_enveloping, args.padding, args.receptor_mode)
            generate_vina_refinement_script(scripts_dir, receptor_flex="receptor_grid/receptor_flex.pdbqt" if args.flex_residues else None)
        if args.vina_log_dir and args.vina_log_dir.exists():
            vina_tables = parse_vina_results(args.vina_log_dir, outdir / "03_vina_refinement")
        if args.complex_dir and args.complex_dir.exists():
            post_tables = analyze_complexes(args.complex_dir, outdir / "04_post_refinement_analysis", args.engine, args.ligand_resname, args.recursive_complexes, "*.pdb", args.contact_cutoff, args.ser_resid, args.his_resid, args.acid_resid)
        final_screening_report(outdir / "05_integrated_screen", diffdock_tables, vina_tables, post_tables, args.engine)
        report = outdir / "RUN_SUMMARY.md"
        report.write_text(f"""# Multi-stage ensemble docking run summary

Engine: `{args.engine}`

Generated modules:

1. `01_diffdock_pose_selection/` — DiffDock-L/SMINA clustering and selected poses, if a CSV/SDF directory was provided.
2. `02_scripts/` — ligand preparation, receptor-grid preparation, and AutoDock Vina refinement scripts.
3. `03_vina_refinement/` — parsed Vina log summaries, if logs were provided.
4. `04_post_refinement_analysis/` — contact analysis and CEH catalytic geometry, if complex PDBs were provided.
5. `05_integrated_screen/` — final integrated substrate-specificity ranking.

External commands are not executed automatically. Run the generated shell scripts after checking paths and grid parameters.
""")
        print(f"Workflow outputs written to {outdir}")


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    command_main(args)


if __name__ == "__main__":
    main()
