#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
protein_model_validation_selection_pipeline.py

Integrated protein structural validation and representative-model selection pipeline.

The pipeline can work in two modes:
  1) table mode: score and select models from an existing CSV/XLSX metrics table.
  2) structure-assisted mode: optionally clean/convert structures and run external validation or
     similarity tools when they are available on PATH.

Selection engines:
  - general: fold-aware docking/readiness-style scoring for general protein models.
  - hydrolase: CEH/alpha-beta-hydrolase-specific scoring emphasizing catalytic pocket RMSD,
    catalytic-distance preservation, and motif-region confidence.

External tools are optional and are only called when the corresponding flags are supplied:
  MolProbity/phenix.molprobity, CaBLAM/phenix.cablam, voronota-js-voromqa, TMalign, Foldseek.
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
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except Exception:  # pragma: no cover
    plt = None

try:
    from Bio.PDB import MMCIFParser, PDBIO, PDBParser
except Exception:  # pragma: no cover
    MMCIFParser = PDBIO = PDBParser = None

# ---------------------------------------------------------------------------
# Constants and column normalization
# ---------------------------------------------------------------------------

STRUCTURE_EXTS = (".pdb", ".cif", ".mmcif")

COLUMN_ALIASES = {
    # names used across uploaded scripts/tables
    "model": ["model", "Model", "File", "file", "Protein_ID", "Short_name", "name", "query"],
    "short_name": ["Short_name", "short_name", "File", "model", "Protein_ID"],
    "family": ["family", "Family", "template_family"],
    "molprobity_score": ["molprobity_score", "MolProbity", "molprobity"],
    "clashscore": ["clashscore", "clash_score"],
    "rama_favored": ["rama_favored", "rama_favored_%", "ramachandran_favored", "rama_favored_percent"],
    "rama_outliers": ["rama_outliers", "rama_outliers_%", "ramachandran_outliers", "rama_outlier_percent"],
    "rotamer_outliers": ["rotamer_outliers", "rotamer_outliers_%"],
    "cbeta_deviations": ["cbeta_deviations", "cbeta_outliers", "cbeta_deviations_count"],
    "rms_bonds": ["rms_bonds", "bond_rmsd"],
    "rms_angles": ["rms_angles", "angle_rmsd"],
    "cablam_disfavored": ["cablam_disfavored", "cablam_disfavored_%"],
    "cablam_outliers": ["cablam_outliers", "cablam_outliers_%"],
    "cablam_severe": ["cablam_severe", "cablam_severe_%"],
    "cablam_helix": ["cablam_helix", "cablam_helix_%"],
    "voro_dark": ["voro_dark", "voromqa_dark_score", "voro_dark_score"],
    "voro_light": ["voro_light", "voromqa_light_score", "voro_light_score"],
    "voro_energy_norm": ["voro_energy_norm"],
    "voro_clash": ["voro_clash"],
    "voro_mean": ["voro_mean", "perres_mean"],
    "voro_std": ["voro_std"],
    "perres_q10_mean": ["perres_q10_mean"],
    "pocket_score": ["pocket_score", "pocket_mean"],
    "pocket_score_z": ["pocket_score_z"],
    "tm_score": ["TM_score", "tm_score", "alntmscore"],
    "q_score": ["Q_score", "q_score"],
    "cmo": ["CMO", "cmo"],
    "rmsd": ["RMSD", "rmsd", "best_rmsd_global"],
    "global_rmsd_A": ["global_rmsd_A", "best_rmsd_global"],
    "pocket_rmsd_A": ["pocket_rmsd_A", "best_rmsd_pocket"],
    "catalytic_distance_delta_A": ["catalytic_distance_delta_A", "catalytic_delta_A"],
    "model_catalytic_distance_A": ["model_catalytic_distance_A", "catalytic_distance_A"],
    "ref_catalytic_distance_A": ["ref_catalytic_distance_A", "reference_catalytic_distance_A"],
    "global_cluster": ["global_cluster"],
    "pocket_cluster": ["pocket_cluster"],
    "ptm": ["ptm", "pTM"],
    "iptm": ["iptm", "ipTM"],
    "ranking_score": ["ranking_score"],
    "plddt_mean": ["atom_plddts_mean", "mean_plddt", "plddt_mean"],
    "HGG_plddt_mean": ["HGG_plddt_mean", "HGG_like_plddt_mean"],
    "GESAG_plddt_mean": ["GESAG_plddt_mean", "GxSxG_plddt_mean", "GXSXG_plddt_mean"],
    "HGG_plddt_count": ["HGG_plddt_count", "HGG_like_plddt_count"],
    "GESAG_plddt_count": ["GESAG_plddt_count", "GxSxG_plddt_count", "GXSXG_plddt_count"],
}

METRIC_DIRECTIONS = {
    "tm_score": "higher",
    "q_score": "higher",
    "cmo": "higher",
    "rmsd": "lower",
    "global_rmsd_A": "lower",
    "pocket_rmsd_A": "lower",
    "molprobity_score": "lower",
    "clashscore": "lower",
    "rama_favored": "higher",
    "rama_outliers": "lower",
    "rotamer_outliers": "lower",
    "cbeta_deviations": "lower",
    "rms_bonds": "lower",
    "rms_angles": "lower",
    "cablam_disfavored": "lower",
    "cablam_outliers": "lower",
    "cablam_severe": "lower",
    "voro_dark": "higher",
    "voro_light": "infer",
    "voro_energy_norm": "lower",
    "voro_clash": "lower",
    "voro_mean": "higher",
    "voro_std": "lower",
    "perres_q10_mean": "infer",
    "pocket_score": "higher",
    "pocket_score_z": "higher",
    "ptm": "higher",
    "iptm": "higher",
    "ranking_score": "higher",
    "plddt_mean": "higher",
    "HGG_plddt_mean": "higher",
    "GESAG_plddt_mean": "higher",
}

FOLD_HINTS = {
    "6w5v": "NPC2",
    "npc2": "NPC2",
    "3gkj": "NPC1-NTD",
    "npc1": "NPC1-NTD",
    "ntd": "NPC1-NTD",
    "4tw0": "CD36RP",
    "cd36": "CD36RP",
    "sr-bi": "CD36RP",
    "scarb": "CD36RP",
    "3zwq": "alpha/beta hydrolase",
    "4m0e": "alpha/beta hydrolase",
    "4n5h": "alpha/beta hydrolase",
    "ceh": "alpha/beta hydrolase",
    "lipe": "alpha/beta hydrolase",
    "esterase": "alpha/beta hydrolase",
    "hydrolase": "alpha/beta hydrolase",
}

# ---------------------------------------------------------------------------
# Generic helpers
# ---------------------------------------------------------------------------

def eprint(*args: Any) -> None:
    print(*args, file=sys.stderr)


def tool_exists(name: str) -> bool:
    return shutil.which(name) is not None


def find_first_tool(names: Sequence[str]) -> Optional[str]:
    for n in names:
        if shutil.which(n):
            return n
    return None


def run_cmd(cmd: Sequence[str], timeout: Optional[int] = None) -> str:
    proc = subprocess.run(cmd, text=True, capture_output=True, timeout=timeout)
    if proc.returncode != 0:
        raise RuntimeError(
            "Command failed: " + " ".join(map(str, cmd)) + "\n" +
            "STDERR:\n" + proc.stderr[:4000]
        )
    return proc.stdout


def safe_float(x: Any) -> Optional[float]:
    try:
        if pd.isna(x):
            return None
        y = float(x)
        if math.isfinite(y):
            return y
    except Exception:
        pass
    return None


def list_structures(directory: Optional[str]) -> List[Path]:
    if not directory:
        return []
    d = Path(directory)
    if not d.exists():
        return []
    return sorted(p for p in d.iterdir() if p.suffix.lower() in STRUCTURE_EXTS)


def safe_stem(path_or_name: Any) -> str:
    stem = Path(str(path_or_name)).stem
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", stem).strip("_") or "model"


def read_table(path: str, sheet: Optional[str] = None) -> pd.DataFrame:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(path)
    if p.suffix.lower() in (".xlsx", ".xls"):
        return pd.read_excel(p, sheet_name=sheet or 0)
    if p.suffix.lower() in (".tsv", ".tab"):
        return pd.read_csv(p, sep="\t")
    return pd.read_csv(p)


def canonicalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    used = set(out.columns)
    for canonical, aliases in COLUMN_ALIASES.items():
        if canonical in out.columns:
            continue
        for a in aliases:
            if a in out.columns:
                out[canonical] = out[a]
                break
    if "model" not in out.columns:
        out["model"] = [f"model_{i+1}" for i in range(len(out))]
    if "short_name" not in out.columns:
        out["short_name"] = out["model"].map(safe_stem)
    return out


def tag_family_from_name(name: Any, explicit: Optional[str] = None) -> str:
    if explicit and str(explicit).strip() and str(explicit).lower() not in ("nan", "none", "other"):
        val = str(explicit).strip()
        if any(k in val.lower() for k in ["full", "truncated", "configuration", "domain", "refined"]):
            pass
        else:
            return val
    s = str(name).lower()
    for key, fam in FOLD_HINTS.items():
        if key in s:
            return fam
    return "Other"


def robust_percentile_scale(series: pd.Series, direction: str = "higher") -> pd.Series:
    """Map a metric to [0,1] using 5th/95th percentiles. Missing values become NaN."""
    x = pd.to_numeric(series, errors="coerce")
    finite = x[np.isfinite(x)]
    if finite.empty:
        return pd.Series(np.nan, index=series.index, dtype=float)
    p5, p95 = np.nanpercentile(finite, [5, 95])
    if not np.isfinite(p5) or not np.isfinite(p95) or p95 == p5:
        return pd.Series(0.5, index=series.index, dtype=float).where(x.notna(), np.nan)
    clipped = x.clip(p5, p95)
    if direction == "lower":
        score = (p95 - clipped) / (p95 - p5)
    else:
        score = (clipped - p5) / (p95 - p5)
    return score.clip(0, 1)


def scale_abs_closer_zero(series: pd.Series) -> pd.Series:
    x = pd.to_numeric(series, errors="coerce").abs()
    finite = x[np.isfinite(x)]
    if finite.empty:
        return pd.Series(np.nan, index=series.index, dtype=float)
    p95 = np.nanpercentile(finite, 95)
    if not np.isfinite(p95) or p95 == 0:
        return pd.Series(1.0, index=series.index, dtype=float).where(x.notna(), np.nan)
    return (1 - (x / p95)).clip(0, 1)


def infer_voro_direction(df: pd.DataFrame, voro_col: str = "voro_light", ref_col: str = "molprobity_score") -> str:
    if voro_col not in df.columns or ref_col not in df.columns:
        return "higher"
    v = pd.to_numeric(df[voro_col], errors="coerce")
    r = pd.to_numeric(df[ref_col], errors="coerce")
    m = v.notna() & r.notna()
    if m.sum() < 5:
        return "higher"
    corr = np.corrcoef(v[m], r[m])[0, 1]
    # If higher VoroMQA light score tracks worse MolProbity, treat lower as better.
    return "lower" if np.isfinite(corr) and corr > 0 else "higher"


def weighted_component(scores: Dict[str, pd.Series], weights: Dict[str, float], index: pd.Index) -> pd.Series:
    available = [k for k in weights if k in scores and scores[k].notna().any()]
    if not available:
        return pd.Series(np.nan, index=index, dtype=float)
    total = float(sum(weights[k] for k in available))
    out = pd.Series(0.0, index=index, dtype=float)
    w_used = pd.Series(0.0, index=index, dtype=float)
    for k in available:
        s = scores[k]
        w = weights[k] / total
        out = out.add(s.fillna(0) * w, fill_value=0)
        w_used = w_used.add(s.notna().astype(float) * w, fill_value=0)
    out = out.where(w_used > 0, np.nan)
    return out / w_used.where(w_used > 0, np.nan)

# ---------------------------------------------------------------------------
# Structure conversion and optional validation wrappers
# ---------------------------------------------------------------------------

def convert_cif_to_pdb(cif_path: Path, out_dir: Path) -> Optional[Path]:
    if MMCIFParser is None or PDBIO is None:
        eprint("[WARN] Biopython Bio.PDB is unavailable; cannot convert CIF/mmCIF to PDB.")
        return None
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{safe_stem(cif_path)}.pdb"
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure(safe_stem(cif_path), str(cif_path))
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(out_path))
    return out_path


def clean_pdb_file(inp: Path, out: Path, keep_hetatm: bool = False, remove_hydrogens: bool = True,
                   first_model_only: bool = True, altloc_policy: str = "A") -> Path:
    out.parent.mkdir(parents=True, exist_ok=True)
    altloc_policy = altloc_policy.lower()

    def keep_altloc(line: str) -> bool:
        alt = line[16:17].strip()
        if altloc_policy == "all":
            return True
        if altloc_policy == "blank":
            return alt == ""
        return alt in ("", "A")

    with inp.open("r", errors="ignore") as fin, out.open("w") as fout:
        for line in fin:
            rec = line[:6].strip()
            if first_model_only and rec == "ENDMDL":
                break
            if rec in ("ATOM", "HETATM"):
                if rec == "HETATM" and not keep_hetatm:
                    continue
                if remove_hydrogens and line[76:78].strip() == "H":
                    continue
                if not keep_altloc(line):
                    continue
                fout.write(line)
            elif rec == "TER":
                fout.write(line)
            elif rec in ("MODEL", "END"):
                fout.write(line)
    return out


def parse_molprobity(text: str) -> Dict[str, Optional[float]]:
    def extract(pattern: str, cast: Any = float) -> Optional[Any]:
        m = re.search(pattern, text, flags=re.IGNORECASE)
        if not m:
            return None
        try:
            return cast(m.group(1))
        except Exception:
            return None
    return {
        "molprobity_score": extract(r"MolProbity score\s*=\s*([\d.]+)"),
        "clashscore": extract(r"Clashscore\s*=\s*([\d.]+)"),
        "rama_favored": extract(r"Favored\s*[:=]\s*([\d.]+)\s*%"),
        "rama_outliers": extract(r"Outliers\s*[:=]\s*([\d.]+)\s*%"),
        "rotamer_outliers": extract(r"Rotamer.*?Outliers\s*[:=]\s*([\d.]+)\s*%"),
        "cbeta_deviations": extract(r"C-beta deviations\s*=\s*(\d+)", int),
        "rms_bonds": extract(r"RMS\(bonds\)\s*=\s*([\d.]+)"),
        "rms_angles": extract(r"RMS\(angles\)\s*=\s*([\d.]+)"),
    }


def parse_cablam(text: str) -> Dict[str, Optional[float]]:
    data = {"cablam_disfavored": None, "cablam_outliers": None, "cablam_severe": None, "cablam_helix": None}
    for line in text.splitlines():
        m = re.search(r"\(([\d.]+)%\)", line)
        if not m:
            continue
        value = float(m.group(1))
        low = line.lower()
        if "disfavored conformations" in low:
            data["cablam_disfavored"] = value
        elif "outlier conformations" in low and "severe" not in low:
            data["cablam_outliers"] = value
        elif "severe ca geometry outliers" in low:
            data["cablam_severe"] = value
        elif "helix-like" in low:
            data["cablam_helix"] = value
    return data


def run_voromqa_global(pdb: Path) -> Dict[str, Optional[float]]:
    cmd = find_first_tool(["voronota-js-voromqa"])
    if not cmd:
        raise RuntimeError("voronota-js-voromqa not found in PATH")
    out = run_cmd([cmd, "--input", str(pdb)])
    nums = [float(x) for x in out.split() if re.fullmatch(r"-?[0-9.]+", x)]
    return {
        "voro_dark": nums[0] if len(nums) > 0 else None,
        "voro_light": nums[1] if len(nums) > 1 else None,
    }


def validate_structures_external(structures: List[Path], outdir: Path) -> pd.DataFrame:
    """Run optional external validation tools. Missing tools are recorded, not fatal."""
    rows = []
    clean_dir = outdir / "cleaned_pdb"
    pdb_dir = outdir / "pdb_converted"
    for path in structures:
        row: Dict[str, Any] = {"model": path.name, "source_path": str(path)}
        try:
            if path.suffix.lower() in (".cif", ".mmcif"):
                converted = convert_cif_to_pdb(path, pdb_dir)
                if converted is None:
                    row["status"] = "conversion_failed"
                    rows.append(row)
                    continue
                pdb = converted
            else:
                pdb = path
            clean = clean_pdb_file(pdb, clean_dir / pdb.name)
            row["clean_pdb"] = str(clean)
            mp_cmd = find_first_tool(["molprobity", "phenix.molprobity"])
            if mp_cmd:
                try:
                    row.update(parse_molprobity(run_cmd([mp_cmd, str(clean)])))
                except Exception as e:
                    row["molprobity_error"] = str(e)
            else:
                row["molprobity_error"] = "not_found"
            cb_cmd = find_first_tool(["cablam", "phenix.cablam"])
            if cb_cmd:
                try:
                    row.update(parse_cablam(run_cmd([cb_cmd, str(clean)])))
                except Exception as e:
                    row["cablam_error"] = str(e)
            else:
                row["cablam_error"] = "not_found"
            if tool_exists("voronota-js-voromqa"):
                try:
                    row.update(run_voromqa_global(clean))
                except Exception as e:
                    row["voromqa_error"] = str(e)
            else:
                row["voromqa_error"] = "not_found"
            row["status"] = "ok"
        except Exception as e:
            row["status"] = "failed"
            row["error"] = str(e)
        rows.append(row)
    return pd.DataFrame(rows)

# ---------------------------------------------------------------------------
# Optional Foldseek / TMalign wrappers
# ---------------------------------------------------------------------------

def extract_best_foldseek_hits(results_tsv: Path, out_tsv: Path) -> pd.DataFrame:
    cols = ["query", "target", "evalue", "bits", "alntmscore", "rmsd", "pident"]
    df = pd.read_csv(results_tsv, sep="\t", header=None, names=cols)
    df["alntmscore"] = pd.to_numeric(df["alntmscore"], errors="coerce")
    best = df.sort_values(["query", "alntmscore"], ascending=[True, False]).groupby("query", as_index=False).first()
    best.to_csv(out_tsv, sep="\t", index=False)
    return best


def run_foldseek_search(query_dir: str, target_db: str, outdir: Path, threads: int = 8,
                        evalue: str = "1e-3") -> Optional[pd.DataFrame]:
    if not tool_exists("foldseek"):
        eprint("[WARN] Foldseek not found; skipping structural template search.")
        return None
    outdir.mkdir(parents=True, exist_ok=True)
    tmp = outdir / "tmp"
    tmp.mkdir(exist_ok=True)
    query_db = outdir / "query_db"
    result_db = outdir / "result_db"
    results_tsv = outdir / "foldseek_results.tsv"
    best_tsv = outdir / "foldseek_best_hits.tsv"
    run_cmd(["foldseek", "createdb", query_dir, str(query_db)])
    run_cmd(["foldseek", "search", str(query_db), target_db, str(result_db), str(tmp),
             "--threads", str(threads), "-e", str(evalue), "--alignment-type", "2", "-a"])
    run_cmd(["foldseek", "convertalis", str(query_db), target_db, str(result_db), str(results_tsv),
             "--format-output", "query,target,evalue,bits,alntmscore,rmsd,pident"])
    return extract_best_foldseek_hits(results_tsv, best_tsv)


def run_tmalign_pair(model: Path, template: Path) -> Dict[str, Optional[float]]:
    if not tool_exists("TMalign"):
        raise RuntimeError("TMalign not found in PATH")
    out = run_cmd(["TMalign", str(model), str(template)])
    rmsd = None
    alnlen = None
    len1 = None
    len2 = None
    tm1 = None
    m = re.search(r"RMSD=\s*([\d.]+)", out)
    if m:
        rmsd = float(m.group(1))
    m = re.search(r"Aligned length=\s*(\d+)", out)
    if m:
        alnlen = int(m.group(1))
    m = re.search(r"Length of Chain_1:\s*(\d+)", out)
    if m:
        len1 = int(m.group(1))
    m = re.search(r"Length of Chain_2:\s*(\d+)", out)
    if m:
        len2 = int(m.group(1))
    tms = re.findall(r"TM-score=\s*([\d.]+)", out)
    if tms:
        tm1 = float(tms[0])
    q = None
    if rmsd is not None and alnlen and len1 and len2:
        q = (alnlen / math.sqrt(len1 * len2)) * (1.0 / (1.0 + (rmsd ** 2 / 3.0 ** 2)))
    return {"rmsd": rmsd, "tm_score": tm1, "q_score": q, "aligned_length": alnlen}


def compare_models_to_templates(models_dir: str, templates_dir: str, outdir: Path) -> Optional[pd.DataFrame]:
    if not tool_exists("TMalign"):
        eprint("[WARN] TMalign not found; skipping model-template comparison.")
        return None
    models = list_structures(models_dir)
    templates = list_structures(templates_dir)
    rows = []
    for m in models:
        for t in templates:
            try:
                r = run_tmalign_pair(m, t)
                rows.append({"model": m.name, "template": t.name, **r})
            except Exception as e:
                rows.append({"model": m.name, "template": t.name, "tmalign_error": str(e)})
    df = pd.DataFrame(rows)
    if not df.empty:
        df.to_csv(outdir / "model_template_tmalign.tsv", sep="\t", index=False)
        best = df.copy()
        best["tm_score"] = pd.to_numeric(best.get("tm_score"), errors="coerce")
        best = best.sort_values(["model", "tm_score", "rmsd"], ascending=[True, False, True]).groupby("model", as_index=False).first()
        best.to_csv(outdir / "best_template_per_model.tsv", sep="\t", index=False)
    return df

# ---------------------------------------------------------------------------
# Selection engines
# ---------------------------------------------------------------------------

@dataclass
class SelectionConfig:
    engine: str = "auto"
    family: Optional[str] = None
    top_n: int = 10
    global_cluster_col: str = "global_cluster"
    pocket_cluster_col: str = "pocket_cluster"
    high_threshold: float = 75.0
    pass_threshold: float = 60.0
    borderline_threshold: float = 45.0


def make_scaled_metrics(df: pd.DataFrame) -> Tuple[Dict[str, pd.Series], Dict[str, str]]:
    directions = dict(METRIC_DIRECTIONS)
    voro_dir = infer_voro_direction(df, "voro_light", "molprobity_score")
    directions["voro_light"] = voro_dir
    directions["perres_q10_mean"] = voro_dir
    scaled: Dict[str, pd.Series] = {}
    for col, direction in directions.items():
        if col not in df.columns:
            continue
        if col == "catalytic_distance_delta_A":
            scaled[col] = scale_abs_closer_zero(df[col])
        else:
            scaled[col] = robust_percentile_scale(df[col], "lower" if direction == "lower" else "higher")
    return scaled, {"voro_direction_inferred": voro_dir}


def validation_gate_general(df: pd.DataFrame) -> pd.Series:
    gate = pd.Series(1.0, index=df.index, dtype=float)
    rules = [
        ("molprobity_score", 3.0, "le"),
        ("clashscore", 20.0, "le"),
        ("rama_outliers", 2.0, "le"),
        ("cablam_severe", 2.0, "le"),
    ]
    for col, thr, direction in rules:
        if col in df.columns:
            x = pd.to_numeric(df[col], errors="coerce")
            ok = (x <= thr) if direction == "le" else (x >= thr)
            gate *= (ok | x.isna()).astype(float)
    return gate


def validation_gate_hydrolase(df: pd.DataFrame) -> pd.Series:
    gate = pd.Series(1.0, index=df.index, dtype=float)
    rules = [
        ("molprobity_score", 2.5, "le"),
        ("clashscore", 15.0, "le"),
        ("rama_outliers", 1.0, "le"),
        ("cablam_severe", 1.0, "le"),
    ]
    for col, thr, direction in rules:
        if col in df.columns:
            x = pd.to_numeric(df[col], errors="coerce")
            ok = (x <= thr) if direction == "le" else (x >= thr)
            gate *= (ok | x.isna()).astype(float)
    return gate


def compute_general_scores(df: pd.DataFrame, cfg: SelectionConfig) -> pd.DataFrame:
    out = df.copy()
    scaled, meta = make_scaled_metrics(out)
    out["voro_direction_inferred"] = meta["voro_direction_inferred"]

    geometry_w = {
        "clashscore": 0.20,
        "rama_favored": 0.20,
        "rotamer_outliers": 0.10,
        "molprobity_score": 0.20,
        "cbeta_deviations": 0.10,
        "rama_outliers": 0.10,
        "rms_bonds": 0.05,
        "rms_angles": 0.05,
    }
    backbone_w = {"cablam_disfavored": 0.20, "cablam_outliers": 0.30, "cablam_severe": 0.50}
    topology_w = {"tm_score": 0.55, "q_score": 0.25, "cmo": 0.10, "rmsd": 0.10}
    packing_w = {
        "voro_dark": 0.20,
        "voro_light": 0.20,
        "voro_mean": 0.20,
        "voro_energy_norm": 0.10,
        "voro_clash": 0.10,
        "voro_std": 0.10,
        "pocket_score": 0.05,
        "pocket_score_z": 0.05,
    }
    confidence_w = {"plddt_mean": 0.40, "ptm": 0.25, "iptm": 0.25, "ranking_score": 0.10}

    out["S_geometry"] = weighted_component(scaled, geometry_w, out.index)
    out["S_backbone"] = weighted_component(scaled, backbone_w, out.index)
    out["S_topology"] = weighted_component(scaled, topology_w, out.index)
    out["S_packing"] = weighted_component(scaled, packing_w, out.index)
    out["S_confidence"] = weighted_component(scaled, confidence_w, out.index)

    comp = {
        "S_geometry": 0.28,
        "S_backbone": 0.22,
        "S_topology": 0.20,
        "S_packing": 0.20,
        "S_confidence": 0.10,
    }
    comp_scores = {k: out[k] for k in comp if k in out.columns}
    out["validation_gate"] = validation_gate_general(out)
    raw = weighted_component(comp_scores, comp, out.index)
    out["SelectionScore_raw"] = raw
    out["SelectionScore"] = 100.0 * out["validation_gate"] * raw.fillna(0)
    out["SelectionEngine"] = "general"
    return assign_decisions(out, cfg)


def compute_hydrolase_scores(df: pd.DataFrame, cfg: SelectionConfig) -> pd.DataFrame:
    out = df.copy()
    scaled, meta = make_scaled_metrics(out)
    out["voro_direction_inferred"] = meta["voro_direction_inferred"]
    # Catalytic distance is special: absolute delta closest to zero is best.
    if "catalytic_distance_delta_A" in out.columns:
        scaled["catalytic_distance_delta_A"] = scale_abs_closer_zero(out["catalytic_distance_delta_A"])

    validation_w = {
        "molprobity_score": 0.20,
        "clashscore": 0.10,
        "rama_outliers": 0.10,
        "rotamer_outliers": 0.05,
        "cbeta_deviations": 0.03,
        "cablam_severe": 0.20,
        "cablam_outliers": 0.10,
        "cablam_disfavored": 0.05,
        "voro_light": 0.12,
        "perres_q10_mean": 0.05,
    }
    ensemble_w = {
        "pocket_rmsd_A": 0.45,
        "catalytic_distance_delta_A": 0.35,
        "global_rmsd_A": 0.20,
    }
    functional_w = {
        "HGG_plddt_mean": 0.30,
        "GESAG_plddt_mean": 0.30,
        "plddt_mean": 0.20,
        "ptm": 0.10,
        "iptm": 0.10,
    }

    out["S_validation"] = weighted_component(scaled, validation_w, out.index)
    out["S_ensemble"] = weighted_component(scaled, ensemble_w, out.index)
    out["S_functional_site"] = weighted_component(scaled, functional_w, out.index)

    # Motif completeness modifier. Missing motif hits should not silently receive full confidence.
    motif_counts = []
    for col in ["HGG_plddt_count", "GESAG_plddt_count"]:
        if col in out.columns:
            motif_counts.append(pd.to_numeric(out[col], errors="coerce").fillna(0) > 0)
    if motif_counts:
        motif_ok = motif_counts[0].astype(float)
        for m in motif_counts[1:]:
            motif_ok += m.astype(float)
        motif_modifier = (0.75 + 0.25 * (motif_ok / len(motif_counts))).clip(0.75, 1.0)
    else:
        motif_modifier = pd.Series(0.90, index=out.index, dtype=float)
    out["motif_completeness_modifier"] = motif_modifier

    comp = {"S_validation": 0.35, "S_ensemble": 0.35, "S_functional_site": 0.30}
    comp_scores = {k: out[k] for k in comp if k in out.columns}
    out["validation_gate"] = validation_gate_hydrolase(out)
    raw = weighted_component(comp_scores, comp, out.index)
    out["SelectionScore_raw"] = raw
    out["SelectionScore"] = 100.0 * out["validation_gate"] * out["motif_completeness_modifier"] * raw.fillna(0)
    out["SelectionEngine"] = "hydrolase"
    return assign_decisions(out, cfg)


def assign_decisions(df: pd.DataFrame, cfg: SelectionConfig) -> pd.DataFrame:
    out = df.copy()
    score = pd.to_numeric(out["SelectionScore"], errors="coerce").fillna(0)
    decisions = []
    for x in score:
        if x >= cfg.high_threshold:
            decisions.append("High-confidence representative")
        elif x >= cfg.pass_threshold:
            decisions.append("Representative candidate")
        elif x >= cfg.borderline_threshold:
            decisions.append("Borderline; inspect/refine")
        else:
            decisions.append("Reject or remodel")
    out["SelectionDecision"] = decisions
    out["Rank"] = score.rank(ascending=False, method="dense").astype(int)
    return out.sort_values(["SelectionScore", "short_name"], ascending=[False, True])


def resolve_engine(df: pd.DataFrame, requested: str) -> str:
    requested = requested.lower()
    if requested in ("general", "hydrolase"):
        return requested
    # auto-detect hydrolase/CEH when catalytic geometry or motif-level columns are present.
    hydrolase_cols = {"pocket_rmsd_A", "catalytic_distance_delta_A", "HGG_plddt_mean", "GESAG_plddt_mean"}
    if any(c in df.columns for c in hydrolase_cols):
        return "hydrolase"
    return "general"


def select_representatives(scored: pd.DataFrame, cfg: SelectionConfig) -> pd.DataFrame:
    out_rows = []
    engine = str(scored["SelectionEngine"].iloc[0]) if not scored.empty and "SelectionEngine" in scored else cfg.engine
    # Always keep overall top models.
    top = scored.head(cfg.top_n).copy()
    top["RepresentativeClass"] = "overall_top"
    out_rows.append(top)

    # Add cluster representatives. Hydrolases should use pocket cluster first.
    cluster_col = cfg.pocket_cluster_col if engine == "hydrolase" else cfg.global_cluster_col
    if cluster_col in scored.columns:
        cluster_reps = (scored.sort_values("SelectionScore", ascending=False)
                        .groupby(cluster_col, dropna=False, as_index=False)
                        .head(1)
                        .copy())
        cluster_reps["RepresentativeClass"] = f"best_per_{cluster_col}"
        out_rows.append(cluster_reps)

    # Add family representatives when a family column is available.
    if "family_resolved" in scored.columns:
        fam_reps = (scored.sort_values("SelectionScore", ascending=False)
                    .groupby("family_resolved", dropna=False, as_index=False)
                    .head(1)
                    .copy())
        fam_reps["RepresentativeClass"] = "best_per_family"
        out_rows.append(fam_reps)

    reps = pd.concat(out_rows, ignore_index=True) if out_rows else pd.DataFrame()
    if reps.empty:
        return reps
    # Preserve multiple classes per same model.
    key = "model" if "model" in reps.columns else "short_name"
    reps["RepresentativeClass"] = reps.groupby(key)["RepresentativeClass"].transform(lambda x: ";".join(sorted(set(map(str, x)))))
    reps = reps.drop_duplicates(subset=[key]).sort_values("SelectionScore", ascending=False)
    return reps

# ---------------------------------------------------------------------------
# Plotting and reports
# ---------------------------------------------------------------------------

def save_score_plots(scored: pd.DataFrame, outdir: Path, top_n: int = 20) -> None:
    if plt is None or scored.empty:
        return
    plot_dir = outdir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    view = scored.sort_values("SelectionScore", ascending=False).head(top_n)
    labels = view.get("short_name", view.get("model", pd.Series(range(len(view))))).astype(str).tolist()
    scores = pd.to_numeric(view["SelectionScore"], errors="coerce").fillna(0).to_numpy()

    fig, ax = plt.subplots(figsize=(max(8, 0.45 * len(labels)), 5), dpi=200)
    ax.bar(range(len(labels)), scores)
    ax.set_ylabel("Selection score")
    ax.set_ylim(0, 100)
    ax.set_title("Top structural model candidates")
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=75, ha="right", fontsize=8)
    fig.tight_layout()
    fig.savefig(plot_dir / "top_selection_scores.png")
    plt.close(fig)

    # Scatter for hydrolase geometry when available.
    if {"pocket_rmsd_A", "catalytic_distance_delta_A", "SelectionScore"}.issubset(scored.columns):
        x = pd.to_numeric(scored["pocket_rmsd_A"], errors="coerce")
        y = pd.to_numeric(scored["catalytic_distance_delta_A"], errors="coerce").abs()
        c = pd.to_numeric(scored["SelectionScore"], errors="coerce")
        m = x.notna() & y.notna() & c.notna()
        if m.sum() >= 2:
            fig, ax = plt.subplots(figsize=(6, 5), dpi=200)
            sc = ax.scatter(x[m], y[m], c=c[m], s=30)
            ax.set_xlabel("Pocket RMSD (Å)")
            ax.set_ylabel("|Catalytic-distance delta| (Å)")
            ax.set_title("Hydrolase pocket geometry vs selection score")
            cb = fig.colorbar(sc, ax=ax)
            cb.set_label("Selection score")
            fig.tight_layout()
            fig.savefig(plot_dir / "hydrolase_pocket_geometry_selection.png")
            plt.close(fig)


def write_markdown_report(scored: pd.DataFrame, reps: pd.DataFrame, outdir: Path, cfg: SelectionConfig) -> Path:
    out = outdir / "validation_selection_report.md"
    engine = scored["SelectionEngine"].iloc[0] if not scored.empty and "SelectionEngine" in scored else cfg.engine
    lines = [
        "# Protein Structural Validation and Representative Model Selection Report",
        "",
        f"Selection engine: **{engine}**",
        f"Models scored: **{len(scored)}**",
        f"Representatives selected: **{len(reps)}**",
        "",
        "## Decision summary",
        "",
    ]
    if "SelectionDecision" in scored.columns:
        counts = scored["SelectionDecision"].value_counts().to_dict()
        for k, v in counts.items():
            lines.append(f"- {k}: {v}")
    lines += ["", "## Top representative candidates", ""]
    cols = [c for c in ["Rank", "short_name", "model", "family_resolved", "SelectionScore", "SelectionDecision", "RepresentativeClass"] if c in reps.columns]
    if cols and not reps.empty:
        lines.append("| " + " | ".join(cols) + " |")
        lines.append("| " + " | ".join(["---"] * len(cols)) + " |")
        for _, r in reps.head(20).iterrows():
            vals = []
            for c in cols:
                v = r.get(c, "")
                if isinstance(v, (float, np.floating)):
                    vals.append(f"{float(v):.3f}")
                else:
                    vals.append(str(v))
            lines.append("| " + " | ".join(vals) + " |")
    else:
        lines.append("No representatives selected.")
    lines += [
        "",
        "## Interpretation notes",
        "",
        "- `SelectionScore` is scaled from 0 to 100 after metric normalization and validation gating.",
        "- The general engine balances stereochemistry, backbone quality, fold/template similarity, packing, and model-confidence fields.",
        "- The hydrolase engine uses stricter validation gates and adds catalytic pocket RMSD, catalytic-distance preservation, and motif-level pLDDT.",
        "- Missing metrics are ignored and available component weights are renormalized; therefore, compare models most confidently when they share the same input metric set.",
        "",
    ]
    out.write_text("\n".join(lines), encoding="utf-8")
    return out


def write_outputs(scored: pd.DataFrame, reps: pd.DataFrame, outdir: Path, cfg: SelectionConfig) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    scored.to_csv(outdir / "all_models_scored.tsv", sep="\t", index=False)
    reps.to_csv(outdir / "representative_models.tsv", sep="\t", index=False)
    try:
        with pd.ExcelWriter(outdir / "structural_validation_selection.xlsx", engine="openpyxl") as writer:
            scored.to_excel(writer, sheet_name="ALL_MODELS_SCORED", index=False)
            reps.to_excel(writer, sheet_name="REPRESENTATIVES", index=False)
            score_cols = [c for c in scored.columns if c.startswith("S_") or c in ["SelectionScore", "validation_gate", "SelectionDecision", "SelectionEngine"]]
            if score_cols:
                scored[["short_name"] + [c for c in score_cols if c != "short_name"]].to_excel(writer, sheet_name="SCORE_COMPONENTS", index=False)
    except Exception as e:
        eprint(f"[WARN] Excel workbook could not be written: {e}")
    save_score_plots(scored, outdir, top_n=max(cfg.top_n, 20))
    write_markdown_report(scored, reps, outdir, cfg)

# ---------------------------------------------------------------------------
# Main workflow
# ---------------------------------------------------------------------------

def merge_metric_tables(primary: Optional[pd.DataFrame], extras: List[pd.DataFrame]) -> pd.DataFrame:
    tables = []
    if primary is not None and not primary.empty:
        tables.append(primary)
    tables.extend([x for x in extras if x is not None and not x.empty])
    if not tables:
        return pd.DataFrame()
    # Combine by model where possible; otherwise concatenate.
    normalized = [canonicalize_columns(t) for t in tables]
    if len(normalized) == 1:
        return normalized[0]
    merged = normalized[0]
    for t in normalized[1:]:
        common_key = "model" if "model" in merged.columns and "model" in t.columns else None
        if common_key:
            # Avoid overwriting existing columns with _y unless needed.
            merged = merged.merge(t, on=common_key, how="outer", suffixes=("", "_extra"))
            for c in list(merged.columns):
                if c.endswith("_extra"):
                    base = c[:-6]
                    if base not in merged.columns:
                        merged[base] = merged[c]
                    else:
                        merged[base] = merged[base].combine_first(merged[c])
                    merged = merged.drop(columns=[c])
        else:
            merged = pd.concat([merged, t], ignore_index=True, sort=False)
    return merged


def run_pipeline(args: argparse.Namespace) -> Tuple[pd.DataFrame, pd.DataFrame]:
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    input_table = None
    if args.metrics_table:
        input_table = read_table(args.metrics_table, args.sheet)
        input_table["input_table"] = os.path.basename(args.metrics_table)

    extra_tables: List[pd.DataFrame] = []

    # Optional structure discovery and validation.
    structures = list_structures(args.models_dir)
    if structures:
        # Always create converted/cleaned PDBs when requested or when validation uses them.
        if args.convert_cif:
            conv_rows = []
            conv_dir = outdir / "pdb_converted"
            for s in structures:
                if s.suffix.lower() in (".cif", ".mmcif"):
                    try:
                        p = convert_cif_to_pdb(s, conv_dir)
                        conv_rows.append({"model": s.name, "converted_pdb": str(p) if p else ""})
                    except Exception as e:
                        conv_rows.append({"model": s.name, "conversion_error": str(e)})
            if conv_rows:
                extra_tables.append(pd.DataFrame(conv_rows))
        if args.run_external_validation:
            extra_tables.append(validate_structures_external(structures, outdir / "external_validation"))

    if args.run_tmalign and args.models_dir and args.templates_dir:
        tdf = compare_models_to_templates(args.models_dir, args.templates_dir, outdir / "tmalign")
        if tdf is not None and not tdf.empty:
            # Use best template summary as scoring input.
            best = (tdf.sort_values(["model", "tm_score", "rmsd"], ascending=[True, False, True])
                    .groupby("model", as_index=False).first())
            extra_tables.append(best)

    if args.run_foldseek and args.models_dir and args.foldseek_db:
        fdf = run_foldseek_search(args.models_dir, args.foldseek_db, outdir / "foldseek",
                                  threads=args.threads, evalue=args.foldseek_evalue)
        if fdf is not None and not fdf.empty:
            fdf = fdf.rename(columns={"query": "model", "alntmscore": "tm_score"})
            extra_tables.append(fdf)

    df = merge_metric_tables(input_table, extra_tables)
    if df.empty:
        raise SystemExit("No model metrics available. Provide --metrics-table and/or --models-dir with external-tool flags.")

    df = canonicalize_columns(df)
    df["family_resolved"] = [tag_family_from_name(n, f if "family" in df.columns else args.family)
                              for n, f in zip(df.get("short_name", df["model"]), df.get("family", pd.Series([args.family] * len(df))))]
    if args.family:
        df["family_resolved"] = args.family

    cfg = SelectionConfig(
        engine=args.engine,
        family=args.family,
        top_n=args.top_n,
        high_threshold=args.high_threshold,
        pass_threshold=args.pass_threshold,
        borderline_threshold=args.borderline_threshold,
    )
    engine = resolve_engine(df, args.engine)
    cfg.engine = engine
    if engine == "hydrolase":
        scored = compute_hydrolase_scores(df, cfg)
    else:
        scored = compute_general_scores(df, cfg)
    reps = select_representatives(scored, cfg)
    write_outputs(scored, reps, outdir, cfg)
    return scored, reps


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Integrated structural validation and representative model-selection pipeline."
    )
    p.add_argument("--metrics-table", default=None, help="Existing CSV/TSV/XLSX metrics table to score.")
    p.add_argument("--sheet", default=None, help="Excel sheet name/index for --metrics-table.")
    p.add_argument("--models-dir", default=None, help="Directory containing model structures (.pdb/.cif/.mmcif).")
    p.add_argument("--templates-dir", default=None, help="Directory containing template/reference structures for TMalign.")
    p.add_argument("--outdir", default="validation_selection_results", help="Output directory.")
    p.add_argument("--engine", choices=["auto", "general", "hydrolase"], default="auto",
                   help="Selection engine. auto switches to hydrolase when catalytic/pocket metrics are present.")
    p.add_argument("--family", default=None, help="Optional fixed family label for all rows.")
    p.add_argument("--top-n", type=int, default=10, help="Number of overall top models to retain as representatives.")
    p.add_argument("--convert-cif", action="store_true", help="Convert CIF/mmCIF models to PDB using Biopython.")
    p.add_argument("--run-external-validation", action="store_true",
                   help="Run MolProbity/CaBLAM/VoroMQA if available on PATH.")
    p.add_argument("--run-tmalign", action="store_true", help="Run TMalign model-template comparison if available.")
    p.add_argument("--run-foldseek", action="store_true", help="Run Foldseek template search if available.")
    p.add_argument("--foldseek-db", default=None, help="Foldseek target database path/name.")
    p.add_argument("--foldseek-evalue", default="1e-3", help="Foldseek E-value threshold.")
    p.add_argument("--threads", type=int, default=8, help="Threads for Foldseek.")
    p.add_argument("--high-threshold", type=float, default=75.0)
    p.add_argument("--pass-threshold", type=float, default=60.0)
    p.add_argument("--borderline-threshold", type=float, default=45.0)
    return p


def main() -> None:
    args = build_arg_parser().parse_args()
    scored, reps = run_pipeline(args)
    print(f"Scored models: {len(scored)}")
    print(f"Representatives: {len(reps)}")
    print(f"Wrote outputs to: {args.outdir}")


if __name__ == "__main__":
    main()
