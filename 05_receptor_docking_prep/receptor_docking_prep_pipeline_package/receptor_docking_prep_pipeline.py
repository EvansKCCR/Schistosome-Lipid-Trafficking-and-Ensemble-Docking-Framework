#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Integrated receptor preparation, MD validation, and representative-structure selection pipeline.

Designed for preparing docking receptors from predicted or experimental protein models.
Includes a CEH-like / alpha-beta hydrolase workflow with catalytic-distance checks.

Core modules
------------
1) refinement-script generation for restrained GROMACS relaxation
2) MD validation command generation for RMSD, RMSF, Rg, SASA, DSSP, and CEH distances
3) XVG parsing and comparative plotting
4) model/frame scoring and representative selection
5) docking-receptor extraction command generation

The pipeline writes shell scripts by default. Use --execute only when GROMACS and
other external tools are available and you want the commands run immediately.
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
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# =============================================================================
# General helpers
# =============================================================================

STRUCTURE_EXTS = (".pdb", ".cif", ".mmcif", ".gro")
PDB_EXTS = (".pdb",)


def safe_mkdir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def write_text(path: Path, text: str, executable: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    if executable:
        path.chmod(path.stat().st_mode | 0o111)


def run_cmd(cmd: Sequence[str], cwd: Optional[Path] = None) -> None:
    print("[RUN]", " ".join(str(x) for x in cmd))
    subprocess.run(cmd, check=True, cwd=str(cwd) if cwd else None)


def list_files(directory: Path, suffixes: Tuple[str, ...]) -> List[Path]:
    if not directory.exists():
        return []
    return sorted([p for p in directory.rglob("*") if p.is_file() and p.suffix.lower() in suffixes])


def stem_clean(path: Path) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", path.stem)


def to_numeric_series(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce")


def first_col(df: pd.DataFrame, candidates: Sequence[str]) -> Optional[str]:
    lower_to_original = {c.lower(): c for c in df.columns}
    for c in candidates:
        if c in df.columns:
            return c
        if c.lower() in lower_to_original:
            return lower_to_original[c.lower()]
    return None


# =============================================================================
# XVG parsing and statistics
# =============================================================================


def read_xvg(path: Path) -> pd.DataFrame:
    """Read GROMACS XVG-like file into a DataFrame.

    Handles files with legends, single y column, or multiple y columns. Comments and
    Grace metadata lines are skipped. The first column is named time_ps for time
    series unless the caller later treats it as residue index.
    """
    legends: List[str] = []
    rows: List[List[float]] = []

    with path.open("r", errors="ignore") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith("@"):
                if "legend" in line and '"' in line:
                    try:
                        legends.append(line.split('"')[1])
                    except Exception:
                        pass
                continue
            if line.startswith("#"):
                continue
            parts = line.split()
            try:
                rows.append([float(x) for x in parts])
            except ValueError:
                continue

    if not rows:
        return pd.DataFrame()

    width = max(len(r) for r in rows)
    padded = [r + [np.nan] * (width - len(r)) for r in rows]
    df = pd.DataFrame(padded)

    y_count = max(0, df.shape[1] - 1)
    if len(legends) >= y_count:
        y_names = legends[:y_count]
    elif y_count == 1:
        y_names = [path.stem]
    else:
        y_names = [legends[i] if i < len(legends) else f"series_{i+1}" for i in range(y_count)]

    df.columns = ["x"] + y_names
    return df


def metric_from_filename(path: Path) -> str:
    name = path.name.lower()
    if "rmsf" in name:
        return "rmsf"
    if "rmsd" in name or name.startswith("rms"):
        return "rmsd"
    if "gyr" in name or "rg" in name or "radius" in name:
        return "rg"
    if "sasa" in name or "area" in name:
        return "sasa"
    if "dssp" in name:
        return "dssp"
    if "distance" in name or any(x in name for x in ["ser", "his", "glu", "catalytic"]):
        return "distance"
    return "other"


def infer_system_name(path: Path, root: Optional[Path] = None) -> str:
    if root:
        try:
            rel = path.relative_to(root)
            if len(rel.parts) > 1:
                return rel.parts[0]
        except Exception:
            pass
    return path.parent.name if path.parent.name not in (".", "") else path.stem


def summarize_xvg(path: Path, metric: Optional[str] = None, system: Optional[str] = None) -> Dict[str, Any]:
    metric = metric or metric_from_filename(path)
    df = read_xvg(path)
    row: Dict[str, Any] = {
        "source_file": str(path),
        "system": system or infer_system_name(path),
        "metric": metric,
        "n_points": int(df.shape[0]) if not df.empty else 0,
    }
    if df.empty or df.shape[1] < 2:
        return row

    # For RMSF, x usually means residue. For other time series, convert ps to ns.
    y_cols = [c for c in df.columns if c != "x"]
    values = df[y_cols].to_numpy(dtype=float).ravel()
    values = values[np.isfinite(values)]
    if values.size == 0:
        return row

    row.update({
        "mean": float(np.nanmean(values)),
        "median": float(np.nanmedian(values)),
        "std": float(np.nanstd(values)),
        "min": float(np.nanmin(values)),
        "max": float(np.nanmax(values)),
        "q10": float(np.nanpercentile(values, 10)),
        "q90": float(np.nanpercentile(values, 90)),
    })

    if metric in {"rmsd", "rg", "sasa", "distance", "dssp"} and df.shape[0] >= 5:
        last = df.tail(max(3, int(0.2 * df.shape[0])))[y_cols].to_numpy(dtype=float).ravel()
        last = last[np.isfinite(last)]
        first = df.head(max(3, int(0.2 * df.shape[0])))[y_cols].to_numpy(dtype=float).ravel()
        first = first[np.isfinite(first)]
        if last.size:
            row["last20_mean"] = float(np.nanmean(last))
            row["last20_std"] = float(np.nanstd(last))
        if last.size and first.size:
            row["last20_minus_first20"] = float(np.nanmean(last) - np.nanmean(first))

    return row


# =============================================================================
# Plotting
# =============================================================================


def label_for_metric(metric: str) -> Tuple[str, str, str]:
    m = metric.lower()
    if m == "rmsd":
        return "Time (ns)", "RMSD (nm)", "Backbone RMSD"
    if m == "rmsf":
        return "Residue index", "RMSF (nm)", "Residue-wise RMSF"
    if m == "rg":
        return "Time (ns)", "Radius of gyration (nm)", "Radius of gyration"
    if m == "sasa":
        return "Time (ns)", "SASA (nm²)", "Solvent-accessible surface area"
    if m == "dssp":
        return "Time (ns)", "Secondary-structure count", "DSSP secondary-structure content"
    if m == "distance":
        return "Time (ns)", "Distance (nm)", "Catalytic/contact distance"
    return "x", "Value", metric


def plot_metric(paths: Sequence[Path], metric: str, out_png: Path, root: Optional[Path] = None, dpi: int = 600) -> Optional[Path]:
    if not paths:
        return None

    xlab, ylab, title = label_for_metric(metric)
    fig, ax = plt.subplots(figsize=(9.5, 5.2))
    plotted = 0

    for path in paths:
        df = read_xvg(path)
        if df.empty or df.shape[1] < 2:
            continue
        system = infer_system_name(path, root)
        x = df["x"].copy()
        if metric.lower() != "rmsf":
            x = x / 1000.0
        y_cols = [c for c in df.columns if c != "x"]
        if metric.lower() == "dssp":
            for col in y_cols:
                ax.plot(x, df[col], linewidth=1.1, label=f"{system}: {col}")
                plotted += 1
        else:
            # If multiple series in a non-DSSP XVG, plot each; otherwise one clean line per system.
            for col in y_cols:
                label = system if len(y_cols) == 1 else f"{system}: {col}"
                ax.plot(x, df[col], linewidth=1.4, label=label)
                plotted += 1

    if plotted == 0:
        plt.close(fig)
        return None

    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_title(title)
    ax.grid(True, alpha=0.25)
    if plotted <= 18:
        ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8, frameon=False)
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return out_png


def plot_all_metrics(xvg_root: Path, outdir: Path, dpi: int = 600) -> Tuple[pd.DataFrame, List[Path]]:
    xvg_files = sorted(xvg_root.rglob("*.xvg"))
    if not xvg_files:
        return pd.DataFrame(), []

    summary_rows = []
    by_metric: Dict[str, List[Path]] = {}
    for p in xvg_files:
        metric = metric_from_filename(p)
        by_metric.setdefault(metric, []).append(p)
        summary_rows.append(summarize_xvg(p, metric=metric, system=infer_system_name(p, xvg_root)))

    plots_dir = outdir / "plots"
    safe_mkdir(plots_dir)
    made: List[Path] = []
    for metric, paths in sorted(by_metric.items()):
        out_png = plots_dir / f"comparative_{metric}.png"
        result = plot_metric(paths, metric, out_png, root=xvg_root, dpi=dpi)
        if result:
            made.append(result)

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(outdir / "md_xvg_summary.tsv", sep="\t", index=False)
    return summary, made


# =============================================================================
# GROMACS command generation: refinement and MD validation
# =============================================================================


@dataclass
class RefinementConfig:
    forcefield: str = "amber99sb-ildn"
    water: str = "tip3p"
    box_type: str = "dodecahedron"
    box_distance_nm: float = 1.0
    salt_molar: float = 0.15
    position_restraint_force: int = 200
    ntomp: int = 8
    ntmpi: int = 1


def refinement_workflow_script(input_pdb: Path, mdp_dir: Path, outdir: Path, cfg: RefinementConfig) -> str:
    basename = stem_clean(input_pdb)
    return f"""#!/usr/bin/env bash
set -euo pipefail

INPUT=\"{input_pdb.resolve()}\"
MDP_DIR=\"{mdp_dir.resolve()}\"
OUTROOT=\"{outdir.resolve()}\"
BASENAME=\"{basename}\"
WORK=\"${{OUTROOT}}/refinement_work/${{BASENAME}}\"
REFINED=\"${{OUTROOT}}/refined_receptors\"
mkdir -p \"${{WORK}}\" \"${{REFINED}}\"
cd \"${{WORK}}\"

echo \"=== ${{BASENAME}}: topology ===\"
gmx pdb2gmx -f \"${{INPUT}}\" -o processed.gro -water {cfg.water} -ff {cfg.forcefield} -ignh

echo \"=== ${{BASENAME}}: box ===\"
gmx editconf -f processed.gro -o boxed.gro -c -d {cfg.box_distance_nm:.3f} -bt {cfg.box_type}

echo \"=== ${{BASENAME}}: solvation ===\"
gmx solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top

echo \"=== ${{BASENAME}}: ions ===\"
gmx grompp -f \"${{MDP_DIR}}/ions.mdp\" -c solvated.gro -p topol.top -o ions.tpr -maxwarn 1
echo \"SOL\" | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc {cfg.salt_molar:.3f}

echo \"=== ${{BASENAME}}: backbone restraints ===\"
echo \"q\" | gmx make_ndx -f processed.gro -o index.ndx
# Default group 4 is commonly Backbone after pdb2gmx; inspect index.ndx if this differs.
echo \"4\" | gmx genrestr -f processed.gro -n index.ndx -o posre_backbone.itp -fc {cfg.position_restraint_force} {cfg.position_restraint_force} {cfg.position_restraint_force}

if ! grep -q \"POSRES_BB\" topol.top; then
  sed -i '/Include Position restraint file/a \\\n#ifdef POSRES_BB\\n#include "posre_backbone.itp"\\n#endif' topol.top
fi

echo \"=== ${{BASENAME}}: restrained minimization ===\"
gmx grompp -f \"${{MDP_DIR}}/em.mdp\" -c solv_ions.gro -r solv_ions.gro -p topol.top -o em.tpr
gmx mdrun -deffnm em -v -ntmpi {cfg.ntmpi} -ntomp {cfg.ntomp}
grep \"Maximum force\" em.log || true

echo \"=== ${{BASENAME}}: optional restrained short MD ===\"
if [[ -f \"${{MDP_DIR}}/md_relax.mdp\" ]]; then
  gmx grompp -f \"${{MDP_DIR}}/md_relax.mdp\" -c em.gro -r em.gro -p topol.top -o md_relax.tpr -maxwarn 1
  gmx mdrun -deffnm md_relax -v -ntmpi {cfg.ntmpi} -ntomp {cfg.ntomp}
  RELAXED=md_relax.gro
else
  RELAXED=em.gro
fi

echo \"=== ${{BASENAME}}: final restrained minimization ===\"
gmx grompp -f \"${{MDP_DIR}}/final_em.mdp\" -c \"${{RELAXED}}\" -r \"${{RELAXED}}\" -p topol.top -o final_em.tpr -maxwarn 1
gmx mdrun -deffnm final_em -v -ntmpi {cfg.ntmpi} -ntomp {cfg.ntomp}

echo \"=== ${{BASENAME}}: export receptor ===\"
printf \"1\\n1\\n\" | gmx trjconv -f final_em.gro -s final_em.tpr -o \"${{REFINED}}/${{BASENAME}}_refined_receptor.pdb\" -pbc mol -center

echo \"Refined receptor written: ${{REFINED}}/${{BASENAME}}_refined_receptor.pdb\"
"""


def write_refinement_module(input_dir: Path, mdp_dir: Path, outdir: Path, cfg: RefinementConfig, execute: bool = False) -> List[Path]:
    structures = [p for p in sorted(input_dir.glob("*.pdb"))]
    if not structures:
        raise SystemExit(f"No PDB files found in {input_dir}")

    script_dir = outdir / "scripts" / "refinement"
    safe_mkdir(script_dir)
    script_paths = []
    for pdb in structures:
        script = refinement_workflow_script(pdb, mdp_dir, outdir, cfg)
        sp = script_dir / f"refine_{stem_clean(pdb)}.sh"
        write_text(sp, script, executable=True)
        script_paths.append(sp)

    batch_lines = ["#!/usr/bin/env bash", "set -euo pipefail"]
    for sp in script_paths:
        batch_lines.append(f"bash {sp.resolve()}")
    batch = outdir / "scripts" / "run_all_refinement.sh"
    write_text(batch, "\n".join(batch_lines) + "\n", executable=True)

    if execute:
        run_cmd(["bash", str(batch)])

    return [batch] + script_paths


@dataclass
class CEHConfig:
    ser_resid: int = 207
    his_resid: int = 508
    acid_resid: int = 388
    alt_acid_resid: int = 206
    ser_atom: str = "OG"
    his_atom_ser: str = "NE2"
    his_atom_acid: str = "ND1"
    acid_atom1: str = "OE1"
    acid_atom2: str = "OE2"
    alt_acid_atom1: str = "OE1"
    alt_acid_atom2: str = "OE2"


def md_validation_script(md_dir: Path, outdir: Path, mode: str, ceh: CEHConfig) -> str:
    # Uses interactive selections compatible with common GROMACS index groups; users can edit if their groups differ.
    ceh_block = ""
    if mode.lower() in {"ceh", "hydrolase"}:
        ceh_block = f"""

echo \"=== CEH catalytic index and distances ===\"
gmx make_ndx -f md.tpr -o catalytic.ndx <<'EOF'
a {ceh.ser_atom} & r {ceh.ser_resid}
name 17 SER{ceh.ser_resid}_{ceh.ser_atom}
a {ceh.his_atom_ser} & r {ceh.his_resid}
name 18 HIS{ceh.his_resid}_{ceh.his_atom_ser}
a {ceh.his_atom_acid} & r {ceh.his_resid}
name 19 HIS{ceh.his_resid}_{ceh.his_atom_acid}
a {ceh.acid_atom1} & r {ceh.acid_resid}
name 20 ACID{ceh.acid_resid}_{ceh.acid_atom1}
a {ceh.acid_atom2} & r {ceh.acid_resid}
name 21 ACID{ceh.acid_resid}_{ceh.acid_atom2}
a {ceh.alt_acid_atom1} & r {ceh.alt_acid_resid}
name 22 ALTACID{ceh.alt_acid_resid}_{ceh.alt_acid_atom1}
a {ceh.alt_acid_atom2} & r {ceh.alt_acid_resid}
name 23 ALTACID{ceh.alt_acid_resid}_{ceh.alt_acid_atom2}
q
EOF

gmx distance -s md.tpr -f md_fit.xtc -n catalytic.ndx -oall ser{ceh.ser_resid}_his{ceh.his_resid}.xvg -select 'group "SER{ceh.ser_resid}_{ceh.ser_atom}" plus group "HIS{ceh.his_resid}_{ceh.his_atom_ser}"'
gmx distance -s md.tpr -f md_fit.xtc -n catalytic.ndx -oall his{ceh.his_resid}_acid{ceh.acid_resid}_{ceh.acid_atom1}.xvg -select 'group "HIS{ceh.his_resid}_{ceh.his_atom_acid}" plus group "ACID{ceh.acid_resid}_{ceh.acid_atom1}"'
gmx distance -s md.tpr -f md_fit.xtc -n catalytic.ndx -oall his{ceh.his_resid}_acid{ceh.acid_resid}_{ceh.acid_atom2}.xvg -select 'group "HIS{ceh.his_resid}_{ceh.his_atom_acid}" plus group "ACID{ceh.acid_resid}_{ceh.acid_atom2}"'
gmx distance -s md.tpr -f md_fit.xtc -n catalytic.ndx -oall his{ceh.his_resid}_altacid{ceh.alt_acid_resid}_{ceh.alt_acid_atom1}.xvg -select 'group "HIS{ceh.his_resid}_{ceh.his_atom_ser}" plus group "ALTACID{ceh.alt_acid_resid}_{ceh.alt_acid_atom1}"' || true
gmx distance -s md.tpr -f md_fit.xtc -n catalytic.ndx -oall his{ceh.his_resid}_altacid{ceh.alt_acid_resid}_{ceh.alt_acid_atom2}.xvg -select 'group "HIS{ceh.his_resid}_{ceh.his_atom_ser}" plus group "ALTACID{ceh.alt_acid_resid}_{ceh.alt_acid_atom2}"' || true
"""

    return f"""#!/usr/bin/env bash
set -euo pipefail
cd \"{md_dir.resolve()}\"
mkdir -p \"{outdir.resolve()}/md_analysis/{md_dir.name}\"
ANALYSIS_OUT=\"{outdir.resolve()}/md_analysis/{md_dir.name}\"

echo \"=== Trajectory centering ===\"
printf \"Protein\\nSystem\\n\" | gmx trjconv -s md.tpr -f md.xtc -o centered.xtc -pbc mol -center

echo \"=== Protein fitting ===\"
printf \"Protein\\nSystem\\n\" | gmx trjconv -s md.tpr -f centered.xtc -o md_fit.xtc -fit rot+trans

echo \"=== RMSD ===\"
printf \"Backbone\\nBackbone\\n\" | gmx rms -s md.tpr -f md_fit.xtc -o rmsd_backbone.xvg

echo \"=== RMSF ===\"
printf \"C-alpha\\n\" | gmx rmsf -s md.tpr -f md_fit.xtc -o rmsf_ca.xvg -res

echo \"=== Radius of gyration ===\"
printf \"Protein\\n\" | gmx gyrate -s md.tpr -f md_fit.xtc -o gyrate.xvg

echo \"=== SASA ===\"
gmx sasa -s md.tpr -f md_fit.xtc -o sasa_total.xvg -surface 'Protein' -output 'Protein'
gmx sasa -s md.tpr -f md_fit.xtc -or residue_sasa.xvg -surface 'Protein' -output 'Protein' || true

echo \"=== DSSP ===\"
gmx dssp -s md.tpr -f centered.xtc -o dssp.dat -num dssp_count.xvg -dt 10 -sel 'Protein' -hmode dssp -clear || true
{ceh_block}

cp -f *.xvg \"${{ANALYSIS_OUT}}/\" 2>/dev/null || true
cp -f dssp.dat \"${{ANALYSIS_OUT}}/\" 2>/dev/null || true
echo \"MD validation files copied to ${{ANALYSIS_OUT}}\"
"""


def write_md_validation_module(md_root: Path, outdir: Path, mode: str, ceh: CEHConfig, execute: bool = False) -> List[Path]:
    candidates = []
    if (md_root / "md.tpr").exists() and (md_root / "md.xtc").exists():
        candidates = [md_root]
    else:
        for d in sorted([p for p in md_root.iterdir() if p.is_dir()] if md_root.exists() else []):
            if (d / "md.tpr").exists() and (d / "md.xtc").exists():
                candidates.append(d)
    if not candidates:
        raise SystemExit(f"No MD run directories containing md.tpr and md.xtc found in {md_root}")

    script_dir = outdir / "scripts" / "md_validation"
    safe_mkdir(script_dir)
    scripts = []
    for d in candidates:
        sp = script_dir / f"analyze_{d.name}.sh"
        write_text(sp, md_validation_script(d, outdir, mode, ceh), executable=True)
        scripts.append(sp)

    batch = outdir / "scripts" / "run_all_md_validation.sh"
    lines = ["#!/usr/bin/env bash", "set -euo pipefail"] + [f"bash {s.resolve()}" for s in scripts]
    write_text(batch, "\n".join(lines) + "\n", executable=True)
    if execute:
        run_cmd(["bash", str(batch)])
    return [batch] + scripts


# =============================================================================
# Selection engines
# =============================================================================


def robust_scale(values: pd.Series, direction: str = "higher") -> pd.Series:
    x = to_numeric_series(values)
    if x.notna().sum() == 0:
        return pd.Series(np.nan, index=x.index)
    med = x.median()
    x = x.fillna(med)
    q5, q95 = np.nanpercentile(x.to_numpy(dtype=float), [5, 95])
    if not np.isfinite(q5) or not np.isfinite(q95) or q5 == q95:
        return pd.Series(0.5, index=x.index)
    clipped = x.clip(q5, q95)
    s = (clipped - q5) / (q95 - q5)
    if direction.lower().startswith("lower"):
        s = 1.0 - s
    return s.clip(0, 1)


def score_abs_closer_zero(values: pd.Series) -> pd.Series:
    d = to_numeric_series(values).abs()
    if d.notna().sum() == 0:
        return pd.Series(np.nan, index=d.index)
    d = d.fillna(d.median())
    q95 = np.nanpercentile(d.to_numpy(dtype=float), 95)
    if not np.isfinite(q95) or q95 == 0:
        return pd.Series(0.5, index=d.index)
    return (1.0 - d / q95).clip(0, 1)


def weighted_available(scores: Dict[str, pd.Series], weights: Dict[str, float], index: pd.Index) -> pd.Series:
    cols = [k for k in weights if k in scores and scores[k].notna().any()]
    if not cols:
        return pd.Series(0.5, index=index)
    wsum = float(sum(weights[k] for k in cols))
    out = pd.Series(0.0, index=index)
    for k in cols:
        out += (weights[k] / wsum) * scores[k].fillna(0.5)
    return out


def validation_gate(df: pd.DataFrame, mode: str) -> pd.Series:
    gate = pd.Series(1.0, index=df.index)

    # Stricter gate for CEH/hydrolase active-site structures.
    thresholds = {
        "general": {
            "molprobity_score": (3.0, "le"),
            "clashscore": (25.0, "le"),
            "rama_outliers_%": (2.0, "le"),
            "cablam_severe_%": (2.0, "le"),
        },
        "ceh": {
            "molprobity_score": (2.5, "le"),
            "clashscore": (15.0, "le"),
            "rama_outliers_%": (1.0, "le"),
            "cablam_severe_%": (1.0, "le"),
        },
    }
    cfg = thresholds["ceh"] if mode.lower() in {"ceh", "hydrolase"} else thresholds["general"]
    for col, (thr, op) in cfg.items():
        if col not in df.columns:
            continue
        x = to_numeric_series(df[col])
        ok = (x <= thr) if op == "le" else (x >= thr)
        ok = ok | x.isna()
        gate *= ok.astype(float)
    return gate


def add_metric_score(scores: Dict[str, pd.Series], df: pd.DataFrame, aliases: Sequence[str], name: str, direction: str) -> None:
    col = first_col(df, aliases)
    if col is not None:
        scores[name] = robust_scale(df[col], direction)


def select_representatives(metrics_table: Path, outdir: Path, mode: str = "general", top_n: int = 3, sheet: Optional[str] = None) -> pd.DataFrame:
    if metrics_table.suffix.lower() in {".xlsx", ".xls"}:
        df = pd.read_excel(metrics_table, sheet_name=sheet or 0, engine="openpyxl")
    else:
        # Auto-detect delimiter.
        df = pd.read_csv(metrics_table, sep=None, engine="python")

    if df.empty:
        raise SystemExit(f"Empty metrics table: {metrics_table}")

    name_col = first_col(df, ["model", "name", "Short_name", "structure", "system", "file", "source_file"])
    if name_col is None:
        df["model"] = [f"model_{i+1}" for i in range(df.shape[0])]
        name_col = "model"

    gate = validation_gate(df, mode)
    df["validation_gate"] = gate

    S: Dict[str, pd.Series] = {}
    # Validation/stereochemistry
    add_metric_score(S, df, ["molprobity_score"], "molprobity", "lower")
    add_metric_score(S, df, ["clashscore"], "clashscore", "lower")
    add_metric_score(S, df, ["rama_outliers_%", "rama_outliers", "Ramachandran_outliers"], "rama_outliers", "lower")
    add_metric_score(S, df, ["rama_favored_%", "rama_favored"], "rama_favored", "higher")
    add_metric_score(S, df, ["rotamer_outliers_%", "rotamer_outliers"], "rotamer_outliers", "lower")
    add_metric_score(S, df, ["cablam_severe_%", "cablam_severe"], "cablam_severe", "lower")
    add_metric_score(S, df, ["cablam_outliers_%", "cablam_outliers"], "cablam_outliers", "lower")
    add_metric_score(S, df, ["voromqa_light_score", "voro_mean", "voro_dark", "voromqa_dark_score"], "packing", "higher")

    # MD stability and fold agreement
    add_metric_score(S, df, ["rmsd_mean", "backbone_rmsd_mean", "RMSD", "global_rmsd_A", "global_rmsd"], "global_rmsd", "lower")
    add_metric_score(S, df, ["rmsd_last20_std", "rmsd_std", "RMSD_std"], "rmsd_stability", "lower")
    add_metric_score(S, df, ["rmsf_mean", "mean_rmsf", "RMSF_mean"], "rmsf", "lower")
    add_metric_score(S, df, ["rg_std", "gyrate_std", "rg_last20_std", "Rg_std"], "rg_stability", "lower")
    add_metric_score(S, df, ["sasa_std", "sasa_last20_std", "SASA_std"], "sasa_stability", "lower")
    add_metric_score(S, df, ["TM_score", "alntmscore", "tm_score"], "tm_score", "higher")
    add_metric_score(S, df, ["Q_score", "q_score"], "q_score", "higher")
    add_metric_score(S, df, ["CMO", "contact_map_overlap"], "cmo", "higher")

    # CEH/hydrolase-specific columns.
    add_metric_score(S, df, ["pocket_rmsd_A", "pocket_rmsd", "active_site_rmsd_A"], "pocket_rmsd", "lower")
    add_metric_score(S, df, ["HGG_plddt_mean", "HGG_mean_plddt", "hgg_plddt"], "hgg_plddt", "higher")
    add_metric_score(S, df, ["GESAG_plddt_mean", "GxSxG_plddt_mean", "gxsxg_plddt", "nucleophile_loop_plddt"], "nuc_loop_plddt", "higher")
    add_metric_score(S, df, ["distance_ser_his_std", "ser_his_std", "catalytic_distance_std"], "catalytic_distance_stability", "lower")

    delta_col = first_col(df, ["catalytic_distance_delta_A", "distance_delta_A", "ser_his_delta", "active_site_distance_delta_A"])
    if delta_col is not None:
        S["catalytic_distance_delta"] = score_abs_closer_zero(df[delta_col])

    # Components.
    validation_weights = {
        "molprobity": 0.22, "clashscore": 0.16, "rama_outliers": 0.14,
        "rama_favored": 0.08, "rotamer_outliers": 0.08, "cablam_severe": 0.18,
        "cablam_outliers": 0.08, "packing": 0.06,
    }
    md_weights = {
        "global_rmsd": 0.25, "rmsd_stability": 0.30, "rmsf": 0.20,
        "rg_stability": 0.15, "sasa_stability": 0.10,
    }
    fold_weights = {"tm_score": 0.50, "q_score": 0.25, "cmo": 0.25}
    ceh_geometry_weights = {
        "pocket_rmsd": 0.32,
        "catalytic_distance_delta": 0.26,
        "catalytic_distance_stability": 0.16,
        "hgg_plddt": 0.13,
        "nuc_loop_plddt": 0.13,
    }

    df["S_validation"] = weighted_available(S, validation_weights, df.index)
    df["S_md_stability"] = weighted_available(S, md_weights, df.index)
    df["S_fold_agreement"] = weighted_available(S, fold_weights, df.index)

    if mode.lower() in {"ceh", "hydrolase"}:
        df["S_ceh_geometry"] = weighted_available(S, ceh_geometry_weights, df.index)
        df["SelectionScore"] = gate * (
            0.34 * df["S_validation"] +
            0.18 * df["S_md_stability"] +
            0.10 * df["S_fold_agreement"] +
            0.38 * df["S_ceh_geometry"]
        )
        df["selection_engine"] = "CEH/hydrolase catalytic-geometry engine"
    else:
        df["SelectionScore"] = gate * (
            0.45 * df["S_validation"] +
            0.30 * df["S_md_stability"] +
            0.25 * df["S_fold_agreement"]
        )
        df["selection_engine"] = "general receptor docking-readiness engine"

    df["SelectionScore_100"] = 100.0 * df["SelectionScore"]
    df = df.sort_values("SelectionScore", ascending=False).reset_index(drop=True)
    df["selection_rank"] = np.arange(1, len(df) + 1)
    df["representative_selected"] = df["selection_rank"] <= top_n

    safe_mkdir(outdir)
    df.to_csv(outdir / "all_receptors_scored.tsv", sep="\t", index=False)
    df[df["representative_selected"]].to_csv(outdir / "representative_receptors.tsv", sep="\t", index=False)

    with pd.ExcelWriter(outdir / "receptor_selection_workbook.xlsx", engine="openpyxl") as writer:
        df.to_excel(writer, sheet_name="ALL_SCORED", index=False)
        df[df["representative_selected"]].to_excel(writer, sheet_name="REPRESENTATIVES", index=False)

    plot_selection_scores(df, outdir, name_col=name_col, mode=mode)
    write_selection_report(df, outdir, name_col=name_col, mode=mode, top_n=top_n)
    return df


def plot_selection_scores(df: pd.DataFrame, outdir: Path, name_col: str, mode: str) -> None:
    plots_dir = outdir / "plots"
    safe_mkdir(plots_dir)
    top = df.head(min(15, len(df))).copy()
    labels = top[name_col].astype(str).tolist()[::-1]
    scores = top["SelectionScore_100"].to_numpy(dtype=float)[::-1]

    fig, ax = plt.subplots(figsize=(9, max(4, 0.35 * len(top))))
    ax.barh(labels, scores)
    ax.set_xlabel("Selection score")
    ax.set_title("Top receptor models for docking")
    ax.set_xlim(0, 100)
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(plots_dir / "top_receptor_selection_scores.png", dpi=600, bbox_inches="tight")
    plt.close(fig)

    if mode.lower() in {"ceh", "hydrolase"} and "S_ceh_geometry" in df.columns:
        fig, ax = plt.subplots(figsize=(6.5, 5.2))
        ax.scatter(df["S_validation"] * 100, df["S_ceh_geometry"] * 100, s=40, alpha=0.85)
        for _, row in df.head(8).iterrows():
            ax.text(row["S_validation"] * 100, row["S_ceh_geometry"] * 100, str(row[name_col])[:18], fontsize=7)
        ax.set_xlabel("Validation component")
        ax.set_ylabel("CEH catalytic-geometry component")
        ax.set_title("CEH receptor selection: validation vs catalytic geometry")
        ax.set_xlim(0, 100)
        ax.set_ylim(0, 100)
        ax.grid(True, alpha=0.25)
        fig.tight_layout()
        fig.savefig(plots_dir / "ceh_validation_vs_catalytic_geometry.png", dpi=600, bbox_inches="tight")
        plt.close(fig)


def write_selection_report(df: pd.DataFrame, outdir: Path, name_col: str, mode: str, top_n: int) -> None:
    top = df.head(top_n)
    lines = []
    lines.append("# Receptor docking-preparation selection report")
    lines.append("")
    lines.append(f"Selection engine: **{top['selection_engine'].iloc[0]}**")
    lines.append(f"Total candidates scored: {df.shape[0]}")
    lines.append(f"Representatives selected: {top_n}")
    lines.append("")
    lines.append("## Selected representatives")
    lines.append("")
    lines.append("| Rank | Model | Selection score | Validation gate |")
    lines.append("|---:|---|---:|---:|")
    for _, row in top.iterrows():
        lines.append(f"| {int(row['selection_rank'])} | {row[name_col]} | {row['SelectionScore_100']:.2f} | {row['validation_gate']:.0f} |")
    lines.append("")
    lines.append("## Notes")
    if mode.lower() in {"ceh", "hydrolase"}:
        lines.append("The CEH/hydrolase engine prioritizes local catalytic-pocket geometry, catalytic-distance preservation, and motif-region reliability in addition to global validation metrics.")
    else:
        lines.append("The general engine prioritizes stereochemical validation, MD stability, fold/template agreement, and packing quality.")
    write_text(outdir / "receptor_selection_report.md", "\n".join(lines) + "\n")


# =============================================================================
# Representative frame extraction
# =============================================================================


def write_frame_extraction_script(selection_table: Path, md_root: Path, outdir: Path, top_n: int = 3, time_col: Optional[str] = None) -> Path:
    df = pd.read_csv(selection_table, sep=None, engine="python")
    if df.empty:
        raise SystemExit(f"Empty selection table: {selection_table}")
    system_col = first_col(df, ["system", "model", "name", "Short_name"])
    if system_col is None:
        raise SystemExit("Selection table must include system/model/name column for frame extraction")
    time_col = time_col or first_col(df, ["time_ps", "time", "frame_time_ps", "representative_time_ps", "dump_ps"])
    if time_col is None:
        raise SystemExit("Selection table must include time_ps or representative_time_ps for frame extraction")

    top = df.head(top_n)
    lines = ["#!/usr/bin/env bash", "set -euo pipefail", f"OUTDIR=\"{(outdir / 'representative_receptors').resolve()}\"", "mkdir -p \"$OUTDIR\""]
    for _, row in top.iterrows():
        system = str(row[system_col])
        time_ps = float(row[time_col])
        md_dir = md_root / system
        if not md_dir.exists():
            md_dir = md_root
        out_pdb = f"${{OUTDIR}}/{re.sub(r'[^A-Za-z0-9_.-]+', '_', system)}_{int(round(time_ps))}ps_receptor.pdb"
        lines.append(f"# Representative for {system} at {time_ps:.1f} ps")
        lines.append(f"printf \"Protein\\n\" | gmx trjconv -s \"{(md_dir / 'md.tpr').resolve()}\" -f \"{(md_dir / 'md_fit.xtc').resolve()}\" -o \"{out_pdb}\" -dump {time_ps:.3f}")
    script = outdir / "scripts" / "extract_representative_frames.sh"
    write_text(script, "\n".join(lines) + "\n", executable=True)
    return script


# =============================================================================
# Main CLI
# =============================================================================


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Integrated receptor preparation pipeline for docking: refinement, MD validation, plotting, and representative selection."
    )
    p.add_argument("--mode", choices=["general", "ceh", "hydrolase"], default="general",
                   help="Use general or CEH/hydrolase-specific workflow.")
    p.add_argument("--steps", nargs="+", default=["all"],
                   choices=["all", "refine", "md-validation", "plot", "select", "extract-frames"],
                   help="Pipeline modules to run or generate.")
    p.add_argument("--outdir", default="receptor_docking_prep_results", help="Output directory.")

    # Refinement inputs
    p.add_argument("--input-dir", default=None, help="Directory of input PDB receptor structures for refinement.")
    p.add_argument("--mdp-dir", default="mdp", help="Directory containing ions.mdp, em.mdp, md_relax.mdp, final_em.mdp.")
    p.add_argument("--forcefield", default="amber99sb-ildn")
    p.add_argument("--water", default="tip3p")
    p.add_argument("--box-type", default="dodecahedron")
    p.add_argument("--box-distance-nm", type=float, default=1.0)
    p.add_argument("--salt-molar", type=float, default=0.15)
    p.add_argument("--posres-force", type=int, default=200)
    p.add_argument("--ntomp", type=int, default=8)
    p.add_argument("--ntmpi", type=int, default=1)

    # MD validation inputs
    p.add_argument("--md-root", default=None, help="MD run directory or parent directory containing md.tpr/md.xtc runs.")
    p.add_argument("--xvg-root", default=None, help="Directory containing XVG files for plotting/summary.")
    p.add_argument("--plot-dpi", type=int, default=600)

    # CEH options
    p.add_argument("--ceh-ser", type=int, default=207)
    p.add_argument("--ceh-his", type=int, default=508)
    p.add_argument("--ceh-acid", type=int, default=388)
    p.add_argument("--ceh-alt-acid", type=int, default=206)

    # Selection inputs
    p.add_argument("--metrics-table", default=None, help="CSV/TSV/XLSX table containing validation/MD/fold metrics for representative selection.")
    p.add_argument("--sheet", default=None, help="Excel sheet name/index for --metrics-table.")
    p.add_argument("--top-n", type=int, default=3, help="Number of representative receptors to select.")
    p.add_argument("--selection-table", default=None, help="Representative table for frame extraction.")
    p.add_argument("--time-col", default=None, help="Column containing representative frame time in ps for extraction.")

    p.add_argument("--execute", action="store_true", help="Execute generated shell scripts. Default is to write scripts only.")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    safe_mkdir(outdir)

    steps = set(args.steps)
    if "all" in steps:
        steps = {"refine", "md-validation", "plot", "select"}

    manifest: Dict[str, Any] = {"mode": args.mode, "steps": sorted(steps), "outputs": []}

    ceh_cfg = CEHConfig(
        ser_resid=args.ceh_ser,
        his_resid=args.ceh_his,
        acid_resid=args.ceh_acid,
        alt_acid_resid=args.ceh_alt_acid,
    )

    if "refine" in steps:
        if not args.input_dir:
            print("[SKIP] refine: --input-dir was not provided")
        else:
            cfg = RefinementConfig(
                forcefield=args.forcefield,
                water=args.water,
                box_type=args.box_type,
                box_distance_nm=args.box_distance_nm,
                salt_molar=args.salt_molar,
                position_restraint_force=args.posres_force,
                ntomp=args.ntomp,
                ntmpi=args.ntmpi,
            )
            paths = write_refinement_module(Path(args.input_dir), Path(args.mdp_dir), outdir, cfg, execute=args.execute)
            manifest["outputs"].extend([str(p) for p in paths])
            print(f"[OK] Refinement scripts written: {len(paths)}")

    if "md-validation" in steps:
        if not args.md_root:
            print("[SKIP] md-validation: --md-root was not provided")
        else:
            paths = write_md_validation_module(Path(args.md_root), outdir, args.mode, ceh_cfg, execute=args.execute)
            manifest["outputs"].extend([str(p) for p in paths])
            print(f"[OK] MD validation scripts written: {len(paths)}")

    if "plot" in steps:
        xvg_root = Path(args.xvg_root) if args.xvg_root else (outdir / "md_analysis")
        if not xvg_root.exists():
            print(f"[SKIP] plot: XVG root not found: {xvg_root}")
        else:
            summary, plots = plot_all_metrics(xvg_root, outdir, dpi=args.plot_dpi)
            manifest["outputs"].append(str(outdir / "md_xvg_summary.tsv"))
            manifest["outputs"].extend([str(p) for p in plots])
            print(f"[OK] XVG summary rows: {summary.shape[0]}; plots: {len(plots)}")

    if "select" in steps:
        if not args.metrics_table:
            print("[SKIP] select: --metrics-table was not provided")
        else:
            selected = select_representatives(Path(args.metrics_table), outdir, mode=args.mode, top_n=args.top_n, sheet=args.sheet)
            manifest["outputs"].extend([
                str(outdir / "all_receptors_scored.tsv"),
                str(outdir / "representative_receptors.tsv"),
                str(outdir / "receptor_selection_workbook.xlsx"),
                str(outdir / "receptor_selection_report.md"),
            ])
            print(f"[OK] Representative selection complete: {min(args.top_n, selected.shape[0])} selected")

    if "extract-frames" in steps:
        table = Path(args.selection_table) if args.selection_table else (outdir / "representative_receptors.tsv")
        if not table.exists():
            print(f"[SKIP] extract-frames: selection table not found: {table}")
        elif not args.md_root:
            print("[SKIP] extract-frames: --md-root was not provided")
        else:
            script = write_frame_extraction_script(table, Path(args.md_root), outdir, top_n=args.top_n, time_col=args.time_col)
            manifest["outputs"].append(str(script))
            print(f"[OK] Frame extraction script written: {script}")
            if args.execute:
                run_cmd(["bash", str(script)])

    write_text(outdir / "pipeline_manifest.json", json.dumps(manifest, indent=2) + "\n")
    print(f"[DONE] Manifest written: {outdir / 'pipeline_manifest.json'}")


if __name__ == "__main__":
    main()
