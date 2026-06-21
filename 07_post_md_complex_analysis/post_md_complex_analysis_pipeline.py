#!/usr/bin/env python3
"""
Post-MD receptor-ligand complex analysis pipeline.

Modules
-------
1) workbook
   Comparative analysis from multi-sheet Excel workbooks generated from GROMACS/XVG outputs.
   Supports general receptor-ligand complexes and CEH-like catalytic-state analysis.

2) contacts
   Direct residue-contact occupancy analysis from trajectory/topology using MDAnalysis.

3) xvg-collect
   Collect GROMACS XVG outputs from multiple systems into a standardized workbook.

4) make-gromacs-scripts
   Generate reproducible GROMACS analysis shell scripts for RMSD, ligand RMSD, Rg, RMSF,
   SASA, DSSP, contacts, ligand-pocket COM distance, and CEH catalytic distances.

The CEH-like engine emphasizes catalytic persistence, including Ser207-His508,
His508-acidic residue, Ser207-ligand carbonyl, and Gly128/Gly129 oxyanion-hole distances.
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
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple, Any

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    import MDAnalysis as mda
    from MDAnalysis.analysis.distances import distance_array
    HAS_MDA = True
except Exception:
    mda = None
    distance_array = None
    HAS_MDA = False

# =============================================================================
# Defaults
# =============================================================================

INDEX_COLUMNS = {
    "Frame", "Time (ps)", "Time_ps", "Time_ns", "Residue#", "Residue_Number",
    "Residue_Name", "Residue_ID", "Atom", "Index"
}

DEFAULT_CEH_CATALYTIC_DISTANCES = [
    "ser207_his508",
    "his508_glu388_oe2",
    "ser207_carbonylC_distance",
    "gly128N_O2_distance",
    "gly129N_O2_distance",
]

DEFAULT_CONTACT_COLUMNS = [
    "Residue_Number",
    "Residue_Name",
    "Occupancy_Percent",
    "Minimum_Distance",
    "Maximum_Distance",
]

CEH_KEY_CONTACT_RESIDUES = [127, 128, 129, 130, 206, 207, 508]
CEH_POCKET_RESIDUES = [128, 129, 130, 206, 207, 242, 349, 391, 508]

PLOT_METRIC_SHEETS = {
    "C_alpha_RMSD": {"ylabel": "Cα RMSD (nm)", "title": "Cα RMSD", "stem": "01_Calpha_RMSD_timeseries"},
    "Lig_RMSD": {"ylabel": "Ligand RMSD (nm)", "title": "Ligand RMSD", "stem": "02_ligand_RMSD_timeseries"},
    "rGyr": {"ylabel": "Radius of gyration (nm)", "title": "Radius of gyration", "stem": "03_rGyr_timeseries"},
    "SASA": {"ylabel": "SASA (nm²)", "title": "Solvent-accessible surface area", "stem": "04_SASA_timeseries"},
    "Protein_SASA": {"ylabel": "SASA (nm²)", "title": "Protein SASA", "stem": "04_Protein_SASA_timeseries"},
    "Pocket_COM": {"ylabel": "COM distance (nm)", "title": "Ligand-pocket COM distance", "stem": "05_pocket_COM_distance"},
}

# =============================================================================
# General utilities
# =============================================================================

def clean_column_name(value: object) -> str:
    text = str(value).strip()
    text = re.sub(r"\.\d+$", "", text)
    text = text.replace("Å", "A")
    return text


def safe_name(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(text)).strip("_") or "unnamed"


def system_label(system: str) -> str:
    s = str(system)
    # CEH_CE20_5 -> CE20:5, RP4_CE20_5 -> RP4-CE20:5
    m = re.search(r"(RP\d+|CEH|NPC\d*)[_-]?CE(\d+)[_:](\d+)", s, flags=re.I)
    if not m:
        m = re.search(r"(RP\d+|CEH|NPC\d*)[_-]?CE(\d+)_(\d+)", s, flags=re.I)
    if m:
        return f"{m.group(1).upper()}–CE{m.group(2)}:{m.group(3)}"
    m = re.search(r"CE(\d+)[_:](\d+)", s, flags=re.I)
    if not m:
        m = re.search(r"CE(\d+)_(\d+)", s, flags=re.I)
    if m:
        return f"CE{m.group(1)}:{m.group(2)}"
    return s


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def build_color_map(systems: Iterable[str]) -> Dict[str, tuple]:
    unique = list(dict.fromkeys([str(s) for s in systems]))
    palette = list(plt.cm.tab10.colors) + list(plt.cm.Set2.colors) + list(plt.cm.Dark2.colors) + list(plt.cm.Paired.colors)
    return {system: palette[i % len(palette)] for i, system in enumerate(unique)}


def get_system_color(system: str, color_map: Dict[str, tuple]) -> tuple:
    return color_map.get(str(system), plt.cm.tab10.colors[0])


def save_figure(fig: plt.Figure, figdir: Path, stem: str, dpi: int, formats: Iterable[str]) -> None:
    for fmt in formats:
        fmt = fmt.strip().lower()
        if not fmt:
            continue
        kwargs = {"bbox_inches": "tight"}
        if fmt in {"png", "jpg", "jpeg", "tif", "tiff"}:
            kwargs["dpi"] = dpi
        fig.savefig(figdir / f"{stem}.{fmt}", **kwargs)
    plt.close(fig)


def numeric_cols(df: pd.DataFrame, exclude: Iterable[str] = INDEX_COLUMNS) -> List[str]:
    excluded = set(str(x) for x in exclude)
    return [c for c in df.columns if str(c) not in excluded and not str(c).startswith("Unnamed")]


def robust_numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def stats_dict(values: pd.Series) -> Dict[str, Any]:
    x = pd.to_numeric(values, errors="coerce").dropna()
    if x.empty:
        return {"N": 0, "Mean": np.nan, "SD": np.nan, "SEM": np.nan, "Median": np.nan, "Q1": np.nan, "Q3": np.nan, "Min": np.nan, "Max": np.nan}
    return {
        "N": int(x.size),
        "Mean": float(x.mean()),
        "SD": float(x.std(ddof=1)) if x.size > 1 else 0.0,
        "SEM": float(x.sem(ddof=1)) if x.size > 1 else 0.0,
        "Median": float(x.median()),
        "Q1": float(x.quantile(0.25)),
        "Q3": float(x.quantile(0.75)),
        "Min": float(x.min()),
        "Max": float(x.max()),
    }

# =============================================================================
# XVG parsing and collection
# =============================================================================

def read_xvg(path: Path) -> pd.DataFrame:
    legends: List[str] = []
    data: List[List[float]] = []
    with open(path, "r", errors="ignore") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith("@") and "legend" in line:
                m = re.search(r'"([^"]+)"', line)
                if m:
                    legends.append(m.group(1))
                continue
            if line.startswith(("#", "@")):
                continue
            parts = line.split()
            try:
                data.append([float(x) for x in parts])
            except ValueError:
                continue
    if not data:
        return pd.DataFrame()
    df = pd.DataFrame(data)
    if len(legends) >= df.shape[1] - 1:
        cols = ["Time_ps"] + legends[:df.shape[1] - 1]
    elif df.shape[1] == 2:
        cols = ["Time_ps", path.stem]
    else:
        cols = ["Time_ps"] + [f"value_{i}" for i in range(1, df.shape[1])]
    df.columns = cols
    df["Time_ns"] = df["Time_ps"] / 1000.0
    return df


def collect_xvg_to_workbook(input_dir: Path, out_xlsx: Path, pattern: str = "*.xvg", recursive: bool = True) -> None:
    iterator = input_dir.rglob(pattern) if recursive else input_dir.glob(pattern)
    files = sorted(p for p in iterator if p.is_file())
    if not files:
        raise SystemExit(f"No XVG files matching {pattern!r} found under {input_dir}")
    with pd.ExcelWriter(out_xlsx, engine="xlsxwriter") as writer:
        index_rows = []
        for path in files:
            df = read_xvg(path)
            if df.empty:
                continue
            # Path convention: system/metric.xvg or metric_system.xvg
            rel = path.relative_to(input_dir)
            sheet = safe_name("_".join(rel.with_suffix("").parts))[:31]
            df.to_excel(writer, sheet_name=sheet, index=False)
            index_rows.append({"sheet": sheet, "file": str(rel), "rows": len(df), "columns": ";".join(map(str, df.columns))})
        pd.DataFrame(index_rows).to_excel(writer, sheet_name="XVG_index", index=False)

# =============================================================================
# Workbook loading and comparative analysis
# =============================================================================

def available_sheets(workbook: Path) -> List[str]:
    return pd.ExcelFile(workbook).sheet_names


def load_time_series(workbook: Path, sheet: str) -> Tuple[pd.DataFrame, List[str]]:
    df = pd.read_excel(workbook, sheet_name=sheet)
    df.columns = [clean_column_name(c) for c in df.columns]
    if "Time_ns" not in df.columns:
        if "Time (ps)" in df.columns:
            df["Time_ns"] = pd.to_numeric(df["Time (ps)"], errors="coerce") / 1000.0
        elif "Time_ps" in df.columns:
            df["Time_ns"] = pd.to_numeric(df["Time_ps"], errors="coerce") / 1000.0
        elif "Frame" in df.columns:
            df["Time_ns"] = pd.to_numeric(df["Frame"], errors="coerce")
    systems = numeric_cols(df)
    for system in systems:
        df[system] = pd.to_numeric(df[system], errors="coerce")
    return df, systems


def load_rmsf(workbook: Path, sheet: str = "RMSF") -> Tuple[pd.DataFrame, List[str]]:
    df = pd.read_excel(workbook, sheet_name=sheet)
    df.columns = [clean_column_name(c) for c in df.columns]
    residue_col = "Residue#"
    if residue_col not in df.columns:
        for candidate in ("Residue_Number", "residue", "Residue"):
            if candidate in df.columns:
                df = df.rename(columns={candidate: "Residue#"})
                break
    systems = [c for c in df.columns if c != "Residue#" and not str(c).startswith("Unnamed")]
    for system in systems:
        df[system] = pd.to_numeric(df[system], errors="coerce")
    return df, systems


def load_block_sheet(workbook: Path, sheet_name: str, required_cols: List[str]) -> pd.DataFrame:
    raw = pd.read_excel(workbook, sheet_name=sheet_name, header=[0, 1])
    blocks = []
    for system in raw.columns.get_level_values(0).unique():
        system = str(system)
        if system.startswith("Unnamed"):
            continue
        block = raw[system].copy()
        block = block.dropna(axis=1, how="all")
        block.columns = [clean_column_name(c) for c in block.columns]
        block = block.loc[:, ~pd.Index(block.columns).duplicated()]
        # Normalize distance units if workbook uses A suffix.
        rename = {}
        for c in block.columns:
            if c == "Minimum_Distance_A": rename[c] = "Minimum_Distance"
            if c == "Maximum_Distance_A": rename[c] = "Maximum_Distance"
        if rename:
            block = block.rename(columns=rename)
        available = [c for c in required_cols if c in block.columns]
        if len(available) < 2:
            continue
        block = block[available].copy()
        first_required = required_cols[0]
        block = block.dropna(subset=[first_required])
        block.insert(0, "System", system)
        block.insert(1, "Ligand", system_label(system))
        blocks.append(block)
    if not blocks:
        raise ValueError(f"No valid blocks found in sheet {sheet_name!r}.")
    return pd.concat(blocks, ignore_index=True)


def load_catalytic_state(workbook: Path, distances: List[str]) -> pd.DataFrame:
    required = ["Frame", "Time (ps)"] + distances
    # If columns do not exactly exist, try all default and available subset.
    df = load_block_sheet(workbook, "catalytic_state", required)
    if "Time_ns" not in df.columns:
        df["Time_ns"] = pd.to_numeric(df.get("Time (ps)"), errors="coerce") / 1000.0
    for col in distances:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def load_contacts(workbook: Path, sheet_name: str = "P_L_contact") -> pd.DataFrame:
    df = load_block_sheet(workbook, sheet_name, DEFAULT_CONTACT_COLUMNS)
    for col in ["Residue_Number", "Occupancy_Percent", "Minimum_Distance", "Maximum_Distance"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    df["Residue_Number"] = df["Residue_Number"].astype("Int64")
    df["Residue_ID"] = df["Residue_Name"].astype(str) + df["Residue_Number"].astype(str)
    df["Distance_Range"] = df["Maximum_Distance"] - df["Minimum_Distance"]
    df["Distance_Midpoint"] = (df["Maximum_Distance"] + df["Minimum_Distance"]) / 2.0
    return df.dropna(subset=["Residue_Number", "Occupancy_Percent"])


def summarize_time_series(df: pd.DataFrame, systems: List[str], metric: str, analysis_start_ns: float) -> pd.DataFrame:
    rows = []
    work = df.copy()
    if "Time_ns" in work.columns:
        work = work.loc[pd.to_numeric(work["Time_ns"], errors="coerce") >= analysis_start_ns]
    for system in systems:
        if system not in work.columns:
            continue
        s = stats_dict(work[system])
        if s["N"] == 0:
            continue
        rows.append({"Metric": metric, "System": system, "Ligand": system_label(system), "Analysis_Start_ns": analysis_start_ns, **s})
    return pd.DataFrame(rows)


def summarize_rmsf(df: pd.DataFrame, systems: List[str]) -> pd.DataFrame:
    rows = []
    for system in systems:
        values = pd.to_numeric(df[system], errors="coerce").dropna()
        if values.empty:
            continue
        max_idx = values.idxmax()
        rows.append({
            "System": system,
            "Ligand": system_label(system),
            "Mean_RMSF": float(values.mean()),
            "SD_RMSF": float(values.std(ddof=1)) if len(values) > 1 else 0.0,
            "Median_RMSF": float(values.median()),
            "Max_RMSF": float(values.max()),
            "Residue_at_Max_RMSF": df.loc[max_idx, "Residue#"] if "Residue#" in df.columns else max_idx,
        })
    return pd.DataFrame(rows)


def annotate_ceh_catalytic_states(
    df: pd.DataFrame,
    distances: List[str],
    hbond_cutoff: float,
    attack_cutoff: float,
    oxyanion_cutoff: float,
) -> pd.DataFrame:
    out = df.copy()
    cutoffs = {}
    for col in distances:
        lc = col.lower()
        if "ser" in lc and "his" in lc:
            cutoffs[col] = hbond_cutoff
        elif "his" in lc and ("glu" in lc or "acid" in lc):
            cutoffs[col] = hbond_cutoff
        elif "carbonyl" in lc or "c26" in lc or "attack" in lc:
            cutoffs[col] = attack_cutoff
        elif "gly" in lc or "oxyanion" in lc or "o2" in lc:
            cutoffs[col] = oxyanion_cutoff
        else:
            cutoffs[col] = hbond_cutoff
    ok_cols = []
    for col, cutoff in cutoffs.items():
        if col not in out.columns:
            continue
        ok = f"{col}_ok"
        out[ok] = pd.to_numeric(out[col], errors="coerce") <= cutoff
        ok_cols.append(ok)

    serhis_cols = [c for c in ok_cols if "ser" in c.lower() and "his" in c.lower()]
    hisacid_cols = [c for c in ok_cols if "his" in c.lower() and ("glu" in c.lower() or "acid" in c.lower())]
    attack_cols = [c for c in ok_cols if "carbonyl" in c.lower() or "c26" in c.lower() or "attack" in c.lower()]
    oxy_cols = [c for c in ok_cols if "gly" in c.lower() or "oxyanion" in c.lower() or "o2" in c.lower()]

    out["Triad_contact_ok"] = True
    if serhis_cols:
        out["Triad_contact_ok"] &= out[serhis_cols].any(axis=1)
    if hisacid_cols:
        out["Triad_contact_ok"] &= out[hisacid_cols].any(axis=1)
    out["Attack_distance_ok"] = out[attack_cols].any(axis=1) if attack_cols else False
    out["Oxyanion_any_ok"] = out[oxy_cols].any(axis=1) if oxy_cols else False
    out["Oxyanion_both_ok"] = out[oxy_cols].all(axis=1) if len(oxy_cols) >= 2 else out["Oxyanion_any_ok"]
    out["Partial_catalytic_ok"] = out["Triad_contact_ok"] & out["Attack_distance_ok"] & out["Oxyanion_any_ok"]
    out["Fully_productive_ok"] = out["Triad_contact_ok"] & out["Attack_distance_ok"] & out["Oxyanion_both_ok"]
    out["Catalytic_score_0_to_5"] = out[ok_cols].sum(axis=1) if ok_cols else 0
    out["Catalytic_score_fraction"] = out["Catalytic_score_0_to_5"] / max(len(ok_cols), 1)
    return out


def summarize_catalytic(df: pd.DataFrame, distances: List[str], analysis_start_ns: float) -> Tuple[pd.DataFrame, pd.DataFrame]:
    work = df.copy()
    if "Time_ns" in work.columns:
        work = work.loc[pd.to_numeric(work["Time_ns"], errors="coerce") >= analysis_start_ns]
    distance_rows = []
    for (system, ligand), group in work.groupby(["System", "Ligand"]):
        for metric in distances:
            if metric not in group.columns:
                continue
            s = stats_dict(group[metric])
            if s["N"] == 0:
                continue
            renamed = {f"{k}_nm" if k in ["Mean", "SD", "SEM", "Median", "Q1", "Q3", "Min", "Max"] else k: v for k, v in s.items()}
            distance_rows.append({"System": system, "Ligand": ligand, "Distance": metric, "Analysis_Start_ns": analysis_start_ns, **renamed})

    criteria = [c for c in work.columns if c.endswith("_ok")]
    criteria += ["Triad_contact_ok", "Attack_distance_ok", "Oxyanion_any_ok", "Oxyanion_both_ok", "Partial_catalytic_ok", "Fully_productive_ok"]
    criteria = list(dict.fromkeys([c for c in criteria if c in work.columns]))

    criteria_rows = []
    for (system, ligand), group in work.groupby(["System", "Ligand"]):
        row = {
            "System": system,
            "Ligand": ligand,
            "Analysis_Start_ns": analysis_start_ns,
            "N_frames": int(len(group)),
            "Mean_catalytic_score_0_to_5": float(group["Catalytic_score_0_to_5"].mean()) if "Catalytic_score_0_to_5" in group else np.nan,
            "Mean_catalytic_score_fraction": float(group["Catalytic_score_fraction"].mean()) if "Catalytic_score_fraction" in group else np.nan,
        }
        for criterion in criteria:
            row[f"{criterion}_percent"] = 100.0 * pd.to_numeric(group[criterion], errors="coerce").mean()
        criteria_rows.append(row)
    return pd.DataFrame(distance_rows), pd.DataFrame(criteria_rows)


def summarize_contacts(contacts: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    contact_summary = (
        contacts.groupby(["Residue_ID", "Residue_Name", "Residue_Number"], dropna=False)
        .agg(
            Systems_present=("System", "nunique"),
            Mean_occupancy_percent=("Occupancy_Percent", "mean"),
            Max_occupancy_percent=("Occupancy_Percent", "max"),
            Min_distance=("Minimum_Distance", "min"),
            Max_distance=("Maximum_Distance", "max"),
        )
        .reset_index()
        .sort_values(["Systems_present", "Mean_occupancy_percent"], ascending=[False, False])
    )
    top_contacts = (
        contacts.sort_values(["System", "Occupancy_Percent"], ascending=[True, False])
        .groupby("System", group_keys=False)
        .head(20)
        .copy()
    )
    occupancy_matrix = contacts.pivot_table(index="Residue_ID", columns="Ligand", values="Occupancy_Percent", aggfunc="max").fillna(0.0)
    distance_matrix = contacts.pivot_table(index="Residue_ID", columns="Ligand", values="Maximum_Distance", aggfunc="max")
    return contact_summary, top_contacts, occupancy_matrix, distance_matrix

# =============================================================================
# Plotting
# =============================================================================

def plot_timeseries(df: pd.DataFrame, systems: List[str], y_label: str, title: str, figdir: Path,
                    stem: str, dpi: int, formats: Iterable[str], color_map: Dict[str, tuple]) -> None:
    if "Time_ns" not in df.columns:
        return
    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    for system in systems:
        ax.plot(df["Time_ns"], df[system], label=system_label(system), linewidth=1.2, color=get_system_color(system, color_map))
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel(y_label)
    ax.set_title(title)
    ax.legend(frameon=False, ncols=2)
    ax.grid(True, linewidth=0.3, alpha=0.5)
    save_figure(fig, figdir, stem, dpi, formats)


def plot_rmsf(df: pd.DataFrame, systems: List[str], figdir: Path, dpi: int, formats: Iterable[str], color_map: Dict[str, tuple]) -> None:
    if "Residue#" not in df.columns:
        return
    fig, ax = plt.subplots(figsize=(9.5, 4.8))
    for system in systems:
        ax.plot(df["Residue#"], df[system], label=system_label(system), linewidth=1.0, color=get_system_color(system, color_map))
    ax.set_xlabel("Residue number")
    ax.set_ylabel("RMSF (nm)")
    ax.set_title("Protein RMSF profile")
    ax.legend(frameon=False, ncols=2)
    ax.grid(True, linewidth=0.3, alpha=0.5)
    save_figure(fig, figdir, "06_RMSF_profile", dpi, formats)


def plot_metric_boxplots(metric_frames: Dict[str, Tuple[pd.DataFrame, List[str]]], figdir: Path, dpi: int,
                         formats: Iterable[str], color_map: Dict[str, tuple]) -> None:
    rows = []
    for metric, (df, systems) in metric_frames.items():
        for system in systems:
            if system not in df:
                continue
            for value in pd.to_numeric(df[system], errors="coerce").dropna():
                rows.append({"Metric": metric, "System": system, "Ligand": system_label(system), "Value": value})
    data = pd.DataFrame(rows)
    if data.empty:
        return
    metrics = list(data["Metric"].drop_duplicates())
    fig, axes = plt.subplots(1, len(metrics), figsize=(5 * len(metrics), 4.8), squeeze=False)
    for ax, metric in zip(axes.flat, metrics):
        sub = data.loc[data["Metric"] == metric]
        labels = list(sub["Ligand"].drop_duplicates())
        values = [sub.loc[sub["Ligand"] == label, "Value"].to_numpy() for label in labels]
        box = ax.boxplot(values, tick_labels=labels, showfliers=False, patch_artist=True)
        ligand_to_system = dict(zip(sub["Ligand"], sub["System"]))
        for patch, label in zip(box["boxes"], labels):
            patch.set_facecolor(get_system_color(ligand_to_system.get(label, label), color_map))
            patch.set_alpha(0.6)
        ax.set_title(metric)
        ax.set_ylabel(metric)
        ax.tick_params(axis="x", rotation=35)
        ax.grid(True, axis="y", linewidth=0.3, alpha=0.5)
    fig.suptitle("Distribution of global trajectory metrics", y=1.02)
    save_figure(fig, figdir, "07_global_metric_boxplots", dpi, formats)


def plot_heatmap(matrix: pd.DataFrame, title: str, cbar_label: str, figdir: Path, stem: str, dpi: int, formats: Iterable[str],
                 vmin=None, vmax=None) -> None:
    if matrix is None or matrix.empty:
        return
    matrix = matrix.copy()
    fig_height = max(4.5, min(16, 0.24 * matrix.shape[0]))
    fig_width = max(7.2, min(16, 0.60 * matrix.shape[1] + 4))
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    image = ax.imshow(matrix.to_numpy(dtype=float), aspect="auto", vmin=vmin, vmax=vmax)
    ax.set_xticks(np.arange(matrix.shape[1]))
    ax.set_yticks(np.arange(matrix.shape[0]))
    ax.set_xticklabels(matrix.columns, rotation=45, ha="right")
    ax.set_yticklabels(matrix.index, fontsize=6 if matrix.shape[0] > 35 else 8)
    ax.set_title(title)
    cbar = fig.colorbar(image, ax=ax, shrink=0.75)
    cbar.set_label(cbar_label)
    save_figure(fig, figdir, stem, dpi, formats)


def plot_contacts(contacts: pd.DataFrame, occupancy_matrix: pd.DataFrame, distance_matrix: pd.DataFrame,
                  figdir: Path, top_n: int, dpi: int, formats: Iterable[str], color_map: Dict[str, tuple]) -> None:
    if occupancy_matrix.empty:
        return
    sorted_occ = occupancy_matrix.loc[occupancy_matrix.mean(axis=1).sort_values(ascending=False).index]
    plot_heatmap(sorted_occ, "Protein-ligand contact occupancy", "Occupancy (%)", figdir, "10_contact_occupancy_heatmap", dpi, formats, vmin=0, vmax=100)
    if distance_matrix is not None and not distance_matrix.empty:
        plot_heatmap(distance_matrix.reindex(sorted_occ.index), "Maximum contact distance by residue", "Maximum distance", figdir, "11_contact_max_distance_heatmap", dpi, formats)
    systems = list(contacts["System"].drop_duplicates())
    if not systems:
        return
    ncols = 2
    nrows = math.ceil(len(systems) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(12, max(4.5, 4.2 * nrows)), squeeze=False)
    for ax, system in zip(axes.flat, systems):
        top = contacts.loc[contacts["System"] == system].sort_values("Occupancy_Percent", ascending=False).head(top_n)
        top = top.sort_values("Occupancy_Percent", ascending=True)
        ax.barh(top["Residue_ID"], top["Occupancy_Percent"], color=get_system_color(system, color_map))
        ax.set_title(f"Top contacts: {system_label(system)}")
        ax.set_xlabel("Occupancy (%)")
        ax.set_xlim(0, max(100, top["Occupancy_Percent"].max() * 1.05 if len(top) else 100))
        ax.grid(True, axis="x", linewidth=0.3, alpha=0.5)
    for ax in axes.flat[len(systems):]:
        ax.axis("off")
    fig.suptitle("Top protein-ligand contacts per system", y=1.01)
    save_figure(fig, figdir, "12_top_contacts_by_system", dpi, formats)


def plot_catalytic_distances(df: pd.DataFrame, distances: List[str], figdir: Path, dpi: int,
                             formats: Iterable[str], color_map: Dict[str, tuple]) -> None:
    distances = [d for d in distances if d in df.columns]
    if not distances or "Time_ns" not in df.columns:
        return
    n = len(distances)
    fig, axes = plt.subplots(n, 1, figsize=(9.2, 2.15 * n), sharex=True)
    if n == 1:
        axes = [axes]
    for ax, distance in zip(axes, distances):
        for system, group in df.groupby("System"):
            ax.plot(group["Time_ns"], group[distance], label=system_label(system), linewidth=1.0, color=get_system_color(system, color_map))
        ax.set_ylabel("nm")
        ax.set_title(distance)
        ax.grid(True, linewidth=0.3, alpha=0.5)
    axes[-1].set_xlabel("Time (ns)")
    axes[0].legend(frameon=False, ncols=4, loc="upper center", bbox_to_anchor=(0.5, 1.55))
    fig.suptitle("CEH catalytic geometry distance profiles", y=1.02)
    save_figure(fig, figdir, "08_CEH_catalytic_distance_profiles", dpi, formats)


def plot_catalytic_criteria(criteria_summary: pd.DataFrame, figdir: Path, dpi: int,
                            formats: Iterable[str], color_map: Dict[str, tuple]) -> None:
    if criteria_summary is None or criteria_summary.empty:
        return
    crit_cols = [c for c in criteria_summary.columns if c.endswith("_percent")]
    if not crit_cols:
        return
    matrix = criteria_summary.set_index("Ligand")[crit_cols]
    matrix.columns = [c.replace("_percent", "").replace("_ok", "").replace("_", " ") for c in matrix.columns]
    plot_heatmap(matrix, "CEH catalytic criterion occupancy", "Frames satisfying criterion (%)", figdir, "09_CEH_catalytic_criteria_heatmap", dpi, formats, vmin=0, vmax=100)
    score_col = "Mean_catalytic_score_0_to_5" if "Mean_catalytic_score_0_to_5" in criteria_summary else "Mean_catalytic_score_fraction"
    score = criteria_summary.sort_values(score_col, ascending=True)
    fig, ax = plt.subplots(figsize=(7, 4.2))
    ax.barh(score["Ligand"], score[score_col], color=[get_system_color(s, color_map) for s in score["System"]])
    ax.set_xlabel(score_col.replace("_", " "))
    ax.set_title("Relative CEH catalytic organization score")
    ax.grid(True, axis="x", linewidth=0.3, alpha=0.5)
    save_figure(fig, figdir, "09b_CEH_catalytic_score_bar", dpi, formats)


def plot_dssp_like(df: pd.DataFrame, systems: List[str], figdir: Path, dpi: int, formats: Iterable[str], color_map: Dict[str, tuple]) -> None:
    # Supports a single DSSP_count sheet where columns are secondary-structure classes, or multi-system DSSP sheets.
    if df.empty or "Time_ns" not in df.columns:
        return
    plot_cols = [c for c in df.columns if c not in INDEX_COLUMNS]
    if not plot_cols:
        return
    fig, ax = plt.subplots(figsize=(10, 5))
    for col in plot_cols:
        ax.plot(df["Time_ns"], df[col], linewidth=1.2, label=str(col))
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Secondary structure count")
    ax.set_title("DSSP secondary-structure content over time")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    save_figure(fig, figdir, "13_DSSP_secondary_structure", dpi, formats)

# =============================================================================
# Direct MDAnalysis contact analysis
# =============================================================================

def residue_contact_analysis(topology: Path, trajectory: Path, outdir: Path, ligand_resname: str = "UNL",
                             cutoff_A: float = 4.5, stride: int = 20, protein_selection: str = "protein") -> pd.DataFrame:
    if not HAS_MDA:
        raise SystemExit("MDAnalysis is required for direct trajectory contact analysis. Install with: conda install -c conda-forge mdanalysis")
    ensure_dir(outdir)
    u = mda.Universe(str(topology), str(trajectory))
    ligand = u.select_atoms(f"resname {ligand_resname}")
    if len(ligand) == 0:
        raise SystemExit(f"No ligand atoms found with resname {ligand_resname!r}")
    protein_residues = u.select_atoms(protein_selection).residues
    contact_frames: Dict[Tuple[int, str], int] = {}
    min_distances: Dict[Tuple[int, str], float] = {}
    max_distances: Dict[Tuple[int, str], float] = {}
    nframes = 0
    for ts in u.trajectory[::max(1, stride)]:
        nframes += 1
        for res in protein_residues:
            dmat = distance_array(res.atoms.positions, ligand.positions)
            dmin = float(np.min(dmat))
            key = (int(res.resid), str(res.resname))
            if key not in contact_frames:
                contact_frames[key] = 0
                min_distances[key] = dmin
                max_distances[key] = dmin
            min_distances[key] = min(min_distances[key], dmin)
            max_distances[key] = max(max_distances[key], dmin)
            if dmin <= cutoff_A:
                contact_frames[key] += 1
    rows = []
    for (resid, resname), count in contact_frames.items():
        occupancy = 100.0 * count / max(nframes, 1)
        if occupancy > 0:
            rows.append({
                "Residue_Number": resid,
                "Residue_Name": resname,
                "Residue_ID": f"{resname}{resid}",
                "Occupancy_Percent": round(occupancy, 2),
                "Minimum_Distance_A": round(min_distances[(resid, resname)], 3),
                "Maximum_Distance_A": round(max_distances[(resid, resname)], 3),
                "Frames_Analyzed": nframes,
                "Cutoff_A": cutoff_A,
            })
    df = pd.DataFrame(rows).sort_values("Occupancy_Percent", ascending=False)
    df.to_csv(outdir / "protein_ligand_residue_contacts.csv", index=False)
    # Plot top contacts
    if not df.empty:
        top = df.head(25).sort_values("Occupancy_Percent", ascending=True)
        fig, ax = plt.subplots(figsize=(8, max(4, 0.28 * len(top))))
        ax.barh(top["Residue_ID"], top["Occupancy_Percent"])
        ax.set_xlabel("Contact occupancy (%)")
        ax.set_title("Top protein-ligand residue contacts")
        ax.grid(True, axis="x", linewidth=0.3, alpha=0.5)
        save_figure(fig, outdir, "protein_ligand_top_contacts", 300, ["png"])
    return df

# =============================================================================
# GROMACS analysis script generation
# =============================================================================

def write_gromacs_analysis_script(
    out_script: Path,
    mode: str,
    tpr: str,
    xtc: str,
    index: str,
    ligand_group: str,
    ligand_resname: str,
    pocket_residues: Optional[List[int]] = None,
    key_residues: Optional[List[int]] = None,
    carbonyl_atom: str = "C26",
    oxyanion_oxygen: str = "O2",
) -> None:
    mode = mode.lower()
    pocket_residues = pocket_residues or (CEH_POCKET_RESIDUES if mode == "ceh" else [])
    key_residues = key_residues or (CEH_KEY_CONTACT_RESIDUES if mode == "ceh" else [])
    pocket_expr = " | ".join(f"r {r}" for r in pocket_residues)
    key_expr = " | ".join(f"r {r}" for r in key_residues)
    lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        "mkdir -p analysis_xvg frames",
        f"TPR=\"{tpr}\"",
        f"XTC=\"{xtc}\"",
        f"INDEX=\"{index}\"",
        "",
        "echo 'Step 1: center and fit trajectory'",
        "printf 'Protein\\nSystem\\n' | gmx trjconv -s \"$TPR\" -f \"$XTC\" -o centered.xtc -pbc mol -center",
        "printf 'Protein\\nSystem\\n' | gmx trjconv -s \"$TPR\" -f centered.xtc -o md_fit.xtc -fit rot+trans",
        "",
        "echo 'Step 2: global stability metrics'",
        "printf 'Backbone\\nBackbone\\n' | gmx rms -s \"$TPR\" -f md_fit.xtc -o analysis_xvg/rmsd_backbone.xvg",
        f"printf 'Protein\\n{ligand_group}\\n' | gmx rms -s \"$TPR\" -f md_fit.xtc -n \"$INDEX\" -o analysis_xvg/ligand_rmsd.xvg || true",
        "printf 'Protein\\n' | gmx gyrate -s \"$TPR\" -f md_fit.xtc -o analysis_xvg/gyrate.xvg",
        "printf 'C-alpha\\n' | gmx rmsf -s \"$TPR\" -f md_fit.xtc -o analysis_xvg/rmsf_ca.xvg -res",
        "printf 'C-alpha\\n' | gmx rmsf -s \"$TPR\" -f md_fit.xtc -res -oq analysis_xvg/rmsf_bfactor.pdb || true",
        "",
        "echo 'Step 3: surface and secondary-structure metrics'",
        "gmx sasa -s \"$TPR\" -f md_fit.xtc -o analysis_xvg/sasa_total.xvg -surface 'Protein' -output 'Protein' || true",
        "gmx sasa -s \"$TPR\" -f md_fit.xtc -or analysis_xvg/residue_sasa.xvg -surface 'Protein' -output 'Protein' || true",
        "gmx dssp -s \"$TPR\" -f centered.xtc -o analysis_xvg/dssp.dat -num analysis_xvg/dssp_count.xvg -sel 'Protein' -hmode dssp -clear || true",
        "",
        "echo 'Step 4: protein-ligand interaction metrics'",
        f"printf 'Protein\\n{ligand_group}\\n' | gmx mindist -s \"$TPR\" -f md_fit.xtc -n \"$INDEX\" -od analysis_xvg/contacts.xvg -on analysis_xvg/contacts_num.xvg || true",
        "",
    ]
    if pocket_expr:
        lines += [
            "echo 'Creating pocket index'",
            "cat > pocket_commands.ndx <<'EOF'",
            pocket_expr,
            "name 21 Pocket",
            "q",
            "EOF",
            "gmx make_ndx -f \"$TPR\" -o pocket.ndx < pocket_commands.ndx",
            f"gmx distance -s \"$TPR\" -f md_fit.xtc -n pocket.ndx -select 'com of group \"{ligand_group}\" plus com of group \"Pocket\"' -oall analysis_xvg/ligand_pocket_com_distance.xvg || true",
            f"printf 'Pocket\\n{ligand_group}\\n' | gmx mindist -s \"$TPR\" -f md_fit.xtc -n pocket.ndx -respertime -or analysis_xvg/residue_mindist.xvg || true",
            "",
        ]
    if mode == "ceh":
        lines += [
            "echo 'Step 5: CEH catalytic-state distances'",
            "cat > catalytic_commands.ndx <<'EOF'",
            "a OG & r 207",
            "name 21 SER207_OG",
            "a NE2 & r 508",
            "name 22 HIS508_NE2",
            "a ND1 & r 508",
            "name 23 HIS508_ND1",
            "a OE2 & r 388",
            "name 24 GLU388_OE2",
            "a OE1 & r 388",
            "name 25 GLU388_OE1",
            "a OE1 & r 206",
            "name 26 GLU206_OE1",
            "a OE2 & r 206",
            "name 27 GLU206_OE2",
            f"a {carbonyl_atom} & r {ligand_resname}",
            "name 28 carbonylC",
            "r 128 & a N",
            "name 29 GLY128_N",
            "r 129 & a N",
            "name 30 GLY129_N",
            f"r {ligand_resname} & a {oxyanion_oxygen}",
            "name 31 ligand_O2",
            key_expr,
            "name 32 Key_contact",
            "q",
            "EOF",
            "gmx make_ndx -f \"$TPR\" -o catalytic.ndx < catalytic_commands.ndx",
            "gmx distance -s \"$TPR\" -f md_fit.xtc -n catalytic.ndx -oall analysis_xvg/ser207_his508.xvg -select 'group \"SER207_OG\" plus group \"HIS508_NE2\"' || true",
            "gmx distance -s \"$TPR\" -f md_fit.xtc -n catalytic.ndx -oall analysis_xvg/his508_glu388_oe2.xvg -select 'group \"HIS508_ND1\" plus group \"GLU388_OE2\"' || true",
            "gmx distance -s \"$TPR\" -f md_fit.xtc -n catalytic.ndx -oall analysis_xvg/his508_glu206_oe2.xvg -select 'group \"HIS508_ND1\" plus group \"GLU206_OE2\"' || true",
            "gmx distance -s \"$TPR\" -f md_fit.xtc -n catalytic.ndx -oall analysis_xvg/ser207_carbonylC_distance.xvg -select 'group \"SER207_OG\" plus group \"carbonylC\"' || true",
            "gmx distance -s \"$TPR\" -f md_fit.xtc -n catalytic.ndx -oall analysis_xvg/gly128N_O2_distance.xvg -select 'group \"GLY128_N\" plus group \"ligand_O2\"' || true",
            "gmx distance -s \"$TPR\" -f md_fit.xtc -n catalytic.ndx -oall analysis_xvg/gly129N_O2_distance.xvg -select 'group \"GLY129_N\" plus group \"ligand_O2\"' || true",
            f"gmx distance -s \"$TPR\" -f md_fit.xtc -n catalytic.ndx -select 'com of group \"{ligand_group}\" plus com of group \"Key_contact\"' -oall analysis_xvg/catalytic_contact_distance.xvg || true",
            f"printf 'Key_contact\\n{ligand_group}\\n' | gmx mindist -s \"$TPR\" -f md_fit.xtc -n catalytic.ndx -respertime -or analysis_xvg/catalytic_contact_mindist.xvg || true",
            "",
            "echo 'Step 6: representative protein-ligand frame extraction'",
            "cat > analysis_commands.ndx <<'EOF'",
            f"Protein | {ligand_group}",
            "name 33 Protein_Ligand",
            "q",
            "EOF",
            "gmx make_ndx -f \"$TPR\" -o analysis.ndx < analysis_commands.ndx || true",
            "printf 'Protein\\nProtein_Ligand\\n' | gmx trjconv -f \"$XTC\" -s \"$TPR\" -n analysis.ndx -pbc mol -center -ur compact -o fixed.xtc || true",
            "printf 'Protein_Ligand\\n' | gmx trjconv -f fixed.xtc -s \"$TPR\" -n analysis.ndx -sep -dt 50 -o frames/frame_.pdb || true",
        ]
    lines += ["echo 'Post-MD GROMACS analysis script complete.'", ""]
    out_script.write_text("\n".join(lines))
    os.chmod(out_script, 0o755)

# =============================================================================
# Workbook command
# =============================================================================

def run_workbook_analysis(args: argparse.Namespace) -> None:
    workbook = Path(args.input)
    outdir = ensure_dir(Path(args.outdir))
    figdir = ensure_dir(outdir / "figures")
    formats = [f.strip() for f in args.formats.split(",") if f.strip()]

    plt.rcParams.update({"font.size": 9, "axes.titlesize": 11, "axes.labelsize": 10, "legend.fontsize": 8, "figure.dpi": 120})

    sheets = available_sheets(workbook)
    metric_frames: Dict[str, Tuple[pd.DataFrame, List[str]]] = {}
    all_systems: List[str] = []
    summary_tables: Dict[str, pd.DataFrame] = {}

    for sheet, meta in PLOT_METRIC_SHEETS.items():
        if sheet in sheets:
            df, systems = load_time_series(workbook, sheet)
            metric_frames[sheet] = (df, systems)
            all_systems += systems

    if "RMSF" in sheets:
        rmsf_df, rmsf_systems = load_rmsf(workbook)
        all_systems += rmsf_systems
    else:
        rmsf_df, rmsf_systems = pd.DataFrame(), []

    catalytic = None
    catalytic_distances = None
    catalytic_criteria = None
    if args.mode.lower() == "ceh" and "catalytic_state" in sheets:
        requested = [d.strip() for d in args.catalytic_distances.split(",") if d.strip()]
        catalytic = load_catalytic_state(workbook, requested)
        available_distances = [d for d in requested if d in catalytic.columns]
        catalytic = annotate_ceh_catalytic_states(catalytic, available_distances, args.hbond_cutoff, args.attack_cutoff, args.oxyanion_cutoff)
        catalytic_distances, catalytic_criteria = summarize_catalytic(catalytic, available_distances, args.analysis_start_ns)
        all_systems += list(catalytic["System"].drop_duplicates())

    contacts = None
    contact_summary = top_contacts = occupancy_matrix = distance_matrix = pd.DataFrame()
    if "P_L_contact" in sheets:
        contacts = load_contacts(workbook)
        contact_summary, top_contacts, occupancy_matrix, distance_matrix = summarize_contacts(contacts)
        all_systems += list(contacts["System"].drop_duplicates())

    color_map = build_color_map(all_systems)

    metric_summary_parts = []
    for sheet, (df, systems) in metric_frames.items():
        meta = PLOT_METRIC_SHEETS.get(sheet, {"ylabel": sheet, "title": sheet, "stem": safe_name(sheet)})
        plot_timeseries(df, systems, meta["ylabel"], meta["title"], figdir, meta["stem"], args.dpi, formats, color_map)
        metric_summary_parts.append(summarize_time_series(df, systems, sheet, args.analysis_start_ns))

    if rmsf_systems:
        plot_rmsf(rmsf_df, rmsf_systems, figdir, args.dpi, formats, color_map)
        summary_tables["RMSF_summary"] = summarize_rmsf(rmsf_df, rmsf_systems)

    if metric_frames:
        # Keep core metrics for boxplots to avoid overly wide panels.
        box_metrics = {k: v for k, v in metric_frames.items() if k in ["C_alpha_RMSD", "Lig_RMSD", "rGyr", "SASA"]}
        plot_metric_boxplots(box_metrics, figdir, args.dpi, formats, color_map)

    if catalytic is not None:
        plot_catalytic_distances(catalytic, [d for d in args.catalytic_distances.split(",") if d.strip()], figdir, args.dpi, formats, color_map)
        plot_catalytic_criteria(catalytic_criteria, figdir, args.dpi, formats, color_map)
        summary_tables["CEH_catalytic_distances"] = catalytic_distances
        summary_tables["CEH_catalytic_criteria"] = catalytic_criteria
        summary_tables["CEH_catalytic_frame_states"] = catalytic

    if contacts is not None:
        plot_contacts(contacts, occupancy_matrix, distance_matrix, figdir, args.top_n_contacts, args.dpi, formats, color_map)
        summary_tables["Contacts_long"] = contacts
        summary_tables["Contact_summary"] = contact_summary
        summary_tables["Top_contacts"] = top_contacts
        summary_tables["Contact_occupancy_matrix"] = occupancy_matrix.reset_index()
        summary_tables["Contact_distance_matrix"] = distance_matrix.reset_index()

    if "DSSP_count" in sheets:
        dssp, dssp_cols = load_time_series(workbook, "DSSP_count")
        plot_dssp_like(dssp, dssp_cols, figdir, args.dpi, formats, color_map)
        summary_tables["DSSP_count"] = dssp

    metric_summary = pd.concat(metric_summary_parts, ignore_index=True) if metric_summary_parts else pd.DataFrame()
    if not metric_summary.empty:
        summary_tables["Trajectory_metrics"] = metric_summary

    color_key = pd.DataFrame([{"System": s, "Ligand": system_label(s), "Color_RGB": str(color_map[s])} for s in color_map])
    summary_tables["Color_key"] = color_key

    parameter_table = pd.DataFrame([
        {"Parameter": "mode", "Value": args.mode, "Meaning": "general or ceh analysis mode"},
        {"Parameter": "analysis_start_ns", "Value": args.analysis_start_ns, "Meaning": "Frames before this time are excluded from summary statistics"},
        {"Parameter": "hbond_cutoff_nm", "Value": args.hbond_cutoff, "Meaning": "CEH Ser-His and His-acid distance cutoff"},
        {"Parameter": "attack_cutoff_nm", "Value": args.attack_cutoff, "Meaning": "CEH Ser207-to-ligand carbonyl distance cutoff"},
        {"Parameter": "oxyanion_cutoff_nm", "Value": args.oxyanion_cutoff, "Meaning": "CEH Gly128/Gly129-to-ligand oxygen distance cutoff"},
    ])
    summary_tables = {"ReadMe_parameters": parameter_table, **summary_tables}

    # Comparative index
    comparative = pd.DataFrame()
    if not metric_summary.empty:
        comparative = metric_summary.pivot_table(index=["System", "Ligand"], columns="Metric", values="Mean", aggfunc="first").reset_index()
        if catalytic_criteria is not None and not catalytic_criteria.empty:
            keep = ["System", "Mean_catalytic_score_0_to_5", "Mean_catalytic_score_fraction", "Fully_productive_ok_percent", "Partial_catalytic_ok_percent"]
            keep = [c for c in keep if c in catalytic_criteria.columns]
            comparative = comparative.merge(catalytic_criteria[keep], on="System", how="left")
        summary_tables["Comparative_index"] = comparative

    out_xlsx = outdir / "post_MD_complex_analysis_summary.xlsx"
    with pd.ExcelWriter(out_xlsx, engine="xlsxwriter") as writer:
        for sheet_name, table in summary_tables.items():
            if table is None or not isinstance(table, pd.DataFrame) or table.empty:
                continue
            safe_sheet = safe_name(sheet_name)[:31]
            table.to_excel(writer, sheet_name=safe_sheet, index=False)
            ws = writer.sheets[safe_sheet]
            ws.freeze_panes(1, 0)
            for i, col in enumerate(table.columns):
                try:
                    width = min(max(12, len(str(col)) + 2, int(table[col].astype(str).str.len().quantile(0.95)) + 2), 42)
                except Exception:
                    width = min(max(12, len(str(col)) + 2), 42)
                ws.set_column(i, i, width)

    # Markdown report
    report = outdir / "post_MD_complex_analysis_report.md"
    lines = [
        "# Post-MD receptor-ligand complex analysis report",
        "",
        f"Input workbook: `{workbook}`",
        f"Mode: `{args.mode}`",
        f"Summary workbook: `{out_xlsx.name}`",
        f"Figures directory: `figures/`",
        "",
        "## Outputs",
        "",
        "- `post_MD_complex_analysis_summary.xlsx`: summary tables and long-form analysis outputs.",
        "- `figures/`: comparative RMSD, ligand RMSD, Rg, SASA, RMSF, contact, DSSP, and CEH catalytic plots when available.",
    ]
    if comparative is not None and not comparative.empty:
        lines += ["", "## Comparative index preview", "", comparative.head(20).to_markdown(index=False)]
    if catalytic_criteria is not None and not catalytic_criteria.empty:
        lines += ["", "## CEH catalytic criteria preview", "", catalytic_criteria.head(20).to_markdown(index=False)]
    report.write_text("\n".join(lines))

    print(f"Analysis complete. Summary workbook: {out_xlsx}")
    print(f"Report: {report}")
    print(f"Figures: {figdir}")

# =============================================================================
# CLI
# =============================================================================

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Integrated post-MD analysis pipeline for receptor-ligand complexes.")
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("workbook", help="Analyze a multi-sheet Excel workbook of MD metrics.")
    p.add_argument("--input", required=True, help="Input trajectory workbook (.xlsx).")
    p.add_argument("--outdir", default="post_MD_analysis_output")
    p.add_argument("--mode", choices=["general", "ceh"], default="general")
    p.add_argument("--analysis-start-ns", type=float, default=0.0)
    p.add_argument("--hbond-cutoff", type=float, default=0.35, help="CEH H-bond/triad cutoff in nm.")
    p.add_argument("--attack-cutoff", type=float, default=0.45, help="CEH Ser-to-carbonyl cutoff in nm.")
    p.add_argument("--oxyanion-cutoff", type=float, default=0.35, help="CEH oxyanion cutoff in nm.")
    p.add_argument("--catalytic-distances", default=",".join(DEFAULT_CEH_CATALYTIC_DISTANCES))
    p.add_argument("--top-n-contacts", type=int, default=20)
    p.add_argument("--dpi", type=int, default=300)
    p.add_argument("--formats", default="png,svg")
    p.set_defaults(func=run_workbook_analysis)

    p = sub.add_parser("contacts", help="Run direct residue contact occupancy analysis from trajectory/topology.")
    p.add_argument("--topology", required=True)
    p.add_argument("--trajectory", required=True)
    p.add_argument("--outdir", default="contact_analysis")
    p.add_argument("--ligand-resname", default="UNL")
    p.add_argument("--cutoff-A", type=float, default=4.5)
    p.add_argument("--stride", type=int, default=20)
    p.add_argument("--protein-selection", default="protein")
    def _contacts(args):
        residue_contact_analysis(Path(args.topology), Path(args.trajectory), Path(args.outdir), args.ligand_resname, args.cutoff_A, args.stride, args.protein_selection)
    p.set_defaults(func=_contacts)

    p = sub.add_parser("xvg-collect", help="Collect GROMACS XVG files into a workbook.")
    p.add_argument("--input-dir", required=True)
    p.add_argument("--out", default="xvg_collected_workbook.xlsx")
    p.add_argument("--pattern", default="*.xvg")
    p.add_argument("--no-recursive", action="store_true")
    def _xvg(args):
        collect_xvg_to_workbook(Path(args.input_dir), Path(args.out), args.pattern, not args.no_recursive)
        print(f"Wrote: {args.out}")
    p.set_defaults(func=_xvg)

    p = sub.add_parser("make-gromacs-scripts", help="Generate a GROMACS post-MD analysis shell script.")
    p.add_argument("--out", default="run_post_md_gromacs_analysis.sh")
    p.add_argument("--mode", choices=["general", "ceh"], default="general")
    p.add_argument("--tpr", default="md.tpr")
    p.add_argument("--xtc", default="md.xtc")
    p.add_argument("--index", default="index.ndx")
    p.add_argument("--ligand-group", default="UNL")
    p.add_argument("--ligand-resname", default="UNL")
    p.add_argument("--pocket-residues", default="", help="Comma-separated residue numbers for ligand-pocket COM and residue mindist.")
    p.add_argument("--key-residues", default="", help="Comma-separated key catalytic/contact residues for CEH mode.")
    p.add_argument("--carbonyl-atom", default="C26")
    p.add_argument("--oxyanion-oxygen", default="O2")
    def _make(args):
        pocket = [int(x) for x in re.split(r"[,\s]+", args.pocket_residues.strip()) if x] if args.pocket_residues.strip() else None
        key = [int(x) for x in re.split(r"[,\s]+", args.key_residues.strip()) if x] if args.key_residues.strip() else None
        write_gromacs_analysis_script(Path(args.out), args.mode, args.tpr, args.xtc, args.index, args.ligand_group, args.ligand_resname, pocket, key, args.carbonyl_atom, args.oxyanion_oxygen)
        print(f"Wrote: {args.out}")
    p.set_defaults(func=_make)

    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
