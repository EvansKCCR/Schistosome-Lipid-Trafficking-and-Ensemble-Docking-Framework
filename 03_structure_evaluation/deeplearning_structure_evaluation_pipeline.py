#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
deeplearning_structure_evaluation_pipeline.py

Integrated evaluator for next-generation deep-learning structure predictions,
including AlphaFold3, Boltz/Boltz-2, Chai-1-like, and generic mmCIF/PDB +
confidence JSON/NPZ outputs.

Core evaluations
----------------
1. Discover and group coordinate files (.cif/.mmcif/.pdb), confidence JSONs,
   and NPZ confidence arrays by model/job name.
2. Parse JSON confidence files recursively for AF3/Boltz-like metrics:
   pLDDT, PAE, PDE, contact probabilities, pTM, ipTM, ranking/confidence score,
   clash/disorder flags, chain-pair matrices, and affinity outputs.
3. Parse mmCIF/PDB coordinates to compute residue, chain, and structure pLDDT
   from the B-factor / _atom_site.B_iso_or_equiv field.
4. Scan catalytic or functional motifs in each chain sequence and report
   motif-level pLDDT using 1-based sequence positions and residue labels.
5. Optionally convert mmCIF files to PDB using Biopython.
6. Generate TSV tables, optional Excel workbook, and publication-friendly plots.

Author: Evans Asamoah Adu workflow integration, generated with ChatGPT.
"""

from __future__ import annotations

import argparse
import csv
import glob
import hashlib
import io
import json
import math
import os
import re
import shlex
import sys
from collections import defaultdict
from pathlib import Path
from statistics import mean, median
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

# Optional packages. The pipeline still writes TSVs if these are absent.
try:
    import numpy as np  # type: ignore
except Exception:  # pragma: no cover
    np = None

try:
    import pandas as pd  # type: ignore
except Exception:  # pragma: no cover
    pd = None

try:
    import matplotlib.pyplot as plt  # type: ignore
except Exception:  # pragma: no cover
    plt = None

try:
    from Bio.PDB import MMCIFParser, PDBIO  # type: ignore
except Exception:  # pragma: no cover
    MMCIFParser = None
    PDBIO = None

AA3_TO_AA1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "SEC": "U", "PYL": "O", "ASX": "B", "GLX": "Z",
}

# Recursively searched JSON aliases. These include AF3-style and Boltz-style names
# but remain conservative to avoid treating unrelated fields as pLDDT/PAE.
JSON_METRIC_ALIASES: Dict[str, List[str]] = {
    "pae": [
        "pae", "pae_matrix", "predicted_aligned_error", "PAE",
        "PredictedAlignedError", "chain_pair_pae_min", "pair_pae_min",
    ],
    "pde": [
        "pde", "pde_matrix", "predicted_distance_error", "pair_pde",
        "predicted_distance_errors",
    ],
    "plddt": [
        "plddt", "pLDDT", "atom_plddt", "atom_plddts", "per_atom_plddt",
        "per_residue_plddt", "local_confidence", "plddts",
    ],
    "contact_prob": [
        "contact_prob", "contact_probability", "contact_probs", "predicted_contacts",
        "contact_matrix", "contact_probabilities",
    ],
    "ptm": ["ptm", "pTM", "ptm_score", "chain_ptm", "aggregate_ptm"],
    "iptm": ["iptm", "ipTM", "iptm_score", "chain_iptm", "interface_ptm", "chain_pair_iptm"],
    "ranking_score": ["ranking_score", "ranking_confidence", "confidence_score", "score"],
    "fraction_disordered": ["fraction_disordered", "disorder", "disordered_fraction"],
    "has_clash": ["has_clash", "clash", "has_clashes"],
    "num_recycles": ["num_recycles", "recycles"],
    "affinity_pred_value": ["affinity_pred_value", "affinity", "predicted_affinity"],
    "affinity_probability_binary": [
        "affinity_probability_binary", "binder_probability", "binding_probability",
        "affinity_prob", "probability_binary",
    ],
}

PAIRWISE_ALIASES: Dict[str, List[str]] = {
    "chain_pair_iptm": ["chain_pair_iptm", "pair_iptm"],
    "chain_pair_pae_min": ["chain_pair_pae_min", "pair_pae_min"],
    "chain_pair_pde": ["chain_pair_pde", "pair_pde"],
}

QUALITY_BANDS = [
    (90.0, "very_high"),
    (70.0, "confident"),
    (50.0, "low"),
    (-math.inf, "very_low"),
]


def is_finite_number(value: Any) -> bool:
    try:
        x = float(value)
    except (TypeError, ValueError):
        return False
    return math.isfinite(x)


def to_float(value: Any) -> Optional[float]:
    if is_finite_number(value):
        return float(value)
    return None


def flatten_numeric(value: Any) -> List[float]:
    out: List[float] = []

    def walk(v: Any) -> None:
        if v is None:
            return
        if np is not None and isinstance(v, np.ndarray):
            for item in v.ravel().tolist():
                walk(item)
        elif isinstance(v, dict):
            # Do not flatten arbitrary dicts under a metric; these are usually metadata.
            return
        elif isinstance(v, (list, tuple)):
            for item in v:
                walk(item)
        else:
            x = to_float(v)
            if x is not None:
                out.append(x)

    walk(value)
    return out


def stat_dict(values: Sequence[float]) -> Dict[str, Any]:
    vals = [float(v) for v in values if is_finite_number(v)]
    if not vals:
        return {"mean": "", "median": "", "min": "", "max": "", "count": 0}
    return {
        "mean": mean(vals),
        "median": median(vals),
        "min": min(vals),
        "max": max(vals),
        "count": len(vals),
    }


def fmt(value: Any, ndigits: int = 4) -> str:
    if value is None or value == "":
        return ""
    if isinstance(value, bool):
        return str(value)
    if isinstance(value, (int, float)) and math.isfinite(float(value)):
        return f"{float(value):.{ndigits}f}"
    return str(value)


def recursive_find(data: Any, aliases: Sequence[str]) -> List[Any]:
    aliases_l = {a.lower() for a in aliases}
    found: List[Any] = []

    def walk(obj: Any) -> None:
        if isinstance(obj, dict):
            for k, v in obj.items():
                if isinstance(k, str) and k.lower() in aliases_l:
                    found.append(v)
                walk(v)
        elif isinstance(obj, list):
            for item in obj:
                walk(item)

    walk(data)
    return found


def first_scalar(data: Any, aliases: Sequence[str]) -> Any:
    for raw in recursive_find(data, aliases):
        vals = flatten_numeric(raw)
        if len(vals) == 1:
            return vals[0]
        if isinstance(raw, (str, bool, int, float)):
            return raw
    return ""


def load_json(path: Path) -> Optional[Any]:
    try:
        with path.open("r", encoding="utf-8") as handle:
            return json.load(handle)
    except Exception as exc:
        print(f"[WARN] Could not parse JSON {path}: {exc}", file=sys.stderr)
        return None


def json_metric_rows(json_path: Path, model_id: str, platform: str) -> List[Dict[str, Any]]:
    data = load_json(json_path)
    if data is None:
        return []
    rows: List[Dict[str, Any]] = []
    for metric, aliases in JSON_METRIC_ALIASES.items():
        values: List[float] = []
        raw_hits = recursive_find(data, aliases)
        for hit in raw_hits:
            values.extend(flatten_numeric(hit))
        s = stat_dict(values)
        rows.append({
            "model_id": model_id,
            "platform": platform,
            "source_file": str(json_path),
            "metric": metric,
            "mean": s["mean"],
            "median": s["median"],
            "min": s["min"],
            "max": s["max"],
            "count": s["count"],
        })
    return rows


def pairwise_json_rows(json_path: Path, model_id: str, platform: str) -> List[Dict[str, Any]]:
    data = load_json(json_path)
    if data is None:
        return []
    rows: List[Dict[str, Any]] = []
    for metric, aliases in PAIRWISE_ALIASES.items():
        for matrix in recursive_find(data, aliases):
            if not isinstance(matrix, list) or not matrix:
                continue
            for i, row in enumerate(matrix, start=1):
                if not isinstance(row, list):
                    continue
                for j, value in enumerate(row, start=1):
                    x = to_float(value)
                    if x is None:
                        continue
                    rows.append({
                        "model_id": model_id,
                        "platform": platform,
                        "source_file": str(json_path),
                        "metric": metric,
                        "row_chain_index": i,
                        "col_chain_index": j,
                        "value": x,
                    })
    return rows


def npz_metric_rows(npz_path: Path, model_id: str, platform: str) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    if np is None:
        rows.append({
            "model_id": model_id, "platform": platform, "source_file": str(npz_path),
            "array_name": "", "metric_guess": "", "mean": "", "median": "", "min": "",
            "max": "", "count": 0, "status": "numpy_not_installed",
        })
        return rows
    try:
        with np.load(str(npz_path), allow_pickle=False) as data:
            for key in data.files:
                arr = data[key]
                values = flatten_numeric(arr)
                s = stat_dict(values)
                metric_guess = guess_metric_from_name(key)
                rows.append({
                    "model_id": model_id,
                    "platform": platform,
                    "source_file": str(npz_path),
                    "array_name": key,
                    "metric_guess": metric_guess,
                    "mean": s["mean"],
                    "median": s["median"],
                    "min": s["min"],
                    "max": s["max"],
                    "count": s["count"],
                    "status": "ok",
                })
    except Exception as exc:
        rows.append({
            "model_id": model_id, "platform": platform, "source_file": str(npz_path),
            "array_name": "", "metric_guess": "", "mean": "", "median": "", "min": "",
            "max": "", "count": 0, "status": f"parse_error: {exc}",
        })
    return rows


def guess_metric_from_name(name: str) -> str:
    n = name.lower()
    if "pae" in n:
        return "pae"
    if "pde" in n:
        return "pde"
    if "plddt" in n:
        return "plddt"
    if "contact" in n:
        return "contact_prob"
    return "unknown"


def split_cif_row(row: str) -> List[str]:
    try:
        return shlex.split(row, posix=True)
    except Exception:
        return row.split()


def first_existing_column(tags: Sequence[str], names: Sequence[str]) -> Optional[int]:
    for name in names:
        try:
            return list(tags).index(name)
        except ValueError:
            continue
    return None


def iter_atom_site_rows(cif_path: Path) -> Iterable[Tuple[List[str], List[str]]]:
    with io.open(cif_path, "r", encoding="utf-8") as handle:
        lines = handle.readlines()

    i = 0
    while i < len(lines):
        if lines[i].strip() != "loop_":
            i += 1
            continue

        tags: List[str] = []
        j = i + 1
        while j < len(lines) and lines[j].strip().startswith("_"):
            tags.append(lines[j].strip())
            j += 1

        if not any(tag.startswith("_atom_site.") for tag in tags):
            i = j
            continue

        k = j
        while k < len(lines):
            row = lines[k].strip()
            if (not row) or row.startswith("#") or row == "loop_" or row.startswith("data_") or row.startswith("_"):
                break
            parts = split_cif_row(row)
            if len(parts) >= len(tags):
                yield tags, parts
            k += 1
        i = k


def residue_sort_key(row: Dict[str, Any]) -> Tuple[str, int, Any]:
    r = row.get("residue_sort")
    if isinstance(r, int):
        return str(row.get("chain", "")), 0, r
    return str(row.get("chain", "")), 1, str(r)


def extract_cif_residues(cif_path: Path) -> List[Dict[str, Any]]:
    residue_atoms: Dict[Tuple[str, str, Any, str], List[float]] = defaultdict(list)
    for tags, parts in iter_atom_site_rows(cif_path):
        idx_chain = first_existing_column(tags, ["_atom_site.auth_asym_id", "_atom_site.label_asym_id"])
        idx_resid = first_existing_column(tags, ["_atom_site.auth_seq_id", "_atom_site.label_seq_id"])
        idx_comp = first_existing_column(tags, ["_atom_site.auth_comp_id", "_atom_site.label_comp_id"])
        idx_b = first_existing_column(tags, ["_atom_site.B_iso_or_equiv"])
        if idx_chain is None or idx_resid is None or idx_comp is None or idx_b is None:
            continue
        chain = parts[idx_chain]
        resid_raw = parts[idx_resid]
        comp3 = parts[idx_comp].upper()
        b = to_float(parts[idx_b])
        if b is None:
            continue
        try:
            resid_sort: Any = int(resid_raw)
        except ValueError:
            resid_sort = resid_raw
        residue_atoms[(chain, resid_raw, resid_sort, comp3)].append(b)

    rows: List[Dict[str, Any]] = []
    for (chain, resid_raw, resid_sort, comp3), vals in residue_atoms.items():
        rows.append({
            "chain": chain,
            "residue_number": resid_raw,
            "residue_sort": resid_sort,
            "residue_name_3": comp3,
            "residue_name_1": AA3_TO_AA1.get(comp3, "X"),
            "mean_plddt": mean(vals),
            "median_plddt": median(vals),
            "atom_count": len(vals),
            "atom_plddts": vals,
        })
    return sorted(rows, key=residue_sort_key)


def extract_pdb_residues(pdb_path: Path) -> List[Dict[str, Any]]:
    residue_atoms: Dict[Tuple[str, str, Any, str], List[float]] = defaultdict(list)
    try:
        with pdb_path.open("r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if not line.startswith(("ATOM  ", "HETATM")):
                    continue
                comp3 = line[17:20].strip().upper()
                chain = line[21].strip() or "_"
                resid_raw = line[22:26].strip()
                ins = line[26].strip()
                if ins:
                    resid_raw = f"{resid_raw}{ins}"
                try:
                    b = float(line[60:66].strip())
                except ValueError:
                    continue
                try:
                    resid_sort: Any = int(re.match(r"-?\d+", resid_raw).group(0))  # type: ignore
                except Exception:
                    resid_sort = resid_raw
                residue_atoms[(chain, resid_raw, resid_sort, comp3)].append(b)
    except Exception as exc:
        print(f"[WARN] Could not parse PDB {pdb_path}: {exc}", file=sys.stderr)
        return []

    rows: List[Dict[str, Any]] = []
    for (chain, resid_raw, resid_sort, comp3), vals in residue_atoms.items():
        rows.append({
            "chain": chain,
            "residue_number": resid_raw,
            "residue_sort": resid_sort,
            "residue_name_3": comp3,
            "residue_name_1": AA3_TO_AA1.get(comp3, "X"),
            "mean_plddt": mean(vals),
            "median_plddt": median(vals),
            "atom_count": len(vals),
            "atom_plddts": vals,
        })
    return sorted(rows, key=residue_sort_key)


def extract_structure_residues(path: Path) -> List[Dict[str, Any]]:
    ext = path.suffix.lower()
    if ext in {".cif", ".mmcif"}:
        return extract_cif_residues(path)
    if ext == ".pdb":
        return extract_pdb_residues(path)
    return []


def chain_groups(residues: Sequence[Dict[str, Any]]) -> Dict[str, List[Dict[str, Any]]]:
    groups: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
    for row in residues:
        groups[str(row["chain"])].append(row)
    return dict(groups)


def structure_summary_rows(path: Path, model_id: str, platform: str, residues: Sequence[Dict[str, Any]]) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    per_residue: List[Dict[str, Any]] = []
    for r in residues:
        per_residue.append({
            "model_id": model_id,
            "platform": platform,
            "source_file": str(path),
            "chain": r["chain"],
            "residue_number": r["residue_number"],
            "residue_name_3": r["residue_name_3"],
            "residue_name_1": r["residue_name_1"],
            "mean_plddt": r["mean_plddt"],
            "median_plddt": r["median_plddt"],
            "quality_band": quality_band(r["mean_plddt"]),
            "atom_count": r["atom_count"],
        })

    per_chain: List[Dict[str, Any]] = []
    for chain, rows in chain_groups(residues).items():
        vals = [r["mean_plddt"] for r in rows]
        s = stat_dict(vals)
        per_chain.append({
            "model_id": model_id,
            "platform": platform,
            "source_file": str(path),
            "chain": chain,
            "mean_plddt": s["mean"],
            "median_plddt": s["median"],
            "min_plddt": s["min"],
            "max_plddt": s["max"],
            "quality_band": quality_band(s["mean"]),
            "residue_count": len(rows),
            "atom_count": sum(int(r["atom_count"]) for r in rows),
        })
    return per_residue, per_chain


def quality_band(plddt: Any) -> str:
    x = to_float(plddt)
    if x is None:
        return "unknown"
    for threshold, label in QUALITY_BANDS:
        if x >= threshold:
            return label
    return "unknown"


def motif_hits_exact(sequence: str, motif: str) -> List[Tuple[int, int, str]]:
    motif_u = motif.upper().strip()
    seq_u = sequence.upper()
    hits: List[Tuple[int, int, str]] = []
    start = 0
    while motif_u:
        pos = seq_u.find(motif_u, start)
        if pos < 0:
            break
        hits.append((pos, pos + len(motif_u), sequence[pos:pos + len(motif_u)]))
        start = pos + 1
    return hits


def motif_hits_regex(sequence: str, pattern: str) -> List[Tuple[int, int, str]]:
    hits: List[Tuple[int, int, str]] = []
    try:
        for m in re.finditer(pattern, sequence, flags=re.IGNORECASE):
            if m.end() > m.start():
                hits.append((m.start(), m.end(), m.group(0)))
    except re.error as exc:
        print(f"[WARN] Invalid motif regex {pattern!r}: {exc}", file=sys.stderr)
    return hits


def motif_rows_for_structure(path: Path, model_id: str, platform: str, residues: Sequence[Dict[str, Any]], motifs: Sequence[str], motif_mode: str) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for chain, chain_residues in chain_groups(residues).items():
        seq = "".join(r["residue_name_1"] for r in chain_residues)
        for motif in motifs:
            hits = motif_hits_regex(seq, motif) if motif_mode == "regex" else motif_hits_exact(seq, motif)
            for hit_no, (start0, end0, match_seq) in enumerate(hits, start=1):
                hit_res = chain_residues[start0:end0]
                atom_vals: List[float] = []
                for r in hit_res:
                    atom_vals.extend(r.get("atom_plddts", []))
                s = stat_dict(atom_vals)
                residue_numbers = ";".join(str(r["residue_number"]) for r in hit_res)
                residue_labels = ";".join(f"{r['residue_name_1']}{r['residue_number']}" for r in hit_res)
                rows.append({
                    "model_id": model_id,
                    "platform": platform,
                    "source_file": str(path),
                    "chain": chain,
                    "motif_query": motif,
                    "motif_match": match_seq,
                    "motif_mode": motif_mode,
                    "hit_number": hit_no,
                    "sequence_start_1based": start0 + 1,
                    "sequence_end_1based": end0,
                    "residue_numbers": residue_numbers,
                    "residue_labels": residue_labels,
                    "residue_count": len(hit_res),
                    "atom_count": s["count"],
                    "mean_plddt": s["mean"],
                    "median_plddt": s["median"],
                    "min_plddt": s["min"],
                    "max_plddt": s["max"],
                    "quality_band": quality_band(s["mean"]),
                })
    return rows


def normalize_stem(path: Path) -> str:
    stem = path.stem.lower()
    # Boltz-style and AF3-style filename normalization.
    prefixes = ["fold_", "confidence_", "affinity_", "pae_", "pde_", "plddt_"]
    for prefix in prefixes:
        if stem.startswith(prefix):
            stem = stem[len(prefix):]
    suffix_patterns = [
        r"_model_\d+$", r"_model$", r"_rank_\d+$", r"_seed_\d+$",
        r"_confidences$", r"_summary_confidences$", r"_full_data_\d+$",
        r"_data_\d+$", r"_full_data$", r"_data$", r"_prediction$",
        r"_scores$", r"_score$", r"_\d+$",
    ]
    changed = True
    while changed:
        changed = False
        for pat in suffix_patterns:
            new = re.sub(pat, "", stem)
            if new != stem and new:
                stem = new
                changed = True
    stem = re.sub(r"[^a-z0-9]+", "_", stem).strip("_")
    return stem or path.stem.lower()


def stable_model_id(group_key: str, paths: Sequence[Path]) -> str:
    digest_src = "|".join(sorted(str(p) for p in paths)) or group_key
    digest = hashlib.md5(digest_src.encode("utf-8")).hexdigest()[:6]
    base = re.sub(r"[^A-Za-z0-9_.-]+", "_", group_key)[:60] or "model"
    return f"{base}_{digest}"


def detect_platform(paths: Sequence[Path], forced: str = "auto") -> str:
    if forced and forced != "auto":
        return forced
    text = " ".join(str(p).lower() for p in paths)
    names = " ".join(p.name.lower() for p in paths)
    if "boltz" in text or "affinity_" in names or "confidence_" in names:
        return "boltz_or_boltz2"
    if "alphafold" in text or "af3" in text or "summary_confidences" in names or "_confidences" in names or "fold_" in names:
        return "alphafold3"
    if "chai" in text:
        return "chai_or_chai1"
    return "generic"


def discover_files(inputs: Sequence[str], recursive: bool = True) -> List[Path]:
    allowed = {".cif", ".mmcif", ".pdb", ".json", ".npz"}
    files: List[Path] = []
    for item in inputs:
        # Support shell-style glob patterns passed quoted.
        matches = glob.glob(item, recursive=recursive)
        if matches:
            for m in matches:
                p = Path(m)
                if p.is_file() and p.suffix.lower() in allowed:
                    files.append(p.resolve())
                elif p.is_dir():
                    files.extend(iter_supported_files(p, recursive=recursive))
            continue
        p = Path(item)
        if p.is_file() and p.suffix.lower() in allowed:
            files.append(p.resolve())
        elif p.is_dir():
            files.extend(iter_supported_files(p, recursive=recursive))
    return sorted(set(files), key=lambda x: str(x))


def iter_supported_files(directory: Path, recursive: bool = True) -> Iterable[Path]:
    patterns = ["*.cif", "*.mmcif", "*.pdb", "*.json", "*.npz"]
    for pat in patterns:
        yield from directory.rglob(pat) if recursive else directory.glob(pat)


def group_files(files: Sequence[Path]) -> Dict[str, Dict[str, List[Path]]]:
    groups: Dict[str, Dict[str, List[Path]]] = defaultdict(lambda: {"structures": [], "jsons": [], "npzs": []})
    for p in files:
        key = normalize_stem(p)
        ext = p.suffix.lower()
        if ext in {".cif", ".mmcif", ".pdb"}:
            groups[key]["structures"].append(p)
        elif ext == ".json":
            groups[key]["jsons"].append(p)
        elif ext == ".npz":
            groups[key]["npzs"].append(p)
    return groups


def write_tsv(path: Path, rows: Sequence[Dict[str, Any]], fieldnames: Optional[List[str]] = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fields: List[str] = []
        for row in rows:
            for key in row.keys():
                if key not in fields:
                    fields.append(key)
        fieldnames = fields
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            clean = {k: fmt(v) if isinstance(v, float) else v for k, v in row.items()}
            writer.writerow(clean)


def read_score_from_rows(rows: Sequence[Dict[str, Any]], metric: str, field: str = "mean") -> Optional[float]:
    vals: List[float] = []
    for r in rows:
        if r.get("metric") == metric:
            x = to_float(r.get(field))
            if x is not None:
                vals.append(x)
    if not vals:
        return None
    return mean(vals)


def summarise_model(model_id: str, platform: str, group_key: str, group: Dict[str, List[Path]],
                    chain_rows: Sequence[Dict[str, Any]], motif_rows: Sequence[Dict[str, Any]],
                    json_rows: Sequence[Dict[str, Any]], npz_rows: Sequence[Dict[str, Any]]) -> Dict[str, Any]:
    plddt_vals = [to_float(r.get("mean_plddt")) for r in chain_rows]
    plddt_vals = [v for v in plddt_vals if v is not None]
    motif_vals = [to_float(r.get("mean_plddt")) for r in motif_rows]
    motif_vals = [v for v in motif_vals if v is not None]

    structure_mean = mean(plddt_vals) if plddt_vals else None
    motif_mean = mean(motif_vals) if motif_vals else None
    json_pae = read_score_from_rows(json_rows, "pae", "mean")
    json_plddt = read_score_from_rows(json_rows, "plddt", "mean")
    json_iptm = read_score_from_rows(json_rows, "iptm", "mean")
    json_ptm = read_score_from_rows(json_rows, "ptm", "mean")
    ranking_score = read_score_from_rows(json_rows, "ranking_score", "mean")
    affinity = read_score_from_rows(json_rows, "affinity_pred_value", "mean")
    affinity_prob = read_score_from_rows(json_rows, "affinity_probability_binary", "mean")

    heuristic = integrated_quality_score(
        structure_plddt=structure_mean,
        motif_plddt=motif_mean,
        json_plddt=json_plddt,
        pae=json_pae,
        iptm=json_iptm,
        ptm=json_ptm,
        ranking_score=ranking_score,
    )

    return {
        "model_id": model_id,
        "group_key": group_key,
        "platform": platform,
        "n_structure_files": len(group.get("structures", [])),
        "n_json_files": len(group.get("jsons", [])),
        "n_npz_files": len(group.get("npzs", [])),
        "structure_files": ";".join(str(p) for p in group.get("structures", [])),
        "json_files": ";".join(str(p) for p in group.get("jsons", [])),
        "npz_files": ";".join(str(p) for p in group.get("npzs", [])),
        "structure_mean_plddt": structure_mean,
        "structure_quality_band": quality_band(structure_mean),
        "motif_mean_plddt": motif_mean,
        "motif_quality_band": quality_band(motif_mean),
        "json_mean_plddt": json_plddt,
        "json_mean_pae": json_pae,
        "json_mean_iptm": json_iptm,
        "json_mean_ptm": json_ptm,
        "ranking_score": ranking_score,
        "affinity_pred_value": affinity,
        "affinity_probability_binary": affinity_prob,
        "integrated_quality_score_0_100": heuristic,
        "integrated_rank_note": "heuristic; prioritize experimental context and visual inspection",
    }


def integrated_quality_score(structure_plddt: Optional[float], motif_plddt: Optional[float],
                             json_plddt: Optional[float], pae: Optional[float], iptm: Optional[float],
                             ptm: Optional[float], ranking_score: Optional[float]) -> Optional[float]:
    components: List[Tuple[float, float]] = []
    # Normalize pLDDT-like values to 0-100.
    if structure_plddt is not None:
        components.append((0.35, max(0.0, min(100.0, structure_plddt))))
    if motif_plddt is not None:
        components.append((0.25, max(0.0, min(100.0, motif_plddt))))
    elif json_plddt is not None:
        components.append((0.20, max(0.0, min(100.0, json_plddt))))
    if iptm is not None:
        components.append((0.15, max(0.0, min(1.0, iptm)) * 100.0))
    if ptm is not None:
        components.append((0.10, max(0.0, min(1.0, ptm)) * 100.0))
    if ranking_score is not None:
        # Some tools already use 0-1, others 0-100. Infer scale safely.
        score = ranking_score * 100.0 if ranking_score <= 1.5 else ranking_score
        components.append((0.10, max(0.0, min(100.0, score))))
    if pae is not None:
        # PAE is better when lower. Penalize only softly; max useful range here is 0-30 Å.
        pae_score = max(0.0, min(100.0, 100.0 - (pae / 30.0 * 100.0)))
        components.append((0.10, pae_score))
    if not components:
        return None
    wsum = sum(w for w, _ in components)
    return sum(w * x for w, x in components) / wsum


def write_excel(out_xlsx: Path, sheet_rows: Dict[str, Sequence[Dict[str, Any]]]) -> None:
    if pd is None:
        print("[WARN] pandas not installed; skipping Excel workbook.", file=sys.stderr)
        return
    out_xlsx.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(out_xlsx, engine="openpyxl") as writer:  # type: ignore
        for name, rows in sheet_rows.items():
            safe = safe_sheet_name(name)
            if not rows:
                df = pd.DataFrame()
            else:
                df = pd.DataFrame(list(rows))
            df.to_excel(writer, sheet_name=safe, index=False)


def safe_sheet_name(name: str) -> str:
    return re.sub(r"[\\/*?:\[\]]", "_", name)[:31] or "Sheet"


def plot_plddt_profile(model_id: str, residues: Sequence[Dict[str, Any]], out_png: Path) -> None:
    if plt is None or not residues:
        return
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(10, 3.2))
    offset = 0
    x_all: List[int] = []
    y_all: List[float] = []
    chain_boundaries: List[Tuple[int, str]] = []
    for chain, rows in chain_groups(residues).items():
        xs = list(range(offset + 1, offset + len(rows) + 1))
        ys = [float(r["mean_plddt"]) for r in rows]
        ax.plot(xs, ys, linewidth=1.2, label=f"Chain {chain}")
        x_all.extend(xs)
        y_all.extend(ys)
        chain_boundaries.append((offset + len(rows), chain))
        offset += len(rows)
    ax.axhline(90, linestyle="--", linewidth=0.8)
    ax.axhline(70, linestyle="--", linewidth=0.8)
    ax.axhline(50, linestyle="--", linewidth=0.8)
    ax.set_ylim(0, 100)
    ax.set_xlabel("Residue index across concatenated chains")
    ax.set_ylabel("Mean pLDDT")
    ax.set_title(f"pLDDT profile: {model_id}")
    if len(chain_boundaries) <= 12:
        for boundary, chain in chain_boundaries[:-1]:
            ax.axvline(boundary + 0.5, linewidth=0.5)
    if len(chain_boundaries) <= 8:
        ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(out_png, dpi=300)
    plt.close(fig)


def plot_model_summary(summary_rows: Sequence[Dict[str, Any]], out_png: Path) -> None:
    if plt is None or not summary_rows:
        return
    rows = [r for r in summary_rows if to_float(r.get("integrated_quality_score_0_100")) is not None]
    if not rows:
        return
    rows = sorted(rows, key=lambda r: float(r["integrated_quality_score_0_100"]), reverse=True)[:30]
    labels = [str(r["model_id"])[:35] for r in rows]
    values = [float(r["integrated_quality_score_0_100"]) for r in rows]
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(max(8, len(rows) * 0.35), 4.5))
    ax.bar(range(len(values)), values)
    ax.set_ylim(0, 100)
    ax.set_ylabel("Integrated quality score")
    ax.set_title("Top model quality summary")
    ax.set_xticks(range(len(values)))
    ax.set_xticklabels(labels, rotation=75, ha="right", fontsize=7)
    fig.tight_layout()
    fig.savefig(out_png, dpi=300)
    plt.close(fig)


def convert_cif_to_pdb(cif_path: Path, out_dir: Path) -> Optional[Path]:
    if MMCIFParser is None or PDBIO is None:
        print("[WARN] Biopython not installed; cannot convert mmCIF to PDB.", file=sys.stderr)
        return None
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{cif_path.stem}.pdb"
    try:
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure(cif_path.stem, str(cif_path))
        io_obj = PDBIO()
        io_obj.set_structure(structure)
        io_obj.save(str(out_path))
        return out_path
    except Exception as exc:
        print(f"[WARN] Failed to convert {cif_path} to PDB: {exc}", file=sys.stderr)
        return None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Evaluate AlphaFold3, Boltz/Boltz-2, and generic deep-learning predicted structures."
    )
    parser.add_argument("inputs", nargs="+", help="Input directories, files, or glob patterns containing .cif/.mmcif/.pdb/.json/.npz outputs.")
    parser.add_argument("--outdir", default="dl_structure_eval_results", help="Output directory.")
    parser.add_argument("--platform", default="auto", choices=["auto", "alphafold3", "boltz2", "boltz_or_boltz2", "chai_or_chai1", "generic"], help="Force platform label or use auto-detection.")
    parser.add_argument("--motifs", default="HGG,GESAG", help="Comma-separated motifs or regex patterns to evaluate. Default: HGG,GESAG")
    parser.add_argument("--motif-mode", choices=["exact", "regex"], default="exact", help="Treat --motifs as exact peptide strings or regular expressions.")
    parser.add_argument("--recursive", action="store_true", default=True, help="Recursively search input directories. Default: true.")
    parser.add_argument("--no-recursive", dest="recursive", action="store_false", help="Do not recursively search input directories.")
    parser.add_argument("--convert-pdb", action="store_true", help="Convert discovered mmCIF files to PDB using Biopython.")
    parser.add_argument("--write-excel", action="store_true", help="Write integrated Excel workbook in addition to TSV tables.")
    parser.add_argument("--no-plots", action="store_true", help="Skip pLDDT profile and summary plots.")
    parser.add_argument("--fail-on-empty", action="store_true", help="Raise an error when no supported files are discovered.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir).resolve()
    table_dir = outdir / "tables"
    plot_dir = outdir / "plots"
    pdb_dir = outdir / "pdb"
    table_dir.mkdir(parents=True, exist_ok=True)

    motifs = [m.strip() for m in args.motifs.split(",") if m.strip()]
    files = discover_files(args.inputs, recursive=args.recursive)
    if not files:
        msg = "No supported files found. Expected .cif, .mmcif, .pdb, .json, or .npz."
        if args.fail_on_empty:
            raise SystemExit(msg)
        print(f"[WARN] {msg}", file=sys.stderr)

    groups = group_files(files)

    all_model_summary: List[Dict[str, Any]] = []
    all_residue: List[Dict[str, Any]] = []
    all_chain: List[Dict[str, Any]] = []
    all_motif: List[Dict[str, Any]] = []
    all_json: List[Dict[str, Any]] = []
    all_pairwise: List[Dict[str, Any]] = []
    all_npz: List[Dict[str, Any]] = []
    all_conversions: List[Dict[str, Any]] = []

    for group_key, group in sorted(groups.items()):
        paths = group["structures"] + group["jsons"] + group["npzs"]
        model_id = stable_model_id(group_key, paths)
        platform = detect_platform(paths, forced=args.platform)
        model_residue: List[Dict[str, Any]] = []
        model_chain: List[Dict[str, Any]] = []
        model_motif: List[Dict[str, Any]] = []
        model_json: List[Dict[str, Any]] = []
        model_pairwise: List[Dict[str, Any]] = []
        model_npz: List[Dict[str, Any]] = []

        for struct_path in group["structures"]:
            residues = extract_structure_residues(struct_path)
            if not residues:
                print(f"[WARN] No residue pLDDT parsed from {struct_path}", file=sys.stderr)
            per_res, per_chain = structure_summary_rows(struct_path, model_id, platform, residues)
            motif_rows = motif_rows_for_structure(struct_path, model_id, platform, residues, motifs, args.motif_mode)
            model_residue.extend(per_res)
            model_chain.extend(per_chain)
            model_motif.extend(motif_rows)
            all_residue.extend(per_res)
            all_chain.extend(per_chain)
            all_motif.extend(motif_rows)
            if not args.no_plots and residues:
                plot_plddt_profile(model_id, residues, plot_dir / f"{model_id}_{struct_path.stem}_plddt_profile.png")
            if args.convert_pdb and struct_path.suffix.lower() in {".cif", ".mmcif"}:
                out_pdb = convert_cif_to_pdb(struct_path, pdb_dir)
                all_conversions.append({
                    "model_id": model_id,
                    "source_cif": str(struct_path),
                    "output_pdb": str(out_pdb) if out_pdb else "",
                    "status": "ok" if out_pdb else "failed_or_biopython_missing",
                })

        for json_path in group["jsons"]:
            rows = json_metric_rows(json_path, model_id, platform)
            pair_rows = pairwise_json_rows(json_path, model_id, platform)
            model_json.extend(rows)
            model_pairwise.extend(pair_rows)
            all_json.extend(rows)
            all_pairwise.extend(pair_rows)

        for npz_path in group["npzs"]:
            rows = npz_metric_rows(npz_path, model_id, platform)
            model_npz.extend(rows)
            all_npz.extend(rows)

        summary = summarise_model(model_id, platform, group_key, group, model_chain, model_motif, model_json, model_npz)
        all_model_summary.append(summary)

    # Sort summary by integrated score when available.
    all_model_summary = sorted(
        all_model_summary,
        key=lambda r: (to_float(r.get("integrated_quality_score_0_100")) is not None,
                       to_float(r.get("integrated_quality_score_0_100")) or -1),
        reverse=True,
    )

    write_tsv(table_dir / "model_summary.tsv", all_model_summary)
    write_tsv(table_dir / "per_residue_plddt.tsv", all_residue)
    write_tsv(table_dir / "per_chain_plddt.tsv", all_chain)
    write_tsv(table_dir / "motif_plddt.tsv", all_motif)
    write_tsv(table_dir / "json_metrics_long.tsv", all_json)
    write_tsv(table_dir / "json_pairwise_long.tsv", all_pairwise)
    write_tsv(table_dir / "npz_metrics_long.tsv", all_npz)
    write_tsv(table_dir / "cif_to_pdb_conversions.tsv", all_conversions)

    if not args.no_plots:
        plot_model_summary(all_model_summary, plot_dir / "integrated_quality_summary.png")

    if args.write_excel:
        write_excel(outdir / "deep_learning_structure_evaluation.xlsx", {
            "model_summary": all_model_summary,
            "per_chain_plddt": all_chain,
            "motif_plddt": all_motif,
            "json_metrics_long": all_json,
            "json_pairwise_long": all_pairwise,
            "npz_metrics_long": all_npz,
            "per_residue_plddt": all_residue,
            "cif_to_pdb_conversions": all_conversions,
        })

    print(f"[DONE] Evaluated {len(all_model_summary)} model group(s).")
    print(f"[DONE] Tables: {table_dir}")
    if not args.no_plots:
        print(f"[DONE] Plots:  {plot_dir}")
    if args.write_excel:
        print(f"[DONE] Excel:  {outdir / 'deep_learning_structure_evaluation.xlsx'}")


if __name__ == "__main__":
    main()
