#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Integrated phylogeny, motif annotation, topology parsing, and hydrolase triage
for Schistosoma lipid-metabolism candidate proteins.

Inputs
------
- One or more family-specific FASTA files from the Pfam/domain extraction step.
- Optional curated motif library as XLSX/CSV/TSV.
- Optional SignalP-6.0 TSV output.
- Optional DeepTMHMM Markdown, GFF3, or 3-line topology output.
- Optional CCTOP TSV/text output.

Main outputs
------------
- Per-family MAFFT alignments, IQ-TREE ML trees, sequence frequency matrices,
  conservation summaries, tree plots, and sequence logos.
- Motif hit tables with 1-based positions and matched residues.
- Integrated topology/signal-peptide/GPI-like feature table.
- Hydrolase classification table and strip diagrams.
- Combined Excel workbook and optional PDF evidence report.
"""

from __future__ import annotations

import argparse
import csv
import html
import math
import os
import re
import shutil
import subprocess
import sys
import warnings
from collections import defaultdict, Counter
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import pandas as pd

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
except Exception as exc:  # pragma: no cover
    raise RuntimeError("matplotlib is required for figures. Install with: pip install matplotlib") from exc

try:
    from Bio import AlignIO, Phylo, SeqIO
except Exception as exc:  # pragma: no cover
    raise RuntimeError("Biopython is required. Install with: pip install biopython") from exc

try:
    import logomaker  # optional, but recommended
except Exception:  # pragma: no cover
    logomaker = None

try:
    from reportlab.lib import colors
    from reportlab.lib.pagesizes import landscape, letter
    from reportlab.lib.styles import getSampleStyleSheet
    from reportlab.platypus import Image as RLImage
    from reportlab.platypus import Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle
except Exception:  # pragma: no cover
    colors = None
    landscape = None
    letter = None
    getSampleStyleSheet = None
    RLImage = None
    Paragraph = None
    SimpleDocTemplate = None
    Spacer = None
    Table = None
    TableStyle = None


# ---------------------------------------------------------------------------
# Explicit thresholds and default heuristic parameters
# ---------------------------------------------------------------------------
AA20 = list("ACDEFGHIKLMNPQRSTVWY")

KYTE_DOOLITTLE = {
    "I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5, "M": 1.9,
    "A": 1.8, "G": -0.4, "T": -0.7, "S": -0.8, "W": -0.9, "Y": -1.3,
    "P": -1.6, "H": -3.2, "E": -3.5, "Q": -3.5, "D": -3.5, "N": -3.5,
    "K": -3.9, "R": -4.5,
}

HYDRO = set("IVLFCMWY")
SMALL = set("ASGNTQDE")
POSITIVE = set("KR")

DEFAULT_KD_TM_WINDOW = 19
DEFAULT_KD_TM_THRESHOLD = 1.60
DEFAULT_SIGNAL_WINDOW = 19
DEFAULT_SIGNAL_H_REGION_THRESHOLD = 1.60
DEFAULT_N_TERM_SCAN = 70
DEFAULT_GPI_TAIL_SCAN = 50
DEFAULT_GPI_HYDROPHOBIC_FRACTION = 0.70
DEFAULT_GPI_TAIL_MIN = 18
DEFAULT_GPI_TAIL_MAX = 25
DEFAULT_SUPPORT_DISPLAY = 70.0
DEFAULT_COLLAPSE_SUPPORT = 50.0

BUILTIN_MOTIFS = [
    {
        "Family": "ABHD/HSL",
        "Motif_ID": "HYDROLASE_GxSxG",
        "Consensus_Motif": "GxSxG",
        "Regex": r"G[A-Z]S[A-Z]G",
        "Functional_Region": "α/β-hydrolase nucleophile loop",
        "Conservation": "diagnostic",
    },
    {
        "Family": "ABHD/HSL",
        "Motif_ID": "HYDROLASE_HGG_like",
        "Consensus_Motif": "HGG/HGGG/HGS/HGA",
        "Regex": r"HGGG?|HGS|HGA",
        "Functional_Region": "HGG-like oxyanion motif",
        "Conservation": "diagnostic",
    },
    {
        "Family": "Secretory/ecto",
        "Motif_ID": "N_GLYCO",
        "Consensus_Motif": "N-x-S/T, x != P",
        "Regex": r"N[^P][ST]",
        "Functional_Region": "putative N-linked glycosylation sequon",
        "Conservation": "context-dependent",
    },
]


@dataclass
class SeqRecordLite:
    family: str
    seq_id: str
    description: str
    sequence: str
    fasta_path: str


@dataclass
class MotifHit:
    family: str
    seq_id: str
    motif_family: str
    motif_id: str
    consensus: str
    regex: str
    functional_region: str
    conservation: str
    start_1based: int
    end_1based: int
    matched_sequence: str


@dataclass
class TopologyCall:
    has_sp: Optional[bool] = None
    sp_cut: Optional[int] = None
    signal_range: Optional[Tuple[int, int]] = None
    tm_ranges: Optional[List[Tuple[int, int]]] = None
    source: str = "none"
    label: str = ""

    @property
    def tm_count(self) -> int:
        return len(self.tm_ranges or [])


# ---------------------------------------------------------------------------
# FASTA utilities
# ---------------------------------------------------------------------------

def sanitize_name(name: str) -> str:
    name = re.sub(r"\.(fa|faa|fasta|fas)$", "", name, flags=re.I)
    name = re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("_")
    return name or "family"


def parse_fasta(path: Path, family: Optional[str] = None) -> List[SeqRecordLite]:
    fam = family or sanitize_name(path.name)
    records: List[SeqRecordLite] = []
    for rec in SeqIO.parse(str(path), "fasta"):
        seq = str(rec.seq).upper().replace(" ", "").replace("*", "")
        if not seq:
            continue
        records.append(SeqRecordLite(fam, rec.id, rec.description, seq, str(path)))
    return records


def collect_fasta_paths(args: argparse.Namespace) -> List[Path]:
    paths: List[Path] = []
    if args.fasta:
        for p in args.fasta:
            paths.append(Path(p))
    if args.input_dir:
        for ext in ("*.fa", "*.faa", "*.fasta", "*.fas"):
            paths.extend(sorted(Path(args.input_dir).glob(ext)))
    seen = []
    for p in paths:
        if p.exists() and p.is_file() and p not in seen:
            seen.append(p)
    if not seen:
        raise FileNotFoundError("No FASTA files found. Use --fasta or --input-dir.")
    return seen


# ---------------------------------------------------------------------------
# Motif library parsing and motif scanning
# ---------------------------------------------------------------------------

def infer_motif_family(motif_id: str, current_section: str = "") -> str:
    mid = (motif_id or "").strip()
    sec = (current_section or "").strip().strip(":")
    if mid:
        prefix = re.split(r"[_\-\s]", mid)[0].strip()
        if prefix:
            if prefix.upper().startswith("CD36"):
                return "CD36RP"
            if prefix.upper().startswith("NPC1"):
                return "NPC1"
            if prefix.upper().startswith("NPC2"):
                return "NPC2"
            if prefix.upper().startswith(("HSL", "ABHD", "HYDROLASE")):
                return "ABHD/HSL"
            return prefix
    if sec:
        if "NPC1" in sec.upper():
            return "NPC1"
        if "NPC2" in sec.upper():
            return "NPC2"
        if "CD36" in sec.upper():
            return "CD36RP"
        if "HYDROLASE" in sec.upper() or "HSL" in sec.upper():
            return "ABHD/HSL"
        return sec
    return "Unspecified"


def normalize_motif_df(df: pd.DataFrame) -> pd.DataFrame:
    """Normalize either a tidy table or an Excel sheet with section rows."""
    raw = df.copy()
    raw = raw.dropna(how="all")
    if raw.empty:
        return pd.DataFrame(columns=["Family", "Motif_ID", "Consensus_Motif", "Regex", "Functional_Region", "Conservation"])

    # If columns are unnamed integers, look for the header row containing "Motif ID".
    if not any(str(c).lower().replace(" ", "_") in {"motif_id", "motifid"} for c in raw.columns):
        header_idx = None
        for idx, row in raw.iterrows():
            vals = [str(x).strip().lower() for x in row.tolist()]
            if "motif id" in vals or "motif_id" in vals:
                header_idx = idx
                break
        if header_idx is not None:
            header_vals = [str(x).strip() if not pd.isna(x) else f"col{i}" for i, x in enumerate(raw.loc[header_idx].tolist())]
            raw = raw.loc[header_idx + 1:].copy()
            raw.columns = header_vals[:len(raw.columns)]

    def pick_col(candidates: Sequence[str]) -> Optional[str]:
        norm = {str(c).strip().lower().replace(" ", "_").replace("-", "_"): c for c in raw.columns}
        for cand in candidates:
            if cand in norm:
                return norm[cand]
        return None

    motif_col = pick_col(["motif_id", "motifid", "motif"])
    cons_col = pick_col(["consensus_motif", "consensus", "sequence_motif"])
    regex_col = pick_col(["regex_compatible_pattern", "regex", "pattern", "regex_pattern"])
    region_col = pick_col(["likely_functional_region", "functional_region", "region"])
    conservation_col = pick_col(["conservation", "conservation_level"])
    family_col = pick_col(["family", "motif_family", "protein_family"])

    rows = []
    current_section = ""
    for _, row in raw.iterrows():
        vals = ["" if pd.isna(x) else str(x).strip() for x in row.tolist()]
        first = vals[0] if vals else ""
        second_nonempty = sum(bool(v) for v in vals[1:])
        if first and second_nonempty == 0 and not re.search(r"_M\d+|HYDROLASE|GLYCO", first, flags=re.I):
            current_section = first
            continue

        motif_id = ("" if motif_col is None or pd.isna(row.get(motif_col, "")) else str(row.get(motif_col)).strip())
        regex = ("" if regex_col is None or pd.isna(row.get(regex_col, "")) else str(row.get(regex_col)).strip())
        if not motif_id or not regex or motif_id.lower() in {"motif id", "motif_id"}:
            continue

        consensus = "" if cons_col is None or pd.isna(row.get(cons_col, "")) else str(row.get(cons_col)).strip()
        region = "" if region_col is None or pd.isna(row.get(region_col, "")) else str(row.get(region_col)).strip()
        conservation = "" if conservation_col is None or pd.isna(row.get(conservation_col, "")) else str(row.get(conservation_col)).strip()
        motif_family = "" if family_col is None or pd.isna(row.get(family_col, "")) else str(row.get(family_col)).strip()
        if not motif_family:
            motif_family = infer_motif_family(motif_id, current_section)
        rows.append({
            "Family": motif_family,
            "Motif_ID": motif_id,
            "Consensus_Motif": consensus,
            "Regex": regex.replace(" ", ""),
            "Functional_Region": region,
            "Conservation": conservation,
        })
    out = pd.DataFrame(rows)
    if out.empty:
        return pd.DataFrame(columns=["Family", "Motif_ID", "Consensus_Motif", "Regex", "Functional_Region", "Conservation"])
    return out.drop_duplicates()


def load_motif_library(path: Optional[Path], include_builtin: bool = True) -> pd.DataFrame:
    frames: List[pd.DataFrame] = []
    if path:
        if not path.exists():
            raise FileNotFoundError(f"Motif library not found: {path}")
        suffix = path.suffix.lower()
        if suffix in {".xlsx", ".xls"}:
            book = pd.read_excel(path, sheet_name=None, header=None, engine=None)
            for _, df in book.items():
                frames.append(normalize_motif_df(df))
        elif suffix in {".tsv", ".txt"}:
            frames.append(normalize_motif_df(pd.read_csv(path, sep="\t")))
        elif suffix == ".csv":
            frames.append(normalize_motif_df(pd.read_csv(path)))
        else:
            raise ValueError("Motif library must be .xlsx, .xls, .tsv, .txt, or .csv")
    if include_builtin:
        frames.append(pd.DataFrame(BUILTIN_MOTIFS))
    if not frames:
        return pd.DataFrame(columns=["Family", "Motif_ID", "Consensus_Motif", "Regex", "Functional_Region", "Conservation"])
    motif_df = pd.concat(frames, ignore_index=True).drop_duplicates(subset=["Motif_ID", "Regex"])
    motif_df = motif_df.sort_values(["Family", "Motif_ID"], kind="stable")
    return motif_df


def regex_search_all(pattern: str, seq: str) -> Iterable[re.Match]:
    """Find overlapping regex matches by advancing one residue at a time."""
    try:
        compiled = re.compile(pattern, flags=re.I)
    except re.error as exc:
        warnings.warn(f"Skipping invalid regex {pattern!r}: {exc}")
        return []
    matches = []
    for i in range(len(seq)):
        m = compiled.match(seq, i)
        if m:
            matches.append(m)
    return matches


def scan_motifs(records: Sequence[SeqRecordLite], motif_df: pd.DataFrame) -> List[MotifHit]:
    hits: List[MotifHit] = []
    for rec in records:
        for _, motif in motif_df.iterrows():
            pattern = str(motif.get("Regex", "")).strip()
            if not pattern:
                continue
            for m in regex_search_all(pattern, rec.sequence):
                hits.append(MotifHit(
                    family=rec.family,
                    seq_id=rec.seq_id,
                    motif_family=str(motif.get("Family", "")),
                    motif_id=str(motif.get("Motif_ID", "")),
                    consensus=str(motif.get("Consensus_Motif", "")),
                    regex=pattern,
                    functional_region=str(motif.get("Functional_Region", "")),
                    conservation=str(motif.get("Conservation", "")),
                    start_1based=m.start() + 1,
                    end_1based=m.end(),
                    matched_sequence=m.group(0),
                ))
    return hits


# ---------------------------------------------------------------------------
# Topology and signal peptide parsing
# ---------------------------------------------------------------------------

def kd_values(seq: str) -> List[float]:
    return [KYTE_DOOLITTLE.get(a.upper(), 0.0) for a in seq]


def sliding_average(vals: Sequence[float], window: int) -> List[float]:
    if not vals:
        return []
    if window <= 1:
        return list(vals)
    if len(vals) < window:
        avg = sum(vals) / len(vals)
        return [avg] * len(vals)
    out = []
    csum = sum(vals[:window])
    out.append(csum / window)
    for i in range(window, len(vals)):
        csum += vals[i] - vals[i - window]
        out.append(csum / window)
    out.extend([out[-1]] * (len(vals) - len(out)))
    return out


def merge_ranges(ranges: Sequence[Tuple[int, int]], max_gap: int = 2) -> List[Tuple[int, int]]:
    if not ranges:
        return []
    ranges = sorted((int(s), int(e)) for s, e in ranges if s and e and e >= s)
    merged = [ranges[0]]
    for s, e in ranges[1:]:
        ps, pe = merged[-1]
        if s <= pe + max_gap:
            merged[-1] = (ps, max(pe, e))
        else:
            merged.append((s, e))
    return merged


def kd_tm_segments(seq: str, window: int = DEFAULT_KD_TM_WINDOW, threshold: float = DEFAULT_KD_TM_THRESHOLD) -> List[Tuple[int, int]]:
    vals = kd_values(seq)
    if len(vals) < window:
        return []
    raw = []
    i = 0
    while i <= len(seq) - window:
        avg = sum(vals[i:i + window]) / window
        if avg >= threshold:
            start = i + 1
            j = i + window
            while j < len(seq):
                avg2 = sum(vals[j - window + 1:j + 1]) / window
                if avg2 < threshold:
                    break
                j += 1
            raw.append((start, j))
            i = j
        else:
            i += 1
    return merge_ranges(raw, max_gap=3)


def n_terminal_hydrophobicity(seq: str, n: int = 30) -> float:
    vals = kd_values(seq[:n])
    return sum(vals) / len(vals) if vals else 0.0


def max_kd_window(seq: str, scan_len: int = DEFAULT_N_TERM_SCAN, window: int = DEFAULT_SIGNAL_WINDOW) -> Tuple[float, Optional[Tuple[int, int]]]:
    region = seq[:min(scan_len, len(seq))]
    vals = kd_values(region)
    if len(vals) < window:
        return (sum(vals) / len(vals), (1, len(vals))) if vals else (0.0, None)
    best = (-999.0, None)
    for i in range(0, len(vals) - window + 1):
        avg = sum(vals[i:i + window]) / window
        if avg > best[0]:
            best = (avg, (i + 1, i + window))
    return best


def heuristic_signal_peptide(seq: str, h_thr: float = DEFAULT_SIGNAL_H_REGION_THRESHOLD) -> Tuple[bool, Optional[int], str]:
    """Conservative fallback: N-region positive charge + hydrophobic h-region + plausible cleavage zone."""
    if len(seq) < 25:
        return False, None, "sequence_too_short"
    n_pos = sum(1 for aa in seq[:8] if aa in POSITIVE)
    h_score, h_rng = max_kd_window(seq, scan_len=45, window=14)
    if h_score < h_thr:
        return False, None, "no_hydrophobic_h_region"
    # Try a simple (-3, -1) small-residue cleavage rule between positions 16 and 35.
    for cut in range(16, min(35, len(seq) - 1)):
        minus3 = seq[cut - 3] if cut - 3 >= 0 else ""
        minus1 = seq[cut - 1] if cut - 1 >= 0 else ""
        if minus3 in "ASTGVC" and minus1 in "ASTGVC" and n_pos >= 0:
            return True, cut, f"heuristic_h_score={h_score:.2f};h_region={h_rng}"
    return True, h_rng[1] if h_rng else None, f"heuristic_h_score={h_score:.2f};no_clear_cleavage_rule"


def contiguous_char_ranges(s: str, chars: str) -> List[Tuple[int, int]]:
    ranges = []
    start = None
    char_set = set(chars)
    for i, ch in enumerate(s, start=1):
        if ch in char_set and start is None:
            start = i
        elif ch not in char_set and start is not None:
            ranges.append((start, i - 1))
            start = None
    if start is not None:
        ranges.append((start, len(s)))
    return ranges


def parse_signalp_tsv(path: Optional[Path]) -> Dict[str, TopologyCall]:
    if not path:
        return {}
    if not path.exists():
        raise FileNotFoundError(f"SignalP TSV not found: {path}")
    calls: Dict[str, TopologyCall] = {}
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                parts = re.split(r"\s+", line)
            if len(parts) < 2:
                continue
            sid = parts[0].strip()
            pred = parts[1].strip().upper()
            has_sp = pred in {"SP", "SPI", "SEC/SPI"} or pred.startswith("SP")
            cut = None
            joined = "\t".join(parts)
            m = re.search(r"CS\s*pos:\s*(\d+)\s*-\s*(\d+)", joined, flags=re.I)
            if m:
                cut = int(m.group(2))
            calls[sid] = TopologyCall(has_sp=has_sp, sp_cut=cut, signal_range=(1, cut) if cut else None, source="SignalP", label=pred)
    return calls


def parse_deeptmhmm(path: Optional[Path]) -> Dict[str, TopologyCall]:
    """Parse DeepTMHMM Markdown with 3-line topologies, predicted_topologies.3line, or GFF/GFF3."""
    if not path:
        return {}
    if not path.exists():
        raise FileNotFoundError(f"DeepTMHMM file not found: {path}")
    text = path.read_text(encoding="utf-8", errors="replace")
    calls: Dict[str, TopologyCall] = {}

    # GFF/GFF3 blocks: ID, source, feature, start, end ... or loose ID feature start end.
    for line in text.splitlines():
        if not line or line.startswith("#") or line.startswith("```") or line.startswith("##"):
            continue
        parts = line.split("\t")
        if len(parts) >= 5:
            sid = parts[0].strip()
            feature = parts[2].strip().lower()
            try:
                s, e = int(parts[3]), int(parts[4])
            except Exception:
                continue
        else:
            toks = re.split(r"\s+", line.strip())
            if len(toks) < 4:
                continue
            sid = toks[0]
            feature = toks[1].lower()
            try:
                s, e = int(toks[2]), int(toks[3])
            except Exception:
                continue
        if feature not in {"tmhelix", "transmembrane", "tm", "signal", "signal_peptide", "sp"}:
            continue
        call = calls.setdefault(sid, TopologyCall(tm_ranges=[], source="DeepTMHMM/GFF"))
        if feature in {"tmhelix", "transmembrane", "tm"}:
            call.tm_ranges = (call.tm_ranges or []) + [(s, e)]
        elif feature in {"signal", "signal_peptide", "sp"}:
            call.signal_range = (s, e)
            call.has_sp = True
            call.sp_cut = e

    # 3-line topology format embedded in Markdown or standalone.
    # Header: >ID | SP+TM ; then sequence ; then topology string with S/M/I/O.
    lines = [ln.strip() for ln in text.splitlines() if ln.strip() and ln.strip() != "```"]
    i = 0
    while i < len(lines):
        if lines[i].startswith(">"):
            header = lines[i][1:].strip()
            sid = header.split("|")[0].strip().split()[0]
            label = header.split("|", 1)[1].strip() if "|" in header else ""
            if i + 2 < len(lines):
                seq_line = re.sub(r"\s+", "", lines[i + 1])
                topo_line = re.sub(r"\s+", "", lines[i + 2])
                if seq_line and topo_line and set(topo_line) <= set("SIMOoousimM") and len(topo_line) >= max(5, int(0.5 * len(seq_line))):
                    topo = topo_line.upper()
                    tm = contiguous_char_ranges(topo, "M")
                    sig = contiguous_char_ranges(topo, "S")
                    call = calls.setdefault(sid, TopologyCall(tm_ranges=[], source="DeepTMHMM/3line"))
                    call.tm_ranges = merge_ranges((call.tm_ranges or []) + tm, max_gap=1)
                    if sig:
                        call.signal_range = sig[0]
                        call.has_sp = True
                        call.sp_cut = sig[0][1]
                    elif call.has_sp is None:
                        call.has_sp = "SP" in label.upper()
                    call.label = label
                    call.source = "DeepTMHMM/3line"
                    i += 3
                    continue
        i += 1

    for call in calls.values():
        call.tm_ranges = merge_ranges(call.tm_ranges or [], max_gap=1)
    return calls


def parse_cctop(path: Optional[Path]) -> Dict[str, TopologyCall]:
    """Flexible parser for CCTOP-like text/TSV exports. It extracts ID and TM ranges when present."""
    if not path:
        return {}
    if not path.exists():
        raise FileNotFoundError(f"CCTOP output not found: {path}")
    calls: Dict[str, TopologyCall] = {}
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row or not any(x.strip() for x in row):
                continue
            line = "\t".join(row).strip()
            if line.startswith("#"):
                continue
            # Find ranges like 15-35 or 15..35; collect all plausible TM ranges.
            ranges = [(int(a), int(b)) for a, b in re.findall(r"(\d+)\s*(?:-|\.\.)\s*(\d+)", line)]
            if not ranges:
                continue
            sid = row[0].strip().lstrip(">")
            if sid.lower() in {"id", "seqid", "sequence", "protein"} and len(row) > 1:
                continue
            calls[sid] = TopologyCall(tm_ranges=merge_ranges(ranges), source="CCTOP")
    return calls


def predict_gpi_like(seq: str) -> Tuple[bool, Optional[int], Optional[Tuple[int, int]], str]:
    if len(seq) < DEFAULT_GPI_TAIL_SCAN:
        return False, None, None, "sequence_shorter_than_tail_scan"
    tail = seq[-DEFAULT_GPI_TAIL_SCAN:]
    base = len(seq) - DEFAULT_GPI_TAIL_SCAN
    hydrophobic_tail = None
    for window in range(DEFAULT_GPI_TAIL_MIN, DEFAULT_GPI_TAIL_MAX + 1):
        for i in range(max(0, len(tail) - window - 8), len(tail) - window + 1):
            frag = tail[i:i + window]
            hyd_frac = sum(1 for aa in frag if aa in HYDRO) / window
            if hyd_frac >= DEFAULT_GPI_HYDROPHOBIC_FRACTION:
                hydrophobic_tail = (base + i + 1, base + i + window)
                break
        if hydrophobic_tail:
            break
    if not hydrophobic_tail:
        return False, None, None, "no_c_terminal_hydrophobic_tail"

    h_start = hydrophobic_tail[0]
    omega = None
    # GPI omega site is usually a small residue upstream of C-terminal hydrophobic tail.
    for pos in range(max(3, h_start - 15), max(3, h_start - 4)):
        aa = seq[pos - 1]
        aa_minus2 = seq[pos - 3]
        if aa in SMALL and aa_minus2 in SMALL:
            omega = pos
            break
    if omega is None:
        return False, None, hydrophobic_tail, "tail_present_no_omega_like_site"
    return True, omega, hydrophobic_tail, "small_residue_omega_like_site_plus_hydrophobic_tail"


def integrate_topology(rec: SeqRecordLite,
                       signalp: Dict[str, TopologyCall],
                       deeptmhmm: Dict[str, TopologyCall],
                       cctop: Dict[str, TopologyCall],
                       tm_window: int,
                       tm_threshold: float) -> Dict[str, object]:
    seq = rec.sequence
    sp_call = signalp.get(rec.seq_id)
    tm_call = deeptmhmm.get(rec.seq_id) or cctop.get(rec.seq_id)
    deeptmhmm_call = deeptmhmm.get(rec.seq_id)
    cctop_call = cctop.get(rec.seq_id)

    sp_heur, sp_heur_cut, sp_heur_reason = heuristic_signal_peptide(seq)
    kd_tm = kd_tm_segments(seq, window=tm_window, threshold=tm_threshold)
    nterm_kd = n_terminal_hydrophobicity(seq, n=30)
    max_n_kd, max_n_rng = max_kd_window(seq)

    has_sp = None
    sp_cut = None
    signal_source = "heuristic"
    if sp_call is not None and sp_call.has_sp is not None:
        has_sp = bool(sp_call.has_sp)
        sp_cut = sp_call.sp_cut
        signal_source = sp_call.source
    elif deeptmhmm_call is not None and deeptmhmm_call.has_sp is not None:
        has_sp = bool(deeptmhmm_call.has_sp)
        sp_cut = deeptmhmm_call.sp_cut
        signal_source = deeptmhmm_call.source
    else:
        has_sp = sp_heur
        sp_cut = sp_heur_cut

    if deeptmhmm_call and deeptmhmm_call.tm_ranges is not None:
        tm_ranges = deeptmhmm_call.tm_ranges
        tm_source = deeptmhmm_call.source
    elif cctop_call and cctop_call.tm_ranges is not None:
        tm_ranges = cctop_call.tm_ranges
        tm_source = cctop_call.source
    else:
        tm_ranges = kd_tm
        tm_source = "KD_fallback"

    gpi_like, gpi_omega, gpi_tail, gpi_reason = predict_gpi_like(seq)

    return {
        "SeqID": rec.seq_id,
        "Family": rec.family,
        "Length": len(seq),
        "SignalPeptide": bool(has_sp),
        "SP_Cleavage": sp_cut if sp_cut is not None else "",
        "Signal_Source": signal_source,
        "SignalP_Call": sp_call.label if sp_call else "",
        "DeepTMHMM_Label": deeptmhmm_call.label if deeptmhmm_call else "",
        "TM_Count": len(tm_ranges),
        "TM_Ranges": ";".join(f"{s}-{e}" for s, e in tm_ranges),
        "TM_Source": tm_source,
        "DeepTMHMM_TM_Ranges": ";".join(f"{s}-{e}" for s, e in (deeptmhmm_call.tm_ranges if deeptmhmm_call and deeptmhmm_call.tm_ranges else [])),
        "CCTOP_TM_Ranges": ";".join(f"{s}-{e}" for s, e in (cctop_call.tm_ranges if cctop_call and cctop_call.tm_ranges else [])),
        "KD_TM_Ranges": ";".join(f"{s}-{e}" for s, e in kd_tm),
        "Nterm_KD_30aa": round(nterm_kd, 3),
        "Max_Nterm_KD": round(max_n_kd, 3),
        "Max_Nterm_KD_Range": f"{max_n_rng[0]}-{max_n_rng[1]}" if max_n_rng else "",
        "Heuristic_SP": sp_heur,
        "Heuristic_SP_Cut": sp_heur_cut if sp_heur_cut else "",
        "Heuristic_SP_Reason": sp_heur_reason,
        "GPI_like": gpi_like,
        "GPI_omega": gpi_omega if gpi_omega else "",
        "GPI_tail_range": f"{gpi_tail[0]}-{gpi_tail[1]}" if gpi_tail else "",
        "GPI_reason": gpi_reason,
        "_tm_ranges_obj": tm_ranges,
        "_sp_cut_obj": sp_cut,
        "_gpi_tail_obj": gpi_tail,
        "_gpi_omega_obj": gpi_omega,
    }


# ---------------------------------------------------------------------------
# Hydrolase classification and feature strips
# ---------------------------------------------------------------------------

def sequence_has_motif(motif_hits: Sequence[MotifHit], seq_id: str, motif_id_contains: str) -> bool:
    target = motif_id_contains.upper()
    for h in motif_hits:
        if h.seq_id == seq_id and target in h.motif_id.upper():
            return True
    return False


def hydrolase_classification(seq: str, feature_row: Dict[str, object], motif_hits: Sequence[MotifHit]) -> str:
    seq_id = str(feature_row["SeqID"])
    has_sp = bool(feature_row["SignalPeptide"])
    tm_count = int(feature_row["TM_Count"])
    gpi_like = bool(feature_row["GPI_like"])

    gxs_hits = [h for h in motif_hits if h.seq_id == seq_id and ("GXSXG" in h.motif_id.upper() or "HYDROLASE_GXSXG" in h.motif_id.upper())]
    hgg_hits = [h for h in motif_hits if h.seq_id == seq_id and ("HGG" in h.motif_id.upper() or "OXYANION" in h.functional_region.upper())]
    has_gxs = bool(gxs_hits)
    has_hgg_like = bool(hgg_hits)
    has_hggg = any("HGGG" in h.matched_sequence.upper() for h in hgg_hits)

    if has_sp and gpi_like and tm_count <= 1:
        return "Ecto CEH (secreted / GPI-like)"
    if has_sp and tm_count == 1:
        return "Ecto CEH (type-I single-pass)"
    if has_sp and tm_count == 0:
        return "Ecto CEH (secreted)"
    if (not has_sp) and tm_count == 0:
        if has_hggg and has_gxs:
            return "Cytosolic HSL-like neutral lipase"
        if has_gxs:
            return "Cytosolic α/β-hydrolase (neutral hydrolase candidate)"
    if has_gxs or has_hgg_like:
        return "α/β-hydrolase motif-positive, topology unresolved"
    return "Unassigned"


def family_is_hydrolase(family: str) -> bool:
    fam = family.upper()
    return any(token in fam for token in ["HYDROL", "ABHYD", "ABHD", "HSL", "LIPASE", "CEH", "ESTERASE"])


def draw_feature_strip(rec: SeqRecordLite, feature_row: Dict[str, object], hits: Sequence[MotifHit], out_png: Path) -> None:
    L = len(rec.sequence)
    fig, ax = plt.subplots(figsize=(11, 1.3))
    ax.set_xlim(1, max(1, L))
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.plot([1, L], [0.5, 0.5], color="#8a8a8a", lw=7, alpha=0.35, solid_capstyle="round")

    sp_cut = feature_row.get("_sp_cut_obj")
    if feature_row.get("SignalPeptide") and sp_cut:
        ax.add_patch(plt.Rectangle((1, 0.32), int(sp_cut), 0.36, color="#1b9e77", alpha=0.75, label="Signal peptide"))

    for s, e in feature_row.get("_tm_ranges_obj") or []:
        ax.add_patch(plt.Rectangle((s, 0.24), e - s + 1, 0.52, color="#d95f02", alpha=0.65, label="TM helix"))

    gtail = feature_row.get("_gpi_tail_obj")
    gomega = feature_row.get("_gpi_omega_obj")
    if gtail:
        s, e = gtail
        ax.add_patch(plt.Rectangle((s, 0.15), e - s + 1, 0.70, fill=False, ec="#7570b3", lw=2, label="GPI-like tail"))
    if gomega:
        ax.axvline(int(gomega), ymin=0.12, ymax=0.88, color="#7570b3", lw=2)

    # Plot selected motif classes to avoid overcrowding.
    for h in hits:
        mid = h.motif_id.upper()
        if "GXSXG" in mid:
            ax.plot(h.start_1based, 0.88, marker="v", ms=6, color="#e7298a")
        elif "HGG" in mid:
            ax.plot(h.start_1based, 0.12, marker="^", ms=6, color="#66a61e")
        elif "NPC" in h.motif_family.upper() or "CD36" in h.motif_family.upper():
            ax.plot(h.start_1based, 0.78, marker="o", ms=3.5, color="#377eb8")
        elif "GLYCO" in mid:
            ax.plot(h.start_1based, 0.22, marker="|", ms=7, color="#984ea3")

    ax.text(1, 0.98, f"{rec.seq_id} | {rec.family} | {L} aa", ha="left", va="top", fontsize=8.5)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout(pad=0.2)
    fig.savefig(out_png, dpi=220, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Alignment, conservation, tree processing
# ---------------------------------------------------------------------------

def run_cmd(cmd: Sequence[str], log_path: Path, cwd: Optional[Path] = None, strict: bool = True) -> bool:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with open(log_path, "w", encoding="utf-8") as log:
        log.write("COMMAND: " + " ".join(map(str, cmd)) + "\n\n")
        try:
            subprocess.run(list(map(str, cmd)), cwd=str(cwd) if cwd else None, stdout=log, stderr=subprocess.STDOUT, check=True)
            return True
        except subprocess.CalledProcessError as exc:
            if strict:
                raise
            log.write(f"\nWARNING: command failed with exit code {exc.returncode}\n")
            return False


def write_family_fasta(records: Sequence[SeqRecordLite], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as out:
        for r in records:
            out.write(f">{r.seq_id} {r.description}\n")
            seq = r.sequence
            for i in range(0, len(seq), 70):
                out.write(seq[i:i + 70] + "\n")


def run_mafft(fasta: Path, alignment: Path, threads: int, strict: bool) -> bool:
    mafft = shutil.which("mafft")
    if not mafft:
        msg = "MAFFT not found on PATH; skipping alignment/phylogeny for this family."
        if strict:
            raise RuntimeError(msg)
        warnings.warn(msg)
        return False
    alignment.parent.mkdir(parents=True, exist_ok=True)
    cmd = [mafft, "--auto", "--thread", str(threads), str(fasta)]
    with open(alignment, "w", encoding="utf-8") as out, open(alignment.with_suffix(".mafft.log"), "w", encoding="utf-8") as log:
        log.write("COMMAND: " + " ".join(cmd) + "\n\n")
        try:
            subprocess.run(cmd, stdout=out, stderr=log, check=True)
            return True
        except subprocess.CalledProcessError:
            if strict:
                raise
            warnings.warn(f"MAFFT failed for {fasta}; see {alignment.with_suffix('.mafft.log')}")
            return False


def run_iqtree(alignment: Path, prefix: Path, threads: int, bootstrap: int, alrt: int, strict: bool) -> bool:
    iqtree = shutil.which("iqtree2") or shutil.which("iqtree")
    if not iqtree:
        msg = "IQ-TREE not found on PATH; skipping ML tree inference."
        if strict:
            raise RuntimeError(msg)
        warnings.warn(msg)
        return False
    cmd = [iqtree, "-s", str(alignment), "-m", "MFP", "-bb", str(bootstrap), "-alrt", str(alrt), "-nt", str(threads), "-pre", str(prefix)]
    return run_cmd(cmd, prefix.with_suffix(".iqtree.run.log"), strict=strict)


def frequency_matrix(alignment: Path, out_csv: Path) -> pd.DataFrame:
    aln = AlignIO.read(str(alignment), "fasta")
    rows = []
    for idx in range(aln.get_alignment_length()):
        col = [str(rec.seq[idx]).upper() for rec in aln]
        counts = Counter(aa for aa in col if aa in AA20)
        gap_count = sum(1 for aa in col if aa in {"-", "."})
        total_non_gap = sum(counts.values())
        row = {"Alignment_Position": idx + 1, "Gap_Count": gap_count, "NonGap_Count": total_non_gap}
        for aa in AA20:
            row[aa] = counts.get(aa, 0)
        if total_non_gap:
            max_freq = max(counts.values()) / total_non_gap if counts else 0
            entropy = 0.0
            for count in counts.values():
                p = count / total_non_gap
                entropy -= p * math.log2(p)
            row["Max_AA_Frequency"] = round(max_freq, 4)
            row["Shannon_Entropy"] = round(entropy, 4)
        else:
            row["Max_AA_Frequency"] = 0.0
            row["Shannon_Entropy"] = 0.0
        rows.append(row)
    df = pd.DataFrame(rows)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)
    return df


def make_sequence_logo(alignment: Path, out_pdf: Path, out_png: Optional[Path] = None) -> bool:
    if logomaker is None:
        warnings.warn("logomaker not installed; skipping sequence logo.")
        return False
    aln = AlignIO.read(str(alignment), "fasta")
    matrix = []
    for idx in range(aln.get_alignment_length()):
        counts = {aa: 0 for aa in AA20}
        for rec in aln:
            aa = str(rec.seq[idx]).upper()
            if aa in counts:
                counts[aa] += 1
        matrix.append(counts)
    df = pd.DataFrame(matrix)
    fig_width = max(10, min(28, aln.get_alignment_length() / 12))
    fig, ax = plt.subplots(figsize=(fig_width, 4.0))
    logomaker.Logo(df, ax=ax)
    ax.set_title("Sequence conservation from full-length MSA")
    ax.set_xlabel("Alignment position")
    ax.set_ylabel("Residue count")
    fig.tight_layout()
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_pdf, dpi=600)
    if out_png:
        fig.savefig(out_png, dpi=220)
    plt.close(fig)
    return True


def parse_support_from_clade(clade) -> Optional[float]:
    if getattr(clade, "confidence", None) is not None:
        try:
            return float(clade.confidence)
        except Exception:
            pass
    name = getattr(clade, "name", None)
    if not name:
        return None
    # IQ-TREE often writes SH-aLRT/UFBoot as 89.1/97.
    nums = re.findall(r"\d+(?:\.\d+)?", str(name))
    if not nums:
        return None
    try:
        return float(nums[-1])  # retain UFBoot for display/collapse by default.
    except Exception:
        return None


def process_tree(treefile: Path, out_newick: Path, out_pdf: Path, out_svg: Path,
                 support_display: float, collapse_support: float,
                 clade_keywords: Dict[str, List[str]]) -> bool:
    if not treefile.exists():
        warnings.warn(f"Tree file not found: {treefile}")
        return False
    tree = Phylo.read(str(treefile), "newick")
    try:
        tree.root_at_midpoint()
    except Exception:
        warnings.warn("Midpoint rooting failed; writing unrooted tree layout.")
    tree.ladderize()

    # Store support values, remove numerical internal node names from labels.
    for clade in tree.find_clades():
        support = parse_support_from_clade(clade)
        if support is not None:
            clade.confidence = support
            if getattr(clade, "name", None) and not clade.is_terminal():
                clade.name = None

    # Collapse low confidence branches.
    for clade in list(tree.find_clades(order="postorder")):
        if clade is tree.root or clade.is_terminal():
            continue
        support = parse_support_from_clade(clade)
        if support is not None and support < collapse_support:
            try:
                tree.collapse(clade)
            except Exception:
                pass

    # Color keyword-defined clades if possible.
    palette = ["#d62728", "#1f77b4", "#2ca02c", "#9467bd", "#ff7f0e", "#17becf", "#8c564b"]
    assigned = {}
    for idx, (label, keywords) in enumerate(clade_keywords.items()):
        found = find_best_keyword_clade(tree, keywords)
        if found:
            color_clade(found, palette[idx % len(palette)])
            assigned[label] = found

    def label_func(clade):
        if clade.is_terminal():
            return clade.name
        support = parse_support_from_clade(clade)
        if support is not None and support >= support_display:
            return f"{support:.0f}"
        return None

    n_terms = len(tree.get_terminals())
    fig_h = max(6, min(40, 0.22 * n_terms + 2.5))
    fig, ax = plt.subplots(figsize=(12, fig_h))
    Phylo.draw(tree, axes=ax, do_show=False, label_func=label_func)
    ax.set_title("Maximum likelihood phylogeny")
    ax.set_xlabel("Substitutions per site")
    ax.set_ylabel("")
    annotate_keyword_clades(ax, assigned)
    for line in ax.get_lines():
        line.set_linewidth(1.2)
    fig.tight_layout()
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_pdf, dpi=600)
    fig.savefig(out_svg)
    plt.close(fig)
    Phylo.write(tree, str(out_newick), "newick")
    return True


def find_best_keyword_clade(tree, keywords: Sequence[str]):
    keywords = [k for k in keywords if k]
    if not keywords:
        return None
    terminals = tree.get_terminals()
    matching = [t for t in terminals if any(k.lower() in (t.name or "").lower() for k in keywords)]
    if not matching:
        return None
    if len(matching) == 1:
        return matching[0]
    try:
        return tree.common_ancestor(matching)
    except Exception:
        return None


def color_clade(clade, color: str) -> None:
    clade.color = color
    for child in clade.clades:
        color_clade(child, color)


def annotate_keyword_clades(ax, assigned: Dict[str, object]) -> None:
    label_positions = {text.get_text(): text.get_position() for text in ax.texts}
    for idx, (label, clade) in enumerate(assigned.items()):
        leaves = [t.name for t in clade.get_terminals()]
        ys = [label_positions[name][1] for name in leaves if name in label_positions]
        if not ys:
            continue
        y_mid = sum(ys) / len(ys)
        ax.text(-0.02, y_mid, label, fontsize=9, fontweight="bold", ha="right", va="center", transform=ax.get_yaxis_transform())


def parse_clade_keywords(text: str) -> Dict[str, List[str]]:
    if not text:
        return {
            "HSL/LIPE-like": ["LIPS", "LIPE", "HSL"],
            "RBBP9-like": ["RBBP9"],
            "CES/EST-like": ["CES", "EST", "CEH"],
            "NPC1-like": ["NPC1"],
            "NPC2-like": ["NPC2"],
            "CD36RP-like": ["CD36", "SRB", "SCARB"],
        }
    out: Dict[str, List[str]] = {}
    for block in text.split(";"):
        if not block.strip():
            continue
        if "=" in block:
            label, keys = block.split("=", 1)
            out[label.strip()] = [k.strip() for k in re.split(r"[,|]", keys) if k.strip()]
        else:
            out[block.strip()] = [block.strip()]
    return out


# ---------------------------------------------------------------------------
# Reports and outputs
# ---------------------------------------------------------------------------

def build_pdf_report(pdf_path: Path, feature_df: pd.DataFrame, strip_paths: Sequence[Path], title: str) -> None:
    if SimpleDocTemplate is None:
        warnings.warn("reportlab not installed; skipping PDF report.")
        return
    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    doc = SimpleDocTemplate(str(pdf_path), pagesize=landscape(letter), rightMargin=18, leftMargin=18, topMargin=18, bottomMargin=18)
    styles = getSampleStyleSheet()
    story = [Paragraph(f"<b>{html.escape(title)}</b>", styles["Title"]), Spacer(1, 8)]
    story.append(Paragraph(
        "Integrated evidence: SignalP/heuristic signal peptide, DeepTMHMM/CCTOP/KD topology, "
        "GPI-like tail heuristic, curated motifs, and hydrolase classification.",
        styles["Normal"],
    ))
    story.append(Spacer(1, 10))
    cols = ["SeqID", "Family", "Length", "SignalPeptide", "TM_Count", "GPI_like", "GxSxG_count", "HGG_like_count", "Hydrolase_Classification"]
    display = feature_df[[c for c in cols if c in feature_df.columns]].copy().head(80)
    data = [display.columns.tolist()] + display.astype(str).values.tolist()
    table = Table(data, hAlign="LEFT", repeatRows=1)
    table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#eef2ff")),
        ("GRID", (0, 0), (-1, -1), 0.25, colors.gray),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTSIZE", (0, 0), (-1, -1), 6.7),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
    ]))
    story.append(table)
    story.append(Spacer(1, 12))
    for path in strip_paths[:120]:
        if path.exists():
            story.append(RLImage(str(path), width=720, height=82))
            story.append(Spacer(1, 5))
    doc.build(story)


def write_excel(out_xlsx: Path, feature_df: pd.DataFrame, motif_df: pd.DataFrame, hits_df: pd.DataFrame, family_summary: pd.DataFrame) -> None:
    out_xlsx.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(out_xlsx, engine="openpyxl") as writer:
        feature_df.to_excel(writer, sheet_name="SEQUENCE_FEATURES", index=False)
        hits_df.to_excel(writer, sheet_name="MOTIF_HITS", index=False)
        family_summary.to_excel(writer, sheet_name="FAMILY_SUMMARY", index=False)
        motif_df.to_excel(writer, sheet_name="MOTIF_LIBRARY", index=False)


def summarize_by_family(feature_df: pd.DataFrame, hits_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for fam, sub in feature_df.groupby("Family", dropna=False):
        fam_hits = hits_df[hits_df["family"] == fam] if not hits_df.empty else pd.DataFrame()
        rows.append({
            "Family": fam,
            "Sequences": len(sub),
            "Mean_Length": round(sub["Length"].mean(), 2) if len(sub) else 0,
            "SignalPeptide_n": int(sub["SignalPeptide"].sum()) if "SignalPeptide" in sub else 0,
            "TM_positive_n": int((sub["TM_Count"] > 0).sum()) if "TM_Count" in sub else 0,
            "GPI_like_n": int(sub["GPI_like"].sum()) if "GPI_like" in sub else 0,
            "Motif_hits": len(fam_hits),
            "Unique_motif_positive_sequences": fam_hits["seq_id"].nunique() if not fam_hits.empty else 0,
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Main workflow
# ---------------------------------------------------------------------------

def process_phylogeny_for_family(family: str, records: Sequence[SeqRecordLite], outdir: Path, args: argparse.Namespace) -> Dict[str, str]:
    fam_dir = outdir / "phylogeny" / sanitize_name(family)
    fam_fasta = fam_dir / f"{sanitize_name(family)}.fasta"
    alignment = fam_dir / "alignment.mafft.fasta"
    tree_prefix = fam_dir / "iqtree_ml"
    outputs = {
        "family": family,
        "fasta": str(fam_fasta),
        "alignment": "",
        "treefile": "",
        "processed_tree": "",
        "tree_pdf": "",
        "tree_svg": "",
        "logo_pdf": "",
        "frequency_csv": "",
        "status": "not_run",
    }
    if len(records) < 3:
        outputs["status"] = "skipped_less_than_3_sequences"
        return outputs
    write_family_fasta(records, fam_fasta)
    ok_align = run_mafft(fam_fasta, alignment, args.threads, strict=args.strict_external)
    if not ok_align:
        outputs["status"] = "skipped_mafft_unavailable_or_failed"
        return outputs
    outputs["alignment"] = str(alignment)
    freq_csv = fam_dir / "alignment_frequency_matrix.csv"
    frequency_matrix(alignment, freq_csv)
    outputs["frequency_csv"] = str(freq_csv)
    logo_pdf = fam_dir / "sequence_logo.pdf"
    logo_png = fam_dir / "sequence_logo.png"
    if make_sequence_logo(alignment, logo_pdf, logo_png):
        outputs["logo_pdf"] = str(logo_pdf)
    ok_tree = run_iqtree(alignment, tree_prefix, args.threads, args.bootstrap, args.alrt, strict=args.strict_external)
    if not ok_tree:
        outputs["status"] = "alignment_done_iqtree_unavailable_or_failed"
        return outputs
    treefile = Path(str(tree_prefix) + ".treefile")
    outputs["treefile"] = str(treefile)
    out_newick = fam_dir / "phylogeny_midpoint_ladderized_collapsed.nwk"
    out_pdf = fam_dir / "phylogeny_final.pdf"
    out_svg = fam_dir / "phylogeny_final.svg"
    process_tree(treefile, out_newick, out_pdf, out_svg, args.support_display, args.collapse_support, parse_clade_keywords(args.clade_keywords))
    outputs.update({"processed_tree": str(out_newick), "tree_pdf": str(out_pdf), "tree_svg": str(out_svg), "status": "complete"})
    return outputs


def run_pipeline(args: argparse.Namespace) -> None:
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "tables").mkdir(exist_ok=True)
    (outdir / "strips").mkdir(exist_ok=True)

    fasta_paths = collect_fasta_paths(args)
    all_records: List[SeqRecordLite] = []
    for fp in fasta_paths:
        all_records.extend(parse_fasta(fp))
    if not all_records:
        raise RuntimeError("No valid protein records found in the input FASTA files.")

    motif_df = load_motif_library(Path(args.motif_library) if args.motif_library else None, include_builtin=not args.no_builtin_motifs)
    motif_df.to_csv(outdir / "tables" / "motif_library_normalized.tsv", sep="\t", index=False)

    signalp_calls = parse_signalp_tsv(Path(args.signalp) if args.signalp else None)
    deeptmhmm_calls = parse_deeptmhmm(Path(args.deeptmhmm) if args.deeptmhmm else None)
    cctop_calls = parse_cctop(Path(args.cctop) if args.cctop else None)

    motif_hits = scan_motifs(all_records, motif_df)
    hits_df = pd.DataFrame([asdict(h) for h in motif_hits])
    if hits_df.empty:
        hits_df = pd.DataFrame(columns=[f.name for f in MotifHit.__dataclass_fields__.values()])
    hits_df.to_csv(outdir / "tables" / "all_motif_hits.tsv", sep="\t", index=False)

    hit_map: Dict[str, List[MotifHit]] = defaultdict(list)
    for h in motif_hits:
        hit_map[h.seq_id].append(h)

    feature_rows = []
    strip_paths: List[Path] = []
    for rec in all_records:
        row = integrate_topology(rec, signalp_calls, deeptmhmm_calls, cctop_calls, args.tm_window, args.tm_threshold)
        seq_hits = hit_map.get(rec.seq_id, [])
        row["Motif_Hit_Count"] = len(seq_hits)
        row["Motif_IDs"] = ";".join(sorted({h.motif_id for h in seq_hits}))
        row["GxSxG_count"] = sum(1 for h in seq_hits if "GXSXG" in h.motif_id.upper())
        row["HGG_like_count"] = sum(1 for h in seq_hits if "HGG" in h.motif_id.upper() or "OXYANION" in h.functional_region.upper())
        row["N_glyco_count"] = sum(1 for h in seq_hits if "GLYCO" in h.motif_id.upper())
        row["Hydrolase_Classification"] = hydrolase_classification(rec.sequence, row, seq_hits) if family_is_hydrolase(rec.family) or row["GxSxG_count"] or row["HGG_like_count"] else "Not hydrolase-family candidate"
        row["Description"] = rec.description
        if args.plot_strips:
            strip_path = outdir / "strips" / sanitize_name(rec.family) / f"{sanitize_name(rec.seq_id)}.png"
            draw_feature_strip(rec, row, seq_hits, strip_path)
            row["StripPNG"] = str(strip_path)
            strip_paths.append(strip_path)
        # Drop object columns before tabular export.
        for key in ["_tm_ranges_obj", "_sp_cut_obj", "_gpi_tail_obj", "_gpi_omega_obj"]:
            row.pop(key, None)
        feature_rows.append(row)

    feature_df = pd.DataFrame(feature_rows)
    feature_df.to_csv(outdir / "tables" / "all_sequence_features.tsv", sep="\t", index=False)

    family_summary = summarize_by_family(feature_df, hits_df)
    family_summary.to_csv(outdir / "tables" / "family_summary.tsv", sep="\t", index=False)

    if args.excel:
        write_excel(outdir / "schisto_phylo_motif_results.xlsx", feature_df, motif_df, hits_df, family_summary)

    phylo_rows = []
    if args.run_phylogeny:
        by_family: Dict[str, List[SeqRecordLite]] = defaultdict(list)
        for rec in all_records:
            by_family[rec.family].append(rec)
        for fam, recs in sorted(by_family.items()):
            phylo_rows.append(process_phylogeny_for_family(fam, recs, outdir, args))
        pd.DataFrame(phylo_rows).to_csv(outdir / "tables" / "phylogeny_outputs.tsv", sep="\t", index=False)

    if args.pdf_report:
        build_pdf_report(outdir / "integrated_evidence_report.pdf", feature_df, strip_paths, "Schistosome lipid-protein phylogeny and motif evidence report")

    print("Integrated pipeline complete.")
    print(f"Sequences processed: {len(all_records)}")
    print(f"Motif hits: {len(hits_df)}")
    print(f"Results directory: {outdir.resolve()}")


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Integrated phylogeny, motif annotation, topology parsing, and hydrolase classification pipeline for Schistosoma lipid-metabolism candidates.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--fasta", nargs="*", help="One or more family-specific FASTA files.")
    p.add_argument("--input-dir", help="Directory containing family-specific .fa/.faa/.fasta files.")
    p.add_argument("--motif-library", help="Curated motif library in XLSX/CSV/TSV format. Built-in hydrolase motifs are always added unless --no-builtin-motifs is used.")
    p.add_argument("--signalp", help="SignalP-6.0 TSV output, e.g. prediction_results.txt.")
    p.add_argument("--deeptmhmm", help="DeepTMHMM Markdown, GFF3, or predicted_topologies.3line output.")
    p.add_argument("--cctop", help="Optional CCTOP TSV/text output with sequence IDs and TM ranges.")
    p.add_argument("--outdir", default="schisto_phylo_motif_results", help="Output directory.")
    p.add_argument("--threads", type=int, default=4, help="Threads for MAFFT/IQ-TREE.")
    p.add_argument("--run-phylogeny", action="store_true", help="Run MAFFT and IQ-TREE for each family FASTA.")
    p.add_argument("--bootstrap", type=int, default=1000, help="IQ-TREE UFBoot replicates.")
    p.add_argument("--alrt", type=int, default=1000, help="IQ-TREE SH-aLRT replicates.")
    p.add_argument("--support-display", type=float, default=DEFAULT_SUPPORT_DISPLAY, help="Minimum support value displayed on tree plots.")
    p.add_argument("--collapse-support", type=float, default=DEFAULT_COLLAPSE_SUPPORT, help="Collapse internal branches below this support value.")
    p.add_argument("--clade-keywords", default="", help="Optional clade labels and terminal-keywords: 'Label=key1,key2;Other=key3'.")
    p.add_argument("--strict-external", action="store_true", help="Fail instead of warning when MAFFT/IQ-TREE are unavailable or fail.")
    p.add_argument("--tm-window", type=int, default=DEFAULT_KD_TM_WINDOW, help="Kyte-Doolittle TM fallback window.")
    p.add_argument("--tm-threshold", type=float, default=DEFAULT_KD_TM_THRESHOLD, help="Kyte-Doolittle mean threshold for TM fallback.")
    p.add_argument("--plot-strips", action=argparse.BooleanOptionalAction, default=True, help="Generate per-sequence strip diagrams.")
    p.add_argument("--excel", action=argparse.BooleanOptionalAction, default=True, help="Write combined Excel workbook.")
    p.add_argument("--pdf-report", action="store_true", help="Write PDF evidence report containing summary table and strip diagrams.")
    p.add_argument("--no-builtin-motifs", action="store_true", help="Do not add built-in GxSxG/HGG/N-glyco motifs.")
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_arg_parser().parse_args(argv)
    try:
        run_pipeline(args)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
