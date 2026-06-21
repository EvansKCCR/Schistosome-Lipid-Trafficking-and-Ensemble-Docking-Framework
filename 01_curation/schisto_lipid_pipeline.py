#!/usr/bin/env python3
"""
Schistosome lipid metabolism-associated protein discovery pipeline.

This pipeline screens schistosome proteomes for Pfam domains associated with
lipid binding, transport, and enzymatic processing. It can run hmmscan, parse
HMMER domtblout files, extract candidate protein sequences, and write
family-specific FASTA files plus Excel workbooks.

Author: Adapted from species-specific Sh/Sj/Sm lipid-gene extraction scripts.
"""

from __future__ import annotations

import argparse
import csv
import os
import shutil
import subprocess
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


DEFAULT_DOMAINS: Dict[str, List[str]] = {
    "FABP": ["PF00061"],
    "CD36": ["PF01130"],
    "ABC": ["PF00005"],
    "NPC1": ["PF12349"],
    "NPC2": ["PF02221"],
    "FATP": ["PF02259"],
    "HSL_N": ["PF06350"],
    "Abhydrolase_3": ["PF07859"],
}


@dataclass(frozen=True)
class Sample:
    """A proteome sample/species to process."""

    species: str
    fasta: Path
    domtblout: Optional[Path] = None


def strip_pfam_version(pfam_acc: str) -> str:
    """Return PF accession without release-version suffix, e.g. PF00061.30 -> PF00061."""
    return pfam_acc.split(".", 1)[0]


def safe_name(value: str) -> str:
    """Return a file-system friendly sample/family name."""
    return (
        value.strip()
        .replace(" ", "_")
        .replace("/", "_")
        .replace("\\", "_")
        .replace(":", "_")
    )


def read_samples(samples_tsv: Path) -> List[Sample]:
    """
    Read a tab-separated sample sheet.

    Required columns:
      - species
      - fasta

    Optional column:
      - domtblout
    """
    samples: List[Sample] = []
    with samples_tsv.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"species", "fasta"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Sample sheet is missing required column(s): {', '.join(sorted(missing))}")

        base_dir = samples_tsv.parent
        for row_number, row in enumerate(reader, start=2):
            species = (row.get("species") or "").strip()
            fasta = (row.get("fasta") or "").strip()
            domtblout = (row.get("domtblout") or "").strip()

            if not species or not fasta:
                raise ValueError(f"Row {row_number} must contain non-empty species and fasta values.")

            fasta_path = Path(fasta)
            if not fasta_path.is_absolute():
                fasta_path = base_dir / fasta_path

            domtblout_path: Optional[Path] = None
            if domtblout:
                domtblout_path = Path(domtblout)
                if not domtblout_path.is_absolute():
                    domtblout_path = base_dir / domtblout_path

            samples.append(Sample(species=species, fasta=fasta_path, domtblout=domtblout_path))

    if not samples:
        raise ValueError("Sample sheet contains no samples.")
    return samples


def read_domain_config(domain_tsv: Optional[Path]) -> Dict[str, List[str]]:
    """
    Read domain-family definitions from a two-column TSV: Family<TAB>Pfam.

    Pfam values may be versioned or unversioned. They are normalized internally
    to unversioned PF accessions to make the pipeline robust to Pfam release changes.
    """
    if domain_tsv is None:
        return DEFAULT_DOMAINS

    domains: Dict[str, List[str]] = defaultdict(list)
    with domain_tsv.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"Family", "Pfam"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Domain TSV is missing required column(s): {', '.join(sorted(missing))}")
        for row_number, row in enumerate(reader, start=2):
            family = (row.get("Family") or "").strip()
            pfam = (row.get("Pfam") or "").strip()
            if not family or not pfam:
                raise ValueError(f"Row {row_number} in domain TSV has an empty Family or Pfam field.")
            domains[family].append(strip_pfam_version(pfam))

    if not domains:
        raise ValueError("Domain TSV contains no domain definitions.")
    return dict(domains)


def index_pfam_to_family(domains: Mapping[str, Sequence[str]]) -> Dict[str, List[str]]:
    """Create PF accession -> family/families lookup."""
    lookup: Dict[str, List[str]] = defaultdict(list)
    for family, pfams in domains.items():
        for pfam in pfams:
            normalized = strip_pfam_version(pfam)
            if family not in lookup[normalized]:
                lookup[normalized].append(family)
    return dict(lookup)


def check_inputs(samples: Sequence[Sample], pfam_hmm: Optional[Path], run_hmmscan: bool) -> None:
    """Validate expected input files and executable availability."""
    for sample in samples:
        if not sample.fasta.exists():
            raise FileNotFoundError(f"FASTA file not found for {sample.species}: {sample.fasta}")
        if not run_hmmscan:
            if sample.domtblout is None:
                raise ValueError(
                    f"No domtblout was supplied for {sample.species}. Provide domtblout in the sample sheet "
                    "or run with --run-hmmscan and --pfam-hmm."
                )
            if not sample.domtblout.exists():
                raise FileNotFoundError(f"domtblout file not found for {sample.species}: {sample.domtblout}")

    if run_hmmscan:
        if pfam_hmm is None:
            raise ValueError("--pfam-hmm is required when --run-hmmscan is used.")
        if not pfam_hmm.exists():
            raise FileNotFoundError(f"Pfam HMM database not found: {pfam_hmm}")
        if shutil.which("hmmscan") is None:
            raise EnvironmentError("hmmscan was not found on PATH. Install HMMER or activate the correct environment.")


def run_hmmscan(sample: Sample, pfam_hmm: Path, domtblout: Path, log_file: Path, cpu: int) -> None:
    """Run HMMER hmmscan for one proteome."""
    domtblout.parent.mkdir(parents=True, exist_ok=True)
    log_file.parent.mkdir(parents=True, exist_ok=True)

    command = [
        "hmmscan",
        "--cpu",
        str(cpu),
        "--domtblout",
        str(domtblout),
        str(pfam_hmm),
        str(sample.fasta),
    ]

    with log_file.open("w") as log_handle:
        result = subprocess.run(command, stdout=log_handle, stderr=subprocess.STDOUT, text=True)

    if result.returncode != 0:
        raise RuntimeError(f"hmmscan failed for {sample.species}. See log: {log_file}")


def parse_domtblout(
    domtblout: Path,
    pfam_to_family: Mapping[str, Sequence[str]],
    evalue_cutoff: float,
) -> List[Dict[str, object]]:
    """
    Parse HMMER --domtblout and return lipid-domain hit records.

    For hmmscan domtblout, the relevant columns are:
      parts[1]  = target/Pfam accession
      parts[3]  = query/protein identifier
      parts[12] = independent domain E-value
    """
    records: List[Dict[str, object]] = []

    with domtblout.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith("#") or not line.strip():
                continue

            parts = line.rstrip().split()
            if len(parts) < 13:
                raise ValueError(f"Malformed domtblout line {line_number} in {domtblout}")

            raw_pfam_acc = parts[1]
            pfam_acc = strip_pfam_version(raw_pfam_acc)
            protein_id = parts[3]

            try:
                i_evalue = float(parts[12])
            except ValueError as exc:
                raise ValueError(
                    f"Could not parse independent E-value at line {line_number} in {domtblout}: {parts[12]}"
                ) from exc

            if i_evalue <= evalue_cutoff and pfam_acc in pfam_to_family:
                for family in pfam_to_family[pfam_acc]:
                    records.append(
                        {
                            "Protein_ID": protein_id,
                            "Family": family,
                            "Pfam": pfam_acc,
                            "Pfam_with_version": raw_pfam_acc,
                            "Independent_Evalue": i_evalue,
                        }
                    )

    return records


def load_sequences(fasta: Path) -> Dict[str, SeqRecord]:
    """Load FASTA records by primary sequence ID."""
    seqs = {record.id: record for record in SeqIO.parse(str(fasta), "fasta")}
    if not seqs:
        raise ValueError(f"No sequences were loaded from FASTA: {fasta}")
    return seqs


def annotate_records(records: List[Dict[str, object]], seqs: Mapping[str, SeqRecord], species: str) -> pd.DataFrame:
    """Add species and sequence-length metadata to parsed records."""
    annotated: List[Dict[str, object]] = []
    missing_sequences: List[str] = []

    for record in records:
        protein_id = str(record["Protein_ID"])
        seq = seqs.get(protein_id)
        if seq is None:
            missing_sequences.append(protein_id)
            continue
        annotated.append(
            {
                "Species": species,
                "Protein_ID": protein_id,
                "Family": record["Family"],
                "Pfam": record["Pfam"],
                "Pfam_with_version": record["Pfam_with_version"],
                "Independent_Evalue": record["Independent_Evalue"],
                "Length": len(seq.seq),
            }
        )

    if missing_sequences:
        unique_missing = sorted(set(missing_sequences))
        example = ", ".join(unique_missing[:5])
        raise KeyError(
            f"{len(unique_missing)} protein IDs in domtblout were not found in FASTA for {species}. "
            f"Examples: {example}"
        )

    df = pd.DataFrame(annotated)
    if df.empty:
        return df

    return df.drop_duplicates().sort_values(
        ["Family", "Protein_ID", "Independent_Evalue", "Pfam"], ignore_index=True
    )


def write_family_fastas(df: pd.DataFrame, seqs: Mapping[str, SeqRecord], outdir: Path) -> None:
    """Write one deduplicated FASTA file per lipid-associated family."""
    fasta_dir = outdir / "FASTA"
    fasta_dir.mkdir(parents=True, exist_ok=True)

    for family in sorted(df["Family"].unique()):
        family_df = df[df["Family"] == family]
        protein_ids = sorted(family_df["Protein_ID"].unique())
        output_fasta = fasta_dir / f"{safe_name(family)}.fasta"
        with output_fasta.open("w") as handle:
            for protein_id in protein_ids:
                SeqIO.write(seqs[protein_id], handle, "fasta")


def write_excel_workbook(df: pd.DataFrame, output_xlsx: Path) -> None:
    """Write per-family worksheets and ALL_CANDIDATES worksheet."""
    output_xlsx.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(output_xlsx, engine="openpyxl") as writer:
        for family in sorted(df["Family"].unique()):
            sheet = safe_name(family)[:31]
            df[df["Family"] == family].to_excel(writer, sheet_name=sheet, index=False)
        df.to_excel(writer, sheet_name="ALL_CANDIDATES", index=False)


def summarize(df: pd.DataFrame, output_csv: Path) -> pd.DataFrame:
    """Write a compact count table by species and family."""
    summary = (
        df.groupby(["Species", "Family"], as_index=False)["Protein_ID"]
        .nunique()
        .rename(columns={"Protein_ID": "Candidate_Count"})
        .sort_values(["Species", "Family"], ignore_index=True)
    )
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(output_csv, index=False)
    return summary


def process_sample(
    sample: Sample,
    outdir: Path,
    pfam_to_family: Mapping[str, Sequence[str]],
    evalue_cutoff: float,
    run_scan: bool,
    pfam_hmm: Optional[Path],
    cpu: int,
) -> pd.DataFrame:
    """Run/parse one sample and write species-level outputs."""
    species_safe = safe_name(sample.species)
    sample_outdir = outdir / species_safe
    hmmer_dir = sample_outdir / "hmmer"

    domtblout = sample.domtblout
    if run_scan:
        assert pfam_hmm is not None
        domtblout = hmmer_dir / f"{species_safe}.domtblout.txt"
        log_file = hmmer_dir / f"{species_safe}.hmmscan.log"
        print(f"[hmmscan] {sample.species}")
        run_hmmscan(sample, pfam_hmm, domtblout, log_file, cpu)

    assert domtblout is not None
    print(f"[parse] {sample.species}: {domtblout}")
    seqs = load_sequences(sample.fasta)
    records = parse_domtblout(domtblout, pfam_to_family, evalue_cutoff)
    df = annotate_records(records, seqs, sample.species)

    if df.empty:
        raise RuntimeError(
            f"No lipid Pfam hits detected for {sample.species}. "
            "Check Pfam IDs, HMM file/index, input FASTA, domtblout format, or E-value cutoff."
        )

    write_family_fastas(df, seqs, sample_outdir)
    write_excel_workbook(df, sample_outdir / f"{species_safe}_Lipid_Gene_Candidates.xlsx")
    summarize(df, sample_outdir / f"{species_safe}_candidate_counts_by_family.csv")
    return df


def write_combined_outputs(all_df: pd.DataFrame, outdir: Path) -> None:
    """Write combined multi-species outputs."""
    combined_dir = outdir / "combined"
    combined_dir.mkdir(parents=True, exist_ok=True)

    all_df = all_df.sort_values(
        ["Species", "Family", "Protein_ID", "Independent_Evalue", "Pfam"], ignore_index=True
    )
    all_df.to_csv(combined_dir / "ALL_SPECIES_Lipid_Gene_Candidates.csv", index=False)
    write_excel_workbook(all_df, combined_dir / "ALL_SPECIES_Lipid_Gene_Candidates.xlsx")
    summarize(all_df, combined_dir / "ALL_SPECIES_candidate_counts_by_family.csv")


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Identify lipid metabolism-associated proteins in schistosome proteomes using Pfam hmmscan results."
    )
    parser.add_argument(
        "--samples",
        required=True,
        type=Path,
        help="Tab-separated sample sheet with columns: species, fasta, and optional domtblout.",
    )
    parser.add_argument(
        "--pfam-hmm",
        type=Path,
        default=None,
        help="Path to Pfam-A.hmm. Required with --run-hmmscan.",
    )
    parser.add_argument(
        "--domains",
        type=Path,
        default=None,
        help="Optional TSV with columns Family and Pfam. Defaults to lipid-associated Pfam families.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("lipid_pipeline_results"),
        help="Output directory. Default: lipid_pipeline_results",
    )
    parser.add_argument(
        "--evalue",
        type=float,
        default=1e-5,
        help="Independent domain E-value cutoff. Default: 1e-5",
    )
    parser.add_argument(
        "--cpu",
        type=int,
        default=8,
        help="Number of CPUs for hmmscan. Default: 8",
    )
    parser.add_argument(
        "--run-hmmscan",
        action="store_true",
        help="Run hmmscan before extraction. If omitted, existing domtblout files from the sample sheet are parsed.",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_arg_parser().parse_args(argv)

    samples = read_samples(args.samples)
    domains = read_domain_config(args.domains)
    pfam_to_family = index_pfam_to_family(domains)

    check_inputs(samples, args.pfam_hmm, args.run_hmmscan)
    args.outdir.mkdir(parents=True, exist_ok=True)

    all_sample_frames: List[pd.DataFrame] = []
    for sample in samples:
        sample_df = process_sample(
            sample=sample,
            outdir=args.outdir,
            pfam_to_family=pfam_to_family,
            evalue_cutoff=args.evalue,
            run_scan=args.run_hmmscan,
            pfam_hmm=args.pfam_hmm,
            cpu=args.cpu,
        )
        all_sample_frames.append(sample_df)

    all_df = pd.concat(all_sample_frames, ignore_index=True).drop_duplicates()
    if all_df.empty:
        raise RuntimeError(
            "No lipid Pfam hits detected in any sample. Check domain definitions, HMM database, "
            "threshold, FASTA inputs, and domtblout files."
        )

    write_combined_outputs(all_df, args.outdir)
    print(f"\n✅ Lipid-associated protein discovery complete. Results: {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
