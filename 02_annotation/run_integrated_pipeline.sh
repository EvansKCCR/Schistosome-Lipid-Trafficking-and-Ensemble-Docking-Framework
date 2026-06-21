#!/usr/bin/env bash
set -euo pipefail

# Edit these paths to match your local directory layout.
FASTA_DIR="FASTA"
MOTIF_LIBRARY="Motif_library.xlsx"
SIGNALP_TSV="prediction_results.txt"
DEEPTMHMM_MD="deeptmhmm_results.md"
OUTDIR="results_phylo_motif"
THREADS=8

python schisto_phylo_motif_pipeline.py \
  --input-dir "${FASTA_DIR}" \
  --motif-library "${MOTIF_LIBRARY}" \
  --signalp "${SIGNALP_TSV}" \
  --deeptmhmm "${DEEPTMHMM_MD}" \
  --outdir "${OUTDIR}" \
  --threads "${THREADS}" \
  --run-phylogeny \
  --bootstrap 1000 \
  --alrt 1000 \
  --pdf-report
