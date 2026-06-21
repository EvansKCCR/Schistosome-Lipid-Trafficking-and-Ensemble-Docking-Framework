#!/usr/bin/env bash
set -euo pipefail

# Example: run hmmscan and extract lipid-associated candidates for all samples.
python schisto_lipid_pipeline.py \
  --samples samples.tsv \
  --pfam-hmm Pfam-A.hmm \
  --domains lipid_domains.tsv \
  --outdir lipid_pipeline_results \
  --evalue 1e-5 \
  --cpu 8 \
  --run-hmmscan
