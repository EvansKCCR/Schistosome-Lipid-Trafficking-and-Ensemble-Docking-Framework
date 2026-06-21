#!/usr/bin/env bash
set -euo pipefail

# Example 1: hydrolase / CEH representative model selection
python protein_model_validation_selection_pipeline.py \
  --metrics-table CEH_combined_model_selection.csv \
  --engine hydrolase \
  --family "alpha/beta hydrolase" \
  --outdir ceh_representative_selection \
  --top-n 10

# Example 2: general structural-validation model selection
python protein_model_validation_selection_pipeline.py \
  --metrics-table quality_matrices.xlsx \
  --sheet "Table S6" \
  --engine general \
  --outdir general_model_selection \
  --top-n 10
