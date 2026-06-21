#!/usr/bin/env bash
set -euo pipefail

# General CD36RP-like workbook analysis
python post_md_complex_analysis_pipeline.py workbook \
  --input trajectory_data.xlsx \
  --outdir post_MD_general \
  --mode general \
  --formats png,svg \
  --dpi 600

# CEH-like workbook analysis
python post_md_complex_analysis_pipeline.py workbook \
  --input CEH_trajectory_data.xlsx \
  --outdir post_MD_CEH \
  --mode ceh \
  --hbond-cutoff 0.35 \
  --attack-cutoff 0.45 \
  --oxyanion-cutoff 0.35 \
  --formats png,svg \
  --dpi 600
