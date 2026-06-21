#!/usr/bin/env bash
set -euo pipefail

# Example 1: generate refinement scripts for input receptor PDBs.
python receptor_docking_prep_pipeline.py \
  --steps refine \
  --mode general \
  --input-dir input \
  --mdp-dir mdp \
  --outdir dockprep_results

# Example 2: generate CEH-like MD validation scripts with catalytic-distance checks.
python receptor_docking_prep_pipeline.py \
  --steps md-validation \
  --mode ceh \
  --md-root md_runs \
  --ceh-ser 207 \
  --ceh-his 508 \
  --ceh-acid 388 \
  --ceh-alt-acid 206 \
  --outdir ceh_dockprep_results

# Example 3: summarize/plot existing XVG files.
python receptor_docking_prep_pipeline.py \
  --steps plot \
  --xvg-root ceh_dockprep_results/md_analysis \
  --outdir ceh_dockprep_results

# Example 4: CEH-specific representative receptor selection.
python receptor_docking_prep_pipeline.py \
  --steps select \
  --mode ceh \
  --metrics-table examples/example_metrics_ceh.tsv \
  --top-n 3 \
  --outdir ceh_dockprep_results
