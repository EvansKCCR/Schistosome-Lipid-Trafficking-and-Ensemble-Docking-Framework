#!/usr/bin/env bash
set -euo pipefail
python deeplearning_structure_evaluation_pipeline.py \
  examples/demo_input \
  --outdir demo_results \
  --motifs HGG,GESAG \
  --write-excel
