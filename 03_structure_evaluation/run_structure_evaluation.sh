#!/usr/bin/env bash
set -euo pipefail

# Example run for AF3/Boltz/Boltz-2/generic outputs.
# Replace "predictions" with the folder containing .cif/.mmcif/.pdb/.json/.npz files.
python deeplearning_structure_evaluation_pipeline.py \
  predictions \
  --outdir dl_structure_eval_results \
  --motifs HGG,GESAG \
  --motif-mode exact \
  --write-excel \
  --convert-pdb
