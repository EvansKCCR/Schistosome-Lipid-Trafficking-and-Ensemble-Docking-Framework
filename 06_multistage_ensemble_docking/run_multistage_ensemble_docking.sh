#!/bin/bash
set -euo pipefail

# 1. Cluster DiffDock-L / SMINA pose clouds
python multistage_ensemble_docking_pipeline.py cluster-diffdock \
  --csv DiffDock_L_results.csv \
  --sdf-dir diffdock_sdf \
  --outdir 01_diffdock_pose_selection \
  --rmsd-cutoff 2.0 \
  --top-k-clusters 3 \
  --top-n-poses 10 \
  --drop-zero-scores

# 2. General rigid receptor preparation using prepare_receptor.txt-style --box_enveloping
python multistage_ensemble_docking_pipeline.py make-scripts \
  --outdir 02_scripts \
  --receptor-pdb pdb_fixed/RP1_receptor.pdb \
  --receptor-name RP1 \
  --receptor-mode general-rigid \
  --box-enveloping ligand_pdbqt/rank1_rp1_ce18_1.pdbqt \
  --padding 4 \
  --selected-sdf-dir ../01_diffdock_pose_selection/selected_poses \
  --receptor-rigid receptor_grid/RP1_receptor_grid.pdbqt \
  --box-file receptor_grid/RP1_box.txt

# 3. Run generated scripts manually after checking paths:
# cd 02_scripts
# bash prepare_ligands_meeko.sh
# bash prepare_receptor_grid.sh
# bash run_vina_refinement.sh

# 4. Parse Vina logs
python multistage_ensemble_docking_pipeline.py parse-vina \
  --log-dir 02_scripts/vina_logs \
  --outdir 03_vina_refinement

# 5. General/CD36RP-style contact analysis
python multistage_ensemble_docking_pipeline.py batch-contact-general \
  --project-root . \
  --proteins RP1 RP3 RP4 \
  --ligands CE18_1 CE18_2 CE18_3 CE20_4 CE20_5 CE16_0 \
  --receptor-template "{protein}/receptor.pdb" \
  --workdir-template "{protein}/{ligand}" \
  --vina-glob "*.pdbqt" \
  --freq-threshold 0.5 \
  --outdir 04_general_contact_analysis

# 6. CEH flexible-output analysis example
# python multistage_ensemble_docking_pipeline.py make-flexdock-ceh-script \
#   --outdir 04_ceh_flexdock_scripts \
#   --flex-out rank1_confidence-2.10_flex_out.pdbqt \
#   --receptor receptor_rigid.pdbqt \
#   --ligand-sdf rank1_confidence-2.10.sdf \
#   --results-dir ceh_flexdock_results
# bash 04_ceh_flexdock_scripts/run_flexdock_v3_ceh.sh
