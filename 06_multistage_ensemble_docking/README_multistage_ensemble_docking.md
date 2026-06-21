# Multi-Stage Ensemble Docking and Substrate Specificity Screening Pipeline

This package implements a modular workflow for screening substrate specificity across receptor ensembles using machine-learning-guided global docking followed by empirical local refinement and post-refinement structural analysis.

The updated workflow explicitly includes three post-DiffDock/Vina modules that are essential for this project:

1. **General rigid receptor preparation** using the `prepare_receptor.txt` pattern with `mk_prepare_receptor.py`, `--box_enveloping`, and `--padding`.
2. **General/CD36RP-style Vina post-processing** based on `batch_contact_cd36rp.py`: split Vina modes, build receptor-ligand complexes, and compute contact frequencies automatically.
3. **CEH-like flexible docking engine** based on `flexdock_v3`: parse flexible Vina output ensembles, cluster His508 conformations, cluster ligand poses, select representative modes, and build catalytically interpretable receptor-ligand complexes.

---

## Conceptual workflow

```text
Ligand library + receptor ensemble
        │
        ▼
Stage 1. DiffDock-L blind global docking
        │
        ├─ 100 global pose candidates per receptor-ligand pair
        ├─ SMINA local minimization/rescoring
        └─ pose cloud CSV + SDF pose files
        │
        ▼
Stage 2. DiffDock-L / SMINA pose clustering
        │
        ├─ ligand RMSD clustering when RDKit is available
        ├─ score-space fallback clustering when structures are unavailable
        ├─ major-cluster representative selection
        └─ top-N pose selection for localized refinement
        │
        ▼
Stage 3. Receptor and ligand preparation
        │
        ├─ general rigid receptor preparation using --box_enveloping
        ├─ CEH-like flexible receptor preparation when catalytic flexibility is required
        └─ ligand PDBQT generation using OpenBabel + Meeko
        │
        ▼
Stage 4. AutoDock Vina localized refinement
        │
        ├─ rigid Vina refinement for general receptors
        ├─ optional flexible-residue Vina refinement for CEH-like enzymes
        └─ high-exhaustiveness local resampling of selected global poses
        │
        ▼
Stage 5. Post-refinement analysis
        │
        ├─ general/CD36RP mode splitting, complex construction, and contact analysis
        ├─ Vina log parsing and affinity/RMSD summary
        ├─ CEH flexdock_v3 ensemble clustering and representative-complex construction
        └─ CEH catalytic-competence scoring when applicable
        │
        ▼
Integrated substrate-specificity ranking
```

---

## Installation

```bash
conda create -n ensemble-docking python=3.11 -y
conda activate ensemble-docking
pip install -r requirements_multistage_ensemble_docking.txt
```

Recommended external tools:

```text
AutoDock Vina
OpenBabel
Meeko
RDKit
```

The CEH flexible-docking engine additionally uses the included `flexdock_v3` package and requires:

```text
numpy
pandas
scipy
scikit-learn
rdkit
matplotlib
seaborn
```

---

## Input DiffDock-L / SMINA table

The DiffDock-L/SMINA input CSV should contain at least:

```text
Receptor
Ligand
Prediction
```

Recommended columns:

```text
DiffDock Confidence
SMINA Minimized Affinity
SMINA Minimized RMSD
SMINA Affinity
```

Flexible column aliases such as `receptor`, `ligand`, `prediction`, `confidence`, `smina_minimized_affinity`, and `smina_minimized_rmsd` are accepted.

---

## Step 1: Cluster DiffDock-L / SMINA pose clouds

```bash
python multistage_ensemble_docking_pipeline.py cluster-diffdock \
  --csv DiffDock_L_results.csv \
  --sdf-dir diffdock_sdf \
  --outdir 01_diffdock_pose_selection \
  --rmsd-cutoff 2.0 \
  --top-k-clusters 3 \
  --top-n-poses 10 \
  --drop-zero-scores
```

Outputs:

```text
01_diffdock_pose_selection/diffdock_membership.tsv
01_diffdock_pose_selection/diffdock_clusters.tsv
01_diffdock_pose_selection/diffdock_representatives.tsv
01_diffdock_pose_selection/diffdock_top_poses.tsv
01_diffdock_pose_selection/selected_poses/*.sdf
01_diffdock_pose_selection/diffdock_pose_selection.xlsx
01_diffdock_pose_selection/plots/
```

---

## Step 2A: General rigid receptor preparation

For general rigid receptors such as CD36RPs, use the ligand-enveloping grid command pattern from `prepare_receptor.txt`:

```bash
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
```

This generates a script containing the project-specific rigid receptor command:

```bash
mk_prepare_receptor.py -i pdb_fixed/RP1_receptor.pdb \
  -p receptor_grid/RP1_receptor_grid.pdbqt \
  --box_enveloping ligand_pdbqt/rank1_rp1_ce18_1.pdbqt \
  --padding 4 \
  -v receptor_grid/RP1_box.txt
```

This mode is recommended for standard rigid Vina refinement when a representative ligand pose defines the local search box.

---

## Step 2B: CEH-like flexible receptor preparation

For CEH-like enzymes, flexible receptor preparation may be used when catalytic or gatekeeping side chains must adapt during local refinement. For example, His508 can be prepared as a flexible residue:

```bash
python multistage_ensemble_docking_pipeline.py make-scripts \
  --outdir 02_ceh_scripts \
  --receptor-pdb receptor.pdb \
  --receptor-name receptor \
  --box-center 53.75 45.5 56.5 \
  --box-size 33.5 27.5 24.5 \
  --flex-residues A:508 \
  --receptor-flex receptor_grid/receptor_flex.pdbqt
```

---

## Step 3: Prepare ligands and run AutoDock Vina refinement

After reviewing paths in the generated scripts:

```bash
cd 02_scripts
bash prepare_ligands_meeko.sh
bash prepare_receptor_grid.sh
bash run_vina_refinement.sh
```

Default Vina settings are intentionally stringent for local refinement of a reduced pose set:

```text
exhaustiveness = 64
num_modes = 50
energy_range = 3
seed = 42
```

---

## Step 4: Parse Vina logs

```bash
python multistage_ensemble_docking_pipeline.py parse-vina \
  --log-dir 02_scripts/vina_logs \
  --outdir 03_vina_refinement
```

Outputs:

```text
03_vina_refinement/vina_modes_long.tsv
03_vina_refinement/vina_summary.tsv
03_vina_refinement/vina_refinement_summary.xlsx
03_vina_refinement/plots/top_vina_affinities.png
```

---

## Step 5A: General/CD36RP Vina mode splitting and contact analysis

The `batch-contact-general` module is a parameterized version of `batch_contact_cd36rp.py`. It automatically:

1. finds Vina PDBQT outputs for each protein-ligand pair;
2. splits multi-MODEL Vina PDBQT files into individual pose PDB files using OpenBabel;
3. constructs receptor-ligand complex PDB files;
4. calculates residue contact frequencies using a 4.5 Å default cutoff;
5. writes all-contact, filtered-contact, binding-statistics, and plot outputs.

Example for RP1/RP3/RP4 and CE ligands:

```bash
python multistage_ensemble_docking_pipeline.py batch-contact-general \
  --project-root . \
  --proteins RP1 RP3 RP4 \
  --ligands CE18_1 CE18_2 CE18_3 CE20_4 CE20_5 CE16_0 \
  --receptor-template "{protein}/receptor.pdb" \
  --workdir-template "{protein}/{ligand}" \
  --vina-glob "*.pdbqt" \
  --freq-threshold 0.5 \
  --outdir 04_general_contact_analysis
```

Key outputs:

```text
04_general_contact_analysis/ALL_contact_summary.tsv
04_general_contact_analysis/binding_statistics.tsv
04_general_contact_analysis/*_contacts_all.tsv
04_general_contact_analysis/*_contacts_filtered.tsv
04_general_contact_analysis/general_vina_contact_analysis.xlsx
04_general_contact_analysis/complexes/
04_general_contact_analysis/plots/
```

---

## Step 5B: CEH-like flexible-output ensemble analysis with flexdock_v3

For CEH-like enzymes, flexible Vina output should not be analyzed as a simple set of rigid poses. The bundled `flexdock_v3` engine parses every Vina MODEL, clusters His508 conformations, clusters ligand poses within each His508 cluster, selects the lowest-energy representative from each combined cluster, and builds representative receptor-ligand complexes.

Generate a CEH flexdock script:

```bash
python multistage_ensemble_docking_pipeline.py make-flexdock-ceh-script \
  --outdir 04_ceh_flexdock_scripts \
  --flex-out rank1_confidence-2.10_flex_out.pdbqt \
  --receptor receptor_rigid.pdbqt \
  --ligand-sdf rank1_confidence-2.10.sdf \
  --results-dir ceh_flexdock_results \
  --expected-ligand-atoms 47 \
  --his-cutoff 0.75 \
  --ligand-cutoff 2.0
```

Run the generated script from the package directory:

```bash
bash 04_ceh_flexdock_scripts/run_flexdock_v3_ceh.sh
```

Or run directly if `flexdock_v3` is on `PYTHONPATH`:

```bash
python multistage_ensemble_docking_pipeline.py run-flexdock-ceh \
  --flex-out rank1_confidence-2.10_flex_out.pdbqt \
  --receptor receptor_rigid.pdbqt \
  --ligand-sdf rank1_confidence-2.10.sdf \
  --results-dir ceh_flexdock_results
```

Important flexdock_v3 outputs:

```text
ceh_flexdock_results/pose_summary.csv
ceh_flexdock_results/his508_clusters.csv
ceh_flexdock_results/ligand_clusters.csv
ceh_flexdock_results/representative_modes.csv
ceh_flexdock_results/representative_complexes/complex_manifest.csv
ceh_flexdock_results/representative_complexes/*.pdb
ceh_flexdock_results/representative_complexes/*_ligand.sdf
```

The generated representative complexes can then be analyzed using the CEH catalytic-geometry engine:

```bash
python multistage_ensemble_docking_pipeline.py analyze-complexes \
  --complex-dir ceh_flexdock_results/representative_complexes \
  --engine ceh \
  --ser-resid 207 \
  --his-resid 508 \
  --acid-resid 388 \
  --outdir 05_ceh_catalytic_analysis
```

---

## CEH catalytic-competence engine

CEH-like and α/β-hydrolase poses are scored using active-site geometry rather than docking energy alone.

Default CEH features:

```text
Ser207 OG       nucleophile
His508 NE2/ND1 catalytic histidine
Glu388 OE1/OE2 acidic residue
Gly127-Gly129   HGG-like oxyanion region / backbone donors
```

The CEH engine evaluates:

- ligand carbonyl identification;
- Ser OG to ligand carbonyl carbon distance;
- Bürgi-Dunitz-like attack angle;
- oxyanion-hole contacts;
- Ser-His distance;
- acid-His distance;
- near-attack conformation status;
- catalytic competence score;
- residue contact hotspots.

Score composition:

```text
0.30 × Ser-carbonyl distance score
0.20 × Bürgi-Dunitz angle score
0.30 × oxyanion contact score
0.20 × catalytic triad score
```

Pose classes:

```text
catalytically-competent
near-reactive
binding
distant
```

For CEH-like enzymes, final substrate prioritization should require both favorable Vina affinity and productive catalytic geometry.

---

## Integrated ranking

The final ranking combines all available evidence.

### General engine

```text
Vina affinity          high weight
DiffDock/SMINA pose    moderate weight
Vina RMSD stability    moderate/low weight
contact consistency    supporting evidence
```

### CEH-like engine

```text
CEH catalytic competence    highest weight
Vina affinity               high weight
DiffDock/SMINA pose          supporting weight
Vina RMSD stability          supporting weight
His508/ligand cluster support supporting evidence
```

This prevents CEH-like substrates from being prioritized solely by binding energy when the ligand carbonyl is not positioned for catalysis.

---

## Typical full workflow

```bash
python multistage_ensemble_docking_pipeline.py workflow \
  --diffdock-csv DiffDock_L_results.csv \
  --sdf-dir diffdock_sdf \
  --receptor-pdb pdb_fixed/RP1_receptor.pdb \
  --receptor-mode general-rigid \
  --box-enveloping ligand_pdbqt/rank1_rp1_ce18_1.pdbqt \
  --padding 4 \
  --vina-log-dir vina_logs \
  --complex-dir refined_complexes \
  --engine general \
  --outdir multistage_ensemble_docking
```

For CEH-like flexible docking, run `make-flexdock-ceh-script` or `run-flexdock-ceh` after Vina refinement to analyze flexible-output ensembles and construct representative complexes.

---

## Notes and limitations

1. DiffDock-L and SMINA outputs from Neurosnap are treated as global-placement evidence.
2. Vina is used as a localized empirical refinement step, not as a replacement for global pose exploration.
3. General rigid receptor preparation should use the `--box-enveloping` mode when a selected ligand pose defines the local docking box.
4. General/CD36RP contact analysis requires OpenBabel for automatic PDBQT mode splitting.
5. CEH flexible-output analysis requires the bundled `flexdock_v3` package and RDKit-compatible ligand SDF topology.
6. CEH residue numbering must be adjusted for other hydrolase homologs.
7. Final biological interpretation should consider receptor preparation quality, MD stability, and experimental plausibility.
