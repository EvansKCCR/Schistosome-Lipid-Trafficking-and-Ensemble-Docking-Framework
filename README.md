# Integrated Schistosoma Lipid-Handling Protein Discovery, Structure Validation, Docking, and MD Analysis Pipeline

A modular computational workflow for identifying lipid-handling protein candidates, evaluating predicted structures, preparing docking-ready receptors, screening substrate specificity through ensemble docking, and validating receptor–ligand complexes using molecular dynamics (MD).

The pipeline was developed for **Schistosoma lipid biology**, with dedicated support for:

- CD36-related proteins (CD36RPs)
- NPC1/NPC2-like sterol-binding proteins
- CEH-like / α/β-hydrolase enzymes
- General lipid-binding and lipid-transport protein candidates

The workflow is designed to be **reproducible, modular, scalable, and manuscript-ready**. Each stage can be executed independently or as part of an end-to-end discovery and validation pipeline.

---

## Workflow Overview

```text
Proteomes / HMMER outputs
        │
        ▼
Lipid-domain candidate discovery
        │
        ▼
Topology, motif, and hydrolase classification
        │
        ▼
Alignment, phylogeny, and family annotation
        │
        ▼
Deep-learning structure prediction evaluation
        │
        ▼
Structural validation and representative model selection
        │
        ▼
Docking-ready receptor preparation
        │
        ▼
Multi-stage ensemble docking and substrate screening
        │
        ▼
MD validation of receptor–ligand complexes
        │
        ▼
Post-MD stability, contact, and catalytic-state analysis
```

---

## Main Features

### 1. Lipid Candidate Discovery

- Parses HMMER/Pfam `domtblout` outputs.
- Maps Pfam domains to curated lipid-handling protein families.
- Extracts family-specific FASTA files.
- Generates candidate summary tables and Excel workbooks.

### 2. Topology, Motif, and Functional Annotation

- Integrates SignalP, DeepTMHMM, CCTOP, and hydrophobicity-based fallback logic.
- Detects curated motifs from user-provided motif libraries.
- Supports CD36RP, NPC1/NPC2, CEH-like, and α/β-hydrolase motif classes.
- Includes hydrolase triage for secreted, membrane-associated, and cytosolic candidates.

### 3. Alignment and Phylogenetic Analysis

- Runs MAFFT alignments.
- Supports IQ-TREE model selection and maximum-likelihood phylogeny.
- Generates alignment frequency matrices, sequence logos, and annotated trees.
- Collapses weak branches and annotates clades by family-specific keywords.

### 4. Deep-Learning Structure Evaluation

Supports outputs from:

- AlphaFold3
- Boltz / Boltz-2
- Chai-like predictors
- Generic mmCIF/PDB + JSON/NPZ confidence outputs

Extracted metrics include:

- pLDDT
- PAE
- pTM / ipTM
- chain-pair ipTM
- contact probabilities
- ranking score
- motif-level pLDDT
- per-residue and per-chain pLDDT

### 5. Structural Validation and Representative Model Selection

Integrates:

- MolProbity
- CaBLAM
- VoroMQA
- Foldseek
- TMalign
- RMSD / TM-score / Q-score / contact-map overlap

Includes two model-selection engines:

- **General engine** for lipid-binding and transporter-like proteins.
- **Hydrolase-specific engine** for CEH-like enzymes, prioritizing catalytic pocket geometry.

### 6. Receptor Preparation for Docking

- Performs restrained receptor refinement using GROMACS.
- Supports short restrained MD relaxation.
- Generates docking-ready receptors.
- Produces MD validation scripts for RMSD, RMSF, radius of gyration, SASA, DSSP, and catalytic-distance analysis.
- Supports general rigid receptor preparation and CEH-like flexible receptor preparation.

### 7. Multi-Stage Ensemble Docking

Combines:

- DiffDock-L global pose generation
- SMINA rescoring
- pose clustering and representative selection
- AutoDock Vina localized refinement
- rigid or flexible receptor docking
- automated Vina mode splitting
- receptor–ligand complex construction
- residue contact analysis
- substrate-specificity ranking

### 8. CEH-Like Catalytic Docking Engine

For CEH-like / α/β-hydrolase complexes, the pipeline scores:

- Ser207–ligand carbonyl carbon distance
- Bürgi–Dunitz-like attack angle
- Gly128/Gly129 oxyanion-hole contacts
- Ser207–His508 distance
- His508–Glu388 or His508–Glu206 catalytic relay support
- catalytic competence class

Pose classes include:

- `catalytically-competent`
- `near-reactive`
- `binding`
- `distant`

### 9. Post-MD Receptor–Ligand Complex Analysis

Analyzes:

- protein RMSD
- ligand RMSD
- RMSF
- radius of gyration
- SASA
- DSSP secondary-structure content
- protein–ligand contact occupancy
- ligand–pocket COM distance
- CEH catalytic-state persistence

For CEH-like complexes, the post-MD engine calculates catalytic-frame occupancy and catalytic organization scores.

---

## Repository Structure

```text
.
├── pipelines/
│   ├── schisto_lipid_pipeline.py
│   ├── schisto_phylo_motif_pipeline.py
│   ├── deeplearning_structure_evaluation_pipeline.py
│   ├── protein_model_validation_selection_pipeline.py
│   ├── receptor_docking_prep_pipeline.py
│   ├── multistage_ensemble_docking_pipeline.py
│   └── post_md_complex_analysis_pipeline.py
│
├── configs/
│   ├── samples.tsv
│   ├── lipid_domains.tsv
│   ├── motif_library.tsv
│   ├── motifs_default.tsv
│   ├── ceh_engine_config.example.json
│   └── ceh_post_md_config.example.json
│
├── examples/
│   ├── example_metrics_general.tsv
│   ├── example_metrics_ceh.tsv
│   ├── hydrolase_pocket_residues.example.txt
│   └── run_examples.sh
│
├── scripts/
│   ├── run_integrated_pipeline.sh
│   ├── run_model_validation_selection.sh
│   ├── run_receptor_docking_prep.sh
│   ├── run_multistage_ensemble_docking.sh
│   └── run_post_md_complex_analysis.sh
│
├── docs/
│   ├── methodology_protocol.md
│   ├── README_model_validation_selection.md
│   ├── README_receptor_docking_prep.md
│   ├── README_multistage_ensemble_docking.md
│   └── README_post_md_complex_analysis.md
│
├── requirements.txt
├── environment.yml
├── LICENSE
└── README.md
```

---

## Installation

### Option 1: Conda

```bash
conda create -n schisto-lipid-pipeline python=3.10
conda activate schisto-lipid-pipeline
pip install -r requirements.txt
```

### Option 2: Mamba

```bash
mamba create -n schisto-lipid-pipeline python=3.10
mamba activate schisto-lipid-pipeline
pip install -r requirements.txt
```

### Optional External Tools

Some modules require external structural bioinformatics tools:

```text
HMMER
MAFFT
IQ-TREE
SignalP
DeepTMHMM
CCTOP
Foldseek
TMalign
MolProbity / Phenix MolProbity
CaBLAM
VoroMQA
GROMACS
OpenBabel
Meeko
AutoDock Vina
SMINA
RDKit
PyMOL or ChimeraX
```

Not all tools are required for every module. Each pipeline can be run independently depending on the available inputs and installed software.

---

## Quick Start

### 1. Lipid Candidate Discovery

```bash
python pipelines/schisto_lipid_pipeline.py \
  --samples configs/samples.tsv \
  --domain-map configs/lipid_domains.tsv \
  --outdir results/01_lipid_candidates
```

### 2. Motif, Topology, and Phylogeny Integration

```bash
python pipelines/schisto_phylo_motif_pipeline.py \
  --input-dir results/01_lipid_candidates/family_fastas \
  --motif-library configs/motif_library.tsv \
  --deeptmhmm results/deeptmhmm_results.md \
  --signalp results/signalp_results.tsv \
  --outdir results/02_phylo_motif \
  --run-phylogeny \
  --pdf-report
```

### 3. Deep-Learning Structure Evaluation

```bash
python pipelines/deeplearning_structure_evaluation_pipeline.py \
  --input-dir structures/predicted_models \
  --motifs configs/motifs_default.tsv \
  --outdir results/03_structure_evaluation \
  --excel
```

### 4. Model Validation and Representative Selection

General proteins:

```bash
python pipelines/protein_model_validation_selection_pipeline.py \
  --metrics-table quality_matrices.xlsx \
  --engine general \
  --outdir results/04_model_selection_general
```

CEH-like enzymes:

```bash
python pipelines/protein_model_validation_selection_pipeline.py \
  --metrics-table CEH_combined_model_selection.csv \
  --engine hydrolase \
  --family "alpha/beta hydrolase" \
  --outdir results/04_model_selection_ceh
```

### 5. Receptor Docking Preparation

```bash
python pipelines/receptor_docking_prep_pipeline.py \
  select \
  --metrics-table example_metrics_ceh.tsv \
  --engine ceh \
  --outdir results/05_receptor_dockprep
```

### 6. Multi-Stage Ensemble Docking

```bash
python pipelines/multistage_ensemble_docking_pipeline.py \
  cluster-diffdock \
  --input-dir docking/diffdock_outputs \
  --outdir results/06_docking_clusters
```

```bash
python pipelines/multistage_ensemble_docking_pipeline.py \
  parse-vina \
  --logs-dir docking/vina_logs \
  --outdir results/06_vina_summary
```

CEH-like catalytic pose analysis:

```bash
python pipelines/multistage_ensemble_docking_pipeline.py \
  ceh-analyze \
  --base-dir docking/ceh_complexes \
  --config configs/ceh_engine_config.example.json \
  --outdir results/06_ceh_catalytic_pose_analysis
```

### 7. Post-MD Complex Analysis

```bash
python pipelines/post_md_complex_analysis_pipeline.py \
  workbook \
  --input trajectory_data.xlsx \
  --engine general \
  --outdir results/07_post_md_general
```

CEH-like complexes:

```bash
python pipelines/post_md_complex_analysis_pipeline.py \
  workbook \
  --input CEH_trajectory_data.xlsx \
  --engine ceh \
  --config configs/ceh_post_md_config.example.json \
  --outdir results/07_post_md_ceh
```

---

## Input Requirements

### Sequence Discovery

- Protein FASTA files
- HMMER/Pfam `domtblout` files
- Lipid-domain mapping table

### Motif and Topology Annotation

- Family-specific FASTA files
- Motif library file
- SignalP output
- DeepTMHMM output
- Optional CCTOP output

### Structure Evaluation

- mmCIF, CIF, or PDB files
- AlphaFold3 / Boltz / Chai JSON or NPZ confidence files
- Optional motif table

### Model Validation

- Predicted or refined structures
- MolProbity / CaBLAM / VoroMQA outputs
- Foldseek or TMalign outputs
- Optional pocket-residue definitions

### Docking

- Prepared receptor structures
- Ligand PDBQT files
- DiffDock-L or SMINA pose outputs
- AutoDock Vina logs and pose outputs

### Post-MD Analysis

- GROMACS `.tpr` and `.xtc` files
- GROMACS `.xvg` metric outputs
- Multi-sheet trajectory workbooks
- Optional catalytic-state distance tables

---

## Key Outputs

The integrated workflow generates:

```text
candidate discovery tables
family-specific FASTA files
motif annotation workbooks
phylogenetic trees
structure confidence summaries
per-residue pLDDT profiles
motif-level pLDDT tables
structural validation workbooks
representative model tables
docking-ready receptor files
DiffDock/Vina pose summaries
protein-ligand contact matrices
CEH catalytic pose rankings
MD comparative summary workbooks
post-MD catalytic-state reports
publication-ready plots
```

---

## CEH-Like Enzyme Workflow

CEH-like enzymes are handled using specialized decision logic throughout the pipeline.

### CEH-Specific Evidence Chain

```text
HGG-like motif
      │
GxSxG / GESAG-like nucleophile loop
      │
Ser–His–acid catalytic relay
      │
high motif-region pLDDT
      │
validated α/β-hydrolase fold
      │
docking pose with productive carbonyl placement
      │
oxyanion-hole stabilization
      │
MD-stable catalytic geometry
```

### CEH-Specific Prioritization Criteria

For CEH-like systems, the final interpretation prioritizes:

- catalytic motif conservation
- motif-region structural confidence
- catalytic pocket geometry
- Ser207–His508 distance
- His508–Glu388 or His508–Glu206 distance
- Ser207-to-ligand carbonyl proximity
- Gly128/Gly129 oxyanion-hole engagement
- persistence of catalytic-state criteria during MD

Binding affinity alone is not considered sufficient evidence for catalytic plausibility.

---

## Reproducibility Notes

To support reproducibility:

- Keep all configuration files under `configs/`.
- Record software versions for all external tools.
- Use fixed random seeds for Monte Carlo scoring and docking where possible.
- Store raw and processed outputs separately.
- Preserve all intermediate tables.
- Use consistent receptor, ligand, and model identifiers across modules.
- Keep figure-generation scripts with the processed data used to make each figure.

Recommended command to record Python package versions:

```bash
pip freeze > environment_freeze.txt
```

Recommended command to record GROMACS version:

```bash
gmx --version > gromacs_version.txt
```

---

## Scaling the Workflow

The pipeline can be scaled using:

- GNU Parallel
- SLURM job arrays
- Snakemake
- Nextflow
- CWL
- Docker or Singularity containers

Suggested parallelization levels:

```text
species
protein family
candidate protein
structure model
receptor conformer
ligand
docking replicate
MD replicate
trajectory
```

---

## Recommended Citation Statement

If using this repository, cite the external tools used in your analysis, including HMMER, MAFFT, IQ-TREE, AlphaFold3 or other structure-prediction platforms, Foldseek, TMalign, MolProbity, VoroMQA, GROMACS, AutoDock Vina, SMINA, DiffDock-L, MDAnalysis, and any ligand or force-field tools used for simulation setup.

---

## Limitations

- Structure-prediction confidence does not guarantee functional correctness.
- Docking affinity alone should not be interpreted as substrate specificity.
- For CEH-like enzymes, catalytic geometry must be evaluated separately from binding energy.
- MD trajectory length affects interpretation of stability and catalytic-state persistence.
- Force-field and ligand-parameter quality can strongly influence MD-derived conclusions.
- Some modules depend on external tools that must be installed separately.

---

## License

MIT License

---

## Author

Evans Asamoah Adu  
Structural Bioinformatics / Schistosoma Lipid Biology  
GitHub: `@EvansKCCR`

---

## Repository Status

This repository is under active development. Scripts are modular and intended to support reproducible computational analysis of lipid-handling proteins, docking-ready receptor preparation, substrate-specificity screening, and post-MD receptor–ligand interpretation.
