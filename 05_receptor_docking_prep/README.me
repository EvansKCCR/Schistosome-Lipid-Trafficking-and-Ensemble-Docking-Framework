# Integrated Receptor Preparation Pipeline for Docking

## Purpose

This pipeline prepares protein receptors for docking through a modular workflow that links:

1. **restrained receptor refinement**;
2. **MD validation**;
3. **comparative trajectory plotting**; and
4. **representative receptor/model selection**.

It is designed for predicted and experimentally derived structures, including lipid-handling proteins such as CD36-like receptors, NPC1/NPC2 proteins, and CEH-like / α/β-hydrolase enzymes.

A special **CEH-like enzyme workflow** is included for catalytic-site preservation checks involving the Ser–His–acid catalytic network.

---

## Pipeline architecture

```text
Input receptor structures
        │
        ▼
Module 1: restrained refinement
        │
        ├─ PDB cleanup/topology generation
        ├─ solvation and neutralization
        ├─ backbone-restrained energy minimization
        ├─ optional short restrained MD relaxation
        └─ final restrained minimization
        │
        ▼
Refined receptor candidates
        │
        ▼
Module 2: MD validation
        │
        ├─ trajectory centering and fitting
        ├─ RMSD
        ├─ RMSF
        ├─ radius of gyration
        ├─ SASA
        ├─ DSSP secondary structure
        └─ CEH catalytic-distance analysis when enabled
        │
        ▼
Module 3: comparative plotting
        │
        ├─ comparative RMSD plots
        ├─ comparative RMSF plots
        ├─ comparative radius of gyration plots
        ├─ comparative SASA plots
        ├─ comparative DSSP plots
        └─ catalytic-distance plots
        │
        ▼
Module 4: representative receptor selection
        │
        ├─ general receptor selection engine
        └─ CEH/hydrolase catalytic-geometry engine
        │
        ▼
Docking-ready representative receptors
```

---

## Files in this package

| File | Purpose |
|---|---|
| `receptor_docking_prep_pipeline.py` | Main integrated Python pipeline |
| `run_receptor_docking_prep.sh` | Example execution script |
| `requirements_receptor_docking_prep.txt` | Python dependencies |
| `mdp/ions.mdp` | Ionization parameter file |
| `mdp/em.mdp` | Restrained minimization parameter file |
| `mdp/md_relax.mdp` | Short restrained relaxation MD parameter file |
| `mdp/final_em.mdp` | Final minimization parameter file |
| `examples/example_metrics_general.tsv` | Example input for general selection engine |
| `examples/example_metrics_ceh.tsv` | Example input for CEH/hydrolase selection engine |

---

## Installation

```bash
conda create -n dockprep python=3.10 -y
conda activate dockprep
pip install -r requirements_receptor_docking_prep.txt
```

External tools are only required for the modules that call them:

- GROMACS for refinement and MD validation;
- DSSP support through `gmx dssp` for secondary-structure analysis;
- MolProbity, CaBLAM, VoroMQA, Foldseek, or TMalign if their outputs are used in the selection table.

The Python pipeline can still generate refinement and MD-analysis shell scripts even when these external tools are not installed.

---

## 1. Generate restrained-refinement scripts

Place receptor PDB files in `input/` and run:

```bash
python receptor_docking_prep_pipeline.py \
  --steps refine \
  --input-dir input \
  --mdp-dir mdp \
  --outdir dockprep_results
```

This writes:

```text
dockprep_results/scripts/refinement/refine_<model>.sh
dockprep_results/scripts/run_all_refinement.sh
```

The generated refinement workflow performs topology generation, box construction, solvation, ion addition, backbone position restraint generation, restrained minimization, optional short restrained MD, final minimization, and receptor extraction.

To run the generated scripts directly:

```bash
bash dockprep_results/scripts/run_all_refinement.sh
```

or add `--execute` to the Python command.

---

## 2. Generate MD-validation scripts

Each MD directory should contain:

```text
md.tpr
md.xtc
```

For a general receptor:

```bash
python receptor_docking_prep_pipeline.py \
  --steps md-validation \
  --mode general \
  --md-root md_runs \
  --outdir dockprep_results
```

For a CEH-like enzyme:

```bash
python receptor_docking_prep_pipeline.py \
  --steps md-validation \
  --mode ceh \
  --md-root md_runs \
  --ceh-ser 207 \
  --ceh-his 508 \
  --ceh-acid 388 \
  --ceh-alt-acid 206 \
  --outdir dockprep_results
```

The CEH mode generates GROMACS commands for:

- Ser207 OG – His508 NE2;
- His508 ND1 – Glu388 OE1/OE2;
- optional His508 – Glu206 OE1/OE2 alternative checks.

These defaults can be changed with the `--ceh-*` options.

---

## 3. Generate comparative RMSD, RMSF, Rg, SASA, DSSP, and distance plots

After MD analysis, or if you already have `.xvg` files:

```bash
python receptor_docking_prep_pipeline.py \
  --steps plot \
  --xvg-root dockprep_results/md_analysis \
  --outdir dockprep_results
```

Outputs:

```text
dockprep_results/md_xvg_summary.tsv
dockprep_results/plots/comparative_rmsd.png
dockprep_results/plots/comparative_rmsf.png
dockprep_results/plots/comparative_rg.png
dockprep_results/plots/comparative_sasa.png
dockprep_results/plots/comparative_dssp.png
dockprep_results/plots/comparative_distance.png
```

The DSSP plotting logic generalizes the original `plot_dssp.py` approach by using the XVG legends as secondary-structure labels and plotting each secondary-structure count over time.

---

## 4. Select representative receptors

### General receptor engine

```bash
python receptor_docking_prep_pipeline.py \
  --steps select \
  --mode general \
  --metrics-table metrics_general.tsv \
  --top-n 3 \
  --outdir dockprep_results
```

The general engine prioritizes:

- MolProbity/CaBLAM/VoroMQA validation;
- RMSD/RMSF/Rg/SASA stability;
- fold or template agreement when TM-score, Q-score, or contact-map overlap are available.

### CEH/hydrolase engine

```bash
python receptor_docking_prep_pipeline.py \
  --steps select \
  --mode ceh \
  --metrics-table metrics_ceh.tsv \
  --top-n 3 \
  --outdir ceh_dockprep_results
```

The CEH/hydrolase engine applies stricter validation gates and gives higher priority to:

- catalytic-pocket RMSD;
- catalytic-distance deviation;
- catalytic-distance stability;
- HGG motif confidence;
- GESAG/GxSxG nucleophile-loop confidence;
- global validation and MD stability.

Outputs:

```text
all_receptors_scored.tsv
representative_receptors.tsv
receptor_selection_workbook.xlsx
receptor_selection_report.md
plots/top_receptor_selection_scores.png
plots/ceh_validation_vs_catalytic_geometry.png
```

---

## 5. Extract representative frames from trajectories

If the representative table contains a `time_ps` or `representative_time_ps` column:

```bash
python receptor_docking_prep_pipeline.py \
  --steps extract-frames \
  --md-root md_runs \
  --selection-table dockprep_results/representative_receptors.tsv \
  --top-n 3 \
  --outdir dockprep_results
```

This writes:

```text
dockprep_results/scripts/extract_representative_frames.sh
```

Run the script to extract docking-ready receptor conformations with `gmx trjconv`.

---

## Recommended CEH-like enzyme workflow

For CEH-like enzymes, use the following order:

```bash
# 1. Refine starting models
python receptor_docking_prep_pipeline.py \
  --steps refine \
  --mode ceh \
  --input-dir input \
  --mdp-dir mdp \
  --outdir ceh_dockprep

# 2. Run short MD or production MD externally using the refined models

# 3. Generate MD-validation scripts with CEH catalytic checks
python receptor_docking_prep_pipeline.py \
  --steps md-validation \
  --mode ceh \
  --md-root md_runs \
  --ceh-ser 207 \
  --ceh-his 508 \
  --ceh-acid 388 \
  --ceh-alt-acid 206 \
  --outdir ceh_dockprep

# 4. Plot comparative MD diagnostics
python receptor_docking_prep_pipeline.py \
  --steps plot \
  --xvg-root ceh_dockprep/md_analysis \
  --outdir ceh_dockprep

# 5. Select representative docking receptors with hydrolase-specific scoring
python receptor_docking_prep_pipeline.py \
  --steps select \
  --mode ceh \
  --metrics-table ceh_metrics.tsv \
  --top-n 3 \
  --outdir ceh_dockprep
```

---

## Suggested metrics table columns

The selection engine accepts many column aliases. Useful columns include:

### Validation

```text
model
molprobity_score
clashscore
rama_outliers_%
rama_favored_%
rotamer_outliers_%
cablam_severe_%
cablam_outliers_%
voromqa_light_score
```

### MD stability

```text
rmsd_mean
rmsd_std
rmsf_mean
rg_std
sasa_std
```

### Fold/template agreement

```text
TM_score
Q_score
CMO
```

### CEH/hydrolase-specific

```text
pocket_rmsd_A
catalytic_distance_delta_A
distance_ser_his_std
HGG_plddt_mean
GESAG_plddt_mean
```

---

## Notes for docking use

- Do not dock directly into an unvalidated predicted structure when substantial clashes, CaBLAM outliers, or catalytic-site distortions are present.
- For CEH-like enzymes, global RMSD alone is insufficient; preserve the Ser–His–acid geometry and the HGG/GxSxG motif region.
- For flexible receptors, use a small ensemble of representative conformations rather than a single rigid receptor.
- Always inspect selected receptors in PyMOL/ChimeraX before docking.
