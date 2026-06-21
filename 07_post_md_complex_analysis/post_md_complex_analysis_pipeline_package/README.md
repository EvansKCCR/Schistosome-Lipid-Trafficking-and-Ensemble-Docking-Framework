# Integrated Post-MD Analysis Pipeline for Receptor–Ligand Complexes

This package provides a modular post-MD analysis workflow for receptor–ligand complexes, with a dedicated CEH-like enzyme analysis engine.

The pipeline integrates four analysis layers:

1. **GROMACS analysis script generation**  
   Generates reproducible shell scripts for trajectory centering/fitting, RMSD, ligand RMSD, radius of gyration, RMSF, SASA, DSSP, protein–ligand contacts, ligand–pocket COM distance, and CEH catalytic-distance analysis.

2. **Workbook-based comparative analysis**  
   Reads multi-sheet Excel workbooks derived from GROMACS outputs and creates comparative summaries and figures.

3. **Direct residue contact analysis**  
   Uses MDAnalysis to calculate per-residue ligand-contact occupancy directly from topology and trajectory files.

4. **CEH-like catalytic-state engine**  
   Evaluates catalytic organization in CEH-like/αβ-hydrolase complexes using Ser–His, His–acid, Ser–carbonyl, and Gly128/Gly129 oxyanion-hole distance criteria.

---

## Input options

### A. Multi-sheet Excel workbook

Expected workbook sheets can include:

- `C_alpha_RMSD`
- `Lig_RMSD`
- `rGyr`
- `RMSF`
- `SASA` or `Protein_SASA`
- `DSSP_count`
- `P_L_contact`
- `catalytic_state` for CEH-like complexes

The CEH workbook is expected to contain catalytic-state columns such as:

- `ser207_his508`
- `his508_glu388_oe2`
- `ser207_carbonylC_distance`
- `gly128N_O2_distance`
- `gly129N_O2_distance`

### B. Raw trajectory/contact analysis

For direct residue-contact analysis, provide:

- `md.tpr` or compatible topology
- `md.xtc` trajectory
- ligand residue name, usually `UNL`

### C. XVG files

The `xvg-collect` module converts GROMACS `.xvg` files into a workbook index for downstream summarization.

---

## Installation

```bash
conda create -n post_md python=3.10 -y
conda activate post_md
pip install -r requirements_post_md_complex_analysis.txt
```

For direct trajectory contact analysis:

```bash
conda install -c conda-forge mdanalysis -y
```

External tools for generated scripts:

- GROMACS 2024+
- DSSP support through `gmx dssp`

---

## Quick start

### 1. General receptor–ligand workbook analysis

```bash
python post_md_complex_analysis_pipeline.py workbook \
  --input trajectory_data.xlsx \
  --outdir post_MD_general \
  --mode general \
  --formats png,svg \
  --dpi 600
```

### 2. CEH-like receptor–ligand workbook analysis

```bash
python post_md_complex_analysis_pipeline.py workbook \
  --input CEH_trajectory_data.xlsx \
  --outdir post_MD_CEH \
  --mode ceh \
  --hbond-cutoff 0.35 \
  --attack-cutoff 0.45 \
  --oxyanion-cutoff 0.35 \
  --formats png,svg \
  --dpi 600
```

### 3. Direct residue contact analysis from trajectory

```bash
python post_md_complex_analysis_pipeline.py contacts \
  --topology md.tpr \
  --trajectory md.xtc \
  --ligand-resname UNL \
  --cutoff-A 4.5 \
  --stride 20 \
  --outdir residue_contacts
```

### 4. Generate GROMACS analysis script for general complexes

```bash
python post_md_complex_analysis_pipeline.py make-gromacs-scripts \
  --mode general \
  --tpr pre_1.tpr \
  --xtc final_result.part0001.xtc \
  --index index.ndx \
  --ligand-group UNL \
  --out run_general_post_md_analysis.sh
```

### 5. Generate GROMACS analysis script for CEH-like complexes

```bash
python post_md_complex_analysis_pipeline.py make-gromacs-scripts \
  --mode ceh \
  --tpr mdsim360/md.tpr \
  --xtc mdsim360/md.xtc \
  --index mdsim360/index.ndx \
  --ligand-group UNL \
  --ligand-resname UNL \
  --carbonyl-atom C26 \
  --oxyanion-oxygen O2 \
  --out run_CEH_post_md_analysis.sh
```

### 6. Collect XVG files

```bash
python post_md_complex_analysis_pipeline.py xvg-collect \
  --input-dir analysis_xvg \
  --out xvg_collected_workbook.xlsx
```

---

## CEH-like enzyme analysis

The CEH engine is designed for carboxyesterase-like / αβ-hydrolase complexes where global trajectory stability is insufficient to judge mechanistic plausibility. It emphasizes:

- Ser207–His508 catalytic contact persistence
- His508–Glu388 or His508–Glu206 acidic-residue support
- Ser207 Oγ to ligand carbonyl-carbon approach distance
- Gly128/Gly129 oxyanion-hole stabilization of ligand oxygen
- catalytic contact persistence across MD frames
- protein–ligand residue contact occupancy

Default CEH pocket residues:

```text
128, 129, 130, 206, 207, 242, 349, 391, 508
```

Default CEH key catalytic/contact residues:

```text
127, 128, 129, 130, 206, 207, 508
```

---

## Outputs

Workbook mode writes:

```text
post_MD_complex_analysis_summary.xlsx
post_MD_complex_analysis_report.md
figures/
```

The workbook can contain:

- `ReadMe_parameters`
- `Comparative_index`
- `Trajectory_metrics`
- `RMSF_summary`
- `CEH_catalytic_distances`
- `CEH_catalytic_criteria`
- `CEH_catalytic_frame_states`
- `Contacts_long`
- `Contact_summary`
- `Top_contacts`
- `Contact_occupancy_matrix`
- `Contact_distance_matrix`

Figures can include:

- Cα RMSD time series
- ligand RMSD time series
- radius of gyration time series
- SASA time series, if supplied
- RMSF profile
- global metric boxplots
- CEH catalytic-distance profiles
- CEH catalytic-criteria heatmap
- CEH catalytic-score barplot
- protein–ligand contact heatmaps
- top contact residues per system
- DSSP secondary-structure profiles, if supplied

---

## Suggested interpretation

For general receptor–ligand complexes, prioritize systems with stable receptor RMSD, stable ligand RMSD, conserved compactness, limited binding-site RMSF, and persistent residue-contact networks.

For CEH-like complexes, prioritize systems that simultaneously maintain global stability and catalytic organization. A ligand with favorable binding stability but poor Ser–carbonyl approach, weak oxyanion-hole occupancy, or disrupted Ser–His–acid geometry should be interpreted as binding-competent but not necessarily catalytically competent.

---

## Citation notes

When using the workflow, cite the original software used for simulation and analysis, including GROMACS, MDAnalysis if direct contact analysis is used, and DSSP if secondary-structure analysis is reported.
