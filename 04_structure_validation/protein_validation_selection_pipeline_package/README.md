# Integrated Protein Structural Validation and Representative Model Selection Pipeline

## Overview

This pipeline integrates protein model quality assessment, structural-template comparison, docking-readiness scoring, and representative-model selection into a single reproducible workflow. It was designed for predicted structural ensembles from AlphaFold, AlphaFold3, Boltz/Boltz-2, and related deep-learning platforms, but it can also be used for homology models or experimentally guided models.

The pipeline supports two selection engines:

1. **General structural-validation engine** for lipid-binding, lipid-transport, membrane-associated, and soluble protein models.
2. **Hydrolase-specific selection engine** for cholesteryl-ester hydrolase-like and α/β-hydrolase candidates. This engine gives greater weight to catalytic-pocket geometry, catalytic-distance preservation, and motif-region confidence.

The implementation consolidates logic from the following workflows:

- MolProbity/CaBLAM/VoroMQA validation.
- VoroMQA global and per-residue packing summaries.
- Foldseek template search and best-hit extraction.
- TMalign/TM-score structural comparison.
- Docking-readiness Monte-Carlo-style component weighting.
- CEH/αβ-hydrolase ensemble benchmarking and catalytic-site model selection.

---

## Main script

```bash
protein_model_validation_selection_pipeline.py
```

---

## Input options

The pipeline can be run from an existing metrics table or from model structures plus optional external tools.

### Option A — score a precomputed metrics table

Use this mode when MolProbity, CaBLAM, VoroMQA, Foldseek, TMalign, pLDDT, or hydrolase ensemble metrics have already been collected.

Accepted table formats:

- `.csv`
- `.tsv`
- `.xlsx`

Recommended columns include any of the following:

| Metric class | Example columns |
|---|---|
| Model identifiers | `model`, `File`, `Short_name`, `Protein_ID` |
| Global validation | `molprobity_score`, `clashscore`, `rama_favored_%`, `rama_outliers_%`, `rotamer_outliers_%`, `cbeta_outliers` |
| Backbone validation | `cablam_disfavored_%`, `cablam_outliers_%`, `cablam_severe_%` |
| Packing | `voromqa_dark_score`, `voromqa_light_score`, `perres_mean`, `perres_q10_mean`, `pocket_score` |
| Fold/template similarity | `TM_score`, `Q_score`, `CMO`, `RMSD`, `best_rmsd_global` |
| Hydrolase ensemble metrics | `global_rmsd_A`, `pocket_rmsd_A`, `catalytic_distance_delta_A`, `global_cluster`, `pocket_cluster` |
| Confidence metrics | `atom_plddts_mean`, `ptm`, `iptm`, `ranking_score` |
| Hydrolase motif confidence | `HGG_plddt_mean`, `GESAG_plddt_mean`, `HGG_plddt_count`, `GESAG_plddt_count` |

Column aliases are normalized automatically. For example, `rama_outliers_%` and `rama_outliers` are treated as the same metric.

### Option B — run optional external validation from structures

Accepted coordinate formats:

- `.pdb`
- `.cif`
- `.mmcif`

Optional tools are only called when the corresponding flags are supplied and the tools are available on `PATH`:

- `phenix.molprobity` or `molprobity`
- `phenix.cablam` or `cablam`
- `voronota-js-voromqa`
- `TMalign`
- `foldseek`

---

## Quick start

### 1. Hydrolase/CEH-specific representative selection

```bash
python protein_model_validation_selection_pipeline.py \
  --metrics-table CEH_combined_model_selection.csv \
  --engine hydrolase \
  --family "alpha/beta hydrolase" \
  --outdir ceh_representative_selection \
  --top-n 10
```

This mode applies the hydrolase-specific engine and prioritizes:

- low pocket RMSD,
- catalytic-distance preservation,
- low global RMSD,
- high HGG/GESAG motif pLDDT,
- good MolProbity/CaBLAM geometry,
- acceptable VoroMQA packing.

### 2. General protein model selection

```bash
python protein_model_validation_selection_pipeline.py \
  --metrics-table quality_matrices.xlsx \
  --sheet "Table S6" \
  --engine general \
  --outdir general_model_selection \
  --top-n 10
```

This mode balances:

- stereochemical quality,
- backbone geometry,
- structural similarity to templates,
- packing quality,
- model confidence values when available.

### 3. Use automatic engine selection

```bash
python protein_model_validation_selection_pipeline.py \
  --metrics-table combined_metrics.csv \
  --engine auto \
  --outdir model_selection_auto
```

`auto` switches to the hydrolase engine when catalytic-pocket or motif-confidence columns are detected.

### 4. Run from coordinate files with optional external tools

```bash
python protein_model_validation_selection_pipeline.py \
  --models-dir models \
  --convert-cif \
  --run-external-validation \
  --outdir structural_validation_results
```

### 5. Add TMalign template comparison

```bash
python protein_model_validation_selection_pipeline.py \
  --models-dir models \
  --templates-dir templates \
  --run-tmalign \
  --run-external-validation \
  --outdir validation_with_templates
```

### 6. Add Foldseek best-hit search

```bash
python protein_model_validation_selection_pipeline.py \
  --models-dir models \
  --run-foldseek \
  --foldseek-db pdb \
  --outdir foldseek_validation_selection
```

---

## Outputs

Each run writes the following files into the output directory:

| Output | Description |
|---|---|
| `all_models_scored.tsv` | Complete scored model table with normalized component scores. |
| `representative_models.tsv` | Selected representative models, including overall top models and cluster/family representatives. |
| `structural_validation_selection.xlsx` | Excel workbook containing all scored models, selected representatives, and score components. |
| `validation_selection_report.md` | Human-readable summary of the run and selected models. |
| `plots/top_selection_scores.png` | Bar plot of top-ranked models. |
| `plots/hydrolase_pocket_geometry_selection.png` | Hydrolase-specific pocket geometry plot, generated when relevant columns are available. |
| `external_validation/` | Optional MolProbity/CaBLAM/VoroMQA outputs if external validation is run. |
| `foldseek/` | Optional Foldseek outputs. |
| `tmalign/` | Optional TMalign outputs. |

---

## General selection engine

The general engine calculates component scores from available metrics using robust 5th/95th percentile scaling. Metric weights are renormalized automatically when some columns are absent.

### Components

| Component | Main metrics |
|---|---|
| `S_geometry` | MolProbity score, clashscore, Ramachandran favored/outliers, rotamer outliers, Cβ deviations, bond/angle RMSD |
| `S_backbone` | CaBLAM disfavored, outlier, and severe scores |
| `S_topology` | TM-score, Q-score, CMO, RMSD |
| `S_packing` | VoroMQA dark/light scores, VoroMQA distribution summaries, pocket score |
| `S_confidence` | pLDDT, pTM, ipTM, ranking score |

The final general score is:

```text
SelectionScore = 100 × validation_gate × weighted_mean(
    S_geometry,
    S_backbone,
    S_topology,
    S_packing,
    S_confidence
)
```

Default general validation gate:

```text
MolProbity score ≤ 3.0
Clashscore ≤ 20
Ramachandran outliers ≤ 2.0%
CaBLAM severe ≤ 2.0%
```

Missing values do not automatically fail the gate.

---

## Hydrolase-specific selection engine

The hydrolase engine is stricter because α/β-hydrolase model selection depends not only on global fold quality but also on catalytic-site reliability.

### Hydrolase-specific components

| Component | Main metrics |
|---|---|
| `S_validation` | MolProbity, clashscore, Ramachandran, rotamer, CaBLAM, VoroMQA |
| `S_ensemble` | pocket RMSD, catalytic-distance delta, global RMSD |
| `S_functional_site` | HGG motif pLDDT, GESAG/GxSxG motif pLDDT, whole-model pLDDT, pTM/ipTM |

The hydrolase score is:

```text
SelectionScore = 100 × validation_gate × motif_completeness_modifier × weighted_mean(
    S_validation,
    S_ensemble,
    S_functional_site
)
```

Default hydrolase validation gate:

```text
MolProbity score ≤ 2.5
Clashscore ≤ 15
Ramachandran outliers ≤ 1.0%
CaBLAM severe ≤ 1.0%
```

Hydrolase ensemble weights:

```text
pocket_rmsd_A                  0.45
|catalytic_distance_delta_A|   0.35
global_rmsd_A                  0.20
```

This means a hydrolase model with excellent global geometry but poor catalytic-site preservation will not be over-prioritized.

---

## Representative selection logic

The pipeline selects:

1. the top `N` overall models;
2. the best model per cluster;
3. the best model per detected family when family information is available.

For hydrolases, the cluster representative is selected using `pocket_cluster` by default. For general proteins, it uses `global_cluster` by default.

---

## Decision classes

| Selection score | Decision |
|---:|---|
| ≥ 75 | High-confidence representative |
| 60–74.99 | Representative candidate |
| 45–59.99 | Borderline; inspect/refine |
| < 45 | Reject or remodel |

Thresholds can be changed using:

```bash
--high-threshold 75 --pass-threshold 60 --borderline-threshold 45
```

---

## Recommended workflow order

For a complete structural-model evaluation study, use the following order:

```text
1. Generate predicted structures using AlphaFold3, Boltz-2, or another platform.
2. Extract confidence metrics and motif-level pLDDT using the deep-learning structure evaluation pipeline.
3. Run MolProbity, CaBLAM, VoroMQA, TMalign, and/or Foldseek as needed.
4. For hydrolases, compute ensemble pocket RMSD and catalytic-distance metrics.
5. Run this integrated validation-selection pipeline.
6. Inspect representative_models.tsv and validation_selection_report.md.
7. Select the final representative models for docking, MD simulation, or manuscript figures.
```

---

## Requirements

Minimum Python packages:

```bash
pip install pandas numpy matplotlib openpyxl biopython
```

Optional external tools:

```text
MolProbity / Phenix MolProbity
CaBLAM / Phenix CaBLAM
VoroMQA via voronota-js-voromqa
TMalign
Foldseek
```

---

## Notes

- Missing metrics are ignored and available metric weights are renormalized.
- VoroMQA light-score direction is inferred relative to MolProbity when both are available.
- The hydrolase engine is intentionally different from the general engine and should be used for CEH-like or α/β-hydrolase candidates.
- For non-hydrolase lipid-binding proteins such as CD36RPs, NPC1, and NPC2, the general engine is usually more appropriate unless catalytic-site metrics are explicitly defined.
