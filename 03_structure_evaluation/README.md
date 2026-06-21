# Integrated Pipeline for Evaluating AlphaFold3, Boltz-2, and Related Predicted Structures

This package provides a single command-line workflow for evaluating predicted structures from next-generation deep-learning platforms, including **AlphaFold3**, **Boltz/Boltz-2**, **Chai-1-like**, and other tools that export coordinate structures and confidence files.

The pipeline consolidates the logic from the separate AF3 JSON metric scripts, mmCIF pLDDT extraction scripts, catalytic-motif scoring scripts, confidence workbook generation scripts, and batch CIF-to-PDB conversion workflow into one reusable evaluator.

---

## 1. What the pipeline evaluates

For each predicted model group, the pipeline extracts and summarizes:

1. **Coordinate-derived pLDDT**
   - Reads `.cif`, `.mmcif`, and `.pdb` files.
   - Extracts pLDDT from mmCIF `_atom_site.B_iso_or_equiv` or the PDB B-factor column.
   - Reports per-residue, per-chain, and whole-structure confidence.

2. **Confidence JSON metrics**
   - Recursively scans JSON outputs for common AlphaFold3/Boltz-style keys:
     - `pae`, `predicted_aligned_error`, `chain_pair_pae_min`
     - `plddt`, `pLDDT`, `atom_plddts`, `per_residue_plddt`
     - `contact_prob`, `contact_probability`, `predicted_contacts`
     - `ptm`, `iptm`, `ranking_score`, `has_clash`, `fraction_disordered`
     - Boltz/Boltz-2-style affinity fields such as `affinity_pred_value` and `affinity_probability_binary`

3. **NPZ confidence arrays**
   - Parses `.npz` outputs such as Boltz-style `pae_*.npz`, `pde_*.npz`, and `plddt_*.npz` where NumPy is installed.

4. **Catalytic or functional motif confidence**
   - Scans motifs in reconstructed chain sequences.
   - Reports motif hits using **1-based sequence positions** and residue labels.
   - Calculates motif-level atom-averaged pLDDT.
   - Defaults to CEH/HSL-relevant motifs: `HGG,GESAG`.

5. **Integrated model ranking**
   - Produces a transparent heuristic `integrated_quality_score_0_100` using available structure pLDDT, motif pLDDT, JSON pLDDT, PAE, pTM, ipTM, and ranking score.
   - This score is for prioritization only; final interpretation should still use biological context, active-site geometry, and visual inspection.

6. **Optional CIF to PDB conversion**
   - Converts `.cif`/`.mmcif` files to `.pdb` using Biopython for downstream tools that require PDB format.

---

## 2. External format notes

AlphaFold3 output directories commonly contain coordinate models and confidence JSON files; AlphaFold3 reports pLDDT as a per-atom confidence estimate and provides additional confidence metrics such as PAE, pTM, and ipTM in JSON outputs.

Boltz/Boltz-2 predictions commonly include mmCIF/PDB coordinate outputs, confidence JSONs, optional affinity JSONs, and confidence arrays such as PAE/PDE/pLDDT NPZ files. Boltz-2 also reports affinity-related outputs for protein-ligand evaluation.

Useful documentation:

- AlphaFold3 output format: https://github.com/google-deepmind/alphafold3/blob/main/docs/output.md
- AlphaFold Server confidence guide: https://alphafoldserver.com/guides
- EMBL-EBI AlphaFold3 quality assessment guide: https://www.ebi.ac.uk/training/online/courses/alphafold/alphafold-3-and-alphafold-server/how-to-assess-the-quality-of-alphafold-3-predictions/
- Boltz prediction documentation: https://github.com/jwohlwend/boltz/blob/main/docs/prediction.md
- Boltz GitHub repository: https://github.com/jwohlwend/boltz

---

## 3. Installation

Create an environment:

```bash
conda create -n dl_struct_eval python=3.10 -y
conda activate dl_struct_eval
pip install -r requirements_deeplearning_structure_eval.txt
```

Required for full functionality:

```text
numpy
pandas
matplotlib
biopython
openpyxl
```

The pipeline still writes TSV tables if `pandas/openpyxl` are absent, but Excel export requires them. CIF-to-PDB conversion requires Biopython.

---

## 4. Recommended input layout

The pipeline recursively scans an input folder for:

```text
.cif
.mmcif
.pdb
.json
.npz
```

Example AlphaFold3-like layout:

```text
AF3_outputs/
├── ceh_model.cif
├── ceh_confidences.json
├── ceh_summary_confidences.json
└── ceh_data.json
```

Example Boltz/Boltz-2-like layout:

```text
boltz_results/
└── predictions/
    └── ceh_complex/
        ├── ceh_complex_model_0.cif
        ├── confidence_ceh_complex_model_0.json
        ├── affinity_ceh_complex_model_0.json
        ├── pae_ceh_complex_model_0.npz
        ├── pde_ceh_complex_model_0.npz
        └── plddt_ceh_complex_model_0.npz
```

The grouping system normalizes common suffixes/prefixes such as:

```text
fold_*
confidence_*
affinity_*
pae_*
pde_*
plddt_*
*_model_0
*_confidences
*_summary_confidences
*_full_data_0
*_data_0
```

---

## 5. Basic run

```bash
python deeplearning_structure_evaluation_pipeline.py predictions \
  --outdir dl_structure_eval_results \
  --motifs HGG,GESAG \
  --motif-mode exact \
  --write-excel
```

For CEH/HSL-like hydrolase evaluation:

```bash
python deeplearning_structure_evaluation_pipeline.py AF3_and_Boltz_outputs \
  --outdir CEH_structure_eval \
  --motifs HGG,GESAG \
  --write-excel \
  --convert-pdb
```

For broad serine-hydrolase regex screening:

```bash
python deeplearning_structure_evaluation_pipeline.py predictions \
  --outdir hydrolase_regex_eval \
  --motifs 'G.S.G,HGG,N[^P][ST]' \
  --motif-mode regex \
  --write-excel
```

---

## 6. Output files

The pipeline creates:

```text
dl_structure_eval_results/
├── tables/
│   ├── model_summary.tsv
│   ├── per_residue_plddt.tsv
│   ├── per_chain_plddt.tsv
│   ├── motif_plddt.tsv
│   ├── json_metrics_long.tsv
│   ├── json_pairwise_long.tsv
│   ├── npz_metrics_long.tsv
│   └── cif_to_pdb_conversions.tsv
├── plots/
│   ├── integrated_quality_summary.png
│   └── *_plddt_profile.png
├── pdb/
│   └── converted PDB files, if --convert-pdb is used
└── deep_learning_structure_evaluation.xlsx, if --write-excel is used
```

### Main output tables

#### `model_summary.tsv`
One row per grouped model. Includes:

- platform label
- associated structure/JSON/NPZ files
- structure mean pLDDT
- motif mean pLDDT
- JSON-derived pLDDT, PAE, ipTM, pTM, ranking score
- Boltz-style affinity fields when present
- integrated quality score

#### `per_residue_plddt.tsv`
One row per residue with:

- chain
- residue number
- residue name
- mean pLDDT
- quality band
- atom count

#### `per_chain_plddt.tsv`
One row per chain with:

- mean/median/min/max pLDDT
- quality band
- residue count
- atom count

#### `motif_plddt.tsv`
One row per motif hit with:

- motif query
- matched sequence
- chain
- 1-based sequence start/end
- residue numbers and labels
- mean/median/min/max motif pLDDT
- quality band

#### `json_metrics_long.tsv`
Long-format summary of all recognized JSON metrics.

#### `json_pairwise_long.tsv`
Long-format chain-pair metrics such as chain-pair ipTM and chain-pair minimum PAE.

#### `npz_metrics_long.tsv`
Summary statistics for arrays found in `.npz` confidence files.

---

## 7. Quality-band interpretation

The pipeline uses conventional pLDDT bands:

| pLDDT range | Pipeline label |
|---:|---|
| ≥90 | `very_high` |
| 70–89.99 | `confident` |
| 50–69.99 | `low` |
| <50 | `very_low` |

For active-site evaluation, the **motif-level pLDDT** is often more informative than the whole-structure mean. A model with high global pLDDT but poor catalytic-motif confidence should be treated cautiously for mechanistic interpretation or docking.

---

## 8. Integrated quality score

The `integrated_quality_score_0_100` is a heuristic ranking metric. It combines available metrics using normalized weights:

- structure mean pLDDT
- motif mean pLDDT
- JSON pLDDT when structure-level pLDDT is unavailable
- ipTM
- pTM
- ranking score
- inverse PAE contribution

This score should be used for triage, not as a replacement for structural inspection, active-site geometry validation, clash analysis, or downstream MD/docking assessment.

---

## 9. Demo run

A small synthetic demo is included.

```bash
cd dl_structure_eval_pipeline_package
bash run_demo.sh
```

Expected outputs:

```text
demo_results/tables/model_summary.tsv
demo_results/tables/motif_plddt.tsv
demo_results/deep_learning_structure_evaluation.xlsx
```

---

## 10. Relationship to the original scripts

This integrated pipeline consolidates:

- `AF3_stats.py` logic for recursive JSON extraction of PAE, pLDDT, and contact probabilities.
- `confidence_score_summary.py` logic for AF3-style summary, per-chain, and chain-pair confidence outputs.
- `AF3_CEH_stats_v1.py` logic for pairing JSON and mmCIF files and computing catalytic motif pLDDT.
- `extract_mmcif_plddt.py` logic for residue, chain, structure, motif, and optional matrix-level pLDDT from mmCIF files.
- `batch_cif_to_pdb.py` logic for batch mmCIF-to-PDB conversion.

---

## 11. Troubleshooting

### No files found

Check that the input folder contains `.cif`, `.mmcif`, `.pdb`, `.json`, or `.npz` files.

```bash
find predictions -type f | head
```

### No pLDDT parsed from mmCIF

The coordinate file must contain `_atom_site.B_iso_or_equiv` values. Some tools may not store confidence in the B-factor field. In that case, rely on JSON/NPZ confidence outputs.

### JSON metrics are empty

The JSON file may use a new schema. Add aliases to `JSON_METRIC_ALIASES` in the script.

### Wrong model grouping

If files are grouped incorrectly, rename files to share a common stem, for example:

```text
CEH1_model_0.cif
confidence_CEH1_model_0.json
pae_CEH1_model_0.npz
```

### Excel was not created

Install `pandas` and `openpyxl`, then rerun with `--write-excel`.

### PDB conversion failed

Install Biopython and rerun with `--convert-pdb`.

---

## 12. Citation statement for methods

When using this in a manuscript, cite the upstream prediction platforms and confidence metrics according to the version used. This helper pipeline is a post-processing workflow and does not generate structure predictions directly.
