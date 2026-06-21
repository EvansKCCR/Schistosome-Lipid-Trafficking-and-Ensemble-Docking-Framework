# Schistosome Lipid-Protein Phylogeny and Motif Annotation Pipeline

This repository provides an integrated command-line workflow for downstream analysis of lipid metabolism-associated proteins identified from Pfam/HMMER screening in *Schistosoma mansoni*, *S. haematobium*, and *S. japonicum*.

The pipeline accepts combined or family-specific FASTA files from the Pfam extraction step and integrates:

1. topology and signal peptide evidence from SignalP-6.0, DeepTMHMM, optional CCTOP, and Kyte-Doolittle fallback rules;
2. full-length multiple sequence alignment with MAFFT;
3. maximum-likelihood phylogeny with IQ-TREE/ModelFinder, UFBoot, and SH-aLRT;
4. amino-acid conservation matrices and optional sequence logos;
5. curated motif scanning for CD36RP, NPC2, NPC1, CEH/α/β-hydrolase, and built-in diagnostic motifs;
6. α/β-hydrolase candidate triage into ecto-CEH, type-I membrane CEH, cytosolic HSL-like lipase, cytosolic α/β-hydrolase, or unresolved classes;
7. per-sequence strip diagrams showing signal peptides, TM helices, GPI-like regions, and motif positions.

---

## 1. Repository contents

```text
schisto_phylo_motif_pipeline.py      # main integrated pipeline
README.md                            # this documentation
README.me                            # same content for compatibility with your naming convention
requirements.txt                     # Python dependencies
run_integrated_pipeline.sh           # example command template
motif_library_normalized.tsv         # normalized motif library extracted from Motif_library.xlsx
```

External tools are required only when `--run-phylogeny` is used:

- MAFFT v7 or later
- IQ-TREE 2 preferred, or IQ-TREE 1 compatible binary named `iqtree`

SignalP-6.0, DeepTMHMM, and CCTOP can be run externally through their web/CLI interfaces. Their exported results are then supplied to this pipeline for integrated annotation.

---

## 2. Installation

Create a clean environment:

```bash
conda create -n schisto_phylo_motif python=3.10 -y
conda activate schisto_phylo_motif
pip install -r requirements.txt
```

Install phylogeny tools if you plan to infer trees:

```bash
conda install -c bioconda mafft iqtree -y
```

---

## 3. Required inputs

### 3.1 Family-specific FASTA files

Use the FASTA outputs from the lipid-gene extraction/Pfam screening step. Each file should contain one candidate family, for example:

```text
FASTA/FABP.fasta
FASTA/CD36.fasta
FASTA/NPC1.fasta
FASTA/NPC2.fasta
FASTA/HSL_N.fasta
FASTA/Abhydrolase_3.fasta
```

The filename stem is used as the family name in outputs.

### 3.2 Motif library

The pipeline accepts the uploaded motif library as Excel, CSV, or TSV:

```bash
--motif-library Motif_library.xlsx
```

The parser detects tables with columns equivalent to:

```text
Motif ID | Consensus Motif | Regex-Compatible Pattern | Likely Functional Region | Conservation
```

The supplied library is normalized to `motif_library_normalized.tsv`. Built-in diagnostic motifs are added unless `--no-builtin-motifs` is specified:

```text
HYDROLASE_GxSxG       G[A-Z]S[A-Z]G
HYDROLASE_HGG_like    HGGG?|HGS|HGA
N_GLYCO               N[^P][ST]
```

### 3.3 SignalP-6.0 output

Expected SignalP file: tab-delimited `prediction_results.txt` with sequence ID, prediction, and cleavage-site position.

```bash
--signalp prediction_results.txt
```

When SignalP is absent, the pipeline uses a conservative N-terminal fallback rule based on positive N-region content, an N-terminal hydrophobic h-region, and a simple cleavage-zone rule.

### 3.4 DeepTMHMM output

The pipeline supports both of these DeepTMHMM formats:

1. Markdown report containing the three-line topology block:

```text
>ID | SP+TM
SEQUENCE
SSSSSSOOOOOOOMMMMMMIIIIII
```

2. GFF/GFF3-like records with `signal` and `TMhelix` features.

Use:

```bash
--deeptmhmm deeptmhmm_results.md
```

DeepTMHMM topology is prioritized over the Kyte-Doolittle fallback.

### 3.5 Optional CCTOP output

CCTOP can be supplied as a TSV/text file containing sequence IDs and TM ranges. DeepTMHMM is used first; CCTOP is used as secondary evidence; Kyte-Doolittle is used only when neither external topology source is present.

```bash
--cctop cctop_results.tsv
```

---

## 4. Basic usage

### Annotation-only run

Use this when you already have topology outputs and want motif annotation, hydrolase triage, tables, strip diagrams, and an Excel workbook:

```bash
python schisto_phylo_motif_pipeline.py \
  --input-dir FASTA \
  --motif-library Motif_library.xlsx \
  --signalp prediction_results.txt \
  --deeptmhmm deeptmhmm_results.md \
  --outdir results_phylo_motif \
  --threads 8 \
  --pdf-report
```

### Full phylogeny + annotation run

Use this when MAFFT and IQ-TREE are installed:

```bash
python schisto_phylo_motif_pipeline.py \
  --input-dir FASTA \
  --motif-library Motif_library.xlsx \
  --signalp prediction_results.txt \
  --deeptmhmm deeptmhmm_results.md \
  --outdir results_phylo_motif \
  --threads 8 \
  --run-phylogeny \
  --bootstrap 1000 \
  --alrt 1000 \
  --pdf-report
```

### Multiple FASTA files without an input directory

```bash
python schisto_phylo_motif_pipeline.py \
  --fasta CD36.fasta NPC1.fasta NPC2.fasta Abhydrolase_3.fasta \
  --motif-library Motif_library.xlsx \
  --deeptmhmm deeptmhmm_results.md \
  --outdir results_phylo_motif
```

---

## 5. Outputs

The pipeline creates the following directory structure:

```text
results_phylo_motif/
├── schisto_phylo_motif_results.xlsx
├── integrated_evidence_report.pdf
├── tables/
│   ├── all_sequence_features.tsv
│   ├── all_motif_hits.tsv
│   ├── family_summary.tsv
│   ├── motif_library_normalized.tsv
│   └── phylogeny_outputs.tsv
├── strips/
│   └── <family>/<sequence_id>.png
└── phylogeny/
    └── <family>/
        ├── <family>.fasta
        ├── alignment.mafft.fasta
        ├── alignment_frequency_matrix.csv
        ├── sequence_logo.pdf
        ├── iqtree_ml.treefile
        ├── phylogeny_midpoint_ladderized_collapsed.nwk
        ├── phylogeny_final.pdf
        └── phylogeny_final.svg
```

### `all_sequence_features.tsv`

Main integrated table. Important columns include:

- `SeqID`
- `Family`
- `Length`
- `SignalPeptide`
- `SP_Cleavage`
- `Signal_Source`
- `DeepTMHMM_Label`
- `TM_Count`
- `TM_Ranges`
- `TM_Source`
- `KD_TM_Ranges`
- `Nterm_KD_30aa`
- `Max_Nterm_KD`
- `GPI_like`
- `GPI_omega`
- `Motif_Hit_Count`
- `Motif_IDs`
- `GxSxG_count`
- `HGG_like_count`
- `N_glyco_count`
- `Hydrolase_Classification`

### `all_motif_hits.tsv`

One row per motif occurrence. Motif coordinates are reported as 1-based inclusive positions:

- `family`
- `seq_id`
- `motif_family`
- `motif_id`
- `consensus`
- `regex`
- `functional_region`
- `conservation`
- `start_1based`
- `end_1based`
- `matched_sequence`

### `alignment_frequency_matrix.csv`

One file per family. It contains residue counts for the 20 amino acids at every MSA column, plus:

- gap count;
- non-gap count;
- maximum residue frequency;
- Shannon entropy.

---

## 6. Hydrolase classification rules

The α/β-hydrolase triage module applies explicit rules:

| Classification | Rule |
|---|---|
| Ecto CEH, secreted/GPI-like | Signal peptide present, GPI-like C-terminal region present, and ≤1 TM helix |
| Ecto CEH, type-I single-pass | Signal peptide present and exactly 1 TM helix |
| Ecto CEH, secreted | Signal peptide present and no TM helix |
| Cytosolic HSL-like neutral lipase | No signal peptide, no TM helix, HGGG-like motif, and GxSxG motif |
| Cytosolic α/β-hydrolase candidate | No signal peptide, no TM helix, and GxSxG motif |
| α/β-hydrolase motif-positive, topology unresolved | Hydrolase motifs present, but topology conflicts or is incomplete |
| Unassigned | Criteria above not met |

---

## 7. Heuristic parameters

The following thresholds are defined in the script and can be adjusted by command-line options where relevant:

| Parameter | Default | Meaning |
|---|---:|---|
| `--tm-window` | 19 | Kyte-Doolittle sliding window for TM fallback |
| `--tm-threshold` | 1.60 | Mean KD score for hydrophobic TM segment call |
| N-terminal hydrophobic scan | first 70 aa | used for fallback signal peptide screening |
| signal h-region threshold | 1.60 | minimum mean KD score for N-terminal h-region |
| GPI tail scan | last 50 aa | C-terminal window searched for GPI-like signal |
| GPI hydrophobic fraction | 0.70 | minimum hydrophobic fraction in C-terminal tail |
| GPI tail length | 18-25 aa | hydrophobic tail length searched |
| displayed tree support | ≥70 | internal support label shown on tree |
| collapsed branch support | <50 | low-confidence internal branches collapsed |

---

## 8. Phylogeny details

For each FASTA file/family, the pipeline performs:

1. MAFFT v7 full-length sequence alignment using `--auto`.
2. IQ-TREE maximum-likelihood inference.
3. ModelFinder automatic model selection with `-m MFP`.
4. UFBoot with 1000 replicates by default.
5. SH-aLRT with 1000 replicates by default.
6. Midpoint rooting and ladderization.
7. Support parsing from IQ-TREE internal labels.
8. Display of support values ≥70.
9. Collapse of branches with support <50.
10. Optional keyword-based clade annotation.

You can provide custom clade keywords:

```bash
--clade-keywords 'HSL/LIPE=LIPS,LIPE,HSL;RBBP9=RBBP9;CES-like=CES,EST,CEH;NPC1=NPC1;NPC2=NPC2;CD36RP=CD36,SCARB,SRB'
```

---

## 9. Recommended manuscript wording

The integrated pipeline takes family-specific FASTA files derived from Pfam-based lipid protein screening and performs topology integration, motif annotation, conservation analysis, and phylogenetic inference. Signal peptide and topology evidence are imported from SignalP-6.0, DeepTMHMM, and optionally CCTOP, with Kyte-Doolittle hydropathy rules used only as fallback evidence. Full-length sequences are aligned using MAFFT, and maximum-likelihood phylogenies are inferred with IQ-TREE using ModelFinder model selection, 1000 UFBoot replicates, and 1000 SH-aLRT replicates. Trees are midpoint-rooted, ladderized, and simplified by collapsing internal branches below 50% support; support labels are displayed for nodes with support values ≥70%. Curated regex motif libraries are used to annotate lipid-binding, sterol-trafficking, cysteine-rich, aromatic/hydrophobic, and α/β-hydrolase motifs, with motif positions reported using 1-based coordinates. α/β-hydrolase candidates are further classified using signal peptide state, TM helix count, GPI-like C-terminal features, GxSxG nucleophile loop presence, and HGG-like oxyanion motif status.

---

## 10. Troubleshooting

### No topology calls were assigned

Check whether sequence IDs in the FASTA exactly match IDs in SignalP/DeepTMHMM/CCTOP outputs. The parser uses exact IDs.

### MAFFT or IQ-TREE is missing

Run without `--run-phylogeny` for motif/topology annotation only, or install the tools with Conda:

```bash
conda install -c bioconda mafft iqtree -y
```

### DeepTMHMM Markdown is not parsed

The parser supports both 3-line topology blocks and GFF/GFF3-like records. For Markdown, each record should look like:

```text
>ID | SP
SEQUENCE
SSSSSSSSSOOOOOOOOO
```

### Invalid motif regex

Invalid regular expressions are skipped with a warning. Inspect:

```text
results_phylo_motif/tables/motif_library_normalized.tsv
```

---

## 11. Citation note

When using the pipeline results in a manuscript, cite the upstream tools used for each analysis step: MAFFT, IQ-TREE/ModelFinder/UFBoot, SignalP-6.0, DeepTMHMM, and CCTOP, as applicable.
