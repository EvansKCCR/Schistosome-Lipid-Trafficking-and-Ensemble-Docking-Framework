# Identification of Lipid Metabolism-Associated Proteins in Schistosomes

This repository contains a reproducible pipeline for identifying candidate lipid metabolism-associated proteins in *Schistosoma mansoni*, *Schistosoma haematobium*, and *Schistosoma japonicum* proteomes using Pfam domain annotation from HMMER `hmmscan` output.

The pipeline screens proteomes for Pfam domains associated with lipid binding, lipid transport, sterol handling, fatty-acid transport, and lipid enzymatic processing. It generates family-specific FASTA files, species-level Excel workbooks, and combined multi-species candidate tables.

---

## 1. Pipeline overview

```text
Proteome FASTA files
        │
        ▼
HMMER hmmscan against Pfam-A.hmm
        │
        ▼
HMMER domtblout parsing
        │
        ▼
Independent domain E-value filtering ≤ 1e-5
        │
        ▼
Pfam lipid-domain family classification
        │
        ▼
Candidate protein sequence extraction
        │
        ▼
Family-specific FASTA files + Excel workbooks + summary tables
```

---

## 2. Input data

Proteome datasets for the following schistosome species should be provided in FASTA format:

- *Schistosoma mansoni*
- *Schistosoma haematobium*
- *Schistosoma japonicum*

Proteomes may be downloaded from WormBase ParaSite:

```text
https://parasite.wormbase.org/index.html
```

Protein domain annotation is performed using the Pfam-A HMM database:

```text
Pfam-A.hmm
```

Before running `hmmscan`, ensure that the Pfam HMM database has been indexed:

```bash
hmmpress Pfam-A.hmm
```

---

## 3. Software requirements

### Required command-line software

- HMMER, tested with HMMER v3.1b2 or later
- Python 3.8 or later

### Required Python packages

Install the Python dependencies with:

```bash
pip install -r requirements.txt
```

The required packages are:

```text
biopython>=1.80
pandas>=1.5
openpyxl>=3.1
```

---

## 4. Repository contents

```text
schisto_lipid_pipeline.py   Main reusable Python pipeline
samples.tsv                 Example sample sheet for the three schistosome proteomes
lipid_domains.tsv           Lipid-associated Pfam domain definitions
requirements.txt            Python package requirements
run_pipeline.sh             Example shell script for a full hmmscan + extraction run
README.md                   Pipeline documentation
```

---

## 5. Pfam domain families screened

The default lipid-associated Pfam families are:

| Family | Pfam accession | Functional category |
|---|---:|---|
| FABP | PF00061 | Fatty acid-binding protein |
| CD36 | PF01130 | CD36/fatty acid translocase-like domain |
| ABC | PF00005 | ABC transporter ATP-binding domain |
| NPC1 | PF12349 | NPC1 sterol-sensing domain-associated family |
| NPC2 | PF02221 | Niemann-Pick type C2 lipid-binding domain |
| FATP | PF02259 | Fatty acid transport protein/AMP-binding-associated family |
| HSL_N | PF06350 | Hormone-sensitive lipase N-terminal domain |
| Abhydrolase_3 | PF07859 | Alpha/beta hydrolase fold family |

The pipeline internally normalizes Pfam IDs by removing version suffixes. Therefore, both `PF00061` and `PF00061.30` are treated as the same Pfam family. This makes the workflow more robust across Pfam releases.

---

## 6. Prepare the sample sheet

Create a tab-separated file named `samples.tsv` with the following columns:

```text
species	fasta	domtblout
Schistosoma_mansoni	Smansoni_proteins.fa	domtblout.txt
Schistosoma_haematobium	Shaematobium_proteins.fa	shdomtblout.txt
Schistosoma_japonicum	Sjaponicum_proteins.fa	Sjdomtblout.txt
```

### Column descriptions

| Column | Description |
|---|---|
| `species` | Species or sample name used in output files |
| `fasta` | Protein FASTA file for that species |
| `domtblout` | Existing HMMER `--domtblout` file; optional if using `--run-hmmscan` |

Relative paths are interpreted relative to the location of `samples.tsv`.

---

## 7. Option A: Run the full pipeline, including hmmscan

Use this mode when you have FASTA proteomes and `Pfam-A.hmm`, but have not yet generated `domtblout` files.

```bash
python schisto_lipid_pipeline.py \
  --samples samples.tsv \
  --pfam-hmm Pfam-A.hmm \
  --domains lipid_domains.tsv \
  --outdir lipid_pipeline_results \
  --evalue 1e-5 \
  --cpu 8 \
  --run-hmmscan
```

This executes a command equivalent to the following for each proteome:

```bash
hmmscan --cpu 8 --domtblout <species>.domtblout.txt Pfam-A.hmm <species>_proteins.fa > <species>.hmmscan.log
```

---

## 8. Option B: Parse existing hmmscan domtblout files only

Use this mode if `hmmscan` has already been run and the `domtblout` files are listed in `samples.tsv`.

```bash
python schisto_lipid_pipeline.py \
  --samples samples.tsv \
  --domains lipid_domains.tsv \
  --outdir lipid_pipeline_results \
  --evalue 1e-5
```

In this mode, `--pfam-hmm` and `--run-hmmscan` are not required.

---

## 9. Domain filtering and candidate definition

The pipeline parses HMMER `--domtblout` files and uses the independent domain E-value column for filtering.

A protein is retained as a lipid metabolism-associated candidate if:

1. It contains at least one target Pfam domain listed in `lipid_domains.tsv`.
2. The independent domain E-value is less than or equal to `1e-5`.
3. The protein identifier in the `domtblout` file matches a sequence identifier in the corresponding FASTA file.

Each retained protein is annotated with:

- Species
- Protein ID
- Lipid-related family
- Pfam accession
- Versioned Pfam accession reported by HMMER
- Independent domain E-value
- Protein sequence length

Proteins with multiple qualifying lipid-associated domains may appear in more than one family worksheet. Duplicate records are removed.

---

## 10. Output structure

After a successful run, the output directory will contain species-specific and combined results:

```text
lipid_pipeline_results/
├── Schistosoma_mansoni/
│   ├── FASTA/
│   │   ├── ABC.fasta
│   │   ├── FABP.fasta
│   │   └── ...
│   ├── Schistosoma_mansoni_Lipid_Gene_Candidates.xlsx
│   ├── Schistosoma_mansoni_candidate_counts_by_family.csv
│   └── hmmer/
│       ├── Schistosoma_mansoni.domtblout.txt
│       └── Schistosoma_mansoni.hmmscan.log
├── Schistosoma_haematobium/
│   └── ...
├── Schistosoma_japonicum/
│   └── ...
└── combined/
    ├── ALL_SPECIES_Lipid_Gene_Candidates.csv
    ├── ALL_SPECIES_Lipid_Gene_Candidates.xlsx
    └── ALL_SPECIES_candidate_counts_by_family.csv
```

### Main outputs

| Output | Description |
|---|---|
| `FASTA/<Family>.fasta` | Deduplicated candidate protein sequences for each lipid-related family |
| `<Species>_Lipid_Gene_Candidates.xlsx` | Species-level Excel workbook with per-family sheets and `ALL_CANDIDATES` |
| `<Species>_candidate_counts_by_family.csv` | Count of unique candidate proteins per family |
| `combined/ALL_SPECIES_Lipid_Gene_Candidates.xlsx` | Combined workbook across all species |
| `combined/ALL_SPECIES_Lipid_Gene_Candidates.csv` | Combined candidate table across all species |
| `combined/ALL_SPECIES_candidate_counts_by_family.csv` | Multi-species family count summary |

---

## 11. Validation and error handling

The pipeline stops with an informative error if:

- A FASTA file is missing.
- A required `domtblout` file is missing in parse-only mode.
- `hmmscan` is requested but HMMER is not available on `PATH`.
- `--pfam-hmm` is missing when `--run-hmmscan` is used.
- No significant lipid-associated Pfam hits are detected.
- Protein IDs in `domtblout` do not match the corresponding FASTA sequence IDs.

If no candidates are detected, check:

1. Whether `Pfam-A.hmm` was correctly downloaded and indexed with `hmmpress`.
2. Whether `hmmscan` completed successfully.
3. Whether the correct FASTA file was used for each species.
4. Whether Pfam accession definitions match the current Pfam release.
5. Whether the E-value cutoff is too stringent for the dataset.

---

## 12. Notes on improvements over the original species-specific scripts

The original workflow used one separate script per schistosome species. This pipeline consolidates those scripts into a single reusable workflow with the following improvements:

- One sample sheet controls all species.
- `hmmscan` can be run automatically or skipped if `domtblout` files already exist.
- Pfam accessions are matched without requiring exact version suffixes.
- Species-specific outputs are written to separate directories to prevent overwriting.
- Combined cross-species outputs are generated automatically.
- Protein IDs are validated against the FASTA file before sequence extraction.
- The *S. japonicum* output naming issue from the original species-specific script is avoided by deriving output names from the `species` column.

---

## 13. Reproducible command used in the manuscript methods

A concise methods-compatible command is:

```bash
python schisto_lipid_pipeline.py \
  --samples samples.tsv \
  --pfam-hmm Pfam-A.hmm \
  --domains lipid_domains.tsv \
  --outdir lipid_pipeline_results \
  --evalue 1e-5 \
  --cpu 8 \
  --run-hmmscan
```

Candidate lipid metabolism-associated proteins were defined as proteins containing one or more lipid-related Pfam domains with independent domain E-value ≤ `1e-5`.
