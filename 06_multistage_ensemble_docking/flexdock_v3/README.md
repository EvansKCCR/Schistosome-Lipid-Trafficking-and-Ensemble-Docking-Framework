# flexdock_v3

`flexdock_v3` selects representative CEH-ligand complexes from an AutoDock
Vina flexible-receptor output. It parses every Vina MODEL, clusters the
flexible HIS508 side chain, clusters ligand poses within each HIS508 cluster,
selects the lowest-scoring pose from every combined cluster, and builds
protein-ligand complexes for subsequent MD system preparation.

## Input files

The default paths assume the package is in the same `data_analysis` directory
as:

- `rank1_confidence-2.10_flex_out.pdbqt`
- `receptor_rigid.pdbqt`
- `rank1_confidence-2.10.sdf`

The supplied flexible output contains 50 MODELs. Each MODEL must contain:

- exactly 47 or 49 `UNL` ligand atoms (change when appropriate);
- one flexible HIS508 block;
- HIS508 atoms `CA`, `CB`, `CG`, `ND1`, `CD2`, `CE1`, and `NE2`;
- `REMARK VINA RESULT` score and RMSD fields;
- a complete `REMARK SMILES IDX` permutation.

## Installation

Conda is recommended because it provides consistent RDKit builds on Linux,
WSL, and Windows:

```bash
conda create -n flexdock_v3 -c conda-forge \
  python=3.11 numpy pandas scipy scikit-learn rdkit matplotlib seaborn
conda activate flexdock_v3
```

Alternatively:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r flexdock_v3/requirements.txt
```

## Run the pipeline

From the directory containing the `flexdock_v3` package:

```bash
python -m flexdock_v3.run_pipeline
```

Explicit paths:

```bash
python -m flexdock_v3.run_pipeline \
  --flex-out rank1_confidence-2.10_flex_out.pdbqt \
  --receptor receptor_rigid.pdbqt \
  --ligand-sdf rank1_confidence-2.10.sdf \
  --results-dir results \
  --his-method agglomerative \
  --his-cutoff 0.75 \
  --ligand-cutoff 2.0
```

The ligand cutoff is 2.0 Å by default. HIS508 conformations are translated so
that CA is at the origin, then the ordered side-chain coordinates
`CB, CG, ND1, CD2, CE1, NE2` are compared. No backbone-N chi1 calculation is
used.

## Individual stages

```bash
python -m flexdock_v3.extract_flex_modes rank1_confidence-2.10_flex_out.pdbqt
python -m flexdock_v3.cluster_his508
python -m flexdock_v3.cluster_ligands
python -m flexdock_v3.select_representatives
python -m flexdock_v3.build_complexes receptor_rigid.pdbqt rank1_confidence-2.10.sdf
```

## Outputs

The `results` directory contains:

- `pose_data.pkl`
- `pose_summary.csv`
- `his508_clusters.csv`
- `ligand_clusters.csv`
- `representative_modes.csv`
- `score_distribution.png`
- `his508_clustering.png`
- `ligand_rmsd_clustering.png`
- `flexdock_v3.log`
- `representative_complexes/complex_manifest.csv`
- `representative_complexes/cluster_X_Y.pdb`
- `representative_complexes/cluster_X_Y_ligand.sdf`

The companion ligand SDF files preserve the exact atom order, stereochemistry,
aromaticity, and bond orders from the supplied topology template. The complex
PDB files contain the same SDF-derived atom order and `CONECT` records. PDB
does not have a complete aromatic-bond representation, so retain the companion
SDF during force-field parameterization.

## Complex construction details

The supplied rigid receptor stores the HIS508 backbone atoms while the flexible
Vina MODEL stores CA and side-chain atoms. The builder removes rigid HIS508,
merges its backbone with the selected flexible coordinates, and inserts one
complete HIS508 residue.

Ligand topology is never inferred from PDBQT. The builder:

1. loads the supplied SDF with RDKit;
2. parses the MODEL's `REMARK SMILES` and `REMARK SMILES IDX` records;
3. graph-matches the SMILES to the SDF;
4. determines the index-map direction by exact element validation;
5. transfers docked coordinates into SDF atom order.

The generated structures are starting complexes for MD preparation. They
still require the usual protonation review, ligand parameterization, force
field assignment, solvation, and minimization.
