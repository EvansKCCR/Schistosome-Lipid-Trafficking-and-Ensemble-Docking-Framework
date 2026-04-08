# Schistosome-Lipid-Trafficking-and-Ensemble-Docking-Framework
This repository contains a reproducible computational pipeline for analyzing host lipid acquisition and sterol trafficking in Schistosoma spp. The framework integrates sequence analysis, structural modeling, ensemble docking, clustering, and convergence-based metrics to identify candidate lipid-handling mechanisms.
Key Features
Multiple sequence alignment and phylogenetic reconstruction (MAFFT + IQ-TREE)
Family-aware motif detection and classification
Ensemble docking (DiffDock → SMINA/Vina refinement)
Dual clustering:
Geometry-based (RMSD)
Score-space (affinity + RMSD)
Convergence metrics:
Dominance fraction
Shannon entropy
Structural validation and docking readiness scoring
Core Concept

This work introduces a convergence-based framework showing that:

Binding geometry and energetic stabilization are decoupled in lipid-binding systems.

Installation
git clone https://github.com/EvansKCCR/schistosome-lipid-trafficking.git
cd schistosome-lipid-trafficking
pip install -r requirements.txt
Usage
1. Sequence analysis
python scripts/alignment_phylogeny_v2.py
2. Motif detection and classification
python scripts/schistosome_sterol_transport_pipeline.py --families all
3. Hydrolase classification
python scripts/classify_hydrolases.py input.fasta \
  --signalp_tsv prediction_results.txt \
  --deeptmhmm_md deeptmhmm_results.md
4. Docking analysis
python scripts/cluster_geometry_poses.py
python scripts/cluster_score_space_poses.py
Outputs
Phylogenetic trees (PDF/SVG)
Sequence logos
Motif annotation tables
Docking clustering results
Heatmaps and multi-panel figures
Reproducibility

All analyses are fully script-driven and reproducible from input FASTA to final figures.

Citation

If you use this work, please cite:

Adu EA. Structure-resolved model of host lipid capture and sterol trafficking at the schistosome tegument.
