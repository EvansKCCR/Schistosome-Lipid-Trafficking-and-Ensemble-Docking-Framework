"""Run the complete flexdock_v3 representative-selection workflow."""

from __future__ import annotations

import argparse
from pathlib import Path

from .build_complexes import build_complexes
from .cluster_his508 import cluster_his508
from .cluster_ligands import cluster_ligands
from .common import configure_logging, get_logger
from .extract_flex_modes import extract_flex_modes
from .select_representatives import select_representatives


PACKAGE_DIR = Path(__file__).resolve().parent
DATA_DIR = PACKAGE_DIR.parent


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--flex-out",
        type=Path,
        default=DATA_DIR / "rank1_confidence-2.10_flex_out.pdbqt",
    )
    parser.add_argument(
        "--receptor",
        type=Path,
        default=DATA_DIR / "receptor_rigid.pdbqt",
    )
    parser.add_argument(
        "--ligand-sdf",
        type=Path,
        default=DATA_DIR / "rank1_confidence-2.10.sdf",
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=DATA_DIR / "results",
    )
    parser.add_argument("--expected-ligand-atoms", type=int, default=47)
    parser.add_argument(
        "--his-method",
        choices=("agglomerative", "dbscan"),
        default="agglomerative",
    )
    parser.add_argument(
        "--his-cutoff",
        type=float,
        default=0.75,
        help="HIS508 CA-centered side-chain RMSD cutoff in Angstrom",
    )
    parser.add_argument("--his-min-samples", type=int, default=2)
    parser.add_argument(
        "--ligand-cutoff",
        type=float,
        default=2.0,
        help="Ligand Kabsch heavy-atom RMSD cutoff in Angstrom",
    )
    parser.add_argument("--verbose", action="store_true")
    return parser


def validate_inputs(paths: list[Path]) -> None:
    missing = [path for path in paths if not path.is_file()]
    if missing:
        raise FileNotFoundError(
            "Required input files are missing:\n"
            + "\n".join(f"  - {path}" for path in missing)
        )


def run_pipeline(args: argparse.Namespace) -> None:
    results_dir = args.results_dir.resolve()
    results_dir.mkdir(parents=True, exist_ok=True)
    configure_logging(results_dir / "flexdock_v3.log", args.verbose)
    logger = get_logger()

    flex_out = args.flex_out.resolve()
    receptor = args.receptor.resolve()
    ligand_sdf = args.ligand_sdf.resolve()
    validate_inputs([flex_out, receptor, ligand_sdf])

    logger.info("Starting flexdock_v3 pipeline")
    logger.info("Flexible output: %s", flex_out)
    logger.info("Rigid receptor: %s", receptor)
    logger.info("Ligand topology: %s", ligand_sdf)
    logger.info("Results directory: %s", results_dir)

    logger.info("Stage 1/5: extracting flexible docking modes")
    extract_flex_modes(
        flex_out,
        results_dir,
        expected_ligand_atoms=args.expected_ligand_atoms,
    )

    logger.info("Stage 2/5: clustering HIS508 conformations")
    cluster_his508(
        results_dir / "pose_data.pkl",
        results_dir,
        method=args.his_method,
        cutoff=args.his_cutoff,
        min_samples=args.his_min_samples,
    )

    logger.info("Stage 3/5: clustering ligand poses")
    cluster_ligands(
        results_dir / "pose_data.pkl",
        results_dir / "his508_clusters.csv",
        results_dir,
        cutoff=args.ligand_cutoff,
    )

    logger.info("Stage 4/5: selecting cluster representatives")
    select_representatives(
        results_dir / "ligand_clusters.csv",
        results_dir / "his508_clusters.csv",
        results_dir / "pose_summary.csv",
        results_dir,
    )

    logger.info("Stage 5/5: building representative complexes")
    manifest = build_complexes(
        receptor,
        ligand_sdf,
        results_dir / "pose_data.pkl",
        results_dir / "representative_modes.csv",
        results_dir,
    )

    logger.info(
        "Pipeline complete: %d representative complexes generated",
        len(manifest),
    )


def main() -> None:
    args = build_parser().parse_args()
    run_pipeline(args)


if __name__ == "__main__":
    main()
