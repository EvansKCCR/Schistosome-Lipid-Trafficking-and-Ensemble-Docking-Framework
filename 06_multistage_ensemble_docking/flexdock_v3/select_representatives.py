"""Select the lowest-Vina-score pose from each combined cluster."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from .common import configure_logging, get_logger


def select_representatives(
    ligand_clusters_path: Path,
    his_clusters_path: Path,
    pose_summary_path: Path,
    results_dir: Path,
) -> pd.DataFrame:
    logger = get_logger()
    ligand = pd.read_csv(ligand_clusters_path)
    his = pd.read_csv(his_clusters_path)
    summary = pd.read_csv(pose_summary_path)

    required_ligand = {
        "mode",
        "his508_cluster",
        "ligand_cluster",
        "combined_cluster",
    }
    missing = required_ligand - set(ligand.columns)
    if missing:
        raise ValueError(f"ligand_clusters.csv missing columns: {sorted(missing)}")

    merged = (
        ligand.drop(columns=["vina_score"], errors="ignore")
        .merge(
            his[
                [
                    "mode",
                    "vina_score",
                    "his508_cluster_size",
                    "his508_medoid_mode",
                    "his508_rmsd_to_medoid",
                    "his508_is_medoid",
                ]
            ],
            on="mode",
            how="inner",
            validate="one_to_one",
        )
        .merge(
            summary[["mode", "rmsd_lb", "rmsd_ub"]],
            on="mode",
            how="inner",
            validate="one_to_one",
        )
    )

    if len(merged) != len(ligand):
        raise ValueError("Cluster tables and pose summary contain inconsistent mode sets")

    representatives = (
        merged.sort_values(
            [
                "his508_cluster",
                "ligand_cluster",
                "vina_score",
                "rmsd_lb",
                "mode",
            ],
            ascending=[True, True, True, True, True],
        )
        .groupby(
            ["his508_cluster", "ligand_cluster"],
            sort=True,
            as_index=False,
        )
        .first()
    )
    representatives["representative_rank"] = range(1, len(representatives) + 1)
    representatives = representatives[
        [
            "representative_rank",
            "combined_cluster",
            "his508_cluster",
            "ligand_cluster",
            "mode",
            "vina_score",
            "rmsd_lb",
            "rmsd_ub",
            "his508_cluster_size",
            "ligand_cluster_size",
            "his508_medoid_mode",
            "ligand_medoid_mode",
            "his508_rmsd_to_medoid",
            "ligand_rmsd_to_medoid",
            "his508_is_medoid",
            "ligand_is_medoid",
        ]
    ].sort_values(["his508_cluster", "ligand_cluster"])

    output_path = results_dir / "representative_modes.csv"
    representatives.to_csv(output_path, index=False)
    logger.info(
        "Selected %d representatives from %d poses",
        len(representatives),
        len(ligand),
    )
    logger.info("Wrote %s", output_path)
    return representatives


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--ligand-clusters",
        type=Path,
        default=Path("results/ligand_clusters.csv"),
    )
    parser.add_argument(
        "--his-clusters",
        type=Path,
        default=Path("results/his508_clusters.csv"),
    )
    parser.add_argument(
        "--pose-summary",
        type=Path,
        default=Path("results/pose_summary.csv"),
    )
    parser.add_argument("--results-dir", type=Path, default=Path("results"))
    parser.add_argument("--verbose", action="store_true")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    configure_logging(args.results_dir / "select_representatives.log", args.verbose)
    select_representatives(
        args.ligand_clusters.resolve(),
        args.his_clusters.resolve(),
        args.pose_summary.resolve(),
        args.results_dir.resolve(),
    )


if __name__ == "__main__":
    main()
