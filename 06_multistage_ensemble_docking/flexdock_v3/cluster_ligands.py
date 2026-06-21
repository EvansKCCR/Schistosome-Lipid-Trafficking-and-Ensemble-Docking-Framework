"""Cluster ligand poses by Kabsch-aligned heavy-atom RMSD within HIS508 clusters."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.cluster import AgglomerativeClustering

from .common import (
    cluster_medoids,
    configure_logging,
    condensed_pairwise_rmsd,
    get_logger,
    load_pickle,
    pose_lookup,
    stable_relabel,
)


def ligand_coordinates(pose: dict) -> np.ndarray:
    atoms = sorted(pose["ligand_atoms"], key=lambda atom: int(atom["serial"]))
    coordinates = np.vstack([np.asarray(atom["coord"], dtype=float) for atom in atoms])
    elements = [str(atom["element"]) for atom in atoms]
    if any(element.upper() == "H" for element in elements):
        raise ValueError(f"MODEL {pose['mode']}: ligand unexpectedly contains hydrogen atoms")
    return coordinates


def agglomerative_cutoff(distance_matrix: np.ndarray, cutoff: float) -> np.ndarray:
    if len(distance_matrix) == 1:
        return np.array([1], dtype=int)
    model = AgglomerativeClustering(
        n_clusters=None,
        metric="precomputed",
        linkage="complete",
        distance_threshold=cutoff,
    )
    return model.fit_predict(distance_matrix)


def plot_ligand_clustering(
    all_coordinates: dict[int, np.ndarray],
    assignments: pd.DataFrame,
    output_path: Path,
) -> None:
    order_df = assignments.sort_values(
        ["his508_cluster", "ligand_cluster", "vina_score", "mode"]
    )
    modes = order_df["mode"].astype(int).tolist()
    coordinates = [all_coordinates[mode] for mode in modes]
    distance_matrix = condensed_pairwise_rmsd(coordinates, align=True)

    sns.set_theme(style="white", context="paper")
    figure, axes = plt.subplots(
        1,
        2,
        figsize=(12.5, 5.2),
        gridspec_kw={"width_ratios": [1.6, 0.8]},
    )
    sns.heatmap(
        distance_matrix,
        cmap="mako",
        square=True,
        xticklabels=modes,
        yticklabels=modes,
        cbar_kws={"label": "Kabsch RMSD (Å)", "shrink": 0.78},
        ax=axes[0],
    )
    axes[0].set_title("Ligand heavy-atom RMSD")
    axes[0].set_xlabel("Vina mode")
    axes[0].set_ylabel("Vina mode")
    axes[0].tick_params(axis="x", labelrotation=90, labelsize=6)
    axes[0].tick_params(axis="y", labelrotation=0, labelsize=6)

    cluster_sizes = (
        assignments.groupby(["his508_cluster", "ligand_cluster"], as_index=False)
        .size()
        .rename(columns={"size": "pose_count"})
    )
    cluster_sizes["combined_cluster"] = cluster_sizes.apply(
        lambda row: f"H{int(row.his508_cluster)}-L{int(row.ligand_cluster)}",
        axis=1,
    )
    cluster_sizes = cluster_sizes.sort_values("pose_count", ascending=True)
    sns.barplot(
        data=cluster_sizes,
        x="pose_count",
        y="combined_cluster",
        color="#2A9D8F",
        ax=axes[1],
    )
    axes[1].set_xlabel("Pose count")
    axes[1].set_ylabel("Combined cluster")
    axes[1].set_title("Cluster occupancy")
    figure.tight_layout()
    figure.savefig(output_path, dpi=350, bbox_inches="tight")
    plt.close(figure)


def cluster_ligands(
    pose_data_path: Path,
    his_clusters_path: Path,
    results_dir: Path,
    cutoff: float = 2.0,
) -> pd.DataFrame:
    logger = get_logger()
    pose_data = load_pickle(pose_data_path)
    poses = pose_lookup(pose_data)
    his_assignments = pd.read_csv(his_clusters_path)

    pose_modes = set(poses)
    cluster_modes = set(his_assignments["mode"].astype(int))
    if pose_modes != cluster_modes:
        raise ValueError(
            "HIS508 assignments and pose_data contain different mode sets: "
            f"pose-only={sorted(pose_modes - cluster_modes)}, "
            f"csv-only={sorted(cluster_modes - pose_modes)}"
        )

    all_coordinates = {mode: ligand_coordinates(pose) for mode, pose in poses.items()}
    reference_shape = next(iter(all_coordinates.values())).shape
    for mode, coordinates in all_coordinates.items():
        if coordinates.shape != reference_shape:
            raise ValueError(
                f"MODEL {mode}: ligand coordinate shape {coordinates.shape} "
                f"differs from {reference_shape}"
            )

    rows = []
    for his_cluster, group in his_assignments.groupby("his508_cluster", sort=True):
        group = group.sort_values("mode").reset_index(drop=True)
        modes = group["mode"].astype(int).to_numpy()
        scores = group["vina_score"].astype(float).to_numpy()
        coordinates = [all_coordinates[int(mode)] for mode in modes]
        distance_matrix = condensed_pairwise_rmsd(coordinates, align=True)
        raw_labels = agglomerative_cutoff(distance_matrix, cutoff)
        labels = stable_relabel(raw_labels, scores, modes)
        medoids = cluster_medoids(distance_matrix, labels)

        for index, mode in enumerate(modes):
            ligand_cluster = int(labels[index])
            medoid_index = medoids[ligand_cluster]
            rows.append(
                {
                    "mode": int(mode),
                    "vina_score": float(scores[index]),
                    "his508_cluster": int(his_cluster),
                    "ligand_cluster": ligand_cluster,
                    "combined_cluster": (
                        f"H{int(his_cluster)}_L{ligand_cluster}"
                    ),
                    "ligand_cluster_size": int(
                        np.count_nonzero(labels == ligand_cluster)
                    ),
                    "ligand_medoid_mode": int(modes[medoid_index]),
                    "ligand_rmsd_to_medoid": float(
                        distance_matrix[index, medoid_index]
                    ),
                    "ligand_is_medoid": bool(index == medoid_index),
                    "rmsd_cutoff": float(cutoff),
                }
            )

        np.save(
            results_dir / f"ligand_rmsd_his_cluster_{int(his_cluster)}.npy",
            distance_matrix,
        )
        logger.debug(
            "HIS cluster %s: %d poses -> %d ligand clusters",
            his_cluster,
            len(group),
            len(set(labels)),
        )

    output = pd.DataFrame(rows).sort_values("mode")
    output_path = results_dir / "ligand_clusters.csv"
    output.to_csv(output_path, index=False)
    plot_ligand_clustering(
        all_coordinates,
        output,
        results_dir / "ligand_rmsd_clustering.png",
    )
    logger.info(
        "Ligand clustering produced %d combined clusters at %.2f Å",
        output["combined_cluster"].nunique(),
        cutoff,
    )
    logger.info("Wrote %s", output_path)
    return output


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pose-data", type=Path, default=Path("results/pose_data.pkl"))
    parser.add_argument(
        "--his-clusters",
        type=Path,
        default=Path("results/his508_clusters.csv"),
    )
    parser.add_argument("--results-dir", type=Path, default=Path("results"))
    parser.add_argument("--cutoff", type=float, default=2.0)
    parser.add_argument("--verbose", action="store_true")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    configure_logging(args.results_dir / "cluster_ligands.log", args.verbose)
    cluster_ligands(
        args.pose_data.resolve(),
        args.his_clusters.resolve(),
        args.results_dir.resolve(),
        cutoff=args.cutoff,
    )


if __name__ == "__main__":
    main()
