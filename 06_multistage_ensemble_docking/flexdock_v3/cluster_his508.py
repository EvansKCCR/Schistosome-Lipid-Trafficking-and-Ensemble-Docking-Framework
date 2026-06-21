"""Cluster flexible HIS508 side-chain conformations."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform
from sklearn.cluster import AgglomerativeClustering, DBSCAN
from sklearn.decomposition import PCA

from .common import (
    HIS_SIDECHAIN_ATOMS,
    cluster_medoids,
    configure_logging,
    coordinates_by_name,
    get_logger,
    load_pickle,
    pose_lookup,
    stable_relabel,
)


def ca_centered_sidechain(pose: dict) -> np.ndarray:
    ca = coordinates_by_name(pose["his508_atoms"], ("CA",))[0]
    sidechain = coordinates_by_name(pose["his508_atoms"], HIS_SIDECHAIN_ATOMS)
    return sidechain - ca


def sidechain_rmsd_matrix(features: list[np.ndarray]) -> np.ndarray:
    n_items = len(features)
    matrix = np.zeros((n_items, n_items), dtype=float)
    for i in range(1, n_items):
        for j in range(i):
            delta = features[i] - features[j]
            rmsd = float(np.sqrt(np.mean(np.sum(delta * delta, axis=1))))
            matrix[i, j] = matrix[j, i] = rmsd
    return matrix


def assign_clusters(
    distance_matrix: np.ndarray,
    method: str,
    cutoff: float,
    min_samples: int,
) -> np.ndarray:
    if len(distance_matrix) == 1:
        return np.array([1], dtype=int)

    if method == "agglomerative":
        model = AgglomerativeClustering(
            n_clusters=None,
            metric="precomputed",
            linkage="average",
            distance_threshold=cutoff,
        )
        return model.fit_predict(distance_matrix)

    labels = DBSCAN(
        eps=cutoff,
        min_samples=min_samples,
        metric="precomputed",
    ).fit_predict(distance_matrix)

    # DBSCAN noise points are retained as singleton clusters.
    next_label = int(labels.max()) + 1
    for index in np.flatnonzero(labels < 0):
        labels[index] = next_label
        next_label += 1
    return labels


def plot_his_clustering(
    features: list[np.ndarray],
    modes: np.ndarray,
    labels: np.ndarray,
    scores: np.ndarray,
    distance_matrix: np.ndarray,
    output_path: Path,
    cutoff: float,
) -> None:
    sns.set_theme(style="whitegrid", context="paper")
    flat_features = np.vstack([feature.reshape(1, -1) for feature in features])
    if len(features) >= 2:
        pca_values = PCA(n_components=2).fit_transform(flat_features)
    else:
        pca_values = np.zeros((1, 2), dtype=float)

    figure, axes = plt.subplots(
        1,
        2,
        figsize=(12, 4.8),
        gridspec_kw={"width_ratios": [1.0, 1.35]},
    )
    palette = sns.color_palette("tab10", n_colors=max(1, len(set(labels))))
    plot_df = pd.DataFrame(
        {
            "PC1": pca_values[:, 0],
            "PC2": pca_values[:, 1],
            "cluster": labels.astype(str),
            "mode": modes,
            "score": scores,
        }
    )
    sns.scatterplot(
        data=plot_df,
        x="PC1",
        y="PC2",
        hue="cluster",
        size="score",
        sizes=(35, 95),
        palette=palette,
        edgecolor="white",
        linewidth=0.5,
        ax=axes[0],
    )
    axes[0].set_title("HIS508 side-chain conformations")
    axes[0].set_xlabel("PCA component 1")
    axes[0].set_ylabel("PCA component 2")
    axes[0].legend(title="Cluster / score", fontsize=7, title_fontsize=8)

    if len(features) >= 2:
        condensed = squareform(distance_matrix, checks=False)
        hierarchy = linkage(condensed, method="average")
        dendrogram(
            hierarchy,
            labels=[str(mode) for mode in modes],
            color_threshold=cutoff,
            above_threshold_color="#777777",
            leaf_rotation=90,
            leaf_font_size=7,
            ax=axes[1],
        )
        axes[1].axhline(cutoff, color="#C44536", linestyle="--", linewidth=1.2)
    axes[1].set_title("Average-linkage dendrogram")
    axes[1].set_xlabel("Vina mode")
    axes[1].set_ylabel("CA-centered side-chain RMSD (Å)")
    figure.tight_layout()
    figure.savefig(output_path, dpi=350, bbox_inches="tight")
    plt.close(figure)


def cluster_his508(
    pose_data_path: Path,
    results_dir: Path,
    method: str = "agglomerative",
    cutoff: float = 0.75,
    min_samples: int = 2,
) -> pd.DataFrame:
    logger = get_logger()
    pose_data = load_pickle(pose_data_path)
    poses = sorted(pose_data["poses"], key=lambda pose: int(pose["mode"]))
    features = [ca_centered_sidechain(pose) for pose in poses]
    modes = np.array([int(pose["mode"]) for pose in poses], dtype=int)
    scores = np.array([float(pose["vina_score"]) for pose in poses], dtype=float)

    distance_matrix = sidechain_rmsd_matrix(features)
    raw_labels = assign_clusters(distance_matrix, method, cutoff, min_samples)
    labels = stable_relabel(raw_labels, scores, modes)
    medoids = cluster_medoids(distance_matrix, labels)

    rows = []
    for index, pose in enumerate(poses):
        label = int(labels[index])
        medoid_index = medoids[label]
        rows.append(
            {
                "mode": int(pose["mode"]),
                "vina_score": float(pose["vina_score"]),
                "his508_cluster": label,
                "his508_cluster_size": int(np.count_nonzero(labels == label)),
                "his508_medoid_mode": int(modes[medoid_index]),
                "his508_rmsd_to_medoid": float(distance_matrix[index, medoid_index]),
                "his508_is_medoid": bool(index == medoid_index),
                "method": method,
                "rmsd_cutoff": float(cutoff),
            }
        )
    output = pd.DataFrame(rows).sort_values("mode")
    output_path = results_dir / "his508_clusters.csv"
    output.to_csv(output_path, index=False)
    np.save(results_dir / "his508_rmsd_matrix.npy", distance_matrix)
    plot_his_clustering(
        features,
        modes,
        labels,
        scores,
        distance_matrix,
        results_dir / "his508_clustering.png",
        cutoff,
    )

    logger.info(
        "HIS508 clustering produced %d clusters using %s at %.3f Å",
        output["his508_cluster"].nunique(),
        method,
        cutoff,
    )
    logger.info("Wrote %s", output_path)
    return output


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pose-data", type=Path, default=Path("results/pose_data.pkl"))
    parser.add_argument("--results-dir", type=Path, default=Path("results"))
    parser.add_argument(
        "--method",
        choices=("agglomerative", "dbscan"),
        default="agglomerative",
    )
    parser.add_argument("--cutoff", type=float, default=0.75)
    parser.add_argument("--min-samples", type=int, default=2)
    parser.add_argument("--verbose", action="store_true")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    configure_logging(args.results_dir / "cluster_his508.log", args.verbose)
    cluster_his508(
        args.pose_data.resolve(),
        args.results_dir.resolve(),
        method=args.method,
        cutoff=args.cutoff,
        min_samples=args.min_samples,
    )


if __name__ == "__main__":
    main()
