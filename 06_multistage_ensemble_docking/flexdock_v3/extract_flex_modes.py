"""Extract ligand and flexible HIS508 coordinates from Vina PDBQT models."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from .common import (
    REQUIRED_HIS_ATOMS,
    configure_logging,
    get_logger,
    parse_pdbqt_atom,
    save_pickle,
)


VINA_RESULT_RE = re.compile(
    r"^REMARK VINA RESULT:\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)"
)
MODEL_RE = re.compile(r"^MODEL\s+(\d+)")
SMILES_IDX_RE = re.compile(r"^REMARK SMILES IDX\s+(.+)$")
SMILES_RE = re.compile(r"^REMARK SMILES\s+(\S+)")


def parse_idx_pairs(text: str) -> list[tuple[int, int]]:
    values = [int(value) for value in text.split()]
    if len(values) % 2:
        raise ValueError(f"Odd number of SMILES IDX values: {text!r}")
    return list(zip(values[0::2], values[1::2]))


def validate_pose(pose: dict, expected_ligand_atoms: int) -> None:
    logger = get_logger()
    mode = pose["mode"]
    ligand_atoms = pose["ligand_atoms"]
    his_atoms = pose["his508_atoms"]

    if len(ligand_atoms) != expected_ligand_atoms:
        raise ValueError(
            f"MODEL {mode}: expected {expected_ligand_atoms} ligand atoms, "
            f"found {len(ligand_atoms)}"
        )
    ligand_serials = [int(atom["serial"]) for atom in ligand_atoms]
    if len(set(ligand_serials)) != expected_ligand_atoms:
        raise ValueError(f"MODEL {mode}: ligand atom serials are not unique")

    his_names = {str(atom["name"]).upper() for atom in his_atoms}
    missing = [name for name in REQUIRED_HIS_ATOMS if name not in his_names]
    if missing:
        raise ValueError(f"MODEL {mode}: missing HIS508 atoms {missing}")

    idx_pairs = pose["smiles_idx_pairs"]
    first = {first for first, _ in idx_pairs}
    second = {second for _, second in idx_pairs}
    expected = set(range(1, expected_ligand_atoms + 1))
    if first != expected or second != expected:
        raise ValueError(
            f"MODEL {mode}: SMILES IDX map is not a complete 1..{expected_ligand_atoms} permutation"
        )

    logger.debug(
        "Validated MODEL %d: ligand=%d atoms, HIS508=%d atoms",
        mode,
        len(ligand_atoms),
        len(his_atoms),
    )


def parse_flex_output(
    flex_out: Path,
    expected_ligand_atoms: int = 47,
) -> dict:
    logger = get_logger()
    lines = flex_out.read_text(encoding="utf-8").splitlines()
    poses: list[dict] = []
    current: dict | None = None
    in_flexible_residue = False

    for line_number, line in enumerate(lines, start=1):
        model_match = MODEL_RE.match(line)
        if model_match:
            if current is not None:
                raise ValueError(f"Line {line_number}: MODEL started before prior ENDMDL")
            current = {
                "mode": int(model_match.group(1)),
                "vina_score": None,
                "rmsd_lb": None,
                "rmsd_ub": None,
                "smiles": None,
                "smiles_idx_pairs": [],
                "ligand_atoms": [],
                "his508_atoms": [],
            }
            in_flexible_residue = False
            continue

        if current is None:
            continue

        result_match = VINA_RESULT_RE.match(line)
        if result_match:
            current["vina_score"] = float(result_match.group(1))
            current["rmsd_lb"] = float(result_match.group(2))
            current["rmsd_ub"] = float(result_match.group(3))
            continue

        smiles_match = SMILES_RE.match(line)
        if smiles_match and not line.startswith("REMARK SMILES IDX"):
            current["smiles"] = smiles_match.group(1)
            continue

        idx_match = SMILES_IDX_RE.match(line)
        if idx_match:
            current["smiles_idx_pairs"].extend(parse_idx_pairs(idx_match.group(1)))
            continue

        if line.startswith("BEGIN_RES"):
            fields = line.split()
            if len(fields) < 4:
                raise ValueError(f"Line {line_number}: malformed BEGIN_RES")
            resname, chain, resid = fields[1], fields[2], int(fields[3])
            if resname != "HIS" or resid != 508:
                raise ValueError(
                    f"Line {line_number}: expected flexible HIS508, got {resname} {chain} {resid}"
                )
            current["his508_chain"] = chain
            in_flexible_residue = True
            continue

        if line.startswith("END_RES"):
            in_flexible_residue = False
            continue

        if line.startswith(("ATOM  ", "HETATM")):
            atom = parse_pdbqt_atom(line).as_dict()
            if in_flexible_residue:
                current["his508_atoms"].append(atom)
            elif atom["resname"] == "UNL":
                current["ligand_atoms"].append(atom)
            continue

        if line.startswith("ENDMDL"):
            for key in ("vina_score", "rmsd_lb", "rmsd_ub", "smiles"):
                if current[key] is None:
                    raise ValueError(
                        f"MODEL {current['mode']}: missing required field {key}"
                    )
            validate_pose(current, expected_ligand_atoms)
            poses.append(current)
            current = None
            in_flexible_residue = False

    if current is not None:
        raise ValueError("File ended before ENDMDL")
    if not poses:
        raise ValueError(f"No MODEL records found in {flex_out}")

    modes = [pose["mode"] for pose in poses]
    if len(set(modes)) != len(modes):
        raise ValueError("Duplicate MODEL numbers found")

    reference_smiles = poses[0]["smiles"]
    reference_idx = poses[0]["smiles_idx_pairs"]
    for pose in poses[1:]:
        if pose["smiles"] != reference_smiles:
            raise ValueError(f"MODEL {pose['mode']}: SMILES differs from MODEL 1")
        if pose["smiles_idx_pairs"] != reference_idx:
            raise ValueError(f"MODEL {pose['mode']}: SMILES IDX map differs from MODEL 1")

    logger.info(
        "Parsed %d models from %s; every ligand has %d atoms",
        len(poses),
        flex_out,
        expected_ligand_atoms,
    )
    return {
        "source": str(flex_out.resolve()),
        "expected_ligand_atoms": expected_ligand_atoms,
        "required_his_atoms": list(REQUIRED_HIS_ATOMS),
        "smiles": reference_smiles,
        "smiles_idx_pairs": reference_idx,
        "poses": poses,
    }


def make_score_plot(summary: pd.DataFrame, output_path: Path) -> None:
    sns.set_theme(style="whitegrid", context="paper")
    figure, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    sns.histplot(
        data=summary,
        x="vina_score",
        bins=min(15, max(5, len(summary) // 3)),
        kde=True,
        color="#2A6F97",
        ax=axes[0],
    )
    axes[0].set_xlabel("Vina score (kcal/mol)")
    axes[0].set_ylabel("Pose count")
    axes[0].set_title("Docking score distribution")

    ordered = summary.sort_values("vina_score")
    sns.scatterplot(
        data=ordered,
        x="mode",
        y="vina_score",
        color="#C44536",
        s=38,
        edgecolor="white",
        linewidth=0.4,
        ax=axes[1],
    )
    axes[1].set_xlabel("Vina mode")
    axes[1].set_ylabel("Vina score (kcal/mol)")
    axes[1].set_title("Scores by docking mode")
    figure.tight_layout()
    figure.savefig(output_path, dpi=350, bbox_inches="tight")
    plt.close(figure)


def extract_flex_modes(
    flex_out: Path,
    results_dir: Path,
    expected_ligand_atoms: int = 47,
) -> dict:
    logger = get_logger()
    results_dir.mkdir(parents=True, exist_ok=True)
    pose_data = parse_flex_output(flex_out, expected_ligand_atoms)

    pose_path = results_dir / "pose_data.pkl"
    save_pickle(pose_data, pose_path)

    rows = []
    for pose in pose_data["poses"]:
        rows.append(
            {
                "mode": pose["mode"],
                "vina_score": pose["vina_score"],
                "rmsd_lb": pose["rmsd_lb"],
                "rmsd_ub": pose["rmsd_ub"],
                "ligand_atom_count": len(pose["ligand_atoms"]),
                "his508_atom_count": len(pose["his508_atoms"]),
                "his508_atom_names": ";".join(
                    atom["name"] for atom in pose["his508_atoms"]
                ),
            }
        )
    summary = pd.DataFrame(rows).sort_values("mode")
    summary_path = results_dir / "pose_summary.csv"
    summary.to_csv(summary_path, index=False)
    make_score_plot(summary, results_dir / "score_distribution.png")

    logger.info("Wrote %s", pose_path)
    logger.info("Wrote %s", summary_path)
    logger.info("Wrote %s", results_dir / "score_distribution.png")
    return pose_data


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("flex_out", type=Path)
    parser.add_argument("--results-dir", type=Path, default=Path("results"))
    parser.add_argument("--expected-ligand-atoms", type=int, default=47)
    parser.add_argument("--verbose", action="store_true")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    configure_logging(
        args.results_dir / "extract_flex_modes.log",
        verbose=args.verbose,
    )
    extract_flex_modes(
        args.flex_out.resolve(),
        args.results_dir.resolve(),
        args.expected_ligand_atoms,
    )


if __name__ == "__main__":
    main()
