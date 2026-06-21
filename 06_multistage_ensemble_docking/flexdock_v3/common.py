"""Shared utilities for the flexdock_v3 pipeline."""

from __future__ import annotations

import logging
import pickle
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np


LOGGER_NAME = "flexdock_v3"
REQUIRED_HIS_ATOMS = ("CA", "CB", "CG", "ND1", "CD2", "CE1", "NE2")
HIS_SIDECHAIN_ATOMS = ("CB", "CG", "ND1", "CD2", "CE1", "NE2")


@dataclass(frozen=True)
class PDBQTAtom:
    serial: int
    name: str
    resname: str
    chain: str
    resid: int
    coord: np.ndarray
    occupancy: float
    tempfactor: float
    charge: float | None
    atom_type: str
    element: str
    record: str

    def as_dict(self) -> dict:
        return {
            "serial": self.serial,
            "name": self.name,
            "resname": self.resname,
            "chain": self.chain,
            "resid": self.resid,
            "coord": np.asarray(self.coord, dtype=float),
            "occupancy": self.occupancy,
            "tempfactor": self.tempfactor,
            "charge": self.charge,
            "atom_type": self.atom_type,
            "element": self.element,
            "record": self.record,
        }


def configure_logging(log_path: Path | None = None, verbose: bool = False) -> logging.Logger:
    logger = logging.getLogger(LOGGER_NAME)
    logger.setLevel(logging.DEBUG)
    logger.propagate = False
    logger.handlers.clear()

    formatter = logging.Formatter(
        "%(asctime)s | %(levelname)-8s | %(name)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    console = logging.StreamHandler()
    console.setLevel(logging.DEBUG if verbose else logging.INFO)
    console.setFormatter(formatter)
    logger.addHandler(console)

    if log_path is not None:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_path, mode="w", encoding="utf-8")
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    return logger


def get_logger() -> logging.Logger:
    return logging.getLogger(LOGGER_NAME)


def infer_element(atom_name: str, atom_type: str = "") -> str:
    pdbqt_type = re.sub(r"[^A-Za-z]", "", atom_type).upper()
    if pdbqt_type:
        if pdbqt_type.startswith("CL"):
            return "Cl"
        if pdbqt_type.startswith("BR"):
            return "Br"
        if pdbqt_type == "A":
            return "C"
        if pdbqt_type.startswith("OA"):
            return "O"
        if pdbqt_type.startswith("NA"):
            return "N"
        if pdbqt_type.startswith("SA"):
            return "S"
        if pdbqt_type.startswith("HD"):
            return "H"
        if pdbqt_type[0] in "CNOSPHFIB":
            return pdbqt_type[0].title()

    letters = re.sub(r"[^A-Za-z]", "", atom_name)
    if not letters:
        return ""
    upper = letters.upper()
    if upper.startswith("CL"):
        return "Cl"
    if upper.startswith("BR"):
        return "Br"
    return upper[0].title()


def parse_pdbqt_atom(line: str) -> PDBQTAtom:
    if not line.startswith(("ATOM  ", "HETATM")):
        raise ValueError(f"Not an ATOM/HETATM line: {line!r}")

    tokens = line.split()
    try:
        serial = int(line[6:11])
        name = line[12:16].strip()
        resname = line[17:20].strip()
        chain = line[21:22].strip()
        resid = int(line[22:26])
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        occupancy = float(line[54:60] or 1.0)
        tempfactor = float(line[60:66] or 0.0)
    except (ValueError, IndexError) as exc:
        raise ValueError(f"Malformed PDBQT atom line: {line.rstrip()}") from exc

    atom_type = tokens[-1] if tokens else ""
    charge = None
    if len(tokens) >= 2:
        try:
            charge = float(tokens[-2])
        except ValueError:
            charge = None

    return PDBQTAtom(
        serial=serial,
        name=name,
        resname=resname,
        chain=chain,
        resid=resid,
        coord=np.array([x, y, z], dtype=float),
        occupancy=occupancy,
        tempfactor=tempfactor,
        charge=charge,
        atom_type=atom_type,
        element=infer_element(name, atom_type),
        record=line[:6].strip() or "ATOM",
    )


def atoms_by_name(atoms: Iterable[dict]) -> dict[str, dict]:
    result: dict[str, dict] = {}
    for atom in atoms:
        name = str(atom["name"]).strip().upper()
        if name in result:
            raise ValueError(f"Duplicate atom name {name}")
        result[name] = atom
    return result


def coordinates_by_name(atoms: Iterable[dict], names: Iterable[str]) -> np.ndarray:
    lookup = atoms_by_name(atoms)
    missing = [name for name in names if name not in lookup]
    if missing:
        raise ValueError(f"Missing required atoms: {missing}")
    return np.vstack([np.asarray(lookup[name]["coord"], dtype=float) for name in names])


def kabsch_transform(mobile: np.ndarray, reference: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return rotation and translation mapping mobile onto reference."""
    mobile = np.asarray(mobile, dtype=float)
    reference = np.asarray(reference, dtype=float)
    if mobile.shape != reference.shape or mobile.ndim != 2 or mobile.shape[1] != 3:
        raise ValueError("Kabsch inputs must have matching (N, 3) shapes")

    mobile_center = mobile.mean(axis=0)
    reference_center = reference.mean(axis=0)
    mobile_centered = mobile - mobile_center
    reference_centered = reference - reference_center
    covariance = mobile_centered.T @ reference_centered
    u, _, vt = np.linalg.svd(covariance)
    correction = np.eye(3)
    if np.linalg.det(u @ vt) < 0:
        correction[-1, -1] = -1
    rotation = u @ correction @ vt
    translation = reference_center - mobile_center @ rotation
    return rotation, translation


def kabsch_align(mobile: np.ndarray, reference: np.ndarray) -> np.ndarray:
    rotation, translation = kabsch_transform(mobile, reference)
    return np.asarray(mobile, dtype=float) @ rotation + translation


def kabsch_rmsd(coords_a: np.ndarray, coords_b: np.ndarray) -> float:
    aligned = kabsch_align(coords_a, coords_b)
    return float(np.sqrt(np.mean(np.sum((aligned - coords_b) ** 2, axis=1))))


def condensed_pairwise_rmsd(coordinates: list[np.ndarray], align: bool = True) -> np.ndarray:
    n_items = len(coordinates)
    matrix = np.zeros((n_items, n_items), dtype=float)
    for i in range(1, n_items):
        for j in range(i):
            if align:
                value = kabsch_rmsd(coordinates[i], coordinates[j])
            else:
                delta = coordinates[i] - coordinates[j]
                value = float(np.sqrt(np.mean(np.sum(delta * delta, axis=1))))
            matrix[i, j] = matrix[j, i] = value
    return matrix


def cluster_medoids(distance_matrix: np.ndarray, labels: np.ndarray) -> dict[int, int]:
    medoids: dict[int, int] = {}
    for label in sorted(set(int(value) for value in labels)):
        members = np.flatnonzero(labels == label)
        submatrix = distance_matrix[np.ix_(members, members)]
        local_index = int(np.argmin(submatrix.mean(axis=1)))
        medoids[label] = int(members[local_index])
    return medoids


def stable_relabel(labels: np.ndarray, scores: np.ndarray, modes: np.ndarray) -> np.ndarray:
    """Relabel clusters by best score, then lowest mode."""
    ranking = []
    for label in sorted(set(int(value) for value in labels)):
        members = np.flatnonzero(labels == label)
        ranking.append(
            (
                float(np.min(scores[members])),
                int(np.min(modes[members])),
                label,
            )
        )
    mapping = {
        old_label: new_label
        for new_label, (_, _, old_label) in enumerate(sorted(ranking), start=1)
    }
    return np.array([mapping[int(value)] for value in labels], dtype=int)


def save_pickle(data: object, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as handle:
        pickle.dump(data, handle, protocol=pickle.HIGHEST_PROTOCOL)


def load_pickle(path: Path):
    with path.open("rb") as handle:
        return pickle.load(handle)


def pose_lookup(pose_data: dict) -> dict[int, dict]:
    return {int(pose["mode"]): pose for pose in pose_data["poses"]}
