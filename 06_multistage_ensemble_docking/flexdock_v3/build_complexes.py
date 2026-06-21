"""Build representative protein-HIS508-ligand complexes using SDF ligand topology."""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Geometry import Point3D

from .common import (
    REQUIRED_HIS_ATOMS,
    configure_logging,
    get_logger,
    kabsch_rmsd,
    load_pickle,
    parse_pdbqt_atom,
    pose_lookup,
)


HIS_ATOM_ORDER = (
    "N",
    "H",
    "CA",
    "HA",
    "CB",
    "HB2",
    "HB3",
    "CG",
    "ND1",
    "HD1",
    "CD2",
    "HD2",
    "CE1",
    "HE1",
    "NE2",
    "HE2",
    "C",
    "O",
    "OXT",
)


def read_receptor_atoms(receptor_path: Path) -> tuple[list[dict], list[dict], int]:
    atoms: list[dict] = []
    his508: list[dict] = []
    insertion_index: int | None = None

    for line in receptor_path.read_text(encoding="utf-8").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        atom = parse_pdbqt_atom(line).as_dict()
        is_his508 = atom["resname"] == "HIS" and int(atom["resid"]) == 508
        if is_his508:
            if insertion_index is None:
                insertion_index = len(atoms)
            his508.append(atom)
        else:
            atoms.append(atom)

    if insertion_index is None or not his508:
        raise ValueError("Rigid receptor does not contain HIS508 backbone atoms")
    return atoms, his508, insertion_index


def merge_his508(rigid_atoms: list[dict], flexible_atoms: list[dict]) -> list[dict]:
    merged: dict[str, dict] = {}
    for atom in rigid_atoms:
        merged[str(atom["name"]).upper()] = atom
    for atom in flexible_atoms:
        merged[str(atom["name"]).upper()] = atom

    required_complete = {"N", "CA", "C", "O", *REQUIRED_HIS_ATOMS}
    missing = sorted(required_complete - set(merged))
    if missing:
        raise ValueError(f"Reconstructed HIS508 is missing atoms: {missing}")

    order = {name: index for index, name in enumerate(HIS_ATOM_ORDER)}
    return sorted(
        merged.values(),
        key=lambda atom: (order.get(str(atom["name"]).upper(), 999), atom["name"]),
    )


def load_ligand_template(sdf_path: Path) -> Chem.Mol:
    supplier = Chem.SDMolSupplier(str(sdf_path), removeHs=False, sanitize=True)
    molecules = [mol for mol in supplier if mol is not None]
    if len(molecules) != 1:
        raise ValueError(f"Expected exactly one valid SDF molecule, found {len(molecules)}")
    return molecules[0]


def substructure_matches(template: Chem.Mol, smiles_mol: Chem.Mol) -> list[tuple[int, ...]]:
    matches = list(
        template.GetSubstructMatches(
            smiles_mol,
            useChirality=True,
            uniquify=False,
            maxMatches=10000,
        )
    )
    if not matches:
        matches = list(
            template.GetSubstructMatches(
                smiles_mol,
                useChirality=False,
                uniquify=False,
                maxMatches=10000,
            )
        )
    return matches


def determine_atom_mapping(
    template: Chem.Mol,
    pose: dict,
) -> tuple[dict[int, int], str]:
    """Map PDBQT serial -> SDF atom index using SMILES IDX and graph matching."""
    smiles_mol = Chem.MolFromSmiles(pose["smiles"])
    if smiles_mol is None:
        raise ValueError("RDKit could not parse REMARK SMILES")
    if smiles_mol.GetNumAtoms() != template.GetNumAtoms():
        raise ValueError(
            "SMILES and SDF atom counts differ: "
            f"{smiles_mol.GetNumAtoms()} vs {template.GetNumAtoms()}"
        )

    matches = substructure_matches(template, smiles_mol)
    if not matches:
        raise ValueError("SDF topology does not match REMARK SMILES")

    pdbqt_atoms = {
        int(atom["serial"]): atom
        for atom in pose["ligand_atoms"]
    }
    idx_pairs = [(int(first), int(second)) for first, second in pose["smiles_idx_pairs"]]
    orientations = {
        "first_smiles_second_pdbqt": idx_pairs,
        "first_pdbqt_second_smiles": [
            (second, first) for first, second in idx_pairs
        ],
    }

    template_coordinates = template.GetConformer().GetPositions()
    best_key: tuple[int, float] | None = None
    best_mapping: dict[int, int] | None = None
    best_orientation: str | None = None
    for orientation, pairs in orientations.items():
        for match in matches:
            mapping: dict[int, int] = {}
            element_matches = 0
            valid = True
            for smiles_index, pdbqt_serial in pairs:
                if not 1 <= smiles_index <= len(match):
                    valid = False
                    break
                if pdbqt_serial not in pdbqt_atoms:
                    valid = False
                    break
                sdf_index = int(match[smiles_index - 1])
                if sdf_index in mapping.values() or pdbqt_serial in mapping:
                    valid = False
                    break
                mapping[pdbqt_serial] = sdf_index
                pdbqt_element = str(pdbqt_atoms[pdbqt_serial]["element"]).upper()
                sdf_element = template.GetAtomWithIdx(sdf_index).GetSymbol().upper()
                element_matches += int(pdbqt_element == sdf_element)

            if not valid or len(mapping) != template.GetNumAtoms():
                continue
            mapped_coordinates = np.zeros_like(template_coordinates)
            for pdbqt_serial, sdf_index in mapping.items():
                mapped_coordinates[sdf_index] = np.asarray(
                    pdbqt_atoms[pdbqt_serial]["coord"],
                    dtype=float,
                )
            conformer_rmsd = kabsch_rmsd(mapped_coordinates, template_coordinates)
            candidate_key = (element_matches, -conformer_rmsd)
            if best_key is None or candidate_key > best_key:
                best_key = candidate_key
                best_mapping = mapping
                best_orientation = orientation

    if best_key is None or best_mapping is None or best_orientation is None:
        raise ValueError("Could not construct a complete PDBQT-to-SDF atom mapping")
    if best_key[0] != template.GetNumAtoms():
        raise ValueError(
            f"Best atom mapping matched {best_key[0]}/{template.GetNumAtoms()} elements"
        )
    return best_mapping, best_orientation


def docked_ligand_molecule(
    template: Chem.Mol,
    pose: dict,
    mapping: dict[int, int],
) -> Chem.Mol:
    molecule = Chem.Mol(template)
    molecule.RemoveAllConformers()
    conformer = Chem.Conformer(molecule.GetNumAtoms())
    pdbqt_atoms = {
        int(atom["serial"]): atom
        for atom in pose["ligand_atoms"]
    }
    assigned = set()
    for pdbqt_serial, sdf_index in mapping.items():
        coordinate = np.asarray(pdbqt_atoms[pdbqt_serial]["coord"], dtype=float)
        conformer.SetAtomPosition(
            int(sdf_index),
            Point3D(float(coordinate[0]), float(coordinate[1]), float(coordinate[2])),
        )
        assigned.add(int(sdf_index))
    if assigned != set(range(molecule.GetNumAtoms())):
        raise ValueError("Not every SDF atom received a docked coordinate")
    molecule.AddConformer(conformer, assignId=True)

    element_counts: Counter[str] = Counter()
    for atom in molecule.GetAtoms():
        symbol = atom.GetSymbol()
        element_counts[symbol] += 1
        atom_name = f"{symbol}{element_counts[symbol]}"
        info = Chem.AtomPDBResidueInfo()
        info.SetName(f"{atom_name:>4}"[-4:])
        info.SetResidueName("UNL")
        info.SetResidueNumber(1)
        info.SetChainId("L")
        info.SetIsHeteroAtom(True)
        info.SetOccupancy(1.0)
        info.SetTempFactor(0.0)
        atom.SetMonomerInfo(info)
    return molecule


def format_pdb_atom(
    serial: int,
    atom: dict,
    record: str | None = None,
    resname: str | None = None,
    chain: str | None = None,
    resid: int | None = None,
) -> str:
    atom_name = str(atom["name"]).strip()[:4]
    atom_resname = (resname or str(atom["resname"]).strip() or "UNK")[:3]
    atom_chain = (chain if chain is not None else str(atom.get("chain", "")))[:1]
    atom_resid = int(resid if resid is not None else atom["resid"])
    x, y, z = np.asarray(atom["coord"], dtype=float)
    occupancy = float(atom.get("occupancy", 1.0))
    tempfactor = float(atom.get("tempfactor", 0.0))
    element = str(atom.get("element", ""))[:2].rjust(2)
    record_name = (record or str(atom.get("record", "ATOM"))).upper()
    if record_name not in {"ATOM", "HETATM"}:
        record_name = "ATOM"
    return (
        f"{record_name:<6}{serial:5d} {atom_name:>4} "
        f"{atom_resname:>3} {atom_chain:1}{atom_resid:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}"
        f"{occupancy:6.2f}{tempfactor:6.2f}          {element}"
    )


def ligand_pdb_records(
    molecule: Chem.Mol,
    first_serial: int,
) -> tuple[list[str], list[str]]:
    conformer = molecule.GetConformer()
    atom_lines = []
    serial_by_index = {}
    for atom in molecule.GetAtoms():
        index = atom.GetIdx()
        serial = first_serial + index
        serial_by_index[index] = serial
        info = atom.GetPDBResidueInfo()
        position = conformer.GetAtomPosition(index)
        atom_dict = {
            "name": info.GetName().strip() if info else f"{atom.GetSymbol()}{index + 1}",
            "resname": "UNL",
            "chain": "L",
            "resid": 1,
            "coord": np.array([position.x, position.y, position.z]),
            "occupancy": 1.0,
            "tempfactor": 0.0,
            "element": atom.GetSymbol(),
            "record": "HETATM",
        }
        atom_lines.append(format_pdb_atom(serial, atom_dict, record="HETATM"))

    conect_lines = []
    for atom in molecule.GetAtoms():
        source = serial_by_index[atom.GetIdx()]
        neighbors: list[int] = []
        for bond in atom.GetBonds():
            other = bond.GetOtherAtomIdx(atom.GetIdx())
            order = int(round(bond.GetBondTypeAsDouble()))
            if bond.GetIsAromatic():
                order = 1
            neighbors.extend([serial_by_index[other]] * max(1, order))
        if neighbors:
            for start in range(0, len(neighbors), 4):
                chunk = neighbors[start : start + 4]
                conect_lines.append(
                    f"CONECT{source:5d}" + "".join(f"{neighbor:5d}" for neighbor in chunk)
                )
    return atom_lines, conect_lines


def validate_complex_lines(lines: list[str], expected_ligand_atoms: int) -> None:
    ligand_lines = [
        line
        for line in lines
        if line.startswith("HETATM") and line[17:20].strip() == "UNL"
    ]
    if len(ligand_lines) != expected_ligand_atoms:
        raise ValueError(
            f"Complex has {len(ligand_lines)} ligand atoms; expected {expected_ligand_atoms}"
        )

    his_lines = [
        line
        for line in lines
        if line.startswith("ATOM")
        and line[17:20].strip() == "HIS"
        and int(line[22:26]) == 508
    ]
    his_names = {line[12:16].strip() for line in his_lines}
    required = {"N", "CA", "C", "O", *REQUIRED_HIS_ATOMS}
    missing = sorted(required - his_names)
    if missing:
        raise ValueError(f"Output complex HIS508 is incomplete: {missing}")


def write_complex(
    receptor_atoms: list[dict],
    rigid_his: list[dict],
    his_insert_index: int,
    flexible_his: list[dict],
    ligand: Chem.Mol,
    output_path: Path,
    mode: int,
    vina_score: float,
    mapping_orientation: str,
) -> None:
    complete_his = merge_his508(rigid_his, flexible_his)
    protein_atoms = (
        receptor_atoms[:his_insert_index]
        + complete_his
        + receptor_atoms[his_insert_index:]
    )

    lines = [
        "REMARK 900 FLEXDOCK_V3 REPRESENTATIVE COMPLEX",
        f"REMARK 901 VINA MODE {mode} SCORE {vina_score:.3f}",
        f"REMARK 902 LIGAND COORDINATE MAP {mapping_orientation}",
        "REMARK 903 LIGAND TOPOLOGY FROM SDF; CONECT RECORDS INCLUDE BOND ORDER",
    ]
    serial = 1
    for atom in protein_atoms:
        lines.append(format_pdb_atom(serial, atom, record="ATOM"))
        serial += 1
    lines.append("TER")

    ligand_lines, conect_lines = ligand_pdb_records(ligand, serial)
    lines.extend(ligand_lines)
    lines.extend(conect_lines)
    lines.append("END")
    validate_complex_lines(lines, ligand.GetNumAtoms())
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_complexes(
    receptor_path: Path,
    ligand_sdf_path: Path,
    pose_data_path: Path,
    representatives_path: Path,
    results_dir: Path,
) -> pd.DataFrame:
    logger = get_logger()
    pose_data = load_pickle(pose_data_path)
    poses = pose_lookup(pose_data)
    representatives = pd.read_csv(representatives_path)
    template = load_ligand_template(ligand_sdf_path)

    expected_atoms = int(pose_data["expected_ligand_atoms"])
    if template.GetNumAtoms() != expected_atoms:
        raise ValueError(
            f"SDF contains {template.GetNumAtoms()} atoms; expected {expected_atoms}"
        )

    receptor_atoms, rigid_his, insertion_index = read_receptor_atoms(receptor_path)
    output_dir = results_dir / "representative_complexes"
    output_dir.mkdir(parents=True, exist_ok=True)
    for pattern in ("cluster_*.pdb", "cluster_*_ligand.sdf"):
        for old_output in output_dir.glob(pattern):
            old_output.unlink()
    old_manifest = output_dir / "complex_manifest.csv"
    if old_manifest.exists():
        old_manifest.unlink()

    manifest_rows = []
    mapping_orientation = None
    mapping = None
    for row in representatives.itertuples(index=False):
        mode = int(row.mode)
        if mode not in poses:
            raise ValueError(f"Representative mode {mode} not found in pose_data")
        pose = poses[mode]

        if mapping is None:
            mapping, mapping_orientation = determine_atom_mapping(template, pose)
            logger.info(
                "Validated PDBQT-to-SDF mapping: %s (%d atoms)",
                mapping_orientation,
                len(mapping),
            )

        ligand = docked_ligand_molecule(template, pose, mapping)
        his_cluster = int(row.his508_cluster)
        ligand_cluster = int(row.ligand_cluster)
        stem = f"cluster_{his_cluster}_{ligand_cluster}"
        pdb_path = output_dir / f"{stem}.pdb"
        sdf_path = output_dir / f"{stem}_ligand.sdf"

        write_complex(
            receptor_atoms,
            rigid_his,
            insertion_index,
            pose["his508_atoms"],
            ligand,
            pdb_path,
            mode,
            float(row.vina_score),
            str(mapping_orientation),
        )
        writer = Chem.SDWriter(str(sdf_path))
        ligand.SetProp("VinaMode", str(mode))
        ligand.SetProp("VinaScore", f"{float(row.vina_score):.3f}")
        ligand.SetProp("CombinedCluster", str(row.combined_cluster))
        writer.write(ligand)
        writer.close()

        manifest_rows.append(
            {
                "combined_cluster": row.combined_cluster,
                "his508_cluster": his_cluster,
                "ligand_cluster": ligand_cluster,
                "mode": mode,
                "vina_score": float(row.vina_score),
                "complex_pdb": str(pdb_path.resolve()),
                "ligand_sdf": str(sdf_path.resolve()),
                "ligand_atom_count": ligand.GetNumAtoms(),
                "his508_atom_count": len(merge_his508(rigid_his, pose["his508_atoms"])),
                "mapping_orientation": mapping_orientation,
            }
        )
        logger.debug("Built %s from MODEL %d", pdb_path.name, mode)

    manifest = pd.DataFrame(manifest_rows)
    manifest_path = output_dir / "complex_manifest.csv"
    manifest.to_csv(manifest_path, index=False)
    logger.info("Built %d representative complexes in %s", len(manifest), output_dir)
    logger.info("Wrote %s", manifest_path)
    return manifest


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("receptor_rigid", type=Path)
    parser.add_argument("ligand_sdf", type=Path)
    parser.add_argument("--pose-data", type=Path, default=Path("results/pose_data.pkl"))
    parser.add_argument(
        "--representatives",
        type=Path,
        default=Path("results/representative_modes.csv"),
    )
    parser.add_argument("--results-dir", type=Path, default=Path("results"))
    parser.add_argument("--verbose", action="store_true")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    configure_logging(args.results_dir / "build_complexes.log", args.verbose)
    build_complexes(
        args.receptor_rigid.resolve(),
        args.ligand_sdf.resolve(),
        args.pose_data.resolve(),
        args.representatives.resolve(),
        args.results_dir.resolve(),
    )


if __name__ == "__main__":
    main()
