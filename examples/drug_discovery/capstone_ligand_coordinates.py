"""Reconcile selected ligand coordinates with the Batch 2 identity record.

This is Batch 3 of the capstone.  It records the deposited coordinates,
occupancy, displacement parameters, alternate locations, bond lengths, angles,
aromatic-ring planarity, and close contacts.  It does not add missing atoms or
hydrogens, choose between alternate conformers, assign protonation states, type
force-field atoms, minimise, dock, or estimate affinity.

Run with gemmi available:

    uv run --with gemmi python examples/drug_discovery/capstone_ligand_coordinates.py
"""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

import gemmi

from moleculab import (
    CoordinateAtom,
    load_ligand_identity,
    load_manifest,
    reconcile_ligand_coordinates,
    sha256_file,
    validate_ligand_coordinates,
    write_ligand_coordinates,
)

DEFAULT_MANIFEST = Path("data/structures/2HYY.manifest.json")
DEFAULT_IDENTITY = Path("data/structures/2HYY.ligand-identity.json")
DEFAULT_SOURCE = Path("data/structures/2HYY.cif")
DEFAULT_OUTPUT = Path("data/structures/2HYY.ligand-coordinates.json")


def _argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--identity", type=Path, default=DEFAULT_IDENTITY)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser


def _unquote(value: str) -> str:
    value = value.strip()
    if len(value) >= 2 and value[0] == value[-1] and value[0] in ("'", '"'):
        return value[1:-1]
    return value


def _component_atom_elements(
    block: gemmi.cif.Block,
    component_id: str,
) -> dict[str, str]:
    elements = {}
    for row in block.find(
        "_chem_comp_atom.",
        ["comp_id", "atom_id", "type_symbol"],
    ):
        values = [_unquote(value) for value in row]
        if values[0] == component_id:
            elements[values[1]] = values[2]
    if not elements:
        raise ValueError(
            f"_chem_comp_atom has no records for component {component_id!r}."
        )
    return elements


def _coordinate_atoms(
    block: gemmi.cif.Block,
    selection,
    model_number: int,
    reference_indices: dict[str, int],
) -> tuple[CoordinateAtom, ...]:
    atoms = []
    for row in block.find(
        "_atom_site.",
        [
            "label_comp_id",
            "label_asym_id",
            "auth_asym_id",
            "auth_seq_id",
            "label_seq_id",
            "pdbx_PDB_ins_code",
            "label_atom_id",
            "type_symbol",
            "Cartn_x",
            "Cartn_y",
            "Cartn_z",
            "occupancy",
            "B_iso_or_equiv",
            "label_alt_id",
            "pdbx_PDB_model_num",
        ],
    ):
        values = [_unquote(value) for value in row]
        if values[0] != selection.component_id:
            continue
        if selection.label_asym_id is not None and values[1] != selection.label_asym_id:
            continue
        if values[2] != selection.auth_asym_id:
            continue
        if values[3] != selection.auth_seq_id:
            continue
        if values[14] not in ("?", ".") and int(values[14]) != model_number:
            continue
        atoms.append(
            CoordinateAtom(
                atom_name=values[6],
                element=values[7],
                reference_index=reference_indices.get(values[6]),
                position_angstrom=(
                    float(values[8]),
                    float(values[9]),
                    float(values[10]),
                ),
                occupancy=(
                    None if values[11] in ("?", ".") else float(values[11])
                ),
                b_iso_or_equiv=(
                    None if values[12] in ("?", ".") else float(values[12])
                ),
                alternate_location=(
                    None if values[13] in ("?", ".") else values[13]
                ),
                model_number=int(values[14]) if values[14] not in ("?", ".") else 1,
            )
        )
    if not atoms:
        raise ValueError(
            "the Batch 1 selection did not match any deposited coordinate records."
        )
    return tuple(atoms)


def main() -> None:
    args = _argument_parser().parse_args()
    for path in (args.manifest, args.identity, args.source):
        if not path.exists():
            raise SystemExit(
                f"missing {path}; run the earlier capstone batches first."
            )

    manifest = load_manifest(args.manifest)
    identity = load_ligand_identity(args.identity)
    block = gemmi.cif.read_file(str(args.source)).sole_block()
    reference_indices = dict(identity.mapping.pairs)
    atoms = _coordinate_atoms(
        block,
        identity.selection,
        manifest.selection.model_number,
        reference_indices,
    )
    component_atom_elements = _component_atom_elements(
        block, identity.selection.component_id
    )

    record = reconcile_ligand_coordinates(
        source_manifest_path=str(args.manifest),
        source_manifest_sha256=sha256_file(args.manifest),
        identity_record_path=str(args.identity),
        identity_record_sha256=sha256_file(args.identity),
        accession=manifest.source.accession,
        selection=identity.selection,
        atoms=atoms,
        bonds=identity.deposited.bonds,
        expected_atom_names=[name for name, _ in identity.mapping.pairs],
        component_atom_elements=component_atom_elements,
        reference_stereocenters=identity.reference.stereocenters,
        deposited_stereochemistry=identity.deposited.stereochemistry,
    )
    errors = validate_ligand_coordinates(record)
    if errors:
        for error in errors:
            print(f"coordinate record error: {error}")
        raise SystemExit(2)

    destination = write_ligand_coordinates(args.output, record)
    statuses = Counter(check.status for check in record.validation.checks)
    discrepancies = Counter(item.status for item in record.discrepancies)
    failures = [check for check in record.validation.checks if check.status == "fail"]

    print(f"{record.accession} ligand coordinate reconciliation")
    print(
        f"  selected: {record.selection.component_id} label "
        f"{record.selection.label_asym_id} auth {record.selection.auth_asym_id} "
        f"seq {record.selection.auth_seq_id}"
    )
    print(
        f"  coordinates: {len(record.atoms)} atoms, {len(record.bonds)} bonds, "
        f"{len(record.angles)} angles, {len(record.rings)} aromatic ring systems"
    )
    print(
        f"  missing: {len(record.missing_component_atoms)} component atoms "
        f"({len(record.missing_heavy_atoms)} heavy), "
        f"{len(record.alternate_location_atoms)} alternate-location atoms"
    )
    print(
        f"  geometry: {len(record.close_contacts)} close contacts; "
        f"criteria {record.criteria.bond_length_min_angstrom}-"
        f"{record.criteria.bond_length_max_angstrom} A bonds"
    )
    print(
        f"  validation: {statuses['pass']} pass, {statuses['fail']} fail, "
        f"{statuses['not_run']} not run, {statuses['not_applicable']} not applicable"
    )
    print(
        f"  discrepancies: {discrepancies['resolved']} resolved, "
        f"{discrepancies['open']} open, "
        f"{discrepancies['not_applicable']} not applicable"
    )
    print(f"  record: {destination}")
    print("  scope: no atoms or hydrogens were added; no conformer was chosen")

    if failures:
        for check in failures:
            print(f"  FAILED {check.check_id}: {check.detail}")
        raise SystemExit(1)


if __name__ == "__main__":
    main()
