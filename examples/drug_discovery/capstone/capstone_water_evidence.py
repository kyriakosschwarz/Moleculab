# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Prepare deposited-water evidence for the selected ligand and protein states.

This is Batch 12B of the capstone.  It preserves the identity, coordinates,
occupancy, and B-factor of each crystallographic water near the ligand.  Water
orientation and donor/acceptor roles remain unassigned, and the record reports
only whether the geometry is compatible with a water bridge.

A geometrically plausible bridge is not a chemically classified water-mediated
interaction, and no affinity is computed.

Run with gemmi available:

    uv run --with gemmi python examples/drug_discovery/capstone/capstone_water_evidence.py
"""

from __future__ import annotations

import argparse
from pathlib import Path
from types import SimpleNamespace

import gemmi

from moleculab import (
    WaterAtomRecord,
    assess_water_bridges,
    build_water_evidence_record,
    load_binding_site,
    load_chemical_preparation,
    load_ligand_coordinates,
    load_ligand_identity,
    load_manifest,
    load_protein_preparation,
    sha256_file,
    validate_water_evidence,
    write_water_evidence,
)

DEFAULT_MANIFEST = Path("data/structures/2HYY.manifest.json")
DEFAULT_IDENTITY = Path("data/structures/2HYY.ligand-identity.json")
DEFAULT_LIGAND = Path("data/structures/2HYY.ligand-coordinates.json")
DEFAULT_SITE = Path("data/structures/2HYY.binding-site.json")
DEFAULT_PREPARATION = Path("data/structures/2HYY.chemical-preparation.json")
DEFAULT_PROTEIN = Path("data/structures/2HYY.protein-preparation.json")
DEFAULT_SOURCE = Path("data/structures/2HYY.cif")
DEFAULT_OUTPUT = Path("data/structures/2HYY.water-evidence.json")
BRIDGE_CUTOFF_ANGSTROM = 3.5


def _argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--identity", type=Path, default=DEFAULT_IDENTITY)
    parser.add_argument("--ligand", type=Path, default=DEFAULT_LIGAND)
    parser.add_argument("--site", type=Path, default=DEFAULT_SITE)
    parser.add_argument("--preparation", type=Path, default=DEFAULT_PREPARATION)
    parser.add_argument("--protein", type=Path, default=DEFAULT_PROTEIN)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser


def _unquote(value: str) -> str:
    value = value.strip()
    if len(value) >= 2 and value[0] == value[-1] and value[0] in ("'", '"'):
        return value[1:-1]
    return value


def _none_if_missing(value: str) -> str | None:
    return None if value in ("?", ".") else value


def _atom_rows(block: gemmi.cif.Block):
    return [
        [_unquote(value) for value in row]
        for row in block.find(
            "_atom_site.",
            [
                "group_PDB",
                "label_comp_id",
                "label_asym_id",
                "auth_asym_id",
                "label_seq_id",
                "auth_seq_id",
                "pdbx_PDB_ins_code",
                "label_atom_id",
                "type_symbol",
                "Cartn_x",
                "Cartn_y",
                "Cartn_z",
                "occupancy",
                "B_iso_or_equiv",
                "pdbx_PDB_model_num",
            ],
        )
    ]


def _deposited_waters(
    block: gemmi.cif.Block,
    auth_chain: str,
    model_number: int,
) -> dict[tuple[str, ...], WaterAtomRecord]:
    waters = {}
    for values in _atom_rows(block):
        if values[0] != "HETATM" or values[1] != "HOH":
            continue
        if values[3] != auth_chain:
            continue
        if values[14] not in ("?", ".") and int(values[14]) != model_number:
            continue
        key = (
            values[2],
            values[3],
            values[5],
            values[7],
        )
        waters[key] = WaterAtomRecord(
            water_id=f"{values[2]}:{values[5]}:{values[7]}",
            label_asym_id=values[2],
            auth_asym_id=values[3],
            label_seq_id=_none_if_missing(values[4]),
            auth_seq_id=values[5],
            atom_name=values[7],
            position_angstrom=(
                float(values[9]),
                float(values[10]),
                float(values[11]),
            ),
            occupancy=None if values[12] in ("?", ".") else float(values[12]),
            b_iso_or_equiv=None if values[13] in ("?", ".") else float(values[13]),
            orientation_known=False,
            donor_role=None,
            acceptor_role=None,
            source="deposited HETATM HOH record",
        )
    return waters


def main() -> None:
    args = _argument_parser().parse_args()
    for path in (
        args.manifest,
        args.identity,
        args.ligand,
        args.site,
        args.preparation,
        args.protein,
        args.source,
    ):
        if not path.exists():
            raise SystemExit(f"missing {path}; run the earlier capstone batches first.")

    manifest = load_manifest(args.manifest)
    identity = load_ligand_identity(args.identity)
    ligand = load_ligand_coordinates(args.ligand)
    site = load_binding_site(args.site)
    preparation = load_chemical_preparation(args.preparation)
    protein = load_protein_preparation(args.protein)
    if not (
        identity.accession
        == ligand.accession
        == site.accession
        == preparation.accession
        == protein.accession
        == manifest.source.accession
    ):
        raise SystemExit("capstone records disagree on accession.")

    ligand_state = next(
        state for state in preparation.states if state.state_id == "reference-neutral"
    )
    protein_state = next(
        state for state in protein.states if state.state_id == "standard-templates"
    )

    block = gemmi.cif.read(str(args.source)).sole_block()
    deposited = _deposited_waters(
        block,
        manifest.selection.protein.auth_asym_id,
        manifest.selection.model_number,
    )
    waters = []
    for nearby in site.waters_near_ligand:
        key = (
            nearby.label_asym_id,
            nearby.auth_asym_id,
            nearby.auth_seq_id,
            nearby.atom_name,
        )
        water = deposited.get(key)
        if water is None:
            raise SystemExit(f"missing deposited water record for {key}.")
        waters.append(water)

    ligand_positions = {
        atom.atom_name: atom.position_angstrom for atom in ligand.atoms
    }
    ligand_polar = [
        SimpleNamespace(
            atom_id=atom.atom_id,
            position_angstrom=ligand_positions[atom.atom_id],
        )
        for atom in ligand_state.atoms
        if (atom.is_donor or atom.is_acceptor) and atom.atom_id in ligand_positions
    ]

    contact_atoms = {
        (contact.auth_asym_id, contact.auth_seq_id, contact.protein_atom)
        for contact in site.contacts
    }
    protein_polar = [
        atom
        for atom in protein_state.atoms
        if (atom.is_donor or atom.is_acceptor)
        and (atom.auth_asym_id, atom.auth_seq_id, atom.atom_name) in contact_atoms
    ]

    bridges = assess_water_bridges(
        waters,
        ligand_polar,
        protein_polar,
        BRIDGE_CUTOFF_ANGSTROM,
    )
    record = build_water_evidence_record(
        source_manifest_path=str(args.manifest),
        source_manifest_sha256=sha256_file(args.manifest),
        binding_site_record_path=str(args.site),
        binding_site_record_sha256=sha256_file(args.site),
        ligand_identity_record_path=str(args.identity),
        ligand_identity_record_sha256=sha256_file(args.identity),
        chemical_preparation_record_path=str(args.preparation),
        chemical_preparation_record_sha256=sha256_file(args.preparation),
        protein_preparation_record_path=str(args.protein),
        protein_preparation_record_sha256=sha256_file(args.protein),
        accession=manifest.source.accession,
        protein_selection=manifest.selection.protein,
        ligand_selection=identity.selection,
        ligand_state_id=ligand_state.state_id,
        protein_state_id=protein_state.state_id,
        bridge_cutoff_angstrom=BRIDGE_CUTOFF_ANGSTROM,
        waters=tuple(waters),
        bridges=bridges,
    )
    errors = validate_water_evidence(record)
    if errors:
        for error in errors:
            print(f"water evidence error: {error}")
        raise SystemExit(2)

    destination = write_water_evidence(args.output, record)
    plausible = sum(1 for item in bridges if item.status == "geometrically_plausible")
    not_supported = sum(1 for item in bridges if item.status == "not_supported")
    unknown = sum(1 for item in bridges if item.status == "unknown")
    print(f"{record.accession} water evidence")
    print(
        f"  ligand state: {record.ligand_state_id}; "
        f"protein state: {record.protein_state_id}"
    )
    print(
        f"  deposited waters: {len(record.waters)}; "
        f"bridge cutoff: {record.bridge_cutoff_angstrom:.1f} A"
    )
    print(
        f"  bridge geometry: plausible={plausible}, "
        f"not_supported={not_supported}, unknown={unknown}"
    )
    print("  water roles: unassigned; orientation: unknown")
    print(f"  record: {destination}")
    print(
        "  scope: geometric bridge plausibility only; "
        "no water-mediated classification or affinity"
    )


if __name__ == "__main__":
    main()
