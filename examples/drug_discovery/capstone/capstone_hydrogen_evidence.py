# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Prepare hydrogen evidence for the selected ligand and protein states.

This is Batch 12A of the capstone.  It keeps deposited heavy atoms separate from
generated ligand hydrogens, validates donor-hydrogen connectivity, and evaluates
ligand-donor hydrogen bonds against explicit protein acceptor features.

Protein hydrogens are not generated in this batch, so protein-donor hydrogen
bonds remain unknown.

Run:

    uv run python examples/drug_discovery/capstone/capstone_hydrogen_evidence.py
"""

from __future__ import annotations

import argparse
from pathlib import Path

from moleculab import (
    HydrogenBondCriteria,
    HydrogenPosition,
    assess_ligand_donor_hydrogen_bonds,
    build_hydrogen_evidence_record,
    load_binding_site,
    load_chemical_preparation,
    load_ligand_coordinates,
    load_ligand_identity,
    load_manifest,
    load_protein_preparation,
    sha256_file,
    validate_hydrogen_evidence,
    write_hydrogen_evidence,
)

DEFAULT_MANIFEST = Path("data/structures/2HYY.manifest.json")
DEFAULT_IDENTITY = Path("data/structures/2HYY.ligand-identity.json")
DEFAULT_LIGAND = Path("data/structures/2HYY.ligand-coordinates.json")
DEFAULT_SITE = Path("data/structures/2HYY.binding-site.json")
DEFAULT_PREPARATION = Path("data/structures/2HYY.chemical-preparation.json")
DEFAULT_PROTEIN = Path("data/structures/2HYY.protein-preparation.json")
DEFAULT_OUTPUT = Path("data/structures/2HYY.hydrogen-evidence.json")


def _argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--identity", type=Path, default=DEFAULT_IDENTITY)
    parser.add_argument("--ligand", type=Path, default=DEFAULT_LIGAND)
    parser.add_argument("--site", type=Path, default=DEFAULT_SITE)
    parser.add_argument("--preparation", type=Path, default=DEFAULT_PREPARATION)
    parser.add_argument("--protein", type=Path, default=DEFAULT_PROTEIN)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser


def main() -> None:
    args = _argument_parser().parse_args()
    for path in (
        args.manifest,
        args.identity,
        args.ligand,
        args.site,
        args.preparation,
        args.protein,
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
        identity.accession == ligand.accession == site.accession == preparation.accession
        == protein.accession == manifest.source.accession
    ):
        raise SystemExit("capstone records disagree on accession.")

    ligand_state = next(
        state for state in preparation.states if state.state_id == "reference-neutral"
    )
    protein_state = next(
        state for state in protein.states if state.state_id == "standard-templates"
    )
    positions = {atom.atom_name: atom.position_angstrom for atom in ligand.atoms}
    hydrogens = []
    for atom in ligand_state.atoms:
        if not atom.is_donor:
            continue
        parent = positions.get(atom.atom_id)
        if parent is None:
            continue
        for index, hydrogen in enumerate(atom.hydrogen_positions_nm):
            hydrogen_angstrom = tuple(value * 10.0 for value in hydrogen)
            distance_angstrom = (
                sum((a - b) ** 2 for a, b in zip(parent, hydrogen_angstrom, strict=True))
                ** 0.5
            )
            hydrogens.append(
                HydrogenPosition(
                    hydrogen_id=f"{atom.atom_id}:H{index + 1}",
                    parent_atom_id=atom.atom_id,
                    parent_side="ligand",
                    parent_element=atom.element,
                    position_angstrom=hydrogen_angstrom,
                    parent_position_angstrom=parent,
                    donor_hydrogen_distance_angstrom=distance_angstrom,
                    state_id=ligand_state.state_id,
                    generation_method="RDKit AddHs(addCoords=True)",
                    toolkit=preparation.toolkit,
                    toolkit_version=preparation.toolkit_version,
                    coordinate_source="generated",
                    parent_coordinate_source="deposited",
                    feature_source=atom.feature_source,
                )
            )

    contact_protein_atoms = {
        (contact.auth_asym_id, contact.auth_seq_id, contact.protein_atom)
        for contact in site.contacts
    }
    acceptors = [
        atom
        for atom in protein_state.atoms
        if atom.is_acceptor
        and (atom.auth_asym_id, atom.auth_seq_id, atom.atom_name) in contact_protein_atoms
    ]
    criteria = HydrogenBondCriteria()
    bonds = assess_ligand_donor_hydrogen_bonds(
        hydrogens,
        acceptors,
        criteria,
    )
    record = build_hydrogen_evidence_record(
        source_manifest_path=str(args.manifest),
        source_manifest_sha256=sha256_file(args.manifest),
        ligand_identity_record_path=str(args.identity),
        ligand_identity_record_sha256=sha256_file(args.identity),
        ligand_coordinate_record_path=str(args.ligand),
        ligand_coordinate_record_sha256=sha256_file(args.ligand),
        binding_site_record_path=str(args.site),
        binding_site_record_sha256=sha256_file(args.site),
        chemical_preparation_record_path=str(args.preparation),
        chemical_preparation_record_sha256=sha256_file(args.preparation),
        protein_preparation_record_path=str(args.protein),
        protein_preparation_record_sha256=sha256_file(args.protein),
        accession=manifest.source.accession,
        protein_selection=manifest.selection.protein,
        ligand_selection=identity.selection,
        ligand_state_id=ligand_state.state_id,
        protein_state_id=protein_state.state_id,
        criteria=criteria,
        ligand_hydrogens=tuple(hydrogens),
        protein_hydrogen_status="not_run",
        bonds=bonds,
    )
    errors = validate_hydrogen_evidence(record)
    if errors:
        for error in errors:
            print(f"hydrogen evidence error: {error}")
        raise SystemExit(2)

    destination = write_hydrogen_evidence(args.output, record)
    supported = sum(1 for item in bonds if item.status == "supported")
    not_supported = sum(1 for item in bonds if item.status == "not_supported")
    print(f"{record.accession} hydrogen evidence")
    print(
        f"  ligand state: {record.ligand_state_id}; "
        f"{len(record.ligand_hydrogens)} generated hydrogens"
    )
    print(
        f"  protein state: {record.protein_state_id}; "
        f"protein hydrogens: {record.protein_hydrogen_status}"
    )
    print(
        f"  donor-acceptor candidates: {len(bonds)}; "
        f"supported={supported}, not_supported={not_supported}"
    )
    print(f"  record: {destination}")
    print("  scope: generated hydrogens are not experimental; no affinity is computed")


if __name__ == "__main__":
    main()
