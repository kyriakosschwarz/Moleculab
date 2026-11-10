# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Assess every ligand/protein chemical-state combination on fixed contacts.

This is Batch 11 of the capstone.  It crosses eight ligand states with two
protein-pocket states and produces 16 separate evidence results against the same
experimental contact geometry.  It does not rank combinations, estimate
populations, or infer affinity.

Hydrogen-bond and water-mediated classifications remain unknown because aligned
hydrogen positions and water roles are not available.

Run with gemmi available:

    uv run --with gemmi python examples/drug_discovery/capstone_combination_interactions.py
"""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

import gemmi

from moleculab import (
    AtomFeature,
    InteractionCriteria,
    assess_combination,
    build_combination_record,
    load_binding_site,
    load_chemical_state_enumeration,
    load_combination_record,
    load_ligand_coordinates,
    load_ligand_identity,
    load_manifest,
    load_protein_preparation,
    sha256_file,
    validate_combination_record,
    write_combination_record,
)

DEFAULT_MANIFEST = Path("data/structures/2HYY.manifest.json")
DEFAULT_IDENTITY = Path("data/structures/2HYY.ligand-identity.json")
DEFAULT_LIGAND = Path("data/structures/2HYY.ligand-coordinates.json")
DEFAULT_SITE = Path("data/structures/2HYY.binding-site.json")
DEFAULT_STATES = Path("data/structures/2HYY.chemical-states.json")
DEFAULT_PROTEIN = Path("data/structures/2HYY.protein-preparation.json")
DEFAULT_SOURCE = Path("data/structures/2HYY.cif")
DEFAULT_OUTPUT = Path("data/structures/2HYY.combination-interactions.json")


def _argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--identity", type=Path, default=DEFAULT_IDENTITY)
    parser.add_argument("--ligand", type=Path, default=DEFAULT_LIGAND)
    parser.add_argument("--site", type=Path, default=DEFAULT_SITE)
    parser.add_argument("--states", type=Path, default=DEFAULT_STATES)
    parser.add_argument("--protein", type=Path, default=DEFAULT_PROTEIN)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser


def _unquote(value: str) -> str:
    value = value.strip()
    if len(value) >= 2 and value[0] == value[-1] and value[0] in ("'", '"'):
        return value[1:-1]
    return value


def _water_lookup(block: gemmi.cif.Block, auth_chain: str):
    lookup = {}
    for row in block.find(
        "_atom_site.",
        [
            "group_PDB",
            "label_comp_id",
            "auth_asym_id",
            "auth_seq_id",
            "label_atom_id",
            "Cartn_x",
            "Cartn_y",
            "Cartn_z",
        ],
    ):
        values = [_unquote(value) for value in row]
        if values[0] == "HETATM" and values[1] == "HOH" and values[2] == auth_chain:
            lookup[(values[3], values[4])] = (
                float(values[5]),
                float(values[6]),
                float(values[7]),
            )
    return lookup


def main() -> None:
    args = _argument_parser().parse_args()
    for path in (
        args.manifest,
        args.identity,
        args.ligand,
        args.site,
        args.states,
        args.protein,
        args.source,
    ):
        if not path.exists():
            raise SystemExit(f"missing {path}; run the earlier capstone batches first.")

    manifest = load_manifest(args.manifest)
    identity = load_ligand_identity(args.identity)
    ligand = load_ligand_coordinates(args.ligand)
    site = load_binding_site(args.site)
    states = load_chemical_state_enumeration(args.states)
    protein = load_protein_preparation(args.protein)
    if not (
        identity.accession == ligand.accession == site.accession == states.accession
        == protein.accession == manifest.source.accession
    ):
        raise SystemExit("capstone records disagree on accession.")

    block = gemmi.cif.read_file(str(args.source)).sole_block()
    water_lookup = _water_lookup(block, manifest.selection.protein.auth_asym_id)
    water_features = []
    for water in site.waters_near_ligand:
        position = water_lookup.get((water.auth_seq_id, water.atom_name))
        if position is None:
            continue
        water_features.append(
            AtomFeature(
                atom_id=f"water:{water.auth_seq_id}:{water.atom_name}",
                side="water",
                element="O",
                position_nm=tuple(value * 0.1 for value in position),
                residue_id=f"HOH {water.auth_asym_id} {water.auth_seq_id}",
                is_donor=None,
                is_acceptor=None,
                formal_charge=None,
                is_hydrophobic=None,
                aromatic_ring_id=None,
                is_metal=None,
                hydrogen_positions_nm=(),
                feature_source=None,
            )
        )

    reference_to_deposited = {
        reference_index: deposited_name
        for deposited_name, reference_index in identity.mapping.pairs
    }
    contact_ligand_atoms = {contact.ligand_atom for contact in site.contacts}
    contact_protein_atoms = {
        (contact.auth_asym_id, contact.auth_seq_id, contact.protein_atom)
        for contact in site.contacts
    }
    criteria = InteractionCriteria()
    combinations = tuple(
        assess_combination(
            ligand_state=ligand_state,
            protein_state=protein_state,
            reference_to_deposited=reference_to_deposited,
            ligand_atoms=ligand.atoms,
            protein_atoms=protein_state.atoms,
            contact_ligand_atoms=contact_ligand_atoms,
            contact_protein_atoms=contact_protein_atoms,
            water_features=water_features,
            criteria=criteria,
        )
        for ligand_state in states.states
        for protein_state in protein.states
    )
    record = build_combination_record(
        source_manifest_path=str(args.manifest),
        source_manifest_sha256=sha256_file(args.manifest),
        ligand_identity_record_path=str(args.identity),
        ligand_identity_record_sha256=sha256_file(args.identity),
        ligand_coordinate_record_path=str(args.ligand),
        ligand_coordinate_record_sha256=sha256_file(args.ligand),
        binding_site_record_path=str(args.site),
        binding_site_record_sha256=sha256_file(args.site),
        chemical_state_record_path=str(args.states),
        chemical_state_record_sha256=sha256_file(args.states),
        protein_preparation_record_path=str(args.protein),
        protein_preparation_record_sha256=sha256_file(args.protein),
        accession=manifest.source.accession,
        protein_selection=manifest.selection.protein,
        ligand_selection=identity.selection,
        criteria=criteria,
        ligand_state_ids=tuple(state.state_id for state in states.states),
        protein_state_ids=tuple(state.state_id for state in protein.states),
        combinations=combinations,
    )
    errors = validate_combination_record(record)
    if errors:
        for error in errors:
            print(f"combination record error: {error}")
        raise SystemExit(2)

    destination = write_combination_record(args.output, record)
    print(f"{record.accession} evidence-complete interaction matrix")
    print(
        f"  states: {len(record.ligand_state_ids)} ligand x "
        f"{len(record.protein_state_ids)} protein = {len(record.combinations)} combinations"
    )
    for combination in record.combinations:
        classified_types = Counter(
            item.interaction_type for item in combination.classified_interactions
        )
        statuses = Counter(status for _, status in combination.category_statuses)
        ligand = combination.ligand_feature_summary
        protein = combination.protein_feature_summary
        print(
            f"    {combination.ligand_state_id} x {combination.protein_state_id}: "
            f"classified={len(combination.classified_interactions)} "
            f"{dict(classified_types)}; "
            f"ligand D/A/H/Q={ligand.donor_count}/{ligand.acceptor_count}/"
            f"{ligand.hydrophobic_count}/{ligand.charged_atom_count}; "
            f"protein D/A/H/Q={protein.donor_count}/{protein.acceptor_count}/"
            f"{protein.hydrophobic_count}/{protein.charged_atom_count}; "
            f"unknown={statuses['unknown']}, "
            f"not_applicable={statuses['not_applicable']}"
        )
    print(f"  record: {destination}")
    print("  scope: fixed geometry; no ranking, population, or affinity")


if __name__ == "__main__":
    main()
