"""Interpret one experimental binding site under multiple chemical states.

This is Batch 9 of the capstone.  It maps each enumerated chemical state back to
the deposited ligand atoms and records state-specific ligand features and
evidence fingerprints against the same experimental contacts.  It does not rank
states, estimate populations, infer protein chemistry, or call chemical
eligibility a demonstrated interaction.

Run:

    uv run python examples/drug_discovery/capstone_state_interactions.py
"""

from __future__ import annotations

import argparse
from pathlib import Path

from moleculab import (
    StateInteractionCriteria,
    assess_state_features,
    build_state_interaction_record,
    compare_state_assessment,
    load_binding_site,
    load_chemical_state_enumeration,
    load_ligand_coordinates,
    load_ligand_identity,
    load_manifest,
    load_state_interactions,
    sha256_file,
    validate_state_interactions,
    write_state_interactions,
)

DEFAULT_MANIFEST = Path("data/structures/2HYY.manifest.json")
DEFAULT_IDENTITY = Path("data/structures/2HYY.ligand-identity.json")
DEFAULT_LIGAND = Path("data/structures/2HYY.ligand-coordinates.json")
DEFAULT_SITE = Path("data/structures/2HYY.binding-site.json")
DEFAULT_STATES = Path("data/structures/2HYY.chemical-states.json")
DEFAULT_OUTPUT = Path("data/structures/2HYY.state-interactions.json")


def _argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--identity", type=Path, default=DEFAULT_IDENTITY)
    parser.add_argument("--ligand", type=Path, default=DEFAULT_LIGAND)
    parser.add_argument("--site", type=Path, default=DEFAULT_SITE)
    parser.add_argument("--states", type=Path, default=DEFAULT_STATES)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser


def main() -> None:
    args = _argument_parser().parse_args()
    for path in (args.manifest, args.identity, args.ligand, args.site, args.states):
        if not path.exists():
            raise SystemExit(f"missing {path}; run the earlier capstone batches first.")

    manifest = load_manifest(args.manifest)
    identity = load_ligand_identity(args.identity)
    ligand = load_ligand_coordinates(args.ligand)
    site = load_binding_site(args.site)
    states = load_chemical_state_enumeration(args.states)
    if not (
        identity.accession == ligand.accession == site.accession == states.accession
        == manifest.source.accession
    ):
        raise SystemExit("capstone records disagree on accession.")

    reference_to_deposited = {
        reference_index: deposited_name
        for deposited_name, reference_index in identity.mapping.pairs
    }
    contact_ligand_atoms = {contact.ligand_atom for contact in site.contacts}
    criteria = StateInteractionCriteria(
        contact_cutoff_angstrom=site.criteria.contact_cutoff_angstrom,
        geometry_source="Batch 3 ligand coordinates and Batch 4 binding-site contacts",
        feature_source=(
            "state-specific donor/acceptor/hydrophobic/aromatic features from "
            "the Batch 8 enumeration"
        ),
        classification_limitations=(
            "protein donor/acceptor, formal-charge, hydrophobic, ring, and "
            "aligned hydrogen evidence is not supplied in Batch 9"
        ),
    )
    selected_state_ids = tuple(state.state_id for state in states.states)
    assessments = tuple(
        assess_state_features(
            state=state,
            reference_to_deposited=reference_to_deposited,
            ligand_atoms=ligand.atoms,
            contact_ligand_atoms=contact_ligand_atoms,
            criteria=criteria,
        )
        for state in states.states
    )
    reference = next(
        item for item in assessments if item.state_id == states.reference_state_id
    )
    comparisons = tuple(
        compare_state_assessment(reference, item)
        for item in assessments
        if item.state_id != reference.state_id
    )
    record = build_state_interaction_record(
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
        accession=manifest.source.accession,
        selection=states.selection,
        selected_state_ids=selected_state_ids,
        reference_state_id=reference.state_id,
        criteria=criteria,
        assessments=assessments,
        comparisons=comparisons,
    )
    errors = validate_state_interactions(record)
    if errors:
        for error in errors:
            print(f"state interaction record error: {error}")
        raise SystemExit(2)

    destination = write_state_interactions(args.output, record)
    print(f"{record.accession} state-aware interaction assessment")
    print(f"  selected states: {', '.join(record.selected_state_ids)}")
    print(
        f"  contact atoms: {reference.contact_atom_count}; "
        f"reference ligand features: "
        f"{reference.ligand_feature_summary.donor_count} donors, "
        f"{reference.ligand_feature_summary.acceptor_count} acceptors, "
        f"{reference.ligand_feature_summary.hydrophobic_count} hydrophobic, "
        f"{reference.ligand_feature_summary.aromatic_atom_count} aromatic"
    )
    for comparison in record.comparisons:
        ligand_diffs = ", ".join(
            f"{name} {first}->{second}"
            for name, first, second in comparison.ligand_feature_differences
        )
        contact_diffs = ", ".join(
            f"{name} {first}->{second}"
            for name, first, second in comparison.contact_feature_differences
        )
        print(
            f"    {comparison.state_id}: differing categories="
            f"{', '.join(comparison.differing_categories) or 'none'}; "
            f"ligand feature deltas={ligand_diffs or 'none'}; "
            f"contact feature deltas={contact_diffs or 'none'}; "
            f"unknown_involved={comparison.comparison.unknown_involved}"
        )
    print(f"  record: {destination}")
    print("  scope: same experimental geometry; no state ranking, population, or affinity")


if __name__ == "__main__":
    main()
