"""Reassess interaction conclusions under explicit analysis assumptions.

This is Batch 13 of the capstone.  It varies contact cutoffs, hydrogen-bond
angular criteria, generated-hydrogen orientations, protein/ligand chemical
states, and water-bridge cutoffs independently.  It reports invariant,
assumption-sensitive, unsupported, and unknown conclusions without combining
them into a confidence score or an affinity estimate.

Run from the repository root:

    uv run python examples/drug_discovery/capstone_interaction_robustness.py
"""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path
from types import SimpleNamespace

from moleculab import (
    AtomFeature,
    HydrogenBondCriteria,
    InteractionCriteria,
    assess_combination,
    assess_contact_cutoff_sensitivity,
    assess_hydrogen_angle_sensitivity,
    assess_hydrogen_placement_sensitivity,
    assess_water_cutoff_sensitivity,
    build_interaction_robustness_record,
    contact_atom_sets,
    hydrogen_candidate_geometries,
    load_binding_site,
    load_chemical_state_enumeration,
    load_hydrogen_evidence,
    load_ligand_coordinates,
    load_ligand_identity,
    load_manifest,
    load_protein_preparation,
    load_water_evidence,
    sha256_file,
    summarize_state_sensitivity,
    validate_interaction_robustness,
    write_interaction_robustness,
)

DEFAULT_MANIFEST = Path("data/structures/2HYY.manifest.json")
DEFAULT_IDENTITY = Path("data/structures/2HYY.ligand-identity.json")
DEFAULT_LIGAND = Path("data/structures/2HYY.ligand-coordinates.json")
DEFAULT_SITE = Path("data/structures/2HYY.binding-site.json")
DEFAULT_STATES = Path("data/structures/2HYY.chemical-states.json")
DEFAULT_PROTEIN = Path("data/structures/2HYY.protein-preparation.json")
DEFAULT_HYDROGEN = Path("data/structures/2HYY.hydrogen-evidence.json")
DEFAULT_WATER = Path("data/structures/2HYY.water-evidence.json")
DEFAULT_OUTPUT = Path("data/structures/2HYY.interaction-robustness.json")

CONTACT_CUTOFFS = (3.0, 3.5, 4.0, 5.0)
BASELINE_CONTACT_CUTOFF = 4.0
HYDROGEN_ANGLE_THRESHOLDS = (90.0, 100.0, 110.0, 120.0, 130.0, 140.0, 150.0)
HYDROGEN_PLACEMENT_OFFSETS = (-20.0, -10.0, 0.0, 10.0, 20.0)
BASELINE_HYDROGEN_CRITERIA = HydrogenBondCriteria()
WATER_CUTOFFS = (3.0, 3.25, 3.5, 3.75, 4.0)
BASELINE_WATER_CUTOFF = 3.5


def _argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--identity", type=Path, default=DEFAULT_IDENTITY)
    parser.add_argument("--ligand", type=Path, default=DEFAULT_LIGAND)
    parser.add_argument("--site", type=Path, default=DEFAULT_SITE)
    parser.add_argument("--states", type=Path, default=DEFAULT_STATES)
    parser.add_argument("--protein", type=Path, default=DEFAULT_PROTEIN)
    parser.add_argument("--hydrogen", type=Path, default=DEFAULT_HYDROGEN)
    parser.add_argument("--water", type=Path, default=DEFAULT_WATER)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser


def main() -> None:
    args = _argument_parser().parse_args()
    for path in (
        args.manifest,
        args.identity,
        args.ligand,
        args.site,
        args.states,
        args.protein,
        args.hydrogen,
        args.water,
    ):
        if not path.exists():
            raise SystemExit(f"missing {path}; run the earlier capstone batches first.")

    manifest = load_manifest(args.manifest)
    identity = load_ligand_identity(args.identity)
    ligand = load_ligand_coordinates(args.ligand)
    site = load_binding_site(args.site)
    states = load_chemical_state_enumeration(args.states)
    protein = load_protein_preparation(args.protein)
    hydrogen = load_hydrogen_evidence(args.hydrogen)
    water = load_water_evidence(args.water)
    if not (
        identity.accession
        == ligand.accession
        == site.accession
        == states.accession
        == protein.accession
        == hydrogen.accession
        == water.accession
        == manifest.source.accession
    ):
        raise SystemExit("capstone records disagree on accession.")

    reference_to_deposited = {
        reference_index: deposited_name
        for deposited_name, reference_index in identity.mapping.pairs
    }
    criteria = InteractionCriteria()
    water_features = [
        AtomFeature(
            atom_id=f"water:{item.water_id}",
            side="water",
            element="O",
            position_nm=tuple(value * 0.1 for value in item.position_angstrom),
            residue_id=f"HOH {item.auth_asym_id} {item.auth_seq_id}",
            is_donor=None,
            is_acceptor=None,
            formal_charge=None,
            is_hydrophobic=None,
            aromatic_ring_id=None,
            is_metal=None,
            hydrogen_positions_nm=(),
            feature_source=None,
        )
        for item in water.waters
    ]

    baseline_ligand_contacts, baseline_protein_contacts = contact_atom_sets(
        ligand.atoms,
        next(
            state for state in protein.states if state.state_id == hydrogen.protein_state_id
        ).atoms,
        BASELINE_CONTACT_CUTOFF,
    )
    baseline_combinations = tuple(
        assess_combination(
            ligand_state=ligand_state,
            protein_state=protein_state,
            reference_to_deposited=reference_to_deposited,
            ligand_atoms=ligand.atoms,
            protein_atoms=protein_state.atoms,
            contact_ligand_atoms=baseline_ligand_contacts,
            contact_protein_atoms=baseline_protein_contacts,
            water_features=water_features,
            criteria=criteria,
        )
        for ligand_state in states.states
        for protein_state in protein.states
    )
    state_sensitivity = summarize_state_sensitivity(
        baseline_combinations,
        reference_ligand_state_id=hydrogen.ligand_state_id,
        reference_protein_state_id=hydrogen.protein_state_id,
    )

    contact_results = assess_contact_cutoff_sensitivity(
        ligand_states=states.states,
        protein_states=protein.states,
        reference_to_deposited=reference_to_deposited,
        ligand_atoms=ligand.atoms,
        water_features=water_features,
        criteria=criteria,
        cutoffs_angstrom=CONTACT_CUTOFFS,
    )

    hydrogen_protein_state = next(
        state for state in protein.states if state.state_id == hydrogen.protein_state_id
    )
    contact_protein_atoms = {
        (contact.auth_asym_id, contact.auth_seq_id, contact.protein_atom)
        for contact in site.contacts
    }
    acceptors = [
        atom
        for atom in hydrogen_protein_state.atoms
        if atom.is_acceptor
        and (atom.auth_asym_id, atom.auth_seq_id, atom.atom_name)
        in contact_protein_atoms
    ]
    hydrogen_candidates = hydrogen_candidate_geometries(
        hydrogen.ligand_hydrogens,
        acceptors,
        hydrogen.criteria,
    )
    hydrogen_angle_results = assess_hydrogen_angle_sensitivity(
        hydrogen.ligand_hydrogens,
        acceptors,
        hydrogen.criteria,
        HYDROGEN_ANGLE_THRESHOLDS,
    )
    hydrogen_placement_results = assess_hydrogen_placement_sensitivity(
        hydrogen.ligand_hydrogens,
        acceptors,
        hydrogen.criteria,
        HYDROGEN_PLACEMENT_OFFSETS,
        modest_band_degrees=10.0,
    )

    ligand_state = next(
        state for state in states.states if state.state_id == hydrogen.ligand_state_id
    )
    ligand_positions = {
        atom.atom_name: atom.position_angstrom for atom in ligand.atoms
    }
    ligand_polar_atoms = [
        SimpleNamespace(
            atom_id=atom.atom_id,
            position_angstrom=ligand_positions[atom.atom_id],
        )
        for atom in ligand_state.atoms
        if (atom.is_donor or atom.is_acceptor) and atom.atom_id in ligand_positions
    ]
    protein_polar_atoms = [
        atom
        for atom in hydrogen_protein_state.atoms
        if (atom.is_donor or atom.is_acceptor)
        and (atom.auth_asym_id, atom.auth_seq_id, atom.atom_name)
        in contact_protein_atoms
    ]
    water_cutoff_results = assess_water_cutoff_sensitivity(
        water.waters,
        ligand_polar_atoms,
        protein_polar_atoms,
        WATER_CUTOFFS,
    )

    water_roles_assigned = any(
        item.donor_role is not None
        or item.acceptor_role is not None
        or item.orientation_known
        for item in water.waters
    )
    record = build_interaction_robustness_record(
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
        hydrogen_evidence_record_path=str(args.hydrogen),
        hydrogen_evidence_record_sha256=sha256_file(args.hydrogen),
        water_evidence_record_path=str(args.water),
        water_evidence_record_sha256=sha256_file(args.water),
        accession=manifest.source.accession,
        protein_selection=manifest.selection.protein,
        ligand_selection=identity.selection,
        baseline_contact_cutoff_angstrom=BASELINE_CONTACT_CUTOFF,
        baseline_hydrogen_criteria=BASELINE_HYDROGEN_CRITERIA,
        baseline_water_cutoff_angstrom=BASELINE_WATER_CUTOFF,
        protein_hydrogen_status=hydrogen.protein_hydrogen_status,
        water_roles_assigned=water_roles_assigned,
        contact_cutoffs_angstrom=CONTACT_CUTOFFS,
        hydrogen_angle_thresholds_degrees=HYDROGEN_ANGLE_THRESHOLDS,
        hydrogen_placement_offsets_degrees=HYDROGEN_PLACEMENT_OFFSETS,
        water_cutoffs_angstrom=WATER_CUTOFFS,
        contact_cutoff_results=contact_results,
        hydrogen_candidates=hydrogen_candidates,
        hydrogen_angle_results=hydrogen_angle_results,
        hydrogen_placement_results=hydrogen_placement_results,
        water_cutoff_results=water_cutoff_results,
        state_sensitivity=state_sensitivity,
    )
    errors = validate_interaction_robustness(record)
    if errors:
        for error in errors:
            print(f"interaction robustness error: {error}")
        raise SystemExit(2)

    destination = write_interaction_robustness(args.output, record)
    statuses = Counter(item.status for item in record.conclusions)
    print(f"{record.accession} interaction robustness")
    print(
        "  contact cutoffs: "
        + ", ".join(
            f"{item.cutoff_angstrom:g} A -> "
            f"{item.minimum_classified_per_combination}-"
            f"{item.maximum_classified_per_combination} classified"
            for item in record.contact_cutoff_results
        )
    )
    print(
        "  hydrogen angle support: "
        + ", ".join(
            f"{item.min_angle_degrees:g} deg -> {item.supported_candidate_count}"
            for item in record.hydrogen_angle_results
        )
    )
    print(
        "  hydrogen placement: "
        + ", ".join(
            f"{item.candidate_id} {item.modest_band_status}"
            for item in record.hydrogen_placement_results
        )
    )
    print(
        "  water bridge plausibility: "
        + ", ".join(
            f"{item.cutoff_angstrom:g} A -> {item.plausible_count}"
            for item in record.water_cutoff_results
        )
    )
    print(
        f"  state patterns: ligand={record.state_sensitivity.unique_ligand_pattern_count}, "
        f"protein={record.state_sensitivity.unique_protein_pattern_count}"
    )
    print(
        "  conclusions: "
        + ", ".join(
            f"{status}={statuses[status]}"
            for status in (
                "invariant",
                "assumption_sensitive",
                "unsupported",
                "unknown",
            )
        )
    )
    print(f"  record: {destination}")
    print(
        "  scope: independent sensitivity axes; no confidence score, ranking, "
        "or affinity"
    )


if __name__ == "__main__":
    main()
