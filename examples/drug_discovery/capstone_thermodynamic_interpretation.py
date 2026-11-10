"""Interpret structural evidence without inventing a binding free energy.

This is Batch 14 of the capstone.  It consumes the structural, interaction,
hydrogen, water, combination, and robustness records and writes a thermodynamic
interpretation record.  The record names the missing thermodynamic terms and
preparation/sampling requirements, and leaves the binding free energy explicitly
``not_computed``.

Run from the repository root:

    uv run python examples/drug_discovery/capstone_thermodynamic_interpretation.py
"""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

from moleculab import (
    StructuralEvidenceSummary,
    RobustnessStatusCount,
    build_thermodynamic_interpretation_record,
    load_binding_site,
    load_chemical_state_enumeration,
    load_combination_record,
    load_hydrogen_evidence,
    load_interaction_robustness,
    load_ligand_coordinates,
    load_ligand_identity,
    load_manifest,
    load_protein_preparation,
    load_water_evidence,
    sha256_file,
    validate_thermodynamic_interpretation,
    write_thermodynamic_interpretation,
)

DEFAULT_MANIFEST = Path("data/structures/2HYY.manifest.json")
DEFAULT_IDENTITY = Path("data/structures/2HYY.ligand-identity.json")
DEFAULT_LIGAND = Path("data/structures/2HYY.ligand-coordinates.json")
DEFAULT_SITE = Path("data/structures/2HYY.binding-site.json")
DEFAULT_STATES = Path("data/structures/2HYY.chemical-states.json")
DEFAULT_PROTEIN = Path("data/structures/2HYY.protein-preparation.json")
DEFAULT_HYDROGEN = Path("data/structures/2HYY.hydrogen-evidence.json")
DEFAULT_WATER = Path("data/structures/2HYY.water-evidence.json")
DEFAULT_COMBINATION = Path("data/structures/2HYY.combination-interactions.json")
DEFAULT_ROBUSTNESS = Path("data/structures/2HYY.interaction-robustness.json")
DEFAULT_OUTPUT = Path("data/structures/2HYY.thermodynamic-interpretation.json")


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
    parser.add_argument("--combination", type=Path, default=DEFAULT_COMBINATION)
    parser.add_argument("--robustness", type=Path, default=DEFAULT_ROBUSTNESS)
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
        args.combination,
        args.robustness,
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
    combination = load_combination_record(args.combination)
    robustness = load_interaction_robustness(args.robustness)
    if not (
        identity.accession
        == ligand.accession
        == site.accession
        == states.accession
        == protein.accession
        == hydrogen.accession
        == water.accession
        == combination.accession
        == robustness.accession
        == manifest.source.accession
    ):
        raise SystemExit("capstone records disagree on accession.")

    reference_combination = next(
        item
        for item in combination.combinations
        if item.ligand_state_id == hydrogen.ligand_state_id
        and item.protein_state_id == hydrogen.protein_state_id
    )
    classified_counts = Counter(
        item.interaction_type for item in reference_combination.classified_interactions
    )
    robustness_counts = Counter(item.status for item in robustness.conclusions)

    summary = StructuralEvidenceSummary(
        ligand_state_count=len(combination.ligand_state_ids),
        protein_state_count=len(combination.protein_state_ids),
        combination_count=len(combination.combinations),
        ligand_contact_atom_count=len(
            {contact.ligand_atom for contact in site.contacts}
        ),
        protein_contact_atom_count=len(
            {
                (contact.auth_asym_id, contact.auth_seq_id, contact.protein_atom)
                for contact in site.contacts
            }
        ),
        hydrophobic_pair_count=classified_counts.get("hydrophobic", 0),
        aromatic_pair_count=classified_counts.get("aromatic", 0),
        supported_hydrogen_bond_count=sum(
            1 for bond in hydrogen.bonds if bond.status == "supported"
        ),
        plausible_water_bridge_count=sum(
            1
            for bridge in water.bridges
            if bridge.status == "geometrically_plausible"
        ),
        robustness_conclusion_counts=tuple(
            RobustnessStatusCount(status=status, count=robustness_counts.get(status, 0))
            for status in (
                "invariant",
                "assumption_sensitive",
                "unsupported",
                "unknown",
            )
        ),
    )
    record = build_thermodynamic_interpretation_record(
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
        combination_interaction_record_path=str(args.combination),
        combination_interaction_record_sha256=sha256_file(args.combination),
        interaction_robustness_record_path=str(args.robustness),
        interaction_robustness_record_sha256=sha256_file(args.robustness),
        accession=manifest.source.accession,
        protein_selection=manifest.selection.protein,
        ligand_selection=identity.selection,
        reference_ligand_state_id=hydrogen.ligand_state_id,
        reference_protein_state_id=hydrogen.protein_state_id,
        structural_summary=summary,
    )
    errors = validate_thermodynamic_interpretation(record)
    if errors:
        for error in errors:
            print(f"thermodynamic interpretation error: {error}")
        raise SystemExit(2)

    destination = write_thermodynamic_interpretation(args.output, record)
    requirement_statuses = Counter(item.status for item in record.requirements)
    print(f"{record.accession} thermodynamic interpretation")
    print(
        f"  structural evidence: {summary.combination_count} state combinations, "
        f"{summary.hydrophobic_pair_count} hydrophobic pairs, "
        f"{summary.aromatic_pair_count} aromatic pair, "
        f"{summary.supported_hydrogen_bond_count} supported H-bonds, "
        f"{summary.plausible_water_bridge_count} plausible water bridge"
    )
    print(
        f"  robustness conclusions: "
        + ", ".join(
            f"{item.status}={item.count}"
            for item in summary.robustness_conclusion_counts
        )
    )
    print(
        f"  boundaries: {len(record.boundaries)}; "
        f"missing thermodynamic terms: {len(record.terms)}; "
        f"requirements: {len(record.requirements)} "
        f"{dict(requirement_statuses)}"
    )
    print(
        f"  binding free energy: {record.free_energy_status} "
        f"(value={record.binding_free_energy_kj_per_mol})"
    )
    print(f"  record: {destination}")
    print(
        "  scope: structural evidence is not an affinity; no binding free "
        "energy or state ranking is computed"
    )


if __name__ == "__main__":
    main()
