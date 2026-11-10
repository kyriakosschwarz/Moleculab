# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Tests for state-aware interaction assessment."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import pytest

from moleculab import (
    STATE_INTERACTION_SCHEMA_VERSION,
    ChemicalState,
    CoordinateAtom,
    LigandSelection,
    PreparedAtom,
    StateInteractionCriteria,
    assess_state_features,
    build_state_interaction_record,
    compare_state_assessment,
    load_state_interactions,
    validate_state_interactions,
    write_state_interactions,
)

CHECKSUM = "2" * 64


def selection() -> LigandSelection:
    return LigandSelection(
        component_id="LIG",
        label_asym_id="E",
        auth_asym_id="A",
        instance="E",
        auth_seq_id="600",
        available_instances=("E",),
    )


def prepared_atom(
    atom_id: str,
    reference_index: int,
    *,
    donor: bool | None = True,
    acceptor: bool | None = True,
    hydrophobic: bool | None = True,
    charge: int | None = 0,
    aromatic: bool = True,
) -> PreparedAtom:
    return PreparedAtom(
        atom_id=atom_id,
        element="C",
        formal_charge=charge,
        is_donor=donor,
        is_acceptor=acceptor,
        is_hydrophobic=hydrophobic,
        aromatic_ring_ids=("ring-1",) if aromatic else (),
        hydrogen_count=0,
        hydrogen_positions_nm=(),
        coordinate_source="generated",
        hydrogen_source=None,
        feature_source="test",
    )


def chemical_state(
    state_id: str,
    *,
    kind: str = "reference",
    charge: int = 0,
    donor: bool | None = True,
    acceptor: bool | None = True,
    hydrophobic: bool | None = True,
    aromatic: bool = True,
) -> ChemicalState:
    return ChemicalState(
        state_id=state_id,
        kind=kind,
        label=state_id,
        parent_state_id=None if kind == "reference" else "reference",
        molecular_graph="CC",
        formula="C2H6",
        proton_count=31,
        formal_charge=charge,
        tautomer_of="reference" if kind == "tautomer" else None,
        ph_assumption=None,
        pka_sites=(),
        enumeration_method="test",
        enumeration_parameters=("test",),
        assumptions=("test",),
        atoms=(
            prepared_atom(
                "atom_0",
                0,
                donor=donor,
                acceptor=acceptor,
                hydrophobic=hydrophobic,
                charge=charge,
                aromatic=aromatic,
            ),
            prepared_atom(
                "atom_1",
                1,
                donor=donor,
                acceptor=acceptor,
                hydrophobic=hydrophobic,
                charge=charge,
                aromatic=aromatic,
            ),
        ),
        reference_atom_mapping=(("atom_0", 0), ("atom_1", 1)),
        rings=(),
        water_roles=(),
    )


def ligand_atom(name: str, position: tuple[float, float, float], index: int) -> CoordinateAtom:
    return CoordinateAtom(
        atom_name=name,
        element="C",
        reference_index=index,
        position_angstrom=position,
        occupancy=1.0,
        b_iso_or_equiv=20.0,
        alternate_location=None,
        model_number=1,
    )


def criteria() -> StateInteractionCriteria:
    return StateInteractionCriteria(
        contact_cutoff_angstrom=4.0,
        geometry_source="test experimental geometry",
        feature_source="test features",
        classification_limitations="protein chemistry not supplied",
    )


def assessment(state: ChemicalState):
    return assess_state_features(
        state=state,
        reference_to_deposited={0: "C1", 1: "C2"},
        ligand_atoms=(
            ligand_atom("C1", (0.0, 0.0, 0.0), 0),
            ligand_atom("C2", (1.4, 0.0, 0.0), 1),
        ),
        contact_ligand_atoms={"C1"},
        criteria=criteria(),
    )


def test_state_features_and_contacts_are_assessed_separately():
    result = assessment(chemical_state("reference"))

    assert result.ligand_feature_summary.atom_count == 2
    assert result.contact_feature_summary.atom_count == 1
    assert result.contact_atom_count == 1
    assert result.fingerprint.state_for(
        result.fingerprint.schema.features[0]
    ) in ("present", "absent", "unknown")


def test_protonation_changes_ionic_evidence():
    neutral = assessment(chemical_state("reference", charge=0))
    protonated = assessment(chemical_state("protonation-1", kind="protonation", charge=1))

    assert neutral.fingerprint.state_for(
        next(f for f in neutral.fingerprint.schema.features if f.interaction_type == "ligand_ionic")
    ) == "absent"
    assert protonated.fingerprint.state_for(
        next(f for f in protonated.fingerprint.schema.features if f.interaction_type == "ligand_ionic")
    ) == "present"

    comparison = compare_state_assessment(neutral, protonated)
    assert "ligand_ionic" in comparison.differing_categories
    assert ("charged_atom_count", 0, 2) in comparison.ligand_feature_differences


def test_unknown_chemistry_is_not_absent():
    result = assessment(
        chemical_state(
            "reference",
            donor=None,
            acceptor=None,
            hydrophobic=None,
        )
    )
    for interaction_type in ("ligand_donor_acceptor", "ligand_hydrophobic"):
        feature = next(
            feature
            for feature in result.fingerprint.schema.features
            if feature.interaction_type == interaction_type
        )
        assert result.fingerprint.state_for(feature) == "unknown"
        assert result.fingerprint.assessable_for(feature) is False


def build_record():
    reference = assessment(chemical_state("reference"))
    protonated = assessment(chemical_state("protonation-1", kind="protonation", charge=1))
    comparison = compare_state_assessment(reference, protonated)
    return build_state_interaction_record(
        source_manifest_path="manifest.json",
        source_manifest_sha256=CHECKSUM,
        ligand_identity_record_path="identity.json",
        ligand_identity_record_sha256=CHECKSUM,
        ligand_coordinate_record_path="coordinates.json",
        ligand_coordinate_record_sha256=CHECKSUM,
        binding_site_record_path="site.json",
        binding_site_record_sha256=CHECKSUM,
        chemical_state_record_path="states.json",
        chemical_state_record_sha256=CHECKSUM,
        accession="TEST",
        selection=selection(),
        selected_state_ids=("reference", "protonation-1"),
        reference_state_id="reference",
        criteria=criteria(),
        assessments=(reference, protonated),
        comparisons=(comparison,),
    )


def test_record_round_trip_and_validation(tmp_path):
    result = build_record()
    path = write_state_interactions(tmp_path / "state-interactions.json", result)
    loaded = load_state_interactions(path)

    assert loaded == result
    assert validate_state_interactions(loaded) == []
    assert loaded.schema_version == STATE_INTERACTION_SCHEMA_VERSION


def test_selected_state_must_be_assessed():
    result = build_record()
    bad = replace(result, selected_state_ids=("missing",))
    errors = validate_state_interactions(bad)

    assert any("selected_state_ids" in error for error in errors)


def test_reference_comparison_records_feature_deltas():
    result = build_record()
    comparison = result.comparisons[0]

    assert comparison.state_id == "protonation-1"
    assert "ligand_ionic" in comparison.differing_categories
    assert comparison.ligand_feature_differences


def test_reference_2hyy_state_interaction_record_is_conservative():
    path = (
        Path(__file__).resolve().parents[1]
        / "data"
        / "structures"
        / "2HYY.state-interactions.json"
    )
    result = load_state_interactions(path)

    assert validate_state_interactions(result) == []
    assert result.accession == "2HYY"
    assert len(result.assessments) == 8
    assert len(result.comparisons) == 7
    assert result.reference_state_id == "reference-neutral"
    assert all(item.contact_atom_count == 34 for item in result.assessments)
    reference = next(item for item in result.assessments if item.state_id == "reference-neutral")
    assert reference.ligand_feature_summary.donor_count == 4
    assert reference.ligand_feature_summary.acceptor_count == 4
    classification = next(
        feature
        for feature in reference.fingerprint.schema.features
        if feature.interaction_type == "interaction_classification"
    )
    assert reference.fingerprint.state_for(classification) == "unknown"
    assert not reference.fingerprint.assessable_for(classification)
    protonation = next(
        item for item in result.comparisons if item.state_id == "protonation-1"
    )
    assert "ligand_ionic" in protonation.differing_categories
