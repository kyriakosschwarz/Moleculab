"""Tests for interaction sensitivity and robustness records."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import pytest

from moleculab import (
    ROBUSTNESS_SCHEMA_VERSION,
    HydrogenBondCriteria,
    HydrogenPosition,
    WaterAtomRecord,
    assess_hydrogen_angle_sensitivity,
    assess_hydrogen_placement_sensitivity,
    assess_water_cutoff_sensitivity,
    contact_atom_sets,
    load_interaction_robustness,
    validate_interaction_robustness,
    write_interaction_robustness,
)


def hydrogen(
    hydrogen_id: str,
    donor: tuple[float, float, float],
    position: tuple[float, float, float],
    parent_atom_id: str = "N1",
) -> HydrogenPosition:
    distance = sum(
        (a - b) ** 2 for a, b in zip(position, donor, strict=True)
    ) ** 0.5
    return HydrogenPosition(
        hydrogen_id=hydrogen_id,
        parent_atom_id=parent_atom_id,
        parent_side="ligand",
        parent_element="N",
        position_angstrom=position,
        parent_position_angstrom=donor,
        donor_hydrogen_distance_angstrom=distance,
        state_id="reference-neutral",
        generation_method="test",
        toolkit="test",
        toolkit_version="1.0",
        coordinate_source="generated",
        parent_coordinate_source="deposited",
        feature_source="test",
    )


def acceptor(atom_id: str, position: tuple[float, float, float]):
    return SimpleNamespace(atom_id=atom_id, position_angstrom=position)


def water(position: tuple[float, float, float] = (0.0, 0.0, 0.0)):
    return WaterAtomRecord(
        water_id="A:101:O",
        label_asym_id="A",
        auth_asym_id="A",
        label_seq_id=None,
        auth_seq_id="101",
        atom_name="O",
        position_angstrom=position,
        occupancy=1.0,
        b_iso_or_equiv=25.0,
        orientation_known=False,
        donor_role=None,
        acceptor_role=None,
        source="deposited HETATM HOH record",
    )


def test_contact_atom_sets_are_cutoff_inclusive():
    ligand_atoms = (
        SimpleNamespace(atom_name="C1", position_angstrom=(0.0, 0.0, 0.0)),
        SimpleNamespace(atom_name="C2", position_angstrom=(20.0, 0.0, 0.0)),
    )
    protein_atoms = (
        SimpleNamespace(
            auth_asym_id="A",
            auth_seq_id="1",
            atom_name="N",
            position_angstrom=(3.0, 0.0, 0.0),
        ),
        SimpleNamespace(
            auth_asym_id="A",
            auth_seq_id="1",
            atom_name="O",
            position_angstrom=(5.0, 0.0, 0.0),
        ),
    )

    ligand, protein = contact_atom_sets(ligand_atoms, protein_atoms, 3.0)
    assert ligand == {"C1"}
    assert protein == {("A", "1", "N")}

    ligand, protein = contact_atom_sets(ligand_atoms, protein_atoms, 5.0)
    assert ligand == {"C1"}
    assert protein == {("A", "1", "N"), ("A", "1", "O")}


def test_hydrogen_angle_sensitivity_reports_threshold_dependence():
    hydrogens = (
        hydrogen("D1:H1", (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), "D1"),
        hydrogen(
            "D2:H1",
            (10.0, 0.0, 0.0),
            (10.866, 0.5, 0.0),
            "D2",
        ),
    )
    acceptors = (
        acceptor("A1", (2.0, 0.0, 0.0)),
        acceptor("A2", (12.0, 0.0, 0.0)),
    )
    criteria = HydrogenBondCriteria()

    results = assess_hydrogen_angle_sensitivity(
        hydrogens,
        acceptors,
        criteria,
        (90.0, 120.0, 150.0),
    )

    assert [item.supported_candidate_count for item in results] == [2, 2, 1]
    assert results[0].supported_candidate_ids == ("D1->A1", "D2->A2")
    assert results[2].supported_candidate_ids == ("D1->A1",)


def test_hydrogen_placement_sensitivity_distinguishes_robust_and_sensitive():
    hydrogens = (
        hydrogen("D1:H1", (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), "D1"),
        hydrogen(
            "D2:H1",
            (10.0, 0.0, 0.0),
            (10.866, 0.5, 0.0),
            "D2",
        ),
    )
    acceptors = (
        acceptor("A1", (2.0, 0.0, 0.0)),
        acceptor("A2", (12.0, 0.0, 0.0)),
    )
    criteria = HydrogenBondCriteria()

    results = assess_hydrogen_placement_sensitivity(
        hydrogens,
        acceptors,
        criteria,
        (-20.0, -10.0, 0.0, 10.0, 20.0),
        modest_band_degrees=10.0,
    )

    first, second = results
    assert first.current_status == "supported"
    assert first.modest_band_status == "invariant_supported"
    assert second.modest_band_status == "placement_sensitive"
    assert second.modest_band_supported_offsets_degrees == (-10.0, 0.0)


def test_water_cutoff_sensitivity_preserves_unknown_separation():
    polar_ligand = (SimpleNamespace(atom_id="C1", position_angstrom=(3.0, 0.0, 0.0)),)
    polar_protein = (
        SimpleNamespace(atom_id="A:1:O", position_angstrom=(-3.0, 0.0, 0.0)),
    )

    results = assess_water_cutoff_sensitivity(
        (water(),),
        polar_ligand,
        polar_protein,
        (2.9, 3.0, 3.5),
    )

    assert [item.plausible_count for item in results] == [0, 1, 1]
    assert all(
        item.unknown_count == 0 for item in results
    )


def test_reference_2hyy_robustness_record_is_conservative():
    path = (
        Path(__file__).resolve().parents[1]
        / "data"
        / "structures"
        / "2HYY.interaction-robustness.json"
    )
    result = load_interaction_robustness(path)

    assert validate_interaction_robustness(result) == []
    assert result.schema_version == ROBUSTNESS_SCHEMA_VERSION
    assert result.accession == "2HYY"
    assert len(result.contact_cutoff_results) == 4
    assert len(result.hydrogen_candidates) == 3
    assert len(result.hydrogen_angle_results) == 7
    assert len(result.hydrogen_placement_results) == 3
    assert len(result.water_cutoff_results) == 5
    assert len(result.conclusions) == 17

    contact = {
        item.cutoff_angstrom: item for item in result.contact_cutoff_results
    }
    assert contact[3.0].minimum_classified_per_combination == 1
    assert contact[3.5].minimum_classified_per_combination == 2
    assert contact[4.0].minimum_classified_per_combination == 19
    assert contact[5.0].minimum_classified_per_combination == 19

    hydrophobic = {
        item.cutoff_angstrom: next(
            (
                status.combination_count
                for status in item.category_status_counts
                if status.interaction_type == "hydrophobic"
                and status.status == "classified"
            ),
            0,
        )
        for item in result.contact_cutoff_results
    }
    assert hydrophobic[3.0] == 0
    assert hydrophobic[3.5] == 16
    assert hydrophobic[4.0] == 16

    placement = {
        item.candidate_id: item for item in result.hydrogen_placement_results
    }
    assert placement["N13->A:315:OG1"].modest_band_status == "placement_sensitive"
    assert placement["N21->A:286:OE2"].modest_band_status == "invariant_supported"
    assert placement["N21->A:381:O"].modest_band_status == "invariant_unsupported"

    water = {
        item.cutoff_angstrom: item.plausible_count
        for item in result.water_cutoff_results
    }
    assert water[3.0] == 0
    assert water[3.25] == 1
    assert water[3.5] == 1

    statuses = {item.conclusion_id: item.status for item in result.conclusions}
    assert statuses["contact_aromatic"] == "invariant"
    assert statuses["contact_hydrophobic"] == "assumption_sensitive"
    assert statuses["contact_ionic"] == "unsupported"
    assert statuses["protein_state_pattern"] == "invariant"
    assert statuses["ligand_state_category_counts"] == "invariant"
    assert statuses["ligand_state_specific_pairs"] == "assumption_sensitive"
    assert statuses["water_bridge_geometry"] == "assumption_sensitive"
    assert statuses["water_chemical_classification"] == "unknown"
    assert statuses["protein_donor_hydrogen"] == "unknown"

    assert result.state_sensitivity.unique_protein_pattern_count == 1
    assert result.state_sensitivity.protein_states_differing_from_reference_pattern == ()
    assert result.state_sensitivity.unique_ligand_pattern_count == 2
    assert result.state_sensitivity.ligand_states_matching_reference_pattern == (
        "reference-neutral",
    )
    assert result.water_roles_assigned is False
    assert result.protein_hydrogen_status == "not_run"


def test_conclusion_statuses_have_consistent_support_sets():
    path = (
        Path(__file__).resolve().parents[1]
        / "data"
        / "structures"
        / "2HYY.interaction-robustness.json"
    )
    result = load_interaction_robustness(path)

    for conclusion in result.conclusions:
        if conclusion.status == "invariant":
            assert conclusion.supporting_conditions
            assert conclusion.contradicting_conditions == ()
        elif conclusion.status == "assumption_sensitive":
            assert conclusion.supporting_conditions
            assert conclusion.contradicting_conditions
        elif conclusion.status == "unsupported":
            assert conclusion.supporting_conditions == ()
        else:
            assert conclusion.status == "unknown"
            assert conclusion.limitations


def test_robustness_record_round_trip(tmp_path):
    path = (
        Path(__file__).resolve().parents[1]
        / "data"
        / "structures"
        / "2HYY.interaction-robustness.json"
    )
    result = load_interaction_robustness(path)
    destination = write_interaction_robustness(tmp_path / "robustness.json", result)
    loaded = load_interaction_robustness(destination)

    assert loaded == result
    assert loaded.to_dict()["state_sensitivity"]["unique_ligand_pattern_count"] == 2


def test_invalid_conclusion_support_sets_are_rejected():
    path = (
        Path(__file__).resolve().parents[1]
        / "data"
        / "structures"
        / "2HYY.interaction-robustness.json"
    )
    result = load_interaction_robustness(path)
    original = result.conclusions[0]

    bad_invariant = replace(
        original,
        status="invariant",
        contradicting_conditions=("contradiction",),
    )
    errors = validate_interaction_robustness(
        replace(result, conclusions=(bad_invariant,) + result.conclusions[1:])
    )
    assert any("invariant" in error for error in errors)

    bad_sensitive = replace(
        original,
        status="assumption_sensitive",
        supporting_conditions=(),
    )
    errors = validate_interaction_robustness(
        replace(result, conclusions=(bad_sensitive,) + result.conclusions[1:])
    )
    assert any("assumption-sensitive" in error for error in errors)

    bad_unsupported = replace(
        original,
        status="unsupported",
        supporting_conditions=("support",),
    )
    errors = validate_interaction_robustness(
        replace(result, conclusions=(bad_unsupported,) + result.conclusions[1:])
    )
    assert any("unsupported" in error for error in errors)


def test_robustness_record_contains_no_composite_score():
    path = (
        Path(__file__).resolve().parents[1]
        / "data"
        / "structures"
        / "2HYY.interaction-robustness.json"
    )
    result = load_interaction_robustness(path)
    serialized = result.to_dict()

    assert "confidence" not in serialized
    assert "confidence_score" not in serialized
    assert "affinity" not in serialized
    assert "rank" not in serialized
    assert all(
        "confidence" not in item.to_dict() for item in result.conclusions
    )
