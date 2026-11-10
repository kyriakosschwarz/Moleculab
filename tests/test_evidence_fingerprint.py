"""Tests for three-state evidence fingerprints and interaction summaries."""

from __future__ import annotations

from pathlib import Path

import pytest

from moleculab import (
    EvidenceFeature,
    EvidenceObservation,
    aggregate_evidence,
    build_evidence_fingerprint,
    compare_evidence_fingerprints,
    criteria_fingerprint,
    load_interaction_evidence,
    validate_interaction_evidence,
    write_interaction_evidence,
)


def atom_feature(interaction_type: str = "hydrogen_bond", water_id: str = "") -> EvidenceFeature:
    return EvidenceFeature(
        scope="atom_pair",
        interaction_type=interaction_type,
        protein_residue="LEU A 1",
        protein_atom="N",
        ligand_atom="O1",
        ligand_residue="LIG E 1",
        water_id=water_id,
    )


def category_feature(interaction_type: str) -> EvidenceFeature:
    return EvidenceFeature(scope="category", interaction_type=interaction_type)


def observation(feature, state: str, *, evidence: str = "test observation") -> EvidenceObservation:
    return EvidenceObservation(
        feature=feature,
        state=state,
        assessable=state != "unknown",
        evidence=evidence,
        criteria="explicit test criteria",
        provenance="test fixture",
    )


def fingerprint(features_states):
    observations = [
        observation(feature, state) for feature, state in features_states
    ]
    return build_evidence_fingerprint(
        observations,
        schema_id="test-schema",
        criteria_fingerprint=criteria_fingerprint({"test": "criteria"}),
    )


def test_feature_identity_separates_category_and_atom_pair_scopes():
    category = category_feature("ionic")
    direct = atom_feature("hydrogen_bond")
    water = atom_feature("water_mediated", water_id="W1")

    assert category.key[0] == "category"
    assert direct.key != water.key
    assert water.is_water_mediated
    with pytest.raises(ValueError, match="atom_pair"):
        EvidenceFeature(scope="atom_pair", interaction_type="ionic")


def test_unknown_is_not_encoded_as_absent():
    feature = atom_feature()
    result = fingerprint([(feature, "unknown")])

    assert result.states == ("unknown",)
    assert result.unknown_count == 1
    assert result.absent_count == 0
    assert result.assessable_count == 0


def test_present_absent_and_unknown_coexist():
    present = atom_feature("hydrogen_bond")
    absent = atom_feature("ionic")
    unknown = atom_feature("hydrophobic")
    result = fingerprint(
        [(present, "present"), (absent, "absent"), (unknown, "unknown")]
    )

    assert result.present_count == 1
    assert result.absent_count == 1
    assert result.unknown_count == 1
    assert result.assessable_count == 2


def test_duplicate_observations_merge_conservatively():
    feature = atom_feature()
    result = fingerprint(
        [
            (feature, "unknown"),
            (feature, "absent"),
            (feature, "unknown"),
        ]
    )
    assert result.states == ("unknown",)

    result = fingerprint(
        [
            (feature, "absent"),
            (feature, "present"),
        ]
    )
    assert result.states == ("present",)


def test_occupancy_uses_only_assessable_observations():
    feature = atom_feature()
    first = fingerprint([(feature, "unknown")])
    second = fingerprint([(feature, "present")])

    occupancy = aggregate_evidence([first, second])

    assert occupancy.present_counts == (1,)
    assert occupancy.unknown_counts == (1,)
    assert occupancy.assessable_counts == (1,)
    assert occupancy.present_fractions == (1.0,)


def test_comparison_does_not_treat_unknown_as_absent():
    first_feature = atom_feature("hydrogen_bond")
    second_feature = atom_feature("ionic")
    third_feature = atom_feature("hydrophobic")
    first = fingerprint(
        [
            (first_feature, "present"),
            (second_feature, "absent"),
            (third_feature, "unknown"),
        ]
    )
    second = fingerprint(
        [
            (first_feature, "absent"),
            (second_feature, "present"),
            (third_feature, "unknown"),
        ]
    )

    comparison = compare_evidence_fingerprints(first, second)

    assert comparison.first_only_present == 1
    assert comparison.second_only_present == 1
    assert comparison.unknown_involved == 1
    assert comparison.assessable_both == 2
    assert comparison.both_present == 0
    assert comparison.both_absent == 0


def test_incompatible_schemas_are_rejected():
    first = fingerprint([(atom_feature("hydrogen_bond"), "present")])
    second = build_evidence_fingerprint(
        [observation(atom_feature("hydrogen_bond"), "present")],
        schema_id="test-schema",
        criteria_fingerprint=criteria_fingerprint({"different": "criteria"}),
    )

    with pytest.raises(ValueError, match="incompatible"):
        compare_evidence_fingerprints(first, second)


def test_reference_2hyy_evidence_fingerprint_is_unassessable():
    path = (
        Path(__file__).resolve().parents[1]
        / "data"
        / "structures"
        / "2HYY.evidence-fingerprint.json"
    )
    result = load_interaction_evidence(path)

    assert validate_interaction_evidence(result) == []
    assert result.accession == "2HYY"
    assert len(result.schema.features) == 6
    assert all(feature.scope == "category" for feature in result.schema.features)
    assert result.fingerprint.present_count == 0
    assert result.fingerprint.absent_count == 0
    assert result.fingerprint.unknown_count == 6
    assert result.fingerprint.assessable_count == 0
    statuses = dict(result.category_statuses)
    assert statuses["hydrogen_bond"] == "unknown"
    assert statuses["metal_coordination"] == "not_applicable"


def test_interaction_evidence_record_round_trip(tmp_path):
    source = Path(__file__).resolve().parents[1] / "data" / "structures" / "2HYY.interactions.json"
    from moleculab import build_interaction_evidence, load_interactions

    interactions = load_interactions(source)
    record = build_interaction_evidence(
        interaction_record=interactions,
        source_manifest_path="manifest.json",
        source_manifest_sha256="a" * 64,
        ligand_identity_record_path="identity.json",
        ligand_identity_record_sha256="b" * 64,
        binding_site_record_path="site.json",
        binding_site_record_sha256="c" * 64,
        interaction_record_path="interactions.json",
        interaction_record_sha256="d" * 64,
    )
    path = write_interaction_evidence(tmp_path / "evidence.json", record)
    loaded = load_interaction_evidence(path)

    assert loaded == record
    assert loaded.fingerprint.schema == loaded.schema


def test_invalid_category_status_is_rejected():
    path = (
        Path(__file__).resolve().parents[1]
        / "data"
        / "structures"
        / "2HYY.evidence-fingerprint.json"
    )
    result = load_interaction_evidence(path)
    from dataclasses import replace

    bad = replace(
        result,
        category_statuses=(("hydrogen_bond", "assumed"),),
    )

    errors = validate_interaction_evidence(bad)

    assert any("invalid status" in error for error in errors)
    assert any("cover every interaction type" in error for error in errors)
