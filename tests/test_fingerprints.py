# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Tests for interaction fingerprints and trajectory occupancy."""

import numpy as np
import pytest

from moleculab import (
    FingerprintSchema,
    HydrogenBond,
    HydrogenBondCriteria,
    InteractionFeature,
    InteractionFingerprint,
    build_fingerprint,
    calculate_occupancy,
)

CRITERIA = HydrogenBondCriteria(max_distance=0.35, min_angle_degrees=120.0)


def make_hbond(
    protein_residue: str = "ASP A 123",
    ligand_atom: str = "N7",
    ligand_residue: str = "LIG B 1",
    *,
    criteria: HydrogenBondCriteria = CRITERIA,
    protein_index: int = 0,
    ligand_index: int = 0,
) -> HydrogenBond:
    return HydrogenBond(
        protein_index=protein_index,
        ligand_index=ligand_index,
        donor_index=0,
        acceptor_index=0,
        donor_is_protein=True,
        protein_residue=protein_residue,
        ligand_residue=ligand_residue,
        protein_label="OD1",
        ligand_label=ligand_atom,
        donor_label="OD1",
        acceptor_label=ligand_atom,
        distance=0.29,
        angle_degrees=165.0,
        hydrogen_position=np.array([0.0, 0.0, 0.0]),
        hydrogen_index=0,
        criteria=criteria,
    )


def test_empty_interactions_produce_a_valid_empty_fingerprint():
    fingerprint = build_fingerprint([])

    assert len(fingerprint.schema) == 0
    assert fingerprint.bits.shape == (0,)
    assert fingerprint.counts.shape == (0,)
    assert fingerprint.criteria is None
    assert fingerprint.present == ()


def test_duplicate_observations_do_not_inflate_binary_features():
    bond = make_hbond()

    fingerprint = build_fingerprint([bond, bond])

    assert len(fingerprint.schema) == 1
    assert fingerprint.bits.tolist() == [True]
    assert fingerprint.counts.tolist() == [2]


def test_schema_is_deterministic_from_interaction_order():
    bonds = [
        make_hbond("PHE A 87", "O2"),
        make_hbond("ASP A 123", "N7"),
        make_hbond("ASP A 123", "N7"),
    ]

    first = FingerprintSchema.from_interactions(bonds)
    second = FingerprintSchema.from_interactions(list(reversed(bonds)))

    assert first == second
    assert [feature.key for feature in first.features] == sorted(
        feature.key for feature in first.features
    )


def test_atom_ordering_does_not_change_feature_identity():
    first = build_fingerprint(
        [make_hbond(protein_index=0, ligand_index=0)]
    )
    reordered = build_fingerprint(
        [make_hbond(protein_index=7, ligand_index=3)]
    )

    assert first.schema == reordered.schema
    assert first.bits.tolist() == reordered.bits.tolist()
    assert first.present == reordered.present


def test_occupancy_zero_one_and_intermediate_values():
    feature_a = make_hbond("ASP A 123", "N7")
    feature_b = make_hbond("PHE A 87", "O2")
    feature_c = make_hbond("LYS A 45", "N1")
    schema = FingerprintSchema.from_interactions([feature_a, feature_b, feature_c])

    frames = [
        build_fingerprint([feature_a, feature_b], schema),
        build_fingerprint([feature_a], schema),
        build_fingerprint([], schema),
    ]
    result = calculate_occupancy(frames)

    occupancy_by_feature = {
        feature: value
        for feature, value in zip(result.features, result.occupancy, strict=True)
    }
    counts_by_feature = {
        feature: value
        for feature, value in zip(result.features, result.counts, strict=True)
    }
    feature_a_id = InteractionFeature.from_hydrogen_bond(feature_a)
    feature_b_id = InteractionFeature.from_hydrogen_bond(feature_b)
    feature_c_id = InteractionFeature.from_hydrogen_bond(feature_c)

    assert occupancy_by_feature[feature_a_id] == pytest.approx(2.0 / 3.0)
    assert occupancy_by_feature[feature_b_id] == pytest.approx(1.0 / 3.0)
    assert occupancy_by_feature[feature_c_id] == pytest.approx(0.0)
    assert counts_by_feature[feature_a_id] == 2
    assert counts_by_feature[feature_b_id] == 1
    assert counts_by_feature[feature_c_id] == 0
    assert result.n_frames == 3
    assert result.n_valid_frames == 3
    assert result.n_missing_frames == 0


def test_missing_frames_are_not_treated_as_absent_interactions():
    bond = make_hbond()
    schema = FingerprintSchema.from_interactions([bond])
    frames = [
        build_fingerprint([bond], schema),
        None,
        build_fingerprint([], schema),
        None,
    ]

    result = calculate_occupancy(frames)

    assert result.occupancy.tolist() == [0.5]
    assert result.counts.tolist() == [1]
    assert result.n_frames == 4
    assert result.n_valid_frames == 2
    assert result.n_missing_frames == 2


def test_all_missing_frames_are_rejected():
    with pytest.raises(ValueError, match="at least one frame"):
        calculate_occupancy([None, None])


def test_empty_frame_list_is_rejected():
    with pytest.raises(ValueError, match="at least one frame"):
        calculate_occupancy([])


def test_residue_chain_and_insertion_code_remain_distinct():
    features = {
        InteractionFeature.from_hydrogen_bond(make_hbond("ASP A 123")),
        InteractionFeature.from_hydrogen_bond(make_hbond("ASP B 123")),
        InteractionFeature.from_hydrogen_bond(make_hbond("ASP A 123A")),
    }

    assert len(features) == 3


def test_results_retain_schema_and_criteria():
    bond = make_hbond()
    fingerprint = build_fingerprint([bond])
    result = calculate_occupancy([fingerprint])

    assert fingerprint.criteria == CRITERIA
    assert result.criteria == CRITERIA
    assert result.schema == fingerprint.schema
    assert result.features == fingerprint.features


def test_mixed_criteria_in_one_frame_are_rejected():
    other = HydrogenBondCriteria(max_distance=0.30, min_angle_degrees=130.0)

    with pytest.raises(ValueError, match="same criteria"):
        build_fingerprint([make_hbond(), make_hbond(criteria=other)])


def test_different_schemas_across_frames_are_rejected():
    first = build_fingerprint([make_hbond("ASP A 123")])
    second = build_fingerprint([make_hbond("PHE A 87", "O2")])

    with pytest.raises(ValueError, match="same fingerprint schema"):
        calculate_occupancy([first, second])


def test_missing_stable_residue_identifier_is_rejected():
    bond = make_hbond(protein_residue="")

    with pytest.raises(ValueError, match="protein_residue"):
        build_fingerprint([bond])


def test_fingerprint_validates_bits_and_counts():
    schema = FingerprintSchema.from_interactions([make_hbond()])

    with pytest.raises(ValueError, match="bits"):
        InteractionFingerprint(schema=schema, bits=[True, False], counts=[1])
    with pytest.raises(ValueError, match="counts"):
        InteractionFingerprint(schema=schema, bits=[True], counts=[-1])
    with pytest.raises(ValueError, match="counts > 0"):
        InteractionFingerprint(schema=schema, bits=[False], counts=[1])


def test_duplicate_counts_do_not_change_frame_occupancy():
    bond = make_hbond()
    schema = FingerprintSchema.from_interactions([bond])
    duplicated = build_fingerprint([bond, bond], schema)

    result = calculate_occupancy([duplicated])

    assert duplicated.counts.tolist() == [2]
    assert result.occupancy.tolist() == [1.0]
    assert result.counts.tolist() == [1]
