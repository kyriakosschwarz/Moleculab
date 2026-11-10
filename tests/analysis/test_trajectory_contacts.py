# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Tests for geometric trajectory contacts, occupancy, and lifetimes."""

import numpy as np
import pytest

from moleculab import (
    ContactLifetimeResult,
    ContactOccupancyResult,
    ResidueContactMatrix,
    aggregate_contacts_by_residue,
    calculate_contact_lifetimes,
    calculate_contact_matrix,
    calculate_contact_occupancy,
)


def test_contact_matrix_uses_minimum_pairwise_distances():
    protein = np.array([[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]])
    ligand = np.array([[0.2, 0.0, 0.0]])

    contacts = calculate_contact_matrix(protein, ligand, cutoff=0.3)

    assert contacts.tolist() == [[True], [False]]


def test_contact_cutoff_is_strict_by_default():
    protein = np.array([[0.3, 0.0, 0.0]])
    ligand = np.array([[0.0, 0.0, 0.0]])

    assert calculate_contact_matrix(protein, ligand, 0.3).tolist() == [[False]]
    assert calculate_contact_matrix(
        protein,
        ligand,
        0.3,
        inclusive=True,
    ).tolist() == [[True]]


def test_contact_matrix_is_translation_and_rotation_invariant():
    protein = np.array([[0.0, 0.0, 0.0], [0.4, 0.0, 0.0]])
    ligand = np.array([[0.2, 0.0, 0.0], [0.6, 0.0, 0.0]])
    rotation = np.array(
        [
            [0.0, -1.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    shift = np.array([1.0, -2.0, 3.0])

    original = calculate_contact_matrix(protein, ligand, 0.25)
    moved = calculate_contact_matrix(
        protein @ rotation.T + shift,
        ligand @ rotation.T + shift,
        0.25,
    )

    np.testing.assert_array_equal(original, moved)


def test_contact_matrix_validation():
    with pytest.raises(ValueError, match="protein_positions"):
        calculate_contact_matrix(np.zeros((2, 2)), np.zeros((1, 3)), 0.3)
    with pytest.raises(ValueError, match="ligand_positions"):
        calculate_contact_matrix(np.zeros((2, 3)), np.zeros((1, 2)), 0.3)
    with pytest.raises(ValueError, match="cutoff"):
        calculate_contact_matrix(np.zeros((2, 3)), np.zeros((1, 3)), 0.0)
    with pytest.raises(ValueError, match="at least one atom"):
        calculate_contact_matrix(np.zeros((2, 3)), np.empty((0, 3)), 0.3)


def test_residue_aggregation_uses_any_contacting_atom():
    contact_matrix = np.array(
        [
            [True, False],
            [False, False],
            [False, True],
            [False, False],
        ]
    )
    labels = ["ASP A 123", "ASP A 123", "PHE A 87", "PHE A 87"]

    result = aggregate_contacts_by_residue(contact_matrix, labels)

    assert isinstance(result, ResidueContactMatrix)
    assert result.residue_ids == ("ASP A 123", "PHE A 87")
    assert result.contacts.tolist() == [
        [True, False],
        [False, True],
    ]


def test_residue_aggregation_validation():
    with pytest.raises(ValueError, match="two-dimensional"):
        aggregate_contacts_by_residue(np.zeros((2,)), ["A", "B"])
    with pytest.raises(ValueError, match="one residue ID"):
        aggregate_contacts_by_residue(np.zeros((2, 1)), ["A"])


def test_contact_occupancy_uses_valid_frames_only():
    frame_present = np.array([[True, False]])
    frame_absent = np.array([[False, False]])
    frames = [frame_present, None, frame_absent, frame_present]

    result = calculate_contact_occupancy(frames)

    assert isinstance(result, ContactOccupancyResult)
    np.testing.assert_allclose(result.occupancy, [[2.0 / 3.0, 0.0]])
    assert result.counts.tolist() == [[2, 0]]
    assert result.n_frames == 4
    assert result.n_valid_frames == 3
    assert result.n_missing_frames == 1


def test_contact_occupancy_validation():
    with pytest.raises(ValueError, match="at least one frame"):
        calculate_contact_occupancy([])
    with pytest.raises(ValueError, match="at least one contact frame"):
        calculate_contact_occupancy([None, None])
    with pytest.raises(ValueError, match="two-dimensional"):
        calculate_contact_occupancy([np.zeros((2,))])
    with pytest.raises(ValueError, match="same shape"):
        calculate_contact_occupancy([np.zeros((1, 2)), np.zeros((2, 1))])


def test_contact_lifetime_run_lengths():
    presence = np.array([True, True, False, True, True, True, False, False, True])

    result = calculate_contact_lifetimes(presence, frame_spacing=0.1)

    assert isinstance(result, ContactLifetimeResult)
    assert result.run_lengths_frames == (2, 3, 1)
    assert result.n_runs == 3
    assert result.mean_lifetime_frames == pytest.approx(2.0)
    assert result.max_lifetime_frames == 3
    assert result.lifetimes == pytest.approx((0.2, 0.3, 0.1))
    assert result.mean_lifetime == pytest.approx(0.2)
    assert result.observed_fraction == pytest.approx(6.0 / 9.0)


def test_same_occupancy_can_have_very_different_lifetimes():
    continuous = np.array([True] * 6 + [False] * 6)
    intermittent = np.array([True, False] * 6, dtype=bool)

    occupancy = calculate_contact_occupancy(
        [
            np.array([[continuous[index], intermittent[index]]], dtype=bool)
            for index in range(len(continuous))
        ]
    )
    continuous_lifetime = calculate_contact_lifetimes(continuous)
    intermittent_lifetime = calculate_contact_lifetimes(intermittent)

    assert np.mean(continuous) == np.mean(intermittent) == 0.5
    np.testing.assert_allclose(occupancy.occupancy, [[0.5, 0.5]])
    assert continuous_lifetime.run_lengths_frames == (6,)
    assert intermittent_lifetime.run_lengths_frames == (1, 1, 1, 1, 1, 1)
    assert continuous_lifetime.mean_lifetime_frames == 6.0
    assert intermittent_lifetime.mean_lifetime_frames == 1.0


def test_contact_lifetime_all_false_and_all_true():
    no_contact = calculate_contact_lifetimes(np.zeros(5, dtype=bool))
    always_contact = calculate_contact_lifetimes(np.ones(5, dtype=bool))

    assert no_contact.run_lengths_frames == ()
    assert no_contact.mean_lifetime_frames == 0.0
    assert no_contact.max_lifetime_frames == 0
    assert no_contact.lifetimes is None
    assert no_contact.mean_lifetime is None
    assert always_contact.run_lengths_frames == (5,)
    assert always_contact.mean_lifetime_frames == 5.0


def test_missing_frames_break_observed_contact_runs():
    presence = [True, True, None, True, False]

    lifetime = calculate_contact_lifetimes(presence, frame_spacing=0.1)
    frames = [
        None if value is None else np.array([[value]], dtype=bool)
        for value in presence
    ]
    occupancy = calculate_contact_occupancy(frames)

    assert lifetime.run_lengths_frames == (2, 1)
    assert lifetime.n_frames == 5
    assert lifetime.n_valid_frames == 4
    assert lifetime.n_missing_frames == 1
    assert lifetime.observed_fraction == pytest.approx(0.75)
    assert lifetime.mean_lifetime_frames == pytest.approx(1.5)
    np.testing.assert_allclose(occupancy.occupancy, [[0.75]])
    assert occupancy.n_missing_frames == 1


def test_all_missing_lifetime_frames_are_rejected():
    with pytest.raises(ValueError, match="at least one contact frame"):
        calculate_contact_lifetimes([None, None])


def test_contact_lifetime_validation():
    with pytest.raises(ValueError, match="at least one frame"):
        calculate_contact_lifetimes(np.array([], dtype=bool))
    with pytest.raises(ValueError, match="frame_spacing"):
        calculate_contact_lifetimes(np.array([True]), frame_spacing=0.0)
