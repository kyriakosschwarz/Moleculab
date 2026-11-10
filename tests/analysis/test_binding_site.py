# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Tests for parser-independent binding-site geometry."""

import numpy as np
import pytest

from moleculab import (
    BindingSiteGeometry,
    binding_site_geometry,
    group_indices_by_residue,
    minimum_distances,
    pocket_geometry,
    select_pocket_atoms,
)


def test_minimum_distances_use_the_closest_ligand_atom():
    protein = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 0.0]])
    ligand = np.array([[3.0, 4.0, 0.0], [0.0, 2.0, 0.0]])

    distances = minimum_distances(protein, ligand)

    np.testing.assert_allclose(distances, [2.0, np.sqrt(2.0)])


def test_minimum_distances_handle_an_empty_protein():
    distances = minimum_distances(
        np.empty((0, 3)),
        np.array([[0.0, 0.0, 0.0]]),
    )

    assert distances.shape == (0,)


def test_minimum_distances_reject_an_empty_ligand():
    with pytest.raises(ValueError, match="at least one atom"):
        minimum_distances(
            np.array([[0.0, 0.0, 0.0]]),
            np.empty((0, 3)),
        )


def test_minimum_distances_validate_shapes_and_values():
    with pytest.raises(ValueError, match="protein_positions"):
        minimum_distances(np.zeros((2, 2)), np.zeros((1, 3)))
    with pytest.raises(ValueError, match="ligand_positions"):
        minimum_distances(np.zeros((2, 3)), np.zeros((1, 2)))
    with pytest.raises(ValueError, match="finite"):
        minimum_distances(
            np.array([[0.0, 0.0, np.nan]]),
            np.zeros((1, 3)),
        )


def test_pocket_selection_is_strict_by_default():
    protein = np.array([[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
    ligand = np.array([[0.0, 0.0, 0.0]])

    strict = select_pocket_atoms(protein, ligand, 1.0)
    inclusive = select_pocket_atoms(protein, ligand, 1.0, inclusive=True)

    assert strict.tolist() == []
    assert inclusive.tolist() == [0]


def test_pocket_selection_returns_sorted_unique_indices():
    protein = np.array(
        [
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.5, 0.0, 0.0],
        ]
    )
    ligand = np.array([[0.0, 0.0, 0.0]])

    indices = select_pocket_atoms(protein, ligand, 1.5)

    assert indices.tolist() == [0, 2]
    assert np.unique(indices).size == indices.size


def test_pocket_selection_validates_the_cutoff():
    protein = np.array([[1.0, 0.0, 0.0]])
    ligand = np.array([[0.0, 0.0, 0.0]])

    for cutoff in [0.0, -1.0, np.nan, np.inf]:
        with pytest.raises(ValueError, match="cutoff"):
            select_pocket_atoms(protein, ligand, cutoff)


def test_distances_are_translation_invariant():
    protein = np.array([[0.1, 0.2, 0.3], [0.4, 0.0, 0.1]])
    ligand = np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]])
    shift = np.array([10.0, -5.0, 2.0])

    original = minimum_distances(protein, ligand)
    shifted = minimum_distances(protein + shift, ligand + shift)

    np.testing.assert_allclose(original, shifted)
    assert select_pocket_atoms(protein, ligand, 0.5).tolist() == (
        select_pocket_atoms(protein + shift, ligand + shift, 0.5).tolist()
    )


def test_distances_are_rotation_invariant():
    protein = np.array([[0.1, 0.2, 0.3], [0.4, 0.0, 0.1]])
    ligand = np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]])
    rotation = np.array(
        [
            [0.0, -1.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
        ]
    )

    original = minimum_distances(protein, ligand)
    rotated = minimum_distances(protein @ rotation.T, ligand @ rotation.T)

    np.testing.assert_allclose(original, rotated)
    assert select_pocket_atoms(protein, ligand, 0.5).tolist() == (
        select_pocket_atoms(protein @ rotation.T, ligand @ rotation.T, 0.5).tolist()
    )


def test_pocket_geometry_reports_centre_and_bounding_box():
    protein = np.array(
        [
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 4.0, 0.0],
            [10.0, 10.0, 10.0],
        ]
    )

    geometry = pocket_geometry(protein, [0, 1, 2], cutoff=5.0)

    assert isinstance(geometry, BindingSiteGeometry)
    assert geometry.n_atoms == 3
    np.testing.assert_allclose(geometry.center, [2.0 / 3.0, 4.0 / 3.0, 0.0])
    np.testing.assert_allclose(geometry.lower, [0.0, 0.0, 0.0])
    np.testing.assert_allclose(geometry.upper, [2.0, 4.0, 0.0])
    np.testing.assert_allclose(geometry.dimensions, [2.0, 4.0, 0.0])
    assert geometry.cutoff == 5.0


def test_pocket_geometry_rejects_an_empty_selection():
    with pytest.raises(ValueError, match="at least one atom"):
        pocket_geometry(np.zeros((3, 3)), [])


def test_pocket_geometry_rejects_invalid_indices():
    protein = np.zeros((3, 3))

    with pytest.raises(ValueError, match="out-of-range"):
        pocket_geometry(protein, [0, 3])
    with pytest.raises(ValueError, match="unique"):
        pocket_geometry(protein, [1, 1])
    with pytest.raises(ValueError, match="integers"):
        pocket_geometry(protein, [0.5])


def test_residue_grouping_rejects_non_integer_indices():
    with pytest.raises(ValueError, match="integers"):
        group_indices_by_residue([0.5], ["A"])
    with pytest.raises(ValueError, match="unique"):
        group_indices_by_residue([0, 0], ["A", "A"])


def test_binding_site_geometry_composes_selection_and_geometry():
    protein = np.array(
        [
            [0.1, 0.0, 0.0],
            [0.0, 0.2, 0.0],
            [0.0, 0.0, 0.3],
            [2.0, 0.0, 0.0],
        ]
    )
    ligand = np.array([[0.0, 0.0, 0.0]])

    geometry = binding_site_geometry(protein, ligand, cutoff=0.5)

    assert geometry.atom_indices.tolist() == [0, 1, 2]
    assert geometry.n_atoms == 3
    np.testing.assert_allclose(geometry.center, [0.1 / 3.0, 0.2 / 3.0, 0.1])


def test_binding_site_geometry_reports_an_empty_selection():
    protein = np.array([[1.0, 0.0, 0.0]])
    ligand = np.array([[0.0, 0.0, 0.0]])

    with pytest.raises(ValueError, match="no protein atoms"):
        binding_site_geometry(protein, ligand, cutoff=0.5)


def test_residue_grouping_keeps_order_and_removes_duplicate_residues():
    grouped = group_indices_by_residue(
        [0, 3, 5, 7],
        ["A", "B", "A", "C"],
    )

    assert grouped == {"A": (0, 5), "B": (3,), "C": (7,)}


def test_residue_grouping_validates_lengths():
    with pytest.raises(ValueError, match="same length"):
        group_indices_by_residue([0, 1], ["A"])


def test_no_periodic_minimum_image_is_applied():
    protein = np.array([[0.95, 0.0, 0.0]])
    ligand = np.array([[0.05, 0.0, 0.0]])

    assert minimum_distances(protein, ligand)[0] == pytest.approx(0.9)
    assert select_pocket_atoms(protein, ligand, 0.5).tolist() == []


def test_distance_scale_changes_by_the_coordinate_scale_factor():
    protein_nm = np.array([[0.1, 0.0, 0.0], [0.2, 0.0, 0.0]])
    ligand_nm = np.array([[0.0, 0.0, 0.0]])

    distances_nm = minimum_distances(protein_nm, ligand_nm)
    distances_angstrom = minimum_distances(protein_nm * 10.0, ligand_nm * 10.0)

    np.testing.assert_allclose(distances_angstrom, distances_nm * 10.0)
