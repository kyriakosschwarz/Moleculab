# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Tests for dipole moments from charges and positions."""

import numpy as np
import pytest

from moleculab import DEBYE_PER_E_NM, DipoleMoment, dipole_moment


def test_simple_charge_separation_dipole():
    charges = np.array([1.0, -1.0])
    positions = np.array([[0.1, 0.0, 0.0], [-0.1, 0.0, 0.0]])

    result = dipole_moment(charges, positions)

    np.testing.assert_allclose(result.vector_e_nm, [0.2, 0.0, 0.0])
    assert result.magnitude_e_nm == pytest.approx(0.2)
    assert result.magnitude_debye == pytest.approx(0.2 * DEBYE_PER_E_NM)


def test_symmetric_charges_have_zero_dipole():
    charges = np.array([1.0, 1.0, -1.0, -1.0])
    positions = np.array(
        [
            [0.1, 0.0, 0.0],
            [-0.1, 0.0, 0.0],
            [0.0, 0.1, 0.0],
            [0.0, -0.1, 0.0],
        ]
    )

    result = dipole_moment(charges, positions)

    np.testing.assert_allclose(result.vector_e_nm, [0.0, 0.0, 0.0], atol=1e-15)
    assert result.magnitude_e_nm == pytest.approx(0.0, abs=1e-15)


def test_neutral_dipole_is_origin_independent():
    charges = np.array([1.0, -1.0])
    positions = np.array([[0.1, 0.0, 0.0], [-0.1, 0.0, 0.0]])
    shift = np.array([5.0, -3.0, 2.0])

    original = dipole_moment(charges, positions)
    shifted = dipole_moment(charges, positions + shift)

    np.testing.assert_allclose(shifted.vector_e_nm, original.vector_e_nm)


def test_charged_system_requires_explicit_origin():
    charges = np.array([1.0])
    positions = np.array([[0.1, 0.0, 0.0]])

    with pytest.raises(ValueError, match="explicit origin"):
        dipole_moment(charges, positions)

    result = dipole_moment(charges, positions, origin=[0.0, 0.0, 0.0])
    np.testing.assert_allclose(result.vector_e_nm, [0.1, 0.0, 0.0])
    assert result.magnitude_e_nm == pytest.approx(0.1)


def test_dipole_conversion_factor():
    assert DEBYE_PER_E_NM == pytest.approx(48.032, rel=1e-4)


def test_dipole_validation():
    with pytest.raises(ValueError, match="at least one"):
        dipole_moment(np.array([]), np.empty((0, 3)))
    with pytest.raises(ValueError, match="shape"):
        dipole_moment(np.array([1.0]), np.zeros((1, 2)))
    with pytest.raises(ValueError, match="same length"):
        dipole_moment(np.array([1.0, -1.0]), np.zeros((1, 3)))
    with pytest.raises(ValueError, match="finite"):
        dipole_moment(np.array([np.nan]), np.zeros((1, 3)))
    with pytest.raises(ValueError, match="origin"):
        dipole_moment(np.array([1.0, -1.0]), np.zeros((2, 3)), origin=[0.0, 0.0])
    with pytest.raises(ValueError, match="neutral_tolerance"):
        dipole_moment(np.array([1.0, -1.0]), np.zeros((2, 3)), neutral_tolerance=-1.0)


def test_dipole_vector_rotates_with_coordinates():
    charges = np.array([1.0, -1.0])
    positions = np.array([[0.1, 0.0, 0.0], [-0.1, 0.0, 0.0]])
    rotation = np.array(
        [
            [0.0, -1.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
        ]
    )

    original = dipole_moment(charges, positions)
    rotated = dipole_moment(charges, positions @ rotation.T)

    np.testing.assert_allclose(
        rotated.vector_e_nm,
        rotation @ original.vector_e_nm,
    )
    assert rotated.magnitude_e_nm == pytest.approx(original.magnitude_e_nm)


def test_dipole_vectors_convert_to_debye():
    charges = np.array([1.0, -1.0])
    positions = np.array([[0.1, 0.0, 0.0], [-0.1, 0.0, 0.0]])

    result = dipole_moment(charges, positions)

    assert isinstance(result, DipoleMoment)
    np.testing.assert_allclose(
        result.vector_debye,
        result.vector_e_nm * DEBYE_PER_E_NM,
    )
