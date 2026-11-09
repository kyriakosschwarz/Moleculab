"""Tests for bond and angle potential models."""

import numpy as np
import pytest

from moleculab import HarmonicAnglePotential, HarmonicBondPotential


def numerical_gradient(energy_function, position, h=1e-6):
    """Central-difference gradient of a scalar function of a 3-vector."""

    gradient = np.empty(3)
    for axis in range(3):
        plus = np.array(position, dtype=float)
        minus = np.array(position, dtype=float)
        plus[axis] += h
        minus[axis] -= h
        gradient[axis] = (energy_function(plus) - energy_function(minus)) / (2.0 * h)
    return gradient


def test_harmonic_bond_energy_is_zero_at_equilibrium():
    potential = HarmonicBondPotential(stiffness=1000.0, length0=0.15)

    assert potential.energy([0.15, 0.0, 0.0]) == pytest.approx(0.0)


def test_harmonic_bond_energy_is_half_k_stretch_squared():
    potential = HarmonicBondPotential(stiffness=1000.0, length0=0.15)

    displacement = [0.15 + 0.01, 0.0, 0.0]
    expected = 0.5 * 1000.0 * 0.01**2
    assert potential.energy(displacement) == pytest.approx(expected)


def test_stretched_bond_pulls_atom_i_toward_atom_j():
    potential = HarmonicBondPotential(stiffness=1000.0, length0=0.1)

    # displacement = r_i - r_j points from j to i (+x), so a stretched bond
    # must pull i in -x.
    force = potential.force([0.2, 0.0, 0.0])

    np.testing.assert_allclose(force, [-100.0, 0.0, 0.0])


def test_compressed_bond_pushes_atom_i_away_from_atom_j():
    potential = HarmonicBondPotential(stiffness=1000.0, length0=0.1)

    force = potential.force([0.05, 0.0, 0.0])

    assert force[0] > 0.0


def test_bond_force_is_negative_energy_gradient():
    potential = HarmonicBondPotential(stiffness=1000.0, length0=0.15)
    displacement = np.array([0.18, -0.02, 0.01])

    numerical_force = -numerical_gradient(potential.energy, displacement)

    np.testing.assert_allclose(potential.force(displacement), numerical_force, atol=1e-6)


def test_bond_parameters_must_be_positive():
    with pytest.raises(ValueError, match="stiffness"):
        HarmonicBondPotential(stiffness=0.0, length0=0.15)
    with pytest.raises(ValueError, match="length0"):
        HarmonicBondPotential(stiffness=1000.0, length0=0.0)


def test_angle_energy_is_zero_at_equilibrium():
    potential = HarmonicAnglePotential(stiffness=100.0, angle0=np.pi / 2.0)
    r_ji = np.array([1.0, 0.0, 0.0])
    r_jk = np.array([0.0, 1.0, 0.0])

    assert potential.energy(r_ji, r_jk) == pytest.approx(0.0)
    for force in potential.forces(r_ji, r_jk):
        np.testing.assert_allclose(force, np.zeros(3), atol=1e-12)


def test_angle_energy_matches_harmonic_formula():
    potential = HarmonicAnglePotential(stiffness=100.0, angle0=np.pi / 3.0)
    r_ji = np.array([1.0, 0.0, 0.0])
    r_jk = np.array([0.0, 1.0, 0.0])

    expected = 0.5 * 100.0 * (np.pi / 2.0 - np.pi / 3.0) ** 2
    assert potential.energy(r_ji, r_jk) == pytest.approx(expected)


def test_angle_forces_match_numerical_position_gradient():
    potential = HarmonicAnglePotential(stiffness=100.0, angle0=np.pi / 3.0)
    positions = np.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ]
    )

    def energy_at(coords):
        r_ji = coords[0] - coords[1]
        r_jk = coords[2] - coords[1]
        return potential.energy(r_ji, r_jk)

    numerical_forces = np.empty_like(positions)
    h = 1e-6
    for atom in range(3):
        for axis in range(3):
            plus = positions.copy()
            minus = positions.copy()
            plus[atom, axis] += h
            minus[atom, axis] -= h
            numerical_forces[atom, axis] = -(
                energy_at(plus) - energy_at(minus)
            ) / (2.0 * h)

    r_ji = positions[0] - positions[1]
    r_jk = positions[2] - positions[1]
    analytical_forces = np.array(potential.forces(r_ji, r_jk))

    np.testing.assert_allclose(analytical_forces, numerical_forces, atol=1e-6)
    np.testing.assert_allclose(analytical_forces.sum(axis=0), np.zeros(3), atol=1e-12)


def test_angle_force_is_undefined_for_linear_geometry():
    potential = HarmonicAnglePotential(stiffness=100.0, angle0=np.pi / 2.0)

    with pytest.raises(ValueError, match="linear"):
        potential.forces([1.0, 0.0, 0.0], [-1.0, 0.0, 0.0])


def test_angle0_must_be_between_zero_and_pi():
    with pytest.raises(ValueError, match="angle0"):
        HarmonicAnglePotential(stiffness=100.0, angle0=-0.1)
    with pytest.raises(ValueError, match="angle0"):
        HarmonicAnglePotential(stiffness=100.0, angle0=np.pi + 0.1)
