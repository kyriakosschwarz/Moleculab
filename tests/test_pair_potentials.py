"""Tests for Lennard-Jones and Coulomb pair potentials."""

import numpy as np
import pytest

from moleculab import CoulombPotential, LennardJonesPotential
from moleculab.units import COULOMB_CONSTANT


def numerical_force(potential, displacement, h=1e-6):
    """Central-difference negative gradient of a pair potential."""

    force = np.empty(3)
    for axis in range(3):
        plus = np.array(displacement, dtype=float)
        minus = np.array(displacement, dtype=float)
        plus[axis] += h
        minus[axis] -= h
        force[axis] = -(potential.energy(plus) - potential.energy(minus)) / (2.0 * h)
    return force


def test_lennard_jones_energy_is_zero_at_sigma():
    potential = LennardJonesPotential(epsilon=1.0, sigma=1.0)

    assert potential.energy([1.0, 0.0, 0.0]) == pytest.approx(0.0, abs=1e-12)


def test_lennard_jones_has_a_well_of_depth_epsilon():
    potential = LennardJonesPotential(epsilon=2.0, sigma=1.0)
    minimum_distance = 2.0 ** (1.0 / 6.0)

    assert potential.energy([minimum_distance, 0.0, 0.0]) == pytest.approx(-2.0)
    np.testing.assert_allclose(
        potential.force([minimum_distance, 0.0, 0.0]),
        [0.0, 0.0, 0.0],
        atol=1e-12,
    )


def test_lennard_jones_repels_at_short_range():
    potential = LennardJonesPotential(epsilon=1.0, sigma=1.0)

    force = potential.force([1.0, 0.0, 0.0])

    # F = 24*epsilon/sigma at r = sigma, pointing away from the other particle.
    np.testing.assert_allclose(force, [24.0, 0.0, 0.0])


def test_lennard_jones_force_is_negative_energy_gradient():
    potential = LennardJonesPotential(epsilon=0.997, sigma=0.34)
    displacement = np.array([0.45, -0.12, 0.08])

    np.testing.assert_allclose(
        potential.force(displacement),
        numerical_force(potential, displacement),
        atol=1e-6,
    )


def test_lennard_jones_is_symmetric_under_particle_swap():
    potential = LennardJonesPotential(epsilon=1.0, sigma=1.0)
    displacement = np.array([1.3, -0.2, 0.4])

    assert potential.energy(displacement) == pytest.approx(
        potential.energy(-displacement)
    )
    np.testing.assert_allclose(
        potential.force(-displacement),
        -potential.force(displacement),
    )


def test_lennard_jones_parameters_must_be_positive():
    with pytest.raises(ValueError, match="epsilon"):
        LennardJonesPotential(epsilon=0.0, sigma=1.0)
    with pytest.raises(ValueError, match="sigma"):
        LennardJonesPotential(epsilon=1.0, sigma=-0.34)


def test_pair_potentials_reject_zero_separation():
    lj = LennardJonesPotential(epsilon=1.0, sigma=1.0)
    coulomb = CoulombPotential(charge_i=1.0, charge_j=1.0)

    with pytest.raises(ValueError, match="displacement"):
        lj.energy([0.0, 0.0, 0.0])
    with pytest.raises(ValueError, match="displacement"):
        coulomb.energy([0.0, 0.0, 0.0])


def test_coulomb_energy_of_two_unit_charges_at_one_nanometre():
    potential = CoulombPotential(charge_i=1.0, charge_j=1.0)

    assert potential.energy([1.0, 0.0, 0.0]) == pytest.approx(COULOMB_CONSTANT)


def test_coulomb_like_charges_repel():
    potential = CoulombPotential(charge_i=1.0, charge_j=1.0)

    np.testing.assert_allclose(
        potential.force([1.0, 0.0, 0.0]),
        [COULOMB_CONSTANT, 0.0, 0.0],
    )


def test_coulomb_opposite_charges_attract():
    potential = CoulombPotential(charge_i=1.0, charge_j=-1.0)

    assert potential.energy([1.0, 0.0, 0.0]) == pytest.approx(-COULOMB_CONSTANT)
    np.testing.assert_allclose(
        potential.force([1.0, 0.0, 0.0]),
        [-COULOMB_CONSTANT, 0.0, 0.0],
    )


def test_coulomb_force_is_negative_energy_gradient():
    potential = CoulombPotential(charge_i=1.5, charge_j=-0.75)
    displacement = np.array([0.7, 0.2, -0.1])

    np.testing.assert_allclose(
        potential.force(displacement),
        numerical_force(potential, displacement),
        atol=1e-6,
    )


def test_coulomb_charges_must_be_finite():
    with pytest.raises(ValueError, match="charge_i"):
        CoulombPotential(charge_i=float("nan"), charge_j=1.0)
    with pytest.raises(ValueError, match="charge_j"):
        CoulombPotential(charge_i=1.0, charge_j=float("inf"))
