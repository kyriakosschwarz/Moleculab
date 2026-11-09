"""Tests for potential-energy models.

The important test here is not merely that ``force`` returns an array.  It is
that the returned force is the negative derivative of the same energy
function.  A potential whose energy and force disagree can run without an
exception while still producing physically wrong dynamics.
"""

import numpy as np
import pytest

from moleculab import HarmonicPotential


def test_energy_is_zero_at_equilibrium():
    potential = HarmonicPotential(stiffness=100.0, center=[0.0, 0.0, 0.0])

    assert potential.energy([0.0, 0.0, 0.0]) == pytest.approx(0.0)


def test_energy_matches_half_k_distance_squared():
    potential = HarmonicPotential(stiffness=100.0, center=[0.0, 0.0, 0.0])

    # 0.5 * 100 * 0.1**2 = 0.5 kJ/mol
    assert potential.energy([0.1, 0.0, 0.0]) == pytest.approx(0.5)


def test_force_points_toward_equilibrium():
    potential = HarmonicPotential(stiffness=100.0, center=[0.0, 0.0, 0.0])

    np.testing.assert_allclose(potential.force([0.1, 0.0, 0.0]), [-10.0, 0.0, 0.0])


def test_force_uses_displacement_from_center():
    potential = HarmonicPotential(stiffness=50.0, center=[0.2, 0.0, 0.0])

    # F = -k (r - r0) = -50 * (0.4 - 0.2, 0, 0)
    np.testing.assert_allclose(potential.force([0.4, 0.0, 0.0]), [-10.0, 0.0, 0.0])


def test_force_is_negative_energy_gradient():
    potential = HarmonicPotential(stiffness=100.0, center=[0.0, 0.0, 0.0])
    position = np.array([0.1, -0.05, 0.02])
    h = 1e-6

    numerical_force = np.empty(3)
    for axis in range(3):
        plus = position.copy()
        minus = position.copy()
        plus[axis] += h
        minus[axis] -= h
        numerical_force[axis] = -(potential.energy(plus) - potential.energy(minus)) / (
            2.0 * h
        )

    np.testing.assert_allclose(
        potential.force(position),
        numerical_force,
        atol=1e-6,
    )


@pytest.mark.parametrize("stiffness", [0.0, -100.0, float("nan")])
def test_stiffness_must_be_finite_and_positive(stiffness):
    with pytest.raises(ValueError, match="stiffness"):
        HarmonicPotential(stiffness=stiffness)


def test_center_is_copied():
    center = np.array([0.0, 0.0, 0.0])
    potential = HarmonicPotential(stiffness=10.0, center=center)

    center[0] = 5.0

    np.testing.assert_allclose(potential.center, [0.0, 0.0, 0.0])
