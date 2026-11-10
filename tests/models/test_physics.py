# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Tests for the general physics equations."""

import numpy as np
import pytest

from moleculab import (
    HarmonicPotential,
    Particle,
    acceleration,
    force_from_potential,
    kinetic_energy,
    potential_energy,
    total_energy,
)


def make_particle():
    return Particle(
        mass=12.0,
        position=[0.1, 0.0, 0.0],
        velocity=[0.1, 0.0, 0.0],
    )


def make_potential():
    return HarmonicPotential(stiffness=100.0, center=[0.0, 0.0, 0.0])


def test_kinetic_energy_matches_half_m_v_squared():
    particle = make_particle()

    # 0.5 * 12 * 0.1**2 = 0.06 kJ/mol
    assert kinetic_energy(particle) == pytest.approx(0.06)


def test_potential_energy_delegates_to_the_potential():
    particle = make_particle()

    assert potential_energy(particle, make_potential()) == pytest.approx(0.5)


def test_total_energy_is_kinetic_plus_potential():
    particle = make_particle()

    # 0.06 + 0.5 = 0.56 kJ/mol
    assert total_energy(particle, make_potential()) == pytest.approx(0.56)


def test_acceleration_follows_newtons_second_law():
    particle = make_particle()

    # a = F / m = (-10, 0, 0) / 12
    np.testing.assert_allclose(
        acceleration(particle, [-10.0, 0.0, 0.0]),
        [-10.0 / 12.0, 0.0, 0.0],
    )


def test_force_from_potential_delegates_to_the_potential():
    particle = make_particle()

    np.testing.assert_allclose(
        force_from_potential(particle, make_potential()),
        [-10.0, 0.0, 0.0],
    )
