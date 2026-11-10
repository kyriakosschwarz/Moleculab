# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Tests for temperature calculation and the Berendsen thermostat."""

import numpy as np
import pytest

from moleculab import (
    Bond,
    HarmonicBondPotential,
    Particle,
    System,
    Topology,
    berendsen_thermostat,
    degrees_of_freedom,
    kinetic_temperature,
)


def two_particle_system(velocity_1=(0.1, 0.0, 0.0), velocity_2=(-0.1, 0.0, 0.0)):
    return System(
        particles=[
            Particle(1.0, [0.0, 0.0, 0.0], velocity_1),
            Particle(1.0, [1.0, 0.0, 0.0], velocity_2),
        ]
    )


def test_degrees_of_freedom_removes_center_of_mass_translation():
    system = two_particle_system()

    assert degrees_of_freedom(system) == 3
    assert degrees_of_freedom(system, remove_com=False) == 6


def test_single_particle_keeps_all_three_degrees_of_freedom():
    system = System(particles=[Particle(12.0, [0.0, 0.0, 0.0], [0.1, 0.0, 0.0])])

    assert degrees_of_freedom(system) == 3


def test_kinetic_temperature_matches_the_equipartition_relation():
    system = System(particles=[Particle(12.0, [0.0, 0.0, 0.0], [0.1, 0.0, 0.0])])

    # K = 0.5 * 12 * 0.1^2 = 0.06 kJ/mol
    # T = 2K / (3R)
    expected = 2.0 * 0.06 / (3.0 * 0.008314462618)
    assert kinetic_temperature(system) == pytest.approx(expected)


def test_thermostat_at_target_temperature_does_not_change_velocities():
    system = two_particle_system()
    target = kinetic_temperature(system)

    updated = berendsen_thermostat(system, target, coupling_time=0.1, dt=0.001)

    np.testing.assert_allclose(updated.velocities, system.velocities)


def test_thermostat_heats_a_cold_system():
    system = two_particle_system()
    target = 2.0 * kinetic_temperature(system)

    updated = berendsen_thermostat(system, target, coupling_time=0.1, dt=0.01)

    assert kinetic_temperature(updated) > kinetic_temperature(system)


def test_thermostat_cools_a_hot_system():
    system = two_particle_system()
    target = 0.5 * kinetic_temperature(system)

    updated = berendsen_thermostat(system, target, coupling_time=0.1, dt=0.01)

    assert kinetic_temperature(updated) < kinetic_temperature(system)


def test_zero_temperature_thermostat_stops_all_motion():
    system = two_particle_system()

    updated = berendsen_thermostat(system, 0.0, coupling_time=0.1, dt=0.01)

    np.testing.assert_allclose(updated.velocities, np.zeros((2, 3)))


def test_thermostat_preserves_positions_topology_and_potentials():
    topology = Topology(
        bonds=[Bond(0, 1, HarmonicBondPotential(stiffness=1000.0, length0=1.0))]
    )
    system = System(
        particles=two_particle_system().particles,
        topology=topology,
    )

    updated = berendsen_thermostat(system, 300.0, coupling_time=0.1, dt=0.001)

    np.testing.assert_allclose(updated.positions, system.positions)
    assert updated.topology is topology


def test_thermostat_validates_inputs():
    system = two_particle_system()

    with pytest.raises(ValueError, match="temperature"):
        berendsen_thermostat(system, -1.0, coupling_time=0.1, dt=0.001)
    with pytest.raises(ValueError, match="coupling_time"):
        berendsen_thermostat(system, 300.0, coupling_time=0.0, dt=0.001)
    with pytest.raises(ValueError, match="dt"):
        berendsen_thermostat(system, 300.0, coupling_time=0.1, dt=0.0)


def test_zero_kinetic_energy_cannot_be_heated_by_scaling():
    system = two_particle_system(velocity_1=(0.0, 0.0, 0.0), velocity_2=(0.0, 0.0, 0.0))

    with pytest.raises(ValueError, match="zero kinetic energy"):
        berendsen_thermostat(system, 300.0, coupling_time=0.1, dt=0.001)
