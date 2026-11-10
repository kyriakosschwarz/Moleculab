# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Tests for multi-particle velocity Verlet integration."""

import numpy as np
import pytest

from moleculab import (
    LennardJonesPotential,
    Particle,
    System,
    simulate_velocity_verlet_system,
    velocity_verlet_step_system,
)


def free_system():
    return System(
        particles=[
            Particle(1.0, [0.0, 0.0, 0.0], [0.1, 0.0, 0.0]),
            Particle(2.0, [1.0, 0.0, 0.0], [-0.05, 0.0, 0.0]),
        ]
    )


def argon_pair_system():
    sigma = 0.34
    epsilon = 0.997
    equilibrium = 2.0 ** (1.0 / 6.0) * sigma
    return System(
        particles=[
            Particle(39.948, [-equilibrium / 2.0, 0.0, 0.0], [0.0, 0.01, 0.0]),
            Particle(39.948, [equilibrium / 2.0, 0.0, 0.0], [0.0, -0.01, 0.0]),
        ],
        pair_potentials=[LennardJonesPotential(epsilon=epsilon, sigma=sigma)],
    )


def test_system_step_does_not_mutate_input():
    system = free_system()
    positions_before = system.positions.copy()
    velocities_before = system.velocities.copy()

    velocity_verlet_step_system(system, dt=0.001)

    np.testing.assert_allclose(system.positions, positions_before)
    np.testing.assert_allclose(system.velocities, velocities_before)


def test_free_particles_move_with_constant_velocity():
    system = free_system()
    dt = 0.01

    trajectory = simulate_velocity_verlet_system(system, dt=dt, steps=10)

    np.testing.assert_allclose(
        trajectory.position[-1],
        system.positions + system.velocities * dt * 10,
        atol=1e-12,
    )
    np.testing.assert_allclose(
        trajectory.velocity[-1],
        system.velocities,
        atol=1e-12,
    )


def test_system_trajectory_shapes():
    system = argon_pair_system()

    trajectory = simulate_velocity_verlet_system(system, dt=0.001, steps=5)

    assert len(trajectory) == 6
    assert trajectory.time.shape == (6,)
    assert trajectory.position.shape == (6, 2, 3)
    assert trajectory.velocity.shape == (6, 2, 3)
    assert trajectory.total_energy.shape == (6,)


def test_system_trajectory_starts_at_the_initial_state():
    system = argon_pair_system()

    trajectory = simulate_velocity_verlet_system(system, dt=0.001, steps=3)

    np.testing.assert_allclose(trajectory.position[0], system.positions)
    np.testing.assert_allclose(trajectory.velocity[0], system.velocities)
    assert trajectory.total_energy[0] == pytest.approx(system.total_energy())


def test_velocity_verlet_conserves_energy_for_a_lennard_jones_pair():
    system = argon_pair_system()

    trajectory = simulate_velocity_verlet_system(system, dt=0.001, steps=1000)

    initial_energy = trajectory.total_energy[0]
    maximum_drift = np.max(np.abs(trajectory.total_energy - initial_energy))
    relative_drift = maximum_drift / abs(initial_energy)
    assert relative_drift < 1e-5


def test_pair_forces_conserve_total_momentum():
    system = argon_pair_system()

    trajectory = simulate_velocity_verlet_system(system, dt=0.001, steps=200)

    masses = np.array([particle.mass for particle in system.particles])[:, np.newaxis]
    initial_momentum = np.sum(masses * trajectory.velocity[0], axis=0)
    final_momentum = np.sum(masses * trajectory.velocity[-1], axis=0)
    np.testing.assert_allclose(final_momentum, initial_momentum, atol=1e-12)


def test_system_integrator_validates_timestep_and_steps():
    system = free_system()

    with pytest.raises(ValueError, match="dt"):
        simulate_velocity_verlet_system(system, dt=0.0, steps=1)
    with pytest.raises(ValueError, match="steps"):
        simulate_velocity_verlet_system(system, dt=0.001, steps=-1)
