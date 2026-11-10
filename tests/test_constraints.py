# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Tests for SHAKE/RATTLE constraints and rigid-water integration."""

import numpy as np
import pytest

from moleculab import (
    DistanceConstraint,
    Particle,
    System,
    assert_constraints_satisfied,
    berendsen_barostat,
    build_water_box,
    constraint_error,
    minimize_energy,
    pressure_bar,
    rattle_velocities,
    shake_positions,
    simulate_velocity_verlet_system,
    velocity_verlet_step_system,
)


def test_distance_constraint_validation():
    with pytest.raises(ValueError, match="different"):
        DistanceConstraint(0, 0, 1.0)
    with pytest.raises(ValueError, match="non-negative"):
        DistanceConstraint(-1, 0, 1.0)
    with pytest.raises(ValueError, match="distance"):
        DistanceConstraint(0, 1, 0.0)


def test_shake_corrects_a_single_stretched_bond():
    positions = np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
    masses = np.array([1.0, 1.0])
    constraints = [DistanceConstraint(0, 1, 1.0)]

    result = shake_positions(positions, masses, constraints)

    np.testing.assert_allclose(
        result.positions,
        [[0.5, 0.0, 0.0], [1.5, 0.0, 0.0]],
    )
    assert result.max_error == pytest.approx(0.0, abs=1e-12)
    np.testing.assert_allclose(
        np.mean(result.positions, axis=0),
        np.mean(positions, axis=0),
    )


def test_shake_preserves_the_mass_weighted_center_of_mass():
    positions = np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
    masses = np.array([1.0, 3.0])
    constraints = [DistanceConstraint(0, 1, 1.0)]

    result = shake_positions(positions, masses, constraints)

    before = np.sum(masses[:, None] * positions, axis=0) / np.sum(masses)
    after = np.sum(masses[:, None] * result.positions, axis=0) / np.sum(masses)
    np.testing.assert_allclose(after, before, atol=1e-12)
    assert constraint_error(result.positions, constraints) == pytest.approx(
        0.0, abs=1e-10
    )


def test_rattle_removes_radial_relative_velocity():
    positions = np.array([[0.5, 0.0, 0.0], [1.5, 0.0, 0.0]])
    velocities = np.array([[0.1, 0.2, 0.0], [-0.1, 0.2, 0.0]])
    masses = np.array([1.0, 1.0])
    constraints = [DistanceConstraint(0, 1, 1.0)]

    result = rattle_velocities(positions, velocities, masses, constraints)

    displacement = positions[0] - positions[1]
    relative_velocity = result.velocities[0] - result.velocities[1]
    assert np.dot(displacement, relative_velocity) == pytest.approx(0.0, abs=1e-12)
    np.testing.assert_allclose(
        np.sum(result.velocities, axis=0),
        np.sum(velocities, axis=0),
    )


def test_water_constraints_correct_a_distorted_geometry():
    system = build_water_box(1, box_length=2.0, seed=1)
    positions = system.positions.copy()
    positions[1] += np.array([0.02, 0.01, 0.0])
    masses = np.array([particle.mass for particle in system.particles])

    result = shake_positions(positions, masses, system.constraints)

    assert len(system.constraints) == 3
    assert constraint_error(result.positions, system.constraints) < 1e-9


def test_water_rattle_makes_velocities_constraint_consistent():
    system = build_water_box(1, box_length=2.0, temperature=300.0, seed=2)
    masses = np.array([particle.mass for particle in system.particles])
    velocities = system.velocities.copy()
    velocities[1] += np.array([0.05, 0.0, 0.0])

    result = rattle_velocities(
        system.positions,
        velocities,
        masses,
        system.constraints,
    )

    for constraint in system.constraints:
        displacement = (
            system.positions[constraint.i] - system.positions[constraint.j]
        )
        relative_velocity = (
            result.velocities[constraint.i] - result.velocities[constraint.j]
        )
        assert np.dot(displacement, relative_velocity) == pytest.approx(
            0.0, abs=1e-9
        )


def test_velocity_verlet_keeps_rigid_water_constraints():
    system = build_water_box(
        2,
        box_length=3.0,
        temperature=300.0,
        seed=3,
    )

    stepped = velocity_verlet_step_system(system, dt=1e-4)

    assert stepped.constraints is system.constraints
    assert constraint_error(stepped.positions, system.constraints) < 1e-9
    for constraint in system.constraints:
        displacement = (
            stepped.positions[constraint.i] - stepped.positions[constraint.j]
        )
        relative_velocity = (
            stepped.velocities[constraint.i] - stepped.velocities[constraint.j]
        )
        assert np.dot(displacement, relative_velocity) == pytest.approx(
            0.0, abs=1e-8
        )


def test_trajectory_integrator_handles_constraints():
    system = build_water_box(
        1,
        box_length=2.0,
        temperature=300.0,
        seed=4,
    )

    trajectory = simulate_velocity_verlet_system(system, dt=1e-4, steps=20)

    assert trajectory.position.shape == (21, 3, 3)
    assert constraint_error(trajectory.position[-1], system.constraints) < 1e-9


def test_minimization_preserves_constraints():
    system = build_water_box(1, box_length=2.0, seed=5)

    result = minimize_energy(system, step_size=1e-6, max_steps=5)

    assert constraint_error(result.system.positions, system.constraints) < 1e-9


def test_barostat_preserves_constraints():
    system = build_water_box(1, box_length=2.0, seed=6)

    scaled = berendsen_barostat(
        system,
        target_pressure=pressure_bar(system) + 10.0,
        coupling_time=0.5,
        dt=0.001,
    )

    assert scaled.max_constraint_error() < 1e-9


def test_constraint_guard_accepts_a_satisfied_system():
    system = build_water_box(1, box_length=2.0, seed=7)

    assert system.max_constraint_error() < 1e-9
    assert assert_constraints_satisfied(system) < 1e-9


def test_constraint_guard_rejects_a_forgotten_shake():
    system = build_water_box(1, box_length=2.0, seed=8)
    positions = system.positions.copy()
    positions[1] += np.array([0.02, 0.01, 0.0])
    broken = System(
        particles=[
            Particle(
                particle.mass,
                positions[index],
                particle.velocity,
                atom_type=particle.atom_type,
                charge=particle.charge,
            )
            for index, particle in enumerate(system.particles)
        ],
        pair_potentials=system.pair_potentials,
        topology=system.topology,
        box=system.box,
        nonbonded_model=system.nonbonded_model,
        constraints=system.constraints,
        indexed_pair_potentials=system.indexed_pair_potentials,
    )

    assert broken.max_constraint_error() > 1e-7
    with pytest.raises(RuntimeError, match="constraints"):
        assert_constraints_satisfied(broken)


def test_flexible_water_has_no_constraints():
    system = build_water_box(1, box_length=2.0, rigid=False)

    assert system.constraints == ()


def test_constraint_input_validation():
    positions = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    masses = np.array([1.0, 1.0])
    constraints = [DistanceConstraint(0, 1, 1.0)]

    with pytest.raises(ValueError, match="shape"):
        shake_positions(np.zeros((2, 2)), masses, constraints)
    with pytest.raises(ValueError, match="masses"):
        shake_positions(positions, np.array([1.0]), constraints)
    with pytest.raises(ValueError, match="out of range"):
        shake_positions(positions, masses, [DistanceConstraint(0, 2, 1.0)])
    with pytest.raises(ValueError, match="velocities"):
        rattle_velocities(positions, np.zeros((2, 2)), masses, constraints)
    with pytest.raises(ValueError, match="max_iterations"):
        shake_positions(positions, masses, constraints, max_iterations=0)
