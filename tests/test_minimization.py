"""Tests for steepest-descent energy minimization."""

import numpy as np
import pytest

from moleculab import (
    Angle,
    Bond,
    HarmonicAnglePotential,
    HarmonicBondPotential,
    Particle,
    System,
    Topology,
    minimize_energy,
    steepest_descent_step,
)


def stretched_diatomic():
    particles = [
        Particle(1.0, [-0.7, 0.0, 0.0]),
        Particle(1.0, [0.7, 0.0, 0.0]),
    ]
    topology = Topology(
        bonds=[Bond(0, 1, HarmonicBondPotential(stiffness=1000.0, length0=1.0))]
    )
    return System(particles=particles, topology=topology)


def test_steepest_descent_step_lowers_energy_for_a_small_step():
    system = stretched_diatomic()
    before = system.potential_energy()

    updated = steepest_descent_step(system, step_size=1e-5)

    assert updated.potential_energy() < before
    assert system.potential_energy() == before


def test_minimization_restores_the_equilibrium_bond_length():
    result = minimize_energy(stretched_diatomic(), step_size=1e-5, max_steps=5000)

    final_positions = result.system.positions
    bond_length = float(np.linalg.norm(final_positions[0] - final_positions[1]))

    assert result.converged
    assert bond_length == pytest.approx(1.0, abs=1e-3)
    assert result.system.potential_energy() == pytest.approx(0.0, abs=1e-4)
    assert result.max_forces[-1] <= 1e-3


def test_minimization_history_has_consistent_lengths():
    result = minimize_energy(stretched_diatomic(), step_size=1e-5, max_steps=5000)

    assert result.energies.size == result.max_forces.size
    assert result.step_sizes.size == result.energies.size - 1
    assert result.steps == result.step_sizes.size
    assert np.all(np.diff(result.energies) <= 0.0)


def test_minimization_reduces_bonded_and_angle_energy():
    particles = [
        Particle(12.0, [1.2, 0.0, 0.0]),
        Particle(12.0, [0.0, 0.0, 0.0]),
        Particle(12.0, [0.0, 1.2, 0.0]),
    ]
    topology = Topology(
        bonds=[Bond(0, 1, HarmonicBondPotential(stiffness=1000.0, length0=1.0))],
        angles=[
            Angle(
                0,
                1,
                2,
                HarmonicAnglePotential(stiffness=100.0, angle0=np.pi / 2.0),
            )
        ],
    )
    system = System(particles=particles, topology=topology)
    initial_energy = system.potential_energy()

    result = minimize_energy(system, step_size=1e-5, max_steps=5000)

    assert result.system.potential_energy() < initial_energy
    assert result.system.potential_energy() == pytest.approx(0.0, abs=1e-3)


def test_already_minimized_system_converges_immediately():
    particles = [
        Particle(1.0, [-0.5, 0.0, 0.0]),
        Particle(1.0, [0.5, 0.0, 0.0]),
    ]
    topology = Topology(
        bonds=[Bond(0, 1, HarmonicBondPotential(stiffness=1000.0, length0=1.0))]
    )
    system = System(particles=particles, topology=topology)

    result = minimize_energy(system)

    assert result.converged
    assert result.steps == 0
    assert result.energies.size == 1


def test_minimization_validates_inputs():
    system = stretched_diatomic()

    with pytest.raises(ValueError, match="step_size"):
        minimize_energy(system, step_size=0.0)
    with pytest.raises(ValueError, match="max_steps"):
        minimize_energy(system, max_steps=-1)
    with pytest.raises(ValueError, match="force_tolerance"):
        minimize_energy(system, force_tolerance=-1.0)
