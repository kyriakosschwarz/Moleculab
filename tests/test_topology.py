"""Tests for molecular topology and bonded terms in a system."""

import numpy as np
import pytest

from moleculab import (
    Angle,
    Bond,
    HarmonicAnglePotential,
    HarmonicBondPotential,
    LennardJonesPotential,
    Particle,
    System,
    Topology,
    simulate_velocity_verlet_system,
    velocity_verlet_step_system,
)


def triatomic_system():
    particles = [
        Particle(12.0, [1.0, 0.0, 0.0]),
        Particle(12.0, [0.0, 0.0, 0.0]),
        Particle(12.0, [0.0, 1.0, 0.0]),
    ]
    topology = Topology(
        bonds=[Bond(0, 1, HarmonicBondPotential(stiffness=1000.0, length0=1.0))],
        angles=[Angle(0, 1, 2, HarmonicAnglePotential(stiffness=100.0, angle0=np.pi / 2))],
    )
    return System(particles=particles, topology=topology)


def test_bond_requires_two_distinct_atoms():
    potential = HarmonicBondPotential(stiffness=1000.0, length0=1.0)

    with pytest.raises(ValueError, match="different"):
        Bond(0, 0, potential)
    with pytest.raises(ValueError, match="non-negative"):
        Bond(-1, 0, potential)


def test_angle_requires_three_distinct_atoms():
    potential = HarmonicAnglePotential(stiffness=100.0, angle0=np.pi / 2)

    with pytest.raises(ValueError, match="three different"):
        Angle(0, 1, 1, potential)


def test_topology_validates_atom_indices_against_system_size():
    topology = Topology(
        bonds=[Bond(0, 2, HarmonicBondPotential(stiffness=1000.0, length0=1.0))]
    )

    with pytest.raises(ValueError, match="atom 2"):
        System(particles=[Particle(1.0, [0.0, 0.0, 0.0])], topology=topology)


def test_topology_counts_terms():
    system = triatomic_system()

    assert system.topology.n_bonds == 1
    assert system.topology.n_angles == 1


def test_topology_energy_is_zero_at_equilibrium_geometry():
    system = triatomic_system()

    assert system.bonded_energy() == pytest.approx(0.0)
    np.testing.assert_allclose(system.forces(), np.zeros((3, 3)), atol=1e-12)


def test_bonded_energy_is_sum_of_bond_and_angle_terms():
    particles = [
        Particle(12.0, [1.2, 0.0, 0.0]),
        Particle(12.0, [0.0, 0.0, 0.0]),
        Particle(12.0, [0.0, 1.0, 0.0]),
    ]
    topology = Topology(
        bonds=[Bond(0, 1, HarmonicBondPotential(stiffness=1000.0, length0=1.0))],
        angles=[Angle(0, 1, 2, HarmonicAnglePotential(stiffness=100.0, angle0=np.pi / 3))],
    )
    system = System(particles=particles, topology=topology)

    expected_bond = 0.5 * 1000.0 * (1.2 - 1.0) ** 2
    expected_angle = 0.5 * 100.0 * (np.pi / 2 - np.pi / 3) ** 2
    assert system.bonded_energy() == pytest.approx(expected_bond + expected_angle)
    assert system.potential_energy() == pytest.approx(expected_bond + expected_angle)


def test_system_combines_bonded_and_nonbonded_energy():
    system = triatomic_system()
    system.pair_potentials = (LennardJonesPotential(epsilon=1.0, sigma=1.0),)

    expected = system.topology.energy(system.positions) + system.nonbonded_energy()
    assert system.potential_energy() == pytest.approx(expected)


def test_bonded_forces_sum_to_zero():
    particles = [
        Particle(12.0, [1.2, 0.1, 0.0]),
        Particle(12.0, [0.0, 0.0, 0.0]),
        Particle(12.0, [-0.1, 0.9, 0.1]),
    ]
    topology = Topology(
        bonds=[Bond(0, 1, HarmonicBondPotential(stiffness=1000.0, length0=1.0))],
        angles=[Angle(0, 1, 2, HarmonicAnglePotential(stiffness=100.0, angle0=np.pi / 3))],
    )
    system = System(particles=particles, topology=topology)

    np.testing.assert_allclose(system.forces().sum(axis=0), np.zeros(3), atol=1e-12)


def test_system_step_preserves_topology():
    system = triatomic_system()

    updated = velocity_verlet_step_system(system, dt=0.001)

    assert updated.topology is system.topology
    assert updated.n_particles == 3


def test_velocity_verlet_conserves_energy_for_a_bonded_diatomic():
    particles = [
        Particle(1.0, [-0.5, 0.0, 0.0], [0.0, 0.1, 0.0]),
        Particle(1.0, [0.5, 0.0, 0.0], [0.0, -0.1, 0.0]),
    ]
    topology = Topology(
        bonds=[Bond(0, 1, HarmonicBondPotential(stiffness=1000.0, length0=1.0))]
    )
    system = System(particles=particles, topology=topology)

    trajectory = simulate_velocity_verlet_system(system, dt=0.001, steps=1000)

    initial_energy = trajectory.total_energy[0]
    maximum_drift = np.max(np.abs(trajectory.total_energy - initial_energy))
    assert maximum_drift / abs(initial_energy) < 1e-5
