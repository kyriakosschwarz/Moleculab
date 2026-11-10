"""Tests for proper torsions and the dihedral angle."""

import numpy as np
import pytest

from moleculab import (
    Particle,
    PeriodicTorsionPotential,
    System,
    Topology,
    Torsion,
    dihedral_angle,
)


def torsion_potential():
    return PeriodicTorsionPotential(barrier=2.0, periodicity=3, phase=0.0)


def torsion_positions():
    return np.array(
        [
            [1.0, 0.2, 0.1],
            [0.1, -0.1, 0.0],
            [-0.2, 0.4, 0.2],
            [0.3, 0.8, -0.1],
        ]
    )


def bond_vectors(positions):
    return (
        positions[1] - positions[0],
        positions[2] - positions[1],
        positions[3] - positions[2],
    )


def test_dihedral_angle_known_geometries():
    assert dihedral_angle(
        [1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 1, 1.0],
    ) == pytest.approx(-np.pi / 2.0)
    assert dihedral_angle(
        [1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [1.0, 1.0, 0.0],
    ) == pytest.approx(0.0)
    assert abs(
        dihedral_angle(
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [-1.0, 1.0, 0.0],
        )
    ) == pytest.approx(np.pi)


def test_energy_from_angle_matches_the_threefold_curve():
    potential = torsion_potential()

    assert potential.energy_from_angle(0.0) == pytest.approx(4.0)
    assert potential.energy_from_angle(np.deg2rad(60.0)) == pytest.approx(0.0)
    assert potential.energy_from_angle(np.deg2rad(120.0)) == pytest.approx(4.0)
    assert potential.energy_from_angle(np.deg2rad(180.0)) == pytest.approx(0.0)


def test_periodic_torsion_energy_at_known_angles():
    potential = PeriodicTorsionPotential(barrier=2.0, periodicity=1, phase=0.0)

    # phi = 0 gives U = 2 * (1 + 1) = 4; phi = pi gives U = 0.
    zero_degree_positions = np.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
        ]
    )
    pi_positions = np.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [-1.0, 1.0, 0.0],
        ]
    )

    assert potential.energy(*bond_vectors(zero_degree_positions)) == pytest.approx(
        4.0, abs=1e-12
    )
    assert potential.energy(*bond_vectors(pi_positions)) == pytest.approx(0.0, abs=1e-12)


def test_energy_from_angle_is_periodic_over_a_full_rotation():
    potential = PeriodicTorsionPotential(barrier=2.0, periodicity=3)

    for angle in np.linspace(-np.pi, np.pi, 17):
        assert potential.energy_from_angle(angle) == pytest.approx(
            potential.energy_from_angle(angle + 2.0 * np.pi),
            abs=1e-12,
        )


def test_energy_from_angle_has_the_expected_threefold_barriers():
    potential = PeriodicTorsionPotential(barrier=2.0, periodicity=3)
    angles = np.linspace(0.0, 2.0 * np.pi, 3601)
    energies = np.array(
        [potential.energy_from_angle(angle) for angle in angles]
    )

    assert float(np.max(energies) - np.min(energies)) == pytest.approx(4.0)
    for minimum in [np.deg2rad(60.0), np.deg2rad(180.0), np.deg2rad(300.0)]:
        assert potential.energy_from_angle(minimum) == pytest.approx(0.0, abs=1e-12)
    for maximum in [0.0, np.deg2rad(120.0), np.deg2rad(240.0)]:
        assert potential.energy_from_angle(maximum) == pytest.approx(4.0, abs=1e-12)


def test_generalized_torsion_force_matches_the_energy_derivative():
    potential = PeriodicTorsionPotential(barrier=2.0, periodicity=3, phase=0.3)
    h = 1e-6

    for angle in [np.deg2rad(30.0), np.deg2rad(100.0), np.deg2rad(-40.0)]:
        numerical_force = -(
            potential.energy_from_angle(angle + h)
            - potential.energy_from_angle(angle - h)
        ) / (2.0 * h)
        assert potential.force_from_angle(angle) == pytest.approx(
            numerical_force,
            abs=1e-6,
        )


def test_periodic_torsion_energy_is_periodic():
    potential = PeriodicTorsionPotential(barrier=2.0, periodicity=2, phase=0.3)
    r_ij = np.array([0.0, 1.0, 0.0])
    r_jk = np.array([0.0, 0.0, 1.0])
    r_kl = np.array([1.0, 0.0, 0.0])

    angle = dihedral_angle(
        np.zeros(3),
        r_ij,
        r_ij + r_jk,
        r_ij + r_jk + r_kl,
    )
    shifted = angle + 2.0 * np.pi / potential.periodicity

    def energy_at(phi):
        return potential.barrier * (
            1.0 + np.cos(potential.periodicity * phi - potential.phase)
        )

    assert energy_at(angle) == pytest.approx(energy_at(shifted))


def test_torsion_forces_match_numerical_position_gradient():
    potential = torsion_potential()
    positions = torsion_positions()

    def energy_at(coords):
        r_ij, r_jk, r_kl = bond_vectors(coords)
        return potential.energy(r_ij, r_jk, r_kl)

    numerical_forces = np.empty_like(positions)
    h = 1e-6
    for atom in range(4):
        for axis in range(3):
            plus = positions.copy()
            minus = positions.copy()
            plus[atom, axis] += h
            minus[atom, axis] -= h
            numerical_forces[atom, axis] = -(
                energy_at(plus) - energy_at(minus)
            ) / (2.0 * h)

    r_ij, r_jk, r_kl = bond_vectors(positions)
    analytical_forces = np.array(potential.forces(r_ij, r_jk, r_kl))

    np.testing.assert_allclose(analytical_forces, numerical_forces, atol=1e-6)
    np.testing.assert_allclose(
        analytical_forces.sum(axis=0),
        np.zeros(3),
        atol=1e-12,
    )


def test_torsion_forces_vanish_at_energy_minima():
    # U = k(1 + cos(phi)); with periodicity 1 and phase 0, minima are at
    # phi = pi, where sin(phi) = 0.
    potential = PeriodicTorsionPotential(barrier=2.0, periodicity=1, phase=0.0)
    positions = np.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [-1.0, 1.0, 0.0],
        ]
    )
    r_ij, r_jk, r_kl = bond_vectors(positions)

    for force in potential.forces(r_ij, r_jk, r_kl):
        np.testing.assert_allclose(force, np.zeros(3), atol=1e-12)


def test_torsion_is_undefined_for_collinear_atoms():
    with pytest.raises(ValueError, match="collinear"):
        dihedral_angle(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
        )


def test_torsion_parameters_are_validated():
    with pytest.raises(ValueError, match="barrier"):
        PeriodicTorsionPotential(barrier=0.0, periodicity=1)
    with pytest.raises(ValueError, match="periodicity"):
        PeriodicTorsionPotential(barrier=1.0, periodicity=0)
    with pytest.raises(ValueError, match="periodicity"):
        PeriodicTorsionPotential(barrier=1.0, periodicity=1.5)


def test_topology_sums_and_applies_torsion_terms():
    positions = torsion_positions()
    particles = [Particle(12.0, position) for position in positions]
    potential = torsion_potential()
    topology = Topology(
        torsions=[Torsion(0, 1, 2, 3, potential)],
    )
    system = System(particles=particles, topology=topology)

    r_ij, r_jk, r_kl = bond_vectors(positions)
    assert topology.n_torsions == 1
    assert system.bonded_energy() == pytest.approx(
        potential.energy(r_ij, r_jk, r_kl)
    )
    np.testing.assert_allclose(system.forces(), topology.forces(positions))
    np.testing.assert_allclose(system.forces().sum(axis=0), np.zeros(3), atol=1e-12)


def test_torsion_requires_four_distinct_atoms():
    potential = torsion_potential()

    with pytest.raises(ValueError, match="four different"):
        Torsion(0, 1, 2, 1, potential)
