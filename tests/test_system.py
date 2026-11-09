"""Tests for multi-particle systems and pair interactions."""

import numpy as np
import pytest

from moleculab import (
    CoulombPotential,
    LennardJonesPotential,
    Particle,
    System,
)


def two_particles(position_1, position_2, velocity_1=None, velocity_2=None):
    return System(
        particles=[
            Particle(39.948, position_1, velocity_1 or [0.0, 0.0, 0.0]),
            Particle(39.948, position_2, velocity_2 or [0.0, 0.0, 0.0]),
        ]
    )


def test_system_without_potentials_has_zero_interaction_energy_and_forces():
    system = two_particles([0.0, 0.0, 0.0], [1.0, 0.0, 0.0])

    assert system.potential_energy() == 0.0
    np.testing.assert_allclose(system.forces(), np.zeros((2, 3)))


def test_system_positions_and_velocities_have_expected_shapes():
    system = two_particles(
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        velocity_1=[0.1, 0.0, 0.0],
        velocity_2=[-0.1, 0.0, 0.0],
    )

    assert system.n_particles == 2
    assert system.positions.shape == (2, 3)
    assert system.velocities.shape == (2, 3)
    assert len(system) == 2


def test_system_kinetic_energy_is_the_sum_of_particle_energies():
    system = two_particles(
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        velocity_1=[0.1, 0.0, 0.0],
        velocity_2=[-0.1, 0.0, 0.0],
    )

    expected = 0.5 * 39.948 * (0.1**2) + 0.5 * 39.948 * (0.1**2)
    assert system.kinetic_energy() == pytest.approx(expected)


def test_two_lennard_jones_particles_have_zero_energy_at_sigma():
    sigma = 0.34
    system = System(
        particles=[
            Particle(39.948, [0.0, 0.0, 0.0]),
            Particle(39.948, [sigma, 0.0, 0.0]),
        ],
        pair_potentials=[LennardJonesPotential(epsilon=1.0, sigma=sigma)],
    )

    assert system.potential_energy() == pytest.approx(0.0, abs=1e-12)
    forces = system.forces()
    np.testing.assert_allclose(forces[0], [-24.0 / sigma, 0.0, 0.0])
    np.testing.assert_allclose(forces[1], [24.0 / sigma, 0.0, 0.0])


def test_two_lennard_jones_particles_are_stable_at_the_well_minimum():
    sigma = 0.34
    epsilon = 0.997
    equilibrium = 2.0 ** (1.0 / 6.0) * sigma
    system = System(
        particles=[
            Particle(39.948, [0.0, 0.0, 0.0]),
            Particle(39.948, [equilibrium, 0.0, 0.0]),
        ],
        pair_potentials=[LennardJonesPotential(epsilon=epsilon, sigma=sigma)],
    )

    assert system.potential_energy() == pytest.approx(-epsilon)
    np.testing.assert_allclose(system.forces(), np.zeros((2, 3)), atol=1e-12)


def test_system_sums_multiple_pair_potentials():
    system = System(
        particles=[
            Particle(1.0, [0.0, 0.0, 0.0]),
            Particle(1.0, [1.0, 0.0, 0.0]),
        ],
        pair_potentials=[
            LennardJonesPotential(epsilon=1.0, sigma=1.0),
            CoulombPotential(charge_i=1.0, charge_j=1.0),
        ],
    )

    lj = LennardJonesPotential(epsilon=1.0, sigma=1.0)
    coulomb = CoulombPotential(charge_i=1.0, charge_j=1.0)
    assert system.potential_energy() == pytest.approx(
        lj.energy([1.0, 0.0, 0.0]) + coulomb.energy([1.0, 0.0, 0.0])
    )

    expected_force = lj.force([-1.0, 0.0, 0.0]) + coulomb.force([-1.0, 0.0, 0.0])
    np.testing.assert_allclose(system.forces()[0], expected_force)


def test_pair_forces_obey_newtons_third_law():
    system = System(
        particles=[
            Particle(1.0, [0.0, 0.0, 0.0]),
            Particle(4.0, [0.7, -0.2, 0.1]),
            Particle(2.0, [-0.4, 0.3, -0.2]),
        ],
        pair_potentials=[LennardJonesPotential(epsilon=1.0, sigma=0.5)],
    )

    forces = system.forces()

    np.testing.assert_allclose(forces.sum(axis=0), np.zeros(3), atol=1e-12)


def test_system_requires_at_least_one_particle():
    with pytest.raises(ValueError, match="at least one"):
        System(particles=[])


def test_system_rejects_non_particle_entries():
    with pytest.raises(TypeError, match="Particle"):
        System(particles=[object()])
