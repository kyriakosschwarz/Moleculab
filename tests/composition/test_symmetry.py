# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Symmetry and combined-system validation.

These tests check properties that are deeper than a single force gradient:

* translating every particle does not change the internal energy;
* rotating every particle does not change the internal energy;
* forces rotate with the coordinates;
* the total force from internal interactions is zero;
* wrapping and lattice translations preserve periodic pair interactions.

The exact identities are separated from numerical expectations in the tutorials
and in ``docs/reference/validation.md``.
"""

import numpy as np
import pytest

from moleculab import (
    Angle,
    Bond,
    Box,
    HarmonicAnglePotential,
    HarmonicBondPotential,
    LennardJonesParameters,
    LennardJonesPotential,
    NonbondedModel,
    Particle,
    PeriodicTorsionPotential,
    ReactionFieldElectrostatics,
    System,
    Topology,
    Torsion,
)


def rotation_matrix(axis, angle):
    axis = np.asarray(axis, dtype=float)
    axis = axis / np.linalg.norm(axis)
    cross = np.array(
        [
            [0.0, -axis[2], axis[1]],
            [axis[2], 0.0, -axis[0]],
            [-axis[1], axis[0], 0.0],
        ]
    )
    return (
        np.eye(3)
        + np.sin(angle) * cross
        + (1.0 - np.cos(angle)) * (cross @ cross)
    )


def rebuild(system, positions):
    return System(
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


def molecular_system():
    particles = [
        Particle(12.0, [0.00, 0.00, 0.00], atom_type="C", charge=-0.20),
        Particle(12.0, [0.15, 0.03, -0.01], atom_type="C", charge=0.10),
        Particle(16.0, [0.23, 0.17, 0.02], atom_type="O", charge=-0.40),
        Particle(1.0, [0.19, 0.21, 0.16], atom_type="H", charge=0.20),
        Particle(1.0, [0.40, -0.15, 0.10], atom_type="H", charge=0.10),
    ]
    topology = Topology(
        bonds=[
            Bond(0, 1, HarmonicBondPotential(1000.0, 0.15)),
            Bond(1, 2, HarmonicBondPotential(800.0, 0.14)),
            Bond(2, 3, HarmonicBondPotential(3000.0, 0.10)),
        ],
        angles=[
            Angle(0, 1, 2, HarmonicAnglePotential(100.0, np.deg2rad(109.5))),
            Angle(1, 2, 3, HarmonicAnglePotential(200.0, np.deg2rad(104.5))),
        ],
        torsions=[
            Torsion(
                0,
                1,
                2,
                3,
                PeriodicTorsionPotential(barrier=2.0, periodicity=3),
            )
        ],
    )
    model = NonbondedModel(
        parameters={
            "C": LennardJonesParameters(epsilon=0.5, sigma=0.34),
            "H": LennardJonesParameters(epsilon=0.1, sigma=0.25),
            "O": LennardJonesParameters(epsilon=0.65, sigma=0.30),
        },
        cutoff=2.0,
        electrostatics=ReactionFieldElectrostatics(dielectric=78.3),
    )
    return System(
        particles=particles,
        topology=topology,
        nonbonded_model=model,
    )


def test_translation_preserves_internal_energy_and_forces():
    system = molecular_system()
    shift = np.array([0.37, -0.42, 0.19])
    translated = rebuild(system, system.positions + shift)

    assert translated.potential_energy() == pytest.approx(
        system.potential_energy(),
        abs=1e-10,
    )
    np.testing.assert_allclose(
        translated.forces(),
        system.forces(),
        atol=1e-10,
    )


def test_rotation_preserves_internal_energy_and_rotates_forces():
    system = molecular_system()
    matrix = rotation_matrix([1.0, 2.0, 3.0], 0.7)

    np.testing.assert_allclose(matrix @ matrix.T, np.eye(3), atol=1e-12)
    assert np.linalg.det(matrix) == pytest.approx(1.0)

    rotated_positions = system.positions @ matrix.T
    rotated = rebuild(system, rotated_positions)

    assert rotated.potential_energy() == pytest.approx(
        system.potential_energy(),
        abs=1e-10,
    )
    np.testing.assert_allclose(
        rotated.forces(),
        system.forces() @ matrix.T,
        atol=1e-10,
    )


def test_total_internal_force_is_zero():
    system = molecular_system()

    np.testing.assert_allclose(
        np.sum(system.forces(), axis=0),
        np.zeros(3),
        atol=1e-10,
    )


def test_pair_forces_are_equal_and_opposite():
    system = System(
        particles=[
            Particle(12.0, [0.0, 0.0, 0.0]),
            Particle(12.0, [0.4, 0.1, -0.05]),
        ],
        pair_potentials=[
            LennardJonesPotential(0.997, 0.34),
        ],
    )

    forces = system.forces()
    np.testing.assert_allclose(forces[0], -forces[1], atol=1e-12)


def test_wrapping_preserves_periodic_pair_interaction():
    box = Box([2.0, 2.0, 2.0])
    system = System(
        particles=[
            Particle(12.0, [-0.2, 0.5, 0.5]),
            Particle(12.0, [1.7, 0.5, 0.5]),
        ],
        pair_potentials=[LennardJonesPotential(1.0, 0.34)],
        box=box,
    )
    wrapped_positions = np.array(
        [box.wrap(particle.position) for particle in system.particles]
    )
    wrapped = rebuild(system, wrapped_positions)

    assert wrapped.potential_energy() == pytest.approx(
        system.potential_energy(),
        abs=1e-10,
    )
    np.testing.assert_allclose(
        wrapped.forces(),
        system.forces(),
        atol=1e-10,
    )


def test_lattice_translation_preserves_periodic_energy():
    box = Box([2.0, 2.0, 2.0])
    system = System(
        particles=[
            Particle(12.0, [0.2, 0.5, 0.5]),
            Particle(12.0, [1.7, 0.5, 0.5]),
        ],
        pair_potentials=[LennardJonesPotential(1.0, 0.34)],
        box=box,
    )
    translated = rebuild(system, system.positions + np.array([2.0, 0.0, 0.0]))

    assert translated.potential_energy() == pytest.approx(
        system.potential_energy(),
        abs=1e-10,
    )


def test_cubic_box_quarter_turn_preserves_energy_and_rotates_forces():
    box = Box([2.0, 2.0, 2.0])
    system = System(
        particles=[
            Particle(12.0, [0.3, 0.4, 0.5]),
            Particle(12.0, [1.6, 1.2, 1.5]),
        ],
        pair_potentials=[LennardJonesPotential(1.0, 0.34)],
        box=box,
    )
    matrix = rotation_matrix([0.0, 0.0, 1.0], np.pi / 2.0)
    rotated = rebuild(system, system.positions @ matrix.T)

    assert rotated.potential_energy() == pytest.approx(
        system.potential_energy(),
        abs=1e-10,
    )
    np.testing.assert_allclose(
        rotated.forces(),
        system.forces() @ matrix.T,
        atol=1e-10,
    )
