# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Tests for Ewald summation and its typed nonbonded integration."""

import math

import numpy as np
import pytest

from moleculab import (
    Bond,
    Box,
    EwaldElectrostatics,
    HarmonicBondPotential,
    LennardJonesParameters,
    NonbondedModel,
    Particle,
    System,
    Topology,
    units,
)


def ewald_real_space_energy(ewald, positions, charges, box, cutoff):
    total = 0.0
    for i in range(len(positions)):
        for j in range(i + 1, len(positions)):
            displacement = box.minimum_image(positions[i] - positions[j])
            distance = float(np.linalg.norm(displacement))
            if distance < cutoff:
                total += (
                    units.COULOMB_CONSTANT
                    * charges[i]
                    * charges[j]
                    * math.erfc(ewald.alpha * distance)
                    / distance
                )
    return total


def ewald_total_energy(ewald, positions, charges, box, cutoff):
    return ewald_real_space_energy(
        ewald, positions, charges, box, cutoff
    ) + ewald.global_energy(positions, charges, box)


def direct_vacuum_lattice_energy(positions, charges, box, radius):
    """Brute-force spherical real-space lattice sum with vacuum boundaries."""

    lengths = box.lengths
    n_max = int(np.ceil(radius / np.min(lengths))) + 1
    total = 0.0
    for nx in range(-n_max, n_max + 1):
        for ny in range(-n_max, n_max + 1):
            for nz in range(-n_max, n_max + 1):
                shift = np.array([nx, ny, nz]) * lengths
                if np.linalg.norm(shift) > radius + 1e-12:
                    continue
                for i in range(len(positions)):
                    for j in range(len(positions)):
                        if i == j and nx == 0 and ny == 0 and nz == 0:
                            continue
                        distance = np.linalg.norm(
                            positions[i] - positions[j] + shift
                        )
                        total += charges[i] * charges[j] / distance
    return 0.5 * units.COULOMB_CONSTANT * total


def surface_correction(positions, charges, box):
    dipole = np.sum(charges[:, np.newaxis] * positions, axis=0)
    return (
        2.0
        * np.pi
        * units.COULOMB_CONSTANT
        / (3.0 * box.volume)
        * float(np.dot(dipole, dipole))
    )


def test_ewald_real_space_energy_matches_erfc_formula():
    ewald = EwaldElectrostatics(alpha=3.0, reciprocal_cutoff=20.0)
    distance = 0.4

    expected = (
        units.COULOMB_CONSTANT
        * 1.0
        * -0.5
        * math.erfc(ewald.alpha * distance)
        / distance
    )
    assert ewald.energy(1.0, -0.5, distance, 1.0) == pytest.approx(expected)
    assert ewald.energy(1.0, -0.5, 1.1, 1.0) == 0.0


def test_ewald_real_space_force_matches_energy_gradient():
    ewald = EwaldElectrostatics(alpha=3.0, reciprocal_cutoff=20.0)
    displacement = np.array([0.4, 0.1, -0.05])
    cutoff = 1.0

    def energy_at(vector):
        distance = float(np.linalg.norm(vector))
        return ewald.energy(1.0, -0.5, distance, cutoff)

    numerical = np.empty(3)
    h = 1e-6
    for axis in range(3):
        plus = displacement.copy()
        minus = displacement.copy()
        plus[axis] += h
        minus[axis] -= h
        numerical[axis] = -(energy_at(plus) - energy_at(minus)) / (2.0 * h)

    np.testing.assert_allclose(
        ewald.force(1.0, -0.5, displacement, cutoff),
        numerical,
        atol=1e-6,
    )


def test_ewald_global_force_matches_energy_gradient():
    box = Box([2.0, 2.0, 2.0])
    positions = np.array([[0.0, 0.0, 0.0], [0.7, 0.2, 0.1]])
    charges = np.array([1.0, -1.0])
    ewald = EwaldElectrostatics(alpha=6.0, reciprocal_cutoff=40.0)

    analytical = ewald.global_forces(positions, charges, box)
    numerical = np.zeros_like(positions)
    h = 1e-6
    for atom in range(2):
        for axis in range(3):
            plus = positions.copy()
            minus = positions.copy()
            plus[atom, axis] += h
            minus[atom, axis] -= h
            numerical[atom, axis] = -(
                ewald.global_energy(plus, charges, box)
                - ewald.global_energy(minus, charges, box)
            ) / (2.0 * h)

    np.testing.assert_allclose(analytical, numerical, atol=2e-6)
    np.testing.assert_allclose(analytical.sum(axis=0), np.zeros(3), atol=1e-10)


def test_ewald_total_matches_direct_lattice_sum_with_surface_correction():
    box = Box([2.0, 2.0, 2.0])
    positions = np.array([[0.0, 0.0, 0.0], [0.7, 0.2, 0.1]])
    charges = np.array([1.0, -1.0])
    ewald = EwaldElectrostatics(alpha=6.0, reciprocal_cutoff=40.0)

    reference = direct_vacuum_lattice_energy(
        positions, charges, box, radius=12.0
    ) - surface_correction(positions, charges, box)
    ewald_energy = ewald_total_energy(
        ewald, positions, charges, box, cutoff=1.0
    )

    assert ewald_energy == pytest.approx(reference, abs=0.05)


def test_ewald_energy_converges_with_reciprocal_cutoff():
    box = Box([2.0, 2.0, 2.0])
    positions = np.array([[0.0, 0.0, 0.0], [0.7, 0.2, 0.1]])
    charges = np.array([1.0, -1.0])

    energies = []
    for reciprocal_cutoff in [20.0, 30.0, 40.0, 50.0]:
        ewald = EwaldElectrostatics(
            alpha=6.0,
            reciprocal_cutoff=reciprocal_cutoff,
        )
        energies.append(
            ewald_total_energy(ewald, positions, charges, box, cutoff=1.0)
        )

    reference = energies[-1]
    errors = [abs(energy - reference) for energy in energies]
    assert errors[0] > errors[1] > errors[2]
    assert errors[2] < 0.1


def test_ewald_charge_correction_is_zero_for_neutral_system():
    box = Box([2.0, 2.0, 2.0])
    ewald = EwaldElectrostatics(alpha=4.0, reciprocal_cutoff=20.0)

    assert ewald.charge_correction([1.0, -1.0], box) == pytest.approx(0.0)
    assert ewald.charge_correction([1.0, 1.0], box) < 0.0


def test_nonbonded_model_and_system_integrate_ewald():
    ewald = EwaldElectrostatics(alpha=6.0, reciprocal_cutoff=40.0)
    model = NonbondedModel(
        parameters={"I": LennardJonesParameters(epsilon=0.0, sigma=0.2)},
        cutoff=1.0,
        electrostatics=ewald,
    )
    positions = np.array([[0.0, 0.0, 0.0], [0.7, 0.2, 0.1]])
    charges = [1.0, -1.0]
    system = System(
        particles=[
            Particle(
                mass=1.0,
                position=positions[index],
                atom_type="I",
                charge=charges[index],
            )
            for index in range(2)
        ],
        nonbonded_model=model,
        box=Box([2.0, 2.0, 2.0]),
    )

    expected = ewald_total_energy(
        ewald,
        positions,
        np.array(charges),
        system.box,
        cutoff=1.0,
    )

    assert system.potential_energy() == pytest.approx(expected)
    np.testing.assert_allclose(system.forces().sum(axis=0), np.zeros(3), atol=1e-10)


def test_ewald_requires_a_box():
    ewald = EwaldElectrostatics(alpha=4.0, reciprocal_cutoff=20.0)
    model = NonbondedModel(
        parameters={"I": LennardJonesParameters(1.0, 0.3)},
        cutoff=1.0,
        electrostatics=ewald,
    )

    with pytest.raises(ValueError, match="periodic box"):
        System(
            particles=[Particle(1.0, [0.0, 0.0, 0.0], atom_type="I", charge=1.0)],
            nonbonded_model=model,
        )


def test_ewald_exclusions_match_self_image_reference():
    ewald = EwaldElectrostatics(alpha=6.0, reciprocal_cutoff=40.0)
    box = Box([2.0, 2.0, 2.0])
    positions = np.array([[0.0, 0.0, 0.0], [0.7, 0.2, 0.1]])
    charges = np.array([1.0, -1.0])
    model = NonbondedModel(
        parameters={"I": LennardJonesParameters(epsilon=0.0, sigma=0.2)},
        cutoff=1.0,
        electrostatics=ewald,
        exclusions=frozenset({(0, 1)}),
    )
    system = System(
        particles=[
            Particle(
                1.0,
                position,
                atom_type="I",
                charge=charges[index],
            )
            for index, position in enumerate(positions)
        ],
        nonbonded_model=model,
        box=box,
    )

    # With the pair excluded, only the self-image energy of each charge remains.
    reference = 0.0
    for index in range(2):
        single_charge = np.array([charges[index]])
        single_position = np.array([positions[index]])
        reference += ewald.global_energy(
            single_position, single_charge, box
        ) - ewald.charge_correction(single_charge, box)

    assert system.potential_energy() == pytest.approx(reference, abs=1e-9)


def test_ewald_exclusion_forces_match_energy_gradient():
    ewald = EwaldElectrostatics(alpha=6.0, reciprocal_cutoff=40.0)
    box = Box([2.0, 2.0, 2.0])
    positions = np.array([[0.0, 0.0, 0.0], [0.7, 0.2, 0.1]])
    charges = np.array([1.0, -1.0])
    model = NonbondedModel(
        parameters={"I": LennardJonesParameters(epsilon=0.0, sigma=0.2)},
        cutoff=1.0,
        electrostatics=ewald,
        exclusions=frozenset({(0, 1)}),
    )

    def make_system(coordinates):
        return System(
            particles=[
                Particle(
                    1.0,
                    coordinates[index],
                    atom_type="I",
                    charge=charges[index],
                )
                for index in range(2)
            ],
            nonbonded_model=model,
            box=box,
        )

    system = make_system(positions)
    numerical = np.zeros_like(positions)
    h = 1e-6
    for atom in range(2):
        for axis in range(3):
            plus = positions.copy()
            minus = positions.copy()
            plus[atom, axis] += h
            minus[atom, axis] -= h
            numerical[atom, axis] = -(
                make_system(plus).potential_energy()
                - make_system(minus).potential_energy()
            ) / (2.0 * h)

    np.testing.assert_allclose(system.forces(), numerical, atol=2e-6)


def test_topology_exclusions_are_applied_to_ewald():
    ewald = EwaldElectrostatics(alpha=6.0, reciprocal_cutoff=40.0)
    box = Box([2.0, 2.0, 2.0])
    positions = np.array([[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]])
    charges = np.array([1.0, -1.0])
    model = NonbondedModel(
        parameters={"I": LennardJonesParameters(epsilon=0.0, sigma=0.2)},
        cutoff=1.0,
        electrostatics=ewald,
    )
    topology = Topology(
        bonds=[Bond(0, 1, HarmonicBondPotential(1000.0, 0.5))]
    )
    system = System(
        particles=[
            Particle(1.0, positions[index], atom_type="I", charge=charges[index])
            for index in range(2)
        ],
        topology=topology,
        nonbonded_model=model,
        box=box,
    )

    reference = 0.0
    for index in range(2):
        single_charge = np.array([charges[index]])
        single_position = np.array([positions[index]])
        reference += ewald.global_energy(
            single_position, single_charge, box
        ) - ewald.charge_correction(single_charge, box)

    assert system.bonded_energy() == pytest.approx(0.0)
    assert system.potential_energy() == pytest.approx(reference, abs=1e-9)


def test_ewald_validates_parameters_and_arrays():
    with pytest.raises(ValueError, match="alpha"):
        EwaldElectrostatics(alpha=0.0, reciprocal_cutoff=10.0)
    with pytest.raises(ValueError, match="reciprocal_cutoff"):
        EwaldElectrostatics(alpha=1.0, reciprocal_cutoff=0.0)

    ewald = EwaldElectrostatics(alpha=3.0, reciprocal_cutoff=10.0)
    with pytest.raises(ValueError, match="positions"):
        ewald.reciprocal_energy(np.zeros((2, 2)), np.zeros(2), Box([1, 1, 1]))
    with pytest.raises(ValueError, match="charges"):
        ewald.reciprocal_energy(np.zeros((2, 3)), np.zeros(3), Box([1, 1, 1]))
    with pytest.raises(TypeError, match="Box"):
        ewald.reciprocal_energy(np.zeros((2, 3)), np.zeros(2), None)
