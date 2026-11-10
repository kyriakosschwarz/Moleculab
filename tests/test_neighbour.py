# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Tests for the optional cell-list neighbour generator.

The cell list is a candidate screen, not a force calculation.  Its contract is:

* every true pair within the cutoff is present;
* the list is duplicate-free and uses ``i < j``;
* filtering the candidates by the exact minimum-image distance reproduces the
  brute-force pair set;
* periodic wrapping does not lose pairs across box faces.
"""

import numpy as np
import pytest

from moleculab import (
    Box,
    EwaldElectrostatics,
    LennardJonesParameters,
    LennardJonesPotential,
    NonbondedModel,
    Particle,
    ReactionFieldElectrostatics,
    System,
    build_neighbour_list,
)


def brute_force_pairs(positions, cutoff, box=None):
    pairs = set()
    for i in range(len(positions)):
        for j in range(i + 1, len(positions)):
            displacement = positions[i] - positions[j]
            if box is not None:
                displacement = box.minimum_image(displacement)
            if np.linalg.norm(displacement) < cutoff:
                pairs.add((i, j))
    return pairs


def filtered_pairs(positions, neighbours, cutoff, box=None):
    pairs = set()
    for i, j in neighbours:
        displacement = positions[i] - positions[j]
        if box is not None:
            displacement = box.minimum_image(displacement)
        if np.linalg.norm(displacement) < cutoff:
            pairs.add((i, j))
    return pairs


def test_nonperiodic_candidates_cover_brute_force():
    rng = np.random.default_rng(11)
    positions = rng.uniform(0.0, 5.0, size=(60, 3))
    cutoff = 0.7

    neighbours = build_neighbour_list(positions, cutoff)

    brute = brute_force_pairs(positions, cutoff)
    candidates = set(neighbours)
    assert brute <= candidates
    assert filtered_pairs(positions, neighbours, cutoff) == brute
    assert len(neighbours) == len(candidates)
    assert all(i < j for i, j in neighbours)
    assert list(neighbours) == sorted(neighbours)


def test_periodic_candidates_cover_brute_force_across_faces():
    box = Box([2.0, 2.0, 2.0])
    rng = np.random.default_rng(12)
    positions = rng.uniform(-1.0, 3.0, size=(80, 3))
    cutoff = 0.6

    neighbours = build_neighbour_list(positions, cutoff, box=box)

    brute = brute_force_pairs(positions, cutoff, box=box)
    candidates = set(neighbours)
    assert brute <= candidates
    assert filtered_pairs(positions, neighbours, cutoff, box=box) == brute
    assert len(neighbours) == len(candidates)
    assert neighbours.periodic is True


def test_periodic_pair_across_the_boundary_is_found():
    box = Box([10.0, 10.0, 10.0])
    positions = np.array([[0.05, 5.0, 5.0], [9.95, 5.0, 5.0]])

    neighbours = build_neighbour_list(positions, 0.2, box=box)

    assert (0, 1) in neighbours
    assert filtered_pairs(positions, neighbours, 0.2, box=box) == {(0, 1)}


def test_cell_list_is_a_screen_not_a_filter():
    positions = np.array([[0.0, 0.0, 0.0], [1.9, 0.0, 0.0]])
    cutoff = 2.0

    neighbours = build_neighbour_list(positions, cutoff)

    # Both particles are within the cutoff, so the pair must be present.
    assert (0, 1) in neighbours


def test_cell_list_scales_better_than_all_pairs_for_separated_points():
    rng = np.random.default_rng(13)
    positions = rng.uniform(0.0, 20.0, size=(200, 3))
    cutoff = 0.5

    neighbours = build_neighbour_list(positions, cutoff)

    total_pairs = len(positions) * (len(positions) - 1) // 2
    assert len(neighbours) < total_pairs
    assert filtered_pairs(positions, neighbours, cutoff) == brute_force_pairs(
        positions, cutoff
    )


def test_periodic_cutoff_must_fit_the_box():
    positions = np.array([[0.0, 0.0, 0.0], [0.1, 0.0, 0.0]])
    box = Box([1.0, 1.0, 1.0])

    with pytest.raises(ValueError, match="half the shortest box edge"):
        build_neighbour_list(positions, 0.6, box=box)


def test_neighbour_input_validation():
    positions = np.array([[0.0, 0.0, 0.0], [0.1, 0.0, 0.0]])

    with pytest.raises(ValueError, match="shape"):
        build_neighbour_list(np.zeros((2, 2)), 0.2)
    with pytest.raises(ValueError, match="cutoff"):
        build_neighbour_list(positions, 0.0)
    with pytest.raises(ValueError, match="finite"):
        build_neighbour_list(
            np.array([[0.0, 0.0, 0.0], [np.nan, 0.0, 0.0]]),
            0.2,
        )
    with pytest.raises(ValueError, match="at least one"):
        build_neighbour_list(np.empty((0, 3)), 0.2)
    with pytest.raises(TypeError, match="Box"):
        build_neighbour_list(positions, 0.2, box="not a box")


def pair_system():
    rng = np.random.default_rng(31)
    positions = rng.uniform(0.0, 3.0, size=(24, 3))
    return System(
        particles=[Particle(12.0, position) for position in positions],
        pair_potentials=[LennardJonesPotential(epsilon=1.0, sigma=0.3)],
        box=Box([3.0, 3.0, 3.0]),
    )


def test_system_neighbour_list_energy_and_force_parity():
    system = pair_system()
    cutoff = 0.8

    energy_brute = system.potential_energy(cutoff=cutoff)
    energy_neighbour = system.potential_energy(
        cutoff=cutoff,
        use_neighbour_list=True,
    )
    forces_brute = system.forces(cutoff=cutoff)
    forces_neighbour = system.forces(
        cutoff=cutoff,
        use_neighbour_list=True,
    )

    assert energy_neighbour == pytest.approx(energy_brute, abs=1e-12)
    np.testing.assert_allclose(forces_neighbour, forces_brute, atol=1e-12)


def test_system_neighbour_list_requires_cutoff_without_a_model():
    system = pair_system()

    with pytest.raises(ValueError, match="cutoff is required"):
        system.forces(use_neighbour_list=True)


def test_system_neighbour_list_defaults_to_the_model_cutoff():
    model = NonbondedModel(
        parameters={
            "C": LennardJonesParameters(epsilon=0.5, sigma=0.34),
            "H": LennardJonesParameters(epsilon=0.1, sigma=0.25),
        },
        cutoff=0.6,
        electrostatics=ReactionFieldElectrostatics(dielectric=78.3),
    )
    system = System(
        particles=[
            Particle(12.0, [0.0, 0.0, 0.0], atom_type="C", charge=-0.2),
            Particle(1.0, [0.4, 0.0, 0.0], atom_type="H", charge=0.2),
            Particle(1.0, [0.7, 0.2, 0.0], atom_type="H", charge=0.1),
        ],
        nonbonded_model=model,
        box=Box([2.0, 2.0, 2.0]),
    )

    assert system.potential_energy(
        use_neighbour_list=True
    ) == pytest.approx(system.potential_energy(), abs=1e-12)
    np.testing.assert_allclose(
        system.forces(use_neighbour_list=True),
        system.forces(),
        atol=1e-12,
    )


def test_system_neighbour_list_rejects_cutoff_below_model_cutoff():
    model = NonbondedModel(
        parameters={"C": LennardJonesParameters(epsilon=0.5, sigma=0.34)},
        cutoff=0.8,
    )
    system = System(
        particles=[
            Particle(12.0, [0.0, 0.0, 0.0], atom_type="C"),
            Particle(12.0, [0.5, 0.0, 0.0], atom_type="C"),
        ],
        nonbonded_model=model,
        box=Box([2.0, 2.0, 2.0]),
    )

    with pytest.raises(ValueError, match="at least the nonbonded model cutoff"):
        system.forces(cutoff=0.5, use_neighbour_list=True)


def test_indexed_pairs_are_not_screened_by_the_neighbour_cutoff():
    system = System(
        particles=[
            Particle(12.0, [0.0, 0.0, 0.0]),
            Particle(12.0, [1.5, 0.0, 0.0]),
        ],
        indexed_pair_potentials=[
            (0, 1, LennardJonesPotential(epsilon=1.0, sigma=0.3))
        ],
    )

    brute = system.forces()
    screened = system.forces(cutoff=0.2, use_neighbour_list=True)

    np.testing.assert_allclose(screened, brute, atol=1e-12)


def test_ewald_global_terms_are_preserved_by_the_neighbour_path():
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

    assert system.potential_energy(
        use_neighbour_list=True
    ) == pytest.approx(system.potential_energy(), abs=1e-12)
    np.testing.assert_allclose(
        system.forces(use_neighbour_list=True),
        system.forces(),
        atol=1e-12,
    )
