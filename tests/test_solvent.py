"""Tests for explicit solvent models and solvent-box construction."""

import numpy as np
import pytest

from moleculab import (
    ARGON,
    COARSE_GRAINED_WATER,
    SolventModel,
    berendsen_thermostat,
    build_solvent_box,
    kinetic_temperature,
    velocity_verlet_step_system,
)


def minimum_pair_distance(system):
    positions = system.positions
    displacement = positions[:, None, :] - positions[None, :, :]
    displacement = system.box.minimum_image_array(displacement)
    distances = np.linalg.norm(displacement, axis=-1)
    upper = np.triu_indices(len(positions), k=1)
    return float(np.min(distances[upper]))


def test_solvent_model_validates_parameters():
    with pytest.raises(ValueError, match="name"):
        SolventModel("", 18.0, 0.65, 0.31)
    with pytest.raises(ValueError, match="mass"):
        SolventModel("bad", 0.0, 0.65, 0.31)
    with pytest.raises(ValueError, match="epsilon"):
        SolventModel("bad", 18.0, 0.0, 0.31)
    with pytest.raises(ValueError, match="sigma"):
        SolventModel("bad", 18.0, 0.65, -0.31)
    with pytest.raises(ValueError, match="charge"):
        SolventModel("bad", 18.0, 0.65, 0.31, float("nan"))


def test_solvent_model_builds_lennard_jones_potential():
    potentials = ARGON.pair_potentials

    assert len(potentials) == 1
    assert potentials[0].epsilon == pytest.approx(ARGON.epsilon)
    assert potentials[0].sigma == pytest.approx(ARGON.sigma)


def test_charged_solvent_adds_coulomb_potential():
    model = SolventModel("ion", 23.0, 0.5, 0.3, charge=1.0)

    assert len(model.pair_potentials) == 2


def test_build_solvent_box_from_box_length():
    system = build_solvent_box(ARGON, 8, box_length=2.0)

    assert system.n_particles == 8
    np.testing.assert_allclose(system.box.lengths, [2.0, 2.0, 2.0])
    assert system.box.volume == pytest.approx(8.0)
    assert all(particle.mass == pytest.approx(ARGON.mass) for particle in system.particles)
    assert np.all(system.positions >= 0.0)
    assert np.all(system.positions < 2.0)
    np.testing.assert_allclose(system.velocities, np.zeros((8, 3)))


def test_build_solvent_box_from_target_density():
    n_molecules = 27
    density = 5.0

    system = build_solvent_box(ARGON, n_molecules, density=density)

    expected_length = (n_molecules / density) ** (1.0 / 3.0)
    np.testing.assert_allclose(system.box.lengths, [expected_length] * 3)
    assert n_molecules / system.box.volume == pytest.approx(density)


def test_build_solvent_box_requires_exactly_one_size_argument():
    with pytest.raises(ValueError, match="exactly one"):
        build_solvent_box(ARGON, 8)
    with pytest.raises(ValueError, match="exactly one"):
        build_solvent_box(ARGON, 8, box_length=2.0, density=1.0)


def test_lattice_placement_does_not_overlap():
    system = build_solvent_box(ARGON, 64, density=4.0)

    assert minimum_pair_distance(system) >= 0.9 * ARGON.sigma


def test_random_placement_is_reproducible_and_non_overlapping():
    first = build_solvent_box(
        ARGON,
        20,
        box_length=5.0,
        placement="random",
        seed=42,
        min_distance=0.5,
    )
    second = build_solvent_box(
        ARGON,
        20,
        box_length=5.0,
        placement="random",
        seed=42,
        min_distance=0.5,
    )

    np.testing.assert_allclose(first.positions, second.positions)
    assert minimum_pair_distance(first) >= 0.5


def test_temperature_initializes_velocities():
    system = build_solvent_box(
        ARGON,
        64,
        density=3.0,
        temperature=300.0,
        seed=11,
    )

    assert kinetic_temperature(system) == pytest.approx(300.0, rel=0.2)
    np.testing.assert_allclose(np.mean(system.velocities, axis=0), np.zeros(3), atol=1e-12)


def test_solvent_box_integrates_with_thermostat_and_integrator():
    system = build_solvent_box(
        ARGON,
        27,
        density=3.0,
        temperature=300.0,
        seed=5,
    )
    initial_temperature = kinetic_temperature(system)

    heated = berendsen_thermostat(
        system,
        temperature=400.0,
        coupling_time=0.1,
        dt=0.001,
    )
    stepped = velocity_verlet_step_system(heated, dt=0.001)

    assert kinetic_temperature(heated) > initial_temperature
    assert stepped.box is system.box
    assert np.all(np.isfinite(stepped.positions))
    assert np.all(np.isfinite(stepped.velocities))


def test_solvent_builder_validates_count_and_placement():
    with pytest.raises(ValueError, match="n_molecules"):
        build_solvent_box(ARGON, 0, box_length=2.0)
    with pytest.raises(ValueError, match="n_molecules"):
        build_solvent_box(ARGON, 2.5, box_length=2.0)
    with pytest.raises(ValueError, match="placement"):
        build_solvent_box(ARGON, 8, box_length=2.0, placement="hexagonal")


def test_overdense_lattice_is_rejected():
    with pytest.raises(ValueError, match="overlapping"):
        build_solvent_box(ARGON, 125, box_length=1.0)


def test_coarse_grained_water_model_is_available():
    assert COARSE_GRAINED_WATER.mass == pytest.approx(18.015)
    assert COARSE_GRAINED_WATER.sigma > 0.0
