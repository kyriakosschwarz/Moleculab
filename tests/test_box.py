"""Tests for orthorhombic periodic boundaries."""

import numpy as np
import pytest

from moleculab import (
    Box,
    LennardJonesPotential,
    Particle,
    System,
    berendsen_thermostat,
    minimize_energy,
    velocity_verlet_step_system,
)


def test_box_volume_and_repr():
    box = Box([2.0, 3.0, 4.0])

    assert box.volume == pytest.approx(24.0)
    assert "nm" in repr(box)


def test_box_requires_positive_lengths():
    with pytest.raises(ValueError, match="positive"):
        Box([1.0, 0.0, 1.0])
    with pytest.raises(ValueError, match="lengths"):
        Box([1.0, 1.0])


def test_minimum_image_uses_the_shortest_separation():
    box = Box([1.0, 1.0, 1.0])

    np.testing.assert_allclose(box.minimum_image([0.9, 0.0, 0.0]), [-0.1, 0.0, 0.0])
    np.testing.assert_allclose(box.minimum_image([-0.9, 0.0, 0.0]), [0.1, 0.0, 0.0])
    assert box.distance([0.9, 0.0, 0.0]) == pytest.approx(0.1)


def test_wrap_maps_positions_into_the_primary_box():
    box = Box([1.0, 1.0, 1.0])

    np.testing.assert_allclose(box.wrap([1.2, -0.2, 2.1]), [0.2, 0.8, 0.1])
    np.testing.assert_allclose(box.wrap([0.5, 0.5, 0.5]), [0.5, 0.5, 0.5])


def test_periodic_nonbonded_interaction_matches_the_minimum_image_pair():
    sigma = 0.1
    potential = LennardJonesPotential(epsilon=1.0, sigma=sigma)

    periodic = System(
        particles=[
            Particle(1.0, [0.05, 0.0, 0.0]),
            Particle(1.0, [0.95, 0.0, 0.0]),
        ],
        pair_potentials=[potential],
        box=Box([1.0, 1.0, 1.0]),
    )
    direct = System(
        particles=[
            Particle(1.0, [0.1, 0.0, 0.0]),
            Particle(1.0, [0.0, 0.0, 0.0]),
        ],
        pair_potentials=[potential],
    )

    assert periodic.potential_energy() == pytest.approx(direct.potential_energy())
    np.testing.assert_allclose(periodic.forces(), direct.forces())


def test_periodic_system_preserves_box_through_updates():
    system = System(
        particles=[Particle(1.0, [0.1, 0.1, 0.1], [0.1, 0.0, 0.0])],
        box=Box([1.0, 2.0, 3.0]),
    )

    stepped = velocity_verlet_step_system(system, dt=0.001)
    minimized = minimize_energy(system).system
    thermostatted = berendsen_thermostat(
        system, 300.0, coupling_time=0.1, dt=0.001
    )

    assert stepped.box is system.box
    assert minimized.box is system.box
    assert thermostatted.box is system.box


def test_system_rejects_invalid_box():
    with pytest.raises(TypeError, match="Box"):
        System(particles=[Particle(1.0, [0.0, 0.0, 0.0])], box="not a box")
