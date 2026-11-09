"""Tests for the particle representation."""

import numpy as np
import pytest

from moleculab import Particle


def test_particle_stores_mass_position_and_velocity():
    particle = Particle(
        mass=12.0,
        position=[0.1, 0.2, 0.3],
        velocity=[0.4, 0.5, 0.6],
    )

    assert particle.mass == 12.0
    np.testing.assert_allclose(particle.position, [0.1, 0.2, 0.3])
    np.testing.assert_allclose(particle.velocity, [0.4, 0.5, 0.6])


def test_speed_is_the_norm_of_velocity():
    particle = Particle(mass=1.0, position=[0.0, 0.0, 0.0], velocity=[3.0, 4.0, 0.0])

    assert particle.speed == pytest.approx(5.0)


@pytest.mark.parametrize("mass", [0.0, -1.0, float("nan"), float("inf")])
def test_mass_must_be_finite_and_positive(mass):
    with pytest.raises(ValueError, match="mass"):
        Particle(mass=mass, position=[0.0, 0.0, 0.0])


def test_position_must_have_three_components():
    with pytest.raises(ValueError, match="position"):
        Particle(mass=1.0, position=[0.0, 0.0])


def test_velocity_must_have_three_components():
    with pytest.raises(ValueError, match="velocity"):
        Particle(mass=1.0, position=[0.0, 0.0, 0.0], velocity=[0.0, 0.0])


def test_non_finite_coordinates_are_rejected():
    with pytest.raises(ValueError, match="position"):
        Particle(mass=1.0, position=[0.0, np.nan, 0.0])


def test_constructor_copies_input_arrays():
    position = np.array([0.1, 0.0, 0.0])
    velocity = np.array([0.0, 0.1, 0.0])

    particle = Particle(mass=1.0, position=position, velocity=velocity)
    position[0] = 99.0
    velocity[1] = 99.0

    np.testing.assert_allclose(particle.position, [0.1, 0.0, 0.0])
    np.testing.assert_allclose(particle.velocity, [0.0, 0.1, 0.0])


def test_copy_is_independent():
    original = Particle(mass=1.0, position=[1.0, 0.0, 0.0], velocity=[0.0, 1.0, 0.0])
    duplicate = original.copy()

    duplicate.position[0] = 5.0

    assert original.position[0] == 1.0
    assert duplicate.mass == original.mass


def test_atom_type_and_charge_defaults_and_copy():
    default = Particle(mass=1.0, position=[0.0, 0.0, 0.0])
    assert default.atom_type == "default"
    assert default.charge == 0.0

    typed = Particle(
        mass=1.0,
        position=[0.0, 0.0, 0.0],
        atom_type="O",
        charge=-0.8476,
    )
    duplicate = typed.copy()

    assert duplicate.atom_type == "O"
    assert duplicate.charge == pytest.approx(-0.8476)


def test_atom_type_must_be_non_empty_and_charge_finite():
    with pytest.raises(ValueError, match="atom_type"):
        Particle(mass=1.0, position=[0.0, 0.0, 0.0], atom_type=" ")
    with pytest.raises(ValueError, match="charge"):
        Particle(mass=1.0, position=[0.0, 0.0, 0.0], charge=float("nan"))
