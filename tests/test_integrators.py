"""Tests for the velocity Verlet integrator and harmonic oscillator."""

import numpy as np
import pytest

from moleculab import (
    HarmonicPotential,
    Particle,
    simulate_velocity_verlet,
    total_energy,
    velocity_verlet_step,
)


def oscillator():
    """Return a harmonic oscillator with period exactly 2 ps.

    For a harmonic oscillator, omega = sqrt(k/m) and T = 2*pi/omega.
    Choosing k = 12*pi^2 with m = 12 gives omega = pi and T = 2 ps.
    """

    mass = 12.0
    stiffness = 12.0 * np.pi**2
    particle = Particle(mass=mass, position=[0.1, 0.0, 0.0], velocity=[0.0, 0.05, 0.0])
    potential = HarmonicPotential(stiffness=stiffness, center=[0.0, 0.0, 0.0])
    return particle, potential


def test_velocity_verlet_step_does_not_mutate_input():
    particle, potential = oscillator()
    position_before = particle.position.copy()
    velocity_before = particle.velocity.copy()

    velocity_verlet_step(particle, potential, dt=0.001)

    np.testing.assert_allclose(particle.position, position_before)
    np.testing.assert_allclose(particle.velocity, velocity_before)


def test_velocity_verlet_step_preserves_chemical_identity():
    particle = Particle(
        mass=12.0,
        position=[0.1, 0.0, 0.0],
        velocity=[0.0, 0.0, 0.0],
        atom_type="carbon",
        charge=-0.2,
    )
    potential = HarmonicPotential(stiffness=100.0, center=[0.0, 0.0, 0.0])

    updated = velocity_verlet_step(particle, potential, dt=0.001)

    assert updated.atom_type == "carbon"
    assert updated.charge == pytest.approx(-0.2)


def test_one_step_matches_the_velocity_verlet_equations():
    particle = Particle(mass=2.0, position=[0.1, 0.0, 0.0], velocity=[0.0, 0.0, 0.0])
    potential = HarmonicPotential(stiffness=100.0, center=[0.0, 0.0, 0.0])
    dt = 0.01

    acceleration_now = potential.force(particle.position) / particle.mass
    expected_position = (
        particle.position + particle.velocity * dt + 0.5 * acceleration_now * dt**2
    )
    acceleration_next = potential.force(expected_position) / particle.mass
    expected_velocity = particle.velocity + 0.5 * (
        acceleration_now + acceleration_next
    ) * dt

    updated = velocity_verlet_step(particle, potential, dt)

    np.testing.assert_allclose(updated.position, expected_position)
    np.testing.assert_allclose(updated.velocity, expected_velocity)


def test_trajectory_has_one_entry_per_state_including_start():
    particle, potential = oscillator()

    trajectory = simulate_velocity_verlet(particle, potential, dt=0.001, steps=10)

    assert len(trajectory) == 11
    assert trajectory.time.shape == (11,)
    assert trajectory.position.shape == (11, 3)
    assert trajectory.velocity.shape == (11, 3)
    assert trajectory.total_energy.shape == (11,)


def test_trajectory_starts_at_the_initial_state():
    particle, potential = oscillator()

    trajectory = simulate_velocity_verlet(particle, potential, dt=0.001, steps=5)

    np.testing.assert_allclose(trajectory.position[0], particle.position)
    np.testing.assert_allclose(trajectory.velocity[0], particle.velocity)
    assert trajectory.kinetic_energy[0] == pytest.approx(
        total_energy(particle, potential) - potential.energy(particle.position)
    )


def test_velocity_verlet_conserves_energy_for_a_harmonic_oscillator():
    particle, potential = oscillator()

    trajectory = simulate_velocity_verlet(particle, potential, dt=0.001, steps=2000)

    initial_energy = trajectory.total_energy[0]
    maximum_drift = np.max(np.abs(trajectory.total_energy - initial_energy))
    relative_drift = maximum_drift / abs(initial_energy)
    assert relative_drift < 1e-5


def test_harmonic_oscillator_returns_near_its_start_after_one_period():
    particle, potential = oscillator()

    # The period is 2 ps, and dt = 0.001 ps gives exactly 2000 steps.
    trajectory = simulate_velocity_verlet(particle, potential, dt=0.001, steps=2000)

    np.testing.assert_allclose(trajectory.position[-1], particle.position, atol=2e-4)
    np.testing.assert_allclose(trajectory.velocity[-1], particle.velocity, atol=2e-3)


def test_zero_steps_returns_only_the_initial_snapshot():
    particle, potential = oscillator()

    trajectory = simulate_velocity_verlet(particle, potential, dt=0.001, steps=0)

    assert len(trajectory) == 1
    np.testing.assert_allclose(trajectory.position[0], particle.position)


@pytest.mark.parametrize("dt", [0.0, -0.001, float("nan"), float("inf")])
def test_timestep_must_be_finite_and_positive(dt):
    particle, potential = oscillator()

    with pytest.raises(ValueError, match="dt"):
        simulate_velocity_verlet(particle, potential, dt=dt, steps=1)


@pytest.mark.parametrize("steps", [-1, 1.5, "5"])
def test_steps_must_be_a_non_negative_integer(steps):
    particle, potential = oscillator()

    with pytest.raises((TypeError, ValueError), match="steps"):
        simulate_velocity_verlet(particle, potential, dt=0.001, steps=steps)
