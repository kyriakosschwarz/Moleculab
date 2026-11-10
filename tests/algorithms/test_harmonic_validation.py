# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Analytic validation of the velocity Verlet harmonic oscillator.

These tests compare the simulated trajectory with the exact solution,

    x(t) = x0 cos(omega t) + (v0 / omega) sin(omega t),

rather than only checking that energy looks approximately constant.  The
second-order convergence tests make the numerical method itself part of the
contract.
"""

import numpy as np
import pytest

from moleculab import HarmonicPotential, Particle, simulate_velocity_verlet

MASS = 12.0
STIFFNESS = 12.0 * np.pi**2
OMEGA = np.sqrt(STIFFNESS / MASS)
PERIOD = 2.0 * np.pi / OMEGA


def exact_state(time, x0, v0):
    """Return the analytic position and velocity for a 1D oscillator."""

    time = np.asarray(time, dtype=float)
    position = x0 * np.cos(OMEGA * time) + (v0 / OMEGA) * np.sin(OMEGA * time)
    velocity = -x0 * OMEGA * np.sin(OMEGA * time) + v0 * np.cos(OMEGA * time)
    return position, velocity


def run_oscillator(dt, periods=2.0, x0=0.1, v0=0.0):
    particle = Particle(
        mass=MASS,
        position=[x0, 0.0, 0.0],
        velocity=[v0, 0.0, 0.0],
    )
    potential = HarmonicPotential(stiffness=STIFFNESS, center=[0.0, 0.0, 0.0])
    steps = round(periods * PERIOD / dt)
    return simulate_velocity_verlet(particle, potential, dt=dt, steps=steps)


def maximum_errors(dt, periods=2.0):
    trajectory = run_oscillator(dt, periods=periods)
    exact_position, exact_velocity = exact_state(trajectory.time, 0.1, 0.0)
    position_error = np.max(np.abs(trajectory.position[:, 0] - exact_position))
    velocity_error = np.max(np.abs(trajectory.velocity[:, 0] - exact_velocity))
    energy_drift = np.max(
        np.abs(trajectory.total_energy - trajectory.total_energy[0])
    )
    return position_error, velocity_error, energy_drift


def test_general_phase_matches_the_analytic_solution():
    x0, v0 = 0.08, 0.03
    trajectory = run_oscillator(dt=0.001, periods=2.0, x0=x0, v0=v0)
    exact_position, exact_velocity = exact_state(trajectory.time, x0, v0)

    np.testing.assert_allclose(
        trajectory.position[:, 0],
        exact_position,
        atol=1e-6,
    )
    np.testing.assert_allclose(
        trajectory.velocity[:, 0],
        exact_velocity,
        atol=3e-6,
    )
    np.testing.assert_allclose(trajectory.position[:, 1:], 0.0, atol=1e-15)


def test_multidimensional_oscillator_decouples():
    x0, vx0 = 0.08, 0.03
    y0, vy0 = -0.05, 0.02
    particle = Particle(
        mass=MASS,
        position=[x0, y0, 0.0],
        velocity=[vx0, vy0, 0.0],
    )
    potential = HarmonicPotential(stiffness=STIFFNESS, center=[0.0, 0.0, 0.0])
    trajectory = simulate_velocity_verlet(
        particle,
        potential,
        dt=0.001,
        steps=4000,
    )
    exact_x, _ = exact_state(trajectory.time, x0, vx0)
    exact_y, _ = exact_state(trajectory.time, y0, vy0)

    np.testing.assert_allclose(trajectory.position[:, 0], exact_x, atol=1e-6)
    np.testing.assert_allclose(trajectory.position[:, 1], exact_y, atol=1e-6)
    np.testing.assert_allclose(trajectory.position[:, 2], 0.0, atol=1e-15)


def test_position_error_is_second_order_in_the_timestep():
    coarse, _, _ = maximum_errors(dt=0.002)
    fine, _, _ = maximum_errors(dt=0.001)

    assert coarse / fine == pytest.approx(4.0, rel=0.05)


def test_velocity_error_is_second_order_in_the_timestep():
    _, coarse, _ = maximum_errors(dt=0.002)
    _, fine, _ = maximum_errors(dt=0.001)

    assert coarse / fine == pytest.approx(4.0, rel=0.05)


def test_energy_drift_is_second_order_in_the_timestep():
    _, _, coarse = maximum_errors(dt=0.002)
    _, _, fine = maximum_errors(dt=0.001)

    assert coarse / fine == pytest.approx(4.0, rel=0.05)


def test_numerical_period_matches_the_analytic_period():
    trajectory = run_oscillator(dt=0.001, periods=4.0)
    position = trajectory.position[:, 0]
    descending = np.where((position[:-1] > 0.0) & (position[1:] <= 0.0))[0]

    crossings = []
    for index in descending:
        t0, t1 = trajectory.time[index], trajectory.time[index + 1]
        x0, x1 = position[index], position[index + 1]
        crossings.append(t0 - x0 * (t1 - t0) / (x1 - x0))

    measured_period = float(np.mean(np.diff(crossings)))
    assert measured_period == pytest.approx(PERIOD, rel=1e-5)


def test_tiny_timestep_approaches_the_analytic_trajectory():
    position_error, _, _ = maximum_errors(dt=1e-4, periods=1.0)

    assert position_error < 1e-8


def test_numerical_energy_matches_the_analytic_energy():
    x0, v0 = 0.08, 0.03
    trajectory = run_oscillator(dt=0.001, periods=2.0, x0=x0, v0=v0)
    analytic_energy = 0.5 * MASS * v0**2 + 0.5 * STIFFNESS * x0**2

    assert trajectory.total_energy[0] == pytest.approx(analytic_energy)
    np.testing.assert_allclose(
        trajectory.total_energy,
        analytic_energy,
        rtol=1e-5,
    )
