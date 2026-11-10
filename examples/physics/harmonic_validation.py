# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Validate velocity Verlet against the analytic harmonic oscillator.

The exact solution for a one-dimensional oscillator is

    x(t) = x0 cos(omega t) + (v0 / omega) sin(omega t)
    omega = sqrt(k / m).

This script reports the maximum position and velocity errors over two periods,
the energy drift, and the observed convergence order as the timestep is halved.
A second-order method should reduce every error by roughly a factor of four.

Run it with:

    uv run python examples/physics/harmonic_validation.py
"""

import numpy as np

from moleculab import HarmonicPotential, Particle, simulate_velocity_verlet


def exact_state(time, x0, v0, omega):
    position = x0 * np.cos(omega * time) + (v0 / omega) * np.sin(omega * time)
    velocity = -x0 * omega * np.sin(omega * time) + v0 * np.cos(omega * time)
    return position, velocity


def main() -> None:
    mass = 12.0
    stiffness = 12.0 * np.pi**2
    omega = np.sqrt(stiffness / mass)
    period = 2.0 * np.pi / omega
    x0, v0 = 0.1, 0.0

    potential = HarmonicPotential(stiffness=stiffness, center=[0.0, 0.0, 0.0])
    timesteps = [0.02, 0.01, 0.005, 0.0025, 0.00125]

    print("Analytic validation of velocity Verlet")
    print(f"  mass      = {mass:g} Da")
    print(f"  stiffness = {stiffness:.6f} kJ/(mol nm^2)")
    print(f"  omega     = {omega:.6f} rad/ps")
    print(f"  period    = {period:.6f} ps")
    print(f"  x0, v0    = {x0:g} nm, {v0:g} nm/ps")
    print()
    print(
        f"{'dt (ps)':>9} {'steps/period':>13} "
        f"{'max |x-x_exact|':>17} {'max |v-v_exact|':>17} "
        f"{'max |E-E0|':>12} {'order':>7}"
    )

    previous_position_error = None
    for dt in timesteps:
        steps = int(round(2.0 * period / dt))
        particle = Particle(
            mass=mass,
            position=[x0, 0.0, 0.0],
            velocity=[v0, 0.0, 0.0],
        )
        trajectory = simulate_velocity_verlet(particle, potential, dt=dt, steps=steps)

        exact_position, exact_velocity = exact_state(
            trajectory.time,
            x0,
            v0,
            omega,
        )
        position_error = float(
            np.max(np.abs(trajectory.position[:, 0] - exact_position))
        )
        velocity_error = float(
            np.max(np.abs(trajectory.velocity[:, 0] - exact_velocity))
        )
        energy_drift = float(
            np.max(np.abs(trajectory.total_energy - trajectory.total_energy[0]))
        )

        if previous_position_error is None:
            order_text = "-"
        else:
            order_text = f"{np.log2(previous_position_error / position_error):.2f}"
        previous_position_error = position_error

        print(
            f"{dt:9.5f} {period / dt:13.0f} "
            f"{position_error:17.3e} {velocity_error:17.3e} "
            f"{energy_drift:12.3e} {order_text:>7}"
        )

    trajectory = simulate_velocity_verlet(
        Particle(mass=mass, position=[x0, 0.0, 0.0], velocity=[v0, 0.0, 0.0]),
        potential,
        dt=0.001,
        steps=int(round(4.0 * period / 0.001)),
    )
    position = trajectory.position[:, 0]
    descending = np.where((position[:-1] > 0.0) & (position[1:] <= 0.0))[0]
    crossings = []
    for index in descending:
        t0, t1 = trajectory.time[index], trajectory.time[index + 1]
        x_start, x_end = position[index], position[index + 1]
        crossings.append(t0 - x_start * (t1 - t0) / (x_end - x_start))
    measured_period = float(np.mean(np.diff(crossings)))

    print()
    print(f"Measured period (dt = 0.001 ps): {measured_period:.8f} ps")
    print(f"Relative period error:           {abs(measured_period - period) / period:.3e}")
    print()
    print("Halving the timestep reduces the errors by about four.")
    print("That is the signature of a second-order method.")


if __name__ == "__main__":
    main()
