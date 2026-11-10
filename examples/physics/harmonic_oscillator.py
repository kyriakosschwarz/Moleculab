# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Simulate a harmonic oscillator and inspect numerical energy conservation.

The oscillator has mass 12 Da and stiffness 12*pi^2 kJ/(mol nm^2), so

    omega = sqrt(k/m) = pi rad/ps
    T     = 2*pi/omega = 2 ps.

The initial displacement is 0.1 nm and the initial velocity is zero.  The
exact solution is x(t) = 0.1 cos(pi t) nm.

Run it with:

    uv run python examples/physics/harmonic_oscillator.py
"""

import numpy as np

from moleculab import HarmonicPotential, Particle, simulate_velocity_verlet


def main() -> None:
    mass = 12.0
    stiffness = 12.0 * np.pi**2
    particle = Particle(mass=mass, position=[0.1, 0.0, 0.0], velocity=[0.0, 0.0, 0.0])
    potential = HarmonicPotential(stiffness=stiffness, center=[0.0, 0.0, 0.0])

    dt = 0.001
    trajectory = simulate_velocity_verlet(particle, potential, dt=dt, steps=2000)

    print("Harmonic oscillator: mass = 12 Da, k = 12*pi^2 kJ/(mol nm^2)")
    print(f"Time step: {dt} ps; steps: {len(trajectory) - 1}; period: 2.0 ps")
    print()
    print(f"{'time (ps)':>10} {'x (nm)':>12} {'exact x':>12} {'E (kJ/mol)':>14}")
    for index in [0, 250, 500, 750, 1000, 1250, 1500, 1750, 2000]:
        time = trajectory.time[index]
        exact = 0.1 * np.cos(np.pi * time)
        print(
            f"{time:10.3f} "
            f"{trajectory.position[index, 0]:12.6f} "
            f"{exact:12.6f} "
            f"{trajectory.total_energy[index]:14.6f}"
        )

    initial_energy = trajectory.total_energy[0]
    maximum_drift = np.max(np.abs(trajectory.total_energy - initial_energy))
    print()
    print(f"Initial total energy: {initial_energy:.6f} kJ/mol")
    print(f"Maximum energy drift: {maximum_drift:.3e} kJ/mol")
    print()
    print("Velocity Verlet keeps the total energy nearly constant.")
    print("A smaller time step reduces the drift further.")


if __name__ == "__main__":
    main()
