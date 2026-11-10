# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Run a small NVT-like simulation with the Berendsen thermostat.

The two particles form an argon-like Lennard-Jones pair.  The loop alternates
between a velocity Verlet step and a thermostat update.  The thermostat does
not change the forces; it only rescales velocities to move the instantaneous
kinetic temperature toward the target.

Run it with:

    uv run python examples/physics/thermostat.py
"""

import numpy as np

from moleculab import (
    LennardJonesPotential,
    Particle,
    System,
    berendsen_thermostat,
    kinetic_temperature,
    velocity_verlet_step_system,
)


def main() -> None:
    sigma = 0.34
    epsilon = 0.997
    equilibrium = 2.0 ** (1.0 / 6.0) * sigma

    system = System(
        particles=[
            Particle(
                39.948,
                [-equilibrium / 2.0, 0.0, 0.0],
                [0.0, 0.05, 0.0],
            ),
            Particle(
                39.948,
                [equilibrium / 2.0, 0.0, 0.0],
                [0.0, -0.05, 0.0],
            ),
        ],
        pair_potentials=[LennardJonesPotential(epsilon=epsilon, sigma=sigma)],
    )

    target_temperature = 100.0
    coupling_time = 0.1
    dt = 0.001
    steps = 2000

    print("Berendsen thermostat")
    print(f"  target temperature: {target_temperature:.1f} K")
    print(f"  coupling time: {coupling_time:.2f} ps")
    print(f"  initial temperature: {kinetic_temperature(system):.2f} K")
    print()

    temperatures = []
    for _ in range(steps):
        system = velocity_verlet_step_system(system, dt)
        system = berendsen_thermostat(
            system,
            temperature=target_temperature,
            coupling_time=coupling_time,
            dt=dt,
        )
        temperatures.append(kinetic_temperature(system))

    print(f"  final temperature: {kinetic_temperature(system):.2f} K")
    print(f"  mean after warm-up: {np.mean(temperatures[steps // 2:]):.2f} K")
    print()
    print("The thermostat nudges the system toward the target temperature.")
    print("It is not part of the physical potential; it is a sampling algorithm.")


if __name__ == "__main__":
    main()
