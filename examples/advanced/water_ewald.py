# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Run bonded rigid water with Ewald long-range electrostatics.

This example builds two rigid three-site water molecules in a periodic box and
uses Ewald summation for the electrostatic interactions.  The 1-2/1-3
intramolecular pairs are excluded from the real-space pair loop, and the Ewald
module subtracts their reciprocal-space contributions with an exclusion
correction.

Run it with:

    uv run python examples/advanced/water_ewald.py
"""

import numpy as np

from moleculab import (
    EwaldElectrostatics,
    berendsen_thermostat,
    build_water_box,
    kinetic_temperature,
    velocity_verlet_step_system,
)


def maximum_constraint_error(positions, constraints):
    return float(
        max(
            abs(
                np.linalg.norm(positions[c.i] - positions[c.j]) - c.distance
            )
            for c in constraints
        )
    )


def main() -> None:
    ewald = EwaldElectrostatics(alpha=6.0, reciprocal_cutoff=40.0)
    system = build_water_box(
        2,
        box_length=3.0,
        temperature=300.0,
        seed=11,
        electrostatics=ewald,
    )
    wave_vectors, _ = ewald.reciprocal_k_vectors(system.box)

    print("Rigid water with Ewald electrostatics")
    print(f"  particles: {system.n_particles}")
    print(f"  constraints: {len(system.constraints)}")
    print(f"  box: {system.box.lengths} nm")
    print(f"  Ewald alpha: {ewald.alpha:.2f} nm^-1")
    print(f"  reciprocal cutoff: {ewald.reciprocal_cutoff:.2f} nm^-1")
    print(f"  reciprocal vectors: {len(wave_vectors)}")
    print(f"  initial temperature: {kinetic_temperature(system):.2f} K")
    print(f"  initial potential energy: {system.potential_energy():.6f} kJ/mol")
    print()

    dt = 5e-5
    steps = 500
    for _ in range(steps):
        system = velocity_verlet_step_system(system, dt)
        system = berendsen_thermostat(
            system,
            temperature=300.0,
            coupling_time=0.02,
            dt=dt,
        )

    print("Short constrained Ewald trajectory")
    print(f"  simulated time: {steps * dt:.4f} ps")
    print(f"  final temperature: {kinetic_temperature(system):.2f} K")
    print(f"  final potential energy: {system.potential_energy():.6f} kJ/mol")
    print(
        "  maximum constraint error: "
        f"{maximum_constraint_error(system.positions, system.constraints):.3e} nm"
    )
    print()
    print("Bonded water can now use long-range electrostatics.  The exclusion")
    print("correction removes the reciprocal part of each excluded intramolecular")
    print("pair; self-image interactions remain part of the periodic Ewald energy.")


if __name__ == "__main__":
    main()
