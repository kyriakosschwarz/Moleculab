# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Build and simulate a minimal explicit three-site water box.

This example demonstrates the pieces that distinguish atomistic water from
the earlier coarse-grained solvent:

* per-atom types (``O`` and ``H``) and partial charges;
* two O-H bonds and one H-O-H angle per molecule;
* automatic 1-2 and 1-3 nonbonded exclusions;
* typed Lennard-Jones parameters plus reaction-field electrostatics.

The model is flexible and educational, not a production rigid-water model.
A conservative timestep is used because the O-H and angle terms are stiff.

Run it with:

    uv run python examples/advanced/water_box.py
"""

import numpy as np

from moleculab import (
    berendsen_thermostat,
    build_water_box,
    kinetic_temperature,
    max_force_magnitude,
    velocity_verlet_step_system,
)


def main() -> None:
    system = build_water_box(
        8,
        density=5.0,
        temperature=300.0,
        seed=7,
    )

    oxygen_count = sum(p.atom_type == "O" for p in system.particles)
    hydrogen_count = sum(p.atom_type == "H" for p in system.particles)
    excluded = system.topology.nonbonded_exclusions()

    print("Three-site water box")
    print(f"  molecules: {system.n_particles // 3}")
    print(f"  particles: {system.n_particles} (O={oxygen_count}, H={hydrogen_count})")
    print(f"  box: {system.box.lengths} nm")
    print(f"  bonds: {system.topology.n_bonds}")
    print(f"  angles: {system.topology.n_angles}")
    print(f"  nonbonded exclusions: {len(excluded)} intramolecular pairs")
    print(f"  reaction-field cutoff: {system.nonbonded_model.cutoff:.4f} nm")
    print(f"  reaction-field dielectric: 78.3")
    print()
    print(f"  initial temperature: {kinetic_temperature(system):.2f} K")
    print(f"  initial potential energy: {system.potential_energy():.4f} kJ/mol")
    print(f"  initial maximum force: {max_force_magnitude(system):.4f} kJ/(mol nm)")
    print()

    dt = 1e-4
    steps = 1000
    temperatures = []
    for _ in range(steps):
        system = velocity_verlet_step_system(system, dt)
        system = berendsen_thermostat(
            system,
            temperature=300.0,
            coupling_time=0.1,
            dt=dt,
        )
        temperatures.append(kinetic_temperature(system))

    print("Short NVT-like trajectory")
    print(f"  simulated time: {steps * dt:.3f} ps")
    print(f"  final temperature: {kinetic_temperature(system):.2f} K")
    print(f"  mean after warm-up: {np.mean(temperatures[steps // 2:]):.2f} K")
    print(f"  final potential energy: {system.potential_energy():.4f} kJ/mol")
    print(f"  final maximum force: {max_force_magnitude(system):.4f} kJ/(mol nm)")
    print()
    print("The nonbonded model excludes bonded O-H and 1-3 H-H pairs.")
    print("Reaction-field electrostatics handles charges inside the cutoff.")
    print("See examples/advanced/water_ewald.py for bonded water with Ewald")
    print("long-range electrostatics and the exclusion correction.")


if __name__ == "__main__":
    main()
