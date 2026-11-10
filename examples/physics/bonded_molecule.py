# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Simulate a small bonded molecule with bond and angle terms.

This is a water-like triatomic model, but the force-field parameters are
deliberately simple rather than a published water model.  The example shows
how topology turns a list of particles into a molecule:

* two O-H bonds,
* one H-O-H angle,
* and no nonbonded terms yet.

Run it with:

    uv run python examples/physics/bonded_molecule.py
"""

import numpy as np

from moleculab import (
    Angle,
    Bond,
    HarmonicAnglePotential,
    HarmonicBondPotential,
    Particle,
    System,
    Topology,
    simulate_velocity_verlet_system,
)


def angle_degrees(positions, i, j, k):
    r_ji = positions[i] - positions[j]
    r_jk = positions[k] - positions[j]
    cosine = np.dot(r_ji, r_jk) / (np.linalg.norm(r_ji) * np.linalg.norm(r_jk))
    return float(np.rad2deg(np.arccos(np.clip(cosine, -1.0, 1.0))))


def main() -> None:
    # Oxygen at the origin, two hydrogens in the xy plane.
    particles = [
        Particle(15.999, [0.0, 0.0, 0.0]),
        Particle(1.008, [0.11, 0.0, 0.0]),
        Particle(1.008, [-0.03, 0.10, 0.0]),
    ]

    topology = Topology(
        bonds=[
            Bond(0, 1, HarmonicBondPotential(stiffness=3000.0, length0=0.10)),
            Bond(0, 2, HarmonicBondPotential(stiffness=3000.0, length0=0.10)),
        ],
        angles=[
            Angle(
                1,
                0,
                2,
                HarmonicAnglePotential(
                    stiffness=500.0,
                    angle0=np.deg2rad(104.5),
                ),
            )
        ],
    )
    system = System(particles=particles, topology=topology)

    print("Initial geometry")
    print(f"  O-H1 bond: {np.linalg.norm(particles[1].position - particles[0].position):.4f} nm")
    print(f"  O-H2 bond: {np.linalg.norm(particles[2].position - particles[0].position):.4f} nm")
    print(f"  H-O-H angle: {angle_degrees(system.positions, 1, 0, 2):.2f} degrees")
    print(f"  bonded energy: {system.bonded_energy():.6f} kJ/mol")
    print()

    dt = 0.0005
    trajectory = simulate_velocity_verlet_system(system, dt=dt, steps=2000)

    initial_energy = trajectory.total_energy[0]
    maximum_drift = np.max(np.abs(trajectory.total_energy - initial_energy))
    final_positions = trajectory.position[-1]

    print(f"Simulated {len(trajectory) - 1} steps of {dt} ps")
    print(f"  final O-H1 bond: {np.linalg.norm(final_positions[1] - final_positions[0]):.4f} nm")
    print(f"  final O-H2 bond: {np.linalg.norm(final_positions[2] - final_positions[0]):.4f} nm")
    print(f"  final H-O-H angle: {angle_degrees(final_positions, 1, 0, 2):.2f} degrees")
    print(f"  initial total energy: {initial_energy:.6f} kJ/mol")
    print(f"  maximum energy drift: {maximum_drift:.3e} kJ/mol")
    print()
    print("The bonds stretch and the angle bends, but total energy is conserved.")


if __name__ == "__main__":
    main()
