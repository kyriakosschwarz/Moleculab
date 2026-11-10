# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Minimize a distorted triatomic molecule.

Energy minimization is the molecular-mechanics counterpart of letting a
physical model settle into its nearest stable structure.  At the minimum the
net force on every atom is zero.

Run it with:

    uv run python examples/physics/minimization.py
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
    max_force_magnitude,
    minimize_energy,
)


def angle_degrees(positions, i, j, k):
    r_ji = positions[i] - positions[j]
    r_jk = positions[k] - positions[j]
    cosine = np.dot(r_ji, r_jk) / (np.linalg.norm(r_ji) * np.linalg.norm(r_jk))
    return float(np.rad2deg(np.arccos(np.clip(cosine, -1.0, 1.0))))


def main() -> None:
    particles = [
        Particle(12.0, [1.20, 0.0, 0.0]),
        Particle(12.0, [0.0, 0.0, 0.0]),
        Particle(12.0, [0.0, 1.20, 0.0]),
    ]
    topology = Topology(
        bonds=[
            Bond(0, 1, HarmonicBondPotential(stiffness=1000.0, length0=1.0)),
            Bond(1, 2, HarmonicBondPotential(stiffness=1000.0, length0=1.0)),
        ],
        angles=[
            Angle(
                0,
                1,
                2,
                HarmonicAnglePotential(stiffness=100.0, angle0=np.pi / 3.0),
            )
        ],
    )
    system = System(particles=particles, topology=topology)

    print("Before minimization")
    print(f"  bond 0-1: {np.linalg.norm(system.positions[0] - system.positions[1]):.4f} nm")
    print(f"  bond 1-2: {np.linalg.norm(system.positions[1] - system.positions[2]):.4f} nm")
    print(f"  angle: {angle_degrees(system.positions, 0, 1, 2):.2f} degrees")
    print(f"  energy: {system.potential_energy():.6f} kJ/mol")
    print(f"  maximum force: {max_force_magnitude(system):.6f} kJ/(mol nm)")
    print()

    result = minimize_energy(system, step_size=1e-5, max_steps=5000)

    final = result.system
    print("After minimization")
    print(f"  bond 0-1: {np.linalg.norm(final.positions[0] - final.positions[1]):.4f} nm")
    print(f"  bond 1-2: {np.linalg.norm(final.positions[1] - final.positions[2]):.4f} nm")
    print(f"  angle: {angle_degrees(final.positions, 0, 1, 2):.2f} degrees")
    print(f"  energy: {final.potential_energy():.6f} kJ/mol")
    print(f"  maximum force: {max_force_magnitude(final):.6f} kJ/(mol nm)")
    print(f"  accepted steps: {result.steps}")
    print(f"  converged: {result.converged}")
    print()
    print("The bonds reach 1.0 nm and the angle reaches the target 60 degrees.")


if __name__ == "__main__":
    main()
