"""Explore chemical conformation with a four-atom torsion.

This is a small, visible model of the conformational problem behind medicinal
chemistry: a molecule has several rotamers separated by barriers.  The example
does not attempt full chemistry perception; it shows how the torsion term in
the force field produces the energy landscape.

Run it with:

    uv run python examples/chemistry/chemical_structure.py
"""

import numpy as np

from moleculab import (
    Angle,
    Bond,
    HarmonicAnglePotential,
    HarmonicBondPotential,
    Particle,
    PeriodicTorsionPotential,
    System,
    Topology,
    Torsion,
    dihedral_angle,
    minimize_energy,
)


def build_chain() -> tuple[System, PeriodicTorsionPotential]:
    particles = [
        Particle(12.0, [-0.30, 0.00, 0.00]),
        Particle(12.0, [0.00, 0.00, 0.00]),
        Particle(12.0, [0.25, 0.20, 0.00]),
        Particle(12.0, [0.25, 0.20, 0.25]),
    ]
    torsion = PeriodicTorsionPotential(barrier=2.0, periodicity=3)
    topology = Topology(
        bonds=[
            Bond(0, 1, HarmonicBondPotential(stiffness=1000.0, length0=0.15)),
            Bond(1, 2, HarmonicBondPotential(stiffness=1000.0, length0=0.15)),
            Bond(2, 3, HarmonicBondPotential(stiffness=1000.0, length0=0.15)),
        ],
        angles=[
            Angle(
                0,
                1,
                2,
                HarmonicAnglePotential(stiffness=100.0, angle0=np.deg2rad(109.5)),
            ),
            Angle(
                1,
                2,
                3,
                HarmonicAnglePotential(stiffness=100.0, angle0=np.deg2rad(109.5)),
            ),
        ],
        torsions=[Torsion(0, 1, 2, 3, torsion)],
    )
    return System(particles=particles, topology=topology), torsion


def main() -> None:
    system, torsion = build_chain()
    initial_dihedral = np.rad2deg(
        dihedral_angle(*[particle.position for particle in system.particles])
    )
    print("Four-atom torsion model")
    print(f"  initial dihedral: {initial_dihedral:7.2f} deg")
    print(f"  initial energy:   {system.potential_energy():7.3f} kJ/mol")

    result = minimize_energy(system, step_size=1e-5, max_steps=5000)
    final_dihedral = np.rad2deg(
        dihedral_angle(
            *[particle.position for particle in result.system.particles]
        )
    )
    print(f"  final dihedral:   {final_dihedral:7.2f} deg")
    print(f"  final energy:     {result.system.potential_energy():7.3f} kJ/mol")
    print(f"  accepted steps:   {result.steps}")
    print()

    print("Torsion energy at selected dihedrals")
    for angle in [0.0, 60.0, 120.0, 180.0]:
        energy = torsion.energy_from_angle(np.deg2rad(angle))
        print(f"  phi = {angle:5.1f} deg: U = {energy:5.2f} kJ/mol")
    print()

    angles = np.linspace(0.0, 360.0, 361)
    energies = np.array(
        [torsion.energy_from_angle(np.deg2rad(angle)) for angle in angles]
    )
    print("Full rotational scan")
    print(f"  minimum: {np.min(energies):5.2f} kJ/mol")
    print(f"  maximum: {np.max(energies):5.2f} kJ/mol")
    print(f"  barrier: {np.max(energies) - np.min(energies):5.2f} kJ/mol")
    print(f"  U(0) = U(360): {energies[0]:.2f} = {energies[-1]:.2f} kJ/mol")
    print()

    print("Conformational entropy: anti/gauche toy model at 300 K")
    state_energies = np.array([0.0, 0.9, 0.9])
    rt = 0.008314462618 * 300.0
    weights = np.exp(-state_energies / rt)
    probabilities = weights / np.sum(weights)
    entropy_j_per_mol_k = -0.008314462618 * 1000.0 * np.sum(
        probabilities * np.log(probabilities)
    )
    for index, (energy, probability) in enumerate(
        zip(state_energies, probabilities)
    ):
        print(f"  state {index}: E = {energy:4.1f} kJ/mol, p = {probability:.3f}")
    print(f"  conformational entropy = {entropy_j_per_mol_k:.2f} J/(mol K)")
    print()
    print("Minima are staggered; the eclipsed barriers are 2k.")


if __name__ == "__main__":
    main()
