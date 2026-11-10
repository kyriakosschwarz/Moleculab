"""Explore chemical conformation with a four-atom torsion.

This is a small, visible model of the conformational problem behind medicinal
chemistry: a molecule has several rotamers separated by barriers.  The example
does not attempt full chemistry perception; it shows how the torsion term in
the force field produces the energy landscape.

Run it with:

    uv run python examples/chemical_structure.py
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
    print("Minima are staggered; the eclipsed barriers are 2k.")


if __name__ == "__main__":
    main()
