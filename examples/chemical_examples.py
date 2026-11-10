"""Chemical examples beyond the harmonic oscillator.

The harmonic oscillator is the clean validation model.  This example shows
three chemically intuitive landscapes:

* a diatomic Lennard-Jones pair with a dissociation limit;
* a threefold torsional potential with rotamer minima;
* a double well with two stable states and a barrier.

Run it with:

    uv run python examples/chemical_examples.py
"""

import numpy as np

from moleculab import (
    DoubleWellPotential,
    LennardJonesPotential,
    Particle,
    PeriodicTorsionPotential,
    System,
    dihedral_angle,
    simulate_velocity_verlet,
    simulate_velocity_verlet_system,
)


def lennard_jones_pair() -> None:
    potential = LennardJonesPotential(epsilon=1.0, sigma=0.34)
    equilibrium = 2.0 ** (1.0 / 6.0) * 0.34

    print("Diatomic Lennard-Jones pair")
    print(f"  sigma:       {0.34:.3f} nm")
    print(f"  equilibrium: {equilibrium:.4f} nm")
    print(f"  well depth:  {potential.energy([equilibrium, 0.0, 0.0]):.4f} kJ/mol")
    print(f"  bond at 0.6 nm: {potential.energy([0.6, 0.0, 0.0]):.4f} kJ/mol")

    system = System(
        particles=[
            Particle(12.0, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]),
            Particle(12.0, [0.45, 0.0, 0.0], [0.5, 0.0, 0.0]),
        ],
        pair_potentials=[potential],
    )
    trajectory = simulate_velocity_verlet_system(system, dt=0.001, steps=4000)
    separation = np.linalg.norm(
        trajectory.position[:, 1, :] - trajectory.position[:, 0, :],
        axis=1,
    )
    print(f"  initial separation: {separation[0]:.3f} nm")
    print(f"  final separation:   {separation[-1]:.3f} nm")
    print("  A high enough kinetic energy lets the pair dissociate.")
    print()


def torsional_rotamers() -> None:
    potential = PeriodicTorsionPotential(barrier=2.0, periodicity=3)
    print("Threefold torsion: rotamers and barriers")
    for angle in [0.0, 60.0, 120.0, 180.0]:
        r_ij = np.array([0.15, 0.0, 0.0])
        r_jk = np.array([0.0, 0.15, 0.0])
        phi = np.deg2rad(angle)
        r_kl = 0.15 * np.array([-np.cos(phi), 0.0, np.sin(phi)])
        positions = [
            np.zeros(3),
            r_ij,
            r_ij + r_jk,
            r_ij + r_jk + r_kl,
        ]
        actual = np.rad2deg(dihedral_angle(*positions))
        energy = potential.energy(r_ij, r_jk, r_kl)
        print(f"  phi = {actual:5.1f} deg: U = {energy:6.3f} kJ/mol")
    print("  Minima are staggered; the eclipsed barrier is 2k.")
    print()


def double_well_transition() -> None:
    potential = DoubleWellPotential(barrier=1.0, width=0.25)
    particle = Particle(12.0, [0.25, 0.0, 0.0], [0.5, 0.0, 0.0])
    trajectory = simulate_velocity_verlet(particle, potential, dt=0.001, steps=8000)

    print("Double-well reaction coordinate")
    print("  minima at x = +/-0.25 nm; barrier at x = 0")
    print(f"  initial x: {trajectory.position[0, 0]:.3f} nm")
    print(f"  final x:   {trajectory.position[-1, 0]:.3f} nm")
    crossings = sum(
        trajectory.position[index - 1, 0] * trajectory.position[index, 0] < 0.0
        for index in range(1, len(trajectory.time))
    )
    print(f"  barrier crossings: {crossings}")
    print(f"  final energy: {trajectory.total_energy[-1]:.4f} kJ/mol")
    print(
        "  With enough energy the particle crosses the barrier and samples both "
        "wells."
    )
    print()


def main() -> None:
    lennard_jones_pair()
    torsional_rotamers()
    double_well_transition()


if __name__ == "__main__":
    main()
