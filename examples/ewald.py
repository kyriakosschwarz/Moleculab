"""Run Ewald summation for a small periodic ionic system.

The system is a neutral 2x2x2 lattice of +1 and -1 charges.  Lennard-Jones
terms are set to zero so that the total energy is purely electrostatic.

The example shows how the Ewald energy converges as the reciprocal-space
cutoff grows, and how the same model plugs into ``NonbondedModel``/``System``.

Run it with:

    uv run python examples/ewald.py
"""

import numpy as np

from moleculab import (
    Box,
    EwaldElectrostatics,
    LennardJonesParameters,
    NonbondedModel,
    Particle,
    System,
    max_force_magnitude,
)


def ionic_system() -> tuple[np.ndarray, np.ndarray, Box]:
    box_length = 4.0
    positions = []
    charges = []
    for i in range(2):
        for j in range(2):
            for k in range(2):
                positions.append(
                    [
                        (i + 0.5) * box_length / 2.0,
                        (j + 0.5) * box_length / 2.0,
                        (k + 0.5) * box_length / 2.0,
                    ]
                )
                charges.append(1.0 if (i + j + k) % 2 == 0 else -1.0)
    coordinates = np.asarray(positions)
    coordinates[1, 0] += 0.07  # break perfect symmetry so forces are visible
    return (
        coordinates,
        np.asarray(charges),
        Box([box_length, box_length, box_length]),
    )


def make_system(positions, charges, box, alpha, reciprocal_cutoff) -> System:
    ewald = EwaldElectrostatics(
        alpha=alpha,
        reciprocal_cutoff=reciprocal_cutoff,
    )
    model = NonbondedModel(
        parameters={
            "P": LennardJonesParameters(epsilon=0.0, sigma=0.3),
            "N": LennardJonesParameters(epsilon=0.0, sigma=0.3),
        },
        cutoff=float(np.min(box.lengths)) / 2.0,
        electrostatics=ewald,
    )
    particles = [
        Particle(
            mass=1.0,
            position=positions[index],
            atom_type="P" if charges[index] > 0 else "N",
            charge=charges[index],
        )
        for index in range(len(positions))
    ]
    return System(
        particles=particles,
        nonbonded_model=model,
        box=box,
    )


def main() -> None:
    positions, charges, box = ionic_system()

    print("Periodic ionic system for Ewald summation")
    print(f"  particles: {len(positions)}")
    print(f"  total charge: {np.sum(charges):+.1f} e")
    print(f"  box: {box.lengths} nm")
    print(f"  real-space cutoff: {np.min(box.lengths) / 2.0:.2f} nm")
    print()
    print("Convergence with reciprocal-space cutoff (alpha = 4.0 nm^-1)")
    print(f"{'k_max (nm^-1)':>14} {'Ewald energy (kJ/mol)':>24} {'max force':>14}")
    energies = []
    for reciprocal_cutoff in [10.0, 15.0, 20.0, 25.0, 30.0]:
        system = make_system(
            positions,
            charges,
            box,
            alpha=4.0,
            reciprocal_cutoff=reciprocal_cutoff,
        )
        energy = system.potential_energy()
        energies.append(energy)
        print(
            f"{reciprocal_cutoff:14.1f} {energy:24.6f} "
            f"{max_force_magnitude(system):14.6f}"
        )

    print()
    print("Convergence with splitting parameter (k_max = 30 nm^-1)")
    print(f"{'alpha (nm^-1)':>14} {'Ewald energy (kJ/mol)':>24}")
    for alpha in [3.0, 4.0, 5.0, 6.0]:
        system = make_system(
            positions,
            charges,
            box,
            alpha=alpha,
            reciprocal_cutoff=30.0,
        )
        print(
            f"{alpha:14.1f} {system.potential_energy():24.6f}"
        )

    print()
    print("The energy plateaus as more reciprocal vectors are included.")
    print("Ewald forces are added to the real-space pair forces by System.")
    print("Intramolecular exclusions are supported through a reciprocal-space")
    print("correction; see examples/water_ewald.py for bonded water.")


if __name__ == "__main__":
    main()
