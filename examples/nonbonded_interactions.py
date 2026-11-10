"""Compare Lennard-Jones and Coulomb interactions at several distances.

This example makes the shape of the two nonbonded potentials visible.  The
Lennard-Jones potential has a repulsive wall at short range, a minimum at
r = 2^(1/6) sigma, and a weak attractive tail.  The Coulomb potential is
positive for like charges and negative for opposite charges.

Run it with:

    uv run python examples/nonbonded_interactions.py
"""

import numpy as np

from moleculab import CoulombPotential, LennardJonesPotential
from moleculab.foundation.units import COULOMB_CONSTANT


def main() -> None:
    # Argon-like parameters: epsilon ~ 1 kJ/mol, sigma ~ 0.34 nm.
    lj = LennardJonesPotential(epsilon=0.997, sigma=0.34)
    like = CoulombPotential(charge_i=1.0, charge_j=1.0)
    opposite = CoulombPotential(charge_i=1.0, charge_j=-1.0)

    print("Lennard-Jones: epsilon = 0.997 kJ/mol, sigma = 0.34 nm")
    print("Coulomb constant:", f"{COULOMB_CONSTANT:.6f} kJ nm/(mol e^2)")
    print()
    print(
        f"{'r (nm)':>8} {'LJ U':>12} {'LJ F':>14} "
        f"{'q=+1 U':>12} {'q=+1 F':>14} {'q=-1 U':>12}"
    )
    for distance in [0.30, 0.34, 0.3816, 0.45, 0.60, 0.80, 1.00]:
        displacement = np.array([distance, 0.0, 0.0])
        lj_energy = lj.energy(displacement)
        lj_force = lj.force(displacement)[0]
        like_energy = like.energy(displacement)
        like_force = like.force(displacement)[0]
        opposite_energy = opposite.energy(displacement)
        print(
            f"{distance:8.4f} {lj_energy:12.4f} {lj_force:14.4f} "
            f"{like_energy:12.4f} {like_force:14.4f} {opposite_energy:12.4f}"
        )

    print()
    print("The LJ force is zero at the well minimum r = 2^(1/6) sigma.")
    print("Positive force on particle i means repulsion along +x.")


if __name__ == "__main__":
    main()
