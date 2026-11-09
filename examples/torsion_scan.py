"""Scan the energy of a four-atom torsion.

The dihedral angle controls molecular conformation.  This example rotates the
last atom around the central bond and prints the resulting dihedral angle and
periodic torsion energy.

Run it with:

    uv run python examples/torsion_scan.py
"""

import numpy as np

from moleculab import PeriodicTorsionPotential, dihedral_angle


def main() -> None:
    # A chain i-j-k-l.  Rotating atom l around the j-k bond changes phi.
    r_i = np.array([0.0, 0.0, 0.0])
    r_j = np.array([0.15, 0.0, 0.0])
    r_k = np.array([0.15, 0.15, 0.0])
    radius = 0.15

    potential = PeriodicTorsionPotential(barrier=2.0, periodicity=3, phase=0.0)

    print("Periodic torsion: U = 2.0 * (1 + cos(3*phi)) kJ/mol")
    print()
    print(f"{'rotation (deg)':>15} {'phi (deg)':>12} {'U (kJ/mol)':>14}")
    for rotation_degrees in range(0, 361, 30):
        angle = np.deg2rad(rotation_degrees)
        r_l = r_k + np.array(
            [radius * np.cos(angle), 0.0, radius * np.sin(angle)]
        )

        phi = dihedral_angle(r_i, r_j, r_k, r_l)
        r_ij = r_j - r_i
        r_jk = r_k - r_j
        r_kl = r_l - r_k
        energy = potential.energy(r_ij, r_jk, r_kl)
        print(
            f"{rotation_degrees:15.0f} {np.rad2deg(phi):12.2f} {energy:14.6f}"
        )

    print()
    print("The energy repeats every 120 degrees because periodicity = 3.")


if __name__ == "__main__":
    main()
