# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Minimal QM/MM: an electron in a 1D box coupled to MM point charges.

The finite-difference Schrödinger equation is solved for an electron in a
1 nm box.  Two fixed MM charges outside the box then perturb the energy
levels.  The unperturbed finite-difference energies are compared with the
analytic particle-in-a-box result.

Run it with:

    uv run python examples/advanced/qm_mm.py
"""

import numpy as np

from moleculab import ELECTRON_MASS_DA, MMPointCharge, ParticleInBox


def main() -> None:
    box = ParticleInBox(length=1.0, mass=ELECTRON_MASS_DA, n_grid=400)
    result = box.solve(n_states=3)
    exact = box.analytic_energies(3)

    print("Quantum particle in a 1D box")
    print(f"  mass: {box.mass:.8f} Da (electron)")
    print(f"  box length: {box.length:.2f} nm")
    print(f"  grid points: {box.n_grid}")
    print()
    print(f"{'state':>6} {'finite difference':>20} {'analytic':>14} {'error':>12}")
    for index, (energy, reference) in enumerate(zip(result.energies, exact), start=1):
        print(
            f"{index:6d} {energy:20.6f} {reference:14.6f} "
            f"{energy - reference:12.6e}"
        )

    charges = [
        MMPointCharge(position=-0.2, charge=0.5),
        MMPointCharge(position=1.2, charge=-0.5),
    ]
    qmmm = box.solve(
        mm_charges=charges,
        qm_charge=-1.0,
        n_states=3,
    )

    print()
    print("QM/MM perturbation from two fixed point charges")
    print(f"{'state':>6} {'QM/MM energy':>16} {'shift':>14}")
    for index, (energy, baseline) in enumerate(
        zip(qmmm.energies, result.energies),
        start=1,
    ):
        print(f"{index:6d} {energy:16.6f} {energy - baseline:14.6f}")

    print()
    print("MM charges at x = -0.20 nm (+0.5 e) and x = 1.20 nm (-0.5 e).")
    print("The QM particle is an electron (q = -1 e).")
    print("This is a 1D teaching model, not a production QM/MM engine.")


if __name__ == "__main__":
    main()
