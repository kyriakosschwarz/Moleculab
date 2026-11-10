# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Molecular FEP, TI, and MBAR on a soft-core Lennard-Jones pair.

A mobile particle is held near a fixed particle by a harmonic trap.  The
Lennard-Jones interaction is coupled through a soft-core path so that
lambda=0 is the pure trap and lambda=1 is the full LJ pair.

The reference free energy is computed by high-accuracy radial quadrature.
FEP, TI, and MBAR are then compared with that reference.

Run it with:

    uv run python examples/advanced/lj_alchemy.py
"""

from moleculab import SoftCoreLennardJones


def main() -> None:
    pair = SoftCoreLennardJones()
    reference = pair.reference_delta_g()

    print("Soft-core Lennard-Jones alchemical pair")
    print(f"  epsilon: {pair.epsilon:.3f} kJ/mol")
    print(f"  sigma: {pair.sigma:.3f} nm")
    print(f"  trap stiffness: {pair.trap_stiffness:.2f} kJ/(mol nm^2)")
    print(f"  temperature: {pair.temperature:.1f} K")
    print(f"  soft-core alpha: {pair.alpha:.2f}")
    print(f"  reference dG: {reference:.6f} kJ/mol")
    print()

    print("Convergence with number of windows")
    print(f"{'windows':>8} {'FEP dG':>12} {'TI dG':>12} {'MBAR dG':>12} "
          f"{'FEP err':>10} {'TI err':>10} {'MBAR err':>10}")
    for n_windows in [5, 11, 21]:
        fep = pair.fep(n_windows=n_windows, n_samples=20_000, seed=100 + n_windows)
        ti = pair.ti(n_windows=n_windows, n_samples=20_000, seed=200 + n_windows)
        mbar = pair.mbar(n_windows=n_windows, n_samples=20_000, seed=300 + n_windows)
        print(
            f"{n_windows:8d} {fep.delta_g:12.6f} {ti.delta_g:12.6f} "
            f"{mbar.delta_g[-1]:12.6f} "
            f"{fep.delta_g - reference:10.6f} "
            f"{ti.delta_g - reference:10.6f} "
            f"{mbar.delta_g[-1] - reference:10.6f}"
        )

    print()
    print("FEP uses adjacent-window work values; TI integrates <dU/dlambda>;")
    print("MBAR combines all windows self-consistently.")
    print("The soft core keeps dU/dlambda finite at lambda=0.")


if __name__ == "__main__":
    main()
