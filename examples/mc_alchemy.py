"""Multi-particle Monte Carlo FEP, TI, and MBAR for a bonded LJ pair.

Two particles are connected by a harmonic bond and interact through a
soft-core Lennard-Jones pair potential.  Metropolis Monte Carlo moves both
particles while preserving their center of mass.  The resulting samples feed
FEP, TI, and MBAR estimators, which are compared with an independent radial
quadrature reference.

Run it with:

    uv run python examples/mc_alchemy.py
"""

from moleculab import BondedLennardJonesMC


def main() -> None:
    pair = BondedLennardJonesMC()
    reference = pair.reference_delta_g()

    print("Bonded two-particle Monte Carlo alchemy")
    print(f"  epsilon: {pair.epsilon:.3f} kJ/mol")
    print(f"  sigma: {pair.sigma:.3f} nm")
    print(f"  bond length: {pair.bond_length:.3f} nm")
    print(f"  bond stiffness: {pair.bond_stiffness:.1f} kJ/(mol nm^2)")
    print(f"  temperature: {pair.temperature:.1f} K")
    print(f"  reference dG: {reference:.6f} kJ/mol")
    print()

    print(f"{'windows':>8} {'FEP dG':>12} {'TI dG':>12} {'MBAR dG':>12} "
          f"{'FEP err':>10} {'TI err':>10} {'MBAR err':>10}")
    for n_windows in [5, 11, 21]:
        fep = pair.fep(n_windows=n_windows, n_samples=1000, seed=100 + n_windows)
        ti = pair.ti(n_windows=n_windows, n_samples=1000, seed=200 + n_windows)
        mbar = pair.mbar(n_windows=n_windows, n_samples=1000, seed=300 + n_windows)
        print(
            f"{n_windows:8d} {fep.delta_g:12.6f} {ti.delta_g:12.6f} "
            f"{mbar.delta_g[-1]:12.6f} "
            f"{fep.delta_g - reference:10.6f} "
            f"{ti.delta_g - reference:10.6f} "
            f"{mbar.delta_g[-1] - reference:10.6f}"
        )

    print()
    print("Sampling uses actual Metropolis moves on a two-particle System.")
    print("The soft core keeps the LJ endpoint finite; the reference is")
    print("an independent radial quadrature of the same relative-coordinate")
    print("partition function.")


if __name__ == "__main__":
    main()
