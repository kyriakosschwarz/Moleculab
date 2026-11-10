# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Learn free-energy perturbation on an exactly solvable harmonic system.

The transformation changes the stiffness of a three-dimensional harmonic
oscillator from k0 to k1.  Its exact free-energy difference is

    dG_exact = (3/2) R T ln(k1 / k0).

FEP estimates the same quantity from samples of the k0 state:

    dG = -R T ln < exp(-(U1 - U0)/(R T)) >_0.

Run it with:

    uv run python examples/advanced/free_energy_perturbation.py
"""

from moleculab import (
    HarmonicPotential,
    fep_harmonic_oscillator,
    harmonic_free_energy_difference,
)


def main() -> None:
    temperature = 300.0
    potential0 = HarmonicPotential(stiffness=100.0, center=[0.0, 0.0, 0.0])
    potential1 = HarmonicPotential(stiffness=200.0, center=[0.0, 0.0, 0.0])

    exact = harmonic_free_energy_difference(potential0, potential1, temperature)
    print("Harmonic free-energy perturbation")
    print(f"  temperature: {temperature:.1f} K")
    print(f"  k0: {potential0.stiffness:.1f} kJ/(mol nm^2)")
    print(f"  k1: {potential1.stiffness:.1f} kJ/(mol nm^2)")
    print(f"  exact dG: {exact:.6f} kJ/mol")
    print()
    print("FEP estimate = -RT ln <exp(-dU/RT)> over samples from state 0")
    print()
    print(f"{'samples':>10} {'dG estimate':>14} {'std error':>12} {'error':>12}")
    for n_samples in [1_000, 10_000, 100_000, 1_000_000]:
        result = fep_harmonic_oscillator(
            potential0,
            potential1,
            temperature,
            n_samples=n_samples,
            seed=2024,
        )
        error = result.delta_g - exact
        print(
            f"{n_samples:10d} "
            f"{result.delta_g:14.6f} "
            f"{result.standard_error:12.6f} "
            f"{error:12.6f}"
        )

    print()
    print("More samples reduce the statistical error, but FEP variance also")
    print("depends on how different the two states are.  Large perturbations")
    print("need better overlap or bidirectional methods such as MBAR.")


if __name__ == "__main__":
    main()
