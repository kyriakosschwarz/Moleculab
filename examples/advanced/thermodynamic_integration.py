"""Learn thermodynamic integration on an exactly solvable harmonic system.

The alchemical path linearly changes the stiffness of a three-dimensional
harmonic oscillator from k0 to k1:

    U_lambda = (1 - lambda) U0 + lambda U1.

Thermodynamic integration uses

    dG = integral_0^1 <dU/dlambda>_lambda dlambda.

The exact result is (3/2) R T ln(k1/k0).

Run it with:

    uv run python examples/advanced/thermodynamic_integration.py
"""

import numpy as np

from moleculab import (
    HarmonicPotential,
    harmonic_free_energy_difference,
    thermodynamic_integration_harmonic,
)


def main() -> None:
    temperature = 300.0
    potential0 = HarmonicPotential(stiffness=100.0, center=[0.0, 0.0, 0.0])
    potential1 = HarmonicPotential(stiffness=200.0, center=[0.0, 0.0, 0.0])

    exact = harmonic_free_energy_difference(potential0, potential1, temperature)
    print("Harmonic thermodynamic integration")
    print(f"  temperature: {temperature:.1f} K")
    print(f"  k0: {potential0.stiffness:.1f} kJ/(mol nm^2)")
    print(f"  k1: {potential1.stiffness:.1f} kJ/(mol nm^2)")
    print(f"  exact dG: {exact:.6f} kJ/mol")
    print()
    print(f"{'windows':>8} {'samples/window':>16} {'TI dG':>14} "
          f"{'std error':>12} {'error':>12}")
    for n_windows in [5, 11, 21, 41]:
        result = thermodynamic_integration_harmonic(
            potential0,
            potential1,
            temperature,
            n_windows=n_windows,
            n_samples=20_000,
            seed=2024,
            method="simpson",
        )
        print(
            f"{n_windows:8d} {20_000:16d} "
            f"{result.delta_g:14.6f} "
            f"{result.standard_error:12.6f} "
            f"{result.delta_g - exact:12.6f}"
        )

    print()
    print("Selected <dU/dlambda> values at 21 windows")
    result = thermodynamic_integration_harmonic(
        potential0,
        potential1,
        temperature,
        n_windows=21,
        n_samples=20_000,
        seed=2024,
        method="simpson",
    )
    for index in [0, 5, 10, 15, 20]:
        print(
            f"  lambda={result.lambda_values[index]:.2f}  "
            f"<dU/dlambda>={result.derivatives[index]:.4f}  "
            f"+/- {result.standard_errors[index]:.4f}"
        )

    print()
    print("Simpson's rule converges quickly for a smooth alchemical path.")
    print("For molecular transformations, many lambda windows are needed")
    print("because the derivative can change sharply.")


if __name__ == "__main__":
    main()
