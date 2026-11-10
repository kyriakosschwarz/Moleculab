"""Use MBAR to combine samples from several alchemical harmonic states.

The alchemical path changes the harmonic stiffness through k = 100, 150, 200
kJ/(mol nm^2).  The exact relative free energies are

    dG_i = (3/2) R T ln(k_i / k_0).

MBAR estimates all states simultaneously from samples drawn from every state.

Run it with:

    uv run python examples/advanced/mbar.py
"""

import numpy as np

from moleculab import (
    HarmonicPotential,
    harmonic_free_energy_difference,
    mbar_harmonic,
)
from moleculab.foundation.units import BOLTZMANN_CONSTANT


def main() -> None:
    temperature = 300.0
    stiffnesses = [100.0, 150.0, 200.0]
    potentials = [
        HarmonicPotential(stiffness=k, center=[0.0, 0.0, 0.0])
        for k in stiffnesses
    ]
    rt = BOLTZMANN_CONSTANT * temperature
    exact = np.array(
        [1.5 * rt * np.log(k / stiffnesses[0]) for k in stiffnesses]
    )

    print("MBAR harmonic alchemical path")
    print(f"  temperature: {temperature:.1f} K")
    print(f"  stiffnesses: {stiffnesses}")
    print("  exact dG:", np.round(exact, 6))
    print()

    print(f"{'samples/state':>14} {'iterations':>12} "
          f"{'dG endpoint':>14} {'error':>12}")
    for n_samples in [1_000, 5_000, 20_000, 50_000]:
        result = mbar_harmonic(
            potentials,
            temperature,
            n_samples=n_samples,
            seed=2024,
        )
        endpoint_error = result.delta_g[-1] - exact[-1]
        print(
            f"{n_samples:14d} {result.iterations:12d} "
            f"{result.delta_g[-1]:14.6f} {endpoint_error:12.6f}"
        )

    result = mbar_harmonic(
        potentials,
        temperature,
        n_samples=20_000,
        seed=2024,
    )
    print()
    print("Final MBAR free energies")
    print(f"{'state':>6} {'k':>8} {'MBAR dG':>14} {'exact dG':>14} {'error':>12}")
    for index, (value, reference) in enumerate(zip(result.delta_g, exact)):
        print(
            f"{index:6d} {stiffnesses[index]:8.1f} "
            f"{value:14.6f} {reference:14.6f} {value - reference:12.6f}"
        )

    print()
    print("MBAR uses all windows together, so it remains accurate even when")
    print("neighbouring states have limited direct overlap.")


if __name__ == "__main__":
    main()
