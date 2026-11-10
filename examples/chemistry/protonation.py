"""Protonation and tautomer populations from acid-base thermodynamics.

This example uses pKa values as inputs. It does not predict pKa. Run it with:

    uv run python examples/chemistry/protonation.py
"""

from __future__ import annotations

import numpy as np

from moleculab import (
    boltzmann_populations,
    deprotonated_fraction,
    mean_net_charge,
    mean_protonation,
    protonated_fraction,
    protonation_microstates,
    tautomer_populations,
    units,
)


def main() -> None:
    ph = 7.4
    print(f"Protonation fractions at pH {ph}")
    print(f"{'species':>14} {'pKa':>5} {'deprotonated':>13} {'protonated':>11}")
    for name, pka in [
        ("acetic acid", 4.76),
        ("ammonium", 9.25),
        ("imidazole", 6.00),
    ]:
        print(
            f"{name:>14} {pka:5.2f} "
            f"{deprotonated_fraction(pka, ph):13.4f} "
            f"{protonated_fraction(pka, ph):11.4f}"
        )
    print()

    print("Acetic acid titration (pKa 4.76)")
    print(f"{'pH':>5} {'f(A-)':>8} {'f(HA)':>8}")
    for value in [2.0, 3.76, 4.76, 5.76, 7.4]:
        print(
            f"{value:5.2f} {deprotonated_fraction(4.76, value):8.4f} "
            f"{protonated_fraction(4.76, value):8.4f}"
        )
    print()

    pkas = [4.0, 9.0]
    protonated_charges = [0.0, 1.0]
    deprotonated_charges = [-1.0, 0.0]
    microstates = protonation_microstates(pkas, ph=7.0)
    print("Two independent ionizable sites at pH 7.0")
    print(f"{'state':>7} {'probability':>12} {'charge':>7}")
    for state, probability in enumerate(microstates):
        label = "".join(
            "H" if (state >> site) & 1 else "D" for site in range(len(pkas))
        )
        charge = sum(
            (
                protonated_charges[site]
                if (state >> site) & 1
                else deprotonated_charges[site]
            )
            for site in range(len(pkas))
        )
        print(f"{label:>7} {probability:12.6f} {charge:7.1f}")

    net_charge = mean_net_charge(
        pkas,
        ph=7.0,
        protonated_charges=protonated_charges,
        deprotonated_charges=deprotonated_charges,
    )
    site_fractions = mean_protonation(pkas, ph=7.0)
    print(
        "  mean protonation: "
        + ", ".join(f"{value:.4f}" for value in site_fractions)
    )
    print(f"  mean net charge:  {net_charge:.4f}")
    print()

    temperature = 298.15
    rt = units.MOLAR_GAS_CONSTANT * temperature
    print(f"Tautomer populations at {temperature:.2f} K (RT = {rt:.4f} kJ/mol)")
    print(f"{'delta G':>8} {'p(A)':>8} {'p(B)':>8}")
    for delta_g in [0.0, -2.0, 2.0, 5.0]:
        populations = tautomer_populations(delta_g, temperature)
        print(f"{delta_g:8.1f} {populations[0]:8.4f} {populations[1]:8.4f}")
    print()

    print("Coupled two-site microstates at 298.15 K")
    coupled_pkas = [4.0, 9.0]
    coupled_ph = 7.0
    coupling = 2.0
    rt_ln10 = rt * np.log(10.0)
    state_energies = []
    for state in range(1 << 2):
        energy = sum(
            (coupled_ph - coupled_pkas[site]) * rt_ln10
            for site in range(2)
            if (state >> site) & 1
        )
        if state == 3:
            energy += coupling
        state_energies.append(energy)
    populations = boltzmann_populations(state_energies, temperature)
    print("  state energies include a +2.0 kJ/mol HH coupling")
    print(f"{'state':>7} {'energy':>8} {'probability':>12}")
    for state, (energy, probability) in enumerate(
        zip(state_energies, populations)
    ):
        label = "".join("H" if (state >> site) & 1 else "D" for site in range(2))
        print(f"{label:>7} {energy:8.2f} {probability:12.6f}")
    print()
    print("pH changes populations of discrete states; it does not scale")
    print("the partial charges of one fixed molecular graph.")


if __name__ == "__main__":
    main()
