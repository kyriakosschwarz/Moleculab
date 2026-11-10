"""Illustrate why contact counts do not determine binding free energy.

The numbers are hypothetical components, not a force field or a free-energy
calculation. The point is the bookkeeping: interaction, desolvation, and
entropy terms can reverse a contact-count ranking.

Run:

    uv run python examples/drug_discovery/structural_observations_thermodynamics.py
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from moleculab import boltzmann_populations, units


@dataclass(frozen=True)
class ToyBindingComponents:
    """Hypothetical binding free-energy components in kJ/mol."""

    label: str
    persistent_contacts: int
    interaction: float
    desolvation: float
    ligand_entropy: float
    protein_entropy: float
    water_reorganisation: float

    @property
    def total(self) -> float:
        return (
            self.interaction
            + self.desolvation
            + self.ligand_entropy
            + self.protein_entropy
            + self.water_reorganisation
        )


def main() -> None:
    ligand_a = ToyBindingComponents(
        label="A",
        persistent_contacts=4,
        interaction=-45.0,
        desolvation=20.0,
        ligand_entropy=12.0,
        protein_entropy=8.0,
        water_reorganisation=5.0,
    )
    ligand_b = ToyBindingComponents(
        label="B",
        persistent_contacts=2,
        interaction=-30.0,
        desolvation=8.0,
        ligand_entropy=4.0,
        protein_entropy=3.0,
        water_reorganisation=2.0,
    )

    print("Hypothetical binding components (kJ/mol)")
    print(
        f"{'ligand':>6} {'contacts':>8} {'interaction':>11} {'desolvation':>11} "
        f"{'lig. entropy':>12} {'prot. entropy':>13} {'water':>7} {'total':>7}"
    )
    for ligand in [ligand_a, ligand_b]:
        print(
            f"{ligand.label:>6} {ligand.persistent_contacts:8d} "
            f"{ligand.interaction:11.1f} {ligand.desolvation:11.1f} "
            f"{ligand.ligand_entropy:12.1f} {ligand.protein_entropy:13.1f} "
            f"{ligand.water_reorganisation:7.1f} {ligand.total:7.1f}"
        )
    print()
    print("Contact count alone would rank A above B.")
    print("The hypothetical free-energy components rank B above A.")
    print("This is the point: contact count is not a thermodynamic ranking.")
    print()

    temperature = 298.15
    rt = units.MOLAR_GAS_CONSTANT * temperature
    # State A: energy 0, degeneracy 1. State B: energy 5 kJ/mol, degeneracy 100.
    effective_energies = np.array([0.0, 5.0 - rt * np.log(100.0)])
    populations = boltzmann_populations(effective_energies, temperature)
    print("Degeneracy and population")
    print("  state A: E = 0 kJ/mol, degeneracy 1")
    print("  state B: E = 5 kJ/mol, degeneracy 100")
    print(f"  population A = {populations[0]:.3f}")
    print(f"  population B = {populations[1]:.3f}")
    print("  The higher-energy state can dominate when it has many more states.")
    print()
    print("Trajectory occupancy is a sampled indicator average, not a free energy.")
    print("A persistent contact can coexist with an unfavourable binding free energy.")

    assert ligand_a.persistent_contacts > ligand_b.persistent_contacts
    assert ligand_b.total < ligand_a.total
    assert populations[1] > populations[0]


if __name__ == "__main__":
    main()
