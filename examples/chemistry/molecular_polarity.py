# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Formal charge, partial charge, and dipole moment with RDKit.

Run with RDKit available:

    uv run --with rdkit python examples/chemistry/molecular_polarity.py
"""

from __future__ import annotations

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

from moleculab import dipole_moment


def formal_charge_report(name: str, smiles: str) -> None:
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"Invalid SMILES for {name}.")
    formal = [atom.GetFormalCharge() for atom in molecule.GetAtoms()]
    print(
        f"  {name:>12}: net={Chem.GetFormalCharge(molecule):+d} "
        f"atom formal charges={formal} "
        f"formula={rdMolDescriptors.CalcMolFormula(molecule)}"
    )


def main() -> None:
    print("Formal charge and protonation state")
    formal_charge_report("acetic acid", "CC(=O)O")
    formal_charge_report("acetate", "CC(=O)[O-]")
    formal_charge_report("methylamine", "CN")
    formal_charge_report("methylammonium", "C[NH3+]")
    print("  Formal charge is integer bookkeeping, not a partial-charge model.")
    print()

    print("Water: Gasteiger partial charges")
    water = Chem.AddHs(Chem.MolFromSmiles("O"))
    AllChem.EmbedMolecule(water, randomSeed=7)
    AllChem.MMFFOptimizeMolecule(water)
    AllChem.ComputeGasteigerCharges(water)
    charges = np.array(
        [float(atom.GetProp("_GasteigerCharge")) for atom in water.GetAtoms()]
    )
    conformer = water.GetConformer()
    positions_nm = np.array(
        [
            [
                conformer.GetAtomPosition(index).x * 0.1,
                conformer.GetAtomPosition(index).y * 0.1,
                conformer.GetAtomPosition(index).z * 0.1,
            ]
            for index in range(water.GetNumAtoms())
        ]
    )
    for atom, charge in zip(water.GetAtoms(), charges, strict=True):
        print(f"  {atom.GetSymbol():>2}: partial charge = {charge:+.3f} e")
    print(f"  net partial charge: {charges.sum():+.6f} e")
    dipole = dipole_moment(charges, positions_nm)
    print(
        "  dipole vector: "
        f"({dipole.vector_e_nm[0]:+.4f}, {dipole.vector_e_nm[1]:+.4f}, "
        f"{dipole.vector_e_nm[2]:+.4f}) e nm"
    )
    print(f"  dipole magnitude: {dipole.magnitude_e_nm:.4f} e nm")
    print(f"  dipole magnitude: {dipole.magnitude_debye:.3f} Debye")
    print("  Gasteiger charges are a model, not an experimental measurement.")
    print("  The dipole also depends on the conformer and charge model.")
    print()

    print("Force-field connection")
    print("  Coulomb energy uses the chosen q_i values:")
    print("    U_ij = k_e q_i q_j / r_ij")
    print("  Different charge models can produce different q_i and different")
    print("  electrostatics, even for the same molecular graph.")

    assert Chem.GetFormalCharge(Chem.MolFromSmiles("CC(=O)O")) == 0
    assert Chem.GetFormalCharge(Chem.MolFromSmiles("CC(=O)[O-]")) == -1
    assert Chem.GetFormalCharge(Chem.MolFromSmiles("C[NH3+]")) == 1
    assert abs(charges.sum()) < 1e-6
    assert dipole.magnitude_debye > 0.0


if __name__ == "__main__":
    main()
