# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Reproducible imatinib chemical-record and provenance anchor.

This script calculates only RDKit descriptors and chemical identity. It does
not download a structure, run docking, or estimate a binding free energy.

Run with RDKit available:

    uv run --with rdkit python examples/drug_discovery/capstone_imatinib.py
"""

from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import Descriptors, FindMolChiralCenters, rdMolDescriptors

IMATINIB_SMILES = (
    "CC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCN(CC3)C)"
    "NC4=NC=CC(=N4)C5=CN=CC=C5"
)


def main() -> None:
    molecule = Chem.MolFromSmiles(IMATINIB_SMILES)
    if molecule is None:
        raise ValueError("Imatinib SMILES failed to parse.")
    Chem.AssignStereochemistry(molecule, cleanIt=True, force=True)

    print("Imatinib chemical record")
    print(f"  input SMILES: {IMATINIB_SMILES}")
    print(f"  canonical SMILES: {Chem.MolToSmiles(molecule)}")
    print(f"  formula: {rdMolDescriptors.CalcMolFormula(molecule)}")
    print(f"  average mass: {Descriptors.MolWt(molecule):.4f} g/mol")
    print(f"  monoisotopic mass: {Descriptors.ExactMolWt(molecule):.6f} g/mol")
    print(f"  cLogP: {Descriptors.MolLogP(molecule):.3f}")
    print(f"  TPSA: {Descriptors.TPSA(molecule):.2f} Angstrom^2")
    print(f"  H-bond donors: {Descriptors.NumHDonors(molecule)}")
    print(f"  H-bond acceptors: {Descriptors.NumHAcceptors(molecule)}")
    print(f"  rotatable bonds: {Descriptors.NumRotatableBonds(molecule)}")
    print(f"  formal charge: {Chem.GetFormalCharge(molecule):+d}")
    print(
        "  stereocentres (including unassigned): "
        f"{FindMolChiralCenters(molecule, includeUnassigned=True)}"
    )
    print("  interpretation: no assigned or detected stereocentres in this representation;")
    print("  this is not a claim that stereochemistry is irrelevant to binding.")
    print()

    print("Provenance ledger")
    print("  experimental: none loaded by this script")
    print("  calculated: formula, mass, cLogP, TPSA, HBD/HBA, rotatable bonds")
    print("  predicted: none")
    print("  hypothesised: none")
    print("  structure source: not bundled; use a documented PDB entry")
    print("  docking/free energy: not computed by this script")
    print()

    print("Next steps in the capstone")
    print("  1. download and inspect PDB 2HYY or 1IEP with provenance recorded")
    print("  2. identify the ATP pocket and ligand copy")
    print("  3. choose protonation and tautomer states explicitly")
    print("  4. analyse contacts, H-bonds, and poses with the Track C tools")
    print("  5. keep measured, calculated, predicted, and hypothesised claims separate")

    assert rdMolDescriptors.CalcMolFormula(molecule) == "C29H31N7O"
    assert 493.0 < Descriptors.MolWt(molecule) < 494.0
    assert Chem.GetFormalCharge(molecule) == 0


if __name__ == "__main__":
    main()
