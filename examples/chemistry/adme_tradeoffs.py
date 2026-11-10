# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Physicochemical trade-offs and LogP/LogD for several drugs.

Run with RDKit available:

    uv run --with rdkit python examples/chemistry/adme_tradeoffs.py
"""

from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import Descriptors

from moleculab import logd_from_logp

MOLECULES = {
    "aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "ibuprofen": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "atenolol": "CC(C)NCC(O)COc1ccc(CC(N)=O)cc1",
    "caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "imatinib": "CC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCN(CC3)C)NC4=NC=CC(=N4)C5=CN=CC=C5",
}


def main() -> None:
    molecules = {
        name: Chem.MolFromSmiles(smiles) for name, smiles in MOLECULES.items()
    }
    for name, molecule in molecules.items():
        if molecule is None:
            raise ValueError(f"Invalid SMILES for {name}.")

    print("Calculated physicochemical descriptors")
    print(
        f"{'molecule':>10} {'MW':>8} {'cLogP':>7} {'TPSA':>7} "
        f"{'HBD':>4} {'HBA':>4} {'rot':>4} {'charge':>6}"
    )
    for name, molecule in molecules.items():
        print(
            f"{name:>10} {Descriptors.MolWt(molecule):8.2f} "
            f"{Descriptors.MolLogP(molecule):7.2f} "
            f"{Descriptors.TPSA(molecule):7.1f} "
            f"{Descriptors.NumHDonors(molecule):4d} "
            f"{Descriptors.NumHAcceptors(molecule):4d} "
            f"{Descriptors.NumRotatableBonds(molecule):4d} "
            f"{Chem.GetFormalCharge(molecule):6d}"
        )
    print()

    ibuprofen = molecules["ibuprofen"]
    atenolol = molecules["atenolol"]
    ibu_logp = Descriptors.MolLogP(ibuprofen)
    aten_logp = Descriptors.MolLogP(atenolol)
    ibu_logd = logd_from_logp(ibu_logp, 4.9, 7.4, ionisation="acid")
    aten_logd = logd_from_logp(aten_logp, 9.6, 7.4, ionisation="base")

    print("LogP versus LogD at pH 7.4")
    print(
        f"  ibuprofen (acid, pKa 4.9): "
        f"LogP={ibu_logp:.3f}, LogD={ibu_logd:.3f}"
    )
    print(
        f"  atenolol (base, pKa 9.6): "
        f"LogP={aten_logp:.3f}, LogD={aten_logd:.3f}"
    )
    print("  LogP describes the neutral form; LogD includes ionisation at pH.")
    print("  The pKa values are inputs, not predictions.")
    print()

    print("Interpretation cautions")
    print("  Ibuprofen is more lipophilic and has lower TPSA than atenolol,")
    print("  but its acid is ionised at pH 7.4, so LogD is lower than LogP.")
    print("  Atenolol is polar and its base is protonated at pH 7.4;")
    print("  its LogD is negative, which is consistent with low passive permeability.")
    print("  Imatinib has high MW and TPSA but is still an oral drug;")
    print("  descriptor thresholds are not proof of bioavailability.")
    print()
    print("Measured ADME outcomes, calculated descriptors, and predictive models")
    print("are different categories and should not be collapsed into one score.")

    assert ibu_logd < ibu_logp
    assert aten_logd < aten_logp
    assert Descriptors.MolWt(molecules["imatinib"]) > 400.0


if __name__ == "__main__":
    main()
