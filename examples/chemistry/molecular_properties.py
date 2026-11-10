"""Molecular descriptors and chemical similarity with RDKit.

Run with RDKit available:

    uv run --with rdkit python examples/chemistry/molecular_properties.py
"""

from __future__ import annotations

try:
    from rdkit import Chem
    from rdkit.Chem import DataStructs, Descriptors, rdFingerprintGenerator
    from rdkit.Chem import rdMolDescriptors
except ImportError as exc:  # pragma: no cover - depends on the environment
    raise SystemExit(
        "RDKit is required for this example. "
        "Run it with: uv run --with rdkit python examples/chemistry/molecular_properties.py"
    ) from exc

import numpy as np

MOLECULES = {
    "aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "ibuprofen": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "naproxen": "COc1ccc2cc(ccc2c1)[C@@H](C)C(=O)O",
    "acetaminophen": "CC(=O)Nc1ccc(O)cc1",
    "caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
}


def main() -> None:
    molecules = {
        name: Chem.MolFromSmiles(smiles) for name, smiles in MOLECULES.items()
    }
    for name, molecule in molecules.items():
        if molecule is None:
            raise ValueError(f"Invalid SMILES for {name}.")

    print("Basic descriptors")
    print(
        f"{'molecule':>14} {'formula':>10} {'MW':>8} {'exact':>8} "
        f"{'HBD':>4} {'HBA':>4} {'TPSA':>7} {'rot':>4} {'charge':>6}"
    )
    for name, molecule in molecules.items():
        print(
            f"{name:>14} {rdMolDescriptors.CalcMolFormula(molecule):>10} "
            f"{Descriptors.MolWt(molecule):8.2f} "
            f"{Descriptors.ExactMolWt(molecule):8.4f} "
            f"{Descriptors.NumHDonors(molecule):4d} "
            f"{Descriptors.NumHAcceptors(molecule):4d} "
            f"{Descriptors.TPSA(molecule):7.1f} "
            f"{Descriptors.NumRotatableBonds(molecule):4d} "
            f"{Chem.GetFormalCharge(molecule):6d}"
        )
    print()

    # A hand-checkable Tanimoto example.
    vector_a = np.array([1, 0, 1, 1, 0, 0, 1], dtype=bool)
    vector_b = np.array([1, 1, 0, 1, 0, 1, 1], dtype=bool)
    shared = int(np.sum(vector_a & vector_b))
    union = int(np.sum(vector_a | vector_b))
    print("Manual Tanimoto example")
    print(f"  shared bits: {shared}")
    print(f"  union bits:  {union}")
    print(f"  Tanimoto:    {shared / union:.3f}")
    print()

    generator = rdFingerprintGenerator.GetMorganGenerator(
        radius=2,
        fpSize=2048,
    )
    fingerprints = {
        name: generator.GetFingerprint(molecule)
        for name, molecule in molecules.items()
    }

    names = list(molecules)
    print("Morgan fingerprint Tanimoto similarity (radius 2, 2048 bits)")
    print(f"{'':>14}" + "".join(f"{name:>14}" for name in names))
    for first in names:
        row = "".join(
            f"{DataStructs.TanimotoSimilarity(fingerprints[first], fingerprints[second]):14.3f}"
            for second in names
        )
        print(f"{first:>14}{row}")
    print()
    print("Structural similarity is not biological equivalence.")
    print("Descriptors are models and screening aids, not measured activities.")


if __name__ == "__main__":
    main()
