"""Stereochemistry, E/Z geometry, and conformational versus configurational isomerism.

Run with RDKit available:

    uv run --with rdkit python examples/chemistry/stereochemistry.py
"""

from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    FindMolChiralCenters,
    rdMolDescriptors,
    rdMolTransforms,
)


def describe(name: str, smiles: str) -> None:
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"Invalid SMILES for {name}: {smiles}")
    Chem.AssignStereochemistry(molecule, cleanIt=True, force=True)
    centers = FindMolChiralCenters(molecule, includeUnassigned=True)
    double_bond_stereo = [
        (info.centeredOn, info.descriptor.name)
        for info in Chem.FindPotentialStereo(molecule)
        if info.type.name == "Bond_Double"
    ]
    print(f"  {name:>8}: {Chem.MolToSmiles(molecule)}")
    print(
        f"            formula={rdMolDescriptors.CalcMolFormula(molecule)} "
        f"centers={centers}"
    )
    if double_bond_stereo:
        print(f"            double-bond stereo={double_bond_stereo}")
    print(f"            InChIKey={Chem.MolToInchiKey(molecule)}")


def main() -> None:
    print("Lactic acid enantiomers")
    describe("R", "C[C@@H](O)C(=O)O")
    describe("S", "C[C@H](O)C(=O)O")
    print("  Same formula and connectivity; different spatial configuration.")
    print()

    print("Tartaric acid stereoisomers")
    describe("R,R", "O[C@@H](C(=O)O)[C@@H](O)C(=O)O")
    describe("S,S", "O[C@H](C(=O)O)[C@H](O)C(=O)O")
    describe("meso", "O[C@@H](C(=O)O)[C@H](O)C(=O)O")
    print("  R,R and S,S are enantiomers; meso R,S is a diastereomer.")
    print()

    print("2-Butene E/Z double-bond stereochemistry")
    describe("E", "C/C=C/C")
    describe("Z", "C/C=C\\C")
    print("  E and Z are configurational isomers, not conformers.")
    print()

    print("Butane conformers: same graph, different geometry")
    butane_keys = []
    for seed in [7, 21]:
        molecule = Chem.AddHs(Chem.MolFromSmiles("CCCC"))
        AllChem.EmbedMolecule(molecule, randomSeed=seed)
        conformer = molecule.GetConformer()
        dihedral = rdMolTransforms.GetDihedralDeg(conformer, 0, 1, 2, 3)
        butane_keys.append(Chem.MolToInchiKey(Chem.RemoveHs(molecule)))
        print(f"  seed {seed}: dihedral C-C-C-C = {dihedral:7.2f} degrees")
    print("  These are conformers, not configurational stereoisomers.")
    print()

    print("Validation checks")
    r_lactic = Chem.MolFromSmiles("C[C@@H](O)C(=O)O")
    s_lactic = Chem.MolFromSmiles("C[C@H](O)C(=O)O")
    for molecule in (r_lactic, s_lactic):
        Chem.AssignStereochemistry(molecule, cleanIt=True, force=True)
    r_centers = FindMolChiralCenters(r_lactic, includeUnassigned=True)
    s_centers = FindMolChiralCenters(s_lactic, includeUnassigned=True)
    print(f"  CIP labels after assignment: R={r_centers}, S={s_centers}")
    print(
        "  same achiral descriptors: "
        f"MW={Descriptors.MolWt(r_lactic):.4f}, "
        f"TPSA={Descriptors.TPSA(r_lactic):.2f}, "
        f"HBD={Descriptors.NumHDonors(r_lactic)}, "
        f"HBA={Descriptors.NumHAcceptors(r_lactic)}"
    )
    print(
        "  enantiomer InChIKeys differ: "
        f"{Chem.MolToInchiKey(r_lactic) != Chem.MolToInchiKey(s_lactic)}"
    )
    for name, molecule in [("R", r_lactic), ("S", s_lactic)]:
        canonical = Chem.MolToSmiles(molecule)
        roundtrip = Chem.MolFromSmiles(canonical)
        preserved = Chem.MolToInchiKey(molecule) == Chem.MolToInchiKey(roundtrip)
        print(f"  {name} SMILES round trip preserves configuration: {preserved}")
    unspecified = Chem.MolFromSmiles("CC(O)C(=O)O")
    print(
        "  unspecified stereocentre: "
        f"{FindMolChiralCenters(unspecified, includeUnassigned=True)}"
    )
    print("  unspecified is not the same as racemic.")
    print(f"  butane configurational identity invariant to seed: {len(set(butane_keys)) == 1}")

    assert r_centers == [(1, "R")]
    assert s_centers == [(1, "S")]
    assert Chem.MolToInchiKey(r_lactic) != Chem.MolToInchiKey(s_lactic)
    assert len(set(butane_keys)) == 1


if __name__ == "__main__":
    main()
