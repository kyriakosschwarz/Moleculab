"""Inspect the hierarchy, uncertainty, and units of a protein structure.

The example uses a small synthetic PDB excerpt so the output is deterministic.
It shows what a structure file contains before any docking or energy
calculation: chains, residues, alternate locations, waters, hetero groups, and
the Angstrom-to-nanometre boundary.

Run with gemmi available:

    uv run --with gemmi python examples/drug_discovery/protein_structure.py
"""

from __future__ import annotations

import numpy as np

try:
    import gemmi
except ImportError as exc:  # pragma: no cover - depends on the environment
    raise SystemExit(
        "gemmi is required for this example. "
        "Run it with: uv run --with gemmi python examples/drug_discovery/protein_structure.py"
    ) from exc

# A teaching excerpt, not a complete experimental structure. Residue 2 is
# intentionally absent, residue 3 is incomplete, and ALA 1 has two CB
# alternate locations.
PDB_TEXT = """\
ATOM      1  N   ALA A   1      11.104   6.134  -6.504  1.00 20.00           N
ATOM      2  CA  ALA A   1      11.639   6.071  -5.147  1.00 20.00           C
ATOM      3  C   ALA A   1      13.149   6.233  -5.164  1.00 20.00           C
ATOM      4  O   ALA A   1      13.750   5.240  -4.737  1.00 20.00           O
ATOM      5  CB  ALA A   1      11.086   4.878  -4.360  0.60 30.00           C
ATOM      6  CB AALA A   1      10.900   4.900  -4.500  0.40 30.00           C
ATOM      7  N   GLY A   3      14.000   7.000  -5.500  1.00 25.00           N
HETATM    8  O   HOH A 101      10.000   5.000  -8.000  1.00 40.00           O
HETATM    9  C1  LIG B   1       2.000   3.000   4.000  1.00 15.00           C
END
"""


def main() -> None:
    structure = gemmi.read_pdb_string(PDB_TEXT)
    model = structure[0]

    chains = list(model)
    atoms = [atom for chain in model for residue in chain for atom in residue]
    polymer = [
        (chain.name, residue)
        for chain in model
        for residue in chain
        if residue.het_flag == "A"
    ]
    waters = [
        (chain.name, residue)
        for chain in model
        for residue in chain
        if residue.name == "HOH"
    ]
    other_hetero = [
        (chain.name, residue)
        for chain in model
        for residue in chain
        if residue.het_flag == "H" and residue.name != "HOH"
    ]

    print("Structure summary")
    print(f"  source format: {structure.input_format.name}")
    print(f"  models: {len(structure)}")
    print(f"  chains: {', '.join(chain.name for chain in chains)}")
    print(f"  atom records: {len(atoms)}")
    print(f"  polymer residues: {len(polymer)}")
    print(f"  waters: {len(waters)}")
    print(f"  other hetero residues: {len(other_hetero)}")
    print()

    print("Polymer residues")
    for chain_name, residue in polymer:
        count = len(residue)
        noun = "record" if count == 1 else "records"
        print(
            f"  {residue.name} {chain_name} {residue.seqid.num}: "
            f"{count} atom {noun}"
        )

    chain_a_numbers = [
        residue.seqid.num
        for chain_name, residue in polymer
        if chain_name == "A"
    ]
    for left, right in zip(chain_a_numbers, chain_a_numbers[1:]):
        if right != left + 1:
            print(f"  residue-number gap in chain A: {left} -> {right}")
    print()

    print("Alternate locations")
    for chain_name, residue in polymer:
        atoms_by_name: dict[str, list] = {}
        for atom in residue:
            atoms_by_name.setdefault(atom.name, []).append(atom)
        for atom_name, group in atoms_by_name.items():
            if len(group) > 1:
                for atom in group:
                    altloc = "" if atom.altloc == "\x00" else atom.altloc
                    print(
                        f"  {residue.name} {chain_name} {residue.seqid.num} "
                        f"{atom_name}: altloc {altloc!r} occupancy {atom.occ:.2f}"
                    )
    print()

    first = atoms[0]
    position_angstrom = np.array([first.pos.x, first.pos.y, first.pos.z])
    position_nm = position_angstrom * 0.1
    print("Coordinate boundary")
    print(
        f"  first atom {first.name}: "
        f"({position_angstrom[0]:.3f}, {position_angstrom[1]:.3f}, "
        f"{position_angstrom[2]:.3f}) Angstrom"
    )
    print(
        f"  converted: ({position_nm[0]:.4f}, {position_nm[1]:.4f}, "
        f"{position_nm[2]:.4f}) nm"
    )
    print()
    print("Provenance: experimental teaching excerpt (declared by the caller)")
    print("Missing atoms are left missing; the example does not rebuild them.")

    np.testing.assert_allclose(position_nm, [1.1104, 0.6134, -0.6504], atol=1e-12)
    assert len(atoms) == 9
    assert [chain.name for chain in chains] == ["A", "B"]


if __name__ == "__main__":
    main()
