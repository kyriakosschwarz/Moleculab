# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Hydrogen-bond donors, acceptors, and geometric criteria.

Run with RDKit available:

    uv run --with rdkit python examples/chemistry/hydrogen_bonding.py
"""

from __future__ import annotations

import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors

from moleculab import (
    HydrogenBondCriteria,
    InteractionAtom,
    detect_hydrogen_bonds,
)


def main() -> None:
    print("Donor and acceptor counts (RDKit)")
    for name, smiles in [
        ("ethanol", "CCO"),
        ("acetic acid", "CC(=O)O"),
        ("acetate", "CC(=O)[O-]"),
        ("methylamine", "CN"),
        ("methylammonium", "C[NH3+]"),
        ("acetamide", "CC(=O)N"),
    ]:
        molecule = Chem.MolFromSmiles(smiles)
        print(
            f"  {name:>15}: HBD={Descriptors.NumHDonors(molecule)} "
            f"HBA={Descriptors.NumHAcceptors(molecule)}"
        )
    print("  Protonation changes donor/acceptor status and molecular charge.")
    print()

    criteria = HydrogenBondCriteria(max_distance=0.35, min_angle_degrees=120.0)
    donor = InteractionAtom(
        position=[0.0, 0.0, 0.0],
        element="N",
        name="N1",
        residue_id="DONOR",
        is_donor=True,
        hydrogen_positions=([0.1, 0.0, 0.0],),
    )
    linear_acceptor = InteractionAtom(
        position=[0.30, 0.0, 0.0],
        element="O",
        name="O1",
        residue_id="ACCEPTOR",
        is_acceptor=True,
    )
    bad_angle_acceptor = InteractionAtom(
        position=[0.05, 0.20, 0.0],
        element="O",
        name="O2",
        residue_id="ACCEPTOR",
        is_acceptor=True,
    )
    no_hydrogen_donor = InteractionAtom(
        position=[0.0, 0.0, 0.0],
        element="N",
        name="N2",
        residue_id="DONOR",
        is_donor=True,
    )

    linear = detect_hydrogen_bonds([donor], [linear_acceptor], criteria)
    bad = detect_hydrogen_bonds([donor], [bad_angle_acceptor], criteria)
    no_h = detect_hydrogen_bonds([no_hydrogen_donor], [linear_acceptor], criteria)

    print("Geometric criteria")
    print("  linear donor-acceptor contact:")
    for bond in linear:
        print(f"    distance={bond.distance:.3f} nm angle={bond.angle_degrees:.1f} degrees")
    print(f"  short contact with unsuitable angle: {len(bad)} H-bond(s)")
    print(f"  donor without hydrogen positions: {len(no_h)} H-bond(s)")
    print()
    print("A geometric H-bond is not automatically a favourable binding term.")
    print("Solvent competition, desolvation, and entropy also matter.")

    assert len(linear) == 1
    assert linear[0].angle_degrees == 180.0
    assert bad == ()
    assert no_h == ()


if __name__ == "__main__":
    main()
