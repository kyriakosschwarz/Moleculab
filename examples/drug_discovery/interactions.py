# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Hydrogen-bond geometry with prepared chemical features.

Run:

    uv run python examples/drug_discovery/interactions.py
"""

from __future__ import annotations

import numpy as np

from moleculab import (
    HydrogenBondCriteria,
    InteractionAtom,
    detect_hydrogen_bonds,
)


def main() -> None:
    criteria = HydrogenBondCriteria(max_distance=0.35, min_angle_degrees=120.0)

    protein_donor = InteractionAtom(
        position=[0.0, 0.0, 0.0],
        element="N",
        name="NZ",
        residue_id="LYS A 1",
        is_donor=True,
        formal_charge=1,
        hydrogen_positions=([0.1, 0.0, 0.0],),
    )
    ligand_acceptor = InteractionAtom(
        position=[0.30, 0.0, 0.0],
        element="O",
        name="O1",
        residue_id="LIG B 1",
        is_acceptor=True,
        formal_charge=-1,
    )

    bonds = detect_hydrogen_bonds(
        [protein_donor],
        [ligand_acceptor],
        criteria,
    )
    print("Hydrogen-bond detection")
    for bond in bonds:
        print(f"  type: {bond.interaction_type}")
        print(f"  protein: {bond.protein_residue} {bond.protein_label}")
        print(f"  ligand:  {bond.ligand_residue} {bond.ligand_label}")
        print(f"  donor is protein: {bond.donor_is_protein}")
        print(f"  donor-acceptor distance: {bond.distance:.3f} nm")
        print(f"                        = {bond.distance * 10.0:.3f} Angstrom")
        print(f"  D-H...A angle: {bond.angle_degrees:.1f} degrees")
    print()

    bad_angle = InteractionAtom(
        position=[0.05, 0.20, 0.0],
        element="O",
        name="O2",
        residue_id="LIG B 1",
        is_acceptor=True,
    )
    print("Short distance, unsuitable angle")
    print(f"  interactions: {len(detect_hydrogen_bonds([protein_donor], [bad_angle], criteria))}")
    print()

    no_hydrogen = InteractionAtom(
        position=[0.0, 0.0, 0.0],
        element="N",
        name="NZ",
        residue_id="LYS A 1",
        is_donor=True,
    )
    print("Donor without hydrogen positions")
    print(
        "  interactions: "
        f"{len(detect_hydrogen_bonds([no_hydrogen], [ligand_acceptor], criteria))}"
    )
    print()

    protein_acceptor = InteractionAtom(
        position=[0.30, 0.0, 0.0],
        element="O",
        name="OD1",
        residue_id="ASP A 2",
        is_acceptor=True,
        formal_charge=-1,
    )
    ligand_donor = InteractionAtom(
        position=[0.0, 0.0, 0.0],
        element="N",
        name="N1",
        residue_id="LIG B 1",
        is_donor=True,
        formal_charge=1,
        hydrogen_positions=([0.1, 0.0, 0.0],),
    )
    reverse = detect_hydrogen_bonds([protein_acceptor], [ligand_donor], criteria)
    print("Ligand donor to protein acceptor")
    if reverse:
        print(f"  donor is protein: {reverse[0].donor_is_protein}")
        print(f"  distance: {reverse[0].distance:.3f} nm")
        print(f"  angle: {reverse[0].angle_degrees:.1f} degrees")

    np.testing.assert_allclose(bonds[0].distance, 0.30)
    assert bonds[0].angle_degrees == 180.0


if __name__ == "__main__":
    main()
