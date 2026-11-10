# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Binding-site geometry from Cartesian coordinates.

This example is parser-independent: the coordinates are already in nm and the
residue labels are supplied by the caller.

Run:

    uv run python examples/drug_discovery/binding_site.py
"""

from __future__ import annotations

import numpy as np

from moleculab import (
    binding_site_geometry,
    group_indices_by_residue,
    minimum_distances,
    select_pocket_atoms,
)


def main() -> None:
    ligand = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.2, 0.0, 0.0],
        ]
    )
    protein = np.array(
        [
            [0.10, 0.00, 0.00],
            [0.00, 0.20, 0.00],
            [0.00, 0.00, 0.35],
            [0.60, 0.00, 0.00],
            [1.20, 0.00, 0.00],
        ]
    )
    residue_labels = ["ALA A 1", "ALA A 1", "GLY A 2", "SER A 3", "VAL A 4"]

    distances = minimum_distances(protein, ligand)
    print("Atom minimum distances (nm)")
    for label, distance in zip(residue_labels, distances, strict=True):
        print(f"  {label:>8}: {distance:.3f}")

    cutoff = 0.5
    geometry = binding_site_geometry(protein, ligand, cutoff=cutoff)
    print()
    print(f"Pocket selection (strict, cutoff {cutoff:.3f} nm)")
    print(f"  selected atom indices: {geometry.atom_indices.tolist()}")
    print(f"  selected atoms: {geometry.n_atoms}")
    grouped = group_indices_by_residue(
        geometry.atom_indices,
        [residue_labels[index] for index in geometry.atom_indices],
    )
    print(f"  residues: {', '.join(grouped)}")
    print(
        "  centre: "
        f"({geometry.center[0]:.3f}, {geometry.center[1]:.3f}, "
        f"{geometry.center[2]:.3f}) nm"
    )
    print(
        "  box dimensions: "
        f"({geometry.dimensions[0]:.3f}, {geometry.dimensions[1]:.3f}, "
        f"{geometry.dimensions[2]:.3f}) nm"
    )
    print()

    edge_atom = np.array([[0.50, 0.0, 0.0]])
    edge_ligand = np.array([[0.0, 0.0, 0.0]])
    strict = select_pocket_atoms(edge_atom, edge_ligand, cutoff)
    inclusive = select_pocket_atoms(edge_atom, edge_ligand, cutoff, inclusive=True)
    print("Cutoff boundary")
    print(f"  strict <  {cutoff:.3f}: {strict.tolist()}")
    print(f"  inclusive <= {cutoff:.3f}: {inclusive.tolist()}")
    print()

    wrapped_protein = np.array([[0.95, 0.0, 0.0]])
    wrapped_ligand = np.array([[0.05, 0.0, 0.0]])
    print("Periodic-boundary note")
    print(
        "  raw distance across a 1 nm box edge: "
        f"{minimum_distances(wrapped_protein, wrapped_ligand)[0]:.3f} nm"
    )
    print("  no minimum image is applied by this geometry layer.")


if __name__ == "__main__":
    main()
