"""Geometric contacts, occupancy, and lifetimes from synthetic frames.

Run:

    uv run python examples/trajectory_contacts.py
"""

from __future__ import annotations

import numpy as np

from moleculab import (
    aggregate_contacts_by_residue,
    calculate_contact_lifetimes,
    calculate_contact_matrix,
    calculate_contact_occupancy,
)


def main() -> None:
    protein = np.array(
        [
            [0.10, 0.00, 0.00],
            [0.00, 0.20, 0.00],
            [0.00, 0.00, 0.40],
            [1.00, 0.00, 0.00],
        ]
    )
    ligand = np.array(
        [
            [0.00, 0.00, 0.00],
            [0.00, 0.00, 0.50],
        ]
    )
    residue_labels = ["ASP A 123", "ASP A 123", "PHE A 87", "LYS A 45"]

    contacts = calculate_contact_matrix(protein, ligand, cutoff=0.5)
    residue_contacts = aggregate_contacts_by_residue(contacts, residue_labels)

    print("Atom contact matrix (rows = protein atoms, columns = ligand atoms)")
    print(contacts.astype(int))
    print()
    print("Residue contact matrix")
    for residue, row in zip(
        residue_contacts.residue_ids,
        residue_contacts.contacts,
        strict=True,
    ):
        print(f"  {residue}: {row.astype(int).tolist()}")
    print()

    continuous = np.array([True] * 6 + [False] * 6)
    intermittent = np.array([True, False] * 6, dtype=bool)
    frames = [
        np.array([[continuous[index], intermittent[index]]], dtype=bool)
        for index in range(len(continuous))
    ]
    occupancy = calculate_contact_occupancy(frames)
    print("Sampled-frame occupancy")
    print(f"  continuous contact: {occupancy.occupancy[0, 0]:.3f}")
    print(f"  intermittent contact: {occupancy.occupancy[0, 1]:.3f}")
    print(f"  valid frames: {occupancy.n_valid_frames}/{occupancy.n_frames}")
    print()

    spacing = 0.1
    continuous_lifetime = calculate_contact_lifetimes(
        continuous,
        frame_spacing=spacing,
    )
    intermittent_lifetime = calculate_contact_lifetimes(
        intermittent,
        frame_spacing=spacing,
    )
    print(f"Contact lifetimes (frame spacing = {spacing} ns)")
    print(
        "  continuous: runs "
        f"{continuous_lifetime.run_lengths_frames}, "
        f"mean {continuous_lifetime.mean_lifetime:.3f} ns"
    )
    print(
        "  intermittent: runs "
        f"{intermittent_lifetime.run_lengths_frames}, "
        f"mean {intermittent_lifetime.mean_lifetime:.3f} ns"
    )
    print()
    print("Same occupancy, very different temporal behaviour.")
    print("Contact persistence is not a binding free energy.")

    assert contacts.tolist() == [
        [True, False],
        [True, False],
        [True, True],
        [False, False],
    ]
    assert residue_contacts.residue_ids == ("ASP A 123", "PHE A 87", "LYS A 45")
    assert occupancy.occupancy.tolist() == [[0.5, 0.5]]
    assert np.isclose(continuous_lifetime.mean_lifetime, 0.6)
    assert np.isclose(intermittent_lifetime.mean_lifetime, 0.1)


if __name__ == "__main__":
    main()
