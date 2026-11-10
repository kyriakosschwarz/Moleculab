# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Interaction fingerprints and occupancy from synthetic frames.

The frames are constructed explicitly so every occupancy can be checked by
hand. They do not come from MD.

Run:

    uv run python examples/drug_discovery/interaction_fingerprints.py
"""

from __future__ import annotations

import numpy as np

from moleculab import (
    FingerprintSchema,
    HydrogenBond,
    HydrogenBondCriteria,
    InteractionFeature,
    build_fingerprint,
    calculate_occupancy,
)

CRITERIA = HydrogenBondCriteria(max_distance=0.35, min_angle_degrees=120.0)


def make_hbond(residue: str, ligand_atom: str) -> HydrogenBond:
    return HydrogenBond(
        protein_index=0,
        ligand_index=0,
        donor_index=0,
        acceptor_index=0,
        donor_is_protein=True,
        protein_residue=residue,
        ligand_residue="LIG B 1",
        protein_label="OD1",
        ligand_label=ligand_atom,
        donor_label="OD1",
        acceptor_label=ligand_atom,
        distance=0.29,
        angle_degrees=165.0,
        hydrogen_position=np.array([0.0, 0.0, 0.0]),
        hydrogen_index=0,
        criteria=CRITERIA,
    )


def main() -> None:
    first = make_hbond("ASP A 123", "N7")
    second = make_hbond("PHE A 87", "O2")
    third = make_hbond("LYS A 45", "N1")
    schema = FingerprintSchema.from_interactions([first, second, third])

    print("Feature schema")
    for index, feature in enumerate(schema.features):
        print(
            f"  {index}: {feature.interaction_type} | "
            f"{feature.protein_residue} | {feature.ligand_residue} | "
            f"{feature.ligand_atom_id}"
        )
    print()

    # Eight valid frames followed by two missing frames.
    frames = [
        build_fingerprint([first, second], schema),
        build_fingerprint([first], schema),
        build_fingerprint([first], schema),
        build_fingerprint([first], schema),
        build_fingerprint([first], schema),
        build_fingerprint([first], schema),
        build_fingerprint([], schema),
        build_fingerprint([second], schema),
        None,
        None,
    ]

    print("Frame fingerprints (1 = present)")
    for frame_number, frame in enumerate(frames, start=1):
        if frame is None:
            print(f"  frame {frame_number:2d}: missing")
            continue
        values = " ".join("1" if bit else "0" for bit in frame.bits)
        print(f"  frame {frame_number:2d}: {values}")
    print()

    result = calculate_occupancy(frames)
    print("Sampled-frame occupancy")
    for feature, occupancy, count in zip(
        result.features,
        result.occupancy,
        result.counts,
        strict=True,
    ):
        print(
            f"  {feature.protein_residue} / {feature.ligand_atom_id}: "
            f"{occupancy:.3f} ({count}/{result.n_valid_frames} valid frames)"
        )
    print(
        f"  total frames: {result.n_frames}; valid: {result.n_valid_frames}; "
        f"missing: {result.n_missing_frames}"
    )
    assert result.criteria == CRITERIA
    print(
        "  criteria: "
        f"distance <= {result.criteria.max_distance:.2f} nm, "
        f"angle >= {result.criteria.min_angle_degrees:.0f} degrees"
    )
    print()
    print(
        "Sampled-frame occupancy describes the analysed frames, "
        "not an equilibrium probability"
    )

    # The hand-calculated expectations: A in 6/8 valid frames, B in 2/8, C in 0.
    values = dict(zip(result.features, result.occupancy, strict=True))
    assert values[InteractionFeature.from_hydrogen_bond(first)] == 0.75
    assert values[InteractionFeature.from_hydrogen_bond(second)] == 0.25
    assert values[InteractionFeature.from_hydrogen_bond(third)] == 0.0


if __name__ == "__main__":
    main()
