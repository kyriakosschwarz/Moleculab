# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Analyse an external docking pose without running a docking search.

Run:

    uv run python examples/drug_discovery/docking_pose_analysis.py
"""

from __future__ import annotations

import numpy as np

from moleculab import (
    LigandPose,
    atom_id_mapping,
    binding_site_geometry,
    evaluate_pose_recovery,
    pose_rmsd,
    symmetry_aware_rmsd,
)


def main() -> None:
    reference = LigandPose(
        positions=np.array(
            [
                [0.00, 0.0, 0.0],
                [0.20, 0.0, 0.0],
                [-0.20, 0.0, 0.0],
            ]
        ),
        atom_ids=("C1", "O1", "O2"),
        source="reference pose",
    )
    pose = LigandPose(
        positions=np.array(
            [
                [0.02, 0.0, 0.0],
                [-0.20, 0.0, 0.0],
                [0.20, 0.0, 0.0],
            ]
        ),
        atom_ids=("C1", "O1", "O2"),
        score=-8.1,
        score_unit="arbitrary",
        score_source="example ranker",
        source="docked pose",
    )

    mapping = atom_id_mapping(pose, reference)
    identity_rmsd = pose_rmsd(pose.positions, reference.positions, mapping)
    swapped = [(0, 0), (1, 2), (2, 1)]
    symmetric = symmetry_aware_rmsd(
        pose.positions,
        reference.positions,
        [mapping, swapped],
    )

    print("Pose comparison")
    print(f"  atom mapping by ID: {mapping}")
    print(f"  RMSD with that mapping: {identity_rmsd:.4f} nm")
    print(
        "  symmetry-aware RMSD: "
        f"{symmetric.rmsd:.4f} nm "
        f"(mapping {symmetric.mapping_index + 1} of {symmetric.n_mappings})"
    )
    print()

    recovery = evaluate_pose_recovery(
        pose,
        reference,
        [mapping, swapped],
        threshold=0.05,
    )
    print("Pose recovery")
    print(f"  RMSD: {recovery.rmsd:.4f} nm")
    print(f"  threshold: {recovery.threshold:.3f} nm")
    print(f"  passed: {recovery.passed}")
    print()

    print("Score provenance")
    print(f"  score: {pose.score}")
    print(f"  unit: {pose.score_unit}")
    print(f"  source: {pose.score_source}")
    print("  This is a docking score, not a binding free energy.")
    print()

    receptor = np.array(
        [
            [0.30, 0.10, 0.00],
            [-0.30, 0.10, 0.00],
            [0.00, 0.30, 0.10],
            [0.00, -0.30, -0.10],
            [0.80, 0.00, 0.00],
        ]
    )
    site = binding_site_geometry(receptor, pose.positions, cutoff=0.5)
    print("Pocket around the docked pose")
    print(f"  selected receptor atoms: {site.atom_indices.tolist()}")
    print(
        "  centre: "
        f"({site.center[0]:.3f}, {site.center[1]:.3f}, "
        f"{site.center[2]:.3f}) nm"
    )
    print()
    print("Symmetry-aware RMSD uses an explicit mapping; no atom is guessed.")

    assert identity_rmsd > 0.3
    assert symmetric.rmsd < 0.02
    assert recovery.passed is True


if __name__ == "__main__":
    main()
