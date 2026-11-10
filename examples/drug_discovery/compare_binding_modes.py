# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Compare binding modes by pose geometry and interaction patterns.

Run:

    uv run python examples/drug_discovery/compare_binding_modes.py
"""

from __future__ import annotations

import numpy as np

from moleculab import (
    FingerprintSchema,
    HydrogenBond,
    HydrogenBondCriteria,
    LigandPose,
    atom_id_mapping,
    build_fingerprint,
    cluster_poses,
    compare_binding_modes,
    compare_interaction_patterns,
    compare_pose_geometry,
)

CRITERIA = HydrogenBondCriteria()


def make_pose(positions, label):
    return LigandPose(
        positions=np.asarray(positions, dtype=float),
        atom_ids=("C1", "O1"),
        metadata={"label": label},
    )


def make_hbond(residue, ligand_atom):
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


def make_fingerprint(*bonds):
    schema = FingerprintSchema.from_interactions(bonds)
    return build_fingerprint(bonds, schema)


def main() -> None:
    pose_a = make_pose([[0.0, 0.0, 0.0], [0.1, 0.0, 0.0]], "A")
    pose_b_same = make_pose([[0.0, 0.0, 0.0], [0.1, 0.0, 0.0]], "B")
    pose_c_shifted = make_pose([[0.5, 0.0, 0.0], [0.6, 0.0, 0.0]], "C")

    fingerprint_a = make_fingerprint(make_hbond("ASP A 123", "N7"))
    fingerprint_b = make_fingerprint(make_hbond("PHE A 87", "O2"))
    fingerprint_c = make_fingerprint(make_hbond("ASP A 123", "N9"))

    same_pose = compare_binding_modes(
        pose_a,
        pose_b_same,
        [atom_id_mapping(pose_a, pose_b_same)],
        fingerprint_a,
        fingerprint_b,
    )
    print("Same geometry, different interactions")
    print(f"  RMSD: {same_pose.geometry.rmsd:.3f} nm")
    print(f"  exact feature Jaccard: {same_pose.interactions.jaccard:.3f}")
    print()

    shifted_pose = compare_binding_modes(
        pose_a,
        pose_c_shifted,
        [atom_id_mapping(pose_a, pose_c_shifted)],
        fingerprint_a,
        fingerprint_c,
        include_ligand_atom=False,
    )
    print("Different geometry, same residue-level pattern")
    print(f"  RMSD: {shifted_pose.geometry.rmsd:.3f} nm")
    print(f"  residue-pattern Jaccard: {shifted_pose.interactions.jaccard:.3f}")
    print()

    residue_comparison = compare_interaction_patterns(
        fingerprint_a,
        fingerprint_c,
        include_ligand_atom=False,
    )
    print("Residue-level shared features")
    print(f"  shared: {residue_comparison.shared}")
    print(f"  only in A: {residue_comparison.only_in_a}")
    print(f"  only in B: {residue_comparison.only_in_b}")
    print()

    clustering = cluster_poses(
        [pose_a, pose_b_same, pose_c_shifted],
        threshold=0.2,
    )
    print("Pose clustering (single linkage)")
    print(f"  threshold: {clustering.threshold:.2f} nm")
    for index, cluster in enumerate(clustering.clusters):
        labels = [
            [pose_a, pose_b_same, pose_c_shifted][member].metadata["label"]
            for member in cluster.member_indices
        ]
        print(f"  cluster {index}: {labels}")
    print()

    geometry = compare_pose_geometry(
        pose_a,
        pose_c_shifted,
        [atom_id_mapping(pose_a, pose_c_shifted)],
    )
    print("Pose geometry")
    print(f"  RMSD: {geometry.rmsd:.3f} nm")
    print(f"  centroid distance: {geometry.centroid_distance:.3f} nm")
    print()
    print("A shared pose and a shared interaction pattern are different questions.")

    assert same_pose.geometry.rmsd == 0.0
    assert same_pose.interactions.jaccard == 0.0
    assert shifted_pose.geometry.rmsd == 0.5
    assert shifted_pose.interactions.jaccard == 1.0
    assert clustering.n_clusters == 2


if __name__ == "__main__":
    main()
