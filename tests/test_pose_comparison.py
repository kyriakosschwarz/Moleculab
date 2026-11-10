# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Tests for pose clustering and binding-mode comparison."""

import numpy as np
import pytest

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
    pose_distance,
)

CRITERIA = HydrogenBondCriteria()


def make_pose(positions, atom_ids=None, label="pose"):
    coordinates = np.asarray(positions, dtype=float)
    if atom_ids is None:
        atom_ids = tuple(f"C{index + 1}" for index in range(coordinates.shape[0]))
    return LigandPose(
        positions=coordinates,
        atom_ids=tuple(atom_ids),
        metadata={"label": label},
    )


def make_hbond(
    protein_residue: str,
    ligand_atom: str,
    ligand_residue: str = "LIG B 1",
) -> HydrogenBond:
    return HydrogenBond(
        protein_index=0,
        ligand_index=0,
        donor_index=0,
        acceptor_index=0,
        donor_is_protein=True,
        protein_residue=protein_residue,
        ligand_residue=ligand_residue,
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


def test_pose_geometry_comparison_reports_rmsd_and_centroid_shift():
    pose_a = make_pose([[0.0, 0.0, 0.0], [0.1, 0.0, 0.0]], atom_ids=("C1", "O1"))
    pose_b = make_pose([[0.1, 0.0, 0.0], [0.2, 0.0, 0.0]], atom_ids=("C1", "O1"))
    mapping = atom_id_mapping(pose_a, pose_b)

    comparison = compare_pose_geometry(pose_a, pose_b, [mapping])

    assert comparison.rmsd == pytest.approx(0.1)
    assert comparison.centroid_distance == pytest.approx(0.1)
    np.testing.assert_allclose(comparison.centroid_displacement, [-0.1, 0.0, 0.0])


def test_geometry_comparison_uses_symmetry_aware_mapping():
    reference = make_pose(
        [[0.0, 0.0, 0.0], [0.2, 0.0, 0.0], [-0.2, 0.0, 0.0]],
        atom_ids=("C1", "O1", "O2"),
    )
    swapped = make_pose(
        [[0.0, 0.0, 0.0], [-0.2, 0.0, 0.0], [0.2, 0.0, 0.0]],
        atom_ids=("C1", "O1", "O2"),
    )
    identity = atom_id_mapping(swapped, reference)
    oxygen_swap = [(0, 0), (1, 2), (2, 1)]

    comparison = compare_pose_geometry(
        swapped,
        reference,
        [identity, oxygen_swap],
    )

    assert comparison.rmsd == pytest.approx(0.0)
    assert comparison.mapping_index == 1


def test_exact_feature_comparison_reports_shared_and_unique_features():
    first = make_fingerprint(
        make_hbond("ASP A 123", "N7"),
        make_hbond("PHE A 87", "O2"),
    )
    second = make_fingerprint(
        make_hbond("ASP A 123", "N7"),
        make_hbond("LYS A 45", "N1"),
    )

    comparison = compare_interaction_patterns(first, second)

    assert comparison.level == "feature"
    assert comparison.jaccard == pytest.approx(1.0 / 3.0)
    assert len(comparison.shared) == 1
    assert len(comparison.only_in_a) == 1
    assert len(comparison.only_in_b) == 1


def test_residue_patterns_ignore_ligand_atom_identity():
    first = make_fingerprint(make_hbond("ASP A 123", "N7"))
    second = make_fingerprint(make_hbond("ASP A 123", "N9"))

    exact = compare_interaction_patterns(first, second)
    residue = compare_interaction_patterns(first, second, include_ligand_atom=False)

    assert exact.jaccard == 0.0
    assert residue.jaccard == 1.0
    assert residue.level == "residue_pattern"


def test_low_rmsd_poses_can_have_different_interaction_patterns():
    pose_a = make_pose([[0.0, 0.0, 0.0], [0.1, 0.0, 0.0]], atom_ids=("C1", "O1"))
    pose_b = make_pose([[0.0, 0.0, 0.0], [0.1, 0.0, 0.0]], atom_ids=("C1", "O1"))
    first = make_fingerprint(make_hbond("ASP A 123", "N7"))
    second = make_fingerprint(make_hbond("PHE A 87", "O2"))

    result = compare_binding_modes(
        pose_a,
        pose_b,
        [atom_id_mapping(pose_a, pose_b)],
        first,
        second,
    )

    assert result.geometry.rmsd == 0.0
    assert result.interactions.jaccard == 0.0


def test_similar_interaction_patterns_can_have_different_geometry():
    pose_a = make_pose([[0.0, 0.0, 0.0], [0.1, 0.0, 0.0]], atom_ids=("C1", "O1"))
    pose_b = make_pose([[0.5, 0.0, 0.0], [0.6, 0.0, 0.0]], atom_ids=("C1", "O1"))
    first = make_fingerprint(make_hbond("ASP A 123", "N7"))
    second = make_fingerprint(make_hbond("ASP A 123", "N9"))

    result = compare_binding_modes(
        pose_a,
        pose_b,
        [atom_id_mapping(pose_a, pose_b)],
        first,
        second,
        include_ligand_atom=False,
    )

    assert result.geometry.rmsd == pytest.approx(0.5)
    assert result.interactions.jaccard == 1.0


def test_cluster_poses_groups_by_distance():
    poses = [
        make_pose([[0.0, 0.0, 0.0], [0.1, 0.0, 0.0]], label="A"),
        make_pose([[0.1, 0.0, 0.0], [0.2, 0.0, 0.0]], label="B"),
        make_pose([[1.0, 0.0, 0.0], [1.1, 0.0, 0.0]], label="C"),
    ]

    result = cluster_poses(poses, threshold=0.2)

    assert result.n_clusters == 2
    assert result.clusters[0].member_indices == (0, 1)
    assert result.clusters[1].member_indices == (2,)
    assert result.labels.tolist() == [0, 0, 1]


def test_clustering_partition_is_independent_of_input_order():
    poses = [
        make_pose([[0.0, 0.0, 0.0], [0.1, 0.0, 0.0]], label="A"),
        make_pose([[0.1, 0.0, 0.0], [0.2, 0.0, 0.0]], label="B"),
        make_pose([[1.0, 0.0, 0.0], [1.1, 0.0, 0.0]], label="C"),
    ]

    def partition(result, pose_list):
        return {
            frozenset(
                pose_list[index].metadata["label"]
                for index in cluster.member_indices
            )
            for cluster in result.clusters
        }

    original = partition(cluster_poses(poses, threshold=0.2), poses)
    reordered_poses = list(reversed(poses))
    reordered = partition(
        cluster_poses(reordered_poses, threshold=0.2),
        reordered_poses,
    )

    assert original == reordered


def test_symmetry_requires_an_explicit_distance_function():
    reference = make_pose(
        [[0.0, 0.0, 0.0], [0.2, 0.0, 0.0], [-0.2, 0.0, 0.0]],
        atom_ids=("C1", "O1", "O2"),
        label="reference",
    )
    swapped = make_pose(
        [[0.0, 0.0, 0.0], [-0.2, 0.0, 0.0], [0.2, 0.0, 0.0]],
        atom_ids=("C1", "O1", "O2"),
        label="swapped",
    )

    default = cluster_poses([reference, swapped], threshold=0.05)
    assert default.n_clusters == 2

    def symmetric_distance(first, second):
        return pose_distance(
            first,
            second,
            atom_mappings=[
                atom_id_mapping(first, second),
                [(0, 0), (1, 2), (2, 1)],
            ],
        )

    explicit = cluster_poses(
        [reference, swapped],
        threshold=0.05,
        distance_function=symmetric_distance,
    )
    assert explicit.n_clusters == 1


def test_default_pose_distance_requires_matching_atom_ids():
    first = make_pose([[0.0, 0.0, 0.0]], atom_ids=("C1",))
    second = make_pose([[0.0, 0.0, 0.0]], atom_ids=("C2",))

    with pytest.raises(ValueError, match="atom IDs do not match"):
        pose_distance(first, second)


def test_single_linkage_can_chain_through_intermediate_poses():
    poses = [
        make_pose([[0.0, 0.0, 0.0], [0.1, 0.0, 0.0]], label="A"),
        make_pose([[0.6, 0.0, 0.0], [0.7, 0.0, 0.0]], label="B"),
        make_pose([[1.2, 0.0, 0.0], [1.3, 0.0, 0.0]], label="C"),
    ]

    result = cluster_poses(poses, threshold=0.65)

    assert result.n_clusters == 1
    assert result.clusters[0].member_indices == (0, 1, 2)


def test_cluster_threshold_is_inclusive():
    poses = [
        make_pose([[0.0, 0.0, 0.0], [0.1, 0.0, 0.0]], label="A"),
        make_pose([[0.2, 0.0, 0.0], [0.3, 0.0, 0.0]], label="B"),
    ]

    result = cluster_poses(poses, threshold=0.2)

    assert result.n_clusters == 1
    assert result.clusters[0].member_indices == (0, 1)


def test_duplicate_poses_have_stable_cluster_assignments():
    first = make_pose([[0.0, 0.0, 0.0], [0.1, 0.0, 0.0]], label="A")
    second = make_pose([[0.0, 0.0, 0.0], [0.1, 0.0, 0.0]], label="B")

    result = cluster_poses([first, second, first], threshold=0.01)

    assert result.n_clusters == 1
    assert result.clusters[0].member_indices == (0, 1, 2)
    assert result.labels.tolist() == [0, 0, 0]


def test_missing_fingerprints_are_not_treated_as_empty_fingerprints():
    fingerprint = make_fingerprint(make_hbond("ASP A 123", "N7"))

    with pytest.raises(TypeError, match="InteractionFingerprint"):
        compare_interaction_patterns(None, fingerprint)
    with pytest.raises(TypeError, match="InteractionFingerprint"):
        compare_interaction_patterns(fingerprint, None)


def test_cluster_validation():
    pose = make_pose([[0.0, 0.0, 0.0]], atom_ids=("C1",))

    with pytest.raises(ValueError, match="at least one"):
        cluster_poses([], threshold=0.1)
    with pytest.raises(ValueError, match="threshold"):
        cluster_poses([pose], threshold=0.0)
    with pytest.raises(TypeError, match="distance_function"):
        cluster_poses([pose], threshold=0.1, distance_function="not callable")


def test_single_pose_forms_one_cluster():
    pose = make_pose([[0.0, 0.0, 0.0]], atom_ids=("C1",))

    result = cluster_poses([pose], threshold=0.1)

    assert result.n_clusters == 1
    assert result.labels.tolist() == [0]
    assert result.distance_matrix.shape == (1, 1)
