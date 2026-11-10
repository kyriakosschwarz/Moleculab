"""Tests for external docking pose records, RMSD, and pose recovery."""

import numpy as np
import pytest

from moleculab import (
    LigandPose,
    PoseRecovery,
    atom_id_mapping,
    evaluate_pose_recovery,
    pose_rmsd,
    symmetry_aware_rmsd,
    validate_atom_mapping,
)


def make_pose(
    positions,
    atom_ids=("C1", "C2", "C3"),
    *,
    score=None,
    score_unit=None,
    score_source="",
) -> LigandPose:
    return LigandPose(
        positions=np.asarray(positions, dtype=float),
        atom_ids=tuple(atom_ids),
        score=score,
        score_unit=score_unit,
        score_source=score_source,
    )


def test_rmsd_is_zero_for_identical_mapped_coordinates():
    positions = np.array([[0.0, 0.0, 0.0], [0.1, 0.0, 0.0]])

    assert pose_rmsd(positions, positions, [(0, 0), (1, 1)]) == 0.0


def test_rmsd_matches_a_hand_calculation():
    pose = np.array([[0.0, 0.0, 0.0], [0.1, 0.0, 0.0]])
    reference = np.array([[0.0, 0.0, 0.0], [0.1, 0.1, 0.0]])

    result = pose_rmsd(pose, reference, [(0, 0), (1, 1)])

    assert result == pytest.approx(np.sqrt(0.005))


def test_rmsd_is_invariant_to_a_consistent_rigid_body_transform():
    pose = np.array([[0.0, 0.0, 0.0], [0.1, 0.0, 0.0], [0.1, 0.2, 0.0]])
    reference = np.array([[0.02, 0.0, 0.0], [0.12, 0.0, 0.0], [0.12, 0.2, 0.0]])
    mapping = [(0, 0), (1, 1), (2, 2)]
    rotation = np.array(
        [
            [0.0, -1.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    shift = np.array([3.0, -2.0, 1.0])

    original = pose_rmsd(pose, reference, mapping)
    moved = pose_rmsd(pose @ rotation.T + shift, reference @ rotation.T + shift, mapping)

    assert moved == pytest.approx(original)


def test_rmsd_is_independent_of_array_order_with_an_explicit_mapping():
    pose = np.array([[0.0, 0.0, 0.0], [0.1, 0.0, 0.0]])
    reference = np.array([[0.1, 0.0, 0.0], [0.2, 0.0, 0.0]])

    direct = pose_rmsd(pose, reference, [(0, 0), (1, 1)])
    reordered = pose_rmsd(
        pose[::-1],
        reference,
        [(0, 1), (1, 0)],
    )

    assert reordered == pytest.approx(direct)


def test_atom_id_mapping_uses_stable_identifiers():
    reference = make_pose(
        [[0.0, 0.0, 0.0], [0.1, 0.0, 0.0]],
        atom_ids=("C1", "O1"),
    )
    pose = make_pose(
        [[0.1, 0.0, 0.0], [0.0, 0.0, 0.0]],
        atom_ids=("O1", "C1"),
    )

    mapping = atom_id_mapping(pose, reference)

    assert mapping == ((0, 1), (1, 0))
    assert pose_rmsd(pose.positions, reference.positions, mapping) == 0.0


def test_atom_id_mapping_rejects_missing_atoms():
    reference = make_pose([[0.0, 0.0, 0.0]], atom_ids=("C1",))
    pose = make_pose([[0.0, 0.0, 0.0]], atom_ids=("C2",))

    with pytest.raises(ValueError, match="atom IDs do not match"):
        atom_id_mapping(pose, reference)


def test_atom_mapping_validation_rejects_invalid_pairs():
    with pytest.raises(ValueError, match="unique"):
        validate_atom_mapping([(0, 0), (0, 1)], 2, 2, require_complete=False)
    with pytest.raises(ValueError, match="unique"):
        validate_atom_mapping([(0, 0), (1, 0)], 2, 2, require_complete=False)
    with pytest.raises(ValueError, match="out-of-range pose"):
        validate_atom_mapping([(2, 0)], 2, 1, require_complete=False)
    with pytest.raises(ValueError, match="out-of-range reference"):
        validate_atom_mapping([(0, 1)], 1, 1, require_complete=False)
    with pytest.raises(ValueError, match="at least one"):
        validate_atom_mapping([], 1, 1, require_complete=False)


def test_partial_mapping_requires_an_explicit_opt_in():
    pose = np.zeros((2, 3))
    reference = np.zeros((3, 3))

    with pytest.raises(ValueError, match="cover every"):
        pose_rmsd(pose, reference, [(0, 0), (1, 1)])
    assert pose_rmsd(
        pose,
        reference,
        [(0, 0), (1, 1)],
        require_complete=False,
    ) == 0.0


def test_symmetry_aware_rmsd_selects_the_best_explicit_mapping():
    reference = make_pose(
        [[0.0, 0.0, 0.0], [0.2, 0.0, 0.0], [0.4, 0.0, 0.0]],
        atom_ids=("C1", "O1", "O2"),
    )
    pose = make_pose(
        [[0.0, 0.0, 0.0], [0.4, 0.0, 0.0], [0.2, 0.0, 0.0]],
        atom_ids=("C1", "O1", "O2"),
    )
    identity = [(0, 0), (1, 1), (2, 2)]
    oxygen_swap = [(0, 0), (1, 2), (2, 1)]

    result = symmetry_aware_rmsd(
        pose.positions,
        reference.positions,
        [identity, oxygen_swap],
    )

    assert result.rmsd == 0.0
    assert result.mapping_index == 1
    assert result.best_mapping == tuple(oxygen_swap)
    assert result.n_mappings == 2


def test_symmetry_aware_rmsd_requires_at_least_one_mapping():
    with pytest.raises(ValueError, match="at least one mapping"):
        symmetry_aware_rmsd(np.zeros((1, 3)), np.zeros((1, 3)), [])


def test_pose_recovery_retains_score_provenance():
    reference = make_pose([[0.0, 0.0, 0.0], [0.1, 0.0, 0.0]], atom_ids=("C1", "O1"))
    pose = make_pose(
        [[0.0, 0.0, 0.0], [0.12, 0.0, 0.0]],
        atom_ids=("C1", "O1"),
        score=-7.2,
        score_unit="arbitrary",
        score_source="example ranker",
    )

    result = evaluate_pose_recovery(
        pose,
        reference,
        [[(0, 0), (1, 1)]],
        threshold=0.05,
    )

    assert isinstance(result, PoseRecovery)
    assert result.rmsd == pytest.approx(np.sqrt(0.0002))
    assert result.passed is True
    assert result.pose_score == -7.2
    assert result.score_unit == "arbitrary"
    assert result.score_source == "example ranker"


def test_pose_recovery_can_fail_the_threshold():
    reference = make_pose([[0.0, 0.0, 0.0]], atom_ids=("C1",))
    pose = make_pose([[0.2, 0.0, 0.0]], atom_ids=("C1",))

    result = evaluate_pose_recovery(
        pose,
        reference,
        [[(0, 0)]],
        threshold=0.1,
    )

    assert result.rmsd == pytest.approx(0.2)
    assert result.passed is False


def test_pose_recovery_validates_threshold():
    pose = make_pose([[0.0, 0.0, 0.0]], atom_ids=("C1",))

    with pytest.raises(ValueError, match="threshold"):
        evaluate_pose_recovery(pose, pose, [[(0, 0)]], threshold=0.0)


def test_ligand_pose_validates_score_provenance():
    with pytest.raises(ValueError, match="score_unit"):
        LigandPose(
            positions=[[0.0, 0.0, 0.0]],
            atom_ids=("C1",),
            score=-7.0,
            score_source="ranker",
        )
    with pytest.raises(ValueError, match="score_source"):
        LigandPose(
            positions=[[0.0, 0.0, 0.0]],
            atom_ids=("C1",),
            score=-7.0,
            score_unit="arbitrary",
        )


def test_ligand_pose_validates_atom_ids_and_shape():
    with pytest.raises(ValueError, match="shape"):
        LigandPose(positions=[[0.0, 0.0]], atom_ids=("C1",))
    with pytest.raises(ValueError, match="one entry"):
        LigandPose(positions=[[0.0, 0.0, 0.0]], atom_ids=("C1", "C2"))
    with pytest.raises(ValueError, match="unique"):
        LigandPose(positions=[[0.0, 0.0, 0.0], [0.1, 0.0, 0.0]], atom_ids=("C1", "C1"))
    with pytest.raises(ValueError, match="finite"):
        LigandPose(positions=[[0.0, 0.0, np.nan]], atom_ids=("C1",))
