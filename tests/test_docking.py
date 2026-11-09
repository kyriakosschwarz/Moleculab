"""Tests for rigid-ligand docking with typed nonbonded scoring."""

import numpy as np
import pytest

from moleculab import (
    Box,
    DockingPose,
    DockingResult,
    LennardJonesParameters,
    NonbondedModel,
    Particle,
    ReactionFieldElectrostatics,
    RigidLigand,
    System,
    dock_rigid_ligand,
    score_pose,
)
from moleculab.units import COULOMB_CONSTANT


def lj_receptor(charge=0.0):
    model = NonbondedModel(
        parameters={"A": LennardJonesParameters(epsilon=1.0, sigma=0.3)},
        cutoff=1.0,
    )
    return System(
        particles=[
            Particle(1.0, [0.0, 0.0, 0.0], atom_type="A", charge=charge)
        ],
        nonbonded_model=model,
    )


def lj_ligand(charge=0.0):
    return RigidLigand(
        [Particle(1.0, [0.0, 0.0, 0.0], atom_type="A", charge=charge)]
    )


def test_rigid_ligand_center_and_local_coordinates():
    ligand = RigidLigand(
        [
            Particle(1.0, [0.0, 0.0, 0.0], atom_type="A"),
            Particle(1.0, [1.0, 0.0, 0.0], atom_type="B"),
        ]
    )

    np.testing.assert_allclose(ligand.center_of_mass, [0.5, 0.0, 0.0])
    np.testing.assert_allclose(
        ligand.local_positions,
        [[-0.5, 0.0, 0.0], [0.5, 0.0, 0.0]],
    )
    assert ligand.n_atoms == 2


def test_rigid_ligand_translation_and_rotation():
    ligand = RigidLigand(
        [
            Particle(1.0, [0.0, 0.0, 0.0], atom_type="A", charge=0.2),
            Particle(1.0, [1.0, 0.0, 0.0], atom_type="B", charge=-0.2),
        ]
    )
    identity = np.eye(3)
    rotation = np.array(
        [
            [0.0, -1.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
        ]
    )

    translated = ligand.transformed_positions([0.5, 0.0, 0.0], identity)
    rotated = ligand.transformed_positions([0.0, 0.0, 0.0], rotation)

    np.testing.assert_allclose(translated, [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    np.testing.assert_allclose(rotated, [[0.0, -0.5, 0.0], [0.0, 0.5, 0.0]])

    particles = ligand.transformed_particles([0.0, 0.0, 0.0], identity)
    assert particles[0].atom_type == "A"
    assert particles[1].charge == pytest.approx(-0.2)


def test_score_pose_matches_lennard_jones_minimum():
    receptor = lj_receptor()
    ligand = lj_ligand()
    minimum_distance = 2.0 ** (1.0 / 6.0) * 0.3

    pose = score_pose(
        receptor,
        ligand,
        [minimum_distance, 0.0, 0.0],
        np.eye(3),
    )

    assert pose.score == pytest.approx(-1.0)
    np.testing.assert_allclose(pose.translation, [minimum_distance, 0.0, 0.0])
    assert pose.positions.shape == (1, 3)


def test_score_pose_at_sigma_is_zero():
    pose = score_pose(lj_receptor(), lj_ligand(), [0.3, 0.0, 0.0], np.eye(3))

    assert pose.score == pytest.approx(0.0, abs=1e-12)


def test_score_pose_beyond_cutoff_is_zero():
    pose = score_pose(lj_receptor(), lj_ligand(), [2.0, 0.0, 0.0], np.eye(3))

    assert pose.score == 0.0


def test_score_pose_includes_reaction_field_coulomb():
    model = NonbondedModel(
        parameters={"A": LennardJonesParameters(epsilon=0.0, sigma=0.3)},
        cutoff=1.0,
        electrostatics=ReactionFieldElectrostatics(dielectric=1.0),
    )
    receptor = System(
        particles=[Particle(1.0, [0.0, 0.0, 0.0], atom_type="A", charge=1.0)],
        nonbonded_model=model,
    )
    ligand = RigidLigand(
        [Particle(1.0, [0.0, 0.0, 0.0], atom_type="A", charge=-1.0)]
    )

    pose = score_pose(receptor, ligand, [0.3, 0.0, 0.0], np.eye(3))

    expected = -COULOMB_CONSTANT * (1.0 / 0.3 - 1.0)
    assert pose.score == pytest.approx(expected)


def test_score_pose_uses_periodic_minimum_image():
    model = NonbondedModel(
        parameters={"A": LennardJonesParameters(epsilon=1.0, sigma=0.3)},
        cutoff=0.5,
    )
    receptor = System(
        particles=[Particle(1.0, [0.95, 0.0, 0.0], atom_type="A")],
        nonbonded_model=model,
        box=Box([1.0, 1.0, 1.0]),
    )
    ligand = lj_ligand()

    # Ligand at 0.05 and receptor at 0.95 are 0.1 nm apart through the boundary.
    pose = score_pose(receptor, ligand, [0.05, 0.0, 0.0], np.eye(3))
    direct = score_pose(lj_receptor(), ligand, [0.1, 0.0, 0.0], np.eye(3))

    assert pose.score == pytest.approx(direct.score)


def test_docking_finds_the_lennard_jones_minimum_distance():
    result = dock_rigid_ligand(
        lj_receptor(),
        lj_ligand(),
        center=[0.0, 0.0, 0.0],
        translation_radius=0.6,
        n_samples=2000,
        seed=42,
    )

    best = result.best
    best_distance = float(np.linalg.norm(best.translation))
    expected_distance = 2.0 ** (1.0 / 6.0) * 0.3

    assert isinstance(result, DockingResult)
    assert isinstance(best, DockingPose)
    assert best_distance == pytest.approx(expected_distance, abs=0.06)
    assert best.score < -0.9
    assert np.all(np.diff(result.scores) >= 0.0)


def test_docking_is_reproducible_with_a_seed():
    first = dock_rigid_ligand(
        lj_receptor(),
        lj_ligand(),
        center=[0.0, 0.0, 0.0],
        translation_radius=0.5,
        n_samples=200,
        seed=7,
    )
    second = dock_rigid_ligand(
        lj_receptor(),
        lj_ligand(),
        center=[0.0, 0.0, 0.0],
        translation_radius=0.5,
        n_samples=200,
        seed=7,
    )

    np.testing.assert_allclose(first.scores, second.scores)
    np.testing.assert_allclose(first.best.translation, second.best.translation)


def test_docking_result_returns_top_poses():
    result = dock_rigid_ligand(
        lj_receptor(),
        lj_ligand(),
        center=[0.0, 0.0, 0.0],
        translation_radius=0.5,
        n_samples=50,
        seed=3,
        max_results=5,
    )

    assert len(result.poses) == 5
    assert len(result.top(2)) == 2
    assert result.best is result.poses[0]


def test_docking_poses_are_orthonormal_rotations():
    result = dock_rigid_ligand(
        lj_receptor(),
        lj_ligand(),
        center=[0.0, 0.0, 0.0],
        translation_radius=0.4,
        n_samples=20,
        seed=11,
    )

    for pose in result.poses:
        np.testing.assert_allclose(pose.rotation.T @ pose.rotation, np.eye(3), atol=1e-10)
        assert np.linalg.det(pose.rotation) == pytest.approx(1.0)


def test_docking_validates_inputs():
    receptor = lj_receptor()
    ligand = lj_ligand()

    with pytest.raises(ValueError, match="nonbonded_model"):
        dock_rigid_ligand(
            System(particles=[Particle(1.0, [0.0, 0.0, 0.0], atom_type="A")]),
            ligand,
        )
    with pytest.raises(ValueError, match="n_samples"):
        dock_rigid_ligand(receptor, ligand, n_samples=0)
    with pytest.raises(ValueError, match="translation_radius"):
        dock_rigid_ligand(receptor, ligand, translation_radius=-0.1)
    with pytest.raises(ValueError, match="max_results"):
        dock_rigid_ligand(receptor, ligand, max_results=0)
    with pytest.raises(ValueError, match="center"):
        dock_rigid_ligand(receptor, ligand, center=[0.0, 0.0])
    with pytest.raises(ValueError, match="orthonormal"):
        score_pose(receptor, ligand, [0.3, 0.0, 0.0], np.ones((3, 3)))
    with pytest.raises(ValueError, match="at least one"):
        RigidLigand([])


def test_docking_validates_ligand_atom_types():
    model = NonbondedModel(
        parameters={"A": LennardJonesParameters(1.0, 0.3)},
        cutoff=1.0,
    )
    receptor = System(
        particles=[Particle(1.0, [0.0, 0.0, 0.0], atom_type="A")],
        nonbonded_model=model,
    )
    ligand = RigidLigand(
        [Particle(1.0, [0.0, 0.0, 0.0], atom_type="Z")]
    )

    with pytest.raises(ValueError, match="Missing"):
        score_pose(receptor, ligand, [0.3, 0.0, 0.0], np.eye(3))
