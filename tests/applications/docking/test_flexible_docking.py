# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Tests for flexible-receptor docking with side-chain strain scoring."""

import numpy as np
import pytest

from moleculab import (
    LennardJonesParameters,
    NonbondedModel,
    Particle,
    RigidLigand,
    SideChain,
    System,
    dock_flexible_receptor,
    dock_rigid_ligand,
    make_receptor_conformer,
    receptor_strain_score,
    score_flexible_pose,
    score_pose,
)


def docking_setup():
    model = NonbondedModel(
        parameters={
            "R": LennardJonesParameters(epsilon=1.0, sigma=0.3),
            "L": LennardJonesParameters(epsilon=1.0, sigma=0.3),
        },
        cutoff=1.0,
    )
    receptor = System(
        particles=[
            Particle(1.0, [0.0, 0.0, 0.0], atom_type="R"),
            Particle(1.0, [0.3, 0.0, 0.1], atom_type="R"),
        ],
        nonbonded_model=model,
    )
    ligand = RigidLigand(
        [Particle(1.0, [0.0, 0.0, 0.0], atom_type="L")]
    )
    return receptor, ligand


def test_side_chain_validation():
    with pytest.raises(ValueError, match="at least one"):
        SideChain(())
    with pytest.raises(ValueError, match="unique"):
        SideChain((0, 0))
    with pytest.raises(ValueError, match="non-negative"):
        SideChain((-1,))
    with pytest.raises(ValueError, match="max_displacement"):
        SideChain((0,), max_displacement=-0.1)
    with pytest.raises(ValueError, match="stiffness"):
        SideChain((0,), stiffness=0.0)


def test_receptor_conformer_translates_only_side_chain():
    receptor, _ = docking_setup()
    chain = SideChain((1,), max_displacement=0.3, stiffness=10.0)
    displacement = np.array([0.1, -0.05, 0.02])

    conformer = make_receptor_conformer(receptor, [chain], [displacement])

    np.testing.assert_allclose(conformer.particles[0].position, [0.0, 0.0, 0.0])
    np.testing.assert_allclose(
        conformer.particles[1].position,
        receptor.particles[1].position + displacement,
    )
    assert conformer.particles[1].atom_type == "R"


def test_strain_score_is_harmonic():
    chain = SideChain((0,), stiffness=10.0)
    displacement = np.array([0.1, 0.2, 0.0])

    expected = 0.5 * 10.0 * (0.1**2 + 0.2**2)
    assert receptor_strain_score([chain], [displacement]) == pytest.approx(expected)


def test_flexible_pose_score_is_interaction_plus_strain():
    receptor, ligand = docking_setup()
    chain = SideChain((1,), max_displacement=0.3, stiffness=2.0)
    displacement = np.array([0.05, 0.0, 0.02])

    pose = score_flexible_pose(
        receptor,
        ligand,
        [chain],
        [displacement],
        translation=[0.3, 0.0, 0.0],
        rotation=np.eye(3),
    )
    conformer = make_receptor_conformer(receptor, [chain], [displacement])
    interaction = score_pose(
        conformer,
        ligand,
        [0.3, 0.0, 0.0],
        np.eye(3),
    ).score
    strain = receptor_strain_score([chain], [displacement])

    assert pose.interaction_score == pytest.approx(interaction)
    assert pose.strain_score == pytest.approx(strain)
    assert pose.score == pytest.approx(interaction + strain)


def test_flexible_docking_relieves_a_side_chain_clash():
    receptor, ligand = docking_setup()
    chain = SideChain((1,), max_displacement=0.3, stiffness=1.0)

    rigid = dock_rigid_ligand(
        receptor,
        ligand,
        center=[0.3, 0.0, 0.0],
        translation_radius=0.05,
        n_samples=100,
        seed=1,
    )
    flexible = dock_flexible_receptor(
        receptor,
        ligand,
        [chain],
        center=[0.3, 0.0, 0.0],
        translation_radius=0.05,
        n_ligand_samples=30,
        n_receptor_samples=80,
        seed=1,
    )

    assert rigid.best.score > 0.0
    assert flexible.best.score < 0.0
    assert flexible.best.score < rigid.best.score
    assert flexible.best.score == pytest.approx(
        flexible.best.interaction_score + flexible.best.strain_score
    )
    assert flexible.best.receptor_displacements[0].shape == (3,)


def test_flexible_docking_is_reproducible_and_sorted():
    receptor, ligand = docking_setup()
    chain = SideChain((1,), max_displacement=0.2, stiffness=2.0)
    kwargs = {
        "center": [0.3, 0.0, 0.0],
        "translation_radius": 0.04,
        "n_ligand_samples": 10,
        "n_receptor_samples": 10,
        "seed": 7,
    }

    first = dock_flexible_receptor(receptor, ligand, [chain], **kwargs)
    second = dock_flexible_receptor(receptor, ligand, [chain], **kwargs)

    np.testing.assert_allclose(first.scores, second.scores)
    assert np.all(np.diff(first.scores) >= 0.0)
    assert first.best is first.poses[0]
    assert len(first.top(3)) == 3


def test_flexible_docking_validation():
    receptor, ligand = docking_setup()
    chain = SideChain((1,))

    with pytest.raises(ValueError, match="nonbonded_model"):
        dock_flexible_receptor(
            System(particles=[Particle(1.0, [0.0, 0.0, 0.0])]),
            ligand,
            [chain],
        )
    with pytest.raises(ValueError, match="outside"):
        dock_flexible_receptor(
            receptor,
            ligand,
            [SideChain((5,))],
        )
    with pytest.raises(ValueError, match="n_ligand_samples"):
        dock_flexible_receptor(
            receptor, ligand, [chain], n_ligand_samples=0
        )
    with pytest.raises(ValueError, match="n_receptor_samples"):
        dock_flexible_receptor(
            receptor, ligand, [chain], n_receptor_samples=0
        )
    with pytest.raises(ValueError, match="max_results"):
        dock_flexible_receptor(
            receptor, ligand, [chain], max_results=0
        )
    with pytest.raises(ValueError, match="translation_radius"):
        dock_flexible_receptor(
            receptor, ligand, [chain], translation_radius=-0.1
        )
    with pytest.raises(ValueError, match="center"):
        dock_flexible_receptor(
            receptor, ligand, [chain], center=[0.0, 0.0]
        )
    with pytest.raises(ValueError, match="length"):
        make_receptor_conformer(receptor, [chain], [])
