"""Tests for hydrogen-bond geometry and the interaction feature record."""

import numpy as np
import pytest

from moleculab import (
    HydrogenBondCriteria,
    InteractionAtom,
    detect_hydrogen_bonds,
)


def linear_donor(hydrogen_x: float = 0.1) -> InteractionAtom:
    return InteractionAtom(
        position=[0.0, 0.0, 0.0],
        element="N",
        name="N1",
        residue_id="LYS A 1",
        is_donor=True,
        hydrogen_positions=([hydrogen_x, 0.0, 0.0],),
    )


def test_linear_hydrogen_bond_is_detected():
    protein = [linear_donor()]
    ligand = [
        InteractionAtom(
            position=[0.30, 0.0, 0.0],
            element="O",
            name="O1",
            residue_id="LIG B 1",
            is_acceptor=True,
            formal_charge=-1,
        )
    ]

    hydrogen_bonds = detect_hydrogen_bonds(protein, ligand)

    assert len(hydrogen_bonds) == 1
    bond = hydrogen_bonds[0]
    assert bond.interaction_type == "hydrogen_bond"
    assert bond.donor_is_protein is True
    assert bond.protein_residue == "LYS A 1"
    assert bond.ligand_residue == "LIG B 1"
    assert bond.distance == pytest.approx(0.30)
    assert bond.angle_degrees == pytest.approx(180.0)
    np.testing.assert_allclose(bond.hydrogen_position, [0.1, 0.0, 0.0])


def test_ligand_donor_to_protein_acceptor_is_detected():
    protein = [
        InteractionAtom(
            position=[0.30, 0.0, 0.0],
            element="O",
            name="O1",
            residue_id="ASP A 2",
            is_acceptor=True,
        )
    ]
    ligand = [linear_donor()]

    hydrogen_bonds = detect_hydrogen_bonds(protein, ligand)

    assert len(hydrogen_bonds) == 1
    assert hydrogen_bonds[0].donor_is_protein is False
    assert hydrogen_bonds[0].protein_residue == "ASP A 2"
    assert hydrogen_bonds[0].ligand_residue == "LYS A 1"


def test_invalid_donor_acceptor_combinations_are_not_detected():
    protein = [linear_donor()]
    ligand_donor = [
        InteractionAtom(
            position=[0.30, 0.0, 0.0],
            element="N",
            is_donor=True,
            hydrogen_positions=([0.4, 0.0, 0.0],),
        )
    ]
    ligand_acceptor = [
        InteractionAtom(position=[0.30, 0.0, 0.0], element="O", is_acceptor=True)
    ]
    protein_acceptor = [
        InteractionAtom(position=[0.0, 0.0, 0.0], element="O", is_acceptor=True)
    ]

    assert detect_hydrogen_bonds(protein, ligand_donor) == ()
    assert detect_hydrogen_bonds(protein_acceptor, ligand_acceptor) == ()


def test_short_distance_with_bad_angle_is_rejected():
    protein = [linear_donor()]
    ligand = [
        InteractionAtom(
            position=[0.05, 0.20, 0.0],
            element="O",
            is_acceptor=True,
        )
    ]

    hydrogen_bonds = detect_hydrogen_bonds(protein, ligand)

    assert hydrogen_bonds == ()


def test_distance_boundary_is_inclusive():
    criteria = HydrogenBondCriteria(max_distance=0.35, min_angle_degrees=120.0)
    ligand_at_boundary = [
        InteractionAtom(position=[0.35, 0.0, 0.0], element="O", is_acceptor=True)
    ]
    ligand_beyond = [
        InteractionAtom(position=[0.3501, 0.0, 0.0], element="O", is_acceptor=True)
    ]

    assert len(detect_hydrogen_bonds([linear_donor()], ligand_at_boundary, criteria)) == 1
    assert detect_hydrogen_bonds([linear_donor()], ligand_beyond, criteria) == ()


def test_angle_boundary_is_inclusive():
    hydrogen = np.array([0.1, 0.0, 0.0])
    acceptor = hydrogen + 0.2 * np.array([0.5, np.sqrt(3.0) / 2.0, 0.0])
    ligand = [
        InteractionAtom(position=acceptor, element="O", is_acceptor=True)
    ]

    at_boundary = HydrogenBondCriteria(max_distance=0.35, min_angle_degrees=120.0)
    above_boundary = HydrogenBondCriteria(
        max_distance=0.35,
        min_angle_degrees=120.1,
    )

    assert len(detect_hydrogen_bonds([linear_donor()], ligand, at_boundary)) == 1
    assert detect_hydrogen_bonds([linear_donor()], ligand, above_boundary) == ()


def test_most_linear_hydrogen_is_selected():
    protein = [
        InteractionAtom(
            position=[0.0, 0.0, 0.0],
            element="N",
            is_donor=True,
            hydrogen_positions=(
                [0.1, 0.0, 0.0],
                [0.0, 0.1, 0.0],
            ),
        )
    ]
    ligand = [
        InteractionAtom(position=[0.3, 0.0, 0.0], element="O", is_acceptor=True)
    ]

    hydrogen_bonds = detect_hydrogen_bonds(protein, ligand)

    assert len(hydrogen_bonds) == 1
    assert hydrogen_bonds[0].angle_degrees == pytest.approx(180.0)
    assert hydrogen_bonds[0].hydrogen_index == 0
    np.testing.assert_allclose(hydrogen_bonds[0].hydrogen_position, [0.1, 0.0, 0.0])


def test_donor_without_hydrogen_positions_is_skipped():
    protein = [
        InteractionAtom(
            position=[0.0, 0.0, 0.0],
            element="N",
            is_donor=True,
        )
    ]
    ligand = [
        InteractionAtom(position=[0.3, 0.0, 0.0], element="O", is_acceptor=True)
    ]

    assert detect_hydrogen_bonds(protein, ligand) == ()


def test_protonation_state_changes_hydrogen_bond_classification():
    acceptor = [
        InteractionAtom(position=[0.3, 0.0, 0.0], element="O", is_acceptor=True)
    ]
    protonated = [
        InteractionAtom(
            position=[0.0, 0.0, 0.0],
            element="N",
            is_donor=True,
            hydrogen_positions=([0.1, 0.0, 0.0],),
        )
    ]
    deprotonated = [
        InteractionAtom(
            position=[0.0, 0.0, 0.0],
            element="N",
            is_donor=False,
        )
    ]

    assert len(detect_hydrogen_bonds(protonated, acceptor)) == 1
    assert detect_hydrogen_bonds(deprotonated, acceptor) == ()


def test_hydrogen_bond_geometry_is_translation_and_rotation_invariant():
    protein = [linear_donor()]
    ligand = [
        InteractionAtom(position=[0.30, 0.0, 0.0], element="O", is_acceptor=True)
    ]
    rotation = np.array(
        [
            [0.0, -1.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    shift = np.array([1.0, -2.0, 3.0])

    original = detect_hydrogen_bonds(protein, ligand)[0]
    moved_protein = [
        InteractionAtom(
            position=protein[0].position @ rotation.T + shift,
            element="N",
            is_donor=True,
            hydrogen_positions=tuple(
                hydrogen @ rotation.T + shift
                for hydrogen in protein[0].hydrogen_positions
            ),
        )
    ]
    moved_ligand = [
        InteractionAtom(
            position=ligand[0].position @ rotation.T + shift,
            element="O",
            is_acceptor=True,
        )
    ]
    moved = detect_hydrogen_bonds(moved_protein, moved_ligand)[0]

    assert moved.distance == pytest.approx(original.distance)
    assert moved.angle_degrees == pytest.approx(original.angle_degrees)


def test_detection_is_reproducible_under_atom_reordering():
    protein = [
        InteractionAtom(
            position=[0.0, 0.0, 0.0],
            element="N",
            name="N1",
            residue_id="LYS A 1",
            is_donor=True,
            hydrogen_positions=([0.1, 0.0, 0.0],),
        ),
        InteractionAtom(
            position=[0.5, 0.0, 0.0],
            element="N",
            name="N2",
            residue_id="LYS A 2",
            is_donor=True,
            hydrogen_positions=([0.6, 0.0, 0.0],),
        ),
    ]
    ligand = [
        InteractionAtom(
            position=[0.3, 0.0, 0.0],
            element="O",
            name="O1",
            residue_id="LIG B 1",
            is_acceptor=True,
        ),
        InteractionAtom(
            position=[0.8, 0.0, 0.0],
            element="O",
            name="O2",
            residue_id="LIG B 2",
            is_acceptor=True,
        ),
    ]

    def signature(bonds):
        return {
            (
                bond.protein_residue,
                bond.ligand_residue,
                round(bond.distance, 12),
                round(bond.angle_degrees, 12),
            )
            for bond in bonds
        }

    assert signature(detect_hydrogen_bonds(protein, ligand)) == signature(
        detect_hydrogen_bonds(list(reversed(protein)), list(reversed(ligand)))
    )


def test_interaction_feature_validation():
    with pytest.raises(ValueError, match="position"):
        InteractionAtom(position=[0.0, 0.0])
    with pytest.raises(ValueError, match="hydrogen_positions"):
        InteractionAtom(position=[0.0, 0.0, 0.0], hydrogen_positions=([0.0, 0.0],))
    with pytest.raises(ValueError, match="max_distance"):
        HydrogenBondCriteria(max_distance=0.0)
    with pytest.raises(ValueError, match="min_angle_degrees"):
        HydrogenBondCriteria(min_angle_degrees=200.0)
    with pytest.raises(TypeError, match="InteractionAtom"):
        detect_hydrogen_bonds(["not an atom"], [])
    with pytest.raises(TypeError, match="criteria"):
        detect_hydrogen_bonds([], [], criteria="not criteria")


def test_empty_atom_lists_return_no_interactions():
    assert detect_hydrogen_bonds([], []) == ()
