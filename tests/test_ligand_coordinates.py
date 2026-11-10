# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Tests for ligand coordinate reconciliation."""

from __future__ import annotations

import math
from dataclasses import replace
from pathlib import Path

import pytest

from moleculab import (
    LIGAND_COORDINATE_SCHEMA_VERSION,
    CoordinateAtom,
    DepositedBond,
    GeometryCriteria,
    LigandSelection,
    load_ligand_coordinates,
    reconcile_ligand_coordinates,
    validate_ligand_coordinates,
    write_ligand_coordinates,
)

CHECKSUM = "c" * 64


def selection() -> LigandSelection:
    return LigandSelection(
        component_id="LIG",
        label_asym_id="A",
        auth_asym_id="A",
        instance="A",
        auth_seq_id="1",
        available_instances=("A",),
    )


def atom(
    name: str,
    element: str,
    position: tuple[float, float, float],
    reference_index: int,
    *,
    occupancy: float | None = 1.0,
    b_factor: float | None = 10.0,
    alternate_location: str | None = None,
) -> CoordinateAtom:
    return CoordinateAtom(
        atom_name=name,
        element=element,
        reference_index=reference_index,
        position_angstrom=position,
        occupancy=occupancy,
        b_iso_or_equiv=b_factor,
        alternate_location=alternate_location,
        model_number=1,
    )


def record(
    *,
    atoms: tuple[CoordinateAtom, ...] | None = None,
    bonds: tuple[DepositedBond, ...] | None = None,
    expected: tuple[str, ...] = ("A1", "A2", "A3"),
    component: dict[str, str] | None = None,
    criteria: GeometryCriteria | None = None,
):
    if atoms is None:
        atoms = (
            atom("A1", "C", (0.0, 0.0, 0.0), 0),
            atom("A2", "C", (1.4, 0.0, 0.0), 1),
            atom("A3", "O", (2.2, 1.2, 0.0), 2),
        )
    if bonds is None:
        bonds = (
            DepositedBond("A1", "A2", "SINGLE"),
            DepositedBond("A2", "A3", "SINGLE"),
        )
    if component is None:
        component = {"A1": "C", "A2": "C", "A3": "O"}
    return reconcile_ligand_coordinates(
        source_manifest_path="manifest.json",
        source_manifest_sha256=CHECKSUM,
        identity_record_path="identity.json",
        identity_record_sha256=CHECKSUM,
        accession="TEST",
        selection=selection(),
        atoms=atoms,
        bonds=bonds,
        expected_atom_names=expected,
        component_atom_elements=component,
        reference_stereocenters=(),
        deposited_stereochemistry=(),
        criteria=criteria,
    )


def status(result, check_id: str) -> str:
    return next(
        check.status for check in result.validation.checks if check.check_id == check_id
    )


def test_complete_coordinate_record_is_valid():
    result = record()

    assert validate_ligand_coordinates(result) == []
    assert result.schema_version == LIGAND_COORDINATE_SCHEMA_VERSION
    assert len(result.atoms) == 3
    assert len(result.bonds) == 2
    assert len(result.angles) == 1
    assert status(result, "coordinates.mapped_atoms_complete") == "pass"
    assert status(result, "geometry.bond_lengths_in_range") == "pass"
    assert status(result, "geometry.angles_in_range") == "pass"
    assert status(result, "geometry.close_contacts") == "pass"
    assert result.missing_heavy_atoms == ()
    assert result.discrepancies == ()


def test_coordinate_record_json_round_trip(tmp_path):
    result = record()
    path = write_ligand_coordinates(tmp_path / "coordinates.json", result)
    loaded = load_ligand_coordinates(path)

    assert loaded == result
    assert loaded.to_dict()["criteria"]["close_contact_min_angstrom"] == 1.5
    assert "least-squares plane" in loaded.to_dict()["conventions"]["plane_fit_method"]
    assert "1-2 and 1-3" in loaded.to_dict()["conventions"]["close_contact_definition"]


def test_missing_coordinate_record_is_reported():
    result = record(expected=("A1", "A2", "A3", "A4"))

    assert status(result, "coordinates.mapped_atoms_complete") == "fail"
    discrepancy = next(
        item for item in result.discrepancies if item.category == "missing_coordinates"
    )
    assert discrepancy.status == "open"
    assert "A4" in discrepancy.detail
    assert validate_ligand_coordinates(result) == []


def test_missing_heavy_atom_is_reported_without_inventing_coordinates():
    component = {"A1": "C", "A2": "C", "A3": "O", "A4": "N"}
    result = record(component=component)

    assert result.missing_heavy_atoms == ("A4",)
    assert status(result, "coordinates.missing_heavy_atoms") == "fail"
    assert any(
        item.category == "missing_heavy_atoms" and item.status == "open"
        for item in result.discrepancies
    )


def test_missing_occupancy_and_b_factor_are_reported():
    atoms = (
        atom("A1", "C", (0.0, 0.0, 0.0), 0, occupancy=None, b_factor=None),
        atom("A2", "C", (1.4, 0.0, 0.0), 1),
        atom("A3", "O", (2.2, 1.2, 0.0), 2),
    )
    result = record(atoms=atoms)

    assert status(result, "coordinates.occupancies_recorded") == "fail"
    assert status(result, "coordinates.b_factors_recorded") == "fail"
    assert any(item.category == "occupancy" for item in result.discrepancies)
    assert any(item.category == "b_factor" for item in result.discrepancies)


def test_alternate_locations_are_preserved_and_geometry_is_withheld():
    atoms = (
        atom("A1", "C", (0.0, 0.0, 0.0), 0, alternate_location="A"),
        atom("A1", "C", (0.1, 0.0, 0.0), 0, alternate_location="B"),
        atom("A2", "C", (1.4, 0.0, 0.0), 1),
        atom("A3", "O", (2.2, 1.2, 0.0), 2),
    )
    result = record(atoms=atoms)

    assert result.alternate_location_atoms == ("A1",)
    assert status(result, "coordinates.alternate_locations_documented") == "pass"
    assert status(result, "geometry.bond_lengths_in_range") == "not_run"
    assert result.bonds == ()
    assert any(
        item.category == "alternate_locations" and item.status == "open"
        for item in result.discrepancies
    )
    assert validate_ligand_coordinates(result) == []


def test_bond_length_out_of_range_is_reported():
    atoms = (
        atom("A1", "C", (0.0, 0.0, 0.0), 0),
        atom("A2", "C", (3.0, 0.0, 0.0), 1),
        atom("A3", "O", (3.8, 1.2, 0.0), 2),
    )
    result = record(atoms=atoms)

    assert status(result, "geometry.bond_lengths_in_range") == "fail"
    assert any(item.category == "bond_length" for item in result.discrepancies)


def test_angle_out_of_range_is_reported():
    atoms = (
        atom("A1", "C", (1.0, 0.0, 0.0), 0),
        atom("A2", "C", (0.0, 0.0, 0.0), 1),
        atom("A3", "O", (0.9848, 0.1736, 0.0), 2),
    )
    result = record(atoms=atoms)

    assert status(result, "geometry.angles_in_range") == "fail"
    assert any(item.category == "bond_angle" for item in result.discrepancies)
    assert result.angles[0].angle_degrees == pytest.approx(10.0, abs=0.1)


def test_ring_planarity_is_checked_from_aromatic_bonds():
    names = ("R1", "R2", "R3", "R4", "R5", "R6")
    angles = [2.0 * math.pi * index / 6 for index in range(6)]
    positions = {
        name: (1.4 * math.cos(angle), 1.4 * math.sin(angle), 0.0)
        for name, angle in zip(names, angles)
    }
    atoms = tuple(
        atom(name, "C", positions[name], index)
        for index, name in enumerate(names)
    )
    bonds = tuple(
        DepositedBond(names[index], names[(index + 1) % 6], "AROMATIC", aromatic=True)
        for index in range(6)
    )
    result = record(
        atoms=atoms,
        bonds=bonds,
        expected=names,
        component={name: "C" for name in names},
    )

    assert len(result.rings) == 1
    assert status(result, "geometry.ring_planarity") == "pass"
    assert result.rings[0].rmsd_angstrom < 0.01

    distorted = tuple(
        replace(atom_value, position_angstrom=(*atom_value.position_angstrom[:2], 1.5))
        if atom_value.atom_name == "R1"
        else atom_value
        for atom_value in atoms
    )
    bad = record(
        atoms=distorted,
        bonds=bonds,
        expected=names,
        component={name: "C" for name in names},
    )

    assert status(bad, "geometry.ring_planarity") == "fail"
    assert any(item.category == "ring_planarity" for item in bad.discrepancies)


def test_nonbonded_close_contact_is_reported():
    atoms = (
        atom("A1", "C", (0.0, 0.0, 0.0), 0),
        atom("A2", "C", (1.5, 0.0, 0.0), 1),
        atom("A3", "C", (0.0, 0.0, 1.0), 2),
        atom("A4", "C", (1.5, 0.0, 1.0), 3),
    )
    bonds = (
        DepositedBond("A1", "A2", "SINGLE"),
        DepositedBond("A3", "A4", "SINGLE"),
    )
    result = record(
        atoms=atoms,
        bonds=bonds,
        expected=("A1", "A2", "A3", "A4"),
        component={"A1": "C", "A2": "C", "A3": "C", "A4": "C"},
    )

    assert status(result, "geometry.close_contacts") == "fail"
    assert any(item.category == "close_contact" for item in result.discrepancies)


def test_reference_stereocentres_are_not_inferred_from_coordinates():
    result = reconcile_ligand_coordinates(
        source_manifest_path="manifest.json",
        source_manifest_sha256=CHECKSUM,
        identity_record_path="identity.json",
        identity_record_sha256=CHECKSUM,
        accession="TEST",
        selection=selection(),
        atoms=(
            atom("A1", "C", (0.0, 0.0, 0.0), 0),
            atom("A2", "C", (1.4, 0.0, 0.0), 1),
            atom("A3", "O", (2.2, 1.2, 0.0), 2),
        ),
        bonds=(
            DepositedBond("A1", "A2", "SINGLE"),
            DepositedBond("A2", "A3", "SINGLE"),
        ),
        expected_atom_names=("A1", "A2", "A3"),
        component_atom_elements={"A1": "C", "A2": "C", "A3": "O"},
        reference_stereocenters=("1:R",),
        deposited_stereochemistry=(),
    )

    assert status(result, "coordinates.stereochemistry") == "not_run"
    assert any(
        item.category == "stereochemistry" and item.status == "open"
        for item in result.discrepancies
    )


def test_validation_rejects_bad_coordinate_record():
    result = record()
    bad_atoms = (replace(result.atoms[0], occupancy=1.5),)
    bad = replace(
        result,
        atoms=bad_atoms,
        criteria=GeometryCriteria(bond_length_min_angstrom=2.0, bond_length_max_angstrom=1.0),
    )

    errors = validate_ligand_coordinates(bad)

    assert any("occupancy must be between" in error for error in errors)
    assert any("bond_length_max_angstrom" in error for error in errors)


def test_reference_2hyy_coordinate_record_is_reproducible():
    path = (
        Path(__file__).resolve().parents[1]
        / "data"
        / "structures"
        / "2HYY.ligand-coordinates.json"
    )
    result = load_ligand_coordinates(path)

    assert validate_ligand_coordinates(result) == []
    assert result.accession == "2HYY"
    assert len(result.atoms) == 37
    assert len(result.bonds) == 41
    assert len(result.angles) == 56
    assert len(result.rings) == 3
    assert result.missing_heavy_atoms == ()
    assert len(result.missing_component_atoms) == 31
    assert result.alternate_location_atoms == ()
    assert status(result, "geometry.bond_lengths_in_range") == "pass"
    assert status(result, "geometry.ring_planarity") == "pass"
    assert status(result, "geometry.close_contacts") == "pass"
    assert any(
        item.category == "missing_hydrogens" and item.status == "not_applicable"
        for item in result.discrepancies
    )
