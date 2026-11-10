# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Tests for separate chemical-state preparation records."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

from moleculab import (
    CHEMICAL_PREPARATION_SCHEMA_VERSION,
    LigandSelection,
    PreparedAtom,
    PreparedRing,
    PreparedState,
    PreparedWaterRole,
    build_chemical_preparation_record,
    load_chemical_preparation,
    validate_chemical_preparation,
    write_chemical_preparation,
)

CHECKSUM = "f" * 64


def selection() -> LigandSelection:
    return LigandSelection(
        component_id="LIG",
        label_asym_id="E",
        auth_asym_id="A",
        instance="E",
        auth_seq_id="600",
        available_instances=("E",),
    )


def atom(
    atom_id: str,
    *,
    formal_charge: int | None = None,
    donor: bool | None = None,
    acceptor: bool | None = None,
    hydrophobic: bool | None = None,
    hydrogens: tuple[tuple[float, float, float], ...] = (),
    coordinate_source: str = "deposited",
    hydrogen_source: str | None = None,
    feature_source: str | None = None,
) -> PreparedAtom:
    return PreparedAtom(
        atom_id=atom_id,
        element="C",
        formal_charge=formal_charge,
        is_donor=donor,
        is_acceptor=acceptor,
        is_hydrophobic=hydrophobic,
        aromatic_ring_ids=(),
        hydrogen_count=len(hydrogens),
        hydrogen_positions_nm=hydrogens,
        coordinate_source=coordinate_source,
        hydrogen_source=hydrogen_source,
        feature_source=feature_source,
    )


def state(state_id: str = "deposited-only") -> PreparedState:
    return PreparedState(
        state_id=state_id,
        label="test state",
        protonation_state="unknown",
        tautomer="unknown",
        formal_charge=None,
        hydrogen_method=None,
        feature_source=None,
        assumptions=("no chemistry was assigned",),
        atoms=(atom("C1"), atom("O1")),
        rings=(),
        water_roles=(
            PreparedWaterRole("W1", "unassigned", "experimental water"),
        ),
    )


def record():
    return build_chemical_preparation_record(
        source_manifest_path="manifest.json",
        source_manifest_sha256=CHECKSUM,
        ligand_identity_record_path="identity.json",
        ligand_identity_record_sha256=CHECKSUM,
        ligand_coordinate_record_path="coordinates.json",
        ligand_coordinate_record_sha256=CHECKSUM,
        binding_site_record_path="site.json",
        binding_site_record_sha256=CHECKSUM,
        accession="TEST",
        selection=selection(),
        toolkit="RDKit",
        toolkit_version="test",
        parameters=("test parameters",),
        states=(state("deposited-only"), state("reference-neutral")),
        selected_state_id=None,
    )


def test_complete_preparation_record_is_valid():
    result = record()

    assert validate_chemical_preparation(result) == []
    assert result.schema_version == CHEMICAL_PREPARATION_SCHEMA_VERSION
    assert len(result.states) == 2
    assert result.selected_state_id is None
    assert any(item.category == "state_selection" for item in result.discrepancies)


def test_preparation_json_round_trip(tmp_path):
    result = record()
    path = write_chemical_preparation(tmp_path / "preparation.json", result)
    loaded = load_chemical_preparation(path)

    assert loaded == result
    assert loaded.to_dict()["states"][0]["state_id"] == "deposited-only"


def test_selected_state_must_exist():
    result = record()
    bad = replace(result, selected_state_id="missing")
    errors = validate_chemical_preparation(bad)

    assert any("selected_state_id" in error for error in errors)


def test_duplicate_state_ids_are_rejected():
    result = record()
    bad = replace(result, states=(state("same"), state("same")))
    errors = validate_chemical_preparation(bad)

    assert any("state_id values must be unique" in error for error in errors)


def test_hydrogen_count_must_match_positions():
    result = record()
    bad_atom = replace(
        result.states[0].atoms[0],
        hydrogen_count=2,
        hydrogen_positions_nm=((0.0, 0.0, 0.0),),
    )
    bad_state = replace(result.states[0], atoms=(bad_atom,))
    bad = replace(result, states=(bad_state, result.states[1]))
    errors = validate_chemical_preparation(bad)

    assert any("hydrogen count does not match" in error for error in errors)


def test_ring_must_reference_known_atoms():
    result = record()
    bad_ring = PreparedRing(
        ring_id="r1",
        atom_ids=("MISSING",),
        centroid_nm=(0.0, 0.0, 0.0),
        normal=(0.0, 0.0, 1.0),
        source="test",
    )
    bad_state = replace(result.states[0], rings=(bad_ring,))
    bad = replace(result, states=(bad_state, result.states[1]))
    errors = validate_chemical_preparation(bad)

    assert any("references an unknown atom" in error for error in errors)


def test_water_role_must_be_explicit():
    result = record()
    bad_water = PreparedWaterRole("W1", "assumed", "test")
    bad_state = replace(result.states[0], water_roles=(bad_water,))
    bad = replace(result, states=(bad_state, result.states[1]))
    errors = validate_chemical_preparation(bad)

    assert any("role is invalid" in error for error in errors)


def test_reference_2hyy_preparation_record_is_reproducible():
    path = (
        Path(__file__).resolve().parents[1]
        / "data"
        / "structures"
        / "2HYY.chemical-preparation.json"
    )
    result = load_chemical_preparation(path)

    assert validate_chemical_preparation(result) == []
    assert result.accession == "2HYY"
    assert len(result.states) == 2
    assert result.selected_state_id is None
    deposited = next(state for state in result.states if state.state_id == "deposited-only")
    reference = next(
        state for state in result.states if state.state_id == "reference-neutral"
    )
    assert all(atom.formal_charge is None for atom in deposited.atoms)
    assert all(atom.is_donor is None for atom in deposited.atoms)
    assert sum(atom.hydrogen_count for atom in reference.atoms) == 31
    assert all(atom.coordinate_source == "deposited" for atom in reference.atoms)
    assert any(atom.hydrogen_source for atom in reference.atoms)
    assert len(reference.rings) == 4
    assert sum(1 for atom in reference.atoms if atom.is_donor) == 4
    assert sum(1 for atom in reference.atoms if atom.is_acceptor) == 4
    assert all(role.role == "unassigned" for role in reference.water_roles)
    assert result.toolkit == "RDKit"
