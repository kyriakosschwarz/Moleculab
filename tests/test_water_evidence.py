"""Tests for deposited-water evidence and geometric bridge plausibility."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import pytest

from moleculab import (
    WATER_EVIDENCE_SCHEMA_VERSION,
    LigandSelection,
    ProteinSelection,
    WaterAtomRecord,
    assess_water_bridges,
    build_water_evidence_record,
    load_water_evidence,
    validate_water_evidence,
    write_water_evidence,
)

CHECKSUM = "6" * 64


def water(
    water_id: str = "A:101:O",
    position: tuple[float, float, float] = (0.0, 0.0, 0.0),
) -> WaterAtomRecord:
    return WaterAtomRecord(
        water_id=water_id,
        label_asym_id="A",
        auth_asym_id="A",
        label_seq_id=None,
        auth_seq_id=water_id.split(":")[1],
        atom_name="O",
        position_angstrom=position,
        occupancy=1.0,
        b_iso_or_equiv=25.0,
        orientation_known=False,
        donor_role=None,
        acceptor_role=None,
        source="deposited HETATM HOH record",
    )


def polar_atom(atom_id: str, position: tuple[float, float, float]):
    return SimpleNamespace(atom_id=atom_id, position_angstrom=position)


def record():
    waters = (water(), water("A:102:O", (10.0, 0.0, 0.0)))
    bridges = assess_water_bridges(
        waters,
        (polar_atom("C1", (3.0, 0.0, 0.0)),),
        (polar_atom("A:1:O", (-3.0, 0.0, 0.0)),),
        3.5,
    )
    return build_water_evidence_record(
        source_manifest_path="manifest.json",
        source_manifest_sha256=CHECKSUM,
        binding_site_record_path="site.json",
        binding_site_record_sha256=CHECKSUM,
        ligand_identity_record_path="identity.json",
        ligand_identity_record_sha256=CHECKSUM,
        chemical_preparation_record_path="preparation.json",
        chemical_preparation_record_sha256=CHECKSUM,
        protein_preparation_record_path="protein.json",
        protein_preparation_record_sha256=CHECKSUM,
        accession="TEST",
        protein_selection=ProteinSelection(label_asym_id="A", auth_asym_id="A"),
        ligand_selection=LigandSelection(
            component_id="LIG",
            label_asym_id="E",
            auth_asym_id="A",
            instance="E",
            auth_seq_id="600",
            available_instances=("E",),
        ),
        ligand_state_id="reference-neutral",
        protein_state_id="standard-templates",
        bridge_cutoff_angstrom=3.5,
        waters=waters,
        bridges=bridges,
    )


def test_bridge_geometry_reports_plausibility_without_chemistry():
    bridges = assess_water_bridges(
        (water(), water("A:102:O", (10.0, 0.0, 0.0))),
        (polar_atom("C1", (3.0, 0.0, 0.0)),),
        (polar_atom("A:1:O", (-3.0, 0.0, 0.0)),),
        3.5,
    )

    assert bridges[0].status == "geometrically_plausible"
    assert bridges[0].chemical_status == "not_assigned"
    assert bridges[0].ligand_distance_angstrom == pytest.approx(3.0)
    assert bridges[0].protein_distance_angstrom == pytest.approx(3.0)
    assert bridges[1].status == "not_supported"
    assert bridges[1].chemical_status == "not_assigned"


def test_missing_polar_evidence_remains_unknown():
    bridges = assess_water_bridges((water(),), (), (), 3.5)

    assert bridges[0].status == "unknown"
    assert bridges[0].ligand_atom_id is None
    assert bridges[0].protein_atom_id is None
    assert bridges[0].chemical_status == "not_assigned"


def test_record_preserves_occupancy_and_b_factor_and_keeps_roles_unassigned():
    result = record()

    assert validate_water_evidence(result) == []
    assert result.schema_version == WATER_EVIDENCE_SCHEMA_VERSION
    assert result.waters[0].occupancy == 1.0
    assert result.waters[0].b_iso_or_equiv == 25.0
    assert result.waters[0].orientation_known is False
    assert result.waters[0].donor_role is None
    assert result.waters[0].acceptor_role is None


def test_record_round_trip(tmp_path):
    result = record()
    path = write_water_evidence(tmp_path / "water.json", result)
    loaded = load_water_evidence(path)

    assert loaded == result
    assert loaded.to_dict()["waters"][0]["water_id"] == "A:101:O"


def test_assigned_water_roles_are_rejected_in_this_batch():
    result = record()
    bad_water = replace(result.waters[0], donor_role="donor")
    bad = replace(result, waters=(bad_water, result.waters[1]))
    errors = validate_water_evidence(bad)

    assert any("roles must be unassigned" in error for error in errors)


def test_chemical_classification_is_rejected_for_geometric_bridges():
    result = record()
    bad_bridge = replace(result.bridges[0], chemical_status="classified")
    bad = replace(result, bridges=(bad_bridge, result.bridges[1]))
    errors = validate_water_evidence(bad)

    assert any("chemical_status must be not_assigned" in error for error in errors)


def test_bridge_must_reference_a_known_water():
    result = record()
    bad_bridge = replace(result.bridges[0], water_id="A:999:O")
    bad = replace(result, bridges=(bad_bridge, result.bridges[1]))
    errors = validate_water_evidence(bad)

    assert any("unknown water" in error for error in errors)


def test_water_record_does_not_contain_affinity_or_ranking():
    serialized = record().to_dict()

    assert "affinity" not in serialized
    assert "score" not in serialized
    assert "rank" not in serialized


def test_reference_2hyy_water_evidence_is_conservative():
    path = (
        Path(__file__).resolve().parents[1]
        / "data"
        / "structures"
        / "2HYY.water-evidence.json"
    )
    result = load_water_evidence(path)

    assert validate_water_evidence(result) == []
    assert result.accession == "2HYY"
    assert len(result.waters) == 4
    assert result.bridge_cutoff_angstrom == 3.5
    assert len(result.bridges) == 4
    assert sum(
        1 for bridge in result.bridges if bridge.status == "geometrically_plausible"
    ) == 1
    assert all(bridge.chemical_status == "not_assigned" for bridge in result.bridges)
    assert all(not water.orientation_known for water in result.waters)
    assert all(water.donor_role is None for water in result.waters)
    assert all(water.acceptor_role is None for water in result.waters)
    assert all(water.occupancy is not None for water in result.waters)
    assert all(water.b_iso_or_equiv is not None for water in result.waters)
