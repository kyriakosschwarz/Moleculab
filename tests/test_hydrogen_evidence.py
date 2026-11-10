"""Tests for hydrogen-position and hydrogen-bond evidence records."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import pytest

from moleculab import (
    HYDROGEN_EVIDENCE_SCHEMA_VERSION,
    HydrogenBondCriteria,
    HydrogenBondEvidence,
    HydrogenPosition,
    LigandSelection,
    ProteinSelection,
    assess_ligand_donor_hydrogen_bonds,
    build_hydrogen_evidence_record,
    load_hydrogen_evidence,
    validate_hydrogen_evidence,
    write_hydrogen_evidence,
)

CHECKSUM = "5" * 64


def position(
    hydrogen_id: str = "N1:H1",
    *,
    hydrogen: tuple[float, float, float] = (1.0, 0.0, 0.0),
    parent: tuple[float, float, float] = (0.0, 0.0, 0.0),
) -> HydrogenPosition:
    distance = sum((a - b) ** 2 for a, b in zip(hydrogen, parent, strict=True)) ** 0.5
    return HydrogenPosition(
        hydrogen_id=hydrogen_id,
        parent_atom_id="N1",
        parent_side="ligand",
        parent_element="N",
        position_angstrom=hydrogen,
        parent_position_angstrom=parent,
        donor_hydrogen_distance_angstrom=distance,
        state_id="reference-neutral",
        generation_method="test AddHs",
        toolkit="test toolkit",
        toolkit_version="1.0",
        coordinate_source="generated",
        parent_coordinate_source="deposited",
        feature_source="test features",
    )


def acceptor(
    atom_id: str = "A:1:O",
    position_angstrom: tuple[float, float, float] = (2.0, 0.0, 0.0),
):
    return SimpleNamespace(atom_id=atom_id, position_angstrom=position_angstrom)


def record():
    criteria = HydrogenBondCriteria(max_distance=0.35, min_angle_degrees=120.0)
    hydrogens = (position(),)
    bonds = assess_ligand_donor_hydrogen_bonds(
        hydrogens,
        (acceptor(),),
        criteria,
    )
    return build_hydrogen_evidence_record(
        source_manifest_path="manifest.json",
        source_manifest_sha256=CHECKSUM,
        ligand_identity_record_path="identity.json",
        ligand_identity_record_sha256=CHECKSUM,
        ligand_coordinate_record_path="coordinates.json",
        ligand_coordinate_record_sha256=CHECKSUM,
        binding_site_record_path="site.json",
        binding_site_record_sha256=CHECKSUM,
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
        criteria=criteria,
        ligand_hydrogens=hydrogens,
        protein_hydrogen_status="not_run",
        bonds=bonds,
    )


def test_aligned_donor_hydrogen_acceptor_geometry_is_supported():
    criteria = HydrogenBondCriteria(max_distance=0.35, min_angle_degrees=120.0)
    bonds = assess_ligand_donor_hydrogen_bonds(
        (position(),),
        (acceptor(),),
        criteria,
    )

    assert len(bonds) == 1
    assert bonds[0].status == "supported"
    assert bonds[0].angle_degrees == pytest.approx(180.0)
    assert bonds[0].donor_acceptor_distance_angstrom == pytest.approx(2.0)
    assert bonds[0].to_dict()["donor_acceptor_distance_nm"] == pytest.approx(0.2)


def test_poor_angular_geometry_is_not_supported():
    criteria = HydrogenBondCriteria(max_distance=0.35, min_angle_degrees=120.0)
    bonds = assess_ligand_donor_hydrogen_bonds(
        (position(),),
        (acceptor(position_angstrom=(0.0, 2.0, 0.0)),),
        criteria,
    )

    assert bonds[0].status == "not_supported"
    assert bonds[0].angle_degrees < 120.0


def test_distance_cutoff_excludes_distant_acceptors():
    criteria = HydrogenBondCriteria(max_distance=0.35, min_angle_degrees=120.0)
    bonds = assess_ligand_donor_hydrogen_bonds(
        (position(),),
        (acceptor(position_angstrom=(10.0, 0.0, 0.0)),),
        criteria,
    )

    assert bonds == ()


def test_record_round_trip_and_validation(tmp_path):
    result = record()
    path = write_hydrogen_evidence(tmp_path / "hydrogen.json", result)
    loaded = load_hydrogen_evidence(path)

    assert loaded == result
    assert validate_hydrogen_evidence(loaded) == []
    assert loaded.schema_version == HYDROGEN_EVIDENCE_SCHEMA_VERSION
    assert loaded.protein_hydrogen_status == "not_run"


def test_generated_and_deposited_provenance_must_remain_distinct():
    result = record()
    bad_hydrogen = replace(result.ligand_hydrogens[0], coordinate_source="deposited")
    bad = replace(result, ligand_hydrogens=(bad_hydrogen,))
    errors = validate_hydrogen_evidence(bad)

    assert any("must be generated" in error for error in errors)


def test_hydrogen_record_does_not_contain_affinity_or_ranking():
    serialized = record().to_dict()

    assert "affinity" not in serialized
    assert "score" not in serialized
    assert "rank" not in serialized


def test_reference_2hyy_hydrogen_evidence_is_conservative():
    path = (
        Path(__file__).resolve().parents[1]
        / "data"
        / "structures"
        / "2HYY.hydrogen-evidence.json"
    )
    result = load_hydrogen_evidence(path)

    assert validate_hydrogen_evidence(result) == []
    assert result.accession == "2HYY"
    assert result.ligand_state_id == "reference-neutral"
    assert result.protein_state_id == "standard-templates"
    assert len(result.ligand_hydrogens) == 2
    assert result.protein_hydrogen_status == "not_run"
    assert len(result.bonds) == 3
    assert sum(1 for bond in result.bonds if bond.status == "supported") == 2
    assert all(
        hydrogen.coordinate_source == "generated"
        and hydrogen.parent_coordinate_source == "deposited"
        for hydrogen in result.ligand_hydrogens
    )
    assert all(bond.donor_side == "ligand" for bond in result.bonds)
    assert all(bond.acceptor_side == "protein" for bond in result.bonds)
    assert all(isinstance(bond, HydrogenBondEvidence) for bond in result.bonds)
