"""Tests for binding-site reconciliation."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import pytest

from moleculab import (
    BINDING_SITE_SCHEMA_VERSION,
    BindingSiteCriteria,
    CoordinateAtom,
    HeteroResidueRecord,
    LigandSelection,
    ProteinAtomRecord,
    ProteinSelection,
    WaterRecord,
    load_binding_site,
    reconcile_binding_site,
    validate_binding_site,
    write_binding_site,
)

CHECKSUM = "d" * 64


def protein_selection() -> ProteinSelection:
    return ProteinSelection(label_asym_id="A", auth_asym_id="A")


def ligand_selection() -> LigandSelection:
    return LigandSelection(
        component_id="LIG",
        label_asym_id="E",
        auth_asym_id="A",
        instance="E",
        auth_seq_id="600",
        available_instances=("E",),
    )


def protein_atom(
    name: str,
    position: tuple[float, float, float],
    *,
    residue_name: str = "ALA",
    auth_seq: str = "10",
    label_seq: str = "1",
    chain: str = "A",
    model: int = 1,
    altloc: str | None = None,
) -> ProteinAtomRecord:
    return ProteinAtomRecord(
        atom_name=name,
        element=name[0],
        residue_name=residue_name,
        label_asym_id=chain,
        auth_asym_id=chain,
        label_seq_id=label_seq,
        auth_seq_id=auth_seq,
        insertion_code=None,
        alternate_location=altloc,
        model_number=model,
        position_angstrom=position,
        occupancy=1.0,
        b_iso_or_equiv=20.0,
    )


def ligand_atom(name: str, position: tuple[float, float, float], index: int) -> CoordinateAtom:
    return CoordinateAtom(
        atom_name=name,
        element=name[0],
        reference_index=index,
        position_angstrom=position,
        occupancy=1.0,
        b_iso_or_equiv=20.0,
        alternate_location=None,
        model_number=1,
    )


def protein_atoms() -> tuple[ProteinAtomRecord, ...]:
    return (
        protein_atom("N", (3.0, 0.0, 0.0), residue_name="LEU", auth_seq="10", label_seq="1"),
        protein_atom("CA", (10.0, 0.0, 0.0), residue_name="LEU", auth_seq="10", label_seq="1"),
        protein_atom("O", (0.0, 3.5, 0.0), residue_name="TYR", auth_seq="11", label_seq="2"),
        protein_atom("CB", (0.0, 5.0, 0.0), residue_name="TYR", auth_seq="11", label_seq="2"),
        protein_atom("N", (8.0, 8.0, 8.0), residue_name="VAL", auth_seq="12", label_seq="3"),
    )


def ligand_atoms() -> tuple[CoordinateAtom, ...]:
    return (
        ligand_atom("C1", (0.0, 0.0, 0.0), 0),
        ligand_atom("C2", (1.4, 0.0, 0.0), 1),
    )


def water(name: str, seq: str, position: tuple[float, float, float]) -> WaterRecord:
    return WaterRecord(
        atom_name="O",
        residue_name="HOH",
        label_asym_id="I",
        auth_asym_id="A",
        label_seq_id=seq,
        auth_seq_id=seq,
        insertion_code=None,
        position_angstrom=position,
        occupancy=1.0,
        b_iso_or_equiv=30.0,
    )


def record():
    return reconcile_binding_site(
        source_manifest_path="manifest.json",
        source_manifest_sha256=CHECKSUM,
        ligand_coordinate_record_path="coordinates.json",
        ligand_coordinate_record_sha256=CHECKSUM,
        accession="TEST",
        model_number=1,
        protein_selection=protein_selection(),
        ligand_selection=ligand_selection(),
        protein_atoms=protein_atoms(),
        ligand_atoms=ligand_atoms(),
        water_atoms=(
            water("O", "101", (0.0, 0.0, 2.0)),
            water("O", "102", (0.0, 0.0, 10.0)),
        ),
        hetero_residues=(
            HeteroResidueRecord(
                residue_name="SO4",
                label_asym_id="A",
                auth_asym_id="A",
                label_seq_id="300",
                auth_seq_id="700",
                insertion_code=None,
                atom_count=5,
            ),
        ),
    )


def status(result, check_id: str) -> str:
    return next(
        check.status for check in result.validation.checks if check.check_id == check_id
    )


def test_complete_binding_site_record_is_valid():
    result = record()

    assert validate_binding_site(result) == []
    assert result.schema_version == BINDING_SITE_SCHEMA_VERSION
    assert result.pocket_atom_count == 2
    assert result.pocket_residue_count == 2
    assert result.protein_atom_count == 5
    assert len(result.contacts) == 4
    assert result.contacts[0].distance_angstrom == 1.6
    assert len(result.waters_near_ligand) == 1
    assert result.water_count_in_selected_chain == 2
    assert len(result.hetero_residues) == 1
    assert result.residue_number_gaps == ()
    assert status(result, "binding_site.pocket_nonempty") == "pass"
    assert status(result, "binding_site.contacts_traceable") == "pass"
    assert status(result, "binding_site.proximity_not_affinity") == "pass"


def test_binding_site_json_round_trip(tmp_path):
    result = record()
    path = write_binding_site(tmp_path / "site.json", result)
    loaded = load_binding_site(path)

    assert loaded == result
    assert loaded.to_dict()["criteria"]["pocket_cutoff_angstrom"] == 4.0


def test_model_mismatch_is_reported():
    atoms = list(protein_atoms())
    atoms[0] = replace(atoms[0], model_number=2)
    result = reconcile_binding_site(
        source_manifest_path="manifest.json",
        source_manifest_sha256=CHECKSUM,
        ligand_coordinate_record_path="coordinates.json",
        ligand_coordinate_record_sha256=CHECKSUM,
        accession="TEST",
        model_number=1,
        protein_selection=protein_selection(),
        ligand_selection=ligand_selection(),
        protein_atoms=tuple(atoms),
        ligand_atoms=ligand_atoms(),
        water_atoms=(),
        hetero_residues=(),
    )

    assert status(result, "binding_site.model_consistent") == "fail"
    assert any(item.category == "model_consistency" for item in result.discrepancies)


def test_chain_mismatch_is_reported():
    atoms = list(protein_atoms())
    atoms[0] = replace(atoms[0], auth_asym_id="B", label_asym_id="B")
    result = reconcile_binding_site(
        source_manifest_path="manifest.json",
        source_manifest_sha256=CHECKSUM,
        ligand_coordinate_record_path="coordinates.json",
        ligand_coordinate_record_sha256=CHECKSUM,
        accession="TEST",
        model_number=1,
        protein_selection=protein_selection(),
        ligand_selection=ligand_selection(),
        protein_atoms=tuple(atoms),
        ligand_atoms=ligand_atoms(),
        water_atoms=(),
        hetero_residues=(),
    )

    assert status(result, "binding_site.protein_chain_matches_manifest") == "fail"
    assert any(item.category == "protein_chain" for item in result.discrepancies)


def test_empty_pocket_is_reported():
    atoms = tuple(
        replace(atom, position_angstrom=(50.0, 50.0, 50.0))
        for atom in protein_atoms()
    )
    result = reconcile_binding_site(
        source_manifest_path="manifest.json",
        source_manifest_sha256=CHECKSUM,
        ligand_coordinate_record_path="coordinates.json",
        ligand_coordinate_record_sha256=CHECKSUM,
        accession="TEST",
        model_number=1,
        protein_selection=protein_selection(),
        ligand_selection=ligand_selection(),
        protein_atoms=atoms,
        ligand_atoms=ligand_atoms(),
        water_atoms=(),
        hetero_residues=(),
    )

    assert result.pocket_atom_count == 0
    assert result.pocket_center_angstrom is None
    assert status(result, "binding_site.pocket_nonempty") == "fail"
    assert any(item.category == "empty_pocket" for item in result.discrepancies)


def test_alternate_locations_are_documented_and_not_selected():
    atoms = list(protein_atoms())
    atoms[0] = replace(atoms[0], alternate_location="A")
    result = reconcile_binding_site(
        source_manifest_path="manifest.json",
        source_manifest_sha256=CHECKSUM,
        ligand_coordinate_record_path="coordinates.json",
        ligand_coordinate_record_sha256=CHECKSUM,
        accession="TEST",
        model_number=1,
        protein_selection=protein_selection(),
        ligand_selection=ligand_selection(),
        protein_atoms=tuple(atoms),
        ligand_atoms=ligand_atoms(),
        water_atoms=(),
        hetero_residues=(),
    )

    assert result.alternate_location_atom_count == 1
    assert status(result, "binding_site.alternate_locations_documented") == "pass"
    assert any(item.category == "alternate_locations" for item in result.discrepancies)


def test_residue_number_gaps_are_reported():
    atoms = list(protein_atoms())
    atoms[-1] = replace(atoms[-1], auth_seq_id="14", label_seq_id="5")
    result = reconcile_binding_site(
        source_manifest_path="manifest.json",
        source_manifest_sha256=CHECKSUM,
        ligand_coordinate_record_path="coordinates.json",
        ligand_coordinate_record_sha256=CHECKSUM,
        accession="TEST",
        model_number=1,
        protein_selection=protein_selection(),
        ligand_selection=ligand_selection(),
        protein_atoms=tuple(atoms),
        ligand_atoms=ligand_atoms(),
        water_atoms=(),
        hetero_residues=(),
    )

    assert "13" in result.residue_number_gaps
    assert status(result, "binding_site.missing_residues") == "not_applicable"
    assert any(item.category == "missing_residues" for item in result.discrepancies)


def test_custom_cutoff_changes_the_pocket():
    result = reconcile_binding_site(
        source_manifest_path="manifest.json",
        source_manifest_sha256=CHECKSUM,
        ligand_coordinate_record_path="coordinates.json",
        ligand_coordinate_record_sha256=CHECKSUM,
        accession="TEST",
        model_number=1,
        protein_selection=protein_selection(),
        ligand_selection=ligand_selection(),
        protein_atoms=protein_atoms(),
        ligand_atoms=ligand_atoms(),
        water_atoms=(),
        hetero_residues=(),
        criteria=BindingSiteCriteria(
            pocket_cutoff_angstrom=2.0,
            pocket_cutoff_nm=0.2,
            water_cutoff_angstrom=6.0,
        ),
    )

    assert result.pocket_atom_count == 1
    assert len(result.contacts) == 1


def test_validation_rejects_bad_site_record():
    result = record()
    bad = replace(
        result,
        pocket_residue_count=999,
        criteria=BindingSiteCriteria(
            pocket_cutoff_angstrom=4.0,
            pocket_cutoff_nm=0.5,
        ),
    )

    errors = validate_binding_site(bad)

    assert any("pocket_residue_count" in error for error in errors)
    assert any("pocket_cutoff_angstrom" in error for error in errors)


def test_reference_2hyy_binding_site_record_is_reproducible():
    path = (
        Path(__file__).resolve().parents[1]
        / "data"
        / "structures"
        / "2HYY.binding-site.json"
    )
    result = load_binding_site(path)

    assert validate_binding_site(result) == []
    assert result.accession == "2HYY"
    assert result.model_number == 1
    assert result.protein_selection.auth_asym_id == "A"
    assert result.ligand_selection.component_id == "STI"
    assert result.pocket_atom_count == 53
    assert result.pocket_residue_count == 23
    assert len(result.contacts) == 112
    assert result.contacts[0].distance_angstrom == pytest.approx(2.723, abs=0.001)
    assert len(result.waters_near_ligand) == 4
    assert result.residue_number_gaps == ("275",)
    assert result.alternate_location_atom_count == 0
    assert status(result, "binding_site.proximity_not_affinity") == "pass"
