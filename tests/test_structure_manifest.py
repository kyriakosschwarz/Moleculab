# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Tests for the experimental-structure provenance manifest."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from moleculab import (
    SCHEMA_VERSION,
    AnalysisRecord,
    ChemicalIdentityRecord,
    DecisionRecord,
    LigandSelection,
    ProcessingRecord,
    ProteinSelection,
    ProvenanceEntry,
    ProvenanceRecord,
    StructureManifest,
    StructureSelection,
    StructureSource,
    ValidationCheck,
    ValidationRecord,
    load_manifest,
    sha256_file,
    validate_manifest,
    write_manifest,
)

CHECKSUM = "a" * 64


def make_manifest() -> StructureManifest:
    return StructureManifest(
        source=StructureSource(
            database="RCSB PDB",
            accession="2HYY",
            record_url="https://www.rcsb.org/structure/2HYY",
            download_url="https://files.rcsb.org/download/2HYY.cif",
            retrieval_date="2026-11-10",
            revision="1.4",
            revision_date="2023-08-30",
            experimental_method="X-RAY DIFFRACTION",
            resolution_angstrom=2.4,
            coordinate_format="mmCIF",
            file_name="2HYY.cif",
            sha256=CHECKSUM,
            file_size_bytes=899620,
            deposited_asymmetric_unit=True,
        ),
        selection=StructureSelection(
            model_number=1,
            protein=ProteinSelection(label_asym_id="A", auth_asym_id="A"),
            ligand=LigandSelection(
                component_id="STI",
                label_asym_id="E",
                auth_asym_id="A",
                instance="E",
                auth_seq_id="600",
                label_seq_id=None,
                insertion_code=None,
                alternate_location=None,
                available_instances=("E", "F", "G", "H"),
            ),
        ),
        processing=ProcessingRecord(
            alternate_location_policy="No selection applied.",
            hydrogen_addition="None.",
            waters="retained",
            ions_cofactors="retained",
            assembly_id=None,
            assembly_transformations_applied=False,
            original_coordinate_frame_preserved=True,
        ),
        decisions=(
            DecisionRecord(
                decision="selected protein chain A and ligand copy E",
                actor="capstone_structure.py",
                rationale="one chain and one ligand copy are analysed.",
            ),
        ),
        chemical_identity=ChemicalIdentityRecord(
            reference="imatinib",
            component_id="STI",
            name="IMATINIB",
            formula="C29H31N7O",
            formal_charge=None,
            stereochemistry=None,
            atom_mapping=None,
        ),
        analysis=AnalysisRecord(
            coordinate_unit="angstrom",
            internal_distance_unit="nm",
        ),
        provenance=ProvenanceRecord(
            records=(
                ProvenanceEntry(
                    "experimental",
                    "coordinates",
                    "PDB 2HYY",
                    actor="RCSB PDB deposition",
                ),
                ProvenanceEntry(
                    "calculated",
                    "checksum",
                    CHECKSUM,
                    actor="capstone_structure.py",
                ),
            )
        ),
        validation=ValidationRecord(
            checks=(
                ValidationCheck(
                    check_id="source.checksum",
                    category="source_integrity",
                    status="pass",
                    subject="2HYY.cif",
                    criterion="SHA-256 digest is recorded",
                    observed=CHECKSUM,
                    expected=CHECKSUM,
                    detail="digest recorded",
                ),
            )
        ),
    )


def test_sha256_file_matches_the_file_bytes(tmp_path):
    source = tmp_path / "payload.txt"
    source.write_bytes(b"moleculab")

    assert sha256_file(source) == (
        "b9e9baa8624ece900280655f5e9b6e38f9a505410fe81f4f14e0cd4634337841"
    )


def test_manifest_json_round_trip(tmp_path):
    manifest = make_manifest()
    path = write_manifest(tmp_path / "manifest.json", manifest)

    loaded = load_manifest(path)

    assert loaded == manifest
    assert loaded.schema_version == SCHEMA_VERSION
    assert loaded.selection.model_number == 1
    assert loaded.selection.ligand.label_asym_id == "E"
    assert loaded.selection.ligand.auth_asym_id == "A"
    assert loaded.selection.ligand.auth_seq_id == "600"
    assert loaded.selection.ligand.label_seq_id is None
    assert loaded.processing.original_coordinate_frame_preserved is True
    assert loaded.decisions[0].actor == "capstone_structure.py"
    assert loaded.validation.checks[0].check_id == "source.checksum"
    assert loaded.chemical_identity.formal_charge is None


def test_manifest_preserves_explicit_nulls(tmp_path):
    manifest = make_manifest()
    path = write_manifest(tmp_path / "manifest.json", manifest)

    payload = json.loads(path.read_text())

    assert payload["schema_version"] == "1.0"
    assert payload["selection"]["model_number"] == 1
    assert payload["selection"]["ligand"]["label_seq_id"] is None
    assert payload["processing"]["assembly_id"] is None
    assert payload["processing"]["translation_applied"] is None
    assert payload["chemical_identity"]["formal_charge"] is None
    assert payload["chemical_identity"]["stereochemistry"] is None
    assert payload["analysis"]["force_field"] is None
    assert "force_field" in payload["analysis"]


def test_complete_manifest_has_no_validation_errors():
    assert validate_manifest(make_manifest()) == []


def test_reference_2hyy_manifest_records_the_selected_copy():
    path = Path(__file__).resolve().parents[1] / "data" / "structures" / "2HYY.manifest.json"
    manifest = load_manifest(path)

    assert validate_manifest(manifest) == []
    assert manifest.source.accession == "2HYY"
    assert manifest.source.revision == "1.4"
    assert manifest.source.resolution_angstrom == pytest.approx(2.4)
    assert manifest.selection.model_number == 1
    assert manifest.selection.protein.auth_asym_id == "A"
    assert manifest.selection.ligand.component_id == "STI"
    assert manifest.selection.ligand.label_asym_id == "E"
    assert manifest.selection.ligand.auth_asym_id == "A"
    assert manifest.selection.ligand.auth_seq_id == "600"
    assert manifest.selection.ligand.available_instances == ("E", "F", "G", "H")
    assert manifest.processing.original_coordinate_frame_preserved is True
    assert manifest.chemical_identity.formula == "C29H31N7O"
    assert all(check.status != "fail" for check in manifest.validation.checks)


def test_validation_reports_missing_source_integrity():
    manifest = make_manifest()
    manifest = StructureManifest(
        source=StructureSource(
            database="RCSB PDB",
            accession="2HYY",
            record_url=None,
            download_url=None,
            retrieval_date=None,
            revision=None,
            revision_date=None,
            experimental_method=None,
            resolution_angstrom=None,
            coordinate_format="mmCIF",
            file_name=None,
            sha256=None,
            file_size_bytes=None,
            deposited_asymmetric_unit=True,
        ),
        selection=manifest.selection,
        processing=manifest.processing,
        decisions=manifest.decisions,
        chemical_identity=manifest.chemical_identity,
        analysis=manifest.analysis,
        provenance=manifest.provenance,
        validation=manifest.validation,
    )

    errors = validate_manifest(manifest)

    assert any("source.download_url" in error for error in errors)
    assert any("source.retrieval_date" in error for error in errors)
    assert any("source.sha256" in error for error in errors)
    assert any("source.file_size_bytes" in error for error in errors)


def test_validation_requires_explicit_selection_labels_and_decisions():
    manifest = make_manifest()
    manifest = StructureManifest(
        source=manifest.source,
        selection=StructureSelection(
            model_number=0,
            protein=ProteinSelection(label_asym_id=None, auth_asym_id=""),
            ligand=LigandSelection(
                component_id="",
                label_asym_id=None,
                auth_asym_id="",
                instance="",
                auth_seq_id="",
                available_instances=("E",),
            ),
        ),
        processing=manifest.processing,
        decisions=(),
        chemical_identity=manifest.chemical_identity,
        analysis=AnalysisRecord(coordinate_unit="nm", internal_distance_unit="angstrom"),
        provenance=ProvenanceRecord(
            records=(ProvenanceEntry("assumed", "pose", "not labelled"),)
        ),
        validation=ValidationRecord(
            checks=(
                ValidationCheck(
                    check_id="",
                    category="analysis",
                    status="unknown",
                    subject="",
                    criterion="",
                    observed=None,
                    expected=None,
                    detail="",
                ),
            )
        ),
    )

    errors = validate_manifest(manifest)

    assert any("selection.model_number" in error for error in errors)
    assert any("selection.protein.auth_asym_id" in error for error in errors)
    assert any("selection.ligand.instance" in error for error in errors)
    assert any("selection.ligand.auth_seq_id" in error for error in errors)
    assert any("available_instances" in error for error in errors)
    assert any("decisions must contain" in error for error in errors)
    assert any("coordinate_unit" in error for error in errors)
    assert any("internal_distance_unit" in error for error in errors)
    assert any("label must be one of" in error for error in errors)
    assert any("status must be one of" in error for error in errors)
    assert any("check_id" in error for error in errors)


def test_validation_rejects_contradictory_coordinate_frame_and_duplicate_checks():
    manifest = make_manifest()
    processing = ProcessingRecord(
        original_coordinate_frame_preserved=True,
        translation_applied=(1.0, 0.0, 0.0),
    )
    duplicate_check = ValidationCheck(
        check_id="source.checksum",
        category="source_integrity",
        status="pass",
        subject="2HYY.cif",
        criterion="SHA-256 digest is recorded",
        observed=CHECKSUM,
        expected=CHECKSUM,
        detail="digest recorded",
    )
    manifest = StructureManifest(
        source=manifest.source,
        selection=manifest.selection,
        processing=processing,
        decisions=manifest.decisions,
        chemical_identity=manifest.chemical_identity,
        analysis=manifest.analysis,
        provenance=manifest.provenance,
        validation=ValidationRecord(
            checks=(manifest.validation.checks[0], duplicate_check)
        ),
    )

    errors = validate_manifest(manifest)

    assert any("original_coordinate_frame_preserved" in error for error in errors)
    assert any("check_id is duplicated" in error for error in errors)


def test_validation_rejects_a_bad_checksum():
    manifest = make_manifest()
    manifest = StructureManifest(
        source=StructureSource(
            database=manifest.source.database,
            accession=manifest.source.accession,
            record_url=manifest.source.record_url,
            download_url=manifest.source.download_url,
            retrieval_date=manifest.source.retrieval_date,
            revision=manifest.source.revision,
            revision_date=manifest.source.revision_date,
            experimental_method=manifest.source.experimental_method,
            resolution_angstrom=manifest.source.resolution_angstrom,
            coordinate_format=manifest.source.coordinate_format,
            file_name=manifest.source.file_name,
            sha256="not-a-checksum",
            file_size_bytes=manifest.source.file_size_bytes,
            deposited_asymmetric_unit=manifest.source.deposited_asymmetric_unit,
        ),
        selection=manifest.selection,
        processing=manifest.processing,
        decisions=manifest.decisions,
        chemical_identity=manifest.chemical_identity,
        analysis=manifest.analysis,
        provenance=manifest.provenance,
        validation=manifest.validation,
    )

    with pytest.raises(ValueError):
        StructureManifest.from_dict(
            {
                **manifest.to_dict(),
                "selection": {"protein": {}, "ligand": {}},
            }
        )

    assert any("source.sha256" in error for error in validate_manifest(manifest))
