"""Provenance records for experimental structure ingestion.

The manifest is the contract for the capstone workflow.  It records which
experimental file was used, which protein chain and ligand instance were
selected, what processing happened, and which claims are experimental,
calculated, predicted, or hypothesised.

This module deliberately does not download structures, parse coordinates, or
choose a chain.  Those are explicit caller decisions.  It provides immutable
records, deterministic JSON serialisation, and validation of the record itself.
"""

from __future__ import annotations

import hashlib
import json
import re
from collections.abc import Mapping
from dataclasses import dataclass
from datetime import date
from pathlib import Path
from typing import Any

SCHEMA_VERSION = "1.0"
PROVENANCE_LABELS = ("experimental", "calculated", "predicted", "hypothesised")
COORDINATE_FORMATS = ("mmCIF", "PDB")
VALIDATION_STATUSES = ("pass", "fail", "not_run")

_SHA256_PATTERN = re.compile(r"[0-9a-f]{64}")


def sha256_file(path: str | Path, chunk_size: int = 1 << 20) -> str:
    """Return the SHA-256 checksum of *path* without altering the file."""

    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _required_string(data: Mapping[str, Any], key: str) -> str:
    value = data.get(key)
    if value is None or str(value) == "":
        raise ValueError(f"{key} is required.")
    return str(value)


def _optional_string(data: Mapping[str, Any], key: str) -> str | None:
    value = data.get(key)
    if value is None:
        return None
    return str(value)


def _optional_float(data: Mapping[str, Any], key: str) -> float | None:
    value = data.get(key)
    if value is None:
        return None
    return float(value)


def _required_int(data: Mapping[str, Any], key: str) -> int:
    value = data.get(key)
    if value is None:
        raise ValueError(f"{key} is required.")
    return int(value)


def _required_bool(data: Mapping[str, Any], key: str) -> bool:
    value = data.get(key)
    if not isinstance(value, bool):
        raise ValueError(f"{key} must be a boolean.")
    return value


def _optional_int(data: Mapping[str, Any], key: str) -> int | None:
    value = data.get(key)
    if value is None:
        return None
    return int(value)


def _optional_float_tuple(
    data: Mapping[str, Any],
    key: str,
    *,
    length: int | None = None,
) -> tuple[float, ...] | None:
    value = data.get(key)
    if value is None:
        return None
    if not isinstance(value, (list, tuple)):
        raise ValueError(f"{key} must be a sequence of numbers.")
    result = tuple(float(item) for item in value)
    if length is not None and len(result) != length:
        raise ValueError(f"{key} must contain {length} numbers.")
    return result


def _as_mapping(value: Any, name: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must be a mapping.")
    return value


def _as_list(value: Any, name: str) -> list[Any]:
    if value is None:
        return []
    if not isinstance(value, (list, tuple)):
        raise ValueError(f"{name} must be a list.")
    return list(value)


@dataclass(frozen=True)
class StructureSource:
    """Identity and integrity of the deposited coordinate file."""

    database: str
    accession: str
    record_url: str | None
    download_url: str | None
    retrieval_date: str | None
    revision: str | None
    revision_date: str | None
    experimental_method: str | None
    resolution_angstrom: float | None
    coordinate_format: str
    file_name: str | None
    sha256: str | None
    file_size_bytes: int | None
    deposited_asymmetric_unit: bool
    biological_assembly: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "database": self.database,
            "accession": self.accession,
            "record_url": self.record_url,
            "download_url": self.download_url,
            "retrieval_date": self.retrieval_date,
            "revision": self.revision,
            "revision_date": self.revision_date,
            "experimental_method": self.experimental_method,
            "resolution_angstrom": self.resolution_angstrom,
            "coordinate_format": self.coordinate_format,
            "file_name": self.file_name,
            "sha256": self.sha256,
            "file_size_bytes": self.file_size_bytes,
            "deposited_asymmetric_unit": self.deposited_asymmetric_unit,
            "biological_assembly": self.biological_assembly,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> StructureSource:
        data = _as_mapping(data, "source")
        return cls(
            database=_required_string(data, "database"),
            accession=_required_string(data, "accession"),
            record_url=_optional_string(data, "record_url"),
            download_url=_optional_string(data, "download_url"),
            retrieval_date=_optional_string(data, "retrieval_date"),
            revision=_optional_string(data, "revision"),
            revision_date=_optional_string(data, "revision_date"),
            experimental_method=_optional_string(data, "experimental_method"),
            resolution_angstrom=_optional_float(data, "resolution_angstrom"),
            coordinate_format=_required_string(data, "coordinate_format"),
            file_name=_optional_string(data, "file_name"),
            sha256=_optional_string(data, "sha256"),
            file_size_bytes=_optional_int(data, "file_size_bytes"),
            deposited_asymmetric_unit=_required_bool(data, "deposited_asymmetric_unit"),
            biological_assembly=_optional_string(data, "biological_assembly"),
        )


@dataclass(frozen=True)
class ProteinSelection:
    """One explicitly selected protein chain."""

    label_asym_id: str | None
    auth_asym_id: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "label_asym_id": self.label_asym_id,
            "auth_asym_id": self.auth_asym_id,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> ProteinSelection:
        data = _as_mapping(data, "selection.protein")
        return cls(
            label_asym_id=_optional_string(data, "label_asym_id"),
            auth_asym_id=_required_string(data, "auth_asym_id"),
        )


@dataclass(frozen=True)
class LigandSelection:
    """One explicitly selected ligand instance.

    ``label_asym_id`` is the mmCIF label asymmetry identifier.  It is not
    interchangeable with ``auth_asym_id``, the author chain identifier, and it
    is not the same thing as the ligand component identifier.
    """

    component_id: str
    label_asym_id: str | None
    auth_asym_id: str
    instance: str
    auth_seq_id: str
    label_seq_id: str | None = None
    insertion_code: str | None = None
    alternate_location: str | None = None
    available_instances: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        return {
            "component_id": self.component_id,
            "label_asym_id": self.label_asym_id,
            "auth_asym_id": self.auth_asym_id,
            "instance": self.instance,
            "auth_seq_id": self.auth_seq_id,
            "label_seq_id": self.label_seq_id,
            "insertion_code": self.insertion_code,
            "alternate_location": self.alternate_location,
            "available_instances": list(self.available_instances),
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> LigandSelection:
        data = _as_mapping(data, "selection.ligand")
        return cls(
            component_id=_required_string(data, "component_id"),
            label_asym_id=_optional_string(data, "label_asym_id"),
            auth_asym_id=_required_string(data, "auth_asym_id"),
            instance=_required_string(data, "instance"),
            auth_seq_id=_required_string(data, "auth_seq_id"),
            label_seq_id=_optional_string(data, "label_seq_id"),
            insertion_code=_optional_string(data, "insertion_code"),
            alternate_location=_optional_string(data, "alternate_location"),
            available_instances=tuple(
                str(value) for value in _as_list(data.get("available_instances"), "available_instances")
            ),
        )


@dataclass(frozen=True)
class StructureSelection:
    model_number: int
    protein: ProteinSelection
    ligand: LigandSelection

    def to_dict(self) -> dict[str, Any]:
        return {
            "model_number": self.model_number,
            "protein": self.protein.to_dict(),
            "ligand": self.ligand.to_dict(),
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> StructureSelection:
        data = _as_mapping(data, "selection")
        return cls(
            model_number=_required_int(data, "model_number"),
            protein=ProteinSelection.from_dict(_as_mapping(data.get("protein"), "selection.protein")),
            ligand=LigandSelection.from_dict(_as_mapping(data.get("ligand"), "selection.ligand")),
        )


@dataclass(frozen=True)
class ProcessingRecord:
    """Explicit record of changes, or the absence of changes."""

    coordinate_transforms: tuple[str, ...] = ()
    removed_waters: tuple[str, ...] = ()
    removed_ions: tuple[str, ...] = ()
    removed_cofactors: tuple[str, ...] = ()
    alternate_location_policy: str | None = None
    hydrogen_addition: str | None = None
    waters: str = "retained"
    ions_cofactors: str = "retained"
    missing_atoms_or_residues: str | None = None
    assembly_id: str | None = None
    assembly_transformations_applied: bool = False
    original_coordinate_frame_preserved: bool = True
    translation_applied: tuple[float, ...] | None = None
    rotation_applied: tuple[float, ...] | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "coordinate_transforms": list(self.coordinate_transforms),
            "removed_waters": list(self.removed_waters),
            "removed_ions": list(self.removed_ions),
            "removed_cofactors": list(self.removed_cofactors),
            "alternate_location_policy": self.alternate_location_policy,
            "hydrogen_addition": self.hydrogen_addition,
            "waters": self.waters,
            "ions_cofactors": self.ions_cofactors,
            "missing_atoms_or_residues": self.missing_atoms_or_residues,
            "assembly_id": self.assembly_id,
            "assembly_transformations_applied": self.assembly_transformations_applied,
            "original_coordinate_frame_preserved": self.original_coordinate_frame_preserved,
            "translation_applied": (
                None if self.translation_applied is None else list(self.translation_applied)
            ),
            "rotation_applied": (
                None if self.rotation_applied is None else list(self.rotation_applied)
            ),
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> ProcessingRecord:
        data = _as_mapping(data, "processing")
        return cls(
            coordinate_transforms=tuple(str(value) for value in _as_list(data.get("coordinate_transforms"), "coordinate_transforms")),
            removed_waters=tuple(str(value) for value in _as_list(data.get("removed_waters"), "removed_waters")),
            removed_ions=tuple(str(value) for value in _as_list(data.get("removed_ions"), "removed_ions")),
            removed_cofactors=tuple(str(value) for value in _as_list(data.get("removed_cofactors"), "removed_cofactors")),
            alternate_location_policy=_optional_string(data, "alternate_location_policy"),
            hydrogen_addition=_optional_string(data, "hydrogen_addition"),
            waters=_required_string(data, "waters"),
            ions_cofactors=_required_string(data, "ions_cofactors"),
            missing_atoms_or_residues=_optional_string(data, "missing_atoms_or_residues"),
            assembly_id=_optional_string(data, "assembly_id"),
            assembly_transformations_applied=_required_bool(
                data, "assembly_transformations_applied"
            ),
            original_coordinate_frame_preserved=_required_bool(
                data, "original_coordinate_frame_preserved"
            ),
            translation_applied=_optional_float_tuple(data, "translation_applied", length=3),
            rotation_applied=_optional_float_tuple(data, "rotation_applied", length=9),
        )


@dataclass(frozen=True)
class ChemicalIdentityRecord:
    reference: str
    component_id: str | None
    name: str | None
    formula: str | None
    formal_charge: int | None
    stereochemistry: str | None
    atom_mapping: str | None
    notes: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "reference": self.reference,
            "component_id": self.component_id,
            "name": self.name,
            "formula": self.formula,
            "formal_charge": self.formal_charge,
            "stereochemistry": self.stereochemistry,
            "atom_mapping": self.atom_mapping,
            "notes": self.notes,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> ChemicalIdentityRecord:
        data = _as_mapping(data, "chemical_identity")
        return cls(
            reference=_required_string(data, "reference"),
            component_id=_optional_string(data, "component_id"),
            name=_optional_string(data, "name"),
            formula=_optional_string(data, "formula"),
            formal_charge=_optional_int(data, "formal_charge"),
            stereochemistry=_optional_string(data, "stereochemistry"),
            atom_mapping=_optional_string(data, "atom_mapping"),
            notes=_optional_string(data, "notes"),
        )


@dataclass(frozen=True)
class AnalysisRecord:
    coordinate_unit: str
    internal_distance_unit: str
    force_field: str | None = None
    charge_model: str | None = None
    interaction_criteria: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "coordinate_unit": self.coordinate_unit,
            "internal_distance_unit": self.internal_distance_unit,
            "force_field": self.force_field,
            "charge_model": self.charge_model,
            "interaction_criteria": self.interaction_criteria,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> AnalysisRecord:
        data = _as_mapping(data, "analysis")
        return cls(
            coordinate_unit=_required_string(data, "coordinate_unit"),
            internal_distance_unit=_required_string(data, "internal_distance_unit"),
            force_field=_optional_string(data, "force_field"),
            charge_model=_optional_string(data, "charge_model"),
            interaction_criteria=_optional_string(data, "interaction_criteria"),
        )


@dataclass(frozen=True)
class DecisionRecord:
    """A selection or processing decision, separate from source evidence."""

    decision: str
    actor: str
    rationale: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "decision": self.decision,
            "actor": self.actor,
            "rationale": self.rationale,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> DecisionRecord:
        data = _as_mapping(data, "decision record")
        return cls(
            decision=_required_string(data, "decision"),
            actor=_required_string(data, "actor"),
            rationale=_required_string(data, "rationale"),
        )


@dataclass(frozen=True)
class ProvenanceEntry:
    label: str
    subject: str
    detail: str
    actor: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "label": self.label,
            "subject": self.subject,
            "detail": self.detail,
            "actor": self.actor,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> ProvenanceEntry:
        data = _as_mapping(data, "provenance record")
        return cls(
            label=_required_string(data, "label"),
            subject=_required_string(data, "subject"),
            detail=_required_string(data, "detail"),
            actor=_optional_string(data, "actor"),
        )


@dataclass(frozen=True)
class ProvenanceRecord:
    records: tuple[ProvenanceEntry, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        return {"records": [record.to_dict() for record in self.records]}

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> ProvenanceRecord:
        data = _as_mapping(data, "provenance")
        return cls(
            records=tuple(
                ProvenanceEntry.from_dict(_as_mapping(item, "provenance record"))
                for item in _as_list(data.get("records"), "provenance.records")
            )
        )


@dataclass(frozen=True)
class ValidationCheck:
    check_id: str
    category: str
    status: str
    subject: str
    criterion: str
    observed: str | None
    expected: str | None
    detail: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "check_id": self.check_id,
            "category": self.category,
            "status": self.status,
            "subject": self.subject,
            "criterion": self.criterion,
            "observed": self.observed,
            "expected": self.expected,
            "detail": self.detail,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> ValidationCheck:
        data = _as_mapping(data, "validation check")
        return cls(
            check_id=_required_string(data, "check_id"),
            category=_required_string(data, "category"),
            status=_required_string(data, "status"),
            subject=_required_string(data, "subject"),
            criterion=_required_string(data, "criterion"),
            observed=_optional_string(data, "observed"),
            expected=_optional_string(data, "expected"),
            detail=_required_string(data, "detail"),
        )


@dataclass(frozen=True)
class ValidationRecord:
    checks: tuple[ValidationCheck, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        return {"checks": [check.to_dict() for check in self.checks]}

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> ValidationRecord:
        data = _as_mapping(data, "validation")
        return cls(
            checks=tuple(
                ValidationCheck.from_dict(_as_mapping(item, "validation check"))
                for item in _as_list(data.get("checks"), "validation.checks")
            )
        )


@dataclass(frozen=True)
class StructureManifest:
    source: StructureSource
    selection: StructureSelection
    processing: ProcessingRecord
    decisions: tuple[DecisionRecord, ...]
    chemical_identity: ChemicalIdentityRecord
    analysis: AnalysisRecord
    provenance: ProvenanceRecord
    validation: ValidationRecord
    schema_version: str = SCHEMA_VERSION

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "source": self.source.to_dict(),
            "selection": self.selection.to_dict(),
            "processing": self.processing.to_dict(),
            "decisions": [decision.to_dict() for decision in self.decisions],
            "chemical_identity": self.chemical_identity.to_dict(),
            "analysis": self.analysis.to_dict(),
            "provenance": self.provenance.to_dict(),
            "validation": self.validation.to_dict(),
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> StructureManifest:
        data = _as_mapping(data, "manifest")
        return cls(
            schema_version=_required_string(data, "schema_version"),
            source=StructureSource.from_dict(_as_mapping(data.get("source"), "source")),
            selection=StructureSelection.from_dict(_as_mapping(data.get("selection"), "selection")),
            processing=ProcessingRecord.from_dict(_as_mapping(data.get("processing"), "processing")),
            decisions=tuple(
                DecisionRecord.from_dict(_as_mapping(item, "decision record"))
                for item in _as_list(data.get("decisions"), "decisions")
            ),
            chemical_identity=ChemicalIdentityRecord.from_dict(
                _as_mapping(data.get("chemical_identity"), "chemical_identity")
            ),
            analysis=AnalysisRecord.from_dict(_as_mapping(data.get("analysis"), "analysis")),
            provenance=ProvenanceRecord.from_dict(_as_mapping(data.get("provenance"), "provenance")),
            validation=ValidationRecord.from_dict(_as_mapping(data.get("validation"), "validation")),
        )


def write_manifest(path: str | Path, manifest: StructureManifest) -> Path:
    """Write a deterministic JSON manifest and return its path."""

    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(manifest.to_dict(), indent=2, sort_keys=False) + "\n",
        encoding="utf-8",
    )
    return destination


def load_manifest(path: str | Path) -> StructureManifest:
    """Read a JSON manifest and reconstruct its immutable records."""

    data = json.loads(Path(path).read_text(encoding="utf-8"))
    return StructureManifest.from_dict(data)


def validate_manifest(manifest: StructureManifest) -> list[str]:
    """Return schema and completeness errors for *manifest*.

    An empty list means the record is complete enough to be used as a
    provenance contract.  A validation check that has not run belongs in the
    manifest's ``validation`` record, not in this function's return value.
    """

    errors: list[str] = []
    if manifest.schema_version != SCHEMA_VERSION:
        errors.append(
            f"schema_version must be {SCHEMA_VERSION!r}, got {manifest.schema_version!r}."
        )

    source = manifest.source
    if not source.database:
        errors.append("source.database is required.")
    if not source.accession:
        errors.append("source.accession is required.")
    if not source.record_url:
        errors.append("source.record_url is required.")
    if not source.download_url:
        errors.append("source.download_url is required.")
    if not source.retrieval_date:
        errors.append("source.retrieval_date is required.")
    else:
        try:
            date.fromisoformat(source.retrieval_date)
        except ValueError:
            errors.append("source.retrieval_date must be an ISO date (YYYY-MM-DD).")
    if not source.revision:
        errors.append("source.revision is required.")
    if not source.revision_date:
        errors.append("source.revision_date is required.")
    else:
        try:
            date.fromisoformat(source.revision_date)
        except ValueError:
            errors.append("source.revision_date must be an ISO date (YYYY-MM-DD).")
    if source.coordinate_format not in COORDINATE_FORMATS:
        errors.append(
            f"source.coordinate_format must be one of {COORDINATE_FORMATS}, got "
            f"{source.coordinate_format!r}."
        )
    if not source.file_name:
        errors.append("source.file_name is required.")
    if source.sha256 is None or not _SHA256_PATTERN.fullmatch(source.sha256):
        errors.append("source.sha256 must be a lowercase 64-character SHA-256 digest.")
    if source.file_size_bytes is None or source.file_size_bytes <= 0:
        errors.append("source.file_size_bytes must be positive.")
    if source.resolution_angstrom is not None and source.resolution_angstrom <= 0:
        errors.append("source.resolution_angstrom must be positive when provided.")
    if not isinstance(source.deposited_asymmetric_unit, bool):
        errors.append("source.deposited_asymmetric_unit must be a boolean.")

    selection = manifest.selection
    if selection.model_number <= 0:
        errors.append("selection.model_number must be a positive integer.")
    if not selection.protein.auth_asym_id:
        errors.append("selection.protein.auth_asym_id is required.")
    if not selection.ligand.component_id:
        errors.append("selection.ligand.component_id is required.")
    if not selection.ligand.auth_asym_id:
        errors.append("selection.ligand.auth_asym_id is required.")
    if not selection.ligand.instance:
        errors.append("selection.ligand.instance is required.")
    if not selection.ligand.auth_seq_id:
        errors.append("selection.ligand.auth_seq_id is required.")
    if selection.ligand.available_instances and (
        selection.ligand.instance not in selection.ligand.available_instances
    ):
        errors.append("selection.ligand.instance must appear in available_instances.")

    processing = manifest.processing
    if processing.assembly_transformations_applied and not processing.assembly_id:
        errors.append(
            "processing.assembly_id is required when assembly transformations were applied."
        )
    if not isinstance(processing.assembly_transformations_applied, bool):
        errors.append("processing.assembly_transformations_applied must be a boolean.")
    if not isinstance(processing.original_coordinate_frame_preserved, bool):
        errors.append("processing.original_coordinate_frame_preserved must be a boolean.")
    if processing.translation_applied is not None and len(processing.translation_applied) != 3:
        errors.append("processing.translation_applied must contain three numbers.")
    if processing.rotation_applied is not None and len(processing.rotation_applied) != 9:
        errors.append("processing.rotation_applied must contain nine numbers.")
    if processing.original_coordinate_frame_preserved and (
        processing.translation_applied is not None
        or processing.rotation_applied is not None
    ):
        errors.append(
            "processing.original_coordinate_frame_preserved cannot be true when a "
            "translation or rotation was applied."
        )
    if source.deposited_asymmetric_unit and processing.assembly_transformations_applied:
        errors.append(
            "assembly transformations cannot be marked as applied when the source "
            "is the deposited asymmetric unit."
        )
    if not source.deposited_asymmetric_unit and not processing.assembly_id:
        errors.append(
            "processing.assembly_id is required when the source is not the "
            "deposited asymmetric unit."
        )
    if (
        not processing.original_coordinate_frame_preserved
        and not processing.coordinate_transforms
        and processing.translation_applied is None
        and processing.rotation_applied is None
    ):
        errors.append(
            "a changed coordinate frame must record coordinate_transforms, a "
            "translation, or a rotation."
        )

    if not manifest.decisions:
        errors.append("decisions must contain at least one selection or processing decision.")
    for index, decision in enumerate(manifest.decisions):
        if not decision.decision or not decision.actor or not decision.rationale:
            errors.append(
                f"decisions[{index}] requires decision, actor, and rationale."
            )

    chemical = manifest.chemical_identity
    if not chemical.reference:
        errors.append("chemical_identity.reference is required.")
    if not chemical.formula:
        errors.append("chemical_identity.formula is required.")

    analysis = manifest.analysis
    if analysis.coordinate_unit != "angstrom":
        errors.append("analysis.coordinate_unit must preserve deposited coordinates in angstrom.")
    if analysis.internal_distance_unit != "nm":
        errors.append("analysis.internal_distance_unit must be 'nm'.")

    if not manifest.provenance.records:
        errors.append("provenance.records must contain at least one labelled record.")
    for index, record in enumerate(manifest.provenance.records):
        if record.label not in PROVENANCE_LABELS:
            errors.append(
                f"provenance.records[{index}].label must be one of {PROVENANCE_LABELS}."
            )
        if not record.subject or not record.detail:
            errors.append(f"provenance.records[{index}] requires subject and detail.")
        if record.actor == "":
            errors.append(f"provenance.records[{index}].actor cannot be empty.")

    check_ids: set[str] = set()
    for index, check in enumerate(manifest.validation.checks):
        if check.status not in VALIDATION_STATUSES:
            errors.append(
                f"validation.checks[{index}].status must be one of {VALIDATION_STATUSES}."
            )
        if not check.check_id or not check.category or not check.detail:
            errors.append(
                f"validation.checks[{index}] requires check_id, category, and detail."
            )
        if not check.subject or not check.criterion:
            errors.append(
                f"validation.checks[{index}] requires subject and criterion."
            )
        if check.check_id in check_ids:
            errors.append(f"validation.checks[{index}].check_id is duplicated.")
        check_ids.add(check.check_id)

    return errors
