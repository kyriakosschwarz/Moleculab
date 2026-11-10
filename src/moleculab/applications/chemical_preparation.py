# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Separate chemical-state preparation records.

Preparation is a computational layer that references the experimental and
identity records; it does not modify them.  A record can hold multiple states
(e.g. an explicit reference state and a deposited-only state) and may leave the
selected state unset.

The module accepts prepared atoms, rings, and water roles from a chemistry
adapter.  It does not import RDKit or infer protonation, formal charges, or
features on its own.
"""

from __future__ import annotations

import json
import math
import re
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .structure_manifest import (
    PROVENANCE_LABELS,
    DecisionRecord,
    LigandSelection,
    ProvenanceEntry,
    ProvenanceRecord,
    ValidationCheck,
    ValidationRecord,
)

CHEMICAL_PREPARATION_SCHEMA_VERSION = "1.0"
WATER_ROLES = ("unassigned", "donor", "acceptor", "both", "unknown")
COORDINATE_SOURCES = ("deposited", "generated", "reference")
DISCREPANCY_STATUSES = ("open", "resolved", "not_applicable")
_SHA256_PATTERN = re.compile(r"[0-9a-f]{64}")


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


def _required_string(data: Mapping[str, Any], key: str) -> str:
    value = data.get(key)
    if value is None or str(value) == "":
        raise ValueError(f"{key} is required.")
    return str(value)


def _optional_string(data: Mapping[str, Any], key: str) -> str | None:
    value = data.get(key)
    return None if value is None else str(value)


def _required_int(data: Mapping[str, Any], key: str) -> int:
    value = data.get(key)
    if value is None:
        raise ValueError(f"{key} is required.")
    return int(value)


def _optional_int(data: Mapping[str, Any], key: str) -> int | None:
    value = data.get(key)
    return None if value is None else int(value)


def _required_float(data: Mapping[str, Any], key: str) -> float:
    value = data.get(key)
    if value is None:
        raise ValueError(f"{key} is required.")
    return float(value)


def _optional_bool(data: Mapping[str, Any], key: str) -> bool | None:
    value = data.get(key)
    return None if value is None else bool(value)


def _required_bool(data: Mapping[str, Any], key: str) -> bool:
    value = data.get(key)
    if not isinstance(value, bool):
        raise ValueError(f"{key} must be a boolean.")
    return value


def _vector3(data: Mapping[str, Any], key: str) -> tuple[float, float, float]:
    value = _as_list(data.get(key), key)
    if len(value) != 3:
        raise ValueError(f"{key} must contain three numbers.")
    return (float(value[0]), float(value[1]), float(value[2]))


def _vector_tuple(data: Mapping[str, Any], key: str) -> tuple[tuple[float, ...], ...]:
    return tuple(
        tuple(float(item) for item in _as_list(vector, key))
        for vector in _as_list(data.get(key), key)
    )


def _string_tuple(data: Mapping[str, Any], key: str) -> tuple[str, ...]:
    return tuple(str(item) for item in _as_list(data.get(key), key))


@dataclass(frozen=True)
class PreparedAtom:
    atom_id: str
    element: str
    formal_charge: int | None
    is_donor: bool | None
    is_acceptor: bool | None
    is_hydrophobic: bool | None
    aromatic_ring_ids: tuple[str, ...]
    hydrogen_count: int
    hydrogen_positions_nm: tuple[tuple[float, float, float], ...]
    coordinate_source: str
    hydrogen_source: str | None
    feature_source: str | None

    def to_dict(self) -> dict[str, Any]:
        return {
            "atom_id": self.atom_id,
            "element": self.element,
            "formal_charge": self.formal_charge,
            "is_donor": self.is_donor,
            "is_acceptor": self.is_acceptor,
            "is_hydrophobic": self.is_hydrophobic,
            "aromatic_ring_ids": list(self.aromatic_ring_ids),
            "hydrogen_count": self.hydrogen_count,
            "hydrogen_positions_nm": [
                list(position) for position in self.hydrogen_positions_nm
            ],
            "coordinate_source": self.coordinate_source,
            "hydrogen_source": self.hydrogen_source,
            "feature_source": self.feature_source,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> PreparedAtom:
        data = _as_mapping(data, "prepared atom")
        return cls(
            atom_id=_required_string(data, "atom_id"),
            element=_required_string(data, "element"),
            formal_charge=_optional_int(data, "formal_charge"),
            is_donor=_optional_bool(data, "is_donor"),
            is_acceptor=_optional_bool(data, "is_acceptor"),
            is_hydrophobic=_optional_bool(data, "is_hydrophobic"),
            aromatic_ring_ids=_string_tuple(data, "aromatic_ring_ids"),
            hydrogen_count=_required_int(data, "hydrogen_count"),
            hydrogen_positions_nm=_vector_tuple(data, "hydrogen_positions_nm"),
            coordinate_source=_required_string(data, "coordinate_source"),
            hydrogen_source=_optional_string(data, "hydrogen_source"),
            feature_source=_optional_string(data, "feature_source"),
        )


@dataclass(frozen=True)
class PreparedRing:
    ring_id: str
    atom_ids: tuple[str, ...]
    centroid_nm: tuple[float, float, float]
    normal: tuple[float, float, float]
    source: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "ring_id": self.ring_id,
            "atom_ids": list(self.atom_ids),
            "centroid_nm": list(self.centroid_nm),
            "normal": list(self.normal),
            "source": self.source,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> PreparedRing:
        data = _as_mapping(data, "prepared ring")
        return cls(
            ring_id=_required_string(data, "ring_id"),
            atom_ids=_string_tuple(data, "atom_ids"),
            centroid_nm=_vector3(data, "centroid_nm"),
            normal=_vector3(data, "normal"),
            source=_required_string(data, "source"),
        )


@dataclass(frozen=True)
class PreparedWaterRole:
    water_id: str
    role: str
    source: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "water_id": self.water_id,
            "role": self.role,
            "source": self.source,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> PreparedWaterRole:
        data = _as_mapping(data, "prepared water role")
        return cls(
            water_id=_required_string(data, "water_id"),
            role=_required_string(data, "role"),
            source=_required_string(data, "source"),
        )


@dataclass(frozen=True)
class PreparedState:
    state_id: str
    label: str
    protonation_state: str
    tautomer: str
    formal_charge: int | None
    hydrogen_method: str | None
    feature_source: str | None
    assumptions: tuple[str, ...]
    atoms: tuple[PreparedAtom, ...]
    rings: tuple[PreparedRing, ...]
    water_roles: tuple[PreparedWaterRole, ...]

    def to_dict(self) -> dict[str, Any]:
        return {
            "state_id": self.state_id,
            "label": self.label,
            "protonation_state": self.protonation_state,
            "tautomer": self.tautomer,
            "formal_charge": self.formal_charge,
            "hydrogen_method": self.hydrogen_method,
            "feature_source": self.feature_source,
            "assumptions": list(self.assumptions),
            "atoms": [atom.to_dict() for atom in self.atoms],
            "rings": [ring.to_dict() for ring in self.rings],
            "water_roles": [role.to_dict() for role in self.water_roles],
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> PreparedState:
        data = _as_mapping(data, "prepared state")
        return cls(
            state_id=_required_string(data, "state_id"),
            label=_required_string(data, "label"),
            protonation_state=_required_string(data, "protonation_state"),
            tautomer=_required_string(data, "tautomer"),
            formal_charge=_optional_int(data, "formal_charge"),
            hydrogen_method=_optional_string(data, "hydrogen_method"),
            feature_source=_optional_string(data, "feature_source"),
            assumptions=_string_tuple(data, "assumptions"),
            atoms=tuple(
                PreparedAtom.from_dict(_as_mapping(item, "prepared atom"))
                for item in _as_list(data.get("atoms"), "atoms")
            ),
            rings=tuple(
                PreparedRing.from_dict(_as_mapping(item, "prepared ring"))
                for item in _as_list(data.get("rings"), "rings")
            ),
            water_roles=tuple(
                PreparedWaterRole.from_dict(_as_mapping(item, "prepared water role"))
                for item in _as_list(data.get("water_roles"), "water_roles")
            ),
        )


@dataclass(frozen=True)
class PreparationDiscrepancy:
    category: str
    status: str
    detail: str
    observed: str | None = None
    expected: str | None = None
    resolution: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "category": self.category,
            "status": self.status,
            "detail": self.detail,
            "observed": self.observed,
            "expected": self.expected,
            "resolution": self.resolution,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> PreparationDiscrepancy:
        data = _as_mapping(data, "preparation discrepancy")
        return cls(
            category=_required_string(data, "category"),
            status=_required_string(data, "status"),
            detail=_required_string(data, "detail"),
            observed=_optional_string(data, "observed"),
            expected=_optional_string(data, "expected"),
            resolution=_optional_string(data, "resolution"),
        )


@dataclass(frozen=True)
class ChemicalPreparationRecord:
    source_manifest_path: str
    source_manifest_sha256: str
    ligand_identity_record_path: str
    ligand_identity_record_sha256: str
    ligand_coordinate_record_path: str
    ligand_coordinate_record_sha256: str
    binding_site_record_path: str
    binding_site_record_sha256: str
    accession: str
    selection: LigandSelection
    toolkit: str
    toolkit_version: str
    parameters: tuple[str, ...]
    states: tuple[PreparedState, ...]
    selected_state_id: str | None
    decisions: tuple[DecisionRecord, ...]
    provenance: ProvenanceRecord
    validation: ValidationRecord
    discrepancies: tuple[PreparationDiscrepancy, ...]
    schema_version: str = CHEMICAL_PREPARATION_SCHEMA_VERSION

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "source_manifest_path": self.source_manifest_path,
            "source_manifest_sha256": self.source_manifest_sha256,
            "ligand_identity_record_path": self.ligand_identity_record_path,
            "ligand_identity_record_sha256": self.ligand_identity_record_sha256,
            "ligand_coordinate_record_path": self.ligand_coordinate_record_path,
            "ligand_coordinate_record_sha256": self.ligand_coordinate_record_sha256,
            "binding_site_record_path": self.binding_site_record_path,
            "binding_site_record_sha256": self.binding_site_record_sha256,
            "accession": self.accession,
            "selection": self.selection.to_dict(),
            "toolkit": self.toolkit,
            "toolkit_version": self.toolkit_version,
            "parameters": list(self.parameters),
            "states": [state.to_dict() for state in self.states],
            "selected_state_id": self.selected_state_id,
            "decisions": [decision.to_dict() for decision in self.decisions],
            "provenance": self.provenance.to_dict(),
            "validation": self.validation.to_dict(),
            "discrepancies": [item.to_dict() for item in self.discrepancies],
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> ChemicalPreparationRecord:
        data = _as_mapping(data, "chemical preparation record")
        return cls(
            schema_version=_required_string(data, "schema_version"),
            source_manifest_path=_required_string(data, "source_manifest_path"),
            source_manifest_sha256=_required_string(data, "source_manifest_sha256"),
            ligand_identity_record_path=_required_string(
                data, "ligand_identity_record_path"
            ),
            ligand_identity_record_sha256=_required_string(
                data, "ligand_identity_record_sha256"
            ),
            ligand_coordinate_record_path=_required_string(
                data, "ligand_coordinate_record_path"
            ),
            ligand_coordinate_record_sha256=_required_string(
                data, "ligand_coordinate_record_sha256"
            ),
            binding_site_record_path=_required_string(data, "binding_site_record_path"),
            binding_site_record_sha256=_required_string(
                data, "binding_site_record_sha256"
            ),
            accession=_required_string(data, "accession"),
            selection=LigandSelection.from_dict(
                _as_mapping(data.get("selection"), "selection")
            ),
            toolkit=_required_string(data, "toolkit"),
            toolkit_version=_required_string(data, "toolkit_version"),
            parameters=_string_tuple(data, "parameters"),
            states=tuple(
                PreparedState.from_dict(_as_mapping(item, "prepared state"))
                for item in _as_list(data.get("states"), "states")
            ),
            selected_state_id=_optional_string(data, "selected_state_id"),
            decisions=tuple(
                DecisionRecord.from_dict(_as_mapping(item, "decision record"))
                for item in _as_list(data.get("decisions"), "decisions")
            ),
            provenance=ProvenanceRecord.from_dict(
                _as_mapping(data.get("provenance"), "provenance")
            ),
            validation=ValidationRecord.from_dict(
                _as_mapping(data.get("validation"), "validation")
            ),
            discrepancies=tuple(
                PreparationDiscrepancy.from_dict(
                    _as_mapping(item, "preparation discrepancy")
                )
                for item in _as_list(data.get("discrepancies"), "discrepancies")
            ),
        )


def _check(
    category: str,
    check_id: str,
    ok: bool,
    *,
    subject: str,
    criterion: str,
    observed: str,
    expected: str,
    detail: str,
    not_run: bool = False,
    not_applicable: bool = False,
) -> ValidationCheck:
    if not_run:
        status = "not_run"
    elif not_applicable:
        status = "not_applicable"
    else:
        status = "pass" if ok else "fail"
    return ValidationCheck(
        check_id=check_id,
        category=category,
        status=status,
        subject=subject,
        criterion=criterion,
        observed=observed,
        expected=expected,
        detail=detail,
    )


def build_chemical_preparation_record(
    *,
    source_manifest_path: str,
    source_manifest_sha256: str,
    ligand_identity_record_path: str,
    ligand_identity_record_sha256: str,
    ligand_coordinate_record_path: str,
    ligand_coordinate_record_sha256: str,
    binding_site_record_path: str,
    binding_site_record_sha256: str,
    accession: str,
    selection: LigandSelection,
    toolkit: str,
    toolkit_version: str,
    parameters: tuple[str, ...],
    states: tuple[PreparedState, ...],
    selected_state_id: str | None = None,
) -> ChemicalPreparationRecord:
    """Assemble and validate a preparation record without selecting a state."""

    if not states:
        raise ValueError("states must contain at least one prepared state.")
    state_ids = [state.state_id for state in states]
    if len(set(state_ids)) != len(state_ids):
        raise ValueError("state_id values must be unique.")
    if selected_state_id is not None and selected_state_id not in state_ids:
        raise ValueError("selected_state_id must name one of the states.")

    checks = (
        _check(
            "provenance",
            "preparation.source_records",
            bool(
                source_manifest_path
                and ligand_identity_record_path
                and ligand_coordinate_record_path
                and binding_site_record_path
            ),
            subject=accession,
            criterion="all experimental and identity inputs are recorded",
            observed="four input record paths present",
            expected="manifest, identity, coordinate, and binding-site records",
            detail="preparation references the prior evidence chain",
        ),
        _check(
            "provenance",
            "preparation.toolkit_recorded",
            bool(toolkit and toolkit_version and parameters),
            subject=selection.component_id,
            criterion="toolkit, version, and parameters are explicit",
            observed=f"{toolkit} {toolkit_version}",
            expected="documented preparation method",
            detail="; ".join(parameters),
        ),
        _check(
            "states",
            "preparation.multiple_states_preserved",
            len(states) > 1,
            subject=selection.component_id,
            criterion="plausible states can coexist without silent selection",
            observed=f"{len(states)} states; selected={selected_state_id!r}",
            expected="one or more states; explicit selection or none",
            detail="no state is selected unless the caller says so",
        ),
        _check(
            "states",
            "preparation.unknown_features_preserved",
            any(
                atom.formal_charge is None
                or atom.is_donor is None
                or atom.is_acceptor is None
                or atom.is_hydrophobic is None
                for state in states
                for atom in state.atoms
            ),
            subject=selection.component_id,
            criterion="unsupported chemistry remains unknown",
            observed="unknown flags are retained where supplied",
            expected="no element-based inference",
            detail="unknown remains distinct from false",
        ),
        _check(
            "interpretation",
            "preparation.no_affinity",
            True,
            subject=selection.component_id,
            criterion="preparation is not converted into affinity",
            observed="no affinity or free-energy estimate is produced",
            expected="chemical-state preparation only",
            detail="prepared features are inputs to later analysis",
        ),
    )

    discrepancies: list[PreparationDiscrepancy] = []
    if selected_state_id is None:
        discrepancies.append(
            PreparationDiscrepancy(
                "state_selection",
                "not_applicable",
                "No prepared state was selected.",
                observed="selected_state_id is null",
                expected="explicit selection or coexistence of states",
                resolution="The record preserves all states without choosing one.",
            )
        )
    if len(states) == 1:
        discrepancies.append(
            PreparationDiscrepancy(
                "alternatives",
                "not_applicable",
                "Only one prepared state is present.",
                observed="1 state",
                expected="multiple states when alternatives are plausible",
                resolution="No alternative state was enumerated.",
            )
        )

    decisions = (
        DecisionRecord(
            decision="created a separate preparation record",
            actor="chemical-state preparation",
            rationale=(
                "deposited coordinates, identity, and binding-site records were not "
                "modified"
            ),
        ),
        DecisionRecord(
            decision="preserved multiple states without silent selection",
            actor="chemical-state preparation",
            rationale=(
                "preparation assumptions and generated chemistry remain separate from "
                "experimental observations"
            ),
        ),
        DecisionRecord(
            decision="recorded toolkit, version, parameters, and assumptions",
            actor="chemical-state preparation",
            rationale="prepared features must be reproducible and auditable",
        ),
        DecisionRecord(
            decision="did not convert preparation into affinity",
            actor="chemical-state preparation",
            rationale="prepared chemistry is not a binding free energy",
        ),
    )

    provenance = ProvenanceRecord(
        records=(
            ProvenanceEntry(
                "experimental",
                "input structural records",
                "Preparation references the manifest, identity, coordinate, and binding-site records.",
                actor="RCSB PDB deposition",
            ),
            ProvenanceEntry(
                "calculated",
                "prepared chemical states",
                (
                    f"{len(states)} state(s) from {toolkit} {toolkit_version}; "
                    f"parameters: {'; '.join(parameters)}."
                ),
                actor="chemical-state preparation",
            ),
            ProvenanceEntry(
                "predicted",
                "none",
                "No affinity or free-energy prediction is produced.",
                actor="chemical-state preparation",
            ),
            ProvenanceEntry(
                "hypothesised",
                "none",
                "No mechanism is inferred from a prepared state.",
                actor="chemical-state preparation",
            ),
        )
    )

    return ChemicalPreparationRecord(
        source_manifest_path=source_manifest_path,
        source_manifest_sha256=source_manifest_sha256,
        ligand_identity_record_path=ligand_identity_record_path,
        ligand_identity_record_sha256=ligand_identity_record_sha256,
        ligand_coordinate_record_path=ligand_coordinate_record_path,
        ligand_coordinate_record_sha256=ligand_coordinate_record_sha256,
        binding_site_record_path=binding_site_record_path,
        binding_site_record_sha256=binding_site_record_sha256,
        accession=accession,
        selection=selection,
        toolkit=toolkit,
        toolkit_version=toolkit_version,
        parameters=parameters,
        states=states,
        selected_state_id=selected_state_id,
        decisions=decisions,
        provenance=provenance,
        validation=ValidationRecord(checks=checks),
        discrepancies=tuple(discrepancies),
    )


def write_chemical_preparation(path: str | Path, record: ChemicalPreparationRecord) -> Path:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(record.to_dict(), indent=2, sort_keys=False) + "\n",
        encoding="utf-8",
    )
    return destination


def load_chemical_preparation(path: str | Path) -> ChemicalPreparationRecord:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    return ChemicalPreparationRecord.from_dict(data)


def validate_chemical_preparation(record: ChemicalPreparationRecord) -> list[str]:
    """Return schema and internal-consistency errors for *record*."""

    errors: list[str] = []
    if record.schema_version != CHEMICAL_PREPARATION_SCHEMA_VERSION:
        errors.append(
            f"schema_version must be {CHEMICAL_PREPARATION_SCHEMA_VERSION!r}, got "
            f"{record.schema_version!r}."
        )
    for name in (
        "source_manifest_path",
        "ligand_identity_record_path",
        "ligand_coordinate_record_path",
        "binding_site_record_path",
    ):
        if not getattr(record, name):
            errors.append(f"{name} is required.")
    for name in (
        "source_manifest_sha256",
        "ligand_identity_record_sha256",
        "ligand_coordinate_record_sha256",
        "binding_site_record_sha256",
    ):
        if not _SHA256_PATTERN.fullmatch(getattr(record, name)):
            errors.append(f"{name} must be a lowercase SHA-256 digest.")
    if not record.accession:
        errors.append("accession is required.")
    if not record.toolkit or not record.toolkit_version:
        errors.append("toolkit and toolkit_version are required.")
    if not record.parameters:
        errors.append("parameters must not be empty.")
    if not record.states:
        errors.append("states must not be empty.")

    state_ids = [state.state_id for state in record.states]
    if len(set(state_ids)) != len(state_ids):
        errors.append("state_id values must be unique.")
    if record.selected_state_id is not None and record.selected_state_id not in state_ids:
        errors.append("selected_state_id must name one of the states.")

    for state_index, state in enumerate(record.states):
        if not state.label or not state.protonation_state or not state.tautomer:
            errors.append(f"states[{state_index}] requires label, protonation, and tautomer.")
        if not state.assumptions:
            errors.append(f"states[{state_index}].assumptions must not be empty.")
        if not state.atoms:
            errors.append(f"states[{state_index}].atoms must not be empty.")
        atom_ids = [atom.atom_id for atom in state.atoms]
        if len(set(atom_ids)) != len(atom_ids):
            errors.append(f"states[{state_index}] atom_ids must be unique.")
        atom_id_set = set(atom_ids)
        for atom_index, atom in enumerate(state.atoms):
            if atom.coordinate_source not in COORDINATE_SOURCES:
                errors.append(
                    f"states[{state_index}].atoms[{atom_index}].coordinate_source is invalid."
                )
            if atom.hydrogen_count < 0:
                errors.append(
                    f"states[{state_index}].atoms[{atom_index}].hydrogen_count cannot be negative."
                )
            if atom.hydrogen_count != len(atom.hydrogen_positions_nm):
                errors.append(
                    f"states[{state_index}].atoms[{atom_index}] hydrogen count does not match positions."
                )
            if atom.hydrogen_count and not atom.hydrogen_source:
                errors.append(
                    f"states[{state_index}].atoms[{atom_index}] generated hydrogens require a source."
                )
            if not atom.hydrogen_count and atom.hydrogen_source:
                errors.append(
                    f"states[{state_index}].atoms[{atom_index}] has a hydrogen source but no hydrogens."
                )
            for hydrogen_index, position in enumerate(atom.hydrogen_positions_nm):
                if len(position) != 3 or not all(math.isfinite(value) for value in position):
                    errors.append(
                        f"states[{state_index}].atoms[{atom_index}].hydrogen_positions_nm"
                        f"[{hydrogen_index}] must be a finite three-vector."
                    )
        for ring_index, ring in enumerate(state.rings):
            if not set(ring.atom_ids).issubset(atom_id_set):
                errors.append(
                    f"states[{state_index}].rings[{ring_index}] references an unknown atom."
                )
            normal = ring.normal
            if not all(math.isfinite(value) for value in normal) or not any(normal):
                errors.append(
                    f"states[{state_index}].rings[{ring_index}].normal must be finite and nonzero."
                )
        for water_index, water in enumerate(state.water_roles):
            if water.role not in WATER_ROLES:
                errors.append(
                    f"states[{state_index}].water_roles[{water_index}].role is invalid."
                )

    if not record.decisions:
        errors.append("decisions must contain at least one record.")
    for index, decision in enumerate(record.decisions):
        if not decision.decision or not decision.actor or not decision.rationale:
            errors.append(f"decisions[{index}] requires decision, actor, and rationale.")

    if not record.provenance.records:
        errors.append("provenance.records must not be empty.")
    for index, entry in enumerate(record.provenance.records):
        if entry.label not in PROVENANCE_LABELS:
            errors.append(f"provenance.records[{index}] has an invalid label.")

    check_ids: set[str] = set()
    for index, check in enumerate(record.validation.checks):
        if check.status not in ("pass", "fail", "not_run", "not_applicable"):
            errors.append(f"validation.checks[{index}] has an invalid status.")
        if not check.check_id or not check.category or not check.detail:
            errors.append(f"validation.checks[{index}] requires check_id, category, and detail.")
        if check.check_id in check_ids:
            errors.append(f"validation.checks[{index}].check_id is duplicated.")
        check_ids.add(check.check_id)

    for index, item in enumerate(record.discrepancies):
        if item.status not in DISCREPANCY_STATUSES:
            errors.append(f"discrepancies[{index}] has an invalid status.")
        if not item.category or not item.detail:
            errors.append(f"discrepancies[{index}] requires category and detail.")
        if item.status == "resolved" and not item.resolution:
            errors.append(f"discrepancies[{index}] is resolved but has no resolution.")

    return errors
