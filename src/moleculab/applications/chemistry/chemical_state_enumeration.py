# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Enumeration records for protonation and tautomer states.

This layer stores explicitly enumerated chemical states.  It does not predict
pKa values, estimate populations, or choose a dominant state.  Protonation and
tautomer alternatives are distinguished, each state carries its own molecular
graph, proton count, formal charge, features, and coordinate provenance, and
the original experimental and preparation records are referenced, not modified.
"""

from __future__ import annotations

import json
import re
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .chemical_preparation import PreparedAtom, PreparedRing, PreparedWaterRole
from ..structure.structure_manifest import (
    PROVENANCE_LABELS,
    DecisionRecord,
    LigandSelection,
    ProvenanceEntry,
    ProvenanceRecord,
    ValidationCheck,
    ValidationRecord,
)

CHEMICAL_STATE_SCHEMA_VERSION = "1.0"
STATE_KINDS = ("reference", "protonation", "tautomer")
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


def _optional_float(data: Mapping[str, Any], key: str) -> float | None:
    value = data.get(key)
    return None if value is None else float(value)


def _string_tuple(data: Mapping[str, Any], key: str) -> tuple[str, ...]:
    return tuple(str(item) for item in _as_list(data.get(key), key))


@dataclass(frozen=True)
class IonizableSite:
    site_id: str
    atom_id: str
    site_type: str
    pka: float | None
    pka_source: str | None
    pka_uncertainty: float | None
    note: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "site_id": self.site_id,
            "atom_id": self.atom_id,
            "site_type": self.site_type,
            "pka": self.pka,
            "pka_source": self.pka_source,
            "pka_uncertainty": self.pka_uncertainty,
            "note": self.note,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> IonizableSite:
        data = _as_mapping(data, "ionizable site")
        return cls(
            site_id=_required_string(data, "site_id"),
            atom_id=_required_string(data, "atom_id"),
            site_type=_required_string(data, "site_type"),
            pka=_optional_float(data, "pka"),
            pka_source=_optional_string(data, "pka_source"),
            pka_uncertainty=_optional_float(data, "pka_uncertainty"),
            note=_optional_string(data, "note"),
        )


@dataclass(frozen=True)
class ChemicalState:
    state_id: str
    kind: str
    label: str
    parent_state_id: str | None
    molecular_graph: str
    formula: str
    proton_count: int
    formal_charge: int
    tautomer_of: str | None
    ph_assumption: float | None
    pka_sites: tuple[IonizableSite, ...]
    enumeration_method: str
    enumeration_parameters: tuple[str, ...]
    assumptions: tuple[str, ...]
    atoms: tuple[PreparedAtom, ...]
    reference_atom_mapping: tuple[tuple[str, int], ...]
    rings: tuple[PreparedRing, ...]
    water_roles: tuple[PreparedWaterRole, ...]

    def to_dict(self) -> dict[str, Any]:
        return {
            "state_id": self.state_id,
            "kind": self.kind,
            "label": self.label,
            "parent_state_id": self.parent_state_id,
            "molecular_graph": self.molecular_graph,
            "formula": self.formula,
            "proton_count": self.proton_count,
            "formal_charge": self.formal_charge,
            "tautomer_of": self.tautomer_of,
            "ph_assumption": self.ph_assumption,
            "pka_sites": [site.to_dict() for site in self.pka_sites],
            "enumeration_method": self.enumeration_method,
            "enumeration_parameters": list(self.enumeration_parameters),
            "assumptions": list(self.assumptions),
            "atoms": [atom.to_dict() for atom in self.atoms],
            "reference_atom_mapping": [
                {"atom_id": atom_id, "reference_index": reference_index}
                for atom_id, reference_index in self.reference_atom_mapping
            ],
            "rings": [ring.to_dict() for ring in self.rings],
            "water_roles": [role.to_dict() for role in self.water_roles],
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> ChemicalState:
        data = _as_mapping(data, "chemical state")
        return cls(
            state_id=_required_string(data, "state_id"),
            kind=_required_string(data, "kind"),
            label=_required_string(data, "label"),
            parent_state_id=_optional_string(data, "parent_state_id"),
            molecular_graph=_required_string(data, "molecular_graph"),
            formula=_required_string(data, "formula"),
            proton_count=_required_int(data, "proton_count"),
            formal_charge=_required_int(data, "formal_charge"),
            tautomer_of=_optional_string(data, "tautomer_of"),
            ph_assumption=_optional_float(data, "ph_assumption"),
            pka_sites=tuple(
                IonizableSite.from_dict(_as_mapping(item, "ionizable site"))
                for item in _as_list(data.get("pka_sites"), "pka_sites")
            ),
            enumeration_method=_required_string(data, "enumeration_method"),
            enumeration_parameters=_string_tuple(data, "enumeration_parameters"),
            assumptions=_string_tuple(data, "assumptions"),
            atoms=tuple(
                PreparedAtom.from_dict(_as_mapping(item, "prepared atom"))
                for item in _as_list(data.get("atoms"), "atoms")
            ),
            reference_atom_mapping=tuple(
                (
                    _required_string(_as_mapping(item, "reference atom mapping"), "atom_id"),
                    _required_int(
                        _as_mapping(item, "reference atom mapping"), "reference_index"
                    ),
                )
                for item in _as_list(
                    data.get("reference_atom_mapping"), "reference_atom_mapping"
                )
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
class StateEnumerationDiscrepancy:
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
    def from_dict(cls, data: Mapping[str, Any]) -> StateEnumerationDiscrepancy:
        data = _as_mapping(data, "state enumeration discrepancy")
        return cls(
            category=_required_string(data, "category"),
            status=_required_string(data, "status"),
            detail=_required_string(data, "detail"),
            observed=_optional_string(data, "observed"),
            expected=_optional_string(data, "expected"),
            resolution=_optional_string(data, "resolution"),
        )


@dataclass(frozen=True)
class ChemicalStateEnumerationRecord:
    source_manifest_path: str
    source_manifest_sha256: str
    ligand_identity_record_path: str
    ligand_identity_record_sha256: str
    ligand_coordinate_record_path: str
    ligand_coordinate_record_sha256: str
    binding_site_record_path: str
    binding_site_record_sha256: str
    chemical_preparation_record_path: str
    chemical_preparation_record_sha256: str
    accession: str
    selection: LigandSelection
    reference_state_id: str
    states: tuple[ChemicalState, ...]
    selected_state_id: str | None
    ph_assumption: float | None
    population_model: str | None
    toolkit: str
    toolkit_version: str
    parameters: tuple[str, ...]
    decisions: tuple[DecisionRecord, ...]
    provenance: ProvenanceRecord
    validation: ValidationRecord
    discrepancies: tuple[StateEnumerationDiscrepancy, ...]
    schema_version: str = CHEMICAL_STATE_SCHEMA_VERSION

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
            "chemical_preparation_record_path": self.chemical_preparation_record_path,
            "chemical_preparation_record_sha256": self.chemical_preparation_record_sha256,
            "accession": self.accession,
            "selection": self.selection.to_dict(),
            "reference_state_id": self.reference_state_id,
            "states": [state.to_dict() for state in self.states],
            "selected_state_id": self.selected_state_id,
            "ph_assumption": self.ph_assumption,
            "population_model": self.population_model,
            "toolkit": self.toolkit,
            "toolkit_version": self.toolkit_version,
            "parameters": list(self.parameters),
            "decisions": [decision.to_dict() for decision in self.decisions],
            "provenance": self.provenance.to_dict(),
            "validation": self.validation.to_dict(),
            "discrepancies": [item.to_dict() for item in self.discrepancies],
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> ChemicalStateEnumerationRecord:
        data = _as_mapping(data, "chemical state enumeration record")
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
            chemical_preparation_record_path=_required_string(
                data, "chemical_preparation_record_path"
            ),
            chemical_preparation_record_sha256=_required_string(
                data, "chemical_preparation_record_sha256"
            ),
            accession=_required_string(data, "accession"),
            selection=LigandSelection.from_dict(
                _as_mapping(data.get("selection"), "selection")
            ),
            reference_state_id=_required_string(data, "reference_state_id"),
            states=tuple(
                ChemicalState.from_dict(_as_mapping(item, "chemical state"))
                for item in _as_list(data.get("states"), "states")
            ),
            selected_state_id=_optional_string(data, "selected_state_id"),
            ph_assumption=_optional_float(data, "ph_assumption"),
            population_model=_optional_string(data, "population_model"),
            toolkit=_required_string(data, "toolkit"),
            toolkit_version=_required_string(data, "toolkit_version"),
            parameters=_string_tuple(data, "parameters"),
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
                StateEnumerationDiscrepancy.from_dict(
                    _as_mapping(item, "state enumeration discrepancy")
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
) -> ValidationCheck:
    if not_run:
        status = "not_run"
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


def build_chemical_state_enumeration_record(
    *,
    source_manifest_path: str,
    source_manifest_sha256: str,
    ligand_identity_record_path: str,
    ligand_identity_record_sha256: str,
    ligand_coordinate_record_path: str,
    ligand_coordinate_record_sha256: str,
    binding_site_record_path: str,
    binding_site_record_sha256: str,
    chemical_preparation_record_path: str,
    chemical_preparation_record_sha256: str,
    accession: str,
    selection: LigandSelection,
    reference_state_id: str,
    states: tuple[ChemicalState, ...],
    toolkit: str,
    toolkit_version: str,
    parameters: tuple[str, ...],
    selected_state_id: str | None = None,
    ph_assumption: float | None = None,
    population_model: str | None = None,
) -> ChemicalStateEnumerationRecord:
    """Assemble a state-enumeration record without selecting a state."""

    state_ids = [state.state_id for state in states]
    if not states:
        raise ValueError("states must contain at least one entry.")
    if len(set(state_ids)) != len(state_ids):
        raise ValueError("state_id values must be unique.")
    if reference_state_id not in state_ids:
        raise ValueError("reference_state_id must name one of the states.")
    if selected_state_id is not None and selected_state_id not in state_ids:
        raise ValueError("selected_state_id must name one of the states.")
    if population_model is not None:
        raise ValueError(
            "population estimation is outside this batch; supply no population model"
        )

    checks = (
        _check(
            "provenance",
            "states.source_records",
            bool(
                source_manifest_path
                and ligand_identity_record_path
                and ligand_coordinate_record_path
                and binding_site_record_path
                and chemical_preparation_record_path
            ),
            subject=accession,
            criterion="all experimental and preparation inputs are recorded",
            observed="five input record paths present",
            expected="manifest, identity, coordinates, binding site, preparation",
            detail="the enumeration references, but does not modify, prior records",
        ),
        _check(
            "enumeration",
            "states.reference_present",
            reference_state_id in state_ids,
            subject=selection.component_id,
            criterion="the reference state is present",
            observed=reference_state_id,
            expected="one reference state",
            detail="alternatives are compared against the reference state",
        ),
        _check(
            "enumeration",
            "states.alternatives_present",
            any(state.kind != "reference" for state in states),
            subject=selection.component_id,
            criterion="at least one protonation or tautomer alternative is present",
            observed=f"{len(states)} states",
            expected="reference plus alternatives",
            detail="no population or dominance claim is made",
        ),
        _check(
            "enumeration",
            "states.no_default_selection",
            selected_state_id is None,
            subject=selection.component_id,
            criterion="no state is silently selected",
            observed=f"selected_state_id={selected_state_id!r}",
            expected="null unless an explicit downstream choice is made",
            detail="alternatives coexist in the record",
        ),
        _check(
            "interpretation",
            "states.population_model_absent",
            population_model is None,
            subject=selection.component_id,
            criterion="population estimates require a defensible thermodynamic model",
            observed="no population model supplied",
            expected="enumeration only",
            detail="pH-dependent populations are not estimated in this batch",
        ),
        _check(
            "interpretation",
            "states.no_affinity",
            True,
            subject=selection.component_id,
            criterion="state enumeration is not converted into affinity",
            observed="no affinity or free-energy estimate is produced",
            expected="chemical-state enumeration only",
            detail="prepared states are inputs to later analyses",
        ),
    )

    discrepancies: list[StateEnumerationDiscrepancy] = []
    if selected_state_id is None:
        discrepancies.append(
            StateEnumerationDiscrepancy(
                "state_selection",
                "not_applicable",
                "No prepared state was selected.",
                observed="selected_state_id is null",
                expected="explicit selection or coexistence of alternatives",
                resolution="All enumerated states are preserved.",
            )
        )
    if population_model is None:
        discrepancies.append(
            StateEnumerationDiscrepancy(
                "population",
                "not_applicable",
                "No pH-dependent population model was supplied.",
                observed="no population model",
                expected="explicit thermodynamic model for populations",
                resolution="The record does not rank or weight the states.",
            )
        )

    decisions = (
        DecisionRecord(
            decision="enumerated alternative states without selecting a default",
            actor="chemical-state enumeration",
            rationale="protonation and tautomer alternatives are hypotheses, not observations",
        ),
        DecisionRecord(
            decision="kept protonation and tautomer enumeration separate",
            actor="chemical-state enumeration",
            rationale="tautomers share composition; protonation changes proton count and charge",
        ),
        DecisionRecord(
            decision="did not estimate populations",
            actor="chemical-state enumeration",
            rationale="no pH-dependent thermodynamic model or pKa source was supplied",
        ),
        DecisionRecord(
            decision="did not modify experimental or Batch 7 preparation records",
            actor="chemical-state enumeration",
            rationale="the enumeration is a separate downstream record",
        ),
    )

    provenance = ProvenanceRecord(
        records=(
            ProvenanceEntry(
                "experimental",
                "input structural records",
                "Enumeration references the experimental chain and the Batch 7 preparation record.",
                actor="RCSB PDB deposition",
            ),
            ProvenanceEntry(
                "calculated",
                "enumerated chemical states",
                (
                    f"{len(states)} states from {toolkit} {toolkit_version}; "
                    f"parameters: {'; '.join(parameters)}."
                ),
                actor="chemical-state enumeration",
            ),
            ProvenanceEntry(
                "predicted",
                "none",
                "No affinity, population, or free-energy estimate is produced.",
                actor="chemical-state enumeration",
            ),
            ProvenanceEntry(
                "hypothesised",
                "protonation and tautomer alternatives",
                "Alternative states are hypotheses for later evaluation, not dominant species.",
                actor="chemical-state enumeration",
            ),
        )
    )

    return ChemicalStateEnumerationRecord(
        source_manifest_path=source_manifest_path,
        source_manifest_sha256=source_manifest_sha256,
        ligand_identity_record_path=ligand_identity_record_path,
        ligand_identity_record_sha256=ligand_identity_record_sha256,
        ligand_coordinate_record_path=ligand_coordinate_record_path,
        ligand_coordinate_record_sha256=ligand_coordinate_record_sha256,
        binding_site_record_path=binding_site_record_path,
        binding_site_record_sha256=binding_site_record_sha256,
        chemical_preparation_record_path=chemical_preparation_record_path,
        chemical_preparation_record_sha256=chemical_preparation_record_sha256,
        accession=accession,
        selection=selection,
        reference_state_id=reference_state_id,
        states=states,
        selected_state_id=selected_state_id,
        ph_assumption=ph_assumption,
        population_model=population_model,
        toolkit=toolkit,
        toolkit_version=toolkit_version,
        parameters=parameters,
        decisions=decisions,
        provenance=provenance,
        validation=ValidationRecord(checks=checks),
        discrepancies=tuple(discrepancies),
    )


def write_chemical_state_enumeration(path: str | Path, record: ChemicalStateEnumerationRecord) -> Path:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(record.to_dict(), indent=2, sort_keys=False) + "\n",
        encoding="utf-8",
    )
    return destination


def load_chemical_state_enumeration(path: str | Path) -> ChemicalStateEnumerationRecord:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    return ChemicalStateEnumerationRecord.from_dict(data)


def validate_chemical_state_enumeration(record: ChemicalStateEnumerationRecord) -> list[str]:
    """Return schema and internal-consistency errors for *record*."""

    errors: list[str] = []
    if record.schema_version != CHEMICAL_STATE_SCHEMA_VERSION:
        errors.append(
            f"schema_version must be {CHEMICAL_STATE_SCHEMA_VERSION!r}, got "
            f"{record.schema_version!r}."
        )
    for name in (
        "source_manifest_path",
        "ligand_identity_record_path",
        "ligand_coordinate_record_path",
        "binding_site_record_path",
        "chemical_preparation_record_path",
    ):
        if not getattr(record, name):
            errors.append(f"{name} is required.")
    for name in (
        "source_manifest_sha256",
        "ligand_identity_record_sha256",
        "ligand_coordinate_record_sha256",
        "binding_site_record_sha256",
        "chemical_preparation_record_sha256",
    ):
        if not _SHA256_PATTERN.fullmatch(getattr(record, name)):
            errors.append(f"{name} must be a lowercase SHA-256 digest.")
    if not record.accession:
        errors.append("accession is required.")
    if not record.toolkit or not record.toolkit_version:
        errors.append("toolkit and toolkit_version are required.")
    if not record.parameters:
        errors.append("parameters must not be empty.")
    if record.population_model is not None:
        errors.append("population_model must be null in this enumeration batch.")

    state_ids = [state.state_id for state in record.states]
    if not record.states:
        errors.append("states must not be empty.")
    if len(set(state_ids)) != len(state_ids):
        errors.append("state_id values must be unique.")
    if record.reference_state_id not in state_ids:
        errors.append("reference_state_id must name one of the states.")
    if record.selected_state_id is not None and record.selected_state_id not in state_ids:
        errors.append("selected_state_id must name one of the states.")

    by_id = {state.state_id: state for state in record.states}
    for index, state in enumerate(record.states):
        if state.kind not in STATE_KINDS:
            errors.append(f"states[{index}].kind is invalid.")
        if not state.label or not state.molecular_graph or not state.formula:
            errors.append(f"states[{index}] requires label, graph, and formula.")
        if state.proton_count < 0:
            errors.append(f"states[{index}].proton_count cannot be negative.")
        if not state.assumptions:
            errors.append(f"states[{index}].assumptions must not be empty.")
        if state.kind != "reference":
            if not state.parent_state_id or state.parent_state_id not in by_id:
                errors.append(f"states[{index}].parent_state_id must name a known state.")
            else:
                parent = by_id[state.parent_state_id]
                if state.kind == "tautomer":
                    if state.tautomer_of is None:
                        errors.append(f"states[{index}].tautomer_of is required for a tautomer.")
                    if state.proton_count != parent.proton_count:
                        errors.append(f"states[{index}] tautomer proton count must match its parent.")
                    if state.formula != parent.formula:
                        errors.append(f"states[{index}] tautomer formula must match its parent.")
                if state.kind == "protonation" and state.proton_count == parent.proton_count:
                    errors.append(f"states[{index}] protonation state must change proton count.")
        for site_index, site in enumerate(state.pka_sites):
            if site.pka is not None and not site.pka_source:
                errors.append(
                    f"states[{index}].pka_sites[{site_index}] requires a pKa source."
                )
            if site.pka_uncertainty is not None and site.pka_uncertainty < 0:
                errors.append(
                    f"states[{index}].pka_sites[{site_index}] uncertainty cannot be negative."
                )
        atom_ids = [atom.atom_id for atom in state.atoms]
        if len(set(atom_ids)) != len(atom_ids):
            errors.append(f"states[{index}] atom_ids must be unique.")
        if not state.atoms:
            errors.append(f"states[{index}].atoms must not be empty.")
        mapping_names = [atom_id for atom_id, _ in state.reference_atom_mapping]
        mapping_indices = [
            reference_index for _, reference_index in state.reference_atom_mapping
        ]
        if set(mapping_names) != set(atom_ids):
            errors.append(
                f"states[{index}].reference_atom_mapping must cover every state atom."
            )
        if len(set(mapping_indices)) != len(mapping_indices):
            errors.append(
                f"states[{index}].reference_atom_mapping indices must be unique."
            )
        if any(reference_index < 0 for reference_index in mapping_indices):
            errors.append(
                f"states[{index}].reference_atom_mapping indices cannot be negative."
            )
        for ring_index, ring in enumerate(state.rings):
            if not set(ring.atom_ids).issubset(set(atom_ids)):
                errors.append(
                    f"states[{index}].rings[{ring_index}] references an unknown atom."
                )
            if not any(ring.normal):
                errors.append(f"states[{index}].rings[{ring_index}].normal is zero.")
        for water_index, water in enumerate(state.water_roles):
            if water.role not in ("unassigned", "donor", "acceptor", "both", "unknown"):
                errors.append(
                    f"states[{index}].water_roles[{water_index}].role is invalid."
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
