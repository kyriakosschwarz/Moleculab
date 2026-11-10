"""Separate chemical preparation for selected protein pocket atoms.

This layer is intentionally conservative.  It records residue-template formal
charges and donor/acceptor/hydrophobic/aromatic assignments for the selected
protein pocket, keeps deposited coordinates separate from any future generated
hydrogens, and preserves multiple protonation hypotheses without selecting one.

It does not generate protein hydrogens, assign pH-dependent populations, infer
water/ion/cofactor roles, or estimate affinity.
"""

from __future__ import annotations

import json
import math
import re
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .chemical_preparation import PreparedWaterRole
from .structure_manifest import (
    PROVENANCE_LABELS,
    DecisionRecord,
    ProteinSelection,
    ProvenanceEntry,
    ProvenanceRecord,
    ValidationCheck,
    ValidationRecord,
)

PROTEIN_PREPARATION_SCHEMA_VERSION = "1.1"

_HYDROPHOBIC_CARBONS = {
    "ALA": {"CB"},
    "VAL": {"CB", "CG1", "CG2"},
    "LEU": {"CB", "CG", "CD1", "CD2"},
    "ILE": {"CB", "CG1", "CG2", "CD1"},
    "MET": {"CB", "CG", "CE"},
    "PHE": {"CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "TRP": {"CB", "CG", "CD1", "CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"},
    "PRO": {"CB", "CG", "CD"},
    "TYR": {"CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
}


def residue_atom_features(
    residue_name: str,
    atom_name: str,
    histidine_state: str = "HIE",
) -> tuple[int, bool, bool, bool]:
    """Return ``(formal_charge, is_donor, is_acceptor, is_hydrophobic)``.

    Backbone and side-chain roles are composed, never overwritten.  The values
    are explicit residue-template hypotheses, not experimental observations.
    """

    charge = 0
    donor = atom_name == "N" and residue_name != "PRO"
    acceptor = atom_name == "O"
    if residue_name == "ARG":
        if atom_name in {"NE", "NH1", "NH2"}:
            donor = True
        if atom_name == "CZ":
            charge = 1
    elif residue_name == "ASN":
        if atom_name == "ND2":
            donor = True
        if atom_name == "OD1":
            acceptor = True
    elif residue_name == "ASP":
        if atom_name in {"OD1", "OD2"}:
            acceptor = True
        if atom_name == "OD1":
            charge = -1
    elif residue_name == "CYS":
        if atom_name == "SG":
            donor = acceptor = True
    elif residue_name == "GLN":
        if atom_name == "NE2":
            donor = True
        if atom_name == "OE1":
            acceptor = True
    elif residue_name == "GLU":
        if atom_name in {"OE1", "OE2"}:
            acceptor = True
        if atom_name == "OE1":
            charge = -1
    elif residue_name == "HIS":
        if histidine_state == "HIP":
            if atom_name in {"ND1", "NE2"}:
                donor = True
            if atom_name == "ND1":
                charge = 1
        else:
            if atom_name == "NE2":
                donor = True
            if atom_name == "ND1":
                acceptor = True
    elif residue_name == "LYS":
        if atom_name == "NZ":
            donor = True
            charge = 1
    elif residue_name == "SER":
        if atom_name == "OG":
            donor = acceptor = True
    elif residue_name == "THR":
        if atom_name == "OG1":
            donor = acceptor = True
    elif residue_name == "TRP":
        if atom_name == "NE1":
            donor = True
    elif residue_name == "TYR":
        if atom_name == "OH":
            donor = acceptor = True
    if atom_name == "OXT":
        acceptor = True
        charge = -1
    return charge, donor, acceptor, atom_name in _HYDROPHOBIC_CARBONS.get(
        residue_name, set()
    )


COORDINATE_SOURCES = ("deposited", "generated")
WATER_ROLES = ("unassigned", "donor", "acceptor", "both", "unknown")
HETERO_ROLES = ("unknown", "ion", "cofactor", "ligand")
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


def _string_tuple(data: Mapping[str, Any], key: str) -> tuple[str, ...]:
    return tuple(str(item) for item in _as_list(data.get(key), key))


@dataclass(frozen=True)
class PreparedProteinAtom:
    atom_id: str
    atom_name: str
    element: str
    residue_name: str
    label_asym_id: str
    auth_asym_id: str
    label_seq_id: str | None
    auth_seq_id: str
    insertion_code: str | None
    alternate_location: str | None
    position_angstrom: tuple[float, float, float]
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
            "atom_name": self.atom_name,
            "element": self.element,
            "residue_name": self.residue_name,
            "label_asym_id": self.label_asym_id,
            "auth_asym_id": self.auth_asym_id,
            "label_seq_id": self.label_seq_id,
            "auth_seq_id": self.auth_seq_id,
            "insertion_code": self.insertion_code,
            "alternate_location": self.alternate_location,
            "position_angstrom": list(self.position_angstrom),
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
    def from_dict(cls, data: Mapping[str, Any]) -> PreparedProteinAtom:
        data = _as_mapping(data, "prepared protein atom")
        return cls(
            atom_id=_required_string(data, "atom_id"),
            atom_name=_required_string(data, "atom_name"),
            element=_required_string(data, "element"),
            residue_name=_required_string(data, "residue_name"),
            label_asym_id=_required_string(data, "label_asym_id"),
            auth_asym_id=_required_string(data, "auth_asym_id"),
            label_seq_id=_optional_string(data, "label_seq_id"),
            auth_seq_id=_required_string(data, "auth_seq_id"),
            insertion_code=_optional_string(data, "insertion_code"),
            alternate_location=_optional_string(data, "alternate_location"),
            position_angstrom=_vector3(data, "position_angstrom"),
            formal_charge=_optional_int(data, "formal_charge"),
            is_donor=_optional_bool(data, "is_donor"),
            is_acceptor=_optional_bool(data, "is_acceptor"),
            is_hydrophobic=_optional_bool(data, "is_hydrophobic"),
            aromatic_ring_ids=_string_tuple(data, "aromatic_ring_ids"),
            hydrogen_count=_required_int(data, "hydrogen_count"),
            hydrogen_positions_nm=tuple(
                tuple(float(item) for item in vector)
                for vector in _as_list(data.get("hydrogen_positions_nm"), "hydrogen_positions_nm")
            ),
            coordinate_source=_required_string(data, "coordinate_source"),
            hydrogen_source=_optional_string(data, "hydrogen_source"),
            feature_source=_optional_string(data, "feature_source"),
        )


@dataclass(frozen=True)
class ProteinResidueState:
    residue_id: str
    residue_name: str
    label_asym_id: str
    auth_asym_id: str
    label_seq_id: str | None
    auth_seq_id: str
    insertion_code: str | None
    protonation_state: str
    formal_charge: int | None
    source: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "residue_id": self.residue_id,
            "residue_name": self.residue_name,
            "label_asym_id": self.label_asym_id,
            "auth_asym_id": self.auth_asym_id,
            "label_seq_id": self.label_seq_id,
            "auth_seq_id": self.auth_seq_id,
            "insertion_code": self.insertion_code,
            "protonation_state": self.protonation_state,
            "formal_charge": self.formal_charge,
            "source": self.source,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> ProteinResidueState:
        data = _as_mapping(data, "protein residue state")
        return cls(
            residue_id=_required_string(data, "residue_id"),
            residue_name=_required_string(data, "residue_name"),
            label_asym_id=_required_string(data, "label_asym_id"),
            auth_asym_id=_required_string(data, "auth_asym_id"),
            label_seq_id=_optional_string(data, "label_seq_id"),
            auth_seq_id=_required_string(data, "auth_seq_id"),
            insertion_code=_optional_string(data, "insertion_code"),
            protonation_state=_required_string(data, "protonation_state"),
            formal_charge=_optional_int(data, "formal_charge"),
            source=_required_string(data, "source"),
        )


@dataclass(frozen=True)
class PreparedProteinRing:
    ring_id: str
    residue_id: str
    atom_ids: tuple[str, ...]
    centroid_nm: tuple[float, float, float]
    normal: tuple[float, float, float]
    source: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "ring_id": self.ring_id,
            "residue_id": self.residue_id,
            "atom_ids": list(self.atom_ids),
            "centroid_nm": list(self.centroid_nm),
            "normal": list(self.normal),
            "source": self.source,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> PreparedProteinRing:
        data = _as_mapping(data, "prepared protein ring")
        return cls(
            ring_id=_required_string(data, "ring_id"),
            residue_id=_required_string(data, "residue_id"),
            atom_ids=_string_tuple(data, "atom_ids"),
            centroid_nm=_vector3(data, "centroid_nm"),
            normal=_vector3(data, "normal"),
            source=_required_string(data, "source"),
        )


@dataclass(frozen=True)
class ProteinHeteroRole:
    hetero_id: str
    residue_name: str
    role: str
    source: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "hetero_id": self.hetero_id,
            "residue_name": self.residue_name,
            "role": self.role,
            "source": self.source,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> ProteinHeteroRole:
        data = _as_mapping(data, "protein hetero role")
        return cls(
            hetero_id=_required_string(data, "hetero_id"),
            residue_name=_required_string(data, "residue_name"),
            role=_required_string(data, "role"),
            source=_required_string(data, "source"),
        )


@dataclass(frozen=True)
class PreparedProteinState:
    state_id: str
    label: str
    protonation_model: str
    hydrogen_method: str | None
    assumptions: tuple[str, ...]
    atoms: tuple[PreparedProteinAtom, ...]
    residue_states: tuple[ProteinResidueState, ...]
    rings: tuple[PreparedProteinRing, ...]
    water_roles: tuple[PreparedWaterRole, ...]
    hetero_roles: tuple[ProteinHeteroRole, ...]

    def to_dict(self) -> dict[str, Any]:
        return {
            "state_id": self.state_id,
            "label": self.label,
            "protonation_model": self.protonation_model,
            "hydrogen_method": self.hydrogen_method,
            "assumptions": list(self.assumptions),
            "atoms": [atom.to_dict() for atom in self.atoms],
            "residue_states": [residue.to_dict() for residue in self.residue_states],
            "rings": [ring.to_dict() for ring in self.rings],
            "water_roles": [role.to_dict() for role in self.water_roles],
            "hetero_roles": [role.to_dict() for role in self.hetero_roles],
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> PreparedProteinState:
        data = _as_mapping(data, "prepared protein state")
        return cls(
            state_id=_required_string(data, "state_id"),
            label=_required_string(data, "label"),
            protonation_model=_required_string(data, "protonation_model"),
            hydrogen_method=_optional_string(data, "hydrogen_method"),
            assumptions=_string_tuple(data, "assumptions"),
            atoms=tuple(
                PreparedProteinAtom.from_dict(_as_mapping(item, "prepared protein atom"))
                for item in _as_list(data.get("atoms"), "atoms")
            ),
            residue_states=tuple(
                ProteinResidueState.from_dict(_as_mapping(item, "protein residue state"))
                for item in _as_list(data.get("residue_states"), "residue_states")
            ),
            rings=tuple(
                PreparedProteinRing.from_dict(_as_mapping(item, "prepared protein ring"))
                for item in _as_list(data.get("rings"), "rings")
            ),
            water_roles=tuple(
                PreparedWaterRole.from_dict(_as_mapping(item, "prepared water role"))
                for item in _as_list(data.get("water_roles"), "water_roles")
            ),
            hetero_roles=tuple(
                ProteinHeteroRole.from_dict(_as_mapping(item, "protein hetero role"))
                for item in _as_list(data.get("hetero_roles"), "hetero_roles")
            ),
        )


@dataclass(frozen=True)
class ProteinPreparationDiscrepancy:
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
    def from_dict(cls, data: Mapping[str, Any]) -> ProteinPreparationDiscrepancy:
        data = _as_mapping(data, "protein preparation discrepancy")
        return cls(
            category=_required_string(data, "category"),
            status=_required_string(data, "status"),
            detail=_required_string(data, "detail"),
            observed=_optional_string(data, "observed"),
            expected=_optional_string(data, "expected"),
            resolution=_optional_string(data, "resolution"),
        )


@dataclass(frozen=True)
class ProteinPreparationRecord:
    source_manifest_path: str
    source_manifest_sha256: str
    binding_site_record_path: str
    binding_site_record_sha256: str
    accession: str
    selection: ProteinSelection
    scope: str
    toolkit: str
    toolkit_version: str
    parameters: tuple[str, ...]
    states: tuple[PreparedProteinState, ...]
    selected_state_id: str | None
    decisions: tuple[DecisionRecord, ...]
    provenance: ProvenanceRecord
    validation: ValidationRecord
    discrepancies: tuple[ProteinPreparationDiscrepancy, ...]
    schema_version: str = PROTEIN_PREPARATION_SCHEMA_VERSION

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "source_manifest_path": self.source_manifest_path,
            "source_manifest_sha256": self.source_manifest_sha256,
            "binding_site_record_path": self.binding_site_record_path,
            "binding_site_record_sha256": self.binding_site_record_sha256,
            "accession": self.accession,
            "selection": self.selection.to_dict(),
            "scope": self.scope,
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
    def from_dict(cls, data: Mapping[str, Any]) -> ProteinPreparationRecord:
        data = _as_mapping(data, "protein preparation record")
        return cls(
            schema_version=_required_string(data, "schema_version"),
            source_manifest_path=_required_string(data, "source_manifest_path"),
            source_manifest_sha256=_required_string(data, "source_manifest_sha256"),
            binding_site_record_path=_required_string(data, "binding_site_record_path"),
            binding_site_record_sha256=_required_string(
                data, "binding_site_record_sha256"
            ),
            accession=_required_string(data, "accession"),
            selection=ProteinSelection.from_dict(
                _as_mapping(data.get("selection"), "selection")
            ),
            scope=_required_string(data, "scope"),
            toolkit=_required_string(data, "toolkit"),
            toolkit_version=_required_string(data, "toolkit_version"),
            parameters=_string_tuple(data, "parameters"),
            states=tuple(
                PreparedProteinState.from_dict(_as_mapping(item, "prepared protein state"))
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
                ProteinPreparationDiscrepancy.from_dict(
                    _as_mapping(item, "protein preparation discrepancy")
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


def build_protein_preparation_record(
    *,
    source_manifest_path: str,
    source_manifest_sha256: str,
    binding_site_record_path: str,
    binding_site_record_sha256: str,
    accession: str,
    selection: ProteinSelection,
    scope: str,
    toolkit: str,
    toolkit_version: str,
    parameters: tuple[str, ...],
    states: tuple[PreparedProteinState, ...],
    selected_state_id: str | None = None,
) -> ProteinPreparationRecord:
    if not states:
        raise ValueError("states must contain at least one prepared protein state.")
    state_ids = [state.state_id for state in states]
    if len(set(state_ids)) != len(state_ids):
        raise ValueError("state_id values must be unique.")
    if selected_state_id is not None and selected_state_id not in state_ids:
        raise ValueError("selected_state_id must name one of the states.")

    checks = (
        _check(
            "provenance",
            "protein_preparation.source_records",
            bool(source_manifest_path and binding_site_record_path),
            subject=accession,
            criterion="manifest and binding-site records are referenced",
            observed="two input record paths present",
            expected="manifest and binding-site record",
            detail="protein preparation does not replace experimental coordinates",
        ),
        _check(
            "scope",
            "protein_preparation.scope_explicit",
            bool(scope),
            subject=selection.auth_asym_id,
            criterion="the prepared protein scope is explicit",
            observed=scope,
            expected="documented atom/residue scope",
            detail="the record does not imply whole-protein preparation",
        ),
        _check(
            "states",
            "protein_preparation.states_unselected",
            selected_state_id is None,
            subject=selection.auth_asym_id,
            criterion="protonation states are not silently selected",
            observed=f"selected_state_id={selected_state_id!r}",
            expected="null unless explicitly chosen downstream",
            detail="alternative residue protonation hypotheses coexist",
        ),
        _check(
            "chemistry",
            "protein_preparation.hydrogens_not_experimental",
            all(atom.hydrogen_source is None for state in states for atom in state.atoms),
            subject=selection.auth_asym_id,
            criterion="no generated hydrogens are mistaken for experimental atoms",
            observed="no hydrogens generated",
            expected="generated hydrogens labelled separately if added",
            detail="donor/acceptor templates do not create experimental hydrogens",
        ),
        _check(
            "chemistry",
            "protein_preparation.roles_unassigned",
            all(
                role.role == "unassigned"
                for state in states
                for role in state.water_roles
            ),
            subject=selection.auth_asym_id,
            criterion="experimental waters remain chemically unassigned",
            observed="waters preserved with unassigned roles",
            expected="no automatic water-role assignment",
            detail="water-mediated interactions remain unproven",
        ),
        _check(
            "interpretation",
            "protein_preparation.no_affinity",
            True,
            subject=selection.auth_asym_id,
            criterion="protein preparation is not converted into affinity",
            observed="no affinity or free-energy estimate is produced",
            expected="chemical preparation only",
            detail="prepared features are inputs to later analyses",
        ),
    )

    discrepancies: list[ProteinPreparationDiscrepancy] = []
    if selected_state_id is None:
        discrepancies.append(
            ProteinPreparationDiscrepancy(
                "state_selection",
                "not_applicable",
                "No protein protonation state was selected.",
                observed="selected_state_id is null",
                expected="explicit selection or coexistence of alternatives",
                resolution="All explicit template states are preserved.",
            )
        )
    discrepancies.append(
        ProteinPreparationDiscrepancy(
            "hydrogen_coordinates",
            "not_applicable",
            "No protein hydrogen coordinates were generated.",
            observed="0 generated protein hydrogens",
            expected="aligned hydrogens required for angular H-bond geometry",
            resolution="Donor/acceptor templates are available, but H-bond geometry remains unproven.",
        )
    )

    decisions = (
        DecisionRecord(
            decision="prepared only the selected pocket atoms",
            actor="protein chemical preparation",
            rationale=scope,
        ),
        DecisionRecord(
            decision="used explicit residue templates for charges and features",
            actor="protein chemical preparation",
            rationale="template assignments are preparation hypotheses, not experimental observations",
        ),
        DecisionRecord(
            decision="kept multiple protonation hypotheses unselected",
            actor="protein chemical preparation",
            rationale="local protein protonation is not established by the experimental record",
        ),
        DecisionRecord(
            decision="did not generate protein hydrogens or assign water roles",
            actor="protein chemical preparation",
            rationale="hydrogen geometry and water-mediated interactions remain unresolved",
        ),
    )

    provenance = ProvenanceRecord(
        records=(
            ProvenanceEntry(
                "experimental",
                "protein pocket coordinates",
                "Coordinates come from the selected experimental chain and binding-site record.",
                actor="RCSB PDB deposition",
            ),
            ProvenanceEntry(
                "calculated",
                "protein template assignments",
                (
                    f"{len(states)} template state(s) from {toolkit} {toolkit_version}; "
                    f"parameters: {'; '.join(parameters)}."
                ),
                actor="protein chemical preparation",
            ),
            ProvenanceEntry(
                "predicted",
                "none",
                "No affinity, population, or free-energy estimate is produced.",
                actor="protein chemical preparation",
            ),
            ProvenanceEntry(
                "hypothesised",
                "protonation and feature assignments",
                "Template-based residue states are hypotheses for later analysis.",
                actor="protein chemical preparation",
            ),
        )
    )

    return ProteinPreparationRecord(
        source_manifest_path=source_manifest_path,
        source_manifest_sha256=source_manifest_sha256,
        binding_site_record_path=binding_site_record_path,
        binding_site_record_sha256=binding_site_record_sha256,
        accession=accession,
        selection=selection,
        scope=scope,
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


def write_protein_preparation(path: str | Path, record: ProteinPreparationRecord) -> Path:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(record.to_dict(), indent=2, sort_keys=False) + "\n",
        encoding="utf-8",
    )
    return destination


def load_protein_preparation(path: str | Path) -> ProteinPreparationRecord:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    return ProteinPreparationRecord.from_dict(data)


def validate_protein_preparation(record: ProteinPreparationRecord) -> list[str]:
    errors: list[str] = []
    if record.schema_version != PROTEIN_PREPARATION_SCHEMA_VERSION:
        errors.append(
            f"schema_version must be {PROTEIN_PREPARATION_SCHEMA_VERSION!r}, got "
            f"{record.schema_version!r}."
        )
    for name in ("source_manifest_path", "binding_site_record_path"):
        if not getattr(record, name):
            errors.append(f"{name} is required.")
    for name in ("source_manifest_sha256", "binding_site_record_sha256"):
        if not _SHA256_PATTERN.fullmatch(getattr(record, name)):
            errors.append(f"{name} must be a lowercase SHA-256 digest.")
    if not record.accession:
        errors.append("accession is required.")
    if not record.scope:
        errors.append("scope is required.")
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
        if not state.label or not state.protonation_model or not state.assumptions:
            errors.append(f"states[{state_index}] requires label, model, and assumptions.")
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
            if len(atom.position_angstrom) != 3 or not all(
                math.isfinite(value) for value in atom.position_angstrom
            ):
                errors.append(
                    f"states[{state_index}].atoms[{atom_index}].position_angstrom must be a finite three-vector."
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
            for hydrogen_index, position in enumerate(atom.hydrogen_positions_nm):
                if len(position) != 3 or not all(
                    math.isfinite(value) for value in position
                ):
                    errors.append(
                        f"states[{state_index}].atoms[{atom_index}]."
                        f"hydrogen_positions_nm[{hydrogen_index}] must be a finite "
                        "three-vector."
                    )
        residue_ids = [residue.residue_id for residue in state.residue_states]
        if len(set(residue_ids)) != len(residue_ids):
            errors.append(f"states[{state_index}] residue_ids must be unique.")
        for ring_index, ring in enumerate(state.rings):
            if not set(ring.atom_ids).issubset(atom_id_set):
                errors.append(
                    f"states[{state_index}].rings[{ring_index}] references an unknown atom."
                )
            if not any(ring.normal):
                errors.append(f"states[{state_index}].rings[{ring_index}].normal is zero.")
        for water_index, water in enumerate(state.water_roles):
            if water.role not in WATER_ROLES:
                errors.append(
                    f"states[{state_index}].water_roles[{water_index}].role is invalid."
                )
        for hetero_index, hetero in enumerate(state.hetero_roles):
            if hetero.role not in HETERO_ROLES:
                errors.append(
                    f"states[{state_index}].hetero_roles[{hetero_index}].role is invalid."
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
