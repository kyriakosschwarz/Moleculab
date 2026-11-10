"""Deposited-water evidence and geometric bridge plausibility.

This layer preserves experimental water identity, coordinates, occupancy, and
B-factors.  It leaves water orientation and donor/acceptor roles unassigned and
reports only whether the geometry is compatible with a water bridge.

A geometrically plausible bridge is not a chemically classified water-mediated
interaction.
"""

from __future__ import annotations

import json
import math
import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .structure_manifest import (
    PROVENANCE_LABELS,
    DecisionRecord,
    LigandSelection,
    ProteinSelection,
    ProvenanceEntry,
    ProvenanceRecord,
    ValidationCheck,
    ValidationRecord,
)

WATER_EVIDENCE_SCHEMA_VERSION = "1.0"
BRIDGE_STATUSES = ("geometrically_plausible", "not_supported", "unknown")
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


def _required_float(data: Mapping[str, Any], key: str) -> float:
    value = data.get(key)
    if value is None:
        raise ValueError(f"{key} is required.")
    return float(value)


def _optional_float(data: Mapping[str, Any], key: str) -> float | None:
    value = data.get(key)
    return None if value is None else float(value)


def _vector3(data: Mapping[str, Any], key: str) -> tuple[float, float, float]:
    value = _as_list(data.get(key), key)
    if len(value) != 3:
        raise ValueError(f"{key} must contain three numbers.")
    return (float(value[0]), float(value[1]), float(value[2]))


@dataclass(frozen=True)
class WaterAtomRecord:
    water_id: str
    label_asym_id: str
    auth_asym_id: str
    label_seq_id: str | None
    auth_seq_id: str
    atom_name: str
    position_angstrom: tuple[float, float, float]
    occupancy: float | None
    b_iso_or_equiv: float | None
    orientation_known: bool
    donor_role: str | None
    acceptor_role: str | None
    source: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "water_id": self.water_id,
            "label_asym_id": self.label_asym_id,
            "auth_asym_id": self.auth_asym_id,
            "label_seq_id": self.label_seq_id,
            "auth_seq_id": self.auth_seq_id,
            "atom_name": self.atom_name,
            "position_angstrom": list(self.position_angstrom),
            "occupancy": self.occupancy,
            "b_iso_or_equiv": self.b_iso_or_equiv,
            "orientation_known": self.orientation_known,
            "donor_role": self.donor_role,
            "acceptor_role": self.acceptor_role,
            "source": self.source,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> WaterAtomRecord:
        data = _as_mapping(data, "water atom record")
        return cls(
            water_id=_required_string(data, "water_id"),
            label_asym_id=_required_string(data, "label_asym_id"),
            auth_asym_id=_required_string(data, "auth_asym_id"),
            label_seq_id=_optional_string(data, "label_seq_id"),
            auth_seq_id=_required_string(data, "auth_seq_id"),
            atom_name=_required_string(data, "atom_name"),
            position_angstrom=_vector3(data, "position_angstrom"),
            occupancy=_optional_float(data, "occupancy"),
            b_iso_or_equiv=_optional_float(data, "b_iso_or_equiv"),
            orientation_known=bool(data.get("orientation_known")),
            donor_role=_optional_string(data, "donor_role"),
            acceptor_role=_optional_string(data, "acceptor_role"),
            source=_required_string(data, "source"),
        )


@dataclass(frozen=True)
class WaterBridgeEvidence:
    water_id: str
    ligand_atom_id: str | None
    protein_atom_id: str | None
    ligand_distance_angstrom: float | None
    protein_distance_angstrom: float | None
    status: str
    chemical_status: str
    criteria: str
    provenance: str
    limitations: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "water_id": self.water_id,
            "ligand_atom_id": self.ligand_atom_id,
            "protein_atom_id": self.protein_atom_id,
            "ligand_distance_angstrom": self.ligand_distance_angstrom,
            "protein_distance_angstrom": self.protein_distance_angstrom,
            "status": self.status,
            "chemical_status": self.chemical_status,
            "criteria": self.criteria,
            "provenance": self.provenance,
            "limitations": self.limitations,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> WaterBridgeEvidence:
        data = _as_mapping(data, "water bridge evidence")
        return cls(
            water_id=_required_string(data, "water_id"),
            ligand_atom_id=_optional_string(data, "ligand_atom_id"),
            protein_atom_id=_optional_string(data, "protein_atom_id"),
            ligand_distance_angstrom=_optional_float(data, "ligand_distance_angstrom"),
            protein_distance_angstrom=_optional_float(data, "protein_distance_angstrom"),
            status=_required_string(data, "status"),
            chemical_status=_required_string(data, "chemical_status"),
            criteria=_required_string(data, "criteria"),
            provenance=_required_string(data, "provenance"),
            limitations=_required_string(data, "limitations"),
        )


@dataclass(frozen=True)
class WaterEvidenceDiscrepancy:
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
    def from_dict(cls, data: Mapping[str, Any]) -> WaterEvidenceDiscrepancy:
        data = _as_mapping(data, "water evidence discrepancy")
        return cls(
            category=_required_string(data, "category"),
            status=_required_string(data, "status"),
            detail=_required_string(data, "detail"),
            observed=_optional_string(data, "observed"),
            expected=_optional_string(data, "expected"),
            resolution=_optional_string(data, "resolution"),
        )


@dataclass(frozen=True)
class WaterEvidenceRecord:
    source_manifest_path: str
    source_manifest_sha256: str
    binding_site_record_path: str
    binding_site_record_sha256: str
    ligand_identity_record_path: str
    ligand_identity_record_sha256: str
    chemical_preparation_record_path: str
    chemical_preparation_record_sha256: str
    protein_preparation_record_path: str
    protein_preparation_record_sha256: str
    accession: str
    protein_selection: ProteinSelection
    ligand_selection: LigandSelection
    ligand_state_id: str
    protein_state_id: str
    bridge_cutoff_angstrom: float
    waters: tuple[WaterAtomRecord, ...]
    bridges: tuple[WaterBridgeEvidence, ...]
    decisions: tuple[DecisionRecord, ...]
    provenance: ProvenanceRecord
    validation: ValidationRecord
    discrepancies: tuple[WaterEvidenceDiscrepancy, ...]
    schema_version: str = WATER_EVIDENCE_SCHEMA_VERSION

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "source_manifest_path": self.source_manifest_path,
            "source_manifest_sha256": self.source_manifest_sha256,
            "binding_site_record_path": self.binding_site_record_path,
            "binding_site_record_sha256": self.binding_site_record_sha256,
            "ligand_identity_record_path": self.ligand_identity_record_path,
            "ligand_identity_record_sha256": self.ligand_identity_record_sha256,
            "chemical_preparation_record_path": self.chemical_preparation_record_path,
            "chemical_preparation_record_sha256": self.chemical_preparation_record_sha256,
            "protein_preparation_record_path": self.protein_preparation_record_path,
            "protein_preparation_record_sha256": self.protein_preparation_record_sha256,
            "accession": self.accession,
            "protein_selection": self.protein_selection.to_dict(),
            "ligand_selection": self.ligand_selection.to_dict(),
            "ligand_state_id": self.ligand_state_id,
            "protein_state_id": self.protein_state_id,
            "bridge_cutoff_angstrom": self.bridge_cutoff_angstrom,
            "waters": [item.to_dict() for item in self.waters],
            "bridges": [item.to_dict() for item in self.bridges],
            "decisions": [decision.to_dict() for decision in self.decisions],
            "provenance": self.provenance.to_dict(),
            "validation": self.validation.to_dict(),
            "discrepancies": [item.to_dict() for item in self.discrepancies],
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> WaterEvidenceRecord:
        data = _as_mapping(data, "water evidence record")
        return cls(
            schema_version=_required_string(data, "schema_version"),
            source_manifest_path=_required_string(data, "source_manifest_path"),
            source_manifest_sha256=_required_string(data, "source_manifest_sha256"),
            binding_site_record_path=_required_string(data, "binding_site_record_path"),
            binding_site_record_sha256=_required_string(
                data, "binding_site_record_sha256"
            ),
            ligand_identity_record_path=_required_string(
                data, "ligand_identity_record_path"
            ),
            ligand_identity_record_sha256=_required_string(
                data, "ligand_identity_record_sha256"
            ),
            chemical_preparation_record_path=_required_string(
                data, "chemical_preparation_record_path"
            ),
            chemical_preparation_record_sha256=_required_string(
                data, "chemical_preparation_record_sha256"
            ),
            protein_preparation_record_path=_required_string(
                data, "protein_preparation_record_path"
            ),
            protein_preparation_record_sha256=_required_string(
                data, "protein_preparation_record_sha256"
            ),
            accession=_required_string(data, "accession"),
            protein_selection=ProteinSelection.from_dict(
                _as_mapping(data.get("protein_selection"), "protein_selection")
            ),
            ligand_selection=LigandSelection.from_dict(
                _as_mapping(data.get("ligand_selection"), "ligand_selection")
            ),
            ligand_state_id=_required_string(data, "ligand_state_id"),
            protein_state_id=_required_string(data, "protein_state_id"),
            bridge_cutoff_angstrom=_required_float(data, "bridge_cutoff_angstrom"),
            waters=tuple(
                WaterAtomRecord.from_dict(_as_mapping(item, "water atom record"))
                for item in _as_list(data.get("waters"), "waters")
            ),
            bridges=tuple(
                WaterBridgeEvidence.from_dict(
                    _as_mapping(item, "water bridge evidence")
                )
                for item in _as_list(data.get("bridges"), "bridges")
            ),
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
                WaterEvidenceDiscrepancy.from_dict(
                    _as_mapping(item, "water evidence discrepancy")
                )
                for item in _as_list(data.get("discrepancies"), "discrepancies")
            ),
        )


def _distance(first, second) -> float:
    return math.dist(
        tuple(float(value) for value in first),
        tuple(float(value) for value in second),
    )


def assess_water_bridges(
    waters: Sequence[WaterAtomRecord],
    ligand_polar_atoms: Sequence[Any],
    protein_polar_atoms: Sequence[Any],
    cutoff_angstrom: float,
) -> tuple[WaterBridgeEvidence, ...]:
    """Report geometric bridge plausibility without assigning water chemistry."""

    records = []
    for water in waters:
        ligand = min(
            ligand_polar_atoms,
            key=lambda atom: _distance(water.position_angstrom, atom.position_angstrom),
            default=None,
        )
        protein = min(
            protein_polar_atoms,
            key=lambda atom: _distance(water.position_angstrom, atom.position_angstrom),
            default=None,
        )
        ligand_distance = (
            None
            if ligand is None
            else _distance(water.position_angstrom, ligand.position_angstrom)
        )
        protein_distance = (
            None
            if protein is None
            else _distance(water.position_angstrom, protein.position_angstrom)
        )
        if ligand_distance is None or protein_distance is None:
            status = "unknown"
            limitations = "ligand or protein polar atom evidence is unavailable"
        elif ligand_distance <= cutoff_angstrom and protein_distance <= cutoff_angstrom:
            status = "geometrically_plausible"
            limitations = "water orientation and chemical roles are unassigned"
        else:
            status = "not_supported"
            limitations = "one side exceeds the bridge distance cutoff"
        records.append(
            WaterBridgeEvidence(
                water_id=water.water_id,
                ligand_atom_id=None if ligand is None else ligand.atom_id,
                protein_atom_id=None if protein is None else protein.atom_id,
                ligand_distance_angstrom=ligand_distance,
                protein_distance_angstrom=protein_distance,
                status=status,
                chemical_status="not_assigned",
                criteria=f"both water contacts <= {cutoff_angstrom} A",
                provenance="experimental water coordinates; template polar features",
                limitations=limitations,
            )
        )
    return tuple(records)


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


def build_water_evidence_record(
    *,
    source_manifest_path: str,
    source_manifest_sha256: str,
    binding_site_record_path: str,
    binding_site_record_sha256: str,
    ligand_identity_record_path: str,
    ligand_identity_record_sha256: str,
    chemical_preparation_record_path: str,
    chemical_preparation_record_sha256: str,
    protein_preparation_record_path: str,
    protein_preparation_record_sha256: str,
    accession: str,
    protein_selection: ProteinSelection,
    ligand_selection: LigandSelection,
    ligand_state_id: str,
    protein_state_id: str,
    bridge_cutoff_angstrom: float,
    waters: tuple[WaterAtomRecord, ...],
    bridges: tuple[WaterBridgeEvidence, ...],
) -> WaterEvidenceRecord:
    checks = (
        _check(
            "provenance",
            "water.source_records",
            True,
            subject=accession,
            criterion="experimental and preparation records are referenced",
            observed="records present",
            expected="manifest, binding site, ligand, and protein preparation",
            detail="water coordinates remain experimental observations",
        ),
        _check(
            "water_observation",
            "water.identity_and_coordinates",
            all(water.water_id and water.auth_seq_id for water in waters),
            subject=accession,
            criterion="every water has a stable residue identity and coordinates",
            observed=f"{len(waters)} waters",
            expected="one identity per deposited water",
            detail="waters are not merged or renumbered",
        ),
        _check(
            "water_chemistry",
            "water.roles_unassigned",
            all(
                water.donor_role is None
                and water.acceptor_role is None
                and not water.orientation_known
                for water in waters
            ),
            subject=accession,
            criterion="water orientation and roles remain unknown",
            observed="unassigned",
            expected="explicit assignment before chemical classification",
            detail="a nearby water is not automatically a donor or acceptor",
        ),
        _check(
            "water_geometry",
            "water.bridge_geometry_only",
            all(
                bridge.chemical_status == "not_assigned" for bridge in bridges
            ),
            subject=accession,
            criterion="bridges are geometric evidence only",
            observed=f"{len(bridges)} geometric bridge records",
            expected="no chemical water-mediated classification",
            detail="geometric plausibility is not a classified interaction",
        ),
        _check(
            "interpretation",
            "water.no_affinity",
            True,
            subject=accession,
            criterion="water evidence is not converted into affinity",
            observed="no affinity or free-energy estimate is produced",
            expected="structural evidence only",
            detail="water displacement thermodynamics is out of scope",
        ),
    )
    discrepancies = (
        WaterEvidenceDiscrepancy(
            "water_roles",
            "not_applicable",
            "Water orientations and chemical roles were not assigned.",
            observed="unassigned",
            expected="explicit water geometry and role evidence",
            resolution="Only geometric bridge plausibility is reported.",
        ),
    )

    decisions = (
        DecisionRecord(
            decision="preserved deposited water identity and coordinates",
            actor="water evidence preparation",
            rationale="water residue identifiers and coordinates are experimental observations",
        ),
        DecisionRecord(
            decision="kept water orientation and roles unassigned",
            actor="water evidence preparation",
            rationale="proximity alone does not establish donor/acceptor roles",
        ),
        DecisionRecord(
            decision="reported bridge geometry without chemical classification",
            actor="water evidence preparation",
            rationale="a geometrically plausible bridge is not a classified interaction",
        ),
        DecisionRecord(
            decision="did not convert water evidence into affinity",
            actor="water evidence preparation",
            rationale="water displacement thermodynamics is not established",
        ),
    )

    provenance = ProvenanceRecord(
        records=(
            ProvenanceEntry(
                "experimental",
                "deposited water coordinates",
                f"{len(waters)} waters preserved with occupancy and B-factors.",
                actor="RCSB PDB deposition",
            ),
            ProvenanceEntry(
                "calculated",
                "geometric bridge plausibility",
                f"{len(bridges)} water bridge records under a {bridge_cutoff_angstrom} A cutoff.",
                actor="water evidence preparation",
            ),
            ProvenanceEntry(
                "predicted",
                "none",
                "No affinity or free-energy estimate is produced.",
                actor="water evidence preparation",
            ),
            ProvenanceEntry(
                "hypothesised",
                "water roles",
                "Water donor/acceptor roles remain hypotheses and are left unassigned.",
                actor="water evidence preparation",
            ),
        )
    )

    return WaterEvidenceRecord(
        source_manifest_path=source_manifest_path,
        source_manifest_sha256=source_manifest_sha256,
        binding_site_record_path=binding_site_record_path,
        binding_site_record_sha256=binding_site_record_sha256,
        ligand_identity_record_path=ligand_identity_record_path,
        ligand_identity_record_sha256=ligand_identity_record_sha256,
        chemical_preparation_record_path=chemical_preparation_record_path,
        chemical_preparation_record_sha256=chemical_preparation_record_sha256,
        protein_preparation_record_path=protein_preparation_record_path,
        protein_preparation_record_sha256=protein_preparation_record_sha256,
        accession=accession,
        protein_selection=protein_selection,
        ligand_selection=ligand_selection,
        ligand_state_id=ligand_state_id,
        protein_state_id=protein_state_id,
        bridge_cutoff_angstrom=bridge_cutoff_angstrom,
        waters=waters,
        bridges=bridges,
        decisions=decisions,
        provenance=provenance,
        validation=ValidationRecord(checks=checks),
        discrepancies=discrepancies,
    )


def write_water_evidence(path: str | Path, record: WaterEvidenceRecord) -> Path:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(record.to_dict(), indent=2, sort_keys=False) + "\n",
        encoding="utf-8",
    )
    return destination


def load_water_evidence(path: str | Path) -> WaterEvidenceRecord:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    return WaterEvidenceRecord.from_dict(data)


def validate_water_evidence(record: WaterEvidenceRecord) -> list[str]:
    errors: list[str] = []
    if record.schema_version != WATER_EVIDENCE_SCHEMA_VERSION:
        errors.append(
            f"schema_version must be {WATER_EVIDENCE_SCHEMA_VERSION!r}, got "
            f"{record.schema_version!r}."
        )
    for name in (
        "source_manifest_path",
        "binding_site_record_path",
        "ligand_identity_record_path",
        "chemical_preparation_record_path",
        "protein_preparation_record_path",
    ):
        if not getattr(record, name):
            errors.append(f"{name} is required.")
    for name in (
        "source_manifest_sha256",
        "binding_site_record_sha256",
        "ligand_identity_record_sha256",
        "chemical_preparation_record_sha256",
        "protein_preparation_record_sha256",
    ):
        if not _SHA256_PATTERN.fullmatch(getattr(record, name)):
            errors.append(f"{name} must be a lowercase SHA-256 digest.")
    if not record.accession:
        errors.append("accession is required.")
    if record.bridge_cutoff_angstrom <= 0:
        errors.append("bridge_cutoff_angstrom must be positive.")

    water_ids = [water.water_id for water in record.waters]
    if len(set(water_ids)) != len(water_ids):
        errors.append("water_id values must be unique.")
    for index, water in enumerate(record.waters):
        if water.donor_role is not None or water.acceptor_role is not None:
            errors.append(f"waters[{index}] roles must be unassigned in this batch.")
        if water.orientation_known:
            errors.append(f"waters[{index}].orientation_known must be false in this batch.")
        if len(water.position_angstrom) != 3 or not all(
            math.isfinite(value) for value in water.position_angstrom
        ):
            errors.append(f"waters[{index}] position must be finite.")
        if water.occupancy is not None and not 0 <= water.occupancy <= 1:
            errors.append(f"waters[{index}] occupancy must be in [0, 1].")

    for index, bridge in enumerate(record.bridges):
        if bridge.status not in BRIDGE_STATUSES:
            errors.append(f"bridges[{index}].status is invalid.")
        if bridge.chemical_status != "not_assigned":
            errors.append(f"bridges[{index}].chemical_status must be not_assigned.")
        if bridge.water_id not in set(water_ids):
            errors.append(f"bridges[{index}] references an unknown water.")

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

    check_ids = set()
    for index, check in enumerate(record.validation.checks):
        if check.status not in ("pass", "fail", "not_run", "not_applicable"):
            errors.append(f"validation.checks[{index}] has an invalid status.")
        if not check.check_id or not check.category or not check.detail:
            errors.append(f"validation.checks[{index}] requires check_id, category, and detail.")
        if check.check_id in check_ids:
            errors.append(f"validation.checks[{index}].check_id is duplicated.")
        check_ids.add(check.check_id)

    for index, item in enumerate(record.discrepancies):
        if item.status not in ("open", "resolved", "not_applicable"):
            errors.append(f"discrepancies[{index}] has an invalid status.")

    return errors
