# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Hydrogen-position and hydrogen-bond evidence records.

This layer preserves the distinction between deposited heavy atoms, generated
hydrogens, and a geometric hydrogen-bond assessment.  It does not generate
hydrogens itself; that is a separate preparation adapter.

For 2HYY, ligand hydrogens come from the Batch 7 reference-neutral state and are
explicitly labelled generated.  Protein hydrogen generation is not performed,
so protein-donor hydrogen bonds remain unavailable.
"""

from __future__ import annotations

import json
import math
import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from ..analysis.interactions import HydrogenBondCriteria
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

HYDROGEN_EVIDENCE_SCHEMA_VERSION = "1.0"
BOND_STATUSES = ("supported", "not_supported", "unknown")
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


def _required_int(data: Mapping[str, Any], key: str) -> int:
    value = data.get(key)
    if value is None:
        raise ValueError(f"{key} is required.")
    return int(value)


def _vector3(data: Mapping[str, Any], key: str) -> tuple[float, float, float]:
    value = _as_list(data.get(key), key)
    if len(value) != 3:
        raise ValueError(f"{key} must contain three numbers.")
    return (float(value[0]), float(value[1]), float(value[2]))


def _string_tuple(data: Mapping[str, Any], key: str) -> tuple[str, ...]:
    return tuple(str(item) for item in _as_list(data.get(key), key))


@dataclass(frozen=True)
class HydrogenPosition:
    hydrogen_id: str
    parent_atom_id: str
    parent_side: str
    parent_element: str
    position_angstrom: tuple[float, float, float]
    parent_position_angstrom: tuple[float, float, float]
    donor_hydrogen_distance_angstrom: float
    state_id: str
    generation_method: str
    toolkit: str
    toolkit_version: str
    coordinate_source: str
    parent_coordinate_source: str
    feature_source: str | None

    def to_dict(self) -> dict[str, Any]:
        return {
            "hydrogen_id": self.hydrogen_id,
            "parent_atom_id": self.parent_atom_id,
            "parent_side": self.parent_side,
            "parent_element": self.parent_element,
            "position_angstrom": list(self.position_angstrom),
            "parent_position_angstrom": list(self.parent_position_angstrom),
            "donor_hydrogen_distance_angstrom": self.donor_hydrogen_distance_angstrom,
            "state_id": self.state_id,
            "generation_method": self.generation_method,
            "toolkit": self.toolkit,
            "toolkit_version": self.toolkit_version,
            "coordinate_source": self.coordinate_source,
            "parent_coordinate_source": self.parent_coordinate_source,
            "feature_source": self.feature_source,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> HydrogenPosition:
        data = _as_mapping(data, "hydrogen position")
        return cls(
            hydrogen_id=_required_string(data, "hydrogen_id"),
            parent_atom_id=_required_string(data, "parent_atom_id"),
            parent_side=_required_string(data, "parent_side"),
            parent_element=_required_string(data, "parent_element"),
            position_angstrom=_vector3(data, "position_angstrom"),
            parent_position_angstrom=_vector3(data, "parent_position_angstrom"),
            donor_hydrogen_distance_angstrom=_required_float(
                data, "donor_hydrogen_distance_angstrom"
            ),
            state_id=_required_string(data, "state_id"),
            generation_method=_required_string(data, "generation_method"),
            toolkit=_required_string(data, "toolkit"),
            toolkit_version=_required_string(data, "toolkit_version"),
            coordinate_source=_required_string(data, "coordinate_source"),
            parent_coordinate_source=_required_string(data, "parent_coordinate_source"),
            feature_source=_optional_string(data, "feature_source"),
        )


@dataclass(frozen=True)
class HydrogenBondEvidence:
    donor_atom_id: str
    hydrogen_id: str
    acceptor_atom_id: str
    donor_side: str
    acceptor_side: str
    donor_acceptor_distance_angstrom: float
    donor_hydrogen_distance_angstrom: float
    angle_degrees: float
    status: str
    criteria: str
    provenance: str
    limitations: str | None

    def to_dict(self) -> dict[str, Any]:
        return {
            "donor_atom_id": self.donor_atom_id,
            "hydrogen_id": self.hydrogen_id,
            "acceptor_atom_id": self.acceptor_atom_id,
            "donor_side": self.donor_side,
            "acceptor_side": self.acceptor_side,
            "donor_acceptor_distance_angstrom": self.donor_acceptor_distance_angstrom,
            "donor_acceptor_distance_nm": self.donor_acceptor_distance_angstrom * 0.1,
            "donor_hydrogen_distance_angstrom": self.donor_hydrogen_distance_angstrom,
            "angle_degrees": self.angle_degrees,
            "status": self.status,
            "criteria": self.criteria,
            "provenance": self.provenance,
            "limitations": self.limitations,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> HydrogenBondEvidence:
        data = _as_mapping(data, "hydrogen bond evidence")
        return cls(
            donor_atom_id=_required_string(data, "donor_atom_id"),
            hydrogen_id=_required_string(data, "hydrogen_id"),
            acceptor_atom_id=_required_string(data, "acceptor_atom_id"),
            donor_side=_required_string(data, "donor_side"),
            acceptor_side=_required_string(data, "acceptor_side"),
            donor_acceptor_distance_angstrom=_required_float(
                data, "donor_acceptor_distance_angstrom"
            ),
            donor_hydrogen_distance_angstrom=_required_float(
                data, "donor_hydrogen_distance_angstrom"
            ),
            angle_degrees=_required_float(data, "angle_degrees"),
            status=_required_string(data, "status"),
            criteria=_required_string(data, "criteria"),
            provenance=_required_string(data, "provenance"),
            limitations=_optional_string(data, "limitations"),
        )


@dataclass(frozen=True)
class HydrogenEvidenceDiscrepancy:
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
    def from_dict(cls, data: Mapping[str, Any]) -> HydrogenEvidenceDiscrepancy:
        data = _as_mapping(data, "hydrogen evidence discrepancy")
        return cls(
            category=_required_string(data, "category"),
            status=_required_string(data, "status"),
            detail=_required_string(data, "detail"),
            observed=_optional_string(data, "observed"),
            expected=_optional_string(data, "expected"),
            resolution=_optional_string(data, "resolution"),
        )


@dataclass(frozen=True)
class HydrogenEvidenceRecord:
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
    protein_preparation_record_path: str
    protein_preparation_record_sha256: str
    accession: str
    protein_selection: ProteinSelection
    ligand_selection: LigandSelection
    ligand_state_id: str
    protein_state_id: str
    criteria: HydrogenBondCriteria
    ligand_hydrogens: tuple[HydrogenPosition, ...]
    protein_hydrogen_status: str
    bonds: tuple[HydrogenBondEvidence, ...]
    decisions: tuple[DecisionRecord, ...]
    provenance: ProvenanceRecord
    validation: ValidationRecord
    discrepancies: tuple[HydrogenEvidenceDiscrepancy, ...]
    schema_version: str = HYDROGEN_EVIDENCE_SCHEMA_VERSION

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
            "protein_preparation_record_path": self.protein_preparation_record_path,
            "protein_preparation_record_sha256": self.protein_preparation_record_sha256,
            "accession": self.accession,
            "protein_selection": self.protein_selection.to_dict(),
            "ligand_selection": self.ligand_selection.to_dict(),
            "ligand_state_id": self.ligand_state_id,
            "protein_state_id": self.protein_state_id,
            "criteria": {
                "max_distance_nm": self.criteria.max_distance,
                "min_angle_degrees": self.criteria.min_angle_degrees,
            },
            "ligand_hydrogens": [item.to_dict() for item in self.ligand_hydrogens],
            "protein_hydrogen_status": self.protein_hydrogen_status,
            "bonds": [item.to_dict() for item in self.bonds],
            "decisions": [decision.to_dict() for decision in self.decisions],
            "provenance": self.provenance.to_dict(),
            "validation": self.validation.to_dict(),
            "discrepancies": [item.to_dict() for item in self.discrepancies],
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> HydrogenEvidenceRecord:
        data = _as_mapping(data, "hydrogen evidence record")
        criteria = _as_mapping(data.get("criteria"), "criteria")
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
            criteria=HydrogenBondCriteria(
                max_distance=_required_float(criteria, "max_distance_nm"),
                min_angle_degrees=_required_float(criteria, "min_angle_degrees"),
            ),
            ligand_hydrogens=tuple(
                HydrogenPosition.from_dict(_as_mapping(item, "hydrogen position"))
                for item in _as_list(data.get("ligand_hydrogens"), "ligand_hydrogens")
            ),
            protein_hydrogen_status=_required_string(
                data, "protein_hydrogen_status"
            ),
            bonds=tuple(
                HydrogenBondEvidence.from_dict(
                    _as_mapping(item, "hydrogen bond evidence")
                )
                for item in _as_list(data.get("bonds"), "bonds")
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
                HydrogenEvidenceDiscrepancy.from_dict(
                    _as_mapping(item, "hydrogen evidence discrepancy")
                )
                for item in _as_list(data.get("discrepancies"), "discrepancies")
            ),
        )


def _distance(first, second) -> float:
    return float(
        math.dist(
            tuple(float(value) for value in first),
            tuple(float(value) for value in second),
        )
    )


def _angle(first, vertex, second) -> float:
    left = [float(first[i]) - float(vertex[i]) for i in range(3)]
    right = [float(second[i]) - float(vertex[i]) for i in range(3)]
    left_norm = math.sqrt(sum(value * value for value in left))
    right_norm = math.sqrt(sum(value * value for value in right))
    if left_norm == 0 or right_norm == 0:
        return 0.0
    cosine = sum(a * b for a, b in zip(left, right, strict=True)) / (
        left_norm * right_norm
    )
    return math.degrees(math.acos(max(-1.0, min(1.0, cosine))))


def assess_ligand_donor_hydrogen_bonds(
    hydrogens: Sequence[HydrogenPosition],
    acceptors: Sequence[Any],
    criteria: HydrogenBondCriteria,
) -> tuple[HydrogenBondEvidence, ...]:
    """Evaluate ligand-donor hydrogens against explicit acceptor features."""

    records = []
    for hydrogen in hydrogens:
        for acceptor in acceptors:
            distance = _distance(
                hydrogen.parent_position_angstrom, acceptor.position_angstrom
            )
            if distance > criteria.max_distance * 10.0:
                continue
            angle = _angle(
                hydrogen.parent_position_angstrom,
                hydrogen.position_angstrom,
                acceptor.position_angstrom,
            )
            supported = angle >= criteria.min_angle_degrees
            records.append(
                HydrogenBondEvidence(
                    donor_atom_id=hydrogen.parent_atom_id,
                    hydrogen_id=hydrogen.hydrogen_id,
                    acceptor_atom_id=acceptor.atom_id,
                    donor_side="ligand",
                    acceptor_side="protein",
                    donor_acceptor_distance_angstrom=distance,
                    donor_hydrogen_distance_angstrom=hydrogen.donor_hydrogen_distance_angstrom,
                    angle_degrees=angle,
                    status="supported" if supported else "not_supported",
                    criteria=(
                        f"D-A <= {criteria.max_distance*10.0:.2f} A; "
                        f"D-H...A >= {criteria.min_angle_degrees:.1f} degrees"
                    ),
                    provenance=(
                        "generated ligand hydrogens; protein acceptor template"
                    ),
                    limitations=(
                        "ligand hydrogens are generated, not experimentally observed"
                    ),
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


def build_hydrogen_evidence_record(
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
    protein_preparation_record_path: str,
    protein_preparation_record_sha256: str,
    accession: str,
    protein_selection: ProteinSelection,
    ligand_selection: LigandSelection,
    ligand_state_id: str,
    protein_state_id: str,
    criteria: HydrogenBondCriteria,
    ligand_hydrogens: tuple[HydrogenPosition, ...],
    protein_hydrogen_status: str,
    bonds: tuple[HydrogenBondEvidence, ...],
) -> HydrogenEvidenceRecord:
    minimum_dh = min(
        (h.donor_hydrogen_distance_angstrom for h in ligand_hydrogens),
        default=0.0,
    )
    maximum_dh = max(
        (h.donor_hydrogen_distance_angstrom for h in ligand_hydrogens),
        default=0.0,
    )
    checks = (
        _check(
            "provenance",
            "hydrogen.source_records",
            True,
            subject=accession,
            criterion="experimental and preparation records are referenced",
            observed="records present",
            expected="manifest, coordinates, binding site, and preparation records",
            detail="generated hydrogens are separate from deposited heavy atoms",
        ),
        _check(
            "hydrogen_provenance",
            "hydrogen.generated_not_experimental",
            all(
                hydrogen.coordinate_source == "generated"
                and hydrogen.parent_coordinate_source == "deposited"
                for hydrogen in ligand_hydrogens
            ),
            subject=ligand_state_id,
            criterion="generated hydrogens have generated provenance",
            observed=f"{len(ligand_hydrogens)} ligand hydrogens",
            expected="generated hydrogens labelled as such",
            detail="parent heavy atoms retain deposited provenance",
        ),
        _check(
            "hydrogen_provenance",
            "hydrogen.donor_connectivity",
            all(
                0.05 <= hydrogen.donor_hydrogen_distance_angstrom <= 0.15
                for hydrogen in ligand_hydrogens
            ),
            subject=ligand_state_id,
            criterion="donor-hydrogen distances are in a covalent range",
            observed=f"{minimum_dh:.3f}-{maximum_dh:.3f} A",
            expected="0.05-0.15 A",
            detail="connectivity is checked before angular geometry",
        ),
        _check(
            "protein_hydrogens",
            "hydrogen.protein_hydrogens_available",
            protein_hydrogen_status == "generated",
            subject=protein_state_id,
            criterion="protein donor hydrogens are available",
            observed=protein_hydrogen_status,
            expected="generated or explicitly unavailable",
            detail="protein-donor hydrogen bonds remain unavailable without hydrogens",
            not_run=protein_hydrogen_status != "generated",
        ),
        _check(
            "interpretation",
            "hydrogen.no_affinity",
            True,
            subject=accession,
            criterion="hydrogen evidence is not converted into affinity",
            observed="no affinity or free-energy estimate is produced",
            expected="geometric hydrogen-bond evidence only",
            detail="generated hydrogen positions are not experimental observations",
        ),
    )
    discrepancies = [
        HydrogenEvidenceDiscrepancy(
            "protein_hydrogens",
            "not_applicable",
            "Protein hydrogens were not generated.",
            observed=protein_hydrogen_status,
            expected="generated protein hydrogens for protein-donor H-bonds",
            resolution="Protein-donor hydrogen bonds remain unknown.",
        )
    ]

    decisions = (
        DecisionRecord(
            decision="kept generated and deposited coordinates separate",
            actor="hydrogen evidence preparation",
            rationale="ligand hydrogens are generated; heavy-atom parent coordinates are deposited",
        ),
        DecisionRecord(
            decision="validated donor-hydrogen connectivity before geometry",
            actor="hydrogen evidence preparation",
            rationale="D-H distances are checked before D-H...A angles",
        ),
        DecisionRecord(
            decision="did not generate protein hydrogens",
            actor="hydrogen evidence preparation",
            rationale="protein-donor hydrogen bonds remain unknown",
        ),
        DecisionRecord(
            decision="did not convert hydrogen-bond geometry into affinity",
            actor="hydrogen evidence preparation",
            rationale="a geometric hydrogen bond is not a binding free energy",
        ),
    )

    provenance = ProvenanceRecord(
        records=(
            ProvenanceEntry(
                "experimental",
                "deposited heavy atoms",
                "Hydrogen parents retain deposited coordinates.",
                actor="RCSB PDB deposition",
            ),
            ProvenanceEntry(
                "calculated",
                "generated ligand hydrogens",
                (
                    f"{len(ligand_hydrogens)} hydrogens from the selected ligand state; "
                    "positions are computationally generated."
                ),
                actor="hydrogen evidence preparation",
            ),
            ProvenanceEntry(
                "predicted",
                "none",
                "No affinity or free-energy estimate is produced.",
                actor="hydrogen evidence preparation",
            ),
            ProvenanceEntry(
                "hypothesised",
                "hydrogen-bond evidence",
                "D-H...A geometry is a structural criterion, not an energy.",
                actor="hydrogen evidence preparation",
            ),
        )
    )

    return HydrogenEvidenceRecord(
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
        protein_preparation_record_path=protein_preparation_record_path,
        protein_preparation_record_sha256=protein_preparation_record_sha256,
        accession=accession,
        protein_selection=protein_selection,
        ligand_selection=ligand_selection,
        ligand_state_id=ligand_state_id,
        protein_state_id=protein_state_id,
        criteria=criteria,
        ligand_hydrogens=ligand_hydrogens,
        protein_hydrogen_status=protein_hydrogen_status,
        bonds=bonds,
        decisions=decisions,
        provenance=provenance,
        validation=ValidationRecord(checks=checks),
        discrepancies=tuple(discrepancies),
    )


def write_hydrogen_evidence(path: str | Path, record: HydrogenEvidenceRecord) -> Path:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(record.to_dict(), indent=2, sort_keys=False) + "\n",
        encoding="utf-8",
    )
    return destination


def load_hydrogen_evidence(path: str | Path) -> HydrogenEvidenceRecord:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    return HydrogenEvidenceRecord.from_dict(data)


def validate_hydrogen_evidence(record: HydrogenEvidenceRecord) -> list[str]:
    errors: list[str] = []
    if record.schema_version != HYDROGEN_EVIDENCE_SCHEMA_VERSION:
        errors.append(
            f"schema_version must be {HYDROGEN_EVIDENCE_SCHEMA_VERSION!r}, got "
            f"{record.schema_version!r}."
        )
    for name in (
        "source_manifest_path",
        "ligand_identity_record_path",
        "ligand_coordinate_record_path",
        "binding_site_record_path",
        "chemical_preparation_record_path",
        "protein_preparation_record_path",
    ):
        if not getattr(record, name):
            errors.append(f"{name} is required.")
    for name in (
        "source_manifest_sha256",
        "ligand_identity_record_sha256",
        "ligand_coordinate_record_sha256",
        "binding_site_record_sha256",
        "chemical_preparation_record_sha256",
        "protein_preparation_record_sha256",
    ):
        if not _SHA256_PATTERN.fullmatch(getattr(record, name)):
            errors.append(f"{name} must be a lowercase SHA-256 digest.")
    if not record.accession:
        errors.append("accession is required.")
    if record.protein_hydrogen_status not in ("generated", "not_run", "unavailable"):
        errors.append("protein_hydrogen_status is invalid.")

    hydrogen_ids = [item.hydrogen_id for item in record.ligand_hydrogens]
    if len(set(hydrogen_ids)) != len(hydrogen_ids):
        errors.append("hydrogen_id values must be unique.")
    for index, hydrogen in enumerate(record.ligand_hydrogens):
        if hydrogen.coordinate_source != "generated":
            errors.append(f"ligand_hydrogens[{index}] must be generated.")
        if hydrogen.parent_coordinate_source != "deposited":
            errors.append(f"ligand_hydrogens[{index}] parent must be deposited.")
        if hydrogen.donor_hydrogen_distance_angstrom <= 0:
            errors.append(f"ligand_hydrogens[{index}] D-H distance must be positive.")
        if len(hydrogen.position_angstrom) != 3 or not all(
            math.isfinite(value) for value in hydrogen.position_angstrom
        ):
            errors.append(f"ligand_hydrogens[{index}] position must be finite.")
    for index, bond in enumerate(record.bonds):
        if bond.status not in BOND_STATUSES:
            errors.append(f"bonds[{index}].status is invalid.")
        if bond.donor_acceptor_distance_angstrom <= 0:
            errors.append(f"bonds[{index}] D-A distance must be positive.")
        if not 0 <= bond.angle_degrees <= 180:
            errors.append(f"bonds[{index}] angle must be in [0, 180].")

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
