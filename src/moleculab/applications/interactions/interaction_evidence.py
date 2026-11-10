# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Evidence-aware fingerprints built from interaction interpretation records.

The application adapter keeps category-level interaction status separate from
pair-level observations.  It records whether a category was classified, assessed
and absent, unknown, or not applicable, and maps that status into the analysis
layer's present/absent/unknown fingerprint without treating unknown as zero.
"""

from __future__ import annotations

import json
import re
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from ...analysis.evidence import (
    EvidenceFeature,
    EvidenceFingerprint,
    EvidenceObservation,
    EvidenceSchema,
    build_evidence_fingerprint,
    criteria_fingerprint,
)
from .interaction_interpretation import (
    INTERACTION_STATUSES,
    InteractionInterpretationRecord,
)
from ..structure.structure_manifest import (
    PROVENANCE_LABELS,
    DecisionRecord,
    ProvenanceEntry,
    ProvenanceRecord,
    ValidationCheck,
    ValidationRecord,
)

INTERACTION_EVIDENCE_SCHEMA_VERSION = "1.0"
INTERACTION_TYPES = (
    "hydrogen_bond",
    "ionic",
    "hydrophobic",
    "aromatic",
    "water_mediated",
    "metal_coordination",
)
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


def _required_int(data: Mapping[str, Any], key: str) -> int:
    value = data.get(key)
    if value is None:
        raise ValueError(f"{key} is required.")
    return int(value)


@dataclass(frozen=True)
class EvidenceDiscrepancy:
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
    def from_dict(cls, data: Mapping[str, Any]) -> EvidenceDiscrepancy:
        data = _as_mapping(data, "evidence discrepancy")
        return cls(
            category=_required_string(data, "category"),
            status=_required_string(data, "status"),
            detail=_required_string(data, "detail"),
            observed=data.get("observed"),
            expected=data.get("expected"),
            resolution=data.get("resolution"),
        )


@dataclass(frozen=True)
class InteractionEvidenceRecord:
    source_manifest_path: str
    source_manifest_sha256: str
    ligand_identity_record_path: str
    ligand_identity_record_sha256: str
    binding_site_record_path: str
    binding_site_record_sha256: str
    interaction_record_path: str
    interaction_record_sha256: str
    accession: str
    schema: EvidenceSchema
    fingerprint: EvidenceFingerprint
    category_statuses: tuple[tuple[str, str], ...]
    decisions: tuple[DecisionRecord, ...]
    provenance: ProvenanceRecord
    validation: ValidationRecord
    discrepancies: tuple[EvidenceDiscrepancy, ...]
    schema_version: str = INTERACTION_EVIDENCE_SCHEMA_VERSION

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "source_manifest_path": self.source_manifest_path,
            "source_manifest_sha256": self.source_manifest_sha256,
            "ligand_identity_record_path": self.ligand_identity_record_path,
            "ligand_identity_record_sha256": self.ligand_identity_record_sha256,
            "binding_site_record_path": self.binding_site_record_path,
            "binding_site_record_sha256": self.binding_site_record_sha256,
            "interaction_record_path": self.interaction_record_path,
            "interaction_record_sha256": self.interaction_record_sha256,
            "accession": self.accession,
            "schema": self.schema.to_dict(),
            "fingerprint": self.fingerprint.to_dict(),
            "category_statuses": [
                {"interaction_type": name, "status": status}
                for name, status in self.category_statuses
            ],
            "decisions": [decision.to_dict() for decision in self.decisions],
            "provenance": self.provenance.to_dict(),
            "validation": self.validation.to_dict(),
            "discrepancies": [item.to_dict() for item in self.discrepancies],
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> InteractionEvidenceRecord:
        data = _as_mapping(data, "interaction evidence record")
        category_statuses = []
        for item in _as_list(data.get("category_statuses"), "category_statuses"):
            item = _as_mapping(item, "category status")
            category_statuses.append(
                (
                    _required_string(item, "interaction_type"),
                    _required_string(item, "status"),
                )
            )
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
            binding_site_record_path=_required_string(data, "binding_site_record_path"),
            binding_site_record_sha256=_required_string(
                data, "binding_site_record_sha256"
            ),
            interaction_record_path=_required_string(data, "interaction_record_path"),
            interaction_record_sha256=_required_string(
                data, "interaction_record_sha256"
            ),
            accession=_required_string(data, "accession"),
            schema=EvidenceSchema.from_dict(_as_mapping(data.get("schema"), "schema")),
            fingerprint=EvidenceFingerprint.from_dict(
                _as_mapping(data.get("fingerprint"), "fingerprint")
            ),
            category_statuses=tuple(category_statuses),
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
                EvidenceDiscrepancy.from_dict(
                    _as_mapping(item, "evidence discrepancy")
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


def build_interaction_evidence(
    *,
    interaction_record: InteractionInterpretationRecord,
    source_manifest_path: str,
    source_manifest_sha256: str,
    ligand_identity_record_path: str,
    ligand_identity_record_sha256: str,
    binding_site_record_path: str,
    binding_site_record_sha256: str,
    interaction_record_path: str,
    interaction_record_sha256: str,
) -> InteractionEvidenceRecord:
    """Summarize category-level interaction evidence as a three-state fingerprint."""

    by_type: dict[str, list] = {name: [] for name in INTERACTION_TYPES}
    for item in interaction_record.interactions:
        by_type.setdefault(item.interaction_type, []).append(item)

    observations: list[EvidenceObservation] = []
    category_statuses: list[tuple[str, str]] = []
    for interaction_type in INTERACTION_TYPES:
        items = by_type.get(interaction_type, [])
        applicable = True
        if not items:
            status = "unknown"
            state = "unknown"
            assessable = False
            evidence = "no interpretation record was supplied for this category"
            criteria = "no category criteria recorded"
        else:
            statuses = {item.status for item in items}
            if "classified" in statuses:
                status = "classified"
                state = "present"
                assessable = True
                classified = [item for item in items if item.status == "classified"]
                evidence = f"{len(classified)} classified interaction record(s)"
            elif "not_classified" in statuses:
                status = "not_classified"
                state = "absent"
                assessable = True
                evidence = "category assessed; no pair satisfied the explicit criteria"
            elif "unknown" in statuses:
                status = "unknown"
                state = "unknown"
                assessable = False
                evidence = items[0].evidence
            elif statuses == {"not_applicable"}:
                status = "not_applicable"
                state = "unknown"
                assessable = False
                applicable = False
                evidence = "not applicable: " + items[0].evidence
            else:
                status = "unknown"
                state = "unknown"
                assessable = False
                evidence = "mixed or unavailable category evidence"
            criteria = items[0].criteria
        category_statuses.append((interaction_type, status))
        feature = EvidenceFeature(scope="category", interaction_type=interaction_type)
        observations.append(
            EvidenceObservation(
                feature=feature,
                state=state,
                assessable=assessable,
                evidence=evidence,
                criteria=criteria,
                provenance="interaction interpretation record",
                applicable=applicable,
            )
        )

    criteria_payload = interaction_record.criteria.to_dict()
    criteria_fp = criteria_fingerprint(criteria_payload)
    criteria_json = json.dumps(criteria_payload, sort_keys=True, separators=(",", ":"))
    schema = EvidenceSchema.from_observations(
        observations,
        schema_id="interaction-category-v1",
        criteria_fingerprint=criteria_fp,
        criteria_json=criteria_json,
    )
    fingerprint = build_evidence_fingerprint(
        observations,
        schema_id="interaction-category-v1",
        criteria_fingerprint=criteria_fp,
        criteria_json=criteria_json,
    )
    if fingerprint.schema != schema:
        raise ValueError("fingerprint schema does not match constructed schema.")

    assessable_count = fingerprint.assessable_count
    checks = (
        _check(
            "provenance",
            "evidence.source_records",
            bool(
                source_manifest_path
                and ligand_identity_record_path
                and binding_site_record_path
                and interaction_record_path
            ),
            subject=interaction_record.accession,
            criterion="all prior records are recorded",
            observed="four record paths present",
            expected="manifest, identity, binding-site, and interaction records",
            detail="the fingerprint is traceable to the structural evidence chain",
        ),
        _check(
            "schema",
            "evidence.schema_compatible",
            True,
            subject=interaction_record.accession,
            criterion="category schema and criteria are explicit",
            observed=f"schema {schema.schema_id}",
            expected="explicit schema identifier and criteria fingerprint",
            detail=f"criteria fingerprint {criteria_fp}",
        ),
        _check(
            "evidence",
            "evidence.assessable",
            assessable_count > 0,
            subject=interaction_record.accession,
            criterion="at least one category is assessable",
            observed=f"{assessable_count} assessable categories",
            expected="one or more, or explicit unknown categories",
            detail="unknown categories are not counted as absent",
            not_run=assessable_count == 0,
        ),
        _check(
            "interpretation",
            "evidence.unknown_not_zero",
            True,
            subject=interaction_record.accession,
            criterion="unknown evidence is preserved rather than encoded as absence",
            observed=f"{fingerprint.unknown_count} unknown categories",
            expected="unknown distinct from absent",
            detail="unknown categories are excluded from assessable denominators",
        ),
        _check(
            "interpretation",
            "evidence.no_affinity",
            True,
            subject=interaction_record.accession,
            criterion="fingerprint is not converted into affinity",
            observed="no affinity or free-energy estimate is produced",
            expected="structural representation only",
            detail="a fingerprint is not a binding free energy",
        ),
    )

    discrepancies: list[EvidenceDiscrepancy] = []
    if assessable_count == 0:
        discrepancies.append(
            EvidenceDiscrepancy(
                "insufficient_evidence",
                "not_applicable",
                "No interaction category is assessable from the available evidence.",
                observed="0 assessable categories",
                expected="explicit chemical features and geometry",
                resolution="Unknown categories remain unknown; no zero was substituted.",
            )
        )

    decisions = (
        DecisionRecord(
            decision="mapped category statuses into present, absent, and unknown states",
            actor="interaction evidence fingerprint",
            rationale=(
                "classified categories are present; assessed categories with no "
                "qualifying pair are absent; missing evidence is unknown"
            ),
        ),
        DecisionRecord(
            decision="preserved category-level status separately from pair-level identity",
            actor="interaction evidence fingerprint",
            rationale=(
                "an unknown category does not imply that every candidate atom pair "
                "was assessed and found absent"
            ),
        ),
        DecisionRecord(
            decision="used assessable categories as the denominator",
            actor="interaction evidence fingerprint",
            rationale="unknown categories are excluded from occupancy denominators",
        ),
        DecisionRecord(
            decision="did not convert the fingerprint into affinity",
            actor="interaction evidence fingerprint",
            rationale="the fingerprint is a structural representation, not a free energy",
        ),
    )

    provenance = ProvenanceRecord(
        records=(
            ProvenanceEntry(
                "experimental",
                "interaction evidence",
                (
                    f"Fingerprint derived from {interaction_record.input_contact_count} "
                    "Batch 4 contacts and the Batch 5 interpretation record."
                ),
                actor="RCSB PDB deposition",
            ),
            ProvenanceEntry(
                "calculated",
                "three-state fingerprint",
                (
                    f"{fingerprint.present_count} present, "
                    f"{fingerprint.absent_count} absent, "
                    f"{fingerprint.unknown_count} unknown categories; "
                    f"{assessable_count} assessable."
                ),
                actor="interaction evidence fingerprint",
            ),
            ProvenanceEntry(
                "predicted",
                "none",
                "No affinity or free-energy prediction is produced.",
                actor="interaction evidence fingerprint",
            ),
            ProvenanceEntry(
                "hypothesised",
                "none",
                "No mechanism is inferred from the fingerprint.",
                actor="interaction evidence fingerprint",
            ),
        )
    )

    return InteractionEvidenceRecord(
        source_manifest_path=source_manifest_path,
        source_manifest_sha256=source_manifest_sha256,
        ligand_identity_record_path=ligand_identity_record_path,
        ligand_identity_record_sha256=ligand_identity_record_sha256,
        binding_site_record_path=binding_site_record_path,
        binding_site_record_sha256=binding_site_record_sha256,
        interaction_record_path=interaction_record_path,
        interaction_record_sha256=interaction_record_sha256,
        accession=interaction_record.accession,
        schema=schema,
        fingerprint=fingerprint,
        category_statuses=tuple(category_statuses),
        decisions=decisions,
        provenance=provenance,
        validation=ValidationRecord(checks=checks),
        discrepancies=tuple(discrepancies),
    )


def write_interaction_evidence(
    path: str | Path,
    record: InteractionEvidenceRecord,
) -> Path:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(record.to_dict(), indent=2, sort_keys=False) + "\n",
        encoding="utf-8",
    )
    return destination


def load_interaction_evidence(path: str | Path) -> InteractionEvidenceRecord:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    return InteractionEvidenceRecord.from_dict(data)


def validate_interaction_evidence(record: InteractionEvidenceRecord) -> list[str]:
    """Return schema and internal-consistency errors for *record*."""

    errors: list[str] = []
    if record.schema_version != INTERACTION_EVIDENCE_SCHEMA_VERSION:
        errors.append(
            f"schema_version must be {INTERACTION_EVIDENCE_SCHEMA_VERSION!r}, got "
            f"{record.schema_version!r}."
        )
    for name in (
        "source_manifest_path",
        "ligand_identity_record_path",
        "binding_site_record_path",
        "interaction_record_path",
    ):
        if not getattr(record, name):
            errors.append(f"{name} is required.")
    for name in (
        "source_manifest_sha256",
        "ligand_identity_record_sha256",
        "binding_site_record_sha256",
        "interaction_record_sha256",
    ):
        if not _SHA256_PATTERN.fullmatch(getattr(record, name)):
            errors.append(f"{name} must be a lowercase SHA-256 digest.")
    if not record.accession:
        errors.append("accession is required.")

    if record.schema.compatibility_key != record.fingerprint.schema.compatibility_key:
        errors.append("schema and fingerprint compatibility keys must match.")
    if record.schema != record.fingerprint.schema:
        errors.append("schema must equal the fingerprint schema.")

    known_types = {name for name, _ in record.category_statuses}
    if len(known_types) != len(record.category_statuses):
        errors.append("category_statuses must contain unique interaction types.")
    for name, status in record.category_statuses:
        if status not in INTERACTION_STATUSES:
            errors.append(f"category_statuses[{name}] has an invalid status.")
    if set(known_types) != set(INTERACTION_TYPES):
        errors.append("category_statuses must cover every interaction type.")

    if len(record.fingerprint.states) != len(record.fingerprint.schema.features):
        errors.append("fingerprint states must match schema features.")
    for feature, state, assessable in zip(
        record.fingerprint.schema.features,
        record.fingerprint.states,
        record.fingerprint.assessable,
        strict=True,
    ):
        if feature.scope != "category":
            errors.append("interaction evidence fingerprints require category features.")
        if state == "unknown" and assessable:
            errors.append("unknown features cannot be assessable.")
        if state in ("present", "absent") and not assessable:
            errors.append("present/absent features must be assessable.")

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
        if item.status not in ("open", "resolved", "not_applicable"):
            errors.append(f"discrepancies[{index}] has an invalid status.")
        if not item.category or not item.detail:
            errors.append(f"discrepancies[{index}] requires category and detail.")
        if item.status == "resolved" and not item.resolution:
            errors.append(f"discrepancies[{index}] is resolved but has no resolution.")

    return errors
