# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Sensitivity and robustness records for structural interaction conclusions.

This layer does not add interaction categories.  It re-evaluates existing
observations under explicit analysis assumptions and reports whether a
conclusion is invariant, assumption-sensitive, unsupported, or unknown.

Each axis is independent.  No composite confidence score, state ranking, or
affinity estimate is produced.
"""

from __future__ import annotations

import json
import math
import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from ...analysis.interactions import HydrogenBondCriteria
from ..chemistry.chemical_state_enumeration import ChemicalState
from .combination_interactions import (
    CATEGORY_ORDER,
    CombinationAssessment,
    assess_combination,
)
from .hydrogen_evidence import (
    HydrogenPosition,
    assess_ligand_donor_hydrogen_bonds,
)
from .interaction_interpretation import InteractionCriteria
from ..structure.ligand_coordinates import CoordinateAtom
from ..chemistry.protein_preparation import PreparedProteinAtom, PreparedProteinState
from ..structure.structure_manifest import (
    PROVENANCE_LABELS,
    DecisionRecord,
    LigandSelection,
    ProteinSelection,
    ProvenanceEntry,
    ProvenanceRecord,
    ValidationCheck,
    ValidationRecord,
)
from .water_evidence import WaterAtomRecord, assess_water_bridges

ROBUSTNESS_SCHEMA_VERSION = "1.0"
CONCLUSION_STATUSES = (
    "invariant",
    "assumption_sensitive",
    "unsupported",
    "unknown",
)
PLACEMENT_STATUSES = (
    "invariant_supported",
    "placement_sensitive",
    "invariant_unsupported",
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


def _required_float(data: Mapping[str, Any], key: str) -> float:
    value = data.get(key)
    if value is None:
        raise ValueError(f"{key} is required.")
    return float(value)


def _required_int(data: Mapping[str, Any], key: str) -> int:
    value = data.get(key)
    if value is None:
        raise ValueError(f"{key} is required.")
    return int(value)


def _optional_string(data: Mapping[str, Any], key: str) -> str | None:
    value = data.get(key)
    return None if value is None else str(value)


def _string_tuple(data: Mapping[str, Any], key: str) -> tuple[str, ...]:
    return tuple(str(item) for item in _as_list(data.get(key), key))


def _float_tuple(data: Mapping[str, Any], key: str) -> tuple[float, ...]:
    return tuple(float(item) for item in _as_list(data.get(key), key))


@dataclass(frozen=True)
class ClassificationCount:
    interaction_type: str
    count: int

    def to_dict(self) -> dict[str, Any]:
        return {
            "interaction_type": self.interaction_type,
            "count": self.count,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> ClassificationCount:
        data = _as_mapping(data, "classification count")
        return cls(
            interaction_type=_required_string(data, "interaction_type"),
            count=_required_int(data, "count"),
        )


@dataclass(frozen=True)
class CategoryStatusCount:
    interaction_type: str
    status: str
    combination_count: int

    def to_dict(self) -> dict[str, Any]:
        return {
            "interaction_type": self.interaction_type,
            "status": self.status,
            "combination_count": self.combination_count,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> CategoryStatusCount:
        data = _as_mapping(data, "category status count")
        return cls(
            interaction_type=_required_string(data, "interaction_type"),
            status=_required_string(data, "status"),
            combination_count=_required_int(data, "combination_count"),
        )


@dataclass(frozen=True)
class ContactCutoffResult:
    cutoff_angstrom: float
    ligand_contact_atom_count: int
    protein_contact_atom_count: int
    combination_count: int
    minimum_classified_per_combination: int
    maximum_classified_per_combination: int
    classified_counts: tuple[ClassificationCount, ...]
    category_status_counts: tuple[CategoryStatusCount, ...]

    def to_dict(self) -> dict[str, Any]:
        return {
            "cutoff_angstrom": self.cutoff_angstrom,
            "ligand_contact_atom_count": self.ligand_contact_atom_count,
            "protein_contact_atom_count": self.protein_contact_atom_count,
            "combination_count": self.combination_count,
            "minimum_classified_per_combination": (
                self.minimum_classified_per_combination
            ),
            "maximum_classified_per_combination": (
                self.maximum_classified_per_combination
            ),
            "classified_counts": [item.to_dict() for item in self.classified_counts],
            "category_status_counts": [
                item.to_dict() for item in self.category_status_counts
            ],
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> ContactCutoffResult:
        data = _as_mapping(data, "contact cutoff result")
        return cls(
            cutoff_angstrom=_required_float(data, "cutoff_angstrom"),
            ligand_contact_atom_count=_required_int(
                data, "ligand_contact_atom_count"
            ),
            protein_contact_atom_count=_required_int(
                data, "protein_contact_atom_count"
            ),
            combination_count=_required_int(data, "combination_count"),
            minimum_classified_per_combination=_required_int(
                data, "minimum_classified_per_combination"
            ),
            maximum_classified_per_combination=_required_int(
                data, "maximum_classified_per_combination"
            ),
            classified_counts=tuple(
                ClassificationCount.from_dict(
                    _as_mapping(item, "classification count")
                )
                for item in _as_list(data.get("classified_counts"), "classified_counts")
            ),
            category_status_counts=tuple(
                CategoryStatusCount.from_dict(
                    _as_mapping(item, "category status count")
                )
                for item in _as_list(
                    data.get("category_status_counts"), "category_status_counts"
                )
            ),
        )


@dataclass(frozen=True)
class HydrogenCandidateGeometry:
    hydrogen_id: str
    donor_atom_id: str
    acceptor_atom_id: str
    donor_acceptor_distance_angstrom: float
    angle_degrees: float
    default_status: str

    @property
    def candidate_id(self) -> str:
        return f"{self.donor_atom_id}->{self.acceptor_atom_id}"

    def to_dict(self) -> dict[str, Any]:
        return {
            "hydrogen_id": self.hydrogen_id,
            "donor_atom_id": self.donor_atom_id,
            "acceptor_atom_id": self.acceptor_atom_id,
            "donor_acceptor_distance_angstrom": self.donor_acceptor_distance_angstrom,
            "angle_degrees": self.angle_degrees,
            "default_status": self.default_status,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> HydrogenCandidateGeometry:
        data = _as_mapping(data, "hydrogen candidate geometry")
        return cls(
            hydrogen_id=_required_string(data, "hydrogen_id"),
            donor_atom_id=_required_string(data, "donor_atom_id"),
            acceptor_atom_id=_required_string(data, "acceptor_atom_id"),
            donor_acceptor_distance_angstrom=_required_float(
                data, "donor_acceptor_distance_angstrom"
            ),
            angle_degrees=_required_float(data, "angle_degrees"),
            default_status=_required_string(data, "default_status"),
        )


@dataclass(frozen=True)
class HydrogenAngleResult:
    min_angle_degrees: float
    supported_candidate_count: int
    supported_candidate_ids: tuple[str, ...]

    def to_dict(self) -> dict[str, Any]:
        return {
            "min_angle_degrees": self.min_angle_degrees,
            "supported_candidate_count": self.supported_candidate_count,
            "supported_candidate_ids": list(self.supported_candidate_ids),
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> HydrogenAngleResult:
        data = _as_mapping(data, "hydrogen angle result")
        return cls(
            min_angle_degrees=_required_float(data, "min_angle_degrees"),
            supported_candidate_count=_required_int(
                data, "supported_candidate_count"
            ),
            supported_candidate_ids=_string_tuple(data, "supported_candidate_ids"),
        )


@dataclass(frozen=True)
class HydrogenPlacementResult:
    hydrogen_id: str
    donor_atom_id: str
    acceptor_atom_id: str
    current_angle_degrees: float
    current_status: str
    tested_offsets_degrees: tuple[float, ...]
    supported_offsets_degrees: tuple[float, ...]
    placement_status: str
    modest_band_offsets_degrees: tuple[float, ...]
    modest_band_supported_offsets_degrees: tuple[float, ...]
    modest_band_status: str

    @property
    def candidate_id(self) -> str:
        return f"{self.donor_atom_id}->{self.acceptor_atom_id}"

    def to_dict(self) -> dict[str, Any]:
        return {
            "hydrogen_id": self.hydrogen_id,
            "donor_atom_id": self.donor_atom_id,
            "acceptor_atom_id": self.acceptor_atom_id,
            "current_angle_degrees": self.current_angle_degrees,
            "current_status": self.current_status,
            "tested_offsets_degrees": list(self.tested_offsets_degrees),
            "supported_offsets_degrees": list(self.supported_offsets_degrees),
            "placement_status": self.placement_status,
            "modest_band_offsets_degrees": list(
                self.modest_band_offsets_degrees
            ),
            "modest_band_supported_offsets_degrees": list(
                self.modest_band_supported_offsets_degrees
            ),
            "modest_band_status": self.modest_band_status,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> HydrogenPlacementResult:
        data = _as_mapping(data, "hydrogen placement result")
        return cls(
            hydrogen_id=_required_string(data, "hydrogen_id"),
            donor_atom_id=_required_string(data, "donor_atom_id"),
            acceptor_atom_id=_required_string(data, "acceptor_atom_id"),
            current_angle_degrees=_required_float(
                data, "current_angle_degrees"
            ),
            current_status=_required_string(data, "current_status"),
            tested_offsets_degrees=_float_tuple(
                data, "tested_offsets_degrees"
            ),
            supported_offsets_degrees=_float_tuple(
                data, "supported_offsets_degrees"
            ),
            placement_status=_required_string(data, "placement_status"),
            modest_band_offsets_degrees=_float_tuple(
                data, "modest_band_offsets_degrees"
            ),
            modest_band_supported_offsets_degrees=_float_tuple(
                data, "modest_band_supported_offsets_degrees"
            ),
            modest_band_status=_required_string(data, "modest_band_status"),
        )


@dataclass(frozen=True)
class WaterBridgeStatus:
    water_id: str
    status: str

    def to_dict(self) -> dict[str, Any]:
        return {"water_id": self.water_id, "status": self.status}

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> WaterBridgeStatus:
        data = _as_mapping(data, "water bridge status")
        return cls(
            water_id=_required_string(data, "water_id"),
            status=_required_string(data, "status"),
        )


@dataclass(frozen=True)
class WaterCutoffResult:
    cutoff_angstrom: float
    plausible_count: int
    not_supported_count: int
    unknown_count: int
    water_statuses: tuple[WaterBridgeStatus, ...]

    def to_dict(self) -> dict[str, Any]:
        return {
            "cutoff_angstrom": self.cutoff_angstrom,
            "plausible_count": self.plausible_count,
            "not_supported_count": self.not_supported_count,
            "unknown_count": self.unknown_count,
            "water_statuses": [item.to_dict() for item in self.water_statuses],
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> WaterCutoffResult:
        data = _as_mapping(data, "water cutoff result")
        return cls(
            cutoff_angstrom=_required_float(data, "cutoff_angstrom"),
            plausible_count=_required_int(data, "plausible_count"),
            not_supported_count=_required_int(data, "not_supported_count"),
            unknown_count=_required_int(data, "unknown_count"),
            water_statuses=tuple(
                WaterBridgeStatus.from_dict(
                    _as_mapping(item, "water bridge status")
                )
                for item in _as_list(data.get("water_statuses"), "water_statuses")
            ),
        )


@dataclass(frozen=True)
class ClassificationCountRange:
    interaction_type: str
    minimum_count: int
    maximum_count: int

    def to_dict(self) -> dict[str, Any]:
        return {
            "interaction_type": self.interaction_type,
            "minimum_count": self.minimum_count,
            "maximum_count": self.maximum_count,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> ClassificationCountRange:
        data = _as_mapping(data, "classification count range")
        return cls(
            interaction_type=_required_string(data, "interaction_type"),
            minimum_count=_required_int(data, "minimum_count"),
            maximum_count=_required_int(data, "maximum_count"),
        )


@dataclass(frozen=True)
class StatePattern:
    state_id: str
    pattern_id: str
    classified_count: int
    classification_counts: tuple[ClassificationCount, ...]

    def to_dict(self) -> dict[str, Any]:
        return {
            "state_id": self.state_id,
            "pattern_id": self.pattern_id,
            "classified_count": self.classified_count,
            "classification_counts": [
                item.to_dict() for item in self.classification_counts
            ],
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> StatePattern:
        data = _as_mapping(data, "state pattern")
        return cls(
            state_id=_required_string(data, "state_id"),
            pattern_id=_required_string(data, "pattern_id"),
            classified_count=_required_int(data, "classified_count"),
            classification_counts=tuple(
                ClassificationCount.from_dict(
                    _as_mapping(item, "classification count")
                )
                for item in _as_list(
                    data.get("classification_counts"), "classification_counts"
                )
            ),
        )


@dataclass(frozen=True)
class StateSensitivitySummary:
    reference_ligand_state_id: str
    reference_protein_state_id: str
    ligand_patterns: tuple[StatePattern, ...]
    protein_patterns: tuple[StatePattern, ...]
    unique_ligand_pattern_count: int
    unique_protein_pattern_count: int
    ligand_states_matching_reference_pattern: tuple[str, ...]
    ligand_states_differing_from_reference_pattern: tuple[str, ...]
    protein_states_matching_reference_pattern: tuple[str, ...]
    protein_states_differing_from_reference_pattern: tuple[str, ...]
    minimum_classified_count: int
    maximum_classified_count: int
    category_count_ranges: tuple[ClassificationCountRange, ...]

    def to_dict(self) -> dict[str, Any]:
        return {
            "reference_ligand_state_id": self.reference_ligand_state_id,
            "reference_protein_state_id": self.reference_protein_state_id,
            "ligand_patterns": [item.to_dict() for item in self.ligand_patterns],
            "protein_patterns": [item.to_dict() for item in self.protein_patterns],
            "unique_ligand_pattern_count": self.unique_ligand_pattern_count,
            "unique_protein_pattern_count": self.unique_protein_pattern_count,
            "ligand_states_matching_reference_pattern": list(
                self.ligand_states_matching_reference_pattern
            ),
            "ligand_states_differing_from_reference_pattern": list(
                self.ligand_states_differing_from_reference_pattern
            ),
            "protein_states_matching_reference_pattern": list(
                self.protein_states_matching_reference_pattern
            ),
            "protein_states_differing_from_reference_pattern": list(
                self.protein_states_differing_from_reference_pattern
            ),
            "minimum_classified_count": self.minimum_classified_count,
            "maximum_classified_count": self.maximum_classified_count,
            "category_count_ranges": [
                item.to_dict() for item in self.category_count_ranges
            ],
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> StateSensitivitySummary:
        data = _as_mapping(data, "state sensitivity summary")
        return cls(
            reference_ligand_state_id=_required_string(
                data, "reference_ligand_state_id"
            ),
            reference_protein_state_id=_required_string(
                data, "reference_protein_state_id"
            ),
            ligand_patterns=tuple(
                StatePattern.from_dict(_as_mapping(item, "state pattern"))
                for item in _as_list(data.get("ligand_patterns"), "ligand_patterns")
            ),
            protein_patterns=tuple(
                StatePattern.from_dict(_as_mapping(item, "state pattern"))
                for item in _as_list(data.get("protein_patterns"), "protein_patterns")
            ),
            unique_ligand_pattern_count=_required_int(
                data, "unique_ligand_pattern_count"
            ),
            unique_protein_pattern_count=_required_int(
                data, "unique_protein_pattern_count"
            ),
            ligand_states_matching_reference_pattern=_string_tuple(
                data, "ligand_states_matching_reference_pattern"
            ),
            ligand_states_differing_from_reference_pattern=_string_tuple(
                data, "ligand_states_differing_from_reference_pattern"
            ),
            protein_states_matching_reference_pattern=_string_tuple(
                data, "protein_states_matching_reference_pattern"
            ),
            protein_states_differing_from_reference_pattern=_string_tuple(
                data, "protein_states_differing_from_reference_pattern"
            ),
            minimum_classified_count=_required_int(
                data, "minimum_classified_count"
            ),
            maximum_classified_count=_required_int(
                data, "maximum_classified_count"
            ),
            category_count_ranges=tuple(
                ClassificationCountRange.from_dict(
                    _as_mapping(item, "classification count range")
                )
                for item in _as_list(
                    data.get("category_count_ranges"), "category_count_ranges"
                )
            ),
        )


@dataclass(frozen=True)
class RobustnessConclusion:
    conclusion_id: str
    axis: str
    statement: str
    status: str
    supporting_conditions: tuple[str, ...]
    contradicting_conditions: tuple[str, ...]
    limitations: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "conclusion_id": self.conclusion_id,
            "axis": self.axis,
            "statement": self.statement,
            "status": self.status,
            "supporting_conditions": list(self.supporting_conditions),
            "contradicting_conditions": list(self.contradicting_conditions),
            "limitations": self.limitations,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> RobustnessConclusion:
        data = _as_mapping(data, "robustness conclusion")
        return cls(
            conclusion_id=_required_string(data, "conclusion_id"),
            axis=_required_string(data, "axis"),
            statement=_required_string(data, "statement"),
            status=_required_string(data, "status"),
            supporting_conditions=_string_tuple(
                data, "supporting_conditions"
            ),
            contradicting_conditions=_string_tuple(
                data, "contradicting_conditions"
            ),
            limitations=_required_string(data, "limitations"),
        )


@dataclass(frozen=True)
class RobustnessDiscrepancy:
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
    def from_dict(cls, data: Mapping[str, Any]) -> RobustnessDiscrepancy:
        data = _as_mapping(data, "robustness discrepancy")
        return cls(
            category=_required_string(data, "category"),
            status=_required_string(data, "status"),
            detail=_required_string(data, "detail"),
            observed=_optional_string(data, "observed"),
            expected=_optional_string(data, "expected"),
            resolution=_optional_string(data, "resolution"),
        )


@dataclass(frozen=True)
class InteractionRobustnessRecord:
    source_manifest_path: str
    source_manifest_sha256: str
    ligand_identity_record_path: str
    ligand_identity_record_sha256: str
    ligand_coordinate_record_path: str
    ligand_coordinate_record_sha256: str
    binding_site_record_path: str
    binding_site_record_sha256: str
    chemical_state_record_path: str
    chemical_state_record_sha256: str
    protein_preparation_record_path: str
    protein_preparation_record_sha256: str
    hydrogen_evidence_record_path: str
    hydrogen_evidence_record_sha256: str
    water_evidence_record_path: str
    water_evidence_record_sha256: str
    accession: str
    protein_selection: ProteinSelection
    ligand_selection: LigandSelection
    baseline_contact_cutoff_angstrom: float
    baseline_hydrogen_criteria: HydrogenBondCriteria
    baseline_water_cutoff_angstrom: float
    protein_hydrogen_status: str
    water_roles_assigned: bool
    contact_cutoffs_angstrom: tuple[float, ...]
    hydrogen_angle_thresholds_degrees: tuple[float, ...]
    hydrogen_placement_offsets_degrees: tuple[float, ...]
    water_cutoffs_angstrom: tuple[float, ...]
    contact_cutoff_results: tuple[ContactCutoffResult, ...]
    hydrogen_candidates: tuple[HydrogenCandidateGeometry, ...]
    hydrogen_angle_results: tuple[HydrogenAngleResult, ...]
    hydrogen_placement_results: tuple[HydrogenPlacementResult, ...]
    water_cutoff_results: tuple[WaterCutoffResult, ...]
    state_sensitivity: StateSensitivitySummary
    conclusions: tuple[RobustnessConclusion, ...]
    decisions: tuple[DecisionRecord, ...]
    provenance: ProvenanceRecord
    validation: ValidationRecord
    discrepancies: tuple[RobustnessDiscrepancy, ...]
    schema_version: str = ROBUSTNESS_SCHEMA_VERSION

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
            "chemical_state_record_path": self.chemical_state_record_path,
            "chemical_state_record_sha256": self.chemical_state_record_sha256,
            "protein_preparation_record_path": self.protein_preparation_record_path,
            "protein_preparation_record_sha256": (
                self.protein_preparation_record_sha256
            ),
            "hydrogen_evidence_record_path": self.hydrogen_evidence_record_path,
            "hydrogen_evidence_record_sha256": self.hydrogen_evidence_record_sha256,
            "water_evidence_record_path": self.water_evidence_record_path,
            "water_evidence_record_sha256": self.water_evidence_record_sha256,
            "accession": self.accession,
            "protein_selection": self.protein_selection.to_dict(),
            "ligand_selection": self.ligand_selection.to_dict(),
            "baseline_contact_cutoff_angstrom": (
                self.baseline_contact_cutoff_angstrom
            ),
            "baseline_hydrogen_criteria": {
                "max_distance_nm": self.baseline_hydrogen_criteria.max_distance,
                "min_angle_degrees": (
                    self.baseline_hydrogen_criteria.min_angle_degrees
                ),
            },
            "baseline_water_cutoff_angstrom": self.baseline_water_cutoff_angstrom,
            "protein_hydrogen_status": self.protein_hydrogen_status,
            "water_roles_assigned": self.water_roles_assigned,
            "contact_cutoffs_angstrom": list(self.contact_cutoffs_angstrom),
            "hydrogen_angle_thresholds_degrees": list(
                self.hydrogen_angle_thresholds_degrees
            ),
            "hydrogen_placement_offsets_degrees": list(
                self.hydrogen_placement_offsets_degrees
            ),
            "water_cutoffs_angstrom": list(self.water_cutoffs_angstrom),
            "contact_cutoff_results": [
                item.to_dict() for item in self.contact_cutoff_results
            ],
            "hydrogen_candidates": [
                item.to_dict() for item in self.hydrogen_candidates
            ],
            "hydrogen_angle_results": [
                item.to_dict() for item in self.hydrogen_angle_results
            ],
            "hydrogen_placement_results": [
                item.to_dict() for item in self.hydrogen_placement_results
            ],
            "water_cutoff_results": [
                item.to_dict() for item in self.water_cutoff_results
            ],
            "state_sensitivity": self.state_sensitivity.to_dict(),
            "conclusions": [item.to_dict() for item in self.conclusions],
            "decisions": [decision.to_dict() for decision in self.decisions],
            "provenance": self.provenance.to_dict(),
            "validation": self.validation.to_dict(),
            "discrepancies": [item.to_dict() for item in self.discrepancies],
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> InteractionRobustnessRecord:
        data = _as_mapping(data, "interaction robustness record")
        criteria = _as_mapping(
            data.get("baseline_hydrogen_criteria"),
            "baseline_hydrogen_criteria",
        )
        return cls(
            schema_version=_required_string(data, "schema_version"),
            source_manifest_path=_required_string(data, "source_manifest_path"),
            source_manifest_sha256=_required_string(
                data, "source_manifest_sha256"
            ),
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
            binding_site_record_path=_required_string(
                data, "binding_site_record_path"
            ),
            binding_site_record_sha256=_required_string(
                data, "binding_site_record_sha256"
            ),
            chemical_state_record_path=_required_string(
                data, "chemical_state_record_path"
            ),
            chemical_state_record_sha256=_required_string(
                data, "chemical_state_record_sha256"
            ),
            protein_preparation_record_path=_required_string(
                data, "protein_preparation_record_path"
            ),
            protein_preparation_record_sha256=_required_string(
                data, "protein_preparation_record_sha256"
            ),
            hydrogen_evidence_record_path=_required_string(
                data, "hydrogen_evidence_record_path"
            ),
            hydrogen_evidence_record_sha256=_required_string(
                data, "hydrogen_evidence_record_sha256"
            ),
            water_evidence_record_path=_required_string(
                data, "water_evidence_record_path"
            ),
            water_evidence_record_sha256=_required_string(
                data, "water_evidence_record_sha256"
            ),
            accession=_required_string(data, "accession"),
            protein_selection=ProteinSelection.from_dict(
                _as_mapping(data.get("protein_selection"), "protein_selection")
            ),
            ligand_selection=LigandSelection.from_dict(
                _as_mapping(data.get("ligand_selection"), "ligand_selection")
            ),
            baseline_contact_cutoff_angstrom=_required_float(
                data, "baseline_contact_cutoff_angstrom"
            ),
            baseline_hydrogen_criteria=HydrogenBondCriteria(
                max_distance=_required_float(criteria, "max_distance_nm"),
                min_angle_degrees=_required_float(
                    criteria, "min_angle_degrees"
                ),
            ),
            baseline_water_cutoff_angstrom=_required_float(
                data, "baseline_water_cutoff_angstrom"
            ),
            protein_hydrogen_status=_required_string(
                data, "protein_hydrogen_status"
            ),
            water_roles_assigned=bool(data.get("water_roles_assigned")),
            contact_cutoffs_angstrom=_float_tuple(
                data, "contact_cutoffs_angstrom"
            ),
            hydrogen_angle_thresholds_degrees=_float_tuple(
                data, "hydrogen_angle_thresholds_degrees"
            ),
            hydrogen_placement_offsets_degrees=_float_tuple(
                data, "hydrogen_placement_offsets_degrees"
            ),
            water_cutoffs_angstrom=_float_tuple(data, "water_cutoffs_angstrom"),
            contact_cutoff_results=tuple(
                ContactCutoffResult.from_dict(
                    _as_mapping(item, "contact cutoff result")
                )
                for item in _as_list(
                    data.get("contact_cutoff_results"),
                    "contact_cutoff_results",
                )
            ),
            hydrogen_candidates=tuple(
                HydrogenCandidateGeometry.from_dict(
                    _as_mapping(item, "hydrogen candidate geometry")
                )
                for item in _as_list(
                    data.get("hydrogen_candidates"), "hydrogen_candidates"
                )
            ),
            hydrogen_angle_results=tuple(
                HydrogenAngleResult.from_dict(
                    _as_mapping(item, "hydrogen angle result")
                )
                for item in _as_list(
                    data.get("hydrogen_angle_results"),
                    "hydrogen_angle_results",
                )
            ),
            hydrogen_placement_results=tuple(
                HydrogenPlacementResult.from_dict(
                    _as_mapping(item, "hydrogen placement result")
                )
                for item in _as_list(
                    data.get("hydrogen_placement_results"),
                    "hydrogen_placement_results",
                )
            ),
            water_cutoff_results=tuple(
                WaterCutoffResult.from_dict(
                    _as_mapping(item, "water cutoff result")
                )
                for item in _as_list(
                    data.get("water_cutoff_results"),
                    "water_cutoff_results",
                )
            ),
            state_sensitivity=StateSensitivitySummary.from_dict(
                _as_mapping(data.get("state_sensitivity"), "state_sensitivity")
            ),
            conclusions=tuple(
                RobustnessConclusion.from_dict(
                    _as_mapping(item, "robustness conclusion")
                )
                for item in _as_list(data.get("conclusions"), "conclusions")
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
                RobustnessDiscrepancy.from_dict(
                    _as_mapping(item, "robustness discrepancy")
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


def contact_atom_sets(
    ligand_atoms: Sequence[CoordinateAtom],
    protein_atoms: Sequence[PreparedProteinAtom],
    cutoff_angstrom: float,
) -> tuple[set[str], set[tuple[str, str, str]]]:
    """Return ligand and protein atom identities within *cutoff_angstrom*."""

    if not math.isfinite(float(cutoff_angstrom)) or cutoff_angstrom <= 0:
        raise ValueError("cutoff_angstrom must be finite and positive.")
    ligand_contact_atoms: set[str] = set()
    protein_contact_atoms: set[tuple[str, str, str]] = set()
    for ligand_atom in ligand_atoms:
        ligand_position = np.asarray(ligand_atom.position_angstrom, dtype=float)
        for protein_atom in protein_atoms:
            protein_position = np.asarray(
                protein_atom.position_angstrom, dtype=float
            )
            if float(np.linalg.norm(ligand_position - protein_position)) <= cutoff_angstrom:
                ligand_contact_atoms.add(ligand_atom.atom_name)
                protein_contact_atoms.add(
                    (
                        protein_atom.auth_asym_id,
                        protein_atom.auth_seq_id,
                        protein_atom.atom_name,
                    )
                )
    return ligand_contact_atoms, protein_contact_atoms


def _classification_counts(
    combinations: Sequence[CombinationAssessment],
) -> tuple[ClassificationCount, ...]:
    counts = {name: 0 for name in CATEGORY_ORDER}
    for combination in combinations:
        for interaction in combination.classified_interactions:
            counts[interaction.interaction_type] = (
                counts.get(interaction.interaction_type, 0) + 1
            )
    return tuple(
        ClassificationCount(name, counts[name]) for name in CATEGORY_ORDER
    )


def _category_status_counts(
    combinations: Sequence[CombinationAssessment],
) -> tuple[CategoryStatusCount, ...]:
    records = []
    for interaction_type in CATEGORY_ORDER:
        statuses: dict[str, int] = {}
        for combination in combinations:
            status = dict(combination.category_statuses)[interaction_type]
            statuses[status] = statuses.get(status, 0) + 1
        for status in sorted(statuses):
            records.append(
                CategoryStatusCount(
                    interaction_type,
                    status,
                    statuses[status],
                )
            )
    return tuple(records)


def summarize_contact_cutoff(
    combinations: Sequence[CombinationAssessment],
    *,
    cutoff_angstrom: float,
    ligand_contact_atom_count: int,
    protein_contact_atom_count: int,
) -> ContactCutoffResult:
    if not combinations:
        raise ValueError("combinations must not be empty.")
    per_combination = [len(item.classified_interactions) for item in combinations]
    return ContactCutoffResult(
        cutoff_angstrom=float(cutoff_angstrom),
        ligand_contact_atom_count=ligand_contact_atom_count,
        protein_contact_atom_count=protein_contact_atom_count,
        combination_count=len(combinations),
        minimum_classified_per_combination=min(per_combination),
        maximum_classified_per_combination=max(per_combination),
        classified_counts=_classification_counts(combinations),
        category_status_counts=_category_status_counts(combinations),
    )


def assess_contact_cutoff_sensitivity(
    *,
    ligand_states: Sequence[ChemicalState],
    protein_states: Sequence[PreparedProteinState],
    reference_to_deposited: Mapping[int, str],
    ligand_atoms: Sequence[CoordinateAtom],
    water_features: Sequence[Any],
    criteria: InteractionCriteria,
    cutoffs_angstrom: Sequence[float],
) -> tuple[ContactCutoffResult, ...]:
    """Reassess every state pair at each contact-atom cutoff."""

    if not ligand_states or not protein_states:
        raise ValueError("ligand_states and protein_states must not be empty.")
    if not cutoffs_angstrom:
        raise ValueError("cutoffs_angstrom must not be empty.")
    reference_atoms = protein_states[0].atoms
    reference_identities = tuple(
        (
            atom.atom_id,
            atom.auth_asym_id,
            atom.auth_seq_id,
            atom.atom_name,
            atom.position_angstrom,
        )
        for atom in reference_atoms
    )
    for state in protein_states[1:]:
        identities = tuple(
            (
                atom.atom_id,
                atom.auth_asym_id,
                atom.auth_seq_id,
                atom.atom_name,
                atom.position_angstrom,
            )
            for atom in state.atoms
        )
        if identities != reference_identities:
            raise ValueError(
                "protein states must share atom identities and coordinates for "
                "a common contact-atom cutoff analysis."
            )

    results = []
    for cutoff in cutoffs_angstrom:
        ligand_contacts, protein_contacts = contact_atom_sets(
            ligand_atoms,
            protein_states[0].atoms,
            cutoff,
        )
        combinations = tuple(
            assess_combination(
                ligand_state=ligand_state,
                protein_state=protein_state,
                reference_to_deposited=reference_to_deposited,
                ligand_atoms=ligand_atoms,
                protein_atoms=protein_state.atoms,
                contact_ligand_atoms=ligand_contacts,
                contact_protein_atoms=protein_contacts,
                water_features=water_features,
                criteria=criteria,
            )
            for ligand_state in ligand_states
            for protein_state in protein_states
        )
        results.append(
            summarize_contact_cutoff(
                combinations,
                cutoff_angstrom=cutoff,
                ligand_contact_atom_count=len(ligand_contacts),
                protein_contact_atom_count=len(protein_contacts),
            )
        )
    return tuple(results)


def hydrogen_candidate_geometries(
    hydrogens: Sequence[HydrogenPosition],
    acceptors: Sequence[Any],
    criteria: HydrogenBondCriteria,
) -> tuple[HydrogenCandidateGeometry, ...]:
    records = []
    for hydrogen in hydrogens:
        for acceptor in acceptors:
            distance = _distance(
                hydrogen.parent_position_angstrom,
                acceptor.position_angstrom,
            )
            if distance > criteria.max_distance * 10.0:
                continue
            angle = _angle(
                hydrogen.parent_position_angstrom,
                hydrogen.position_angstrom,
                acceptor.position_angstrom,
            )
            records.append(
                HydrogenCandidateGeometry(
                    hydrogen_id=hydrogen.hydrogen_id,
                    donor_atom_id=hydrogen.parent_atom_id,
                    acceptor_atom_id=acceptor.atom_id,
                    donor_acceptor_distance_angstrom=distance,
                    angle_degrees=angle,
                    default_status=(
                        "supported"
                        if angle >= criteria.min_angle_degrees
                        else "not_supported"
                    ),
                )
            )
    return tuple(records)


def assess_hydrogen_angle_sensitivity(
    hydrogens: Sequence[HydrogenPosition],
    acceptors: Sequence[Any],
    criteria: HydrogenBondCriteria,
    thresholds_degrees: Sequence[float],
) -> tuple[HydrogenAngleResult, ...]:
    results = []
    for threshold in thresholds_degrees:
        bonds = assess_ligand_donor_hydrogen_bonds(
            hydrogens,
            acceptors,
            HydrogenBondCriteria(
                max_distance=criteria.max_distance,
                min_angle_degrees=float(threshold),
            ),
        )
        supported = tuple(
            f"{bond.donor_atom_id}->{bond.acceptor_atom_id}"
            for bond in bonds
            if bond.status == "supported"
        )
        results.append(
            HydrogenAngleResult(
                min_angle_degrees=float(threshold),
                supported_candidate_count=len(supported),
                supported_candidate_ids=supported,
            )
        )
    return tuple(results)


def _unit(vector: np.ndarray) -> np.ndarray:
    norm = float(np.linalg.norm(vector))
    if norm == 0.0:
        raise ValueError("cannot normalize a zero vector.")
    return vector / norm


def _hydrogen_placement_status(
    supported_offsets: Sequence[float],
    tested_offsets: Sequence[float],
) -> str:
    if all(offset in supported_offsets for offset in tested_offsets):
        return "invariant_supported"
    if not any(offset in supported_offsets for offset in tested_offsets):
        return "invariant_unsupported"
    return "placement_sensitive"


def assess_hydrogen_placement_sensitivity(
    hydrogens: Sequence[HydrogenPosition],
    acceptors: Sequence[Any],
    criteria: HydrogenBondCriteria,
    offsets_degrees: Sequence[float],
    *,
    modest_band_degrees: float = 10.0,
) -> tuple[HydrogenPlacementResult, ...]:
    """Perturb the D-H orientation in the D-H/D-A plane and reassess geometry.

    The perturbation is a sensitivity probe, not a new chemistry-toolkit
    preparation and not an experimental coordinate.
    """

    offsets = tuple(float(offset) for offset in offsets_degrees)
    if not offsets:
        raise ValueError("offsets_degrees must not be empty.")
    if 0.0 not in offsets:
        raise ValueError("offsets_degrees must include 0 for the reference placement.")
    records = []
    for hydrogen in hydrogens:
        donor = np.asarray(hydrogen.parent_position_angstrom, dtype=float)
        current_hydrogen = np.asarray(hydrogen.position_angstrom, dtype=float)
        for acceptor in acceptors:
            acceptor_position = np.asarray(
                acceptor.position_angstrom, dtype=float
            )
            donor_acceptor = acceptor_position - donor
            donor_acceptor_distance = float(np.linalg.norm(donor_acceptor))
            if donor_acceptor_distance > criteria.max_distance * 10.0:
                continue
            donor_hydrogen = current_hydrogen - donor
            bond_length = float(np.linalg.norm(donor_hydrogen))
            if bond_length == 0.0:
                continue
            acceptor_direction = _unit(donor_acceptor)
            hydrogen_direction = _unit(donor_hydrogen)
            cosine = float(
                np.clip(
                    np.dot(acceptor_direction, hydrogen_direction),
                    -1.0,
                    1.0,
                )
            )
            current_theta = math.acos(cosine)
            plane = hydrogen_direction - acceptor_direction * cosine
            if float(np.linalg.norm(plane)) < 1e-9:
                reference = (
                    np.array([1.0, 0.0, 0.0])
                    if abs(acceptor_direction[0]) < 0.9
                    else np.array([0.0, 1.0, 0.0])
                )
                plane = reference - acceptor_direction * np.dot(
                    reference, acceptor_direction
                )
            plane_direction = _unit(plane)

            supported_offsets = []
            for offset in offsets:
                target_theta = min(
                    max(current_theta + math.radians(float(offset)), 0.0),
                    math.pi,
                )
                alternative_direction = (
                    acceptor_direction * math.cos(target_theta)
                    + plane_direction * math.sin(target_theta)
                )
                alternative_hydrogen = donor + alternative_direction * bond_length
                angle = _angle(
                    donor,
                    alternative_hydrogen,
                    acceptor_position,
                )
                if angle >= criteria.min_angle_degrees:
                    supported_offsets.append(float(offset))

            current_angle = _angle(donor, current_hydrogen, acceptor_position)
            current_status = (
                "supported"
                if current_angle >= criteria.min_angle_degrees
                else "not_supported"
            )
            modest_offsets = tuple(
                offset for offset in offsets if abs(offset) <= modest_band_degrees
            )
            modest_supported = tuple(
                float(offset)
                for offset in modest_offsets
                if float(offset) in supported_offsets
            )
            records.append(
                HydrogenPlacementResult(
                    hydrogen_id=hydrogen.hydrogen_id,
                    donor_atom_id=hydrogen.parent_atom_id,
                    acceptor_atom_id=acceptor.atom_id,
                    current_angle_degrees=current_angle,
                    current_status=current_status,
                    tested_offsets_degrees=offsets,
                    supported_offsets_degrees=tuple(supported_offsets),
                    placement_status=_hydrogen_placement_status(
                        supported_offsets,
                        offsets,
                    ),
                    modest_band_offsets_degrees=modest_offsets,
                    modest_band_supported_offsets_degrees=modest_supported,
                    modest_band_status=_hydrogen_placement_status(
                        supported_offsets,
                        modest_offsets,
                    ),
                )
            )
    return tuple(records)


def assess_water_cutoff_sensitivity(
    waters: Sequence[WaterAtomRecord],
    ligand_polar_atoms: Sequence[Any],
    protein_polar_atoms: Sequence[Any],
    cutoffs_angstrom: Sequence[float],
) -> tuple[WaterCutoffResult, ...]:
    if not cutoffs_angstrom:
        raise ValueError("cutoffs_angstrom must not be empty.")
    results = []
    for cutoff in cutoffs_angstrom:
        bridges = assess_water_bridges(
            waters,
            ligand_polar_atoms,
            protein_polar_atoms,
            float(cutoff),
        )
        statuses = tuple(
            WaterBridgeStatus(bridge.water_id, bridge.status)
            for bridge in bridges
        )
        results.append(
            WaterCutoffResult(
                cutoff_angstrom=float(cutoff),
                plausible_count=sum(
                    1 for bridge in bridges if bridge.status == "geometrically_plausible"
                ),
                not_supported_count=sum(
                    1 for bridge in bridges if bridge.status == "not_supported"
                ),
                unknown_count=sum(
                    1 for bridge in bridges if bridge.status == "unknown"
                ),
                water_statuses=statuses,
            )
        )
    return tuple(results)


def _normalized_participants(
    interaction,
) -> tuple[str, ...]:
    participants = []
    for participant in interaction.participants:
        if "ligand:ligand:" in participant:
            participants.append(participant.rsplit(":", 1)[-1])
        else:
            participants.append(participant)
    return tuple(sorted(participants))


def _interaction_signature(
    combination: CombinationAssessment,
) -> frozenset[tuple[str, tuple[str, ...]]]:
    return frozenset(
        (
            interaction.interaction_type,
            _normalized_participants(interaction),
        )
        for interaction in combination.classified_interactions
    )


def _classification_counts_for(
    combination: CombinationAssessment,
) -> tuple[ClassificationCount, ...]:
    counts = {name: 0 for name in CATEGORY_ORDER}
    for interaction in combination.classified_interactions:
        counts[interaction.interaction_type] = (
            counts.get(interaction.interaction_type, 0) + 1
        )
    return tuple(
        ClassificationCount(name, counts[name])
        for name in CATEGORY_ORDER
        if counts.get(name, 0) > 0
    )


def summarize_state_sensitivity(
    combinations: Sequence[CombinationAssessment],
    *,
    reference_ligand_state_id: str = "reference-neutral",
    reference_protein_state_id: str = "standard-templates",
) -> StateSensitivitySummary:
    if not combinations:
        raise ValueError("combinations must not be empty.")
    by_pair = {
        (item.ligand_state_id, item.protein_state_id): item
        for item in combinations
    }
    ligand_state_ids = sorted({item.ligand_state_id for item in combinations})
    protein_state_ids = sorted({item.protein_state_id for item in combinations})
    if (reference_ligand_state_id, reference_protein_state_id) not in by_pair:
        raise ValueError(
            "the reference ligand/protein state pair is not present in the "
            "combinations."
        )

    signatures: dict[tuple[str, str], frozenset] = {}
    pattern_ids: dict[frozenset, str] = {}
    for combination in combinations:
        signature = _interaction_signature(combination)
        signatures[(combination.ligand_state_id, combination.protein_state_id)] = (
            signature
        )
        if signature not in pattern_ids:
            pattern_ids[signature] = f"pattern-{len(pattern_ids) + 1}"

    ligand_patterns = []
    for ligand_state_id in ligand_state_ids:
        combination = by_pair[(ligand_state_id, reference_protein_state_id)]
        signature = signatures[(ligand_state_id, reference_protein_state_id)]
        ligand_patterns.append(
            StatePattern(
                state_id=ligand_state_id,
                pattern_id=pattern_ids[signature],
                classified_count=len(combination.classified_interactions),
                classification_counts=_classification_counts_for(combination),
            )
        )

    protein_patterns = []
    for protein_state_id in protein_state_ids:
        combination = by_pair[(reference_ligand_state_id, protein_state_id)]
        signature = signatures[(reference_ligand_state_id, protein_state_id)]
        protein_patterns.append(
            StatePattern(
                state_id=protein_state_id,
                pattern_id=pattern_ids[signature],
                classified_count=len(combination.classified_interactions),
                classification_counts=_classification_counts_for(combination),
            )
        )

    reference_ligand_pattern = next(
        item.pattern_id
        for item in ligand_patterns
        if item.state_id == reference_ligand_state_id
    )
    reference_protein_pattern = next(
        item.pattern_id
        for item in protein_patterns
        if item.state_id == reference_protein_state_id
    )

    counts_by_type = {name: [] for name in CATEGORY_ORDER}
    per_combination = []
    for combination in combinations:
        per_combination.append(len(combination.classified_interactions))
        classification_counts = dict(
            (item.interaction_type, item.count)
            for item in _classification_counts_for(combination)
        )
        for name in CATEGORY_ORDER:
            counts_by_type[name].append(classification_counts.get(name, 0))

    return StateSensitivitySummary(
        reference_ligand_state_id=reference_ligand_state_id,
        reference_protein_state_id=reference_protein_state_id,
        ligand_patterns=tuple(ligand_patterns),
        protein_patterns=tuple(protein_patterns),
        unique_ligand_pattern_count=len(
            {item.pattern_id for item in ligand_patterns}
        ),
        unique_protein_pattern_count=len(
            {item.pattern_id for item in protein_patterns}
        ),
        ligand_states_matching_reference_pattern=tuple(
            item.state_id
            for item in ligand_patterns
            if item.pattern_id == reference_ligand_pattern
        ),
        ligand_states_differing_from_reference_pattern=tuple(
            item.state_id
            for item in ligand_patterns
            if item.pattern_id != reference_ligand_pattern
        ),
        protein_states_matching_reference_pattern=tuple(
            item.state_id
            for item in protein_patterns
            if item.pattern_id == reference_protein_pattern
        ),
        protein_states_differing_from_reference_pattern=tuple(
            item.state_id
            for item in protein_patterns
            if item.pattern_id != reference_protein_pattern
        ),
        minimum_classified_count=min(per_combination),
        maximum_classified_count=max(per_combination),
        category_count_ranges=tuple(
            ClassificationCountRange(
                name,
                min(counts_by_type[name]),
                max(counts_by_type[name]),
            )
            for name in CATEGORY_ORDER
        ),
    )


def _condition_label(axis: str, value: float, unit: str) -> str:
    return f"{axis}={value:g}_{unit}"


def _split_conditions(
    conditions: Sequence[str],
    predicates: Sequence[bool],
) -> tuple[tuple[str, ...], tuple[str, ...]]:
    support = tuple(
        condition
        for condition, predicate in zip(conditions, predicates, strict=True)
        if predicate
    )
    contradict = tuple(
        condition
        for condition, predicate in zip(conditions, predicates, strict=True)
        if not predicate
    )
    return support, contradict


def _status_from_split(
    support: Sequence[str],
    contradict: Sequence[str],
) -> str:
    if support and not contradict:
        return "invariant"
    if support and contradict:
        return "assumption_sensitive"
    return "unsupported"


def _contact_counts_for(
    result: ContactCutoffResult,
    interaction_type: str,
) -> int:
    return next(
        (
            item.count
            for item in result.classified_counts
            if item.interaction_type == interaction_type
        ),
        0,
    )


def _classified_combination_count(
    result: ContactCutoffResult,
    interaction_type: str,
) -> int:
    return next(
        (
            item.combination_count
            for item in result.category_status_counts
            if item.interaction_type == interaction_type
            and item.status == "classified"
        ),
        0,
    )


def derive_robustness_conclusions(
    *,
    contact_results: Sequence[ContactCutoffResult],
    hydrogen_candidates: Sequence[HydrogenCandidateGeometry],
    hydrogen_angle_results: Sequence[HydrogenAngleResult],
    hydrogen_placement_results: Sequence[HydrogenPlacementResult],
    water_results: Sequence[WaterCutoffResult],
    state_sensitivity: StateSensitivitySummary,
    baseline_contact_cutoff_angstrom: float,
    baseline_hydrogen_criteria: HydrogenBondCriteria,
    baseline_water_cutoff_angstrom: float,
    protein_hydrogen_status: str,
    water_roles_assigned: bool,
) -> tuple[RobustnessConclusion, ...]:
    conclusions: list[RobustnessConclusion] = []

    contact_conditions = tuple(
        _condition_label("contact_cutoff", item.cutoff_angstrom, "A")
        for item in contact_results
    )
    baseline_contact = next(
        item
        for item in contact_results
        if item.cutoff_angstrom == baseline_contact_cutoff_angstrom
    )
    count_predicates = tuple(
        (
            item.minimum_classified_per_combination,
            item.maximum_classified_per_combination,
        )
        == (
            baseline_contact.minimum_classified_per_combination,
            baseline_contact.maximum_classified_per_combination,
        )
        for item in contact_results
    )
    support, contradict = _split_conditions(
        contact_conditions,
        count_predicates,
    )
    counts_text = ", ".join(
        f"{item.cutoff_angstrom:g} A: {item.minimum_classified_per_combination}-"
        f"{item.maximum_classified_per_combination}"
        for item in contact_results
    )
    conclusions.append(
        RobustnessConclusion(
            conclusion_id="contact_classified_count",
            axis="contact_cutoff",
            statement=(
                "The number of classified interactions per state combination "
                f"depends on the contact cutoff ({counts_text})."
            ),
            status=_status_from_split(support, contradict),
            supporting_conditions=support,
            contradicting_conditions=contradict,
            limitations=(
                "The contact cutoff selects candidate atoms; interaction-specific "
                "distance and angle criteria remain fixed."
            ),
        )
    )

    hydrophobic_predicates = tuple(
        _classified_combination_count(item, "hydrophobic") > 0
        for item in contact_results
    )
    support, contradict = _split_conditions(
        contact_conditions,
        hydrophobic_predicates,
    )
    conclusions.append(
        RobustnessConclusion(
            conclusion_id="contact_hydrophobic",
            axis="contact_cutoff",
            statement=(
                "Hydrophobic atom-pair classification is present: "
                + ", ".join(
                    f"{item.cutoff_angstrom:g} A: "
                    f"{_classified_combination_count(item, 'hydrophobic')}/"
                    f"{item.combination_count} combinations"
                    for item in contact_results
                )
                + "."
            ),
            status=_status_from_split(support, contradict),
            supporting_conditions=support,
            contradicting_conditions=contradict,
            limitations=(
                "Hydrophobic assignments are template features, not measured "
                "interaction energies."
            ),
        )
    )

    aromatic_predicates = tuple(
        _classified_combination_count(item, "aromatic") > 0
        for item in contact_results
    )
    support, contradict = _split_conditions(
        contact_conditions,
        aromatic_predicates,
    )
    conclusions.append(
        RobustnessConclusion(
            conclusion_id="contact_aromatic",
            axis="contact_cutoff",
            statement=(
                "An aromatic ring-pair classification is present at every "
                "tested contact cutoff."
            ),
            status=_status_from_split(support, contradict),
            supporting_conditions=support,
            contradicting_conditions=contradict,
            limitations=(
                "The aromatic classification is a geometric ring-centroid "
                "criterion under explicit aromatic templates."
            ),
        )
    )

    ionic_predicates = tuple(
        _classified_combination_count(item, "ionic") > 0
        for item in contact_results
    )
    support, contradict = _split_conditions(contact_conditions, ionic_predicates)
    conclusions.append(
        RobustnessConclusion(
            conclusion_id="contact_ionic",
            axis="contact_cutoff",
            statement="An ionic interaction is supported at the assessed contacts.",
            status=_status_from_split(support, contradict),
            supporting_conditions=support,
            contradicting_conditions=contradict,
            limitations=(
                "Absence of a classified ionic pair is not evidence that "
                "electrostatics does not contribute to binding."
            ),
        )
    )

    angle_conditions = tuple(
        _condition_label("hydrogen_angle", item.min_angle_degrees, "deg")
        for item in hydrogen_angle_results
    )
    baseline_angle_count = next(
        item.supported_candidate_count
        for item in hydrogen_angle_results
        if item.min_angle_degrees == baseline_hydrogen_criteria.min_angle_degrees
    )
    angle_predicates = tuple(
        item.supported_candidate_count == baseline_angle_count
        for item in hydrogen_angle_results
    )
    support, contradict = _split_conditions(angle_conditions, angle_predicates)
    angle_text = ", ".join(
        f"{item.min_angle_degrees:g} deg: {item.supported_candidate_count}"
        for item in hydrogen_angle_results
    )
    conclusions.append(
        RobustnessConclusion(
            conclusion_id="hydrogen_angle_count",
            axis="hydrogen_angle",
            statement=(
                "The number of supported ligand-donor hydrogen-bond candidates "
                f"depends on the angle criterion ({angle_text})."
            ),
            status=_status_from_split(support, contradict),
            supporting_conditions=support,
            contradicting_conditions=contradict,
            limitations=(
                "Candidates use generated ligand hydrogens; protein-donor "
                "hydrogen bonds remain unavailable."
            ),
        )
    )

    for candidate in hydrogen_candidates:
        supported_conditions = tuple(
            condition
            for condition, result in zip(
                angle_conditions, hydrogen_angle_results, strict=True
            )
            if candidate.candidate_id in result.supported_candidate_ids
        )
        unsupported_conditions = tuple(
            condition
            for condition in angle_conditions
            if condition not in supported_conditions
        )
        status = _status_from_split(
            supported_conditions,
            unsupported_conditions,
        )
        conclusions.append(
            RobustnessConclusion(
                conclusion_id=f"hydrogen_angle.{candidate.candidate_id}",
                axis="hydrogen_angle",
                statement=(
                    f"{candidate.candidate_id} is supported at "
                    f"{len(supported_conditions)} of {len(angle_conditions)} "
                    "tested angle thresholds "
                    f"(D-A {candidate.donor_acceptor_distance_angstrom:.2f} A, "
                    f"D-H...A {candidate.angle_degrees:.1f} deg)."
                ),
                status=status,
                supporting_conditions=supported_conditions,
                contradicting_conditions=unsupported_conditions,
                limitations=(
                    "The angular result depends on the generated hydrogen "
                    "position, not an experimentally observed hydrogen."
                ),
            )
        )

    for placement in hydrogen_placement_results:
        modest_offsets = placement.modest_band_offsets_degrees
        modest_conditions = tuple(
            _condition_label("placement_offset", offset, "deg")
            for offset in modest_offsets
        )
        supported_conditions = tuple(
            _condition_label("placement_offset", offset, "deg")
            for offset in placement.modest_band_supported_offsets_degrees
        )
        contradicting_conditions = tuple(
            condition
            for condition in modest_conditions
            if condition not in supported_conditions
        )
        conclusions.append(
            RobustnessConclusion(
                conclusion_id=f"hydrogen_placement.{placement.candidate_id}",
                axis="hydrogen_placement",
                statement=(
                    f"{placement.candidate_id} "
                    + (
                        "remains supported"
                        if placement.modest_band_status == "invariant_supported"
                        else (
                            "remains unsupported"
                            if placement.modest_band_status
                            == "invariant_unsupported"
                            else "changes classification"
                        )
                    )
                    + " under a modest D-H orientation perturbation of "
                    f"+/-{max(abs(offset) for offset in modest_offsets):g} deg."
                ),
                status=_status_from_split(
                    supported_conditions,
                    contradicting_conditions,
                ),
                supporting_conditions=supported_conditions,
                contradicting_conditions=contradicting_conditions,
                limitations=(
                    "The perturbation is a sensitivity probe, not an "
                    "alternative chemistry-toolkit preparation or a calibrated "
                    "uncertainty interval. It samples one D-H/D-A plane, not "
                    "the full orientation space. Outside the modest band, a "
                    "more aligned hypothetical placement can support a "
                    "candidate that is unsupported at the reference placement; "
                    "that is not evidence from the reference structure."
                ),
            )
        )

    water_conditions = tuple(
        _condition_label("water_cutoff", item.cutoff_angstrom, "A")
        for item in water_results
    )
    water_predicates = tuple(item.plausible_count > 0 for item in water_results)
    support, contradict = _split_conditions(water_conditions, water_predicates)
    water_text = ", ".join(
        f"{item.cutoff_angstrom:g} A: {item.plausible_count}"
        for item in water_results
    )
    conclusions.append(
        RobustnessConclusion(
            conclusion_id="water_bridge_geometry",
            axis="water_cutoff",
            statement=(
                "A geometrically plausible water bridge is present under "
                f"some cutoffs ({water_text})."
            ),
            status=_status_from_split(support, contradict),
            supporting_conditions=support,
            contradicting_conditions=contradict,
            limitations=(
                "Geometric plausibility is not a chemically classified "
                "water-mediated interaction."
            ),
        )
    )

    conclusions.append(
        RobustnessConclusion(
            conclusion_id="water_chemical_classification",
            axis="water_cutoff",
            statement=(
                "No water-mediated interaction is chemically classified under "
                "any tested water cutoff."
            ),
            status="unknown" if not water_roles_assigned else "assumption_sensitive",
            supporting_conditions=(),
            contradicting_conditions=water_conditions,
            limitations=(
                "Water orientation and donor/acceptor roles remain unassigned; "
                "a nearby water is not automatically a donor or acceptor."
            ),
        )
    )

    protein_support = tuple(
        f"protein_state={state_id}"
        for state_id in state_sensitivity.protein_states_matching_reference_pattern
    )
    protein_contradict = tuple(
        f"protein_state={state_id}"
        for state_id in state_sensitivity.protein_states_differing_from_reference_pattern
    )
    conclusions.append(
        RobustnessConclusion(
            conclusion_id="protein_state_pattern",
            axis="protein_state",
            statement=(
                "The classified interaction pattern is invariant across the "
                "tested protein histidine-protonation states."
                if not protein_contradict
                else (
                    "The classified interaction pattern changes across the "
                    "tested protein histidine-protonation states."
                )
            ),
            status=_status_from_split(protein_support, protein_contradict),
            supporting_conditions=protein_support,
            contradicting_conditions=protein_contradict,
            limitations=(
                "Only two explicit protein-template states are available; this "
                "does not establish invariance under other protein preparations."
            ),
        )
    )

    ligand_conditions = tuple(
        f"ligand_state={item.state_id}" for item in state_sensitivity.ligand_patterns
    )
    ligand_support = tuple(
        f"ligand_state={state_id}"
        for state_id in state_sensitivity.ligand_states_matching_reference_pattern
    )
    ligand_contradict = tuple(
        f"ligand_state={state_id}"
        for state_id in state_sensitivity.ligand_states_differing_from_reference_pattern
    )
    category_ranges_invariant = all(
        item.minimum_count == item.maximum_count
        for item in state_sensitivity.category_count_ranges
    )
    conclusions.append(
        RobustnessConclusion(
            conclusion_id="ligand_state_category_counts",
            axis="ligand_state",
            statement=(
                "Every ligand/protein state combination classifies the same "
                "number of interactions by category."
                if category_ranges_invariant
                else (
                    "The number of classified interactions by category changes "
                    "across ligand/protein state combinations."
                )
            ),
            status=(
                "invariant"
                if category_ranges_invariant
                else "assumption_sensitive"
            ),
            supporting_conditions=ligand_conditions,
            contradicting_conditions=(
                () if category_ranges_invariant else ligand_contradict
            ),
            limitations=(
                "Equal category counts do not make the specific participants "
                "identical or imply equal binding strength."
            ),
        )
    )
    conclusions.append(
        RobustnessConclusion(
            conclusion_id="ligand_state_specific_pairs",
            axis="ligand_state",
            statement=(
                "The specific classified hydrophobic and aromatic participants "
                "are invariant across ligand states."
                if not ligand_contradict
                else (
                    "The specific classified hydrophobic and aromatic "
                    "participants change across ligand states, even where "
                    "category counts are stable."
                )
            ),
            status=_status_from_split(ligand_support, ligand_contradict),
            supporting_conditions=ligand_support,
            contradicting_conditions=ligand_contradict,
            limitations=(
                "Ligand-state atom and ring identities are state-specific; "
                "the same feature count does not imply the same interaction set."
            ),
        )
    )

    if protein_hydrogen_status != "generated":
        conclusions.append(
            RobustnessConclusion(
                conclusion_id="protein_donor_hydrogen",
                axis="protein_state",
                statement=(
                    "Protein-donor hydrogen bonds remain unknown because "
                    "protein hydrogens were not prepared."
                ),
                status="unknown",
                supporting_conditions=(),
                contradicting_conditions=(),
                limitations=(
                    "The absence of a protein-donor candidate is not evidence "
                    "that no such hydrogen bond exists."
                ),
            )
        )

    return tuple(conclusions)


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


def build_interaction_robustness_record(
    *,
    source_manifest_path: str,
    source_manifest_sha256: str,
    ligand_identity_record_path: str,
    ligand_identity_record_sha256: str,
    ligand_coordinate_record_path: str,
    ligand_coordinate_record_sha256: str,
    binding_site_record_path: str,
    binding_site_record_sha256: str,
    chemical_state_record_path: str,
    chemical_state_record_sha256: str,
    protein_preparation_record_path: str,
    protein_preparation_record_sha256: str,
    hydrogen_evidence_record_path: str,
    hydrogen_evidence_record_sha256: str,
    water_evidence_record_path: str,
    water_evidence_record_sha256: str,
    accession: str,
    protein_selection: ProteinSelection,
    ligand_selection: LigandSelection,
    baseline_contact_cutoff_angstrom: float,
    baseline_hydrogen_criteria: HydrogenBondCriteria,
    baseline_water_cutoff_angstrom: float,
    protein_hydrogen_status: str,
    water_roles_assigned: bool,
    contact_cutoffs_angstrom: tuple[float, ...],
    hydrogen_angle_thresholds_degrees: tuple[float, ...],
    hydrogen_placement_offsets_degrees: tuple[float, ...],
    water_cutoffs_angstrom: tuple[float, ...],
    contact_cutoff_results: tuple[ContactCutoffResult, ...],
    hydrogen_candidates: tuple[HydrogenCandidateGeometry, ...],
    hydrogen_angle_results: tuple[HydrogenAngleResult, ...],
    hydrogen_placement_results: tuple[HydrogenPlacementResult, ...],
    water_cutoff_results: tuple[WaterCutoffResult, ...],
    state_sensitivity: StateSensitivitySummary,
) -> InteractionRobustnessRecord:
    conclusions = derive_robustness_conclusions(
        contact_results=contact_cutoff_results,
        hydrogen_candidates=hydrogen_candidates,
        hydrogen_angle_results=hydrogen_angle_results,
        hydrogen_placement_results=hydrogen_placement_results,
        water_results=water_cutoff_results,
        state_sensitivity=state_sensitivity,
        baseline_contact_cutoff_angstrom=baseline_contact_cutoff_angstrom,
        baseline_hydrogen_criteria=baseline_hydrogen_criteria,
        baseline_water_cutoff_angstrom=baseline_water_cutoff_angstrom,
        protein_hydrogen_status=protein_hydrogen_status,
        water_roles_assigned=water_roles_assigned,
    )

    checks = (
        _check(
            "provenance",
            "robustness.source_records",
            bool(
                source_manifest_path
                and ligand_identity_record_path
                and ligand_coordinate_record_path
                and binding_site_record_path
                and chemical_state_record_path
                and protein_preparation_record_path
                and hydrogen_evidence_record_path
                and water_evidence_record_path
            ),
            subject=accession,
            criterion="all structural and evidence records are referenced",
            observed="eight input records present",
            expected=(
                "manifest, identity, coordinates, binding site, ligand states, "
                "protein states, hydrogen evidence, water evidence"
            ),
            detail="the robustness record references prior records without modifying them",
        ),
        _check(
            "sensitivity",
            "robustness.contact_cutoffs",
            tuple(item.cutoff_angstrom for item in contact_cutoff_results)
            == contact_cutoffs_angstrom,
            subject=accession,
            criterion="contact results cover every tested cutoff",
            observed=f"{len(contact_cutoff_results)} results",
            expected=f"{len(contact_cutoffs_angstrom)} cutoffs",
            detail="each cutoff is an independent condition",
        ),
        _check(
            "sensitivity",
            "robustness.hydrogen_angles",
            tuple(item.min_angle_degrees for item in hydrogen_angle_results)
            == hydrogen_angle_thresholds_degrees,
            subject=accession,
            criterion="hydrogen-angle results cover every tested threshold",
            observed=f"{len(hydrogen_angle_results)} results",
            expected=f"{len(hydrogen_angle_thresholds_degrees)} thresholds",
            detail="each threshold is an independent condition",
        ),
        _check(
            "sensitivity",
            "robustness.hydrogen_placement",
            tuple(
                item.tested_offsets_degrees for item in hydrogen_placement_results
            )
            and all(
                item.tested_offsets_degrees
                == hydrogen_placement_offsets_degrees
                for item in hydrogen_placement_results
            ),
            subject=accession,
            criterion="hydrogen-placement results use the declared perturbation grid",
            observed=f"{len(hydrogen_placement_results)} candidates",
            expected="declared offset grid on every candidate",
            detail="perturbed positions are sensitivity probes, not new experimental coordinates",
        ),
        _check(
            "sensitivity",
            "robustness.water_cutoffs",
            tuple(item.cutoff_angstrom for item in water_cutoff_results)
            == water_cutoffs_angstrom,
            subject=accession,
            criterion="water results cover every tested cutoff",
            observed=f"{len(water_cutoff_results)} results",
            expected=f"{len(water_cutoffs_angstrom)} cutoffs",
            detail="each cutoff is an independent condition",
        ),
        _check(
            "interpretation",
            "robustness.no_composite_score",
            True,
            subject=accession,
            criterion="conclusions are not collapsed into a confidence score",
            observed="separate invariant, sensitive, unsupported, and unknown conclusions",
            expected="no aggregate score or ranking",
            detail="different uncertainties retain different physical meanings",
        ),
        _check(
            "interpretation",
            "robustness.no_affinity",
            True,
            subject=accession,
            criterion="robustness is not converted into affinity",
            observed="no affinity or free-energy estimate is produced",
            expected="assumption-sensitivity evidence only",
            detail="interaction counts are not binding free energies",
        ),
    )

    decisions = (
        DecisionRecord(
            decision="treated each assumption axis independently",
            actor="interaction robustness analysis",
            rationale=(
                "contact, angular, placement, state, and water uncertainties "
                "have different physical meanings"
            ),
        ),
        DecisionRecord(
            decision=(
                "preserved invariant, sensitive, unsupported, and unknown "
                "conclusions separately"
            ),
            actor="interaction robustness analysis",
            rationale=(
                "a robust conclusion is not the same as a high-confidence score"
            ),
        ),
        DecisionRecord(
            decision=(
                "perturbed generated hydrogen orientation only as a sensitivity "
                "probe"
            ),
            actor="interaction robustness analysis",
            rationale=(
                "the perturbations are not experimental coordinates or new "
                "toolkit preparations"
            ),
        ),
        DecisionRecord(
            decision="did not rank states or convert interactions into affinity",
            actor="interaction robustness analysis",
            rationale="sensitivity analysis is not a binding free-energy calculation",
        ),
    )

    provenance = ProvenanceRecord(
        records=(
            ProvenanceEntry(
                "experimental",
                "deposited heavy atoms and waters",
                "All sensitivity conditions reuse the original experimental coordinates.",
                actor="RCSB PDB deposition",
            ),
            ProvenanceEntry(
                "calculated",
                "interaction and geometric sensitivity",
                (
                    "Contact cutoffs, angular criteria, water cutoffs, and state "
                    "combinations are recalculated from prior records."
                ),
                actor="interaction robustness analysis",
            ),
            ProvenanceEntry(
                "predicted",
                "none",
                "No affinity or free-energy estimate is produced.",
                actor="interaction robustness analysis",
            ),
            ProvenanceEntry(
                "hypothesised",
                "hydrogen-placement perturbations",
                (
                    "Perturbed hydrogen orientations are sensitivity probes, "
                    "not experimental or toolkit-generated coordinates."
                ),
                actor="interaction robustness analysis",
            ),
        )
    )

    discrepancies = (
        RobustnessDiscrepancy(
            category="hydrogen_placement",
            status="not_applicable",
            detail="No experimentally observed hydrogen positions are available.",
            observed="generated and perturbed positions",
            expected="experimental hydrogen positions for direct angular evidence",
            resolution="Hydrogen-angle conclusions remain conditioned on generated positions.",
        ),
        RobustnessDiscrepancy(
            category="protein_hydrogens",
            status="not_applicable",
            detail="Protein hydrogens were not prepared.",
            observed=protein_hydrogen_status,
            expected="generated protein hydrogens for protein-donor hydrogen bonds",
            resolution="Protein-donor hydrogen bonds remain unknown.",
        ),
        RobustnessDiscrepancy(
            category="water_chemistry",
            status="not_applicable",
            detail="Water orientation and donor/acceptor roles remain unassigned.",
            observed="unassigned",
            expected="explicit water chemistry for a classified water-mediated interaction",
            resolution="Only geometric bridge plausibility is reported.",
        ),
    )

    return InteractionRobustnessRecord(
        source_manifest_path=source_manifest_path,
        source_manifest_sha256=source_manifest_sha256,
        ligand_identity_record_path=ligand_identity_record_path,
        ligand_identity_record_sha256=ligand_identity_record_sha256,
        ligand_coordinate_record_path=ligand_coordinate_record_path,
        ligand_coordinate_record_sha256=ligand_coordinate_record_sha256,
        binding_site_record_path=binding_site_record_path,
        binding_site_record_sha256=binding_site_record_sha256,
        chemical_state_record_path=chemical_state_record_path,
        chemical_state_record_sha256=chemical_state_record_sha256,
        protein_preparation_record_path=protein_preparation_record_path,
        protein_preparation_record_sha256=protein_preparation_record_sha256,
        hydrogen_evidence_record_path=hydrogen_evidence_record_path,
        hydrogen_evidence_record_sha256=hydrogen_evidence_record_sha256,
        water_evidence_record_path=water_evidence_record_path,
        water_evidence_record_sha256=water_evidence_record_sha256,
        accession=accession,
        protein_selection=protein_selection,
        ligand_selection=ligand_selection,
        baseline_contact_cutoff_angstrom=baseline_contact_cutoff_angstrom,
        baseline_hydrogen_criteria=baseline_hydrogen_criteria,
        baseline_water_cutoff_angstrom=baseline_water_cutoff_angstrom,
        protein_hydrogen_status=protein_hydrogen_status,
        water_roles_assigned=water_roles_assigned,
        contact_cutoffs_angstrom=contact_cutoffs_angstrom,
        hydrogen_angle_thresholds_degrees=hydrogen_angle_thresholds_degrees,
        hydrogen_placement_offsets_degrees=hydrogen_placement_offsets_degrees,
        water_cutoffs_angstrom=water_cutoffs_angstrom,
        contact_cutoff_results=contact_cutoff_results,
        hydrogen_candidates=hydrogen_candidates,
        hydrogen_angle_results=hydrogen_angle_results,
        hydrogen_placement_results=hydrogen_placement_results,
        water_cutoff_results=water_cutoff_results,
        state_sensitivity=state_sensitivity,
        conclusions=conclusions,
        decisions=decisions,
        provenance=provenance,
        validation=ValidationRecord(checks=checks),
        discrepancies=discrepancies,
    )


def write_interaction_robustness(
    path: str | Path,
    record: InteractionRobustnessRecord,
) -> Path:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(record.to_dict(), indent=2, sort_keys=False) + "\n",
        encoding="utf-8",
    )
    return destination


def load_interaction_robustness(
    path: str | Path,
) -> InteractionRobustnessRecord:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    return InteractionRobustnessRecord.from_dict(data)


def validate_interaction_robustness(
    record: InteractionRobustnessRecord,
) -> list[str]:
    errors: list[str] = []
    if record.schema_version != ROBUSTNESS_SCHEMA_VERSION:
        errors.append(
            f"schema_version must be {ROBUSTNESS_SCHEMA_VERSION!r}, got "
            f"{record.schema_version!r}."
        )
    for name in (
        "source_manifest_path",
        "ligand_identity_record_path",
        "ligand_coordinate_record_path",
        "binding_site_record_path",
        "chemical_state_record_path",
        "protein_preparation_record_path",
        "hydrogen_evidence_record_path",
        "water_evidence_record_path",
    ):
        if not getattr(record, name):
            errors.append(f"{name} is required.")
    for name in (
        "source_manifest_sha256",
        "ligand_identity_record_sha256",
        "ligand_coordinate_record_sha256",
        "binding_site_record_sha256",
        "chemical_state_record_sha256",
        "protein_preparation_record_sha256",
        "hydrogen_evidence_record_sha256",
        "water_evidence_record_sha256",
    ):
        if not _SHA256_PATTERN.fullmatch(getattr(record, name)):
            errors.append(f"{name} must be a lowercase SHA-256 digest.")
    if not record.accession:
        errors.append("accession is required.")
    if record.baseline_contact_cutoff_angstrom not in record.contact_cutoffs_angstrom:
        errors.append("baseline contact cutoff must be one of the tested cutoffs.")
    if (
        record.baseline_hydrogen_criteria.min_angle_degrees
        not in record.hydrogen_angle_thresholds_degrees
    ):
        errors.append("baseline hydrogen angle must be one of the tested thresholds.")
    if record.baseline_water_cutoff_angstrom not in record.water_cutoffs_angstrom:
        errors.append("baseline water cutoff must be one of the tested cutoffs.")
    if record.protein_hydrogen_status not in ("generated", "not_run", "unavailable"):
        errors.append("protein_hydrogen_status is invalid.")

    if tuple(item.cutoff_angstrom for item in record.contact_cutoff_results) != (
        record.contact_cutoffs_angstrom
    ):
        errors.append("contact_cutoff_results must cover the tested cutoffs.")
    if tuple(item.min_angle_degrees for item in record.hydrogen_angle_results) != (
        record.hydrogen_angle_thresholds_degrees
    ):
        errors.append("hydrogen_angle_results must cover the tested thresholds.")
    if tuple(item.cutoff_angstrom for item in record.water_cutoff_results) != (
        record.water_cutoffs_angstrom
    ):
        errors.append("water_cutoff_results must cover the tested cutoffs.")
    for index, placement in enumerate(record.hydrogen_placement_results):
        if placement.tested_offsets_degrees != (
            record.hydrogen_placement_offsets_degrees
        ):
            errors.append(
                f"hydrogen_placement_results[{index}] must use the declared offsets."
            )
        if placement.placement_status not in PLACEMENT_STATUSES:
            errors.append(
                f"hydrogen_placement_results[{index}].placement_status is invalid."
            )
        if placement.modest_band_status not in PLACEMENT_STATUSES:
            errors.append(
                f"hydrogen_placement_results[{index}].modest_band_status is invalid."
            )

    conclusion_ids = set()
    for index, conclusion in enumerate(record.conclusions):
        if conclusion.status not in CONCLUSION_STATUSES:
            errors.append(f"conclusions[{index}].status is invalid.")
        if not conclusion.conclusion_id or not conclusion.statement:
            errors.append(
                f"conclusions[{index}] requires conclusion_id and statement."
            )
        if conclusion.conclusion_id in conclusion_ids:
            errors.append(f"conclusions[{index}].conclusion_id is duplicated.")
        conclusion_ids.add(conclusion.conclusion_id)
        if conclusion.status == "invariant" and conclusion.contradicting_conditions:
            errors.append(
                f"conclusions[{index}] is invariant but has contradicting conditions."
            )
        if conclusion.status == "assumption_sensitive" and not (
            conclusion.supporting_conditions
            and conclusion.contradicting_conditions
        ):
            errors.append(
                f"conclusions[{index}] is assumption-sensitive but lacks both "
                "supporting and contradicting conditions."
            )
        if conclusion.status == "unsupported" and conclusion.supporting_conditions:
            errors.append(
                f"conclusions[{index}] is unsupported but has supporting conditions."
            )
        if conclusion.status == "unknown" and not conclusion.limitations:
            errors.append(
                f"conclusions[{index}] is unknown and requires a limitation."
            )

    if not record.decisions:
        errors.append("decisions must contain at least one record.")
    for index, decision in enumerate(record.decisions):
        if not decision.decision or not decision.actor or not decision.rationale:
            errors.append(
                f"decisions[{index}] requires decision, actor, and rationale."
            )
    if not record.provenance.records:
        errors.append("provenance.records must not be empty.")
    for index, entry in enumerate(record.provenance.records):
        if entry.label not in PROVENANCE_LABELS:
            errors.append(f"provenance.records[{index}] has an invalid label.")

    for index, item in enumerate(record.discrepancies):
        if item.status not in ("open", "resolved", "not_applicable"):
            errors.append(f"discrepancies[{index}] has an invalid status.")

    check_ids = set()
    for index, check in enumerate(record.validation.checks):
        if check.status not in ("pass", "fail", "not_run", "not_applicable"):
            errors.append(f"validation.checks[{index}] has an invalid status.")
        if not check.check_id or not check.category or not check.detail:
            errors.append(
                f"validation.checks[{index}] requires check_id, category, and detail."
            )
        if check.check_id in check_ids:
            errors.append(f"validation.checks[{index}].check_id is duplicated.")
        check_ids.add(check.check_id)

    return errors
