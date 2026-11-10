"""State-aware interpretation of the same experimental binding site.

This layer combines the explicit chemical-state enumeration with the
experimental binding-site contacts.  It maps each state's atoms back to the
Batch 2 reference indices and the deposited ligand atoms, then records
state-specific ligand feature summaries and three-state evidence fingerprints.

It does not rank states, estimate populations, infer protein chemistry, or
convert chemical eligibility into affinity.
"""

from __future__ import annotations

import json
import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from ..analysis.evidence import (
    EvidenceComparison,
    EvidenceFeature,
    EvidenceFingerprint,
    EvidenceObservation,
    EvidenceSchema,
    build_evidence_fingerprint,
    compare_evidence_fingerprints,
    criteria_fingerprint,
)
from .chemical_state_enumeration import ChemicalState, ChemicalStateEnumerationRecord
from .ligand_coordinates import CoordinateAtom
from .structure_manifest import (
    PROVENANCE_LABELS,
    DecisionRecord,
    LigandSelection,
    ProvenanceEntry,
    ProvenanceRecord,
    ValidationCheck,
    ValidationRecord,
)

STATE_INTERACTION_SCHEMA_VERSION = "1.0"
CATEGORY_ORDER = (
    "ligand_donor_acceptor",
    "ligand_ionic",
    "ligand_hydrophobic",
    "ligand_aromatic",
    "geometric_contact",
    "interaction_classification",
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


def _required_float(data: Mapping[str, Any], key: str) -> float:
    value = data.get(key)
    if value is None:
        raise ValueError(f"{key} is required.")
    return float(value)


def _string_tuple(data: Mapping[str, Any], key: str) -> tuple[str, ...]:
    return tuple(str(item) for item in _as_list(data.get(key), key))


@dataclass(frozen=True)
class StateInteractionCriteria:
    contact_cutoff_angstrom: float
    geometry_source: str
    feature_source: str
    classification_limitations: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "contact_cutoff_angstrom": self.contact_cutoff_angstrom,
            "geometry_source": self.geometry_source,
            "feature_source": self.feature_source,
            "classification_limitations": self.classification_limitations,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> StateInteractionCriteria:
        data = _as_mapping(data, "state interaction criteria")
        return cls(
            contact_cutoff_angstrom=_required_float(data, "contact_cutoff_angstrom"),
            geometry_source=_required_string(data, "geometry_source"),
            feature_source=_required_string(data, "feature_source"),
            classification_limitations=_required_string(
                data, "classification_limitations"
            ),
        )


@dataclass(frozen=True)
class FeatureSummary:
    atom_count: int
    donor_count: int
    acceptor_count: int
    hydrophobic_count: int
    aromatic_atom_count: int
    charged_atom_count: int
    formal_charge: int | None

    def to_dict(self) -> dict[str, Any]:
        return {
            "atom_count": self.atom_count,
            "donor_count": self.donor_count,
            "acceptor_count": self.acceptor_count,
            "hydrophobic_count": self.hydrophobic_count,
            "aromatic_atom_count": self.aromatic_atom_count,
            "charged_atom_count": self.charged_atom_count,
            "formal_charge": self.formal_charge,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> FeatureSummary:
        data = _as_mapping(data, "feature summary")
        return cls(
            atom_count=_required_int(data, "atom_count"),
            donor_count=_required_int(data, "donor_count"),
            acceptor_count=_required_int(data, "acceptor_count"),
            hydrophobic_count=_required_int(data, "hydrophobic_count"),
            aromatic_atom_count=_required_int(data, "aromatic_atom_count"),
            charged_atom_count=_required_int(data, "charged_atom_count"),
            formal_charge=(
                None if data.get("formal_charge") is None else int(data["formal_charge"])
            ),
        )


@dataclass(frozen=True)
class StateInteractionAssessment:
    state_id: str
    kind: str
    label: str
    formula: str
    proton_count: int
    formal_charge: int
    ligand_feature_summary: FeatureSummary
    contact_feature_summary: FeatureSummary
    contact_atom_count: int
    contact_atom_names: tuple[str, ...]
    donor_contact_atoms: tuple[str, ...]
    acceptor_contact_atoms: tuple[str, ...]
    hydrophobic_contact_atoms: tuple[str, ...]
    aromatic_contact_atoms: tuple[str, ...]
    charged_contact_atoms: tuple[str, ...]
    classified_interaction_participants: tuple[str, ...]
    fingerprint: EvidenceFingerprint

    def to_dict(self) -> dict[str, Any]:
        return {
            "state_id": self.state_id,
            "kind": self.kind,
            "label": self.label,
            "formula": self.formula,
            "proton_count": self.proton_count,
            "formal_charge": self.formal_charge,
            "ligand_feature_summary": self.ligand_feature_summary.to_dict(),
            "contact_feature_summary": self.contact_feature_summary.to_dict(),
            "contact_atom_count": self.contact_atom_count,
            "contact_atom_names": list(self.contact_atom_names),
            "donor_contact_atoms": list(self.donor_contact_atoms),
            "acceptor_contact_atoms": list(self.acceptor_contact_atoms),
            "hydrophobic_contact_atoms": list(self.hydrophobic_contact_atoms),
            "aromatic_contact_atoms": list(self.aromatic_contact_atoms),
            "charged_contact_atoms": list(self.charged_contact_atoms),
            "classified_interaction_participants": list(
                self.classified_interaction_participants
            ),
            "fingerprint": self.fingerprint.to_dict(),
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> StateInteractionAssessment:
        data = _as_mapping(data, "state interaction assessment")
        return cls(
            state_id=_required_string(data, "state_id"),
            kind=_required_string(data, "kind"),
            label=_required_string(data, "label"),
            formula=_required_string(data, "formula"),
            proton_count=_required_int(data, "proton_count"),
            formal_charge=_required_int(data, "formal_charge"),
            ligand_feature_summary=FeatureSummary.from_dict(
                _as_mapping(data.get("ligand_feature_summary"), "ligand_feature_summary")
            ),
            contact_feature_summary=FeatureSummary.from_dict(
                _as_mapping(data.get("contact_feature_summary"), "contact_feature_summary")
            ),
            contact_atom_count=_required_int(data, "contact_atom_count"),
            contact_atom_names=_string_tuple(data, "contact_atom_names"),
            donor_contact_atoms=_string_tuple(data, "donor_contact_atoms"),
            acceptor_contact_atoms=_string_tuple(data, "acceptor_contact_atoms"),
            hydrophobic_contact_atoms=_string_tuple(data, "hydrophobic_contact_atoms"),
            aromatic_contact_atoms=_string_tuple(data, "aromatic_contact_atoms"),
            charged_contact_atoms=_string_tuple(data, "charged_contact_atoms"),
            classified_interaction_participants=_string_tuple(
                data, "classified_interaction_participants"
            ),
            fingerprint=EvidenceFingerprint.from_dict(
                _as_mapping(data.get("fingerprint"), "fingerprint")
            ),
        )


@dataclass(frozen=True)
class StateInteractionComparison:
    reference_state_id: str
    state_id: str
    comparison: EvidenceComparison
    differing_categories: tuple[str, ...]
    ligand_feature_differences: tuple[tuple[str, int, int], ...]
    contact_feature_differences: tuple[tuple[str, int, int], ...]

    def to_dict(self) -> dict[str, Any]:
        return {
            "reference_state_id": self.reference_state_id,
            "state_id": self.state_id,
            "comparison": {
                "schema_id": self.comparison.schema_id,
                "criteria_fingerprint": self.comparison.criteria_fingerprint,
                "both_present": self.comparison.both_present,
                "both_absent": self.comparison.both_absent,
                "first_only_present": self.comparison.first_only_present,
                "second_only_present": self.comparison.second_only_present,
                "unknown_involved": self.comparison.unknown_involved,
                "not_applicable_involved": self.comparison.not_applicable_involved,
                "assessable_both": self.comparison.assessable_both,
            },
            "differing_categories": list(self.differing_categories),
            "ligand_feature_differences": [
                {"feature": feature, "reference_count": reference_count, "state_count": state_count}
                for feature, reference_count, state_count in self.ligand_feature_differences
            ],
            "contact_feature_differences": [
                {"feature": feature, "reference_count": reference_count, "state_count": state_count}
                for feature, reference_count, state_count in self.contact_feature_differences
            ],
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> StateInteractionComparison:
        data = _as_mapping(data, "state interaction comparison")
        comparison = _as_mapping(data.get("comparison"), "comparison")
        return cls(
            reference_state_id=_required_string(data, "reference_state_id"),
            state_id=_required_string(data, "state_id"),
            comparison=EvidenceComparison(
                schema_id=_required_string(comparison, "schema_id"),
                criteria_fingerprint=_required_string(
                    comparison, "criteria_fingerprint"
                ),
                both_present=_required_int(comparison, "both_present"),
                both_absent=_required_int(comparison, "both_absent"),
                first_only_present=_required_int(comparison, "first_only_present"),
                second_only_present=_required_int(comparison, "second_only_present"),
                unknown_involved=_required_int(comparison, "unknown_involved"),
                not_applicable_involved=_required_int(
                    comparison, "not_applicable_involved"
                ),
                assessable_both=_required_int(comparison, "assessable_both"),
            ),
            differing_categories=_string_tuple(data, "differing_categories"),
            ligand_feature_differences=tuple(
                (
                    _required_string(_as_mapping(item, "ligand feature difference"), "feature"),
                    _required_int(
                        _as_mapping(item, "ligand feature difference"), "reference_count"
                    ),
                    _required_int(
                        _as_mapping(item, "ligand feature difference"), "state_count"
                    ),
                )
                for item in _as_list(
                    data.get("ligand_feature_differences"), "ligand_feature_differences"
                )
            ),
            contact_feature_differences=tuple(
                (
                    _required_string(_as_mapping(item, "contact feature difference"), "feature"),
                    _required_int(
                        _as_mapping(item, "contact feature difference"), "reference_count"
                    ),
                    _required_int(
                        _as_mapping(item, "contact feature difference"), "state_count"
                    ),
                )
                for item in _as_list(
                    data.get("contact_feature_differences"), "contact_feature_differences"
                )
            ),
        )


@dataclass(frozen=True)
class StateInteractionRecord:
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
    accession: str
    selection: LigandSelection
    selected_state_ids: tuple[str, ...]
    reference_state_id: str
    criteria: StateInteractionCriteria
    assessments: tuple[StateInteractionAssessment, ...]
    comparisons: tuple[StateInteractionComparison, ...]
    decisions: tuple[DecisionRecord, ...]
    provenance: ProvenanceRecord
    validation: ValidationRecord
    discrepancies: tuple[Any, ...]
    schema_version: str = STATE_INTERACTION_SCHEMA_VERSION

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
            "accession": self.accession,
            "selection": self.selection.to_dict(),
            "selected_state_ids": list(self.selected_state_ids),
            "reference_state_id": self.reference_state_id,
            "criteria": self.criteria.to_dict(),
            "assessments": [item.to_dict() for item in self.assessments],
            "comparisons": [item.to_dict() for item in self.comparisons],
            "decisions": [decision.to_dict() for decision in self.decisions],
            "provenance": self.provenance.to_dict(),
            "validation": self.validation.to_dict(),
            "discrepancies": list(self.discrepancies),
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> StateInteractionRecord:
        data = _as_mapping(data, "state interaction record")
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
            chemical_state_record_path=_required_string(
                data, "chemical_state_record_path"
            ),
            chemical_state_record_sha256=_required_string(
                data, "chemical_state_record_sha256"
            ),
            accession=_required_string(data, "accession"),
            selection=LigandSelection.from_dict(
                _as_mapping(data.get("selection"), "selection")
            ),
            selected_state_ids=_string_tuple(data, "selected_state_ids"),
            reference_state_id=_required_string(data, "reference_state_id"),
            criteria=StateInteractionCriteria.from_dict(
                _as_mapping(data.get("criteria"), "criteria")
            ),
            assessments=tuple(
                StateInteractionAssessment.from_dict(
                    _as_mapping(item, "state interaction assessment")
                )
                for item in _as_list(data.get("assessments"), "assessments")
            ),
            comparisons=tuple(
                StateInteractionComparison.from_dict(
                    _as_mapping(item, "state interaction comparison")
                )
                for item in _as_list(data.get("comparisons"), "comparisons")
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
            discrepancies=tuple(_as_mapping(item, "discrepancy") for item in _as_list(data.get("discrepancies"), "discrepancies")),
        )


def _summary(atoms: Sequence[Any]) -> FeatureSummary:
    aromatic = sum(1 for atom in atoms if atom.aromatic_ring_ids)
    charged = sum(
        1 for atom in atoms if atom.formal_charge is not None and atom.formal_charge != 0
    )
    charges = [atom.formal_charge for atom in atoms]
    total_charge = None if any(charge is None for charge in charges) else sum(charges)
    return FeatureSummary(
        atom_count=len(atoms),
        donor_count=sum(1 for atom in atoms if atom.is_donor),
        acceptor_count=sum(1 for atom in atoms if atom.is_acceptor),
        hydrophobic_count=sum(1 for atom in atoms if atom.is_hydrophobic),
        aromatic_atom_count=aromatic,
        charged_atom_count=charged,
        formal_charge=total_charge,
    )


def assess_state_features(
    *,
    state: ChemicalState,
    reference_to_deposited: Mapping[int, str],
    ligand_atoms: Sequence[CoordinateAtom],
    contact_ligand_atoms: set[str],
    criteria: StateInteractionCriteria,
) -> StateInteractionAssessment:
    """Assess one state's ligand features against the experimental contact set."""

    deposited_atoms = {atom.atom_name: atom for atom in ligand_atoms}
    state_atoms_by_deposited: dict[str, Any] = {}
    for atom in state.atoms:
        reference_index = dict(state.reference_atom_mapping).get(atom.atom_id)
        if reference_index is None:
            continue
        deposited_name = reference_to_deposited.get(reference_index)
        if deposited_name is None:
            continue
        state_atoms_by_deposited[deposited_name] = atom

    mapped_atoms = tuple(
        state_atoms_by_deposited[name]
        for name in deposited_atoms
        if name in state_atoms_by_deposited
    )
    contact_pairs = tuple(
        (name, state_atoms_by_deposited[name])
        for name in contact_ligand_atoms
        if name in state_atoms_by_deposited
    )
    contact_atoms = tuple(atom for _, atom in contact_pairs)
    ligand_summary = _summary(mapped_atoms)
    contact_summary = _summary(contact_atoms)
    contact_atom_names = tuple(name for name, _ in contact_pairs)
    donor_contact_atoms = tuple(
        name for name, atom in contact_pairs if atom.is_donor
    )
    acceptor_contact_atoms = tuple(
        name for name, atom in contact_pairs if atom.is_acceptor
    )
    hydrophobic_contact_atoms = tuple(
        name for name, atom in contact_pairs if atom.is_hydrophobic
    )
    aromatic_contact_atoms = tuple(
        name for name, atom in contact_pairs if atom.aromatic_ring_ids
    )
    charged_contact_atoms = tuple(
        name
        for name, atom in contact_pairs
        if atom.formal_charge is not None and atom.formal_charge != 0
    )

    donor_acceptor_known = all(
        atom.is_donor is not None and atom.is_acceptor is not None
        for atom in mapped_atoms
    )
    donor_acceptor_present = any(
        atom.is_donor or atom.is_acceptor for atom in mapped_atoms
    )
    ionic_known = all(atom.formal_charge is not None for atom in mapped_atoms)
    ionic_present = any(
        atom.formal_charge not in (None, 0) for atom in mapped_atoms
    )
    hydrophobic_known = all(
        atom.is_hydrophobic is not None for atom in mapped_atoms
    )
    hydrophobic_present = any(atom.is_hydrophobic for atom in mapped_atoms)
    aromatic_present = any(atom.aromatic_ring_ids for atom in mapped_atoms)

    observations = [
        EvidenceObservation(
            feature=EvidenceFeature("category", "ligand_donor_acceptor"),
            state=(
                "present"
                if donor_acceptor_known and donor_acceptor_present
                else "absent"
                if donor_acceptor_known
                else "unknown"
            ),
            assessable=donor_acceptor_known,
            evidence=(
                f"{ligand_summary.donor_count} donors, "
                f"{ligand_summary.acceptor_count} acceptors"
            ),
            criteria=criteria.feature_source,
            provenance=state.label,
        ),
        EvidenceObservation(
            feature=EvidenceFeature("category", "ligand_ionic"),
            state=(
                "present"
                if ionic_known and ionic_present
                else "absent"
                if ionic_known
                else "unknown"
            ),
            assessable=ionic_known,
            evidence=(
                f"formal charge {ligand_summary.formal_charge!r}; "
                f"{ligand_summary.charged_atom_count} charged atoms"
            ),
            criteria=criteria.feature_source,
            provenance=state.label,
        ),
        EvidenceObservation(
            feature=EvidenceFeature("category", "ligand_hydrophobic"),
            state=(
                "present"
                if hydrophobic_known and hydrophobic_present
                else "absent"
                if hydrophobic_known
                else "unknown"
            ),
            assessable=hydrophobic_known,
            evidence=f"{ligand_summary.hydrophobic_count} hydrophobic atoms",
            criteria=criteria.feature_source,
            provenance=state.label,
        ),
        EvidenceObservation(
            feature=EvidenceFeature("category", "ligand_aromatic"),
            state="present" if aromatic_present else "absent",
            assessable=True,
            evidence=f"{ligand_summary.aromatic_atom_count} aromatic atoms",
            criteria=criteria.feature_source,
            provenance=state.label,
        ),
        EvidenceObservation(
            feature=EvidenceFeature("category", "geometric_contact"),
            state="present" if contact_atoms else "absent",
            assessable=True,
            evidence=f"{len(contact_atoms)} contact ligand atoms",
            criteria=(
                f"{criteria.geometry_source}; contact cutoff "
                f"{criteria.contact_cutoff_angstrom} A"
            ),
            provenance="Batch 4 binding-site record",
        ),
        EvidenceObservation(
            feature=EvidenceFeature("category", "interaction_classification"),
            state="unknown",
            assessable=False,
            evidence="protein chemical features and aligned hydrogen positions are not supplied",
            criteria=criteria.classification_limitations,
            provenance="state-aware interaction assessment",
        ),
    ]
    criteria_payload = criteria.to_dict()
    fingerprint = build_evidence_fingerprint(
        observations,
        schema_id="state-interaction-v1",
        criteria_fingerprint=criteria_fingerprint(criteria_payload),
        criteria_json=json.dumps(criteria_payload, sort_keys=True, separators=(",", ":")),
    )
    return StateInteractionAssessment(
        state_id=state.state_id,
        kind=state.kind,
        label=state.label,
        formula=state.formula,
        proton_count=state.proton_count,
        formal_charge=state.formal_charge,
        ligand_feature_summary=ligand_summary,
        contact_feature_summary=contact_summary,
        contact_atom_count=len(contact_ligand_atoms),
        contact_atom_names=contact_atom_names,
        donor_contact_atoms=donor_contact_atoms,
        acceptor_contact_atoms=acceptor_contact_atoms,
        hydrophobic_contact_atoms=hydrophobic_contact_atoms,
        aromatic_contact_atoms=aromatic_contact_atoms,
        charged_contact_atoms=charged_contact_atoms,
        classified_interaction_participants=(),
        fingerprint=fingerprint,
    )


def _feature_count_differences(
    first: FeatureSummary,
    second: FeatureSummary,
) -> tuple[tuple[str, int, int], ...]:
    fields = (
        "donor_count",
        "acceptor_count",
        "hydrophobic_count",
        "aromatic_atom_count",
        "charged_atom_count",
    )
    return tuple(
        (name, getattr(first, name), getattr(second, name))
        for name in fields
        if getattr(first, name) != getattr(second, name)
    )


def compare_state_assessment(
    reference: StateInteractionAssessment,
    other: StateInteractionAssessment,
) -> StateInteractionComparison:
    comparison = compare_evidence_fingerprints(reference.fingerprint, other.fingerprint)
    differing = []
    for first, second in zip(
        reference.fingerprint.schema.features,
        other.fingerprint.schema.features,
        strict=True,
    ):
        if (
            reference.fingerprint.state_for(first) != other.fingerprint.state_for(second)
            or reference.fingerprint.applicable_for(first)
            != other.fingerprint.applicable_for(second)
        ):
            differing.append(first.interaction_type)
    return StateInteractionComparison(
        reference_state_id=reference.state_id,
        state_id=other.state_id,
        comparison=comparison,
        differing_categories=tuple(differing),
        ligand_feature_differences=_feature_count_differences(
            reference.ligand_feature_summary, other.ligand_feature_summary
        ),
        contact_feature_differences=_feature_count_differences(
            reference.contact_feature_summary, other.contact_feature_summary
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


def build_state_interaction_record(
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
    accession: str,
    selection: LigandSelection,
    selected_state_ids: tuple[str, ...],
    reference_state_id: str,
    criteria: StateInteractionCriteria,
    assessments: tuple[StateInteractionAssessment, ...],
    comparisons: tuple[StateInteractionComparison, ...],
) -> StateInteractionRecord:
    assessment_ids = [item.state_id for item in assessments]
    if len(set(assessment_ids)) != len(assessment_ids):
        raise ValueError("assessment state_id values must be unique.")
    if reference_state_id not in assessment_ids:
        raise ValueError("reference_state_id must name one of the assessments.")
    if not set(selected_state_ids).issubset(assessment_ids):
        raise ValueError("selected_state_ids must name assessed states.")

    checks = (
        _check(
            "provenance",
            "state_interactions.source_records",
            bool(
                source_manifest_path
                and ligand_identity_record_path
                and binding_site_record_path
                and chemical_state_record_path
            ),
            subject=accession,
            criterion="experimental and state records are referenced",
            observed="four input records present",
            expected="manifest, identity, binding site, and state enumeration",
            detail="states are interpreted against the same experimental geometry",
        ),
        _check(
            "state",
            "state_interactions.explicit_selection",
            bool(selected_state_ids),
            subject=selection.component_id,
            criterion="chemical states are selected explicitly for this analysis",
            observed=", ".join(selected_state_ids),
            expected="one or more named states",
            detail="no default state is chosen by the record builder",
        ),
        _check(
            "mapping",
            "state_interactions.atom_mapping_complete",
            all(
                assessment.contact_atom_count
                == assessment.contact_feature_summary.atom_count
                for assessment in assessments
            ),
            subject=selection.component_id,
            criterion="contact ligand atoms map into each selected state",
            observed=", ".join(
                f"{item.state_id}: {item.contact_atom_count}" for item in assessments
            ),
            expected="all contact atoms mapped",
            detail="state atoms are mapped through the Batch 2 reference indices",
        ),
        _check(
            "interpretation",
            "state_interactions.geometry_separate",
            True,
            subject=selection.component_id,
            criterion="geometric contacts are recorded separately from classification",
            observed="geometric contact category is independent",
            expected="geometry and chemical interpretation remain distinct",
            detail="a contact is not promoted to a classified interaction",
        ),
        _check(
            "interpretation",
            "state_interactions.no_ranking",
            True,
            subject=selection.component_id,
            criterion="states are compared without ranking or population weighting",
            observed="no score, ranking, or affinity is produced",
            expected="traceable state differences only",
            detail="enumeration and geometry do not establish thermodynamic dominance",
        ),
    )

    decisions = (
        DecisionRecord(
            decision="interpreted the same experimental contacts under each selected state",
            actor="state-aware interaction assessment",
            rationale="chemical state changes feature eligibility without changing the deposited geometry",
        ),
        DecisionRecord(
            decision="kept geometric contact and interaction classification separate",
            actor="state-aware interaction assessment",
            rationale="a contact is not a classified interaction",
        ),
        DecisionRecord(
            decision="did not rank states or estimate populations",
            actor="state-aware interaction assessment",
            rationale="no pH-dependent thermodynamic model was supplied",
        ),
    )

    provenance = ProvenanceRecord(
        records=(
            ProvenanceEntry(
                "experimental",
                "binding-site contacts",
                f"{sum(a.contact_atom_count for a in assessments)} state-mapped contact atoms across {len(assessments)} states.",
                actor="RCSB PDB deposition",
            ),
            ProvenanceEntry(
                "calculated",
                "state-specific features",
                (
                    f"{len(assessments)} selected states; features from the enumerated "
                    "chemical-state record."
                ),
                actor="state-aware interaction assessment",
            ),
            ProvenanceEntry(
                "predicted",
                "none",
                "No affinity, population, or free-energy estimate is produced.",
                actor="state-aware interaction assessment",
            ),
            ProvenanceEntry(
                "hypothesised",
                "state alternatives",
                "Differences between states are hypotheses for later evaluation.",
                actor="state-aware interaction assessment",
            ),
        )
    )

    return StateInteractionRecord(
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
        accession=accession,
        selection=selection,
        selected_state_ids=selected_state_ids,
        reference_state_id=reference_state_id,
        criteria=criteria,
        assessments=assessments,
        comparisons=comparisons,
        decisions=decisions,
        provenance=provenance,
        validation=ValidationRecord(checks=checks),
        discrepancies=(),
    )


def write_state_interactions(path: str | Path, record: StateInteractionRecord) -> Path:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(record.to_dict(), indent=2, sort_keys=False) + "\n",
        encoding="utf-8",
    )
    return destination


def load_state_interactions(path: str | Path) -> StateInteractionRecord:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    return StateInteractionRecord.from_dict(data)


def validate_state_interactions(record: StateInteractionRecord) -> list[str]:
    errors: list[str] = []
    if record.schema_version != STATE_INTERACTION_SCHEMA_VERSION:
        errors.append(
            f"schema_version must be {STATE_INTERACTION_SCHEMA_VERSION!r}, got "
            f"{record.schema_version!r}."
        )
    for name in (
        "source_manifest_path",
        "ligand_identity_record_path",
        "ligand_coordinate_record_path",
        "binding_site_record_path",
        "chemical_state_record_path",
    ):
        if not getattr(record, name):
            errors.append(f"{name} is required.")
    for name in (
        "source_manifest_sha256",
        "ligand_identity_record_sha256",
        "ligand_coordinate_record_sha256",
        "binding_site_record_sha256",
        "chemical_state_record_sha256",
    ):
        if not _SHA256_PATTERN.fullmatch(getattr(record, name)):
            errors.append(f"{name} must be a lowercase SHA-256 digest.")
    if not record.accession:
        errors.append("accession is required.")
    if not record.selected_state_ids:
        errors.append("selected_state_ids must not be empty.")

    assessment_ids = [item.state_id for item in record.assessments]
    if len(set(assessment_ids)) != len(assessment_ids):
        errors.append("assessment state_id values must be unique.")
    if record.reference_state_id not in assessment_ids:
        errors.append("reference_state_id must name an assessment.")
    if not set(record.selected_state_ids).issubset(assessment_ids):
        errors.append("selected_state_ids must name assessed states.")
    for index, assessment in enumerate(record.assessments):
        if not assessment.formula or not assessment.label:
            errors.append(f"assessments[{index}] requires label and formula.")
        if assessment.contact_atom_count < 0:
            errors.append(f"assessments[{index}].contact_atom_count cannot be negative.")
        if len(assessment.contact_atom_names) != assessment.contact_feature_summary.atom_count:
            errors.append(
                f"assessments[{index}] mapped contact names must match the contact feature count."
            )
        contact_names = set(assessment.contact_atom_names)
        for field_name in (
            "donor_contact_atoms",
            "acceptor_contact_atoms",
            "hydrophobic_contact_atoms",
            "aromatic_contact_atoms",
            "charged_contact_atoms",
        ):
            if not set(getattr(assessment, field_name)).issubset(contact_names):
                errors.append(
                    f"assessments[{index}].{field_name} must be a subset of contact atoms."
                )
        if len(assessment.fingerprint.schema.features) != len(CATEGORY_ORDER):
            errors.append(f"assessments[{index}] has an unexpected feature count.")

    comparison_ids = {(item.reference_state_id, item.state_id) for item in record.comparisons}
    if len(comparison_ids) != len(record.comparisons):
        errors.append("comparisons must be unique by reference and state.")
    for index, comparison in enumerate(record.comparisons):
        if comparison.reference_state_id not in assessment_ids:
            errors.append(f"comparisons[{index}].reference_state_id is unknown.")
        if comparison.state_id not in assessment_ids:
            errors.append(f"comparisons[{index}].state_id is unknown.")

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

    return errors
