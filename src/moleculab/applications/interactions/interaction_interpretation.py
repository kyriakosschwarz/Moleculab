# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Chemistry-aware interpretation of ligand-protein contacts.

This module consumes explicit atom features and interaction criteria.  It does
not infer donor/acceptor status from elements, assign formal charges, perceive
aromatic rings, or convert proximity into affinity.

Interaction types are classified only when the required evidence is supplied:

- hydrogen bond: donor/acceptor flags and hydrogen positions;
- ionic: explicit formal charges on both partners;
- hydrophobic: explicit hydrophobic flags;
- aromatic: explicit ring centroids and normals on both sides;
- water-mediated: explicit water and polar features;
- metal coordination: explicit metal identity; coordination geometry remains
  unknown unless a future adapter supplies it.

Missing evidence produces ``unknown`` or ``not_applicable`` records, never an
implicit positive classification.
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

from ...analysis.interactions import (
    HydrogenBondCriteria,
    InteractionAtom,
    detect_hydrogen_bonds,
)
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

INTERACTION_SCHEMA_VERSION = "1.0"
INTERACTION_STATUSES = ("classified", "not_classified", "unknown", "not_applicable")
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


def _optional_float(data: Mapping[str, Any], key: str) -> float | None:
    value = data.get(key)
    return None if value is None else float(value)


def _required_bool(data: Mapping[str, Any], key: str) -> bool:
    value = data.get(key)
    if not isinstance(value, bool):
        raise ValueError(f"{key} must be a boolean.")
    return value


def _optional_bool(data: Mapping[str, Any], key: str) -> bool | None:
    value = data.get(key)
    return None if value is None else bool(value)


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
class AtomFeature:
    atom_id: str
    side: str
    element: str
    position_nm: tuple[float, float, float]
    residue_id: str | None = None
    is_donor: bool | None = None
    is_acceptor: bool | None = None
    formal_charge: int | None = None
    is_hydrophobic: bool | None = None
    aromatic_ring_id: str | None = None
    is_metal: bool | None = None
    hydrogen_positions_nm: tuple[tuple[float, float, float], ...] = ()
    feature_source: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "atom_id": self.atom_id,
            "side": self.side,
            "element": self.element,
            "position_nm": list(self.position_nm),
            "residue_id": self.residue_id,
            "is_donor": self.is_donor,
            "is_acceptor": self.is_acceptor,
            "formal_charge": self.formal_charge,
            "is_hydrophobic": self.is_hydrophobic,
            "aromatic_ring_id": self.aromatic_ring_id,
            "is_metal": self.is_metal,
            "hydrogen_positions_nm": [
                list(position) for position in self.hydrogen_positions_nm
            ],
            "feature_source": self.feature_source,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> AtomFeature:
        data = _as_mapping(data, "atom feature")
        return cls(
            atom_id=_required_string(data, "atom_id"),
            side=_required_string(data, "side"),
            element=_required_string(data, "element"),
            position_nm=_vector3(data, "position_nm"),
            residue_id=_optional_string(data, "residue_id"),
            is_donor=_optional_bool(data, "is_donor"),
            is_acceptor=_optional_bool(data, "is_acceptor"),
            formal_charge=_optional_int(data, "formal_charge"),
            is_hydrophobic=_optional_bool(data, "is_hydrophobic"),
            aromatic_ring_id=_optional_string(data, "aromatic_ring_id"),
            is_metal=_optional_bool(data, "is_metal"),
            hydrogen_positions_nm=_vector_tuple(data, "hydrogen_positions_nm"),
            feature_source=_optional_string(data, "feature_source"),
        )


@dataclass(frozen=True)
class AromaticRing:
    ring_id: str
    side: str
    atom_ids: tuple[str, ...]
    centroid_nm: tuple[float, float, float]
    normal: tuple[float, float, float]
    feature_source: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "ring_id": self.ring_id,
            "side": self.side,
            "atom_ids": list(self.atom_ids),
            "centroid_nm": list(self.centroid_nm),
            "normal": list(self.normal),
            "feature_source": self.feature_source,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> AromaticRing:
        data = _as_mapping(data, "aromatic ring")
        return cls(
            ring_id=_required_string(data, "ring_id"),
            side=_required_string(data, "side"),
            atom_ids=_string_tuple(data, "atom_ids"),
            centroid_nm=_vector3(data, "centroid_nm"),
            normal=_vector3(data, "normal"),
            feature_source=_optional_string(data, "feature_source"),
        )


@dataclass(frozen=True)
class InteractionCriteria:
    h_bond_max_distance_nm: float = 0.35
    h_bond_min_angle_degrees: float = 120.0
    ionic_max_distance_nm: float = 0.4
    hydrophobic_max_distance_nm: float = 0.4
    aromatic_centroid_max_distance_nm: float = 0.55
    aromatic_normal_angle_max_degrees: float = 30.0
    water_max_distance_nm: float = 0.35
    metal_max_distance_nm: float = 0.3

    def to_dict(self) -> dict[str, Any]:
        return {
            "h_bond_max_distance_nm": self.h_bond_max_distance_nm,
            "h_bond_min_angle_degrees": self.h_bond_min_angle_degrees,
            "ionic_max_distance_nm": self.ionic_max_distance_nm,
            "hydrophobic_max_distance_nm": self.hydrophobic_max_distance_nm,
            "aromatic_centroid_max_distance_nm": self.aromatic_centroid_max_distance_nm,
            "aromatic_normal_angle_max_degrees": self.aromatic_normal_angle_max_degrees,
            "water_max_distance_nm": self.water_max_distance_nm,
            "metal_max_distance_nm": self.metal_max_distance_nm,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> InteractionCriteria:
        data = _as_mapping(data, "interaction criteria")
        return cls(
            h_bond_max_distance_nm=_required_float(data, "h_bond_max_distance_nm"),
            h_bond_min_angle_degrees=_required_float(
                data, "h_bond_min_angle_degrees"
            ),
            ionic_max_distance_nm=_required_float(data, "ionic_max_distance_nm"),
            hydrophobic_max_distance_nm=_required_float(
                data, "hydrophobic_max_distance_nm"
            ),
            aromatic_centroid_max_distance_nm=_required_float(
                data, "aromatic_centroid_max_distance_nm"
            ),
            aromatic_normal_angle_max_degrees=_required_float(
                data, "aromatic_normal_angle_max_degrees"
            ),
            water_max_distance_nm=_required_float(data, "water_max_distance_nm"),
            metal_max_distance_nm=_required_float(data, "metal_max_distance_nm"),
        )


@dataclass(frozen=True)
class InteractionRecord:
    interaction_type: str
    status: str
    participants: tuple[str, ...]
    distance_nm: float | None = None
    angle_degrees: float | None = None
    protein_residue: str | None = None
    ligand_atom: str | None = None
    water_id: str | None = None
    metal_id: str | None = None
    evidence: str = ""
    criteria: str = ""
    limitations: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "interaction_type": self.interaction_type,
            "status": self.status,
            "participants": list(self.participants),
            "distance_nm": self.distance_nm,
            "angle_degrees": self.angle_degrees,
            "protein_residue": self.protein_residue,
            "ligand_atom": self.ligand_atom,
            "water_id": self.water_id,
            "metal_id": self.metal_id,
            "evidence": self.evidence,
            "criteria": self.criteria,
            "limitations": self.limitations,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> InteractionRecord:
        data = _as_mapping(data, "interaction record")
        return cls(
            interaction_type=_required_string(data, "interaction_type"),
            status=_required_string(data, "status"),
            participants=_string_tuple(data, "participants"),
            distance_nm=_optional_float(data, "distance_nm"),
            angle_degrees=_optional_float(data, "angle_degrees"),
            protein_residue=_optional_string(data, "protein_residue"),
            ligand_atom=_optional_string(data, "ligand_atom"),
            water_id=_optional_string(data, "water_id"),
            metal_id=_optional_string(data, "metal_id"),
            evidence=_required_string(data, "evidence"),
            criteria=_required_string(data, "criteria"),
            limitations=_optional_string(data, "limitations"),
        )


@dataclass(frozen=True)
class InteractionDiscrepancy:
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
    def from_dict(cls, data: Mapping[str, Any]) -> InteractionDiscrepancy:
        data = _as_mapping(data, "interaction discrepancy")
        return cls(
            category=_required_string(data, "category"),
            status=_required_string(data, "status"),
            detail=_required_string(data, "detail"),
            observed=_optional_string(data, "observed"),
            expected=_optional_string(data, "expected"),
            resolution=_optional_string(data, "resolution"),
        )


@dataclass(frozen=True)
class InteractionInterpretationRecord:
    source_manifest_path: str
    source_manifest_sha256: str
    ligand_identity_record_path: str
    ligand_identity_record_sha256: str
    binding_site_record_path: str
    binding_site_record_sha256: str
    accession: str
    model_number: int
    protein_selection: ProteinSelection
    ligand_selection: LigandSelection
    criteria: InteractionCriteria
    input_contact_count: int
    features: tuple[AtomFeature, ...]
    aromatic_rings: tuple[AromaticRing, ...]
    interactions: tuple[InteractionRecord, ...]
    decisions: tuple[DecisionRecord, ...]
    provenance: ProvenanceRecord
    validation: ValidationRecord
    discrepancies: tuple[InteractionDiscrepancy, ...]
    schema_version: str = INTERACTION_SCHEMA_VERSION

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "source_manifest_path": self.source_manifest_path,
            "source_manifest_sha256": self.source_manifest_sha256,
            "ligand_identity_record_path": self.ligand_identity_record_path,
            "ligand_identity_record_sha256": self.ligand_identity_record_sha256,
            "binding_site_record_path": self.binding_site_record_path,
            "binding_site_record_sha256": self.binding_site_record_sha256,
            "accession": self.accession,
            "model_number": self.model_number,
            "protein_selection": self.protein_selection.to_dict(),
            "ligand_selection": self.ligand_selection.to_dict(),
            "criteria": self.criteria.to_dict(),
            "input_contact_count": self.input_contact_count,
            "features": [feature.to_dict() for feature in self.features],
            "aromatic_rings": [ring.to_dict() for ring in self.aromatic_rings],
            "interactions": [interaction.to_dict() for interaction in self.interactions],
            "decisions": [decision.to_dict() for decision in self.decisions],
            "provenance": self.provenance.to_dict(),
            "validation": self.validation.to_dict(),
            "discrepancies": [item.to_dict() for item in self.discrepancies],
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> InteractionInterpretationRecord:
        data = _as_mapping(data, "interaction interpretation record")
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
            accession=_required_string(data, "accession"),
            model_number=_required_int(data, "model_number"),
            protein_selection=ProteinSelection.from_dict(
                _as_mapping(data.get("protein_selection"), "protein_selection")
            ),
            ligand_selection=LigandSelection.from_dict(
                _as_mapping(data.get("ligand_selection"), "ligand_selection")
            ),
            criteria=InteractionCriteria.from_dict(
                _as_mapping(data.get("criteria"), "criteria")
            ),
            input_contact_count=_required_int(data, "input_contact_count"),
            features=tuple(
                AtomFeature.from_dict(_as_mapping(item, "atom feature"))
                for item in _as_list(data.get("features"), "features")
            ),
            aromatic_rings=tuple(
                AromaticRing.from_dict(_as_mapping(item, "aromatic ring"))
                for item in _as_list(data.get("aromatic_rings"), "aromatic_rings")
            ),
            interactions=tuple(
                InteractionRecord.from_dict(_as_mapping(item, "interaction record"))
                for item in _as_list(data.get("interactions"), "interactions")
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
                InteractionDiscrepancy.from_dict(
                    _as_mapping(item, "interaction discrepancy")
                )
                for item in _as_list(data.get("discrepancies"), "discrepancies")
            ),
        )


def _distance(first: tuple[float, float, float], second: tuple[float, float, float]) -> float:
    return float(
        np.linalg.norm(np.asarray(first, dtype=float) - np.asarray(second, dtype=float))
    )


def _angle_between(first: tuple[float, float, float], second: tuple[float, float, float]) -> float:
    left = np.asarray(first, dtype=float)
    right = np.asarray(second, dtype=float)
    denominator = float(np.linalg.norm(left) * np.linalg.norm(right))
    if denominator == 0.0:
        return 0.0
    cosine = float(np.dot(left, right) / denominator)
    return float(math.degrees(math.acos(max(-1.0, min(1.0, cosine)))))


def _feature_label(feature: AtomFeature) -> str:
    if feature.residue_id:
        return f"{feature.residue_id}:{feature.atom_id}"
    return feature.atom_id


def _interaction_atom(feature: AtomFeature) -> InteractionAtom:
    return InteractionAtom(
        position=np.asarray(feature.position_nm, dtype=float),
        element=feature.element,
        atom_type=None,
        residue_id=feature.residue_id,
        name=feature.atom_id,
        is_donor=bool(feature.is_donor),
        is_acceptor=bool(feature.is_acceptor),
        formal_charge=feature.formal_charge or 0,
        is_aromatic=feature.aromatic_ring_id is not None,
        is_hydrophobic=bool(feature.is_hydrophobic),
        hydrogen_positions=tuple(
            np.asarray(position, dtype=float)
            for position in feature.hydrogen_positions_nm
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


def _has_known_polar(features: Sequence[AtomFeature]) -> bool:
    return any(
        feature.is_donor is not None or feature.is_acceptor is not None
        for feature in features
    )


def _has_donor_hydrogens(features: Sequence[AtomFeature]) -> bool:
    return any(feature.is_donor and feature.hydrogen_positions_nm for feature in features)


def _feature_pairs(
    protein_features: Sequence[AtomFeature],
    ligand_features: Sequence[AtomFeature],
    cutoff_nm: float,
) -> list[tuple[AtomFeature, AtomFeature, float]]:
    pairs = []
    for protein in protein_features:
        for ligand in ligand_features:
            distance = _distance(protein.position_nm, ligand.position_nm)
            if distance <= cutoff_nm:
                pairs.append((protein, ligand, distance))
    return pairs


def _classify_hydrogen_bonds(
    protein_features: Sequence[AtomFeature],
    ligand_features: Sequence[AtomFeature],
    criteria: InteractionCriteria,
) -> list[InteractionRecord]:
    if not _has_known_polar(protein_features) or not _has_known_polar(ligand_features):
        return [
            InteractionRecord(
                "hydrogen_bond",
                "unknown",
                (),
                evidence="donor and acceptor assignments are incomplete",
                criteria=(
                    f"D-A <= {criteria.h_bond_max_distance_nm} nm; "
                    f"D-H...A >= {criteria.h_bond_min_angle_degrees} degrees"
                ),
                limitations="donor/acceptor identity was not supplied for both sides",
            )
        ]
    if not (
        _has_donor_hydrogens(protein_features)
        or _has_donor_hydrogens(ligand_features)
    ):
        return [
            InteractionRecord(
                "hydrogen_bond",
                "unknown",
                (),
                evidence="donor flags are present but no donor hydrogen positions are available",
                criteria=(
                    f"D-A <= {criteria.h_bond_max_distance_nm} nm; "
                    f"D-H...A >= {criteria.h_bond_min_angle_degrees} degrees"
                ),
                limitations="angular hydrogen-bond geometry cannot be evaluated",
            )
        ]
    hb_criteria = HydrogenBondCriteria(
        max_distance=criteria.h_bond_max_distance_nm,
        min_angle_degrees=criteria.h_bond_min_angle_degrees,
    )
    bonds = detect_hydrogen_bonds(
        [_interaction_atom(feature) for feature in protein_features],
        [_interaction_atom(feature) for feature in ligand_features],
        hb_criteria,
    )
    if not bonds:
        return [
            InteractionRecord(
                "hydrogen_bond",
                "not_classified",
                (),
                evidence="chemical features and hydrogen positions are present",
                criteria=(
                    f"D-A <= {criteria.h_bond_max_distance_nm} nm; "
                    f"D-H...A >= {criteria.h_bond_min_angle_degrees} degrees"
                ),
                limitations="no pair satisfied the explicit geometry",
            )
        ]
    return [
        InteractionRecord(
            "hydrogen_bond",
            "classified",
            (bond.protein_label, bond.ligand_label),
            distance_nm=bond.distance,
            angle_degrees=bond.angle_degrees,
            protein_residue=bond.protein_residue,
            ligand_atom=bond.ligand_label,
            evidence=(
                f"{'protein' if bond.donor_is_protein else 'ligand'} donor; "
                f"{bond.donor_label}...{bond.acceptor_label}"
            ),
            criteria=(
                f"D-A <= {criteria.h_bond_max_distance_nm} nm; "
                f"D-H...A >= {criteria.h_bond_min_angle_degrees} degrees"
            ),
        )
        for bond in bonds
    ]


def _classify_ionic(
    protein_features: Sequence[AtomFeature],
    ligand_features: Sequence[AtomFeature],
    criteria: InteractionCriteria,
) -> list[InteractionRecord]:
    if any(
        feature.formal_charge is None
        for feature in (*protein_features, *ligand_features)
    ):
        return [
            InteractionRecord(
                "ionic",
                "unknown",
                (),
                evidence="formal charges are not assigned for every candidate atom",
                criteria=f"opposite nonzero charges <= {criteria.ionic_max_distance_nm} nm",
                limitations="no charge model or assignment was supplied",
            )
        ]
    pairs = _feature_pairs(protein_features, ligand_features, criteria.ionic_max_distance_nm)
    records = [
        InteractionRecord(
            "ionic",
            "classified",
            (_feature_label(protein), _feature_label(ligand)),
            distance_nm=distance,
            protein_residue=protein.residue_id,
            ligand_atom=ligand.atom_id,
            evidence=(
                f"formal charges {protein.formal_charge:+d} and "
                f"{ligand.formal_charge:+d}"
            ),
            criteria=f"opposite nonzero charges <= {criteria.ionic_max_distance_nm} nm",
        )
        for protein, ligand, distance in pairs
        if protein.formal_charge and ligand.formal_charge
        and protein.formal_charge * ligand.formal_charge < 0
    ]
    if records:
        return records
    return [
        InteractionRecord(
            "ionic",
            "not_classified",
            (),
            evidence="formal charges are assigned",
            criteria=f"opposite nonzero charges <= {criteria.ionic_max_distance_nm} nm",
            limitations="no opposite-charge pair satisfied the distance criterion",
        )
    ]


def _classify_hydrophobic(
    protein_features: Sequence[AtomFeature],
    ligand_features: Sequence[AtomFeature],
    criteria: InteractionCriteria,
) -> list[InteractionRecord]:
    if any(
        feature.is_hydrophobic is None
        for feature in (*protein_features, *ligand_features)
    ):
        return [
            InteractionRecord(
                "hydrophobic",
                "unknown",
                (),
                evidence="hydrophobic flags are incomplete",
                criteria=f"explicit hydrophobic pairs <= {criteria.hydrophobic_max_distance_nm} nm",
                limitations="hydrophobic atom/residue classification was not supplied",
            )
        ]
    pairs = _feature_pairs(protein_features, ligand_features, criteria.hydrophobic_max_distance_nm)
    records = [
        InteractionRecord(
            "hydrophobic",
            "classified",
            (_feature_label(protein), _feature_label(ligand)),
            distance_nm=distance,
            protein_residue=protein.residue_id,
            ligand_atom=ligand.atom_id,
            evidence="both atoms carry explicit hydrophobic flags",
            criteria=f"explicit hydrophobic pairs <= {criteria.hydrophobic_max_distance_nm} nm",
        )
        for protein, ligand, distance in pairs
        if protein.is_hydrophobic and ligand.is_hydrophobic
    ]
    if records:
        return records
    return [
        InteractionRecord(
            "hydrophobic",
            "not_classified",
            (),
            evidence="hydrophobic flags are assigned",
            criteria=f"explicit hydrophobic pairs <= {criteria.hydrophobic_max_distance_nm} nm",
            limitations="no pair satisfied the distance criterion",
        )
    ]


def _classify_aromatic(
    protein_rings: Sequence[AromaticRing],
    ligand_rings: Sequence[AromaticRing],
    criteria: InteractionCriteria,
) -> list[InteractionRecord]:
    if not protein_rings or not ligand_rings:
        return [
            InteractionRecord(
                "aromatic",
                "unknown",
                (),
                evidence="aromatic ring definitions are incomplete",
                criteria=(
                    f"centroid distance <= {criteria.aromatic_centroid_max_distance_nm} nm; "
                    f"normal angle <= {criteria.aromatic_normal_angle_max_degrees} degrees"
                ),
                limitations="both protein and ligand ring definitions are required",
            )
        ]
    records = []
    for protein in protein_rings:
        for ligand in ligand_rings:
            distance = _distance(protein.centroid_nm, ligand.centroid_nm)
            angle = _angle_between(protein.normal, ligand.normal)
            if (
                distance <= criteria.aromatic_centroid_max_distance_nm
                and min(angle, 180.0 - angle)
                <= criteria.aromatic_normal_angle_max_degrees
            ):
                records.append(
                    InteractionRecord(
                        "aromatic",
                        "classified",
                        (protein.ring_id, ligand.ring_id),
                        distance_nm=distance,
                        angle_degrees=angle,
                        evidence="explicit ring centroids and normals",
                        criteria=(
                            f"centroid distance <= {criteria.aromatic_centroid_max_distance_nm} nm; "
                            f"normal angle <= {criteria.aromatic_normal_angle_max_degrees} degrees"
                        ),
                    )
                )
    if records:
        return records
    return [
        InteractionRecord(
            "aromatic",
            "not_classified",
            (),
            evidence="ring definitions are supplied for both sides",
            criteria=(
                f"centroid distance <= {criteria.aromatic_centroid_max_distance_nm} nm; "
                f"normal angle <= {criteria.aromatic_normal_angle_max_degrees} degrees"
            ),
            limitations="no ring pair satisfied the criteria",
        )
    ]


def _classify_water_mediated(
    protein_features: Sequence[AtomFeature],
    ligand_features: Sequence[AtomFeature],
    water_features: Sequence[AtomFeature],
    criteria: InteractionCriteria,
) -> list[InteractionRecord]:
    if not water_features:
        return [
            InteractionRecord(
                "water_mediated",
                "not_applicable",
                (),
                evidence="no water features supplied",
                criteria=f"water bridges <= {criteria.water_max_distance_nm} nm on both sides",
                limitations="no explicit water record was provided",
            )
        ]
    if (
        not _has_known_polar(protein_features)
        or not _has_known_polar(ligand_features)
        or any(
            water.is_donor is None and water.is_acceptor is None
            for water in water_features
        )
    ):
        return [
            InteractionRecord(
                "water_mediated",
                "unknown",
                (),
                evidence="water or polar atom features are incomplete",
                criteria=f"water bridges <= {criteria.water_max_distance_nm} nm on both sides",
                limitations="donor/acceptor identity is required for all participants",
            )
        ]
    records = []
    for water in water_features:
        for protein in protein_features:
            protein_distance = _distance(water.position_nm, protein.position_nm)
            if protein_distance > criteria.water_max_distance_nm:
                continue
            for ligand in ligand_features:
                ligand_distance = _distance(water.position_nm, ligand.position_nm)
                if ligand_distance > criteria.water_max_distance_nm:
                    continue
                if not (
                    (water.is_donor and protein.is_acceptor and ligand.is_acceptor)
                    or (water.is_acceptor and protein.is_donor and ligand.is_donor)
                ):
                    continue
                records.append(
                    InteractionRecord(
                        "water_mediated",
                        "classified",
                        (
                            _feature_label(protein),
                            water.atom_id,
                            _feature_label(ligand),
                        ),
                        distance_nm=max(protein_distance, ligand_distance),
                        protein_residue=protein.residue_id,
                        ligand_atom=ligand.atom_id,
                        water_id=water.atom_id,
                        evidence=(
                            f"water {water.atom_id} bridges protein and ligand "
                            "with explicit donor/acceptor roles"
                        ),
                        criteria=f"both water contacts <= {criteria.water_max_distance_nm} nm",
                    )
                )
    if records:
        return records
    return [
        InteractionRecord(
            "water_mediated",
            "not_classified",
            (),
            evidence="water and polar features are supplied",
            criteria=f"both water contacts <= {criteria.water_max_distance_nm} nm",
            limitations="no water bridge satisfied the criteria and role assignments",
        )
    ]


def _classify_metal(
    protein_features: Sequence[AtomFeature],
    ligand_features: Sequence[AtomFeature],
    metal_features: Sequence[AtomFeature],
    criteria: InteractionCriteria,
) -> list[InteractionRecord]:
    if not metal_features:
        return [
            InteractionRecord(
                "metal_coordination",
                "not_applicable",
                (),
                evidence="no explicit metal features supplied",
                criteria=f"metal distance <= {criteria.metal_max_distance_nm} nm",
                limitations="no metal identity is available",
            )
        ]
    records = []
    for metal in metal_features:
        for protein in protein_features:
            if _distance(metal.position_nm, protein.position_nm) <= criteria.metal_max_distance_nm:
                records.append(
                    InteractionRecord(
                        "metal_coordination",
                        "unknown",
                        (metal.atom_id, _feature_label(protein)),
                        distance_nm=_distance(metal.position_nm, protein.position_nm),
                        protein_residue=protein.residue_id,
                        metal_id=metal.atom_id,
                        evidence="explicit metal identity and distance",
                        criteria=f"metal distance <= {criteria.metal_max_distance_nm} nm",
                        limitations="coordination geometry and metal identity were not fully supplied",
                    )
                )
        for ligand in ligand_features:
            if _distance(metal.position_nm, ligand.position_nm) <= criteria.metal_max_distance_nm:
                records.append(
                    InteractionRecord(
                        "metal_coordination",
                        "unknown",
                        (metal.atom_id, _feature_label(ligand)),
                        distance_nm=_distance(metal.position_nm, ligand.position_nm),
                        ligand_atom=ligand.atom_id,
                        metal_id=metal.atom_id,
                        evidence="explicit metal identity and distance",
                        criteria=f"metal distance <= {criteria.metal_max_distance_nm} nm",
                        limitations="coordination geometry and metal identity were not fully supplied",
                    )
                )
    if records:
        return records
    return [
        InteractionRecord(
            "metal_coordination",
            "not_classified",
            (),
            evidence="metal features are supplied",
            criteria=f"metal distance <= {criteria.metal_max_distance_nm} nm",
            limitations="no candidate satisfied the distance criterion",
        )
    ]


def interpret_interactions(
    *,
    source_manifest_path: str,
    source_manifest_sha256: str,
    ligand_identity_record_path: str,
    ligand_identity_record_sha256: str,
    binding_site_record_path: str,
    binding_site_record_sha256: str,
    accession: str,
    model_number: int,
    protein_selection: ProteinSelection,
    ligand_selection: LigandSelection,
    input_contact_count: int,
    protein_features: Sequence[AtomFeature],
    ligand_features: Sequence[AtomFeature],
    water_features: Sequence[AtomFeature] = (),
    metal_features: Sequence[AtomFeature] = (),
    aromatic_rings: Sequence[AromaticRing] = (),
    criteria: InteractionCriteria | None = None,
) -> InteractionInterpretationRecord:
    """Classify interactions only where explicit evidence is available."""

    criteria = criteria or InteractionCriteria()
    if not protein_features:
        raise ValueError("protein_features must not be empty.")
    if not ligand_features:
        raise ValueError("ligand_features must not be empty.")

    protein_rings = tuple(ring for ring in aromatic_rings if ring.side == "protein")
    ligand_rings = tuple(ring for ring in aromatic_rings if ring.side == "ligand")
    interactions = tuple(
        [
            *_classify_hydrogen_bonds(protein_features, ligand_features, criteria),
            *_classify_ionic(protein_features, ligand_features, criteria),
            *_classify_hydrophobic(protein_features, ligand_features, criteria),
            *_classify_aromatic(protein_rings, ligand_rings, criteria),
            *_classify_water_mediated(
                protein_features, ligand_features, water_features, criteria
            ),
            *_classify_metal(
                protein_features, ligand_features, metal_features, criteria
            ),
        ]
    )
    classified_count = sum(1 for item in interactions if item.status == "classified")

    checks = (
        _check(
            "provenance",
            "interactions.source_records",
            bool(
                source_manifest_path
                and ligand_identity_record_path
                and binding_site_record_path
            ),
            subject=accession,
            criterion="manifest, identity, and binding-site records are recorded",
            observed="all record paths present",
            expected="three prior records",
            detail="the interpretation is traceable to the structural evidence chain",
        ),
        _check(
            "evidence",
            "interactions.feature_provenance",
            all(feature.feature_source for feature in (*protein_features, *ligand_features)),
            subject=ligand_selection.component_id,
            criterion="every chemical feature has an explicit provenance",
            observed=(
                f"{sum(1 for f in (*protein_features, *ligand_features) if f.feature_source)} "
                f"of {len(protein_features) + len(ligand_features)} feature records"
            ),
            expected="all feature records",
            detail="missing feature provenance is reported, not inferred",
            not_run=not all(
                feature.feature_source
                for feature in (*protein_features, *ligand_features)
            ),
        ),
        _check(
            "classification",
            "interactions.classified",
            classified_count > 0,
            subject=ligand_selection.component_id,
            criterion="at least one interaction satisfies the explicit criteria",
            observed=f"{classified_count} classified interactions",
            expected="one or more, or explicit unknown/not-applicable records",
            detail="unknown chemistry remains unknown",
            not_run=classified_count == 0,
        ),
        _check(
            "interpretation",
            "interactions.proximity_not_affinity",
            True,
            subject=ligand_selection.component_id,
            criterion="interaction classification is not converted into affinity",
            observed="no affinity or free-energy estimate is produced",
            expected="structural and chemical interpretation only",
            detail="interaction records are not binding free energies",
        ),
    )

    discrepancies: list[InteractionDiscrepancy] = []
    if not all(feature.feature_source for feature in (*protein_features, *ligand_features)):
        discrepancies.append(
            InteractionDiscrepancy(
                "feature_provenance",
                "open",
                "Some atom features do not record their chemical provenance.",
                observed="incomplete feature provenance",
                expected="explicit source for every feature",
                resolution="No feature was inferred from an element or from proximity.",
            )
        )
    if classified_count == 0:
        discrepancies.append(
            InteractionDiscrepancy(
                "insufficient_evidence",
                "not_applicable",
                "No interaction type has enough explicit evidence to be classified.",
                observed="0 classified interactions",
                expected="donor/acceptor, charge, hydrophobic, ring, water, or metal evidence",
                resolution="The record preserves the unknown classifications rather than guessing.",
            )
        )

    decisions = (
        DecisionRecord(
            decision="classified interactions only where explicit chemical evidence was supplied",
            actor="interaction interpretation",
            rationale=(
                "donor/acceptor identities, formal charges, hydrophobic flags, ring "
                "geometry, water roles, and metal identity were not inferred from "
                "elements or proximity"
            ),
        ),
        DecisionRecord(
            decision="preserved unknown and not-applicable interaction records",
            actor="interaction interpretation",
            rationale=(
                "missing evidence remains visible rather than being converted into a "
                "positive interaction"
            ),
        ),
        DecisionRecord(
            decision="did not convert interactions into affinity",
            actor="interaction interpretation",
            rationale="interaction geometry is not a binding free energy",
        ),
    )

    provenance = ProvenanceRecord(
        records=(
            ProvenanceEntry(
                "experimental",
                "contact and coordinate records",
                (
                    f"Interpretation consumes {input_contact_count} contacts from "
                    "the binding-site record and the earlier ligand identity and "
                    "coordinate records."
                ),
                actor="RCSB PDB deposition",
            ),
            ProvenanceEntry(
                "calculated",
                "interaction classification",
                (
                    f"{classified_count} classified interactions; "
                    f"{len(interactions) - classified_count} unknown, not classified, "
                    "or not applicable records."
                ),
                actor="interaction interpretation",
            ),
            ProvenanceEntry(
                "predicted",
                "none",
                "No binding affinity or free energy is produced.",
                actor="interaction interpretation",
            ),
            ProvenanceEntry(
                "hypothesised",
                "none",
                "No mechanism is inferred from an interaction record.",
                actor="interaction interpretation",
            ),
        )
    )

    return InteractionInterpretationRecord(
        source_manifest_path=source_manifest_path,
        source_manifest_sha256=source_manifest_sha256,
        ligand_identity_record_path=ligand_identity_record_path,
        ligand_identity_record_sha256=ligand_identity_record_sha256,
        binding_site_record_path=binding_site_record_path,
        binding_site_record_sha256=binding_site_record_sha256,
        accession=accession,
        model_number=model_number,
        protein_selection=protein_selection,
        ligand_selection=ligand_selection,
        criteria=criteria,
        input_contact_count=input_contact_count,
        features=tuple([*protein_features, *ligand_features, *water_features, *metal_features]),
        aromatic_rings=tuple(aromatic_rings),
        interactions=interactions,
        decisions=decisions,
        provenance=provenance,
        validation=ValidationRecord(checks=checks),
        discrepancies=tuple(discrepancies),
    )


def write_interactions(path: str | Path, record: InteractionInterpretationRecord) -> Path:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(record.to_dict(), indent=2, sort_keys=False) + "\n",
        encoding="utf-8",
    )
    return destination


def load_interactions(path: str | Path) -> InteractionInterpretationRecord:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    return InteractionInterpretationRecord.from_dict(data)


def validate_interactions(record: InteractionInterpretationRecord) -> list[str]:
    """Return schema and internal-consistency errors for *record*."""

    errors: list[str] = []
    if record.schema_version != INTERACTION_SCHEMA_VERSION:
        errors.append(
            f"schema_version must be {INTERACTION_SCHEMA_VERSION!r}, got "
            f"{record.schema_version!r}."
        )
    for name in (
        "source_manifest_path",
        "ligand_identity_record_path",
        "binding_site_record_path",
    ):
        if not getattr(record, name):
            errors.append(f"{name} is required.")
    for name in (
        "source_manifest_sha256",
        "ligand_identity_record_sha256",
        "binding_site_record_sha256",
    ):
        if not _SHA256_PATTERN.fullmatch(getattr(record, name)):
            errors.append(f"{name} must be a lowercase SHA-256 digest.")
    if not record.accession:
        errors.append("accession is required.")
    if record.model_number <= 0:
        errors.append("model_number must be positive.")
    if record.input_contact_count < 0:
        errors.append("input_contact_count cannot be negative.")

    criteria = record.criteria
    for name in (
        "h_bond_max_distance_nm",
        "ionic_max_distance_nm",
        "hydrophobic_max_distance_nm",
        "aromatic_centroid_max_distance_nm",
        "water_max_distance_nm",
        "metal_max_distance_nm",
    ):
        if getattr(criteria, name) <= 0:
            errors.append(f"criteria.{name} must be positive.")
    if not 0 < criteria.h_bond_min_angle_degrees <= 180:
        errors.append("criteria.h_bond_min_angle_degrees must be in (0, 180].")
    if not 0 < criteria.aromatic_normal_angle_max_degrees <= 90:
        errors.append(
            "criteria.aromatic_normal_angle_max_degrees must be in (0, 90]."
        )

    feature_ids = [feature.atom_id for feature in record.features]
    if len(set(feature_ids)) != len(feature_ids):
        errors.append("feature atom_id values must be unique.")
    for index, feature in enumerate(record.features):
        if feature.side not in ("protein", "ligand", "water", "metal"):
            errors.append(f"features[{index}].side is invalid.")
        if not feature.element:
            errors.append(f"features[{index}].element is required.")
        position = np.asarray(feature.position_nm, dtype=float)
        if position.shape != (3,) or not np.all(np.isfinite(position)):
            errors.append(f"features[{index}].position_nm must be a finite three-vector.")
        for hydrogen_index, hydrogen in enumerate(feature.hydrogen_positions_nm):
            coordinates = np.asarray(hydrogen, dtype=float)
            if coordinates.shape != (3,) or not np.all(np.isfinite(coordinates)):
                errors.append(
                    f"features[{index}].hydrogen_positions_nm[{hydrogen_index}] "
                    "must be a finite three-vector."
                )

    for index, ring in enumerate(record.aromatic_rings):
        if ring.side not in ("protein", "ligand"):
            errors.append(f"aromatic_rings[{index}].side is invalid.")
        if not ring.atom_ids:
            errors.append(f"aromatic_rings[{index}].atom_ids must not be empty.")
        normal = np.asarray(ring.normal, dtype=float)
        if normal.shape != (3,) or not np.all(np.isfinite(normal)) or np.linalg.norm(normal) == 0:
            errors.append(f"aromatic_rings[{index}].normal must be a finite nonzero vector.")

    for index, interaction in enumerate(record.interactions):
        if interaction.status not in INTERACTION_STATUSES:
            errors.append(f"interactions[{index}].status is invalid.")
        if not interaction.interaction_type or not interaction.evidence:
            errors.append(f"interactions[{index}] requires type and evidence.")
        if interaction.distance_nm is not None and interaction.distance_nm <= 0:
            errors.append(f"interactions[{index}].distance_nm must be positive.")
        if interaction.angle_degrees is not None and not (
            0.0 <= interaction.angle_degrees <= 180.0
        ):
            errors.append(f"interactions[{index}].angle_degrees must be in [0, 180].")

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
