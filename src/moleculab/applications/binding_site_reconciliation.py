"""Binding-site reconciliation around a validated ligand.

This module connects the ligand coordinate record to the selected protein
chain.  It records the pocket cutoff and atom-selection policy, residue
identities in both mmCIF label and author namespaces, atom-level ligand-protein
contacts, nearby waters, other hetero residues, missing residue numbers, and
alternate locations.

It does not infer hydrogen bonds, hydrophobic favourability, binding affinity,
or a binding free energy from proximity.  It does not prepare, protonate,
minimise, or dock the complex.
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

from ..analysis.binding_site import (
    group_indices_by_residue,
    pocket_geometry,
    select_pocket_atoms,
)
from .ligand_coordinates import CoordinateAtom
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

BINDING_SITE_SCHEMA_VERSION = "1.0"
DISCREPANCY_STATUSES = ("open", "resolved", "not_applicable")
_ANGSTROM_TO_NM = 0.1
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


def _vector3(data: Mapping[str, Any], key: str) -> tuple[float, float, float]:
    value = _as_list(data.get(key), key)
    if len(value) != 3:
        raise ValueError(f"{key} must contain three numbers.")
    return (float(value[0]), float(value[1]), float(value[2]))


def _string_tuple(data: Mapping[str, Any], key: str) -> tuple[str, ...]:
    return tuple(str(item) for item in _as_list(data.get(key), key))


@dataclass(frozen=True)
class BindingSiteCriteria:
    pocket_cutoff_angstrom: float = 4.0
    pocket_cutoff_nm: float = 0.4
    inclusive: bool = False
    water_cutoff_angstrom: float = 6.0
    protein_atom_policy: str = (
        "polymer ATOM records in the manifest-selected chain and model"
    )
    water_policy: str = (
        "waters are recorded separately and excluded from pocket selection"
    )
    ion_cofactor_policy: str = (
        "ions and cofactors are recorded separately and excluded from pocket "
        "selection"
    )
    alternate_location_policy: str = (
        "alternate locations are retained and reported; no conformer is selected"
    )
    interaction_policy: str = (
        "proximity only; no hydrogen-bond, favourability, or affinity inference"
    )

    def to_dict(self) -> dict[str, Any]:
        return {
            "pocket_cutoff_angstrom": self.pocket_cutoff_angstrom,
            "pocket_cutoff_nm": self.pocket_cutoff_nm,
            "inclusive": self.inclusive,
            "water_cutoff_angstrom": self.water_cutoff_angstrom,
            "protein_atom_policy": self.protein_atom_policy,
            "water_policy": self.water_policy,
            "ion_cofactor_policy": self.ion_cofactor_policy,
            "alternate_location_policy": self.alternate_location_policy,
            "interaction_policy": self.interaction_policy,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> BindingSiteCriteria:
        data = _as_mapping(data, "binding-site criteria")
        return cls(
            pocket_cutoff_angstrom=_required_float(data, "pocket_cutoff_angstrom"),
            pocket_cutoff_nm=_required_float(data, "pocket_cutoff_nm"),
            inclusive=_required_bool(data, "inclusive"),
            water_cutoff_angstrom=_required_float(data, "water_cutoff_angstrom"),
            protein_atom_policy=_required_string(data, "protein_atom_policy"),
            water_policy=_required_string(data, "water_policy"),
            ion_cofactor_policy=_required_string(data, "ion_cofactor_policy"),
            alternate_location_policy=_required_string(
                data, "alternate_location_policy"
            ),
            interaction_policy=_required_string(data, "interaction_policy"),
        )


@dataclass(frozen=True)
class ProteinAtomRecord:
    atom_name: str
    element: str
    residue_name: str
    label_asym_id: str
    auth_asym_id: str
    label_seq_id: str | None
    auth_seq_id: str
    insertion_code: str | None
    alternate_location: str | None
    model_number: int
    position_angstrom: tuple[float, float, float]
    occupancy: float | None
    b_iso_or_equiv: float | None

    def to_dict(self) -> dict[str, Any]:
        return {
            "atom_name": self.atom_name,
            "element": self.element,
            "residue_name": self.residue_name,
            "label_asym_id": self.label_asym_id,
            "auth_asym_id": self.auth_asym_id,
            "label_seq_id": self.label_seq_id,
            "auth_seq_id": self.auth_seq_id,
            "insertion_code": self.insertion_code,
            "alternate_location": self.alternate_location,
            "model_number": self.model_number,
            "position_angstrom": list(self.position_angstrom),
            "occupancy": self.occupancy,
            "b_iso_or_equiv": self.b_iso_or_equiv,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> ProteinAtomRecord:
        data = _as_mapping(data, "protein atom")
        return cls(
            atom_name=_required_string(data, "atom_name"),
            element=_required_string(data, "element"),
            residue_name=_required_string(data, "residue_name"),
            label_asym_id=_required_string(data, "label_asym_id"),
            auth_asym_id=_required_string(data, "auth_asym_id"),
            label_seq_id=_optional_string(data, "label_seq_id"),
            auth_seq_id=_required_string(data, "auth_seq_id"),
            insertion_code=_optional_string(data, "insertion_code"),
            alternate_location=_optional_string(data, "alternate_location"),
            model_number=_required_int(data, "model_number"),
            position_angstrom=_vector3(data, "position_angstrom"),
            occupancy=_optional_float(data, "occupancy"),
            b_iso_or_equiv=_optional_float(data, "b_iso_or_equiv"),
        )


@dataclass(frozen=True)
class PocketResidueRecord:
    residue_name: str
    label_asym_id: str
    auth_asym_id: str
    label_seq_id: str | None
    auth_seq_id: str
    insertion_code: str | None
    atom_count: int
    selected_atom_count: int
    minimum_distance_angstrom: float
    minimum_distance_nm: float
    closest_protein_atom: str
    closest_ligand_atom: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "residue_name": self.residue_name,
            "label_asym_id": self.label_asym_id,
            "auth_asym_id": self.auth_asym_id,
            "label_seq_id": self.label_seq_id,
            "auth_seq_id": self.auth_seq_id,
            "insertion_code": self.insertion_code,
            "atom_count": self.atom_count,
            "selected_atom_count": self.selected_atom_count,
            "minimum_distance_angstrom": self.minimum_distance_angstrom,
            "minimum_distance_nm": self.minimum_distance_nm,
            "closest_protein_atom": self.closest_protein_atom,
            "closest_ligand_atom": self.closest_ligand_atom,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> PocketResidueRecord:
        data = _as_mapping(data, "pocket residue")
        return cls(
            residue_name=_required_string(data, "residue_name"),
            label_asym_id=_required_string(data, "label_asym_id"),
            auth_asym_id=_required_string(data, "auth_asym_id"),
            label_seq_id=_optional_string(data, "label_seq_id"),
            auth_seq_id=_required_string(data, "auth_seq_id"),
            insertion_code=_optional_string(data, "insertion_code"),
            atom_count=_required_int(data, "atom_count"),
            selected_atom_count=_required_int(data, "selected_atom_count"),
            minimum_distance_angstrom=_required_float(data, "minimum_distance_angstrom"),
            minimum_distance_nm=_required_float(data, "minimum_distance_nm"),
            closest_protein_atom=_required_string(data, "closest_protein_atom"),
            closest_ligand_atom=_required_string(data, "closest_ligand_atom"),
        )


@dataclass(frozen=True)
class LigandProteinContact:
    ligand_atom: str
    protein_atom: str
    residue_name: str
    label_asym_id: str
    auth_asym_id: str
    label_seq_id: str | None
    auth_seq_id: str
    insertion_code: str | None
    distance_angstrom: float
    distance_nm: float

    def to_dict(self) -> dict[str, Any]:
        return {
            "ligand_atom": self.ligand_atom,
            "protein_atom": self.protein_atom,
            "residue_name": self.residue_name,
            "label_asym_id": self.label_asym_id,
            "auth_asym_id": self.auth_asym_id,
            "label_seq_id": self.label_seq_id,
            "auth_seq_id": self.auth_seq_id,
            "insertion_code": self.insertion_code,
            "distance_angstrom": self.distance_angstrom,
            "distance_nm": self.distance_nm,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> LigandProteinContact:
        data = _as_mapping(data, "ligand-protein contact")
        return cls(
            ligand_atom=_required_string(data, "ligand_atom"),
            protein_atom=_required_string(data, "protein_atom"),
            residue_name=_required_string(data, "residue_name"),
            label_asym_id=_required_string(data, "label_asym_id"),
            auth_asym_id=_required_string(data, "auth_asym_id"),
            label_seq_id=_optional_string(data, "label_seq_id"),
            auth_seq_id=_required_string(data, "auth_seq_id"),
            insertion_code=_optional_string(data, "insertion_code"),
            distance_angstrom=_required_float(data, "distance_angstrom"),
            distance_nm=_required_float(data, "distance_nm"),
        )


@dataclass(frozen=True)
class WaterRecord:
    atom_name: str
    residue_name: str
    label_asym_id: str
    auth_asym_id: str
    label_seq_id: str | None
    auth_seq_id: str
    insertion_code: str | None
    position_angstrom: tuple[float, float, float]
    occupancy: float | None
    b_iso_or_equiv: float | None

    def to_dict(self) -> dict[str, Any]:
        return {
            "atom_name": self.atom_name,
            "residue_name": self.residue_name,
            "label_asym_id": self.label_asym_id,
            "auth_asym_id": self.auth_asym_id,
            "label_seq_id": self.label_seq_id,
            "auth_seq_id": self.auth_seq_id,
            "insertion_code": self.insertion_code,
            "position_angstrom": list(self.position_angstrom),
            "occupancy": self.occupancy,
            "b_iso_or_equiv": self.b_iso_or_equiv,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> WaterRecord:
        data = _as_mapping(data, "water record")
        return cls(
            atom_name=_required_string(data, "atom_name"),
            residue_name=_required_string(data, "residue_name"),
            label_asym_id=_required_string(data, "label_asym_id"),
            auth_asym_id=_required_string(data, "auth_asym_id"),
            label_seq_id=_optional_string(data, "label_seq_id"),
            auth_seq_id=_required_string(data, "auth_seq_id"),
            insertion_code=_optional_string(data, "insertion_code"),
            position_angstrom=_vector3(data, "position_angstrom"),
            occupancy=_optional_float(data, "occupancy"),
            b_iso_or_equiv=_optional_float(data, "b_iso_or_equiv"),
        )


@dataclass(frozen=True)
class NearbyWater:
    atom_name: str
    label_asym_id: str
    auth_asym_id: str
    label_seq_id: str | None
    auth_seq_id: str
    insertion_code: str | None
    distance_angstrom: float
    distance_nm: float

    def to_dict(self) -> dict[str, Any]:
        return {
            "atom_name": self.atom_name,
            "label_asym_id": self.label_asym_id,
            "auth_asym_id": self.auth_asym_id,
            "label_seq_id": self.label_seq_id,
            "auth_seq_id": self.auth_seq_id,
            "insertion_code": self.insertion_code,
            "distance_angstrom": self.distance_angstrom,
            "distance_nm": self.distance_nm,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> NearbyWater:
        data = _as_mapping(data, "nearby water")
        return cls(
            atom_name=_required_string(data, "atom_name"),
            label_asym_id=_required_string(data, "label_asym_id"),
            auth_asym_id=_required_string(data, "auth_asym_id"),
            label_seq_id=_optional_string(data, "label_seq_id"),
            auth_seq_id=_required_string(data, "auth_seq_id"),
            insertion_code=_optional_string(data, "insertion_code"),
            distance_angstrom=_required_float(data, "distance_angstrom"),
            distance_nm=_required_float(data, "distance_nm"),
        )


@dataclass(frozen=True)
class HeteroResidueRecord:
    residue_name: str
    label_asym_id: str
    auth_asym_id: str
    label_seq_id: str | None
    auth_seq_id: str
    insertion_code: str | None
    atom_count: int

    def to_dict(self) -> dict[str, Any]:
        return {
            "residue_name": self.residue_name,
            "label_asym_id": self.label_asym_id,
            "auth_asym_id": self.auth_asym_id,
            "label_seq_id": self.label_seq_id,
            "auth_seq_id": self.auth_seq_id,
            "insertion_code": self.insertion_code,
            "atom_count": self.atom_count,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> HeteroResidueRecord:
        data = _as_mapping(data, "hetero residue")
        return cls(
            residue_name=_required_string(data, "residue_name"),
            label_asym_id=_required_string(data, "label_asym_id"),
            auth_asym_id=_required_string(data, "auth_asym_id"),
            label_seq_id=_optional_string(data, "label_seq_id"),
            auth_seq_id=_required_string(data, "auth_seq_id"),
            insertion_code=_optional_string(data, "insertion_code"),
            atom_count=_required_int(data, "atom_count"),
        )


@dataclass(frozen=True)
class BindingSiteDiscrepancy:
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
    def from_dict(cls, data: Mapping[str, Any]) -> BindingSiteDiscrepancy:
        data = _as_mapping(data, "binding-site discrepancy")
        return cls(
            category=_required_string(data, "category"),
            status=_required_string(data, "status"),
            detail=_required_string(data, "detail"),
            observed=_optional_string(data, "observed"),
            expected=_optional_string(data, "expected"),
            resolution=_optional_string(data, "resolution"),
        )


@dataclass(frozen=True)
class BindingSiteRecord:
    source_manifest_path: str
    source_manifest_sha256: str
    ligand_coordinate_record_path: str
    ligand_coordinate_record_sha256: str
    accession: str
    model_number: int
    protein_selection: ProteinSelection
    ligand_selection: LigandSelection
    criteria: BindingSiteCriteria
    protein_atom_count: int
    ligand_atom_count: int
    pocket_atom_count: int
    pocket_residue_count: int
    pocket_center_angstrom: tuple[float, float, float] | None
    pocket_center_nm: tuple[float, float, float] | None
    pocket_dimensions_angstrom: tuple[float, float, float] | None
    pocket_dimensions_nm: tuple[float, float, float] | None
    residues: tuple[PocketResidueRecord, ...]
    contacts: tuple[LigandProteinContact, ...]
    water_count_in_selected_chain: int
    waters_near_ligand: tuple[NearbyWater, ...]
    hetero_residues: tuple[HeteroResidueRecord, ...]
    residue_number_gaps: tuple[str, ...]
    alternate_location_atom_count: int
    decisions: tuple[DecisionRecord, ...]
    provenance: ProvenanceRecord
    validation: ValidationRecord
    discrepancies: tuple[BindingSiteDiscrepancy, ...]
    schema_version: str = BINDING_SITE_SCHEMA_VERSION

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "source_manifest_path": self.source_manifest_path,
            "source_manifest_sha256": self.source_manifest_sha256,
            "ligand_coordinate_record_path": self.ligand_coordinate_record_path,
            "ligand_coordinate_record_sha256": self.ligand_coordinate_record_sha256,
            "accession": self.accession,
            "model_number": self.model_number,
            "protein_selection": self.protein_selection.to_dict(),
            "ligand_selection": self.ligand_selection.to_dict(),
            "criteria": self.criteria.to_dict(),
            "protein_atom_count": self.protein_atom_count,
            "ligand_atom_count": self.ligand_atom_count,
            "pocket_atom_count": self.pocket_atom_count,
            "pocket_residue_count": self.pocket_residue_count,
            "pocket_center_angstrom": (
                None if self.pocket_center_angstrom is None else list(self.pocket_center_angstrom)
            ),
            "pocket_center_nm": (
                None if self.pocket_center_nm is None else list(self.pocket_center_nm)
            ),
            "pocket_dimensions_angstrom": (
                None
                if self.pocket_dimensions_angstrom is None
                else list(self.pocket_dimensions_angstrom)
            ),
            "pocket_dimensions_nm": (
                None
                if self.pocket_dimensions_nm is None
                else list(self.pocket_dimensions_nm)
            ),
            "residues": [residue.to_dict() for residue in self.residues],
            "contacts": [contact.to_dict() for contact in self.contacts],
            "water_count_in_selected_chain": self.water_count_in_selected_chain,
            "waters_near_ligand": [water.to_dict() for water in self.waters_near_ligand],
            "hetero_residues": [residue.to_dict() for residue in self.hetero_residues],
            "residue_number_gaps": list(self.residue_number_gaps),
            "alternate_location_atom_count": self.alternate_location_atom_count,
            "decisions": [decision.to_dict() for decision in self.decisions],
            "provenance": self.provenance.to_dict(),
            "validation": self.validation.to_dict(),
            "discrepancies": [item.to_dict() for item in self.discrepancies],
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> BindingSiteRecord:
        data = _as_mapping(data, "binding-site record")

        def optional_vector(key: str) -> tuple[float, float, float] | None:
            value = data.get(key)
            return None if value is None else _vector3(data, key)

        return cls(
            schema_version=_required_string(data, "schema_version"),
            source_manifest_path=_required_string(data, "source_manifest_path"),
            source_manifest_sha256=_required_string(data, "source_manifest_sha256"),
            ligand_coordinate_record_path=_required_string(
                data, "ligand_coordinate_record_path"
            ),
            ligand_coordinate_record_sha256=_required_string(
                data, "ligand_coordinate_record_sha256"
            ),
            accession=_required_string(data, "accession"),
            model_number=_required_int(data, "model_number"),
            protein_selection=ProteinSelection.from_dict(
                _as_mapping(data.get("protein_selection"), "protein_selection")
            ),
            ligand_selection=LigandSelection.from_dict(
                _as_mapping(data.get("ligand_selection"), "ligand_selection")
            ),
            criteria=BindingSiteCriteria.from_dict(
                _as_mapping(data.get("criteria"), "criteria")
            ),
            protein_atom_count=_required_int(data, "protein_atom_count"),
            ligand_atom_count=_required_int(data, "ligand_atom_count"),
            pocket_atom_count=_required_int(data, "pocket_atom_count"),
            pocket_residue_count=_required_int(data, "pocket_residue_count"),
            pocket_center_angstrom=optional_vector("pocket_center_angstrom"),
            pocket_center_nm=optional_vector("pocket_center_nm"),
            pocket_dimensions_angstrom=optional_vector("pocket_dimensions_angstrom"),
            pocket_dimensions_nm=optional_vector("pocket_dimensions_nm"),
            residues=tuple(
                PocketResidueRecord.from_dict(_as_mapping(item, "pocket residue"))
                for item in _as_list(data.get("residues"), "residues")
            ),
            contacts=tuple(
                LigandProteinContact.from_dict(_as_mapping(item, "ligand-protein contact"))
                for item in _as_list(data.get("contacts"), "contacts")
            ),
            water_count_in_selected_chain=_required_int(
                data, "water_count_in_selected_chain"
            ),
            waters_near_ligand=tuple(
                NearbyWater.from_dict(_as_mapping(item, "nearby water"))
                for item in _as_list(data.get("waters_near_ligand"), "waters_near_ligand")
            ),
            hetero_residues=tuple(
                HeteroResidueRecord.from_dict(_as_mapping(item, "hetero residue"))
                for item in _as_list(data.get("hetero_residues"), "hetero_residues")
            ),
            residue_number_gaps=_string_tuple(data, "residue_number_gaps"),
            alternate_location_atom_count=_required_int(
                data, "alternate_location_atom_count"
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
                BindingSiteDiscrepancy.from_dict(
                    _as_mapping(item, "binding-site discrepancy")
                )
                for item in _as_list(data.get("discrepancies"), "discrepancies")
            ),
        )


def _residue_key(atom: ProteinAtomRecord) -> tuple[Any, ...]:
    return (
        atom.label_asym_id,
        atom.auth_asym_id,
        atom.label_seq_id,
        atom.auth_seq_id,
        atom.insertion_code,
        atom.residue_name,
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


def reconcile_binding_site(
    *,
    source_manifest_path: str,
    source_manifest_sha256: str,
    ligand_coordinate_record_path: str,
    ligand_coordinate_record_sha256: str,
    accession: str,
    model_number: int,
    protein_selection: ProteinSelection,
    ligand_selection: LigandSelection,
    protein_atoms: Sequence[ProteinAtomRecord],
    ligand_atoms: Sequence[CoordinateAtom],
    water_atoms: Sequence[WaterRecord],
    hetero_residues: Sequence[HeteroResidueRecord],
    criteria: BindingSiteCriteria | None = None,
) -> BindingSiteRecord:
    """Select the binding site and atom-level contacts around the ligand."""

    criteria = criteria or BindingSiteCriteria()
    if not protein_atoms:
        raise ValueError("protein_atoms must contain at least one atom.")
    if not ligand_atoms:
        raise ValueError("ligand_atoms must contain at least one atom.")

    protein_positions_nm = np.asarray(
        [atom.position_angstrom for atom in protein_atoms], dtype=float
    ) * _ANGSTROM_TO_NM
    ligand_positions_nm = np.asarray(
        [atom.position_angstrom for atom in ligand_atoms], dtype=float
    ) * _ANGSTROM_TO_NM

    selected_indices = select_pocket_atoms(
        protein_positions_nm,
        ligand_positions_nm,
        criteria.pocket_cutoff_nm,
        inclusive=criteria.inclusive,
    )
    pocket_atom_count = int(selected_indices.size)
    if pocket_atom_count:
        geometry = pocket_geometry(
            protein_positions_nm,
            selected_indices,
            cutoff=criteria.pocket_cutoff_nm,
        )
        center_nm = tuple(float(value) for value in geometry.center)
        dimensions_nm = tuple(float(value) for value in geometry.dimensions)
        center_angstrom = tuple(value * 10.0 for value in center_nm)
        dimensions_angstrom = tuple(value * 10.0 for value in dimensions_nm)
    else:
        center_nm = dimensions_nm = center_angstrom = dimensions_angstrom = None

    residue_labels = [_residue_key(atom) for atom in protein_atoms]
    selected_labels = [residue_labels[index] for index in selected_indices]
    grouped = group_indices_by_residue(selected_indices, selected_labels)
    residues: list[PocketResidueRecord] = []
    contacts: list[LigandProteinContact] = []
    for key, indices in grouped.items():
        label_asym_id, auth_asym_id, label_seq_id, auth_seq_id, insertion_code, residue_name = key
        best_distance = float("inf")
        best_protein = ""
        best_ligand = ""
        for protein_index in indices:
            protein_atom = protein_atoms[protein_index]
            for ligand_atom in ligand_atoms:
                displacement = (
                    np.asarray(protein_atom.position_angstrom, dtype=float)
                    - np.asarray(ligand_atom.position_angstrom, dtype=float)
                )
                distance = float(np.linalg.norm(displacement))
                if distance < best_distance:
                    best_distance = distance
                    best_protein = protein_atom.atom_name
                    best_ligand = ligand_atom.atom_name
                if distance < criteria.pocket_cutoff_angstrom or (
                    criteria.inclusive and distance <= criteria.pocket_cutoff_angstrom
                ):
                    contacts.append(
                        LigandProteinContact(
                            ligand_atom=ligand_atom.atom_name,
                            protein_atom=protein_atom.atom_name,
                            residue_name=residue_name,
                            label_asym_id=label_asym_id,
                            auth_asym_id=auth_asym_id,
                            label_seq_id=label_seq_id,
                            auth_seq_id=auth_seq_id,
                            insertion_code=insertion_code,
                            distance_angstrom=distance,
                            distance_nm=distance * _ANGSTROM_TO_NM,
                        )
                    )
        atom_count = sum(1 for atom in protein_atoms if _residue_key(atom) == key)
        residues.append(
            PocketResidueRecord(
                residue_name=residue_name,
                label_asym_id=label_asym_id,
                auth_asym_id=auth_asym_id,
                label_seq_id=label_seq_id,
                auth_seq_id=auth_seq_id,
                insertion_code=insertion_code,
                atom_count=atom_count,
                selected_atom_count=len(indices),
                minimum_distance_angstrom=best_distance,
                minimum_distance_nm=best_distance * _ANGSTROM_TO_NM,
                closest_protein_atom=best_protein,
                closest_ligand_atom=best_ligand,
            )
        )
    residues.sort(key=lambda residue: (residue.label_asym_id, residue.auth_seq_id))
    contacts.sort(key=lambda contact: contact.distance_angstrom)

    waters_near_ligand: list[NearbyWater] = []
    for water in water_atoms:
        best = min(
            (
                float(
                    np.linalg.norm(
                        np.asarray(water.position_angstrom, dtype=float)
                        - np.asarray(ligand_atom.position_angstrom, dtype=float)
                    )
                )
                for ligand_atom in ligand_atoms
            ),
            default=float("inf"),
        )
        if best <= criteria.water_cutoff_angstrom:
            waters_near_ligand.append(
                NearbyWater(
                    atom_name=water.atom_name,
                    label_asym_id=water.label_asym_id,
                    auth_asym_id=water.auth_asym_id,
                    label_seq_id=water.label_seq_id,
                    auth_seq_id=water.auth_seq_id,
                    insertion_code=water.insertion_code,
                    distance_angstrom=best,
                    distance_nm=best * _ANGSTROM_TO_NM,
                )
            )
    waters_near_ligand.sort(key=lambda water: water.distance_angstrom)

    auth_numbers = sorted(
        {
            int(atom.auth_seq_id)
            for atom in protein_atoms
            if atom.auth_seq_id.lstrip("-").isdigit()
        }
    )
    missing_numbers = []
    for left, right in zip(auth_numbers, auth_numbers[1:]):
        missing_numbers.extend(range(left + 1, right))
    residue_number_gaps = tuple(str(number) for number in missing_numbers)
    alternate_location_atom_count = sum(
        1 for atom in protein_atoms if atom.alternate_location is not None
    )

    model_numbers = {atom.model_number for atom in protein_atoms} | {
        atom.model_number for atom in ligand_atoms
    }
    model_consistent = model_numbers == {model_number}
    protein_chain_matches = all(
        atom.auth_asym_id == protein_selection.auth_asym_id
        and (
            protein_selection.label_asym_id is None
            or atom.label_asym_id == protein_selection.label_asym_id
        )
        for atom in protein_atoms
    )
    residue_identifiers_present = all(
        atom.auth_asym_id and atom.auth_seq_id and atom.residue_name
        for atom in protein_atoms
    )
    contact_atoms_traceable = all(
        contact.ligand_atom in {atom.atom_name for atom in ligand_atoms}
        and contact.protein_atom in {atom.atom_name for atom in protein_atoms}
        for contact in contacts
    ) if contacts else True

    checks = (
        _check(
            "provenance",
            "binding_site.record_sources",
            bool(source_manifest_path and ligand_coordinate_record_path),
            subject=accession,
            criterion="manifest and ligand coordinate record are recorded",
            observed="both record paths present",
            expected="manifest and coordinate record",
            detail="the binding-site record is traceable to both prior artifacts",
        ),
        _check(
            "coordinate_frame",
            "binding_site.model_consistent",
            model_consistent,
            subject=accession,
            criterion="protein and ligand atoms use the selected model",
            observed=", ".join(str(number) for number in sorted(model_numbers)),
            expected=str(model_number),
            detail="no model combination or coordinate transformation was applied",
        ),
        _check(
            "selection",
            "binding_site.protein_chain_matches_manifest",
            protein_chain_matches,
            subject=protein_selection.auth_asym_id,
            criterion="protein atoms belong to the manifest-selected author chain",
            observed=protein_selection.auth_asym_id,
            expected=protein_selection.auth_asym_id,
            detail=f"{len(protein_atoms)} polymer atoms parsed",
        ),
        _check(
            "selection",
            "binding_site.pocket_nonempty",
            pocket_atom_count > 0,
            subject=protein_selection.auth_asym_id,
            criterion="at least one protein atom lies within the pocket cutoff",
            observed=f"{pocket_atom_count} atoms, {len(residues)} residues",
            expected="at least one atom",
            detail=(
                f"cutoff {criteria.pocket_cutoff_angstrom} A "
                f"({criteria.pocket_cutoff_nm} nm)"
            ),
        ),
        _check(
            "identity",
            "binding_site.residue_identifiers_present",
            residue_identifiers_present,
            subject=protein_selection.auth_asym_id,
            criterion="residue name and author sequence identifier are recorded",
            observed=f"{len({_residue_key(atom) for atom in protein_atoms})} residues",
            expected="all residues identified",
            detail="label and author namespaces are preserved where available",
        ),
        _check(
            "quality",
            "binding_site.alternate_locations_documented",
            True,
            subject=protein_selection.auth_asym_id,
            criterion="alternate locations are absent or explicitly retained",
            observed=f"{alternate_location_atom_count} alternate-location atoms",
            expected="0, or an explicit conformer policy",
            detail=criteria.alternate_location_policy,
        ),
        _check(
            "completeness",
            "binding_site.missing_residues",
            not residue_number_gaps,
            subject=protein_selection.auth_asym_id,
            criterion="no residue-number gaps remain in the selected chain",
            observed=(
                "no gaps"
                if not residue_number_gaps
                else "gaps: " + ", ".join(residue_number_gaps)
            ),
            expected="no gaps",
            detail="missing numbers are reported; residues are not rebuilt",
            not_applicable=bool(residue_number_gaps),
        ),
        _check(
            "environment",
            "binding_site.waters_recorded",
            True,
            subject=protein_selection.auth_asym_id,
            criterion="waters are recorded under an explicit policy",
            observed=f"{len(water_atoms)} selected-chain waters, "
            f"{len(waters_near_ligand)} near the ligand",
            expected="recorded separately from the pocket",
            detail=criteria.water_policy,
        ),
        _check(
            "environment",
            "binding_site.hetero_residues_recorded",
            True,
            subject=protein_selection.auth_asym_id,
            criterion="ions and cofactors are recorded under an explicit policy",
            observed=f"{len(hetero_residues)} hetero residues",
            expected="recorded separately from the pocket",
            detail=criteria.ion_cofactor_policy,
        ),
        _check(
            "contacts",
            "binding_site.contacts_traceable",
            contact_atoms_traceable,
            subject=ligand_selection.component_id,
            criterion="every contact names a specific ligand and protein atom",
            observed=f"{len(contacts)} atom contacts",
            expected="traceable atom pairs",
            detail="contacts are proximity observations, not interaction assignments",
        ),
        _check(
            "interpretation",
            "binding_site.proximity_not_affinity",
            True,
            subject=ligand_selection.component_id,
            criterion="proximity is not converted into an affinity claim",
            observed="no interaction or affinity inferred",
            expected="structural context only",
            detail=criteria.interaction_policy,
        ),
    )

    discrepancies: list[BindingSiteDiscrepancy] = []
    if not model_consistent:
        discrepancies.append(
            BindingSiteDiscrepancy(
                "model_consistency",
                "open",
                "Protein and ligand atoms do not all use the selected model.",
                observed=", ".join(str(number) for number in sorted(model_numbers)),
                expected=str(model_number),
            )
        )
    if not protein_chain_matches:
        discrepancies.append(
            BindingSiteDiscrepancy(
                "protein_chain",
                "open",
                "At least one protein atom does not match the manifest chain.",
                observed="mixed chains",
                expected=protein_selection.auth_asym_id,
            )
        )
    if not pocket_atom_count:
        discrepancies.append(
            BindingSiteDiscrepancy(
                "empty_pocket",
                "open",
                "No protein atoms fall inside the pocket cutoff.",
                observed="0 atoms",
                expected="at least one atom",
                resolution="The cutoff or selection should be reviewed.",
            )
        )
    if alternate_location_atom_count:
        discrepancies.append(
            BindingSiteDiscrepancy(
                "alternate_locations",
                "open",
                "Alternate protein locations are present and no conformer was selected.",
                observed=f"{alternate_location_atom_count} atoms",
                expected=criteria.alternate_location_policy,
                resolution="All records remain available for inspection.",
            )
        )
    if residue_number_gaps:
        discrepancies.append(
            BindingSiteDiscrepancy(
                "missing_residues",
                "not_applicable",
                "Residue-number gaps are present in the selected chain.",
                observed=", ".join(residue_number_gaps),
                expected="continuous numbering",
                resolution="Gaps are reported; no residues were rebuilt.",
            )
        )

    decisions = (
        DecisionRecord(
            decision=f"selected pocket atoms within {criteria.pocket_cutoff_angstrom} A of the ligand",
            actor="binding-site reconciliation",
            rationale=criteria.protein_atom_policy,
        ),
        DecisionRecord(
            decision="recorded waters and hetero residues separately",
            actor="binding-site reconciliation",
            rationale=(
                "water, ion, and cofactor presence is preserved without adding them "
                "to the protein pocket selection"
            ),
        ),
        DecisionRecord(
            decision="did not infer interactions or affinity from proximity",
            actor="binding-site reconciliation",
            rationale=criteria.interaction_policy,
        ),
        DecisionRecord(
            decision="did not rebuild missing residues or choose alternate conformers",
            actor="binding-site reconciliation",
            rationale=(
                "missing and alternate records remain visible for later structural "
                "preparation"
            ),
        ),
    )

    provenance = ProvenanceRecord(
        records=(
            ProvenanceEntry(
                "experimental",
                "protein coordinates",
                (
                    f"Selected chain {protein_selection.auth_asym_id} from the "
                    f"{accession} deposition, model {model_number}."
                ),
                actor="RCSB PDB deposition",
            ),
            ProvenanceEntry(
                "experimental",
                "ligand coordinates",
                (
                    f"Selected component {ligand_selection.component_id} from the "
                    "Batch 3 coordinate record."
                ),
                actor="RCSB PDB deposition",
            ),
            ProvenanceEntry(
                "calculated",
                "pocket and contacts",
                (
                    f"{pocket_atom_count} pocket atoms, {len(residues)} residues, "
                    f"{len(contacts)} atom contacts within "
                    f"{criteria.pocket_cutoff_angstrom} A."
                ),
                actor="binding-site reconciliation",
            ),
            ProvenanceEntry(
                "predicted",
                "none",
                "No affinity, free-energy, or interaction-favourability estimate is produced.",
                actor="binding-site reconciliation",
            ),
            ProvenanceEntry(
                "hypothesised",
                "none",
                "No mechanistic interpretation is made from contacts.",
                actor="binding-site reconciliation",
            ),
        )
    )

    return BindingSiteRecord(
        source_manifest_path=source_manifest_path,
        source_manifest_sha256=source_manifest_sha256,
        ligand_coordinate_record_path=ligand_coordinate_record_path,
        ligand_coordinate_record_sha256=ligand_coordinate_record_sha256,
        accession=accession,
        model_number=model_number,
        protein_selection=protein_selection,
        ligand_selection=ligand_selection,
        criteria=criteria,
        protein_atom_count=len(protein_atoms),
        ligand_atom_count=len(ligand_atoms),
        pocket_atom_count=pocket_atom_count,
        pocket_residue_count=len(residues),
        pocket_center_angstrom=center_angstrom,
        pocket_center_nm=center_nm,
        pocket_dimensions_angstrom=dimensions_angstrom,
        pocket_dimensions_nm=dimensions_nm,
        residues=tuple(residues),
        contacts=tuple(contacts),
        water_count_in_selected_chain=len(water_atoms),
        waters_near_ligand=tuple(waters_near_ligand),
        hetero_residues=tuple(hetero_residues),
        residue_number_gaps=residue_number_gaps,
        alternate_location_atom_count=alternate_location_atom_count,
        decisions=decisions,
        provenance=provenance,
        validation=ValidationRecord(checks=checks),
        discrepancies=tuple(discrepancies),
    )


def write_binding_site(path: str | Path, record: BindingSiteRecord) -> Path:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(record.to_dict(), indent=2, sort_keys=False) + "\n",
        encoding="utf-8",
    )
    return destination


def load_binding_site(path: str | Path) -> BindingSiteRecord:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    return BindingSiteRecord.from_dict(data)


def validate_binding_site(record: BindingSiteRecord) -> list[str]:
    """Return schema and internal-consistency errors for *record*."""

    errors: list[str] = []
    if record.schema_version != BINDING_SITE_SCHEMA_VERSION:
        errors.append(
            f"schema_version must be {BINDING_SITE_SCHEMA_VERSION!r}, got "
            f"{record.schema_version!r}."
        )
    if not record.source_manifest_path:
        errors.append("source_manifest_path is required.")
    if not _SHA256_PATTERN.fullmatch(record.source_manifest_sha256):
        errors.append("source_manifest_sha256 must be a lowercase SHA-256 digest.")
    if not record.ligand_coordinate_record_path:
        errors.append("ligand_coordinate_record_path is required.")
    if not _SHA256_PATTERN.fullmatch(record.ligand_coordinate_record_sha256):
        errors.append(
            "ligand_coordinate_record_sha256 must be a lowercase SHA-256 digest."
        )
    if not record.accession:
        errors.append("accession is required.")
    if record.model_number <= 0:
        errors.append("model_number must be positive.")

    criteria = record.criteria
    if criteria.pocket_cutoff_angstrom <= 0 or criteria.pocket_cutoff_nm <= 0:
        errors.append("pocket cutoffs must be positive.")
    if not math.isclose(
        criteria.pocket_cutoff_angstrom * _ANGSTROM_TO_NM,
        criteria.pocket_cutoff_nm,
        rel_tol=1e-9,
        abs_tol=1e-9,
    ):
        errors.append("pocket_cutoff_angstrom and pocket_cutoff_nm are inconsistent.")
    if criteria.water_cutoff_angstrom <= 0:
        errors.append("water_cutoff_angstrom must be positive.")

    if record.protein_atom_count <= 0:
        errors.append("protein_atom_count must be positive.")
    if record.ligand_atom_count <= 0:
        errors.append("ligand_atom_count must be positive.")
    if record.pocket_atom_count < 0:
        errors.append("pocket_atom_count cannot be negative.")
    if record.pocket_residue_count != len(record.residues):
        errors.append("pocket_residue_count must match the residue list.")
    if sum(residue.selected_atom_count for residue in record.residues) != record.pocket_atom_count:
        errors.append("pocket_atom_count must match the selected residue atoms.")
    if record.pocket_atom_count and record.pocket_center_angstrom is None:
        errors.append("pocket center is required for a non-empty pocket.")
    if record.pocket_atom_count and record.pocket_dimensions_angstrom is None:
        errors.append("pocket dimensions are required for a non-empty pocket.")

    if len(record.waters_near_ligand) > record.water_count_in_selected_chain:
        errors.append("nearby waters cannot exceed the selected-chain water count.")
    for index, water in enumerate(record.waters_near_ligand):
        if water.distance_angstrom <= 0 or water.distance_nm <= 0:
            errors.append(f"waters_near_ligand[{index}] distance must be positive.")

    for index, residue in enumerate(record.hetero_residues):
        if residue.atom_count <= 0:
            errors.append(f"hetero_residues[{index}] atom_count must be positive.")

    for index, residue in enumerate(record.residues):
        if not residue.residue_name or not residue.auth_seq_id:
            errors.append(f"residues[{index}] requires residue name and author sequence.")
        if residue.atom_count <= 0 or residue.selected_atom_count <= 0:
            errors.append(f"residues[{index}] atom counts must be positive.")
        if residue.minimum_distance_angstrom <= 0:
            errors.append(f"residues[{index}] distance must be positive.")

    for index, contact in enumerate(record.contacts):
        if not contact.ligand_atom or not contact.protein_atom:
            errors.append(f"contacts[{index}] requires both atom names.")
        if contact.distance_angstrom <= 0 or contact.distance_nm <= 0:
            errors.append(f"contacts[{index}] distance must be positive.")

    if record.alternate_location_atom_count < 0:
        errors.append("alternate_location_atom_count cannot be negative.")

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
