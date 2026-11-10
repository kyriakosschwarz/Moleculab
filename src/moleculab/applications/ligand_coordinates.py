"""Coordinate reconciliation for a chemically identified ligand.

This module connects the explicit atom mapping from a ligand identity record to
the deposited coordinates.  It records every coordinate, occupancy,
displacement parameter, alternate location, bond length, angle, aromatic-ring
planarity, and close contact.  It does not add missing atoms, choose among
alternate conformers, average occupancy, assign hydrogens or protonation states,
type force-field atoms, minimise, or dock.

The caller supplies parsed coordinate records and the identity record's bond
graph.  File parsing and chemistry toolkits stay outside this module.
"""

from __future__ import annotations

import json
import math
import re
from collections import Counter
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from .ligand_identity import DepositedBond
from .structure_manifest import (
    PROVENANCE_LABELS,
    DecisionRecord,
    LigandSelection,
    ProvenanceEntry,
    ProvenanceRecord,
    ValidationCheck,
    ValidationRecord,
)

LIGAND_COORDINATE_SCHEMA_VERSION = "1.0"
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
class GeometryCriteria:
    """Broad sanity bounds, not force-field target values."""

    bond_length_min_angstrom: float = 0.5
    bond_length_max_angstrom: float = 2.5
    angle_min_degrees: float = 30.0
    angle_max_degrees: float = 180.0
    ring_planarity_rmsd_max_angstrom: float = 0.2
    close_contact_min_angstrom: float = 1.5

    def to_dict(self) -> dict[str, Any]:
        return {
            "bond_length_min_angstrom": self.bond_length_min_angstrom,
            "bond_length_max_angstrom": self.bond_length_max_angstrom,
            "angle_min_degrees": self.angle_min_degrees,
            "angle_max_degrees": self.angle_max_degrees,
            "ring_planarity_rmsd_max_angstrom": self.ring_planarity_rmsd_max_angstrom,
            "close_contact_min_angstrom": self.close_contact_min_angstrom,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> GeometryCriteria:
        data = _as_mapping(data, "geometry criteria")
        return cls(
            bond_length_min_angstrom=_required_float(
                data, "bond_length_min_angstrom"
            ),
            bond_length_max_angstrom=_required_float(
                data, "bond_length_max_angstrom"
            ),
            angle_min_degrees=_required_float(data, "angle_min_degrees"),
            angle_max_degrees=_required_float(data, "angle_max_degrees"),
            ring_planarity_rmsd_max_angstrom=_required_float(
                data, "ring_planarity_rmsd_max_angstrom"
            ),
            close_contact_min_angstrom=_required_float(
                data, "close_contact_min_angstrom"
            ),
        )


@dataclass(frozen=True)
class CoordinateAtom:
    atom_name: str
    element: str
    reference_index: int | None
    position_angstrom: tuple[float, float, float]
    occupancy: float | None
    b_iso_or_equiv: float | None
    alternate_location: str | None
    model_number: int

    def to_dict(self) -> dict[str, Any]:
        return {
            "atom_name": self.atom_name,
            "element": self.element,
            "reference_index": self.reference_index,
            "position_angstrom": list(self.position_angstrom),
            "occupancy": self.occupancy,
            "b_iso_or_equiv": self.b_iso_or_equiv,
            "alternate_location": self.alternate_location,
            "model_number": self.model_number,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> CoordinateAtom:
        data = _as_mapping(data, "coordinate atom")
        return cls(
            atom_name=_required_string(data, "atom_name"),
            element=_required_string(data, "element"),
            reference_index=_optional_int(data, "reference_index"),
            position_angstrom=_vector3(data, "position_angstrom"),
            occupancy=_optional_float(data, "occupancy"),
            b_iso_or_equiv=_optional_float(data, "b_iso_or_equiv"),
            alternate_location=_optional_string(data, "alternate_location"),
            model_number=_required_int(data, "model_number"),
        )


@dataclass(frozen=True)
class BondGeometry:
    atom_name_1: str
    atom_name_2: str
    reference_index_1: int | None
    reference_index_2: int | None
    length_angstrom: float
    length_nm: float

    def to_dict(self) -> dict[str, Any]:
        return {
            "atom_name_1": self.atom_name_1,
            "atom_name_2": self.atom_name_2,
            "reference_index_1": self.reference_index_1,
            "reference_index_2": self.reference_index_2,
            "length_angstrom": self.length_angstrom,
            "length_nm": self.length_nm,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> BondGeometry:
        data = _as_mapping(data, "bond geometry")
        return cls(
            atom_name_1=_required_string(data, "atom_name_1"),
            atom_name_2=_required_string(data, "atom_name_2"),
            reference_index_1=_optional_int(data, "reference_index_1"),
            reference_index_2=_optional_int(data, "reference_index_2"),
            length_angstrom=_required_float(data, "length_angstrom"),
            length_nm=_required_float(data, "length_nm"),
        )


@dataclass(frozen=True)
class AngleGeometry:
    atom_name_1: str
    atom_name_2: str
    atom_name_3: str
    angle_degrees: float

    def to_dict(self) -> dict[str, Any]:
        return {
            "atom_name_1": self.atom_name_1,
            "atom_name_2": self.atom_name_2,
            "atom_name_3": self.atom_name_3,
            "angle_degrees": self.angle_degrees,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> AngleGeometry:
        data = _as_mapping(data, "angle geometry")
        return cls(
            atom_name_1=_required_string(data, "atom_name_1"),
            atom_name_2=_required_string(data, "atom_name_2"),
            atom_name_3=_required_string(data, "atom_name_3"),
            angle_degrees=_required_float(data, "angle_degrees"),
        )


@dataclass(frozen=True)
class RingPlanarity:
    atom_names: tuple[str, ...]
    aromatic_bond_count: int
    rmsd_angstrom: float
    max_deviation_angstrom: float

    def to_dict(self) -> dict[str, Any]:
        return {
            "atom_names": list(self.atom_names),
            "aromatic_bond_count": self.aromatic_bond_count,
            "rmsd_angstrom": self.rmsd_angstrom,
            "max_deviation_angstrom": self.max_deviation_angstrom,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> RingPlanarity:
        data = _as_mapping(data, "ring planarity")
        return cls(
            atom_names=_string_tuple(data, "atom_names"),
            aromatic_bond_count=_required_int(data, "aromatic_bond_count"),
            rmsd_angstrom=_required_float(data, "rmsd_angstrom"),
            max_deviation_angstrom=_required_float(data, "max_deviation_angstrom"),
        )


@dataclass(frozen=True)
class CloseContact:
    atom_name_1: str
    atom_name_2: str
    distance_angstrom: float
    distance_nm: float

    def to_dict(self) -> dict[str, Any]:
        return {
            "atom_name_1": self.atom_name_1,
            "atom_name_2": self.atom_name_2,
            "distance_angstrom": self.distance_angstrom,
            "distance_nm": self.distance_nm,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> CloseContact:
        data = _as_mapping(data, "close contact")
        return cls(
            atom_name_1=_required_string(data, "atom_name_1"),
            atom_name_2=_required_string(data, "atom_name_2"),
            distance_angstrom=_required_float(data, "distance_angstrom"),
            distance_nm=_required_float(data, "distance_nm"),
        )


@dataclass(frozen=True)
class CoordinateDiscrepancy:
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
    def from_dict(cls, data: Mapping[str, Any]) -> CoordinateDiscrepancy:
        data = _as_mapping(data, "coordinate discrepancy")
        return cls(
            category=_required_string(data, "category"),
            status=_required_string(data, "status"),
            detail=_required_string(data, "detail"),
            observed=_optional_string(data, "observed"),
            expected=_optional_string(data, "expected"),
            resolution=_optional_string(data, "resolution"),
        )


@dataclass(frozen=True)
class LigandCoordinateRecord:
    source_manifest_path: str
    source_manifest_sha256: str
    identity_record_path: str
    identity_record_sha256: str
    accession: str
    selection: LigandSelection
    criteria: GeometryCriteria
    atoms: tuple[CoordinateAtom, ...]
    bonds: tuple[BondGeometry, ...]
    angles: tuple[AngleGeometry, ...]
    rings: tuple[RingPlanarity, ...]
    close_contacts: tuple[CloseContact, ...]
    missing_component_atoms: tuple[str, ...]
    missing_heavy_atoms: tuple[str, ...]
    alternate_location_atoms: tuple[str, ...]
    decisions: tuple[DecisionRecord, ...]
    provenance: ProvenanceRecord
    validation: ValidationRecord
    discrepancies: tuple[CoordinateDiscrepancy, ...]
    schema_version: str = LIGAND_COORDINATE_SCHEMA_VERSION

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "source_manifest_path": self.source_manifest_path,
            "source_manifest_sha256": self.source_manifest_sha256,
            "identity_record_path": self.identity_record_path,
            "identity_record_sha256": self.identity_record_sha256,
            "accession": self.accession,
            "selection": self.selection.to_dict(),
            "criteria": self.criteria.to_dict(),
            "atoms": [atom.to_dict() for atom in self.atoms],
            "bonds": [bond.to_dict() for bond in self.bonds],
            "angles": [angle.to_dict() for angle in self.angles],
            "rings": [ring.to_dict() for ring in self.rings],
            "close_contacts": [contact.to_dict() for contact in self.close_contacts],
            "missing_component_atoms": list(self.missing_component_atoms),
            "missing_heavy_atoms": list(self.missing_heavy_atoms),
            "alternate_location_atoms": list(self.alternate_location_atoms),
            "decisions": [decision.to_dict() for decision in self.decisions],
            "provenance": self.provenance.to_dict(),
            "validation": self.validation.to_dict(),
            "discrepancies": [item.to_dict() for item in self.discrepancies],
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> LigandCoordinateRecord:
        data = _as_mapping(data, "ligand coordinate record")
        return cls(
            schema_version=_required_string(data, "schema_version"),
            source_manifest_path=_required_string(data, "source_manifest_path"),
            source_manifest_sha256=_required_string(data, "source_manifest_sha256"),
            identity_record_path=_required_string(data, "identity_record_path"),
            identity_record_sha256=_required_string(data, "identity_record_sha256"),
            accession=_required_string(data, "accession"),
            selection=LigandSelection.from_dict(_as_mapping(data.get("selection"), "selection")),
            criteria=GeometryCriteria.from_dict(_as_mapping(data.get("criteria"), "criteria")),
            atoms=tuple(
                CoordinateAtom.from_dict(_as_mapping(item, "coordinate atom"))
                for item in _as_list(data.get("atoms"), "atoms")
            ),
            bonds=tuple(
                BondGeometry.from_dict(_as_mapping(item, "bond geometry"))
                for item in _as_list(data.get("bonds"), "bonds")
            ),
            angles=tuple(
                AngleGeometry.from_dict(_as_mapping(item, "angle geometry"))
                for item in _as_list(data.get("angles"), "angles")
            ),
            rings=tuple(
                RingPlanarity.from_dict(_as_mapping(item, "ring planarity"))
                for item in _as_list(data.get("rings"), "rings")
            ),
            close_contacts=tuple(
                CloseContact.from_dict(_as_mapping(item, "close contact"))
                for item in _as_list(data.get("close_contacts"), "close_contacts")
            ),
            missing_component_atoms=_string_tuple(data, "missing_component_atoms"),
            missing_heavy_atoms=_string_tuple(data, "missing_heavy_atoms"),
            alternate_location_atoms=_string_tuple(data, "alternate_location_atoms"),
            decisions=tuple(
                DecisionRecord.from_dict(_as_mapping(item, "decision record"))
                for item in _as_list(data.get("decisions"), "decisions")
            ),
            provenance=ProvenanceRecord.from_dict(_as_mapping(data.get("provenance"), "provenance")),
            validation=ValidationRecord.from_dict(_as_mapping(data.get("validation"), "validation")),
            discrepancies=tuple(
                CoordinateDiscrepancy.from_dict(_as_mapping(item, "coordinate discrepancy"))
                for item in _as_list(data.get("discrepancies"), "discrepancies")
            ),
        )


def _distance(first: np.ndarray, second: np.ndarray) -> float:
    return float(np.linalg.norm(first - second))


def _angle(first: np.ndarray, center: np.ndarray, second: np.ndarray) -> float:
    left = first - center
    right = second - center
    denominator = np.linalg.norm(left) * np.linalg.norm(right)
    if denominator == 0.0:
        return float("nan")
    cosine = float(np.dot(left, right) / denominator)
    return math.degrees(math.acos(max(-1.0, min(1.0, cosine))))


def _connected_components(
    atoms: set[str],
    edges: Sequence[tuple[str, str]],
) -> list[tuple[tuple[str, ...], int]]:
    adjacency: dict[str, set[str]] = {atom: set() for atom in atoms}
    for first, second in edges:
        adjacency[first].add(second)
        adjacency[second].add(first)
    seen: set[str] = set()
    components: list[tuple[tuple[str, ...], int]] = []
    for start in sorted(atoms):
        if start in seen:
            continue
        stack = [start]
        component = []
        while stack:
            atom = stack.pop()
            if atom in seen:
                continue
            seen.add(atom)
            component.append(atom)
            stack.extend(adjacency[atom] - seen)
        component_set = set(component)
        bond_count = sum(
            1
            for first, second in edges
            if first in component_set and second in component_set
        )
        components.append((tuple(sorted(component)), bond_count))
    return components


def _plane_deviation(points: np.ndarray) -> tuple[float, float]:
    centred = points - points.mean(axis=0)
    _, _, vectors = np.linalg.svd(centred, full_matrices=False)
    normal = vectors[-1]
    distances = np.abs(centred @ normal)
    rmsd = float(np.sqrt(np.mean(distances**2)))
    return rmsd, float(np.max(distances))


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


def reconcile_ligand_coordinates(
    *,
    source_manifest_path: str,
    source_manifest_sha256: str,
    identity_record_path: str,
    identity_record_sha256: str,
    accession: str,
    selection: LigandSelection,
    atoms: Sequence[CoordinateAtom],
    bonds: Sequence[DepositedBond],
    expected_atom_names: Sequence[str],
    component_atom_elements: Mapping[str, str],
    reference_stereocenters: Sequence[str],
    deposited_stereochemistry: Sequence[str],
    criteria: GeometryCriteria | None = None,
) -> LigandCoordinateRecord:
    """Reconcile supplied coordinates with the identity record's atom mapping."""

    criteria = criteria or GeometryCriteria()
    expected_names = tuple(dict.fromkeys(str(name) for name in expected_atom_names))
    component_elements = {
        str(name): str(element) for name, element in component_atom_elements.items()
    }
    atom_names = [atom.atom_name for atom in atoms]
    name_counts = Counter(atom_names)
    alternate_location_atoms = tuple(
        sorted(name for name, count in name_counts.items() if count > 1)
    )
    missing_coordinate_atoms = tuple(
        sorted(set(expected_names) - set(atom_names))
    )
    missing_component_atoms = tuple(
        sorted(set(component_elements) - set(atom_names))
    )
    missing_heavy_atoms = tuple(
        name for name in missing_component_atoms if component_elements[name] != "H"
    )
    ambiguous = bool(alternate_location_atoms)

    if ambiguous:
        atom_by_name = {}
    else:
        atom_by_name = {atom.atom_name: atom for atom in atoms}
        positions = {
            name: np.asarray(atom.position_angstrom, dtype=float)
            for name, atom in atom_by_name.items()
        }

    bond_geometry: list[BondGeometry] = []
    bond_lengths: list[float] = []
    if not ambiguous:
        for bond in bonds:
            first = atom_by_name.get(bond.atom_name_1)
            second = atom_by_name.get(bond.atom_name_2)
            if first is None or second is None:
                continue
            length = _distance(
                positions[bond.atom_name_1], positions[bond.atom_name_2]
            )
            bond_lengths.append(length)
            bond_geometry.append(
                BondGeometry(
                    atom_name_1=bond.atom_name_1,
                    atom_name_2=bond.atom_name_2,
                    reference_index_1=first.reference_index,
                    reference_index_2=second.reference_index,
                    length_angstrom=length,
                    length_nm=length * _ANGSTROM_TO_NM,
                )
            )

    adjacency: dict[str, set[str]] = {}
    if not ambiguous:
        for bond in bonds:
            if bond.atom_name_1 in atom_by_name and bond.atom_name_2 in atom_by_name:
                adjacency.setdefault(bond.atom_name_1, set()).add(bond.atom_name_2)
                adjacency.setdefault(bond.atom_name_2, set()).add(bond.atom_name_1)
    angle_geometry: list[AngleGeometry] = []
    angle_values: list[float] = []
    if not ambiguous:
        for center, neighbours in sorted(adjacency.items()):
            ordered = sorted(neighbours)
            for first_index, first in enumerate(ordered):
                for second in ordered[first_index + 1 :]:
                    angle = _angle(
                        positions[first], positions[center], positions[second]
                    )
                    if math.isfinite(angle):
                        angle_values.append(angle)
                        angle_geometry.append(
                            AngleGeometry(
                                atom_name_1=first,
                                atom_name_2=center,
                                atom_name_3=second,
                                angle_degrees=angle,
                            )
                        )

    aromatic_edges = [
        (bond.atom_name_1, bond.atom_name_2)
        for bond in bonds
        if bond.aromatic
    ]
    aromatic_atoms = {name for edge in aromatic_edges for name in edge}
    ring_planarity: list[RingPlanarity] = []
    ring_rmsds: list[float] = []
    if not ambiguous:
        for component_atoms, aromatic_bond_count in _connected_components(
            aromatic_atoms, aromatic_edges
        ):
            if len(component_atoms) < 3:
                continue
            if not all(
                component in atom_by_name for component in component_atoms
            ):
                continue
            degrees = Counter(
                name
                for edge in aromatic_edges
                if edge[0] in component_atoms and edge[1] in component_atoms
                for name in edge
            )
            if any(degrees[name] < 2 for name in component_atoms):
                continue
            points = np.asarray(
                [positions[name] for name in component_atoms], dtype=float
            )
            rmsd, max_deviation = _plane_deviation(points)
            ring_rmsds.append(rmsd)
            ring_planarity.append(
                RingPlanarity(
                    atom_names=component_atoms,
                    aromatic_bond_count=aromatic_bond_count,
                    rmsd_angstrom=rmsd,
                    max_deviation_angstrom=max_deviation,
                )
            )

    close_contacts: list[CloseContact] = []
    if not ambiguous:
        bonded_pairs = {
            frozenset((bond.atom_name_1, bond.atom_name_2)) for bond in bonds
        }
        neighbour_sets = {
            name: set(neighbours) for name, neighbours in adjacency.items()
        }
        names = sorted(atom_by_name)
        for index, first in enumerate(names):
            for second in names[index + 1 :]:
                pair = frozenset((first, second))
                if pair in bonded_pairs:
                    continue
                if neighbour_sets.get(first, set()) & neighbour_sets.get(second, set()):
                    continue
                distance = _distance(positions[first], positions[second])
                if distance < criteria.close_contact_min_angstrom:
                    close_contacts.append(
                        CloseContact(
                            atom_name_1=first,
                            atom_name_2=second,
                            distance_angstrom=distance,
                            distance_nm=distance * _ANGSTROM_TO_NM,
                        )
                    )
        close_contacts.sort(key=lambda contact: contact.distance_angstrom)

    missing_occupancy = tuple(
        sorted(atom.atom_name for atom in atoms if atom.occupancy is None)
    )
    missing_b_factors = tuple(
        sorted(atom.atom_name for atom in atoms if atom.b_iso_or_equiv is None)
    )
    bond_length_out_of_range = [
        value
        for value in bond_lengths
        if not (
            criteria.bond_length_min_angstrom
            <= value
            <= criteria.bond_length_max_angstrom
        )
    ]
    angle_out_of_range = [
        value
        for value in angle_values
        if not (criteria.angle_min_degrees <= value <= criteria.angle_max_degrees)
    ]
    non_planar_rings = [
        ring
        for ring in ring_planarity
        if ring.rmsd_angstrom > criteria.ring_planarity_rmsd_max_angstrom
    ]

    coordinate_names = set(atom_names)
    all_mapped_atoms_have_coordinates = not missing_coordinate_atoms
    alternate_locations_documented = all(
        atom.atom_name not in alternate_location_atoms
        or atom.alternate_location is not None
        for atom in atoms
    )

    checks = (
        _check(
            "provenance",
            "coordinates.source_records",
            bool(source_manifest_path and identity_record_path),
            subject=accession,
            criterion="source manifest and identity record are recorded",
            observed="both record paths present",
            expected="manifest and identity record",
            detail="the coordinate record is traceable to both prior artifacts",
        ),
        _check(
            "coordinate_mapping",
            "coordinates.mapped_atoms_complete",
            all_mapped_atoms_have_coordinates,
            subject=selection.component_id,
            criterion="every identity-mapped atom has a coordinate record",
            observed=f"{len(missing_coordinate_atoms)} missing coordinate records",
            expected="0 missing",
            detail=(
                "missing: " + ", ".join(missing_coordinate_atoms)
                if missing_coordinate_atoms
                else "all mapped atoms have deposited coordinates"
            ),
        ),
        _check(
            "coordinate_mapping",
            "coordinates.alternate_locations_documented",
            alternate_locations_documented,
            subject=selection.component_id,
            criterion="alternate locations are preserved and not silently chosen",
            observed=(
                f"{len(alternate_location_atoms)} atoms with multiple records"
                if ambiguous
                else "no alternate-location duplicates"
            ),
            expected="none, or explicit records",
            detail=(
                "geometry is not computed while alternate locations remain unresolved"
                if ambiguous
                else "each atom name appears once"
            ),
        ),
        _check(
            "coordinate_quality",
            "coordinates.occupancies_recorded",
            not missing_occupancy,
            subject=selection.component_id,
            criterion="occupancy is recorded for every deposited atom",
            observed=f"{len(atoms) - len(missing_occupancy)} of {len(atoms)}",
            expected="all atoms",
            detail=(
                "missing: " + ", ".join(missing_occupancy)
                if missing_occupancy
                else "occupancies recorded without averaging"
            ),
        ),
        _check(
            "coordinate_quality",
            "coordinates.b_factors_recorded",
            not missing_b_factors,
            subject=selection.component_id,
            criterion="B-factor or displacement parameter is recorded for every atom",
            observed=f"{len(atoms) - len(missing_b_factors)} of {len(atoms)}",
            expected="all atoms",
            detail=(
                "missing: " + ", ".join(missing_b_factors)
                if missing_b_factors
                else "per-atom displacement parameters recorded"
            ),
        ),
        _check(
            "geometry",
            "geometry.bond_lengths_in_range",
            not bond_length_out_of_range and bool(bond_lengths),
            subject=selection.component_id,
            criterion="bond lengths fall within broad sanity bounds",
            observed=(
                f"{len(bond_length_out_of_range)} outside "
                f"{criteria.bond_length_min_angstrom}-{criteria.bond_length_max_angstrom} A"
            ),
            expected="0 outside",
            detail=(
                "broad physical sanity bounds, not force-field target lengths"
            ),
            not_run=ambiguous or not bond_lengths,
        ),
        _check(
            "geometry",
            "geometry.angles_in_range",
            not angle_out_of_range and bool(angle_values),
            subject=selection.component_id,
            criterion="bond angles fall within broad sanity bounds",
            observed=(
                f"{len(angle_out_of_range)} outside "
                f"{criteria.angle_min_degrees}-{criteria.angle_max_degrees} degrees"
            ),
            expected="0 outside",
            detail="broad physical sanity bounds, not force-field target angles",
            not_run=ambiguous or not angle_values,
        ),
        _check(
            "geometry",
            "geometry.ring_planarity",
            not non_planar_rings and bool(ring_planarity),
            subject=selection.component_id,
            criterion="aromatic ring systems are approximately planar",
            observed=f"{len(non_planar_rings)} above RMSD threshold",
            expected=f"RMSD <= {criteria.ring_planarity_rmsd_max_angstrom} A",
            detail=(
                f"{len(ring_planarity)} aromatic ring system(s) assessed"
                if ring_planarity
                else "no aromatic ring systems found"
            ),
            not_run=ambiguous or not ring_planarity,
        ),
        _check(
            "geometry",
            "geometry.close_contacts",
            not close_contacts,
            subject=selection.component_id,
            criterion="no nonbonded 1-4-or-greater contact is unusually short",
            observed=f"{len(close_contacts)} below {criteria.close_contact_min_angstrom} A",
            expected="0 unusually short contacts",
            detail="bonded and 1-3 pairs are excluded",
            not_run=ambiguous,
        ),
        _check(
            "completeness",
            "coordinates.missing_heavy_atoms",
            not missing_heavy_atoms,
            subject=selection.component_id,
            criterion="all heavy component atoms are observed",
            observed=f"{len(missing_heavy_atoms)} missing",
            expected="0 missing heavy atoms",
            detail=(
                "missing: " + ", ".join(missing_heavy_atoms)
                if missing_heavy_atoms
                else "all heavy component atoms have coordinates"
            ),
        ),
        _check(
            "completeness",
            "coordinates.missing_hydrogens",
            not missing_component_atoms,
            subject=selection.component_id,
            criterion="hydrogen coordinates are present",
            observed=f"{len(missing_component_atoms)} component atoms absent",
            expected="hydrogens are not required in this batch",
            detail="missing atoms are reported, not rebuilt",
            not_applicable=bool(missing_component_atoms)
            and not missing_heavy_atoms,
        ),
        _check(
            "stereochemistry",
            "coordinates.stereochemistry",
            False,
            subject=selection.component_id,
            criterion="coordinate-derived chirality matches the reference",
            observed=(
                "no reference stereocentres"
                if not reference_stereocenters
                else "reference stereocentres present; coordinate chirality not computed"
            ),
            expected="not applicable when reference has no stereocentres",
            detail=(
                "no stereochemical comparison is meaningful for this reference"
                if not reference_stereocenters
                else "stereochemical geometry is deferred"
            ),
            not_run=bool(reference_stereocenters),
            not_applicable=not reference_stereocenters,
        ),
    )

    discrepancies: list[CoordinateDiscrepancy] = []
    if missing_coordinate_atoms:
        discrepancies.append(
            CoordinateDiscrepancy(
                "missing_coordinates",
                "open",
                "Mapped atoms without coordinate records: "
                + ", ".join(missing_coordinate_atoms),
                observed=f"{len(missing_coordinate_atoms)} missing",
                expected="0 missing",
            )
        )
    if ambiguous:
        discrepancies.append(
            CoordinateDiscrepancy(
                "alternate_locations",
                "open",
                "Multiple coordinate records exist for at least one atom name.",
                observed=", ".join(alternate_location_atoms),
                expected="one resolved conformer or an explicit alternate-location policy",
                resolution="No conformer was selected and no coordinates were averaged.",
            )
        )
    if missing_occupancy:
        discrepancies.append(
            CoordinateDiscrepancy(
                "occupancy",
                "open",
                "Occupancy is missing for some deposited atoms.",
                observed=", ".join(missing_occupancy),
                expected="occupancy recorded for every atom",
            )
        )
    if missing_b_factors:
        discrepancies.append(
            CoordinateDiscrepancy(
                "b_factor",
                "open",
                "B-factor or displacement parameter is missing for some atoms.",
                observed=", ".join(missing_b_factors),
                expected="B-factor recorded for every atom",
            )
        )
    if bond_length_out_of_range:
        discrepancies.append(
            CoordinateDiscrepancy(
                "bond_length",
                "open",
                f"{len(bond_length_out_of_range)} bond lengths fall outside the broad sanity bounds.",
                observed=(
                    f"{min(bond_length_out_of_range):.3f}-"
                    f"{max(bond_length_out_of_range):.3f} A"
                ),
                expected=(
                    f"{criteria.bond_length_min_angstrom}-"
                    f"{criteria.bond_length_max_angstrom} A"
                ),
            )
        )
    if angle_out_of_range:
        discrepancies.append(
            CoordinateDiscrepancy(
                "bond_angle",
                "open",
                f"{len(angle_out_of_range)} bond angles fall outside the broad sanity bounds.",
                observed=(
                    f"{min(angle_out_of_range):.1f}-"
                    f"{max(angle_out_of_range):.1f} degrees"
                ),
                expected=(
                    f"{criteria.angle_min_degrees}-"
                    f"{criteria.angle_max_degrees} degrees"
                ),
            )
        )
    if non_planar_rings:
        discrepancies.append(
            CoordinateDiscrepancy(
                "ring_planarity",
                "open",
                f"{len(non_planar_rings)} aromatic ring system(s) exceed the planarity threshold.",
                observed=f"max RMSD {max(r.rmsd_angstrom for r in non_planar_rings):.3f} A",
                expected=f"RMSD <= {criteria.ring_planarity_rmsd_max_angstrom} A",
            )
        )
    if close_contacts:
        discrepancies.append(
            CoordinateDiscrepancy(
                "close_contact",
                "open",
                f"{len(close_contacts)} unusually short nonbonded contact(s) found.",
                observed=f"minimum {close_contacts[0].distance_angstrom:.3f} A",
                expected=f">= {criteria.close_contact_min_angstrom} A",
            )
        )
    if missing_heavy_atoms:
        discrepancies.append(
            CoordinateDiscrepancy(
                "missing_heavy_atoms",
                "open",
                "Heavy atoms present in the component dictionary are absent from the coordinates.",
                observed=", ".join(missing_heavy_atoms),
                expected="all heavy atoms observed",
            )
        )
    if missing_component_atoms and not missing_heavy_atoms:
        discrepancies.append(
            CoordinateDiscrepancy(
                "missing_hydrogens",
                "not_applicable",
                "Hydrogen atoms in the component dictionary are absent from the coordinates.",
                observed=f"{len(missing_component_atoms)} component atoms absent",
                expected="heavy-atom coordinate record only",
                resolution="No hydrogens were added.",
            )
        )
    if reference_stereocenters:
        discrepancies.append(
            CoordinateDiscrepancy(
                "stereochemistry",
                "open",
                "Coordinate-derived chirality was not compared with the reference stereocentres.",
                observed="not computed",
                expected=", ".join(reference_stereocenters),
                resolution="Deferred; no stereochemical assignment was applied.",
            )
        )

    decisions = (
        DecisionRecord(
            decision="preserved deposited coordinates in angstrom and recorded nm conversions",
            actor="ligand coordinate reconciliation",
            rationale=(
                "the source frame is retained; lengths and distances are converted "
                "to nm for analysis"
            ),
        ),
        DecisionRecord(
            decision="did not choose or average alternate conformers",
            actor="ligand coordinate reconciliation",
            rationale=(
                "alternate locations are recorded; geometry is withheld when they "
                "remain unresolved"
            ),
        ),
        DecisionRecord(
            decision="did not rebuild missing atoms or hydrogens",
            actor="ligand coordinate reconciliation",
            rationale=(
                "missing component atoms are reported from the deposition rather "
                "than invented"
            ),
        ),
        DecisionRecord(
            decision="used broad geometry sanity bounds, not force-field targets",
            actor="ligand coordinate reconciliation",
            rationale=(
                "the checks detect gross inconsistencies without pretending to be "
                "a parameterised geometry validation"
            ),
        ),
    )

    provenance = ProvenanceRecord(
        records=(
            ProvenanceEntry(
                "experimental",
                "deposited coordinates",
                (
                    f"Selected {selection.component_id} coordinates from the "
                    f"{accession} deposition."
                ),
                actor="RCSB PDB deposition",
            ),
            ProvenanceEntry(
                "calculated",
                "coordinate mapping and geometry",
                (
                    f"{len(atoms)} atom records, {len(bond_geometry)} bonds, "
                    f"{len(angle_geometry)} angles, {len(ring_planarity)} ring systems."
                ),
                actor="ligand coordinate reconciliation",
            ),
            ProvenanceEntry(
                "predicted",
                "none",
                "No docking pose, affinity, or free-energy estimate is produced.",
                actor="ligand coordinate reconciliation",
            ),
            ProvenanceEntry(
                "hypothesised",
                "none",
                "No mechanistic interpretation is made from the coordinates.",
                actor="ligand coordinate reconciliation",
            ),
        )
    )

    return LigandCoordinateRecord(
        source_manifest_path=source_manifest_path,
        source_manifest_sha256=source_manifest_sha256,
        identity_record_path=identity_record_path,
        identity_record_sha256=identity_record_sha256,
        accession=accession,
        selection=selection,
        criteria=criteria,
        atoms=tuple(atoms),
        bonds=tuple(bond_geometry),
        angles=tuple(angle_geometry),
        rings=tuple(ring_planarity),
        close_contacts=tuple(close_contacts),
        missing_component_atoms=missing_component_atoms,
        missing_heavy_atoms=missing_heavy_atoms,
        alternate_location_atoms=alternate_location_atoms,
        decisions=decisions,
        provenance=provenance,
        validation=ValidationRecord(checks=checks),
        discrepancies=tuple(discrepancies),
    )


def write_ligand_coordinates(path: str | Path, record: LigandCoordinateRecord) -> Path:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(record.to_dict(), indent=2, sort_keys=False) + "\n",
        encoding="utf-8",
    )
    return destination


def load_ligand_coordinates(path: str | Path) -> LigandCoordinateRecord:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    return LigandCoordinateRecord.from_dict(data)


def validate_ligand_coordinates(record: LigandCoordinateRecord) -> list[str]:
    """Return schema and internal-consistency errors for *record*."""

    errors: list[str] = []
    if record.schema_version != LIGAND_COORDINATE_SCHEMA_VERSION:
        errors.append(
            f"schema_version must be {LIGAND_COORDINATE_SCHEMA_VERSION!r}, got "
            f"{record.schema_version!r}."
        )
    if not record.source_manifest_path:
        errors.append("source_manifest_path is required.")
    if not _SHA256_PATTERN.fullmatch(record.source_manifest_sha256):
        errors.append("source_manifest_sha256 must be a lowercase SHA-256 digest.")
    if not record.identity_record_path:
        errors.append("identity_record_path is required.")
    if not _SHA256_PATTERN.fullmatch(record.identity_record_sha256):
        errors.append("identity_record_sha256 must be a lowercase SHA-256 digest.")
    if not record.accession:
        errors.append("accession is required.")
    if not record.atoms:
        errors.append("atoms must not be empty.")

    keys = [(atom.atom_name, atom.alternate_location) for atom in record.atoms]
    if len(set(keys)) != len(keys):
        errors.append("atom name and alternate-location pairs must be unique.")
    atom_names = {atom.atom_name for atom in record.atoms}
    for index, atom in enumerate(record.atoms):
        if not (
            math.isfinite(atom.position_angstrom[0])
            and math.isfinite(atom.position_angstrom[1])
            and math.isfinite(atom.position_angstrom[2])
        ):
            errors.append(f"atoms[{index}].position_angstrom must be finite.")
        if atom.occupancy is not None and not 0.0 <= atom.occupancy <= 1.0:
            errors.append(f"atoms[{index}].occupancy must be between 0 and 1.")

    for index, bond in enumerate(record.bonds):
        if bond.atom_name_1 not in atom_names or bond.atom_name_2 not in atom_names:
            errors.append(f"bonds[{index}] references an unknown atom name.")
        if bond.length_angstrom <= 0 or bond.length_nm <= 0:
            errors.append(f"bonds[{index}] length must be positive.")

    for index, angle in enumerate(record.angles):
        if not (
            angle.atom_name_1 in atom_names
            and angle.atom_name_2 in atom_names
            and angle.atom_name_3 in atom_names
        ):
            errors.append(f"angles[{index}] references an unknown atom name.")
        if not 0.0 <= angle.angle_degrees <= 180.0:
            errors.append(f"angles[{index}] must be between 0 and 180 degrees.")

    for index, ring in enumerate(record.rings):
        if not set(ring.atom_names).issubset(atom_names):
            errors.append(f"rings[{index}] references an unknown atom name.")
        if ring.rmsd_angstrom < 0 or ring.max_deviation_angstrom < 0:
            errors.append(f"rings[{index}] deviations cannot be negative.")

    for index, contact in enumerate(record.close_contacts):
        if (
            contact.atom_name_1 not in atom_names
            or contact.atom_name_2 not in atom_names
        ):
            errors.append(f"close_contacts[{index}] references an unknown atom name.")
        if contact.distance_angstrom <= 0 or contact.distance_nm <= 0:
            errors.append(f"close_contacts[{index}] distance must be positive.")

    criteria = record.criteria
    if criteria.bond_length_min_angstrom <= 0:
        errors.append("criteria.bond_length_min_angstrom must be positive.")
    if criteria.bond_length_max_angstrom <= criteria.bond_length_min_angstrom:
        errors.append(
            "criteria.bond_length_max_angstrom must exceed the minimum."
        )
    if not 0 < criteria.angle_min_degrees < criteria.angle_max_degrees <= 180:
        errors.append("criteria angle bounds must satisfy 0 < min < max <= 180.")
    if criteria.ring_planarity_rmsd_max_angstrom <= 0:
        errors.append("criteria.ring_planarity_rmsd_max_angstrom must be positive.")
    if criteria.close_contact_min_angstrom <= 0:
        errors.append("criteria.close_contact_min_angstrom must be positive.")

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
