"""Three-state evidence fingerprints and comparisons.

A conventional binary fingerprint records presence or absence.  That is not
enough when chemical evidence is incomplete: an interaction can be present,
absent, or *unknown* because the features needed to assess it were never
supplied.  Unknown must not be encoded as zero.

This module is deliberately independent of structure parsing and chemistry
toolkits.  Callers supply explicit features and observations.  A feature can be
category-level (for a whole interaction class) or atom-pair-level (for a
specific ligand/protein/water feature).
"""

from __future__ import annotations

import json
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any

EVIDENCE_STATES = ("present", "absent", "unknown")
EVIDENCE_SCOPES = ("category", "atom_pair")


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


def _optional_string(data: Mapping[str, Any], key: str) -> str:
    return str(data.get(key, ""))


@dataclass(frozen=True, order=True)
class EvidenceFeature:
    """Stable identity for one evidence feature.

    A category feature names only the interaction class.  An atom-pair feature
    names the protein residue/atom, ligand atom/residue, and, for
    water-mediated interactions, the water identifier.
    """

    scope: str
    interaction_type: str
    protein_residue: str = ""
    protein_atom: str = ""
    ligand_atom: str = ""
    ligand_residue: str = ""
    water_id: str = ""

    def __post_init__(self) -> None:
        if self.scope not in EVIDENCE_SCOPES:
            raise ValueError(f"scope must be one of {EVIDENCE_SCOPES}.")
        if not self.interaction_type.strip():
            raise ValueError("interaction_type must be a non-empty string.")
        if self.scope == "atom_pair":
            if not self.protein_residue.strip():
                raise ValueError("atom_pair features require protein_residue.")
            if not self.protein_atom.strip():
                raise ValueError("atom_pair features require protein_atom.")
            if not self.ligand_atom.strip():
                raise ValueError("atom_pair features require ligand_atom.")
        elif any(
            value
            for value in (
                self.protein_residue,
                self.protein_atom,
                self.ligand_atom,
                self.ligand_residue,
                self.water_id,
            )
        ):
            raise ValueError("category features must not carry atom identifiers.")

    @property
    def key(self) -> tuple[str, ...]:
        return (
            self.scope,
            self.interaction_type,
            self.protein_residue,
            self.protein_atom,
            self.ligand_atom,
            self.ligand_residue,
            self.water_id,
        )

    @property
    def is_water_mediated(self) -> bool:
        return bool(self.water_id)

    def to_dict(self) -> dict[str, Any]:
        return {
            "scope": self.scope,
            "interaction_type": self.interaction_type,
            "protein_residue": self.protein_residue,
            "protein_atom": self.protein_atom,
            "ligand_atom": self.ligand_atom,
            "ligand_residue": self.ligand_residue,
            "water_id": self.water_id,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> EvidenceFeature:
        data = _as_mapping(data, "evidence feature")
        return cls(
            scope=_required_string(data, "scope"),
            interaction_type=_required_string(data, "interaction_type"),
            protein_residue=_optional_string(data, "protein_residue"),
            protein_atom=_optional_string(data, "protein_atom"),
            ligand_atom=_optional_string(data, "ligand_atom"),
            ligand_residue=_optional_string(data, "ligand_residue"),
            water_id=_optional_string(data, "water_id"),
        )


@dataclass(frozen=True)
class EvidenceObservation:
    feature: EvidenceFeature
    state: str
    assessable: bool
    evidence: str
    criteria: str
    provenance: str
    distance_nm: float | None = None
    angle_degrees: float | None = None

    def __post_init__(self) -> None:
        if self.state not in EVIDENCE_STATES:
            raise ValueError(f"state must be one of {EVIDENCE_STATES}.")
        if self.state in ("present", "absent") and not self.assessable:
            raise ValueError("present/absent observations must be assessable.")
        if self.state == "unknown" and self.assessable:
            raise ValueError("unknown observations must not be assessable.")
        if not self.evidence.strip():
            raise ValueError("evidence must be a non-empty string.")
        if not self.criteria.strip():
            raise ValueError("criteria must be a non-empty string.")
        if not self.provenance.strip():
            raise ValueError("provenance must be a non-empty string.")
        if self.distance_nm is not None and self.distance_nm <= 0:
            raise ValueError("distance_nm must be positive when provided.")
        if self.angle_degrees is not None and not 0.0 <= self.angle_degrees <= 180.0:
            raise ValueError("angle_degrees must be in [0, 180].")

    def to_dict(self) -> dict[str, Any]:
        return {
            "feature": self.feature.to_dict(),
            "state": self.state,
            "assessable": self.assessable,
            "evidence": self.evidence,
            "criteria": self.criteria,
            "provenance": self.provenance,
            "distance_nm": self.distance_nm,
            "angle_degrees": self.angle_degrees,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> EvidenceObservation:
        data = _as_mapping(data, "evidence observation")
        return cls(
            feature=EvidenceFeature.from_dict(_as_mapping(data.get("feature"), "feature")),
            state=_required_string(data, "state"),
            assessable=bool(data.get("assessable")),
            evidence=_required_string(data, "evidence"),
            criteria=_required_string(data, "criteria"),
            provenance=_required_string(data, "provenance"),
            distance_nm=(
                None if data.get("distance_nm") is None else float(data["distance_nm"])
            ),
            angle_degrees=(
                None
                if data.get("angle_degrees") is None
                else float(data["angle_degrees"])
            ),
        )


@dataclass(frozen=True)
class EvidenceSchema:
    schema_id: str
    criteria_fingerprint: str
    features: tuple[EvidenceFeature, ...]

    def __post_init__(self) -> None:
        if not self.schema_id.strip():
            raise ValueError("schema_id must be a non-empty string.")
        if not self.criteria_fingerprint.strip():
            raise ValueError("criteria_fingerprint must be a non-empty string.")
        keys = [feature.key for feature in self.features]
        if len(set(keys)) != len(keys):
            raise ValueError("evidence features must be unique.")
        if list(keys) != sorted(keys):
            raise ValueError("evidence features must be sorted by key.")

    @property
    def compatibility_key(self) -> tuple[str, str, tuple[tuple[str, ...], ...]]:
        return (
            self.schema_id,
            self.criteria_fingerprint,
            tuple(feature.key for feature in self.features),
        )

    @classmethod
    def from_observations(
        cls,
        observations: Sequence[EvidenceObservation],
        *,
        schema_id: str,
        criteria_fingerprint: str,
    ) -> EvidenceSchema:
        features = sorted(
            {observation.feature for observation in observations},
            key=lambda feature: feature.key,
        )
        return cls(schema_id, criteria_fingerprint, tuple(features))

    def index(self, feature: EvidenceFeature) -> int:
        for index, candidate in enumerate(self.features):
            if candidate == feature:
                return index
        raise ValueError(f"feature {feature!r} is not in the schema.")

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_id": self.schema_id,
            "criteria_fingerprint": self.criteria_fingerprint,
            "features": [feature.to_dict() for feature in self.features],
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> EvidenceSchema:
        data = _as_mapping(data, "evidence schema")
        return cls(
            schema_id=_required_string(data, "schema_id"),
            criteria_fingerprint=_required_string(data, "criteria_fingerprint"),
            features=tuple(
                EvidenceFeature.from_dict(_as_mapping(item, "evidence feature"))
                for item in _as_list(data.get("features"), "features")
            ),
        )


@dataclass(frozen=True)
class EvidenceFingerprint:
    schema: EvidenceSchema
    states: tuple[str, ...]
    assessable: tuple[bool, ...]
    observations: tuple[EvidenceObservation, ...] = ()

    def __post_init__(self) -> None:
        if len(self.states) != len(self.schema.features):
            raise ValueError("states must have one entry per schema feature.")
        if len(self.assessable) != len(self.schema.features):
            raise ValueError("assessable must have one entry per schema feature.")
        for state, assessable in zip(self.states, self.assessable):
            if state not in EVIDENCE_STATES:
                raise ValueError(f"invalid evidence state {state!r}.")
            if (state == "unknown") == assessable:
                raise ValueError(
                    "unknown states must be unassessable and present/absent assessable."
                )

    def index(self, feature: EvidenceFeature) -> int:
        return self.schema.index(feature)

    def state_for(self, feature: EvidenceFeature) -> str:
        return self.states[self.index(feature)]

    def assessable_for(self, feature: EvidenceFeature) -> bool:
        return self.assessable[self.index(feature)]

    @property
    def present_count(self) -> int:
        return sum(state == "present" for state in self.states)

    @property
    def absent_count(self) -> int:
        return sum(state == "absent" for state in self.states)

    @property
    def unknown_count(self) -> int:
        return sum(state == "unknown" for state in self.states)

    @property
    def assessable_count(self) -> int:
        return sum(self.assessable)

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": self.schema.to_dict(),
            "states": list(self.states),
            "assessable": list(self.assessable),
            "observations": [observation.to_dict() for observation in self.observations],
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> EvidenceFingerprint:
        data = _as_mapping(data, "evidence fingerprint")
        return cls(
            schema=EvidenceSchema.from_dict(_as_mapping(data.get("schema"), "schema")),
            states=tuple(str(state) for state in _as_list(data.get("states"), "states")),
            assessable=tuple(bool(value) for value in _as_list(data.get("assessable"), "assessable")),
            observations=tuple(
                EvidenceObservation.from_dict(_as_mapping(item, "evidence observation"))
                for item in _as_list(data.get("observations"), "observations")
            ),
        )


@dataclass(frozen=True)
class EvidenceOccupancy:
    schema: EvidenceSchema
    present_counts: tuple[int, ...]
    absent_counts: tuple[int, ...]
    unknown_counts: tuple[int, ...]
    assessable_counts: tuple[int, ...]

    def __post_init__(self) -> None:
        size = len(self.schema.features)
        for name, values in (
            ("present_counts", self.present_counts),
            ("absent_counts", self.absent_counts),
            ("unknown_counts", self.unknown_counts),
            ("assessable_counts", self.assessable_counts),
        ):
            if len(values) != size:
                raise ValueError(f"{name} must have one entry per feature.")
            if any(value < 0 for value in values):
                raise ValueError(f"{name} cannot contain negative values.")
        for present, absent, unknown, assessable in zip(
            self.present_counts,
            self.absent_counts,
            self.unknown_counts,
            self.assessable_counts,
            strict=True,
        ):
            if assessable != present + absent:
                raise ValueError("assessable counts must equal present plus absent.")

    @property
    def present_fractions(self) -> tuple[float | None, ...]:
        return tuple(
            None if assessable == 0 else present / assessable
            for present, assessable in zip(
                self.present_counts, self.assessable_counts, strict=True
            )
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": self.schema.to_dict(),
            "present_counts": list(self.present_counts),
            "absent_counts": list(self.absent_counts),
            "unknown_counts": list(self.unknown_counts),
            "assessable_counts": list(self.assessable_counts),
            "present_fractions": list(self.present_fractions),
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> EvidenceOccupancy:
        data = _as_mapping(data, "evidence occupancy")
        return cls(
            schema=EvidenceSchema.from_dict(_as_mapping(data.get("schema"), "schema")),
            present_counts=tuple(int(value) for value in _as_list(data.get("present_counts"), "present_counts")),
            absent_counts=tuple(int(value) for value in _as_list(data.get("absent_counts"), "absent_counts")),
            unknown_counts=tuple(int(value) for value in _as_list(data.get("unknown_counts"), "unknown_counts")),
            assessable_counts=tuple(int(value) for value in _as_list(data.get("assessable_counts"), "assessable_counts")),
        )


@dataclass(frozen=True)
class EvidenceComparison:
    schema_id: str
    criteria_fingerprint: str
    both_present: int
    both_absent: int
    first_only_present: int
    second_only_present: int
    unknown_involved: int
    assessable_both: int


def build_evidence_fingerprint(
    observations: Sequence[EvidenceObservation],
    *,
    schema_id: str,
    criteria_fingerprint: str,
) -> EvidenceFingerprint:
    """Combine observations into one three-state fingerprint.

    Duplicate feature observations are merged conservatively: any present
    observation makes the feature present; otherwise any unknown observation
    makes it unknown; otherwise it is absent.
    """

    values = tuple(observations)
    if not values:
        raise ValueError("observations must not be empty.")
    schema = EvidenceSchema.from_observations(
        values,
        schema_id=schema_id,
        criteria_fingerprint=criteria_fingerprint,
    )
    by_feature: dict[EvidenceFeature, list[EvidenceObservation]] = {}
    for observation in values:
        by_feature.setdefault(observation.feature, []).append(observation)

    states = []
    assessable = []
    for feature in schema.features:
        feature_observations = by_feature[feature]
        if any(observation.state == "present" for observation in feature_observations):
            state = "present"
        elif any(observation.state == "unknown" for observation in feature_observations):
            state = "unknown"
        else:
            state = "absent"
        states.append(state)
        assessable.append(state != "unknown")

    return EvidenceFingerprint(
        schema=schema,
        states=tuple(states),
        assessable=tuple(assessable),
        observations=values,
    )


def aggregate_evidence(
    fingerprints: Sequence[EvidenceFingerprint],
) -> EvidenceOccupancy:
    """Aggregate compatible fingerprints using assessable denominators."""

    values = tuple(fingerprints)
    if not values:
        raise ValueError("fingerprints must not be empty.")
    schema = values[0].schema
    for fingerprint in values[1:]:
        if fingerprint.schema.compatibility_key != schema.compatibility_key:
            raise ValueError("fingerprints have incompatible evidence schemas.")

    present = [0] * len(schema.features)
    absent = [0] * len(schema.features)
    unknown = [0] * len(schema.features)
    assessable = [0] * len(schema.features)
    for fingerprint in values:
        for index, (state, is_assessable) in enumerate(
            zip(fingerprint.states, fingerprint.assessable, strict=True)
        ):
            if state == "present":
                present[index] += 1
            elif state == "absent":
                absent[index] += 1
            else:
                unknown[index] += 1
            if is_assessable:
                assessable[index] += 1
    return EvidenceOccupancy(
        schema=schema,
        present_counts=tuple(present),
        absent_counts=tuple(absent),
        unknown_counts=tuple(unknown),
        assessable_counts=tuple(assessable),
    )


def compare_evidence_fingerprints(
    first: EvidenceFingerprint,
    second: EvidenceFingerprint,
) -> EvidenceComparison:
    """Compare compatible fingerprints without treating unknown as absent."""

    if first.schema.compatibility_key != second.schema.compatibility_key:
        raise ValueError("fingerprints have incompatible evidence schemas.")
    both_present = both_absent = first_only = second_only = unknown = assessable_both = 0
    for first_state, second_state in zip(
        first.states, second.states, strict=True
    ):
        if first_state == "unknown" or second_state == "unknown":
            unknown += 1
            continue
        assessable_both += 1
        if first_state == "present" and second_state == "present":
            both_present += 1
        elif first_state == "absent" and second_state == "absent":
            both_absent += 1
        elif first_state == "present":
            first_only += 1
        else:
            second_only += 1
    return EvidenceComparison(
        schema_id=first.schema.schema_id,
        criteria_fingerprint=first.schema.criteria_fingerprint,
        both_present=both_present,
        both_absent=both_absent,
        first_only_present=first_only,
        second_only_present=second_only,
        unknown_involved=unknown,
        assessable_both=assessable_both,
    )


def criteria_fingerprint(criteria: Mapping[str, Any]) -> str:
    """Return a stable digest for a criteria mapping."""

    import hashlib

    encoded = json.dumps(criteria, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()
