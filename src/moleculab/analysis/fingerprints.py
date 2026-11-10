"""Interaction fingerprints and trajectory occupancy.

The fingerprint layer consumes interaction records from
:mod:`moleculab.analysis.interactions`. It does not read trajectory files and
does not run MD. Frames are supplied by the caller, so the same functions work
for MD snapshots, experimental ensembles, or docking poses.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

from .interactions import HydrogenBond, HydrogenBondCriteria


@dataclass(frozen=True)
class InteractionFeature:
    """A stable feature identity for a residue-level interaction.

    The identity uses the interaction type, protein residue identifier, ligand
    atom identifier, and ligand residue identifier. These are strings supplied
    by the structure adapter; array positions are not used as molecular
    identity.
    """

    interaction_type: str
    protein_residue: str
    ligand_atom_id: str
    ligand_residue: str = ""

    def __post_init__(self) -> None:
        for name in ("interaction_type", "protein_residue", "ligand_atom_id"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value.strip():
                raise ValueError(f"{name} must be a non-empty string.")
        if not isinstance(self.ligand_residue, str):
            raise TypeError("ligand_residue must be a string.")

    @classmethod
    def from_hydrogen_bond(cls, bond: HydrogenBond) -> InteractionFeature:
        """Build a feature identity from one hydrogen-bond record."""

        return cls(
            interaction_type=bond.interaction_type,
            protein_residue=bond.protein_residue or "",
            ligand_atom_id=bond.ligand_label,
            ligand_residue=bond.ligand_residue or "",
        )

    @property
    def key(self) -> tuple[str, str, str, str]:
        """Sortable feature identity."""

        return (
            self.interaction_type,
            self.protein_residue,
            self.ligand_residue,
            self.ligand_atom_id,
        )


@dataclass(frozen=True)
class FingerprintSchema:
    """Ordered feature definition shared by a set of frames."""

    features: tuple[InteractionFeature, ...]

    def __post_init__(self) -> None:
        if len(set(self.features)) != len(self.features):
            raise ValueError("features must be unique.")

    @classmethod
    def from_interactions(
        cls,
        interactions: Sequence[HydrogenBond],
    ) -> FingerprintSchema:
        """Build a deterministic schema from interaction records."""

        features = sorted(
            {InteractionFeature.from_hydrogen_bond(bond) for bond in interactions},
            key=lambda feature: feature.key,
        )
        return cls(tuple(features))

    def __len__(self) -> int:
        return len(self.features)

    def index(self, feature: InteractionFeature) -> int:
        """Return the position of *feature* in the schema."""

        for index, candidate in enumerate(self.features):
            if candidate == feature:
                return index
        raise ValueError(f"feature {feature!r} is not in the schema.")


@dataclass(frozen=True)
class InteractionFingerprint:
    """Binary and count-based fingerprint for one frame.

    ``bits[j]`` is 1 when feature ``j`` is present at least once. ``counts[j]``
    records how many observations mapped to that feature, so duplicate
    observations do not inflate the binary fingerprint.
    """

    schema: FingerprintSchema
    bits: np.ndarray
    counts: np.ndarray
    criteria: HydrogenBondCriteria | None = None

    def __post_init__(self) -> None:
        bits = np.asarray(self.bits, dtype=bool).reshape(-1)
        counts = np.asarray(self.counts, dtype=int).reshape(-1)
        if bits.shape != (len(self.schema),):
            raise ValueError("bits must have one entry per schema feature.")
        if counts.shape != (len(self.schema),):
            raise ValueError("counts must have one entry per schema feature.")
        if np.any(counts < 0):
            raise ValueError("counts must be non-negative.")
        if not np.array_equal(bits, counts > 0):
            raise ValueError("bits must equal counts > 0.")
        object.__setattr__(self, "bits", bits.copy())
        object.__setattr__(self, "counts", counts.copy())

    @property
    def features(self) -> tuple[InteractionFeature, ...]:
        """Features in schema order."""

        return self.schema.features

    @property
    def present(self) -> tuple[InteractionFeature, ...]:
        """Features present in this frame."""

        return tuple(
            feature
            for feature, bit in zip(self.features, self.bits, strict=True)
            if bit
        )


@dataclass(frozen=True)
class OccupancyResult:
    """Occupancy summary for a set of frame fingerprints."""

    schema: FingerprintSchema
    occupancy: np.ndarray
    counts: np.ndarray
    n_frames: int
    n_valid_frames: int
    n_missing_frames: int
    criteria: HydrogenBondCriteria | None = None

    @property
    def features(self) -> tuple[InteractionFeature, ...]:
        """Features in schema order."""

        return self.schema.features


def _as_hydrogen_bonds(
    interactions: Sequence[HydrogenBond],
) -> tuple[HydrogenBond, ...]:
    bonds = tuple(interactions)
    for bond in bonds:
        if not isinstance(bond, HydrogenBond):
            raise TypeError("interactions must contain HydrogenBond records.")
    return bonds


def _common_criteria(
    bonds: Sequence[HydrogenBond],
) -> HydrogenBondCriteria | None:
    if not bonds:
        return None
    criteria = bonds[0].criteria
    for bond in bonds[1:]:
        if bond.criteria != criteria:
            raise ValueError("all interactions in a frame must use the same criteria.")
    return criteria


def build_fingerprint(
    interactions: Sequence[HydrogenBond],
    schema: FingerprintSchema | None = None,
) -> InteractionFingerprint:
    """Build one binary/count fingerprint from interaction records.

    If *schema* is omitted, the schema is derived deterministically from the
    interactions. Supply a schema explicitly when building fingerprints for
    several frames; otherwise a feature absent from one frame does not appear
    in that frame's schema.
    """

    bonds = _as_hydrogen_bonds(interactions)
    if schema is None:
        schema = FingerprintSchema.from_interactions(bonds)
    elif not isinstance(schema, FingerprintSchema):
        raise TypeError("schema must be a FingerprintSchema instance.")

    counts = np.zeros(len(schema), dtype=int)
    for bond in bonds:
        feature = InteractionFeature.from_hydrogen_bond(bond)
        counts[schema.index(feature)] += 1

    return InteractionFingerprint(
        schema=schema,
        bits=counts > 0,
        counts=counts,
        criteria=_common_criteria(bonds),
    )


def calculate_occupancy(
    fingerprints: Sequence[InteractionFingerprint | None],
) -> OccupancyResult:
    """Return feature occupancy across frame fingerprints.

    A ``None`` frame is treated as **missing or invalid**, not as a frame in
    which every interaction is absent. Occupancy is computed over valid frames
    only:

    .. math::

        O_j = \\frac{1}{N_{\\mathrm{valid}}}\\sum_{t \\in \\mathrm{valid}} f_j(t).

    All valid fingerprints must share the same schema. Detection criteria are
    retained when available; a frame with no interactions may not carry a
    criteria object, so it inherits the common criteria from the other valid
    frames.
    """

    frames = tuple(fingerprints)
    if not frames:
        raise ValueError("fingerprints must contain at least one frame.")
    valid = [frame for frame in frames if frame is not None]
    if not valid:
        raise ValueError("at least one frame must be valid.")
    for frame in valid:
        if not isinstance(frame, InteractionFingerprint):
            raise TypeError("fingerprints must contain InteractionFingerprint records.")

    schema = valid[0].schema
    criteria_values = [frame.criteria for frame in valid if frame.criteria is not None]
    criteria = criteria_values[0] if criteria_values else None
    for frame in valid[1:]:
        if frame.schema != schema:
            raise ValueError("all valid frames must use the same fingerprint schema.")
    for candidate in criteria_values[1:]:
        if candidate != criteria:
            raise ValueError("all valid frames must use the same detection criteria.")

    matrix = np.stack([frame.bits.astype(float) for frame in valid], axis=0)
    occupancy = np.mean(matrix, axis=0)
    counts = np.sum(matrix, axis=0).astype(int)
    return OccupancyResult(
        schema=schema,
        occupancy=occupancy,
        counts=counts,
        n_frames=len(frames),
        n_valid_frames=len(valid),
        n_missing_frames=len(frames) - len(valid),
        criteria=criteria,
    )
