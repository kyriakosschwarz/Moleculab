"""Geometric contacts and their persistence across coordinate frames.

This module works with Cartesian positions and boolean contact matrices. It
does not load trajectories, infer chemical interactions, or assume that a
contact is favourable. Contact occupancy and contact lifetimes are reported as
separate quantities.
"""

from __future__ import annotations

from collections.abc import Hashable, Sequence
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class ResidueContactMatrix:
    """Residue-level contact matrix for one frame."""

    residue_ids: tuple[Hashable, ...]
    contacts: np.ndarray

    def __post_init__(self) -> None:
        contacts = np.asarray(self.contacts, dtype=bool)
        if contacts.ndim != 2:
            raise ValueError("contacts must be a two-dimensional boolean matrix.")
        if contacts.shape[0] != len(self.residue_ids):
            raise ValueError("one residue ID is required per contact-matrix row.")
        object.__setattr__(self, "contacts", contacts.copy())


@dataclass(frozen=True)
class ContactOccupancyResult:
    """Sampled-frame occupancy for a contact matrix."""

    occupancy: np.ndarray
    counts: np.ndarray
    n_frames: int
    n_valid_frames: int
    n_missing_frames: int

    def __post_init__(self) -> None:
        occupancy = np.asarray(self.occupancy, dtype=float)
        counts = np.asarray(self.counts, dtype=int)
        if occupancy.shape != counts.shape:
            raise ValueError("occupancy and counts must have the same shape.")
        if not np.all(np.isfinite(occupancy)):
            raise ValueError("occupancy must contain only finite values.")
        if np.any(occupancy < 0.0) or np.any(occupancy > 1.0):
            raise ValueError("occupancy values must lie in [0, 1].")
        if np.any(counts < 0):
            raise ValueError("counts must be non-negative.")
        object.__setattr__(self, "occupancy", occupancy.copy())
        object.__setattr__(self, "counts", counts.copy())


@dataclass(frozen=True)
class ContactLifetimeResult:
    """Run lengths of one contact and optional sampled durations.

    ``run_lengths_frames`` counts consecutive sampled frames in which the
    contact is present. If ``frame_spacing`` is supplied, ``lifetimes`` is
    ``run_length * frame_spacing``; this is a sampled-duration convention, not
    a physical residence time.
    """

    run_lengths_frames: tuple[int, ...]
    frame_spacing: float | None
    observed_fraction: float
    n_frames: int
    n_missing_frames: int

    @property
    def n_valid_frames(self) -> int:
        """Number of frames with an observed contact state."""

        return self.n_frames - self.n_missing_frames
    @property
    def n_runs(self) -> int:
        """Number of separate contact appearances."""

        return len(self.run_lengths_frames)

    @property
    def mean_lifetime_frames(self) -> float:
        """Mean run length in sampled frames."""

        if not self.run_lengths_frames:
            return 0.0
        return float(np.mean(self.run_lengths_frames))

    @property
    def max_lifetime_frames(self) -> int:
        """Longest run length in sampled frames."""

        return max(self.run_lengths_frames, default=0)

    @property
    def lifetimes(self) -> tuple[float, ...] | None:
        """Run durations in the frame-spacing unit, when spacing is known."""

        if self.frame_spacing is None:
            return None
        return tuple(length * self.frame_spacing for length in self.run_lengths_frames)

    @property
    def mean_lifetime(self) -> float | None:
        """Mean sampled duration when frame spacing is known."""

        lifetimes = self.lifetimes
        if lifetimes is None or not lifetimes:
            return None if lifetimes is None else 0.0
        return float(np.mean(lifetimes))


def _as_positions(name: str, positions: np.ndarray) -> np.ndarray:
    coordinates = np.asarray(positions, dtype=float)
    if coordinates.ndim != 2 or coordinates.shape[1] != 3:
        raise ValueError(f"{name} must have shape (n_atoms, 3).")
    if not np.all(np.isfinite(coordinates)):
        raise ValueError(f"{name} must contain only finite values.")
    return coordinates


def _as_cutoff(cutoff: float) -> float:
    value = float(cutoff)
    if not np.isfinite(value) or value <= 0.0:
        raise ValueError("cutoff must be finite and positive.")
    return value


def calculate_contact_matrix(
    protein_positions: np.ndarray,
    ligand_positions: np.ndarray,
    cutoff: float,
    *,
    inclusive: bool = False,
) -> np.ndarray:
    """Return a boolean protein-atom by ligand-atom contact matrix.

    The default comparison is strict: ``distance < cutoff``. Set
    ``inclusive=True`` for ``distance <= cutoff``. Positions are in nm and no
    periodic minimum image is applied.
    """

    protein = _as_positions("protein_positions", protein_positions)
    ligand = _as_positions("ligand_positions", ligand_positions)
    limit = _as_cutoff(cutoff)
    if ligand.shape[0] == 0:
        raise ValueError("ligand_positions must contain at least one atom.")

    contacts = np.empty((protein.shape[0], ligand.shape[0]), dtype=bool)
    for protein_index, position in enumerate(protein):
        distances = np.linalg.norm(ligand - position, axis=1)
        contacts[protein_index] = (
            distances <= limit if inclusive else distances < limit
        )
    return contacts


def aggregate_contacts_by_residue(
    contact_matrix: np.ndarray,
    protein_residue_ids: Sequence[Hashable],
) -> ResidueContactMatrix:
    """Aggregate atom-level contacts into residue-level contacts.

    A residue contacts a ligand atom when at least one of its protein atoms
    contacts that ligand atom. Residues appear in first-occurrence order.
    """

    contacts = np.asarray(contact_matrix, dtype=bool)
    if contacts.ndim != 2:
        raise ValueError("contact_matrix must be two-dimensional.")
    labels = tuple(protein_residue_ids)
    if len(labels) != contacts.shape[0]:
        raise ValueError("one residue ID is required per contact-matrix row.")

    grouped: dict[Hashable, list[int]] = {}
    for index, label in enumerate(labels):
        grouped.setdefault(label, []).append(index)

    residue_ids = tuple(grouped)
    residue_contacts = np.zeros((len(residue_ids), contacts.shape[1]), dtype=bool)
    for row, label in enumerate(residue_ids):
        residue_contacts[row] = np.any(contacts[grouped[label]], axis=0)
    return ResidueContactMatrix(residue_ids=residue_ids, contacts=residue_contacts)


def calculate_contact_occupancy(
    contact_frames: Sequence[np.ndarray | None],
) -> ContactOccupancyResult:
    """Return sampled-frame contact occupancy over valid frames.

    A ``None`` frame is missing or invalid, not a frame with no contacts.
    Occupancy is computed over valid frames only.
    """

    frames = tuple(contact_frames)
    if not frames:
        raise ValueError("contact_frames must contain at least one frame.")
    valid = [frame for frame in frames if frame is not None]
    if not valid:
        raise ValueError("at least one contact frame must be valid.")

    shape = None
    matrices = []
    for frame in valid:
        matrix = np.asarray(frame, dtype=bool)
        if matrix.ndim != 2:
            raise ValueError("each contact frame must be a two-dimensional matrix.")
        if shape is None:
            shape = matrix.shape
        elif matrix.shape != shape:
            raise ValueError("all contact frames must have the same shape.")
        matrices.append(matrix)

    stack = np.stack(matrices, axis=0).astype(float)
    occupancy = np.mean(stack, axis=0)
    counts = np.sum(stack, axis=0).astype(int)
    return ContactOccupancyResult(
        occupancy=occupancy,
        counts=counts,
        n_frames=len(frames),
        n_valid_frames=len(valid),
        n_missing_frames=len(frames) - len(valid),
    )


def calculate_contact_lifetimes(
    contact_presence: Sequence[bool | None],
    *,
    frame_spacing: float | None = None,
) -> ContactLifetimeResult:
    """Return run lengths for one contact's sampled presence.

    ``contact_presence`` is a one-dimensional sequence of booleans or missing
    values (``None`` or ``NaN``). A missing frame breaks observed continuity: a
    run is not inferred across it. This is conservative; bridging a missing
    frame requires an explicit policy and is not assumed here.

    With ``frame_spacing`` the sampled lifetime is
    ``run_length * frame_spacing``. This convention treats each sampled frame as
    one interval and is not a physical residence time.
    """

    raw = list(contact_presence)
    if not raw:
        raise ValueError("contact_presence must contain at least one frame.")

    parsed: list[bool | None] = []
    missing = 0
    for value in raw:
        if value is None or (
            isinstance(value, (float, np.floating)) and np.isnan(value)
        ):
            parsed.append(None)
            missing += 1
        else:
            parsed.append(bool(value))
    if missing == len(parsed):
        raise ValueError("at least one contact frame must be valid.")

    spacing = None
    if frame_spacing is not None:
        spacing = float(frame_spacing)
        if not np.isfinite(spacing) or spacing <= 0.0:
            raise ValueError("frame_spacing must be finite and positive.")

    runs: list[int] = []
    current = 0
    for present in parsed:
        if present is True:
            current += 1
        else:
            if current:
                runs.append(current)
            current = 0
    if current:
        runs.append(current)

    valid = len(parsed) - missing
    true_count = sum(value is True for value in parsed)
    return ContactLifetimeResult(
        run_lengths_frames=tuple(runs),
        frame_spacing=spacing,
        observed_fraction=true_count / valid,
        n_frames=len(parsed),
        n_missing_frames=missing,
    )
