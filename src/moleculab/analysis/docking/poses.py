# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Ligand pose record for external docking results."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field

import numpy as np


@dataclass(frozen=True)
class LigandPose:
    """Coordinates, atom identifiers, and score provenance for one ligand pose.

    Positions are in nm. ``atom_ids`` must be stable, unique identifiers; array
    order is not used as molecular identity by the mapping helpers.

    If a score is present, both ``score_unit`` and ``score_source`` are
    required. The score is metadata from the producing method, not a
    thermodynamic binding free energy.
    """

    positions: np.ndarray
    atom_ids: tuple[str, ...]
    score: float | None = None
    score_unit: str | None = None
    score_source: str = ""
    source: str = ""
    metadata: Mapping[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        positions = np.asarray(self.positions, dtype=float)
        if positions.ndim != 2 or positions.shape[1] != 3:
            raise ValueError("positions must have shape (n_atoms, 3).")
        if not np.all(np.isfinite(positions)):
            raise ValueError("positions must contain only finite values.")

        atom_ids = tuple(self.atom_ids)
        if len(atom_ids) != positions.shape[0]:
            raise ValueError("atom_ids must have one entry per position.")
        for atom_id in atom_ids:
            if not isinstance(atom_id, str) or not atom_id.strip():
                raise ValueError("atom_ids must be non-empty strings.")
        if len(set(atom_ids)) != len(atom_ids):
            raise ValueError("atom_ids must be unique.")

        if self.score is not None:
            score = float(self.score)
            if not np.isfinite(score):
                raise ValueError("score must be finite when provided.")
            if not self.score_unit:
                raise ValueError("score_unit is required when a score is provided.")
            if not self.score_source:
                raise ValueError("score_source is required when a score is provided.")
            object.__setattr__(self, "score", score)

        object.__setattr__(self, "positions", positions.copy())
        object.__setattr__(self, "atom_ids", atom_ids)
        object.__setattr__(self, "metadata", dict(self.metadata))

    @property
    def n_atoms(self) -> int:
        """Number of ligand atoms in the pose."""

        return int(self.positions.shape[0])
