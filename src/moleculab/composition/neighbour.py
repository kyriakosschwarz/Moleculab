"""Optional cell-list neighbour generation.

The default :class:`~moleculab.composition.system.System` pair loop is deliberately
O(N^2): it is short, readable, and appropriate for the small systems used in
the tutorials.  This module adds an **opt-in** spatial cell list for larger
systems.

The cell list only answers one question:

    Which particle pairs might be within a cutoff?

It may return extra candidate pairs whose true separation is larger than the
cutoff.  The caller still computes the exact minimum-image distance and applies
the potential, exclusions, and any global terms.  This keeps the neighbour list
independent of the force field.

For an orthorhombic periodic box, the cutoff must not exceed half the shortest
box edge.  The cell size is at least the cutoff, so checking the 27 surrounding
cells is sufficient.

Typical use::

    neighbours = build_neighbour_list(system.positions, cutoff, system.box)
    for i, j in neighbours:
        displacement = system.box.minimum_image(
            system.positions[i] - system.positions[j]
        )
        if np.linalg.norm(displacement) < cutoff:
            ...
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product

import numpy as np

from ..foundation.box import Box
from ..foundation.validation import as_positive_float


@dataclass(frozen=True)
class NeighbourList:
    """Candidate particle pairs from a spatial cell list.

    Attributes
    ----------
    pairs:
        Sorted unique ``(i, j)`` pairs with ``i < j``.  Every true pair within
        the cutoff is present; extra candidate pairs are allowed.
    n_cells:
        Number of cells along x, y, and z.
    cell_size:
        Cell edge lengths in nm.
    origin:
        Lower corner of the non-periodic bounding box, or ``(0, 0, 0)`` for a
        periodic box.
    periodic:
        Whether the list was built with a periodic box.
    """

    pairs: tuple[tuple[int, int], ...]
    n_cells: tuple[int, int, int]
    cell_size: tuple[float, float, float]
    origin: tuple[float, float, float]
    periodic: bool

    @property
    def n_pairs(self) -> int:
        """Number of candidate pairs."""

        return len(self.pairs)

    def __len__(self) -> int:
        return len(self.pairs)

    def __iter__(self):
        return iter(self.pairs)


def _validate_positions(positions: np.ndarray) -> np.ndarray:
    coordinates = np.asarray(positions, dtype=float)
    if coordinates.ndim != 2 or coordinates.shape[1] != 3:
        raise ValueError(
            "positions must have shape (n_particles, 3), "
            f"got {coordinates.shape}."
        )
    if coordinates.shape[0] == 0:
        raise ValueError("positions must contain at least one particle.")
    if not np.all(np.isfinite(coordinates)):
        raise ValueError("positions must contain only finite values.")
    return coordinates


def _periodic_cell_count(length: float, cutoff: float) -> int:
    return max(1, int(np.floor(length / cutoff + 1e-12)))


def _nonperiodic_cell_count(extent: float, cutoff: float) -> int:
    if extent <= 0.0:
        return 1
    return max(1, int(np.floor(extent / cutoff + 1e-12)))


def build_neighbour_list(
    positions: np.ndarray,
    cutoff: float,
    box: Box | None = None,
) -> NeighbourList:
    """Return candidate pairs within *cutoff* using a spatial cell list.

    Parameters
    ----------
    positions:
        Particle positions with shape ``(n_particles, 3)`` in nm.
    cutoff:
        Candidate-pair cutoff in nm.  Must be finite and positive.
    box:
        Optional periodic box.  When omitted, a non-periodic bounding box is
        used.

    Returns
    -------
    NeighbourList
        Candidate pairs.  The list is a superset of the true pairs within the
        cutoff; callers must still test the exact distance.

    Notes
    -----
    The cutoff is compared with the **raw minimum-image distance** by the
    caller.  For a periodic box the cell list requires
    ``cutoff <= min(box.lengths) / 2`` so that adjacent cells are sufficient.
    """

    coordinates = _validate_positions(positions)
    cutoff_value = as_positive_float(cutoff, "cutoff")

    if box is None:
        origin = np.min(coordinates, axis=0)
        extents = np.max(coordinates, axis=0) - origin
        n_cells = tuple(
            _nonperiodic_cell_count(float(extent), cutoff_value)
            for extent in extents
        )
        cell_size = tuple(
            max(float(extent) / count, cutoff_value)
            if count == 1
            else float(extent) / count
            for extent, count in zip(extents, n_cells)
        )
        periodic = False
    else:
        if not isinstance(box, Box):
            raise TypeError("box must be a Box or None.")
        lengths = box.lengths
        if cutoff_value > float(np.min(lengths)) / 2.0 + 1e-12:
            raise ValueError(
                "cutoff must not exceed half the shortest box edge "
                f"({float(np.min(lengths)) / 2.0:.6g} nm)."
            )
        origin = np.zeros(3)
        n_cells = tuple(
            _periodic_cell_count(float(length), cutoff_value)
            for length in lengths
        )
        cell_size = tuple(
            float(length) / count for length, count in zip(lengths, n_cells)
        )
        periodic = True

    # Assign every particle to a cell in one NumPy pass.  Cell assignment is
    # vectorized; candidate enumeration below stays explicit so the screening
    # logic is readable.
    cells: dict[tuple[int, int, int], list[int]] = {}
    origin_array = np.asarray(origin, dtype=float)
    cell_size_array = np.asarray(cell_size, dtype=float)
    n_cells_array = np.asarray(n_cells, dtype=np.int64)
    scaled = (coordinates - origin_array) / cell_size_array
    cell_indices = np.floor(scaled).astype(np.int64)
    if periodic:
        cell_indices %= n_cells_array
    else:
        cell_indices = np.clip(cell_indices, 0, n_cells_array - 1)
    for index, cell in enumerate(map(tuple, cell_indices)):
        cells.setdefault(cell, []).append(index)

    # Because the cell edge is at least the cutoff, only the 27 surrounding
    # cells can contain a true pair.  The set removes duplicates created by
    # periodic wrapping.
    offsets = tuple(product((-1, 0, 1), repeat=3))
    pairs: set[tuple[int, int]] = set()
    for cell, members in cells.items():
        for offset in offsets:
            neighbour = tuple(
                (cell[axis] + offset[axis]) % n_cells[axis]
                if periodic
                else cell[axis] + offset[axis]
                for axis in range(3)
            )
            if not periodic and any(
                neighbour[axis] < 0 or neighbour[axis] >= n_cells[axis]
                for axis in range(3)
            ):
                continue
            other = cells.get(neighbour)
            if other is None:
                continue
            for i in members:
                for j in other:
                    if i < j:
                        pairs.add((i, j))
                    elif j < i:
                        pairs.add((j, i))

    return NeighbourList(
        pairs=tuple(sorted(pairs)),
        n_cells=tuple(int(value) for value in n_cells),
        cell_size=tuple(float(value) for value in cell_size),
        origin=tuple(float(value) for value in origin),
        periodic=periodic,
    )
