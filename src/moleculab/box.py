"""Periodic boundary conditions for an orthorhombic simulation box.

A finite simulation box is usually treated as if it were surrounded by copies
of itself.  When a particle leaves through one face, its periodic image enters
through the opposite face.  The **minimum-image convention** then computes the
shortest separation between two particles, including the images:

.. math:: \\mathbf d_{\\min} = \\mathbf d
          - \\mathrm{round}(\\mathbf d / \\mathbf L) \\odot \\mathbf L.

This module currently supports orthorhombic boxes (three independent edge
lengths).  Triclinic boxes and pressure coupling are deliberately left for a
later milestone.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ._validation import as_vector3


@dataclass
class Box:
    """An orthorhombic periodic box with edge lengths in nm.

    Parameters
    ----------
    lengths:
        Three positive edge lengths ``(Lx, Ly, Lz)`` in nm.
    """

    lengths: np.ndarray

    def __post_init__(self) -> None:
        self.lengths = as_vector3(self.lengths, "lengths").copy()
        if np.any(self.lengths <= 0.0):
            raise ValueError("Box lengths must all be positive.")

    @property
    def volume(self) -> float:
        """Box volume in nm^3."""

        return float(np.prod(self.lengths))

    def minimum_image(self, displacement: np.ndarray) -> np.ndarray:
        """Return the shortest periodic image of *displacement*."""

        return self.minimum_image_array(as_vector3(displacement, "displacement"))

    def minimum_image_array(self, displacements: np.ndarray) -> np.ndarray:
        """Return minimum-image vectors for an array of displacements.

        The last axis must have length 3.  Any number of leading axes is
        allowed, which makes this convenient for pair-distance matrices.
        """

        vectors = np.asarray(displacements, dtype=float)
        if vectors.ndim == 0 or vectors.shape[-1] != 3:
            raise ValueError(
                "displacements must have shape (..., 3), "
                f"got {vectors.shape}."
            )
        if not np.all(np.isfinite(vectors)):
            raise ValueError("displacements must contain only finite values.")
        return vectors - np.round(vectors / self.lengths) * self.lengths

    def wrap(self, position: np.ndarray) -> np.ndarray:
        """Wrap a position into the primary box ``[0, L)``."""

        vector = as_vector3(position, "position")
        return vector - np.floor(vector / self.lengths) * self.lengths

    def distance(self, displacement: np.ndarray) -> float:
        """Return the minimum-image distance in nm."""

        return float(np.linalg.norm(self.minimum_image(displacement)))

    def __repr__(self) -> str:
        lengths = np.array2string(self.lengths, precision=3)
        return f"Box(lengths={lengths} nm)"
