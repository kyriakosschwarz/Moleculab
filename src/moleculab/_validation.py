"""Validation helpers shared by the public data classes."""

from __future__ import annotations

import numpy as np


def as_vector3(value: object, name: str) -> np.ndarray:
    """Return *value* as a finite, three-component float array.

    The coordinate system in Moleculab is three-dimensional from the start.
    Failing early on a wrong shape is much easier to debug than a broadcast
    error several equations later.
    """

    try:
        vector = np.asarray(value, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a sequence of three numbers.") from exc

    if vector.shape != (3,):
        raise ValueError(f"{name} must have exactly 3 components, got {vector.shape}.")
    if not np.all(np.isfinite(vector)):
        raise ValueError(f"{name} must contain only finite values.")
    return vector


def as_finite_float(value: object, name: str) -> float:
    """Return *value* as a finite float."""

    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a number.") from exc

    if not np.isfinite(result):
        raise ValueError(f"{name} must be finite.")
    return result


def as_positive_float(value: object, name: str) -> float:
    """Return *value* as a finite, positive float."""

    result = as_finite_float(value, name)
    if result <= 0.0:
        raise ValueError(f"{name} must be positive, got {result}.")
    return result


def as_nonnegative_float(value: object, name: str) -> float:
    """Return *value* as a finite, non-negative float."""

    result = as_finite_float(value, name)
    if result < 0.0:
        raise ValueError(f"{name} must be non-negative, got {result}.")
    return result


def as_displacement(
    value: object,
    name: str = "displacement",
) -> tuple[np.ndarray, float]:
    """Return a validated displacement vector and its length.

    A displacement between two distinct particles must be non-zero; otherwise
    distance-dependent forces are undefined.
    """

    vector = as_vector3(value, name)
    distance = float(np.linalg.norm(vector))
    if distance == 0.0:
        raise ValueError(f"{name} must be non-zero.")
    return vector, distance

