"""Optional compute backends.

Moleculab is NumPy/CPU by default.  The backend selector is deliberately
lazy: PyTorch is imported only when the user explicitly selects it or asks
for ``"auto"``.

The current accelerated path is Ewald reciprocal-space energy/forces.  Other
modules continue to run on the CPU and convert as needed.  This keeps the
default dependency light and avoids making the educational code depend on a
large optional framework.

Examples
--------
>>> from moleculab.backend import set_backend, backend_info
>>> set_backend("cpu")
'cpu'
>>> backend_info()["backend"]
'cpu'
"""

from __future__ import annotations

import os
from collections.abc import Iterator
from contextlib import contextmanager

_BACKEND = os.environ.get("MOLECULAB_BACKEND", "cpu").strip().lower()
_PRECISION = os.environ.get("MOLECULAB_PRECISION", "float64").strip().lower()
_TORCH = None


def _load_torch():
    global _TORCH
    if _TORCH is not None:
        return _TORCH
    try:
        import torch
    except ImportError as exc:
        raise ImportError(
            "PyTorch is not installed. Install the optional acceleration extra "
            "with `uv sync --extra rocm` or `pip install torch`."
        ) from exc
    _TORCH = torch
    return torch


def torch_available() -> bool:
    try:
        _load_torch()
    except ImportError:
        return False
    return True


def available_backends() -> tuple[str, ...]:
    """Return the backend names usable in this environment."""

    return ("cpu", "torch") if torch_available() else ("cpu",)


def torch_device() -> str:
    """Return the best PyTorch device string for this machine."""

    torch = _load_torch()
    if torch.cuda.is_available():
        return "cuda"
    if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        return "mps"
    return "cpu"


def _auto_backend() -> str:
    if not torch_available():
        return "cpu"
    torch = _load_torch()
    if torch.cuda.is_available():
        return "torch"
    if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        return "torch"
    return "cpu"


def set_backend(name: str = "cpu") -> str:
    """Select ``"cpu"``, ``"torch"``, or ``"auto"``.

    ``"auto"`` selects PyTorch only when an accelerator (CUDA/ROCm or Apple
    MPS) is available; otherwise it stays on the CPU.  PyTorch on CPU remains
    available explicitly as ``"torch"`` for testing and benchmarking.
    """

    global _BACKEND
    value = str(name).strip().lower()
    if value == "cpu":
        _BACKEND = "cpu"
    elif value == "auto":
        _BACKEND = _auto_backend()
    elif value == "torch":
        _load_torch()
        _BACKEND = "torch"
    else:
        raise ValueError("Backend must be 'cpu', 'torch', or 'auto'.")
    return _BACKEND


def get_backend() -> str:
    """Return the currently selected backend name."""

    return _BACKEND


def set_precision(name: str = "float64") -> str:
    """Select the PyTorch floating-point precision for accelerated kernels.

    ``"float64"`` is the accuracy-oriented default.  ``"float32"`` is
    substantially faster on RDNA/APU hardware and is usually sufficient for
    MD-scale Ewald calculations.  It is opt-in so correctness tests remain
    strict by default.
    """

    global _PRECISION
    value = str(name).strip().lower()
    if value not in {"float64", "float32"}:
        raise ValueError("Precision must be 'float64' or 'float32'.")
    _PRECISION = value
    return _PRECISION


def get_precision() -> str:
    return _PRECISION


def torch_dtype():
    """Return the selected PyTorch floating-point dtype."""

    torch = _load_torch()
    if _PRECISION == "float32":
        return torch.float32
    return torch.float64


def to_torch(
    array,
    *,
    device: str | None = None,
    dtype=None,
    copy: bool = False,
):
    """Convert a NumPy array to a PyTorch tensor.

    With ``copy=False`` this uses ``torch.from_numpy`` first, which is a
    zero-copy view when the array is C-contiguous and already has the right
    dtype.  It then moves to the selected device with ``non_blocking=True``.
    On APUs the CPU and GPU share physical RAM, so avoiding an intermediate
    NumPy copy is especially useful.
    """

    import numpy as np

    torch = _load_torch()
    values = np.asarray(array)
    if copy:
        tensor = torch.as_tensor(
            values,
            dtype=dtype if dtype is not None else torch_dtype(),
        )
    else:
        contiguous = np.ascontiguousarray(values)
        tensor = torch.from_numpy(contiguous)
        target_dtype = dtype if dtype is not None else torch_dtype()
        if tensor.dtype != target_dtype:
            tensor = tensor.to(target_dtype)
    target_device = device if device is not None else (
        torch_device() if get_backend() == "torch" else "cpu"
    )
    return tensor.to(target_device, non_blocking=True)


def to_numpy(tensor):
    """Return a detached CPU NumPy view/copy of a PyTorch tensor."""

    return tensor.detach().cpu().numpy()


def backend_info() -> dict[str, object]:
    """Return a small diagnostic dictionary."""

    info: dict[str, object] = {
        "backend": _BACKEND,
        "precision": _PRECISION,
        "torch_available": torch_available(),
        "torch_version": None,
        "device": None,
        "accelerator": False,
    }
    if info["torch_available"]:
        torch = _load_torch()
        device = torch_device()
        info["torch_version"] = getattr(torch, "__version__", "unknown")
        info["device"] = device
        info["accelerator"] = device != "cpu"
    return info


@contextmanager
def use_backend(name: str = "cpu") -> Iterator[str]:
    """Temporarily select a backend and restore the previous one."""

    previous = get_backend()
    selected = set_backend(name)
    try:
        yield selected
    finally:
        set_backend(previous)
