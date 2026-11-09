"""Tests for CPU/default and optional PyTorch acceleration backends."""

import numpy as np
import pytest

from moleculab import (
    available_backends,
    backend_info,
    get_backend,
    get_precision,
    set_backend,
    set_precision,
    to_numpy,
    to_torch,
    torch_available,
    use_backend,
)
from moleculab.backend import _load_torch
from moleculab.box import Box
from moleculab.ewald import EwaldElectrostatics


def ionic_system():
    box = Box([2.0, 2.0, 2.0])
    positions = np.array([[0.0, 0.0, 0.0], [0.7, 0.2, 0.1]])
    charges = np.array([1.0, -1.0])
    return box, positions, charges


def test_cpu_is_the_default_backend():
    set_backend("cpu")

    assert get_backend() == "cpu"
    assert "cpu" in available_backends()
    info = backend_info()
    assert info["backend"] == "cpu"
    assert isinstance(info["torch_available"], bool)


def test_backend_context_manager_restores_previous_backend():
    set_backend("cpu")

    with use_backend("cpu") as selected:
        assert selected == "cpu"
        assert get_backend() == "cpu"

    assert get_backend() == "cpu"


def test_auto_backend_selects_available_hardware():
    set_backend("cpu")

    selected = set_backend("auto")

    assert selected in {"cpu", "torch"}
    if not torch_available():
        assert selected == "cpu"
    set_backend("cpu")


def test_invalid_backend_is_rejected():
    with pytest.raises(ValueError, match="Backend"):
        set_backend("quantum")


def test_precision_selection_and_validation():
    set_precision("float64")
    assert get_precision() == "float64"
    assert backend_info()["precision"] == "float64"
    assert set_precision("float32") == "float32"
    assert get_precision() == "float32"
    with pytest.raises(ValueError, match="Precision"):
        set_precision("float16")
    set_precision("float64")


def test_zero_copy_cpu_tensor_when_available():
    if not torch_available():
        pytest.skip("PyTorch is not installed")
    import numpy as np
    from moleculab.backend import _load_torch

    torch = _load_torch()
    array = np.ascontiguousarray(np.arange(12, dtype=np.float64).reshape(3, 4))
    set_precision("float64")
    tensor = to_torch(array, device="cpu")

    assert tensor.data_ptr() == array.ctypes.data
    tensor[0, 0] = 99.0
    assert array[0, 0] == 99.0
    np.testing.assert_allclose(to_numpy(tensor), array)


def test_float32_ewald_close_to_float64_when_torch_available():
    if not torch_available():
        pytest.skip("PyTorch is not installed")

    box, positions, charges = ionic_system()
    ewald = EwaldElectrostatics(alpha=6.0, reciprocal_cutoff=40.0)
    try:
        set_backend("torch")
        set_precision("float64")
        energy64 = ewald.reciprocal_energy(positions, charges, box)
        forces64 = ewald.global_forces(positions, charges, box)
        set_precision("float32")
        energy32 = ewald.reciprocal_energy(positions, charges, box)
        forces32 = ewald.global_forces(positions, charges, box)
    finally:
        set_precision("float64")
        set_backend("cpu")

    assert energy32 == pytest.approx(energy64, rel=1e-4)
    np.testing.assert_allclose(forces32, forces64, rtol=1e-3, atol=1e-3)



def test_explicit_torch_requires_installation():
    if torch_available():
        pytest.skip("PyTorch is installed")
    with pytest.raises(ImportError, match="PyTorch"):
        set_backend("torch")


def test_torch_ewald_matches_cpu_when_available():
    if not torch_available():
        pytest.skip("PyTorch is not installed")

    box, positions, charges = ionic_system()
    ewald = EwaldElectrostatics(alpha=6.0, reciprocal_cutoff=40.0)
    set_backend("cpu")
    cpu_energy = ewald.reciprocal_energy(positions, charges, box)
    cpu_forces = ewald.global_forces(positions, charges, box)
    try:
        set_backend("torch")
        torch_energy = ewald.reciprocal_energy(positions, charges, box)
        torch_forces = ewald.global_forces(positions, charges, box)
    finally:
        set_backend("cpu")

    assert torch_energy == pytest.approx(cpu_energy, rel=1e-10, abs=1e-10)
    np.testing.assert_allclose(torch_forces, cpu_forces, rtol=1e-10, atol=1e-10)
