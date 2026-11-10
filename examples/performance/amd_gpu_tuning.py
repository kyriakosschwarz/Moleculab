# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Tune Molab for an AMD Strix Point APU with ROCm PyTorch.

This example demonstrates the practical optimizations for APUs:

* zero-copy NumPy -> PyTorch handoff on CPU;
* float32 Ewald on ROCm for higher throughput;
* vectorized pairwise distances entirely on the GPU;
* optional bfloat16 autocast for matmul;
* CPU BLAS thread limits for the heterogeneous Zen 5 / Zen 5c layout;
* optional ``torch.compile``.

Run with the optional ROCm extra installed:

    uv run python examples/performance/amd_gpu_tuning.py --backend auto
    uv run python examples/performance/amd_gpu_tuning.py --backend auto --compile

CPU remains the default for the library; this example explicitly selects the
GPU when available.
"""

from __future__ import annotations

import os

# These must be set before NumPy/OpenBLAS/MKL initialize their thread pools.
os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("MKL_NUM_THREADS", "8")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "8")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "8")

import argparse
import time

import numpy as np

from moleculab import (
    EwaldElectrostatics,
    backend_info,
    set_backend,
    set_precision,
    to_torch,
    torch_available,
)
from moleculab.backend import _load_torch, torch_device
from moleculab.foundation.box import Box


def ionic_system():
    box_length = 4.0
    positions = []
    charges = []
    for i in range(2):
        for j in range(2):
            for k in range(2):
                positions.append(
                    [
                        (i + 0.5) * box_length / 2.0,
                        (j + 0.5) * box_length / 2.0,
                        (k + 0.5) * box_length / 2.0,
                    ]
                )
                charges.append(1.0 if (i + j + k) % 2 == 0 else -1.0)
    coordinates = np.asarray(positions)
    coordinates[1, 0] += 0.07
    return coordinates, np.asarray(charges), Box([box_length, box_length, box_length])


def time_ewald(ewald, positions, charges, box, repeats):
    ewald.global_energy(positions, charges, box)
    ewald.global_forces(positions, charges, box)
    start = time.perf_counter()
    for _ in range(repeats):
        energy = ewald.global_energy(positions, charges, box)
        forces = ewald.global_forces(positions, charges, box)
    elapsed = time.perf_counter() - start
    return energy, forces, elapsed


def pairwise_distances(coords):
    torch = _load_torch()
    diff = coords.unsqueeze(1) - coords.unsqueeze(0)
    return torch.linalg.norm(diff, dim=-1)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--backend",
        choices=("cpu", "auto", "torch"),
        default="auto",
    )
    parser.add_argument("--repeats", type=int, default=20)
    parser.add_argument("--compile", action="store_true")
    args = parser.parse_args()

    selected = set_backend(args.backend)
    print("backend:", backend_info())
    print("selected:", selected)
    print("CPU thread limits:", {
        key: os.environ[key]
        for key in (
            "OMP_NUM_THREADS",
            "MKL_NUM_THREADS",
            "OPENBLAS_NUM_THREADS",
        )
    })
    print()

    positions, charges, box = ionic_system()
    ewald = EwaldElectrostatics(alpha=4.0, reciprocal_cutoff=30.0)

    set_precision("float64")
    set_backend("cpu")
    energy_cpu, forces_cpu, cpu_time = time_ewald(
        ewald, positions, charges, box, args.repeats
    )
    print(f"CPU float64: {cpu_time:.6f} s for {args.repeats} energy+force calls")
    print(f"  energy: {energy_cpu:.8f}, force norm: {np.linalg.norm(forces_cpu):.8f}")
    print()

    if not torch_available():
        print("PyTorch is not installed; install with `uv sync --extra rocm`.")
        return
    torch = _load_torch()
    device = torch_device()
    print("PyTorch device:", device)
    if device == "cuda":
        print("GPU:", torch.cuda.get_device_name(0))

    set_backend("torch")
    for precision in ("float64", "float32"):
        set_precision(precision)
        energy, forces, elapsed = time_ewald(
            ewald, positions, charges, box, args.repeats
        )
        print(
            f"ROCm {precision}: {elapsed:.6f} s, "
            f"energy error vs CPU {abs(energy - energy_cpu):.3e}, "
            f"force norm {np.linalg.norm(forces):.8f}"
        )
    set_precision("float32")

    # Zero-copy CPU handoff: torch.from_numpy shares the NumPy buffer.
    array = np.ascontiguousarray(np.arange(12, dtype=np.float32).reshape(3, 4))
    cpu_tensor = to_torch(array, device="cpu", dtype=torch.float32)
    print()
    print("zero-copy CPU buffer shared:", cpu_tensor.data_ptr() == array.ctypes.data)

    # Vectorized O(N^2) pairwise distances on the GPU.
    coords = to_torch(
        np.random.default_rng(0).normal(size=(2048, 3)).astype(np.float32),
        dtype=torch.float32,
    )
    if device == "cuda":
        torch.cuda.synchronize()
    start = time.perf_counter()
    distances = pairwise_distances(coords)
    if device == "cuda":
        torch.cuda.synchronize()
    print(f"GPU pairwise distances 2048 atoms: {time.perf_counter() - start:.6f} s")
    print("distance matrix:", tuple(distances.shape), "device:", distances.device)

    # bfloat16 autocast is useful for bandwidth-bound matmul on RDNA 3.5.
    a = torch.randn(2048, 2048, device=device, dtype=torch.float32)
    b = torch.randn(2048, 2048, device=device, dtype=torch.float32)
    with torch.autocast(device_type="cuda" if device == "cuda" else "cpu", dtype=torch.bfloat16):
        product = a @ b
    print("autocast matmul:", tuple(product.shape), product.dtype)

    if args.compile:
        try:
            compiled = torch.compile(pairwise_distances, backend="inductor")
            distances = compiled(coords)
            print("torch.compile pairwise distances:", tuple(distances.shape))
        except Exception as exc:
            print("torch.compile unavailable:", type(exc).__name__, exc)

    print()
    print("CPU remains Molab's default; select `auto`/`torch` explicitly for ROCm.")
    print("NPUs are not exposed as native PyTorch devices by this backend layer.")


if __name__ == "__main__":
    main()
