# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Benchmark CPU versus optional PyTorch acceleration for Ewald forces.

CPU is the default.  Select ``--backend auto`` to use a CUDA/ROCm or Apple MPS
device when PyTorch is installed and reports one; ``--backend torch`` forces
PyTorch even on CPU.  Without an accelerator, auto stays on CPU.

Run it with:

    uv run python examples/performance/backend_benchmark.py --backend cpu
    uv run python examples/performance/backend_benchmark.py --backend auto

Install the optional framework first with:

    uv sync --extra rocm
"""

import argparse
import time

import numpy as np

from moleculab import (
    EwaldElectrostatics,
    backend_info,
    set_backend,
    torch_available,
    use_backend,
)
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
    return (
        coordinates,
        np.asarray(charges),
        Box([box_length, box_length, box_length]),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--backend",
        choices=("cpu", "auto", "torch"),
        default="cpu",
        help="CPU is the default; auto uses an accelerator if available.",
    )
    parser.add_argument("--repeats", type=int, default=50)
    args = parser.parse_args()

    try:
        selected = set_backend(args.backend)
    except ImportError as exc:
        print(exc)
        print("Continuing with the CPU default.")
        selected = set_backend("cpu")

    print("Backend diagnostics:", backend_info())
    print(f"Selected backend: {selected}")
    print()

    positions, charges, box = ionic_system()
    ewald = EwaldElectrostatics(alpha=4.0, reciprocal_cutoff=30.0)

    backends = ["cpu"]
    if torch_available() and selected != "cpu":
        backends.append("torch")

    for name in backends:
        with use_backend(name):
            # Warm-up
            ewald.global_energy(positions, charges, box)
            ewald.global_forces(positions, charges, box)

            start = time.perf_counter()
            for _ in range(args.repeats):
                energy = ewald.global_energy(positions, charges, box)
                forces = ewald.global_forces(positions, charges, box)
            elapsed = time.perf_counter() - start

            print(f"{name:>6} backend")
            print(f"  repeats: {args.repeats}")
            print(f"  total time: {elapsed:.6f} s")
            print(f"  time per energy+force: {elapsed / args.repeats:.6f} s")
            print(f"  reciprocal energy: {energy:.8f} kJ/mol")
            print(f"  force norm: {np.linalg.norm(forces):.8f}")
            print()

    print("CPU is always available. PyTorch is optional and lazy-loaded.")
    print("Current accelerated path: Ewald reciprocal energy and forces.")
    print("Other modules remain NumPy/CPU unless they explicitly support a backend.")
    print("AMD ROCm devices need a ROCm-enabled PyTorch build; vendor NPUs are")
    print("not supported by this backend layer yet.")


if __name__ == "__main__":
    main()
