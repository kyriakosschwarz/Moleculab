# Performance and acceleration

Performance work is separate from the scientific model. The CPU/NumPy path is
the default; optional acceleration must preserve the same results.

## Where to look

- [Acceleration backends](acceleration.md) — CPU default, optional PyTorch
  Ewald path, ROCm setup, and NPU limitations.
- `examples/performance/backend_benchmark.py` — CPU versus optional backend timing and
  checksum comparison.
- `examples/performance/amd_gpu_tuning.py` — ROCm-oriented zero-copy, vectorized, and
  precision experiments.

## Rules

1. Profile before changing algorithms.
2. Prefer algorithmic reductions, such as the cell list, over micro-optimising
   Python loops.
3. Keep device-specific code in `backend.py` and `models/ewald.py`, not in
   potentials or integrators.
4. Validate accelerated paths against the CPU path.
5. Do not claim NPU support that the backend does not provide.

The CPU path remains the reference implementation.
