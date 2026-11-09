# Acceleration backends and native NPU support

## Default policy

Moleculab is CPU/NumPy by default. The core install is:

```bash
uv sync
```

The heavy framework stack is optional:

```bash
uv sync --extra rocm     # AMD ROCm 7.2 PyTorch, optional
uv sync --extra examples # Matplotlib only
```

`torch` and `triton-rocm` live only in the `rocm` extra and are resolved from
the ROCm 7.2 index declared in `pyproject.toml`. A normal `uv sync` does not
install them.

## What is accelerated now

The Ewald reciprocal-space energy/force path has a PyTorch implementation.
CPU remains the fallback. Select it with:

```python
from moleculab import set_backend, set_precision

set_backend("auto")      # CPU unless CUDA/ROCm/MPS is available
set_precision("float32") # faster ROCm path; float64 remains the default
```

Zero-copy CPU handoff and non-blocking device transfer are available through
`moleculab.backend.to_torch` / `to_numpy`. See
`examples/amd_gpu_tuning.py` for float32 Ewald, vectorized GPU pairwise
distances, bfloat16 autocast, and optional `torch.compile`.

## AMD XDNA NPU: why it is not a `torch.device("npu")`

An NPU is not a drop-in GPU. Making it a native PyTorch device would require
all of the following:

1. **Vendor runtime and driver integration** — XDNA/XRT kernel dispatch, queue
   management, memory allocation, and synchronization.
2. **A PyTorch backend** — either a `PrivateUse1` device backend with kernels,
   or a `torch.compile` lowering path. PyTorch's CUDA/HIP paths do not map
   directly onto NPU tile arrays.
3. **Operator coverage and lowering** — MLIR/AIR-based compilation of the
   required operators, including layout and tiling decisions for the NPU.
4. **Autograd and optimizer support** — only if training is required; NPUs are
   primarily static, low-precision inference engines.
5. **dtype and shape constraints** — NPUs favor fixed shapes and INT8/BF16
   workloads, whereas MD/Ewald calculations use dynamic shapes, reductions,
   complex arithmetic, and often FP64.

That is a substantial backend project, not a configuration flag.

## The practical NPU path

The current industry bridge is:

```text
PyTorch/NumPy model -> ONNX export -> ONNX Runtime -> XDNA/Vitis AI execution provider
```

This is primarily an **inference** path. It bypasses the native PyTorch eager
loop and generally does not support the dynamic, complex-valued kernels used by
the molecular simulation core.

## Is it needed or beneficial?

For Moleculab's current workloads: **no**.

- Ewald, QM/MM diagonalization, force loops, and free-energy sampling are
  small-kernel, dynamic, precision-sensitive workloads.
- The Radeon 890M via ROCm is the appropriate accelerator and is already
  working.
- The NPU is better suited to large, static neural-network inference, such as
  future ML potentials or learned scoring functions.

If Molab later adds a neural surrogate model, an optional ONNX Runtime XDNA
provider could make sense for batch inference. It should remain separate from
the core CPU/NumPy and ROCm paths.
