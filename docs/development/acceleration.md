# Acceleration backends and native NPU support

## Default policy

Moleculab is CPU/NumPy by default. The core install is:

```bash
uv sync
```

The heavy framework stack is optional:

```bash
uv sync --extra rocm                    # AMD ROCm 7.2 PyTorch, optional
uv sync --extra examples                # Matplotlib, Pillow, RDKit, and gemmi
uv sync --extra rocm --extra examples   # both
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
`examples/performance/amd_gpu_tuning.py` for float32 Ewald, vectorized GPU pairwise
distances, bfloat16 autocast, and optional `torch.compile`.

## Python-native performance

Not every bottleneck needs Rust. The first performance layer is algorithmic and
array-centric, and it keeps the physics readable.

| Bottleneck | Python-native strategy | Status |
| --- | --- | --- |
| O(N²) pair screening | cell-list neighbour candidates | implemented as an opt-in `System` path |
| repeated particle attribute access | contiguous arrays from `System` | available for analysis and vectorised work |
| expensive Ewald reciprocal kernel | optional isolated PyTorch backend | implemented, CPU fallback |
| pair evaluation itself | still an explicit Python loop | intentionally readable; future work |

### Algorithmic reduction: the cell list

`composition/neighbour.py` bins particles into cells and proposes candidate pairs. The
`System` force/energy path still computes the exact minimum-image distance and
applies the same potential; the cell list only removes pairs that cannot be
within the cutoff. In the 300-particle example the candidate set is about 12.7%
of all pairs and the measured force calculation is roughly 6x faster. The
default O(N²) loop remains the teaching path.

### Array-centric composition

`System` exposes `positions`, `velocities`, `charges`, and `atom_types` as
NumPy arrays. These are rebuilt from the particle objects, so they are best used
at the boundary of a calculation or for vectorised analysis, not as a claim that
the hot pair loop has no Python overhead.

### Backend isolation

The optional PyTorch path is confined to the Ewald reciprocal kernel. Core
potentials, bonded terms, integrators, and analysis remain NumPy/Python. A
backend change should be invisible to the physics.

### Profile before optimising

Use `time.perf_counter` around the force or energy call, and record the
candidate fraction from the neighbour list. A faster cell list does not help if
the force evaluation or analysis is the real bottleneck. Optimise the measured
hot spot, not the most interesting-looking function.

### Remaining bottlenecks

A representative profile (500 particles, typed Lennard-Jones plus reaction
field, cutoff 1.2 nm, 8 nm box) shows the force path is still dominated by
Python-level per-pair work, not by the cell list itself:

| Bottleneck | Impact | Python-native next step |
| --- | --- | --- |
| per-pair validation and array checks | high | internal unchecked displacements inside `System`, while public APIs keep validation |
| Python pair loop and method dispatch | high | optional vectorized pair kernel; keep the explicit loop as the default |
| integrator rebuilds two `System` objects per step | high for long MD | array-level step with one rebuild per stored frame |
| candidate enumeration using Python sets | medium | vectorized cell-pair generation |
| Ewald CPU loops | high for charged periodic systems | vectorize more NumPy reductions; optional PyTorch path already exists |
| RDF builds an N×N distance matrix per frame | medium | cell-list or chunked pair distances |
| minimisation and barostat copy every particle | medium | in-place array path behind the same tests |
| docking scores poses one at a time | medium | batch receptor-ligand pair evaluation |
| dense QM diagonalisation | high for large grids | sparse or iterative eigensolvers |

The rule is:

1. keep the readable default path;
2. add an opt-in fast path with the same tests;
3. benchmark before and after;
4. never let the fast path change the physics.

The cache of Lorentz-Berthelot parameters, pre-extracted `System` arrays,
squared-cutoff screening, and single-validation `Box.minimum_image` are already
in place.

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
