# Moleculab

An educational molecular physics and simulation library, built from first
principles.  The goal is not to replace GROMACS, OpenMM, or ORCA.  The goal is
to make the physics behind computational drug discovery visible: every
important equation has an implementation, a unit test, and a numerical
experiment.

Moleculab is aimed at a new computational drug discovery (CADD) user who wants
to understand what a molecular simulation engine actually does before treating
it as a black box.

## Repository layout

The source follows the architecture layers:

| Folder | Contents |
| --- | --- |
| `src/moleculab/foundation/` | units, validation, particles, box geometry |
| `src/moleculab/models/` | potentials, bonded terms, topology, nonbonded, Ewald, physics |
| `src/moleculab/composition/` | `System`, neighbour screening, constraints |
| `src/moleculab/algorithms/` | Velocity Verlet, minimisation |
| `src/moleculab/protocols/` | thermostat, barostat |
| `src/moleculab/analysis/` | MSD, RDF |
| `src/moleculab/applications/` | docking, solvent, water, free energy, QM/MM |
| `src/moleculab/backend.py` | optional CPU/ROCm/PyTorch dispatch |

The public API remains the flat `from moleculab import ...` interface.

## Current state

Moleculab now implements a small but coherent molecular-mechanics engine:

- A consistent internal unit system (nm, ps, kJ/mol, Da).
- `Particle`: mass, position, and velocity in three dimensions.
- `HarmonicPotential`: an external harmonic potential.
- Nonbonded pair potentials: Lennard-Jones and Coulomb.
- `System`: many particles, summed nonbonded interactions, bonded topology,
  and an optional periodic box.
- Bonded topology: harmonic bonds, harmonic angles, and periodic torsions.
- General physics functions: kinetic energy, potential energy, total energy,
  and acceleration.
- Velocity Verlet integration for one particle or a full system.
- Steepest-descent energy minimization with backtracking.
- A Berendsen thermostat and kinetic-temperature calculation.
- Orthorhombic periodic boundaries with the minimum-image convention.
- Trajectory analysis: mean-squared displacement and radial distribution
  function, including periodic minimum-image distances.
- Explicit coarse-grained solvent boxes with lattice/random placement,
  density-based sizing, optional thermal velocities, and thermostat support.
- Free-energy perturbation on an analytically solvable harmonic system, with
  standard errors and convergence examples.
- Per-atom types and charges, typed Lennard-Jones parameters, automatic 1-2/1-3
  intramolecular exclusions, and reaction-field electrostatics.
- A three-site water model with SHAKE/RATTLE rigid constraints and a periodic
  water-box builder.  Flexible water is available with `rigid=False`, and the
  builder can use reaction-field or Ewald electrostatics.
- Rigid-ligand docking with translation/rotation sampling and typed
  Lennard-Jones/Coulomb scoring.
- Ewald summation for periodic typed charges, integrated with typed
  nonbonded interactions and `System` forces, including an intramolecular
  exclusion correction for bonded systems.
- Pressure and Berendsen pressure coupling/NPT with kinetic and virial
  contributions and isotropic box/position scaling.
- Thermodynamic integration for an analytic λ-coupled harmonic system, with
  Simpson/trapezoid quadrature and propagated standard errors.
- MBAR over multiple harmonic alchemical windows, validated against exact
  free-energy differences.
- Molecular FEP/TI/MBAR on a soft-core Lennard-Jones pair, validated against
  high-accuracy radial quadrature.
- Multi-particle Monte Carlo FEP/TI/MBAR for a bonded Lennard-Jones pair.
- Solvated-ligand soft-core decoupling in an explicit Ewald water box, with
  FEP/TI/MBAR cross-validation.
- Flexible-receptor docking with explicit side-chain sampling and a harmonic
  receptor-strain calibration term.
- A minimal finite-difference QM/MM particle-in-a-box model coupled to fixed
  MM point charges, validated against analytic and perturbative references.
- Acid-base, protonation, and tautomer population bookkeeping for discrete
  chemical states.
- Parser-independent binding-site geometry: minimum protein-ligand distances,
  cutoff selection, residue aggregation, and pocket centre/bounding box.
- Toolkit-independent interaction features and hydrogen-bond geometry, with
  explicit distance and D-H...A angle criteria.
- Interaction fingerprints, deterministic feature schemas, and sampled-frame
  occupancy with missing-frame handling.
- External docking pose records, explicit atom mapping, symmetry-aware ligand
  RMSD, and pose-recovery thresholds.
- Pose clustering and binding-mode comparison, separating geometry from
  residue-level interaction patterns.
- Reproducible docking benchmark records with valid-case success rates and
  explicit invalid-case reporting.
- CPU-default NumPy execution with an optional lazy PyTorch backend for the
  Ewald reciprocal energy/force path.
- Tests that check physical relationships, not just that code runs.

This is still an educational engine, not a production molecular-dynamics
package.  Solvent models, free-energy methods, docking, and quantum chemistry
remain on the roadmap.

The longer roadmap lives in [`docs/lib.md`](docs/lib.md).  The as-built
architecture is described in [`docs/architecture.md`](docs/architecture.md).

The tutorial entry point is
[`docs/tutorials/README.md`](docs/tutorials/README.md), which maps the three
connected tracks:

- Track A, the physics engine:
  [`tutorial.md`](docs/tutorials/tutorial.md) and the
  [`physics/`](docs/tutorials/physics/README.md) deep dives;
- Track B, chemical structure and state:
  [`chemistry/`](docs/tutorials/chemistry/README.md);
- Track C, structure-based drug discovery:
  [`drug_discovery/`](docs/tutorials/drug_discovery/README.md).

[`docs/tutorials/introduction.md`](docs/tutorials/introduction.md) explains how
the tracks connect to CADD applications, and
[`docs/tutorials/tracks.md`](docs/tutorials/tracks.md) gives the prerequisites
and full curriculum. The tracks meet in the
[`imatinib–BCR-ABL1 capstone`](docs/tutorials/capstone_imatinib.md). The first
implementation notes are in [`docs/energies.md`](docs/energies.md);
[`docs/validation.md`](docs/validation.md) documents the component,
combined-system, numerical, and reproducibility checks, and
[`docs/acceleration.md`](docs/acceleration.md) covers CPU/ROCm/NPU backend
choices.

Companion guides for the tutorial sequence:
[`troubleshooting.md`](docs/tutorials/troubleshooting.md),
[`visualization_prompts.md`](docs/tutorials/visualization_prompts.md), and
[`answer_keys.md`](docs/tutorials/answer_keys.md).

## Install with uv

Moleculab uses a `uv`-managed virtual environment in `.venv`.

```bash
uv sync
```

This installs the package in editable mode together with the development
dependencies (`pytest`).  The core runtime requires Python 3.14+ and
`numpy>=2.5`.

To run the examples, install the optional plotting dependency:

```bash
uv sync --extra examples
```

CPU execution is the default.  For AMD ROCm, install the optional extra:

```bash
uv sync --extra rocm
```

This resolves `torch` and `triton-rocm` from the PyTorch ROCm 7.2 index
declared in `pyproject.toml`; no NVIDIA packages are pulled in.  The extra is
never installed by a plain `uv sync`.

For NVIDIA/CUDA, install a matching PyTorch wheel manually, for example:

```bash
uv pip install torch --index-url https://download.pytorch.org/whl/cu130
```

Then select it at runtime with `set_backend("auto")` or
`MOLECULAB_BACKEND=auto`.  PyTorch is never required for the normal test suite
or the core library.

Activate the environment if you prefer a shell session:

```bash
source .venv/bin/activate
```

## Run the tests

```bash
uv run pytest -v
```

Run one test or one file while developing:

```bash
uv run pytest tests/test_integrators.py::test_velocity_verlet_conserves_energy_for_a_harmonic_oscillator -v
uv run pytest tests/test_potentials.py -v
```

## CPU and optional acceleration

The default backend is NumPy on the CPU.  Select an accelerator explicitly:

```python
from moleculab import backend_info, set_backend

set_backend("auto")   # use CUDA/ROCm or MPS if PyTorch reports one
print(backend_info())
```

Or set the environment variable:

```bash
MOLECULAB_BACKEND=auto uv run python examples/backend_benchmark.py --backend auto
```

`--backend torch` forces PyTorch even on CPU for benchmarking; `auto` stays on
the CPU when no supported accelerator is visible.  The current accelerated
path is Ewald reciprocal-space energy and forces.  Other modules remain
NumPy/CPU and convert at their boundaries.

For APUs, `set_precision("float32")` enables the faster half-width ROCm path,
and `to_torch`/`to_numpy` use a zero-copy `torch.from_numpy` handoff before the
non-blocking device move.  The CPU thread limits and a vectorized GPU pairwise
example are in `examples/amd_gpu_tuning.py`.

### AMD ROCm on this machine

The system ROCm/HIP runtime and the Radeon 890M (`gfx1151`) are present:

```bash
rocminfo        # enumerates the Radeon 890M GPU agent
hipconfig       # HIP_PATH=/usr, ROCM_PATH=/usr, HIP_PLATFORM=amd
```

The optional ROCm extra resolves the matching PyTorch wheel automatically:

```bash
uv sync --extra rocm
```

If you need to install it manually instead:

```bash
uv pip install --index-url https://download.pytorch.org/whl/rocm7.2 \
  "torch==2.14.1+rocm7.2"
```

The manual command is also useful if your installed ROCm is 6.x; substitute a
matching `rocm6.4` wheel/index.

If the device is not detected, set the architecture override before launching:

```bash
export HSA_OVERRIDE_GFX_VERSION=11.5.1
```

Verified on this system:

```text
Python 3.14.8
numpy 2.5.4
pytest 9.2.0
torch 2.14.1+rocm7.2, HIP 7.2.53211
torch.cuda.is_available() -> True
device -> AMD Radeon 890M
backend_info() -> {'backend': 'torch', 'device': 'cuda', 'accelerator': True}
```

Vendor NPUs (AMD XDNA) are present as `/dev/accel0` but are not supported by
this backend layer.  They require an ONNX Runtime/Vitis AI/vendor runtime path
and are inference-oriented; they are not exposed as a native PyTorch device.

## Try it

```bash
uv run python examples/first_particle.py
uv run python examples/harmonic_oscillator.py
uv run python examples/nonbonded_interactions.py
uv run python examples/bonded_molecule.py
uv run python examples/torsion_scan.py
uv run python examples/minimization.py
uv run python examples/thermostat.py
uv run python examples/periodic_box.py
uv run python examples/trajectory_analysis.py
uv run python examples/solvent_box.py
uv run python examples/free_energy_perturbation.py
uv run python examples/water_box.py
uv run python examples/water_ewald.py
uv run python examples/rigid_water.py
uv run python examples/rigid_docking.py
uv run python examples/ewald.py
uv run python examples/barostat.py
uv run python examples/thermodynamic_integration.py
uv run python examples/mbar.py
uv run python examples/lj_alchemy.py
uv run python examples/mc_alchemy.py
uv run python examples/solvated_alchemy.py
uv run python examples/qm_mm.py
uv run python examples/flexible_docking.py
uv run python examples/backend_benchmark.py --backend cpu
uv run python examples/amd_gpu_tuning.py --backend auto
```

- `first_particle.py` reproduces the hand calculation from the first tutorial.
- `harmonic_oscillator.py` integrates the equations of motion and shows energy
  conservation.
- `nonbonded_interactions.py` compares Lennard-Jones and Coulomb energies and
  forces over a range of distances.
- `bonded_molecule.py` simulates a water-like triatomic with bonds and an
  angle.
- `torsion_scan.py` rotates a four-atom chain and prints the periodic torsion
  energy.
- `algorithms/minimization.py` relaxes a distorted molecule to its equilibrium geometry.
- `thermostat.py` runs a short NVT-like loop with the Berendsen thermostat.
- `periodic_box.py` demonstrates minimum-image distances across a box face.
- `trajectory_analysis.py` simulates a periodic Lennard-Jones fluid and
  computes its MSD and radial distribution function.
- `solvent_box.py` builds an explicit coarse-grained water box, runs a short
  thermostatted trajectory, and analyzes its structure.
- `free_energy_perturbation.py` estimates a harmonic free-energy difference and
  shows how the FEP estimate and standard error converge with more samples.
- `water_box.py` builds a three-site water box with per-atom types, charges,
  exclusions, and reaction-field electrostatics, then runs a short trajectory.
- `water_ewald.py` runs bonded rigid water with Ewald electrostatics and the
  intramolecular exclusion correction.
- `rigid_water.py` demonstrates SHAKE position correction and RATTLE velocity
  correction on a rigid-water trajectory.
- `rigid_docking.py` docks a rigid two-atom ligand into a typed receptor pocket
  and prints the best ranked poses.
- `models/ewald.py` computes Ewald energies for a periodic ionic system and shows
  convergence with the reciprocal-space cutoff.
- `barostat.py` runs a short NPT-like argon simulation and reports pressure,
  temperature, and box volume before and after coupling.
- `thermodynamic_integration.py` integrates a harmonic alchemical path and
  shows convergence with the number of λ windows.
- `applications/mbar.py` combines samples from several harmonic alchemical states and
  compares MBAR free energies with the exact result.
- `applications/lj_alchemy.py` compares FEP, TI, and MBAR for a soft-core Lennard-Jones pair
  against a high-accuracy radial-quadrature reference.
- `applications/mc_alchemy.py` runs actual Metropolis MC on a bonded two-particle System and
  compares FEP, TI, and MBAR with the same reference.
- `applications/solvated_alchemy.py` decouples a neutral ligand from an Ewald water box and
  cross-validates FEP, TI, and MBAR.
- `applications/flexible_docking.py` relaxes a clashing receptor side chain and compares the
  calibrated flexible score with rigid-receptor docking.
- `applications/qm_mm.py` solves a finite-difference electron-in-a-box model and perturbs it
  with fixed MM point charges.
- `backend_benchmark.py` benchmarks CPU versus optional PyTorch Ewald
  reciprocal energy/forces.
- `amd_gpu_tuning.py` exercises zero-copy handoff, float32 ROCm Ewald,
  vectorized GPU pairwise distances, bfloat16 autocast, and optional
  `torch.compile`.

A minimal interactive example:

```python
from moleculab import HarmonicPotential, Particle, total_energy

particle = Particle(mass=12.0, position=[0.1, 0.0, 0.0], velocity=[0.1, 0.0, 0.0])
potential = HarmonicPotential(stiffness=100.0, center=[0.0, 0.0, 0.0])

print(total_energy(particle, potential))  # 0.56 kJ/mol
```

## Units

Moleculab uses one internal unit system throughout:

| Quantity    | Unit         |
| ----------- | ------------ |
| Distance    | nm           |
| Time        | ps           |
| Energy      | kJ/mol       |
| Mass        | Da           |
| Force       | kJ/(mol nm)  |
| Temperature | K            |
| Charge      | e            |

With these units, Newton's second law is numerically `a = F / m` without an
extra conversion factor.  This is a deliberate simplification for learning and
is why the tests use forces, masses, and accelerations directly.  The Coulomb
constant is likewise expressed in internal units as
`COULOMB_CONSTANT ≈ 138.935456 kJ nm/(mol e²)`.

## Design principles

1. **Separate physical models from numerical algorithms.**  A potential
   describes an interaction.  A physics function applies a general law.  An
   integrator advances time.  None of them should absorb the others'
   responsibilities.
2. **Energy and force must agree.**  Every potential implements `energy()` and
   `force()`, and tests compare the analytical force with a numerical
   derivative of the energy.
3. **Tests verify physics.**  A function that runs is not necessarily correct.
   The test suite checks equations, units, conservation laws, and limiting
   cases.
4. **Grow in small, tested milestones.**  Add one physical concept at a time,
   understand it, then move on.  The roadmap deliberately delays complex force
   fields and quantum chemistry until the foundations are solid.

## Package layout

```text
src/moleculab/
├── units.py        Unit labels, constants, and unit conversions
├── backend.py      CPU/PyTorch backend selection
├── particles.py    Particle state and validation
├── potentials.py   Nonbonded models: harmonic, Lennard-Jones, Coulomb
├── nonbonded.py    Typed nonbonded model, exclusions, reaction field
├── ewald.py        Periodic Ewald electrostatics
├── bonded.py       Bond, angle, and torsion potential models
├── physics.py      General energy and force equations
├── topology.py     Bonds, angles, torsions, and topology validation
├── system.py       Particles plus interactions and optional periodic box
├── box.py          Orthorhombic periodic boundaries and minimum images
├── integrators.py  Velocity Verlet and trajectory containers
├── minimization.py Steepest-descent energy minimization
├── thermostats.py  Kinetic temperature and the Berendsen thermostat
├── pressure.py     Pressure, virial, and the Berendsen barostat
├── analysis.py     Mean-squared displacement and radial distribution
├── solvent.py      Coarse-grained explicit solvent boxes
├── water.py        Three-site water model
├── constraints.py  SHAKE/RATTLE distance constraints
├── fep.py          Free-energy perturbation for harmonic systems
├── ti.py           Thermodynamic integration for harmonic systems
├── mbar.py         Multistate Bennett acceptance ratio
├── lj_alchemy.py   Soft-core Lennard-Jones FEP/TI/MBAR path
├── mc_alchemy.py   Multi-particle MC FEP/TI/MBAR path
├── solvated_alchemy.py  Solvated-ligand Ewald decoupling
├── qm_mm.py        Finite-difference QM/MM particle in a box
├── flexible_docking.py  Side-chain flexible docking
└── docking.py      Rigid-ligand translation/rotation docking
```

## Roadmap

Moleculab is being built from simple mechanics toward molecular simulation:

- [x] Units and numerical conventions
- [x] Atoms, coordinates, and particle state
- [x] Harmonic potential and energy/force consistency
- [x] Lennard-Jones and Coulomb interactions
- [x] Molecular topology, bonds, angles, and torsions
- [x] Energy minimization and Velocity Verlet dynamics
- [x] Temperature control with a Berendsen thermostat
- [x] Orthorhombic periodic boundaries
- [x] Trajectory analysis: mean-squared displacement and RDF
- [x] Coarse-grained explicit solvent boxes
- [x] Free-energy perturbation on an analytic harmonic system
- [x] Atomistic three-site water with per-atom types and exclusions
- [x] SHAKE/RATTLE rigid-water constraints
- [x] Reaction-field electrostatics
- [x] Rigid-ligand docking with typed nonbonded scoring
- [x] Ewald summation for periodic point charges
- [x] Ewald intramolecular-exclusion correction for bonded systems
- [x] Pressure and Berendsen NPT coupling
- [x] Thermodynamic integration on an analytic harmonic system
- [x] MBAR over analytic harmonic alchemical windows
- [x] Molecular FEP/TI/MBAR on a trapped soft-core Lennard-Jones pair
- [x] Multi-particle Monte Carlo FEP/TI/MBAR for a bonded LJ pair
- [x] Solvated-ligand FEP/TI/MBAR with Ewald water
- [x] Flexible-receptor docking with side-chain strain scoring
- [x] Minimal finite-difference QM/MM particle in a box
- [x] CPU-default backend with optional PyTorch Ewald acceleration
- [ ] Molecular FEP/TI/MBAR for complex systems
- [ ] Trained/empirical binding-affinity scoring functions
- [ ] Real QM/MM, Hartree-Fock, and DFT

The roadmap is a guide, not a race.  Each milestone should have a working
implementation, tests, and a numerical experiment before the next one begins.
