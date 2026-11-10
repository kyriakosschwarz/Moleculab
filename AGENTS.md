# AGENTS.md

## Repository status

- Moleculab is a working **educational** molecular physics/simulation library. Optimize for scientific correctness and readable, inspectable code, not production hardening.
- The code in `src/moleculab/` is now the executable source of truth. `docs/lib.md` describes the long-term architecture/roadmap; `docs/energies.md` preserves the first-milestone derivation and expected numbers. Reconcile docs with code rather than assuming the docs are implemented.
- Implemented: units, typed `Particle` state, `HarmonicPotential`, double-well external potential, Lennard-Jones/Coulomb pair potentials, `System` with pair-indexed potentials, bonds/angles/torsions, Velocity Verlet, analytic validation of the Velocity Verlet harmonic oscillator, an optional cell-list neighbour screen with opt-in `System` parity, SHAKE/RATTLE distance constraints, steepest-descent minimization, Berendsen thermostat, orthorhombic PBC, MSD, RDF, coarse-grained solvent boxes, harmonic FEP/TI/MBAR, molecular soft-core Lennard-Jones FEP/TI/MBAR, multi-particle Monte Carlo FEP/TI/MBAR, solvated-ligand Ewald decoupling, typed nonbonded interactions with 1-2/1-3 exclusions, reaction-field electrostatics, a three-site rigid water model, rigid-ligand docking, flexible-receptor docking with explicit side-chain sampling, Ewald summation with intramolecular-exclusion corrections, pressure/Berendsen NPT coupling, a finite-difference QM/MM particle-in-a-box model, acid-base/protonation and tautomer population bookkeeping, and a CPU-default optional PyTorch Ewald backend.
- Deferred: molecular FEP/TI/MBAR for complex systems, trained/empirical binding-affinity scoring, real QM/MM/DFT, and NPU-specific acceleration. Do not add these without tests and a numerical example.
- Python + NumPy first. Rust/PyO3 is explicitly a later optimization, not part of the current implementation.

## Environment and commands

- The project uses a **uv-managed** `.venv` and `uv.lock`. The development environment currently runs Python 3.14.8, NumPy 2.5.4, pytest 9.2.0, and PyTorch 2.14.1+rocm7.2 on AMD ROCm 7.2/HIP. Prefer current versions when upgrading.

```bash
uv sync                                   # install package + dev deps
uv run pytest -q                          # full suite
uv run pytest tests/test_torsion.py -v     # one test file
uv run pytest tests/test_acid_base.py -v   # protonation populations
uv run pytest tests/test_box.py::test_minimum_image_uses_the_shortest_separation -v
uv run python examples/bonded_molecule.py
uv run python examples/harmonic_oscillator.py
uv run python examples/harmonic_validation.py
uv run python examples/chemical_examples.py
uv run python examples/chemical_structure.py
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
uv run python examples/protonation.py
uv run python examples/flexible_docking.py
uv run python examples/backend_benchmark.py --backend cpu
uv run python examples/neighbour_list.py
uv run python docs/generate_architecture_visuals.py
```

- CPU is the default backend. `moleculab.backend.set_backend("cpu"|"auto"|"torch")` and `MOLECULAB_BACKEND` select it; PyTorch is lazy and optional (`uv sync --extra rocm`). Only the Ewald reciprocal energy/force path is accelerated. `auto` picks a CUDA/ROCm or MPS device only when PyTorch reports one; otherwise CPU. Vendor NPUs are not supported. `set_precision("float32")` enables the faster ROCm path; `to_torch`/`to_numpy` provide zero-copy CPU handoff before non-blocking device moves.
- `pyproject.toml` uses a setuptools `src/` layout, `requires-python >=3.14`, and `numpy>=2.5`. `pytest` is in the `dev` dependency group; `matplotlib` is optional under `examples`; ROCm PyTorch is optional under `rocm`, with `torch`/`triton-rocm` resolved from the PyTorch ROCm 7.2 index.
- Re-export public API from `moleculab/__init__.py`; examples and tests import from `moleculab`, not internal modules.

## Units and constants

- Internal units: distance **nm**, time **ps**, energy **kJ/mol**, mass **Da**, force **kJ/(mol nm)**, temperature **K**, charge **e**.
- With this system, `a = F / m` is numerically correct without conversion factors. Do not add ad-hoc conversion constants inside physics code.
- `BOLTZMANN_CONSTANT = 0.008314462618` is actually the **molar gas constant** `R` in kJ/(mol K), despite the name. Keep the name/value for compatibility; `MOLAR_GAS_CONSTANT` is the descriptive alias, and `BOLTZMANN_CONSTANT_J_PER_K` is the single-particle constant.
- `COULOMB_CONSTANT ≈ 138.935456 kJ nm/(mol e²)` is derived from SI constants in `foundation/units.py`. Charges are in elementary-charge units and distances in nm.

## Architecture rules

- Separate **physical models** from **numerical algorithms**: potentials define interactions; physics functions apply general laws; integrators advance time; thermostats control temperature.
- `Particle` stores mass, position, velocity, `atom_type` (default `"default"`), and `charge` (default `0.0`). All `System` rebuild paths must preserve `atom_type` and `charge`; losing them silently changes a typed nonbonded calculation.
- **Nonbonded** potentials (`HarmonicPotential` is external; `LennardJonesPotential` and `CoulombPotential` are pair potentials) implement `energy`/`force`. Pair displacement is `r_i - r_j`; the force is on particle i and the force on j is its negative.
- `models/nonbonded.py` adds the typed path: `NonbondedModel` combines per-type Lennard-Jones parameters with Lorentz-Berthelot rules, reads charges from particles, and applies `ReactionFieldElectrostatics` or no electrostatics. Its cutoff must not exceed half the shortest box edge. Exclusions come from `Topology.nonbonded_exclusions()` (1-2 bonds and 1-3 angle endpoints) plus any explicit model exclusions.
- **Bonded** potentials live in `models/bonded.py` and are indexed by explicit atom numbers in `models/topology.py`. Bond displacement is `r_j - r_i`; angle vectors are `r_ji = r_i - r_j`, `r_jk = r_k - r_j`; torsion bond vectors are `r_ij = r_j - r_i`, `r_jk = r_k - r_j`, `r_kl = r_l - r_k`.
- Forces must be the negative gradient of the same energy (`F = -∇U`). `tests/test_bonded.py`, `test_pair_potentials.py`, `test_torsion.py`, and the consolidated `tests/test_force_gradients.py` verify this with central finite differences.
- Angle parameters are in **radians** (`angle0`), not degrees. Use `np.deg2rad`.
- `HarmonicPotential` is an **external** one-particle potential, not a chemical bond. Do not use it for bonded topology.
- `System` sums nonbonded `pair_potentials` over all unique pairs and adds `topology` bonded terms. When `System.box` is set, only nonbonded interactions use the minimum-image convention; bonded terms use raw coordinates.
- `analysis/trajectory.py` accepts positions shaped `(n_frames, n_particles, 3)` (or one frame as `(n_particles, 3)`). MSD uses minimum-image displacements when given a box; RDF requires a `Box` for volume normalization and rejects `r_max` beyond half the shortest box edge.
- `applications/solvent.py` provides a **coarse-grained single-site** solvent model, not atomistic water. `build_solvent_box` requires exactly one of `box_length` or `density`, places particles on a lattice or randomly, rejects overlaps below `0.9*sigma`, and can initialize Maxwell-Boltzmann velocities. It returns a regular `System` with a `Box`, so the thermostat/integrator/analysis path is unchanged.
- `applications/fep.py` uses molar units: `RT = BOLTZMANN_CONSTANT * temperature`, not the single-particle `k_B T`. FEP is evaluated in log space for stability, and the harmonic oscillator is the analytic reference (`ΔG = 1.5 RT ln(k1/k0)`). `fep_harmonic_oscillator()` is the fast vectorized path; `fep_between_potentials()` is the generic path for arbitrary `energy(position)` objects.
- `applications/ti.py` implements thermodynamic integration for the linear path `U_lambda = (1-lambda)U0 + lambda U1`. `thermodynamic_integration_harmonic()` samples the exact λ-dependent Gaussian and integrates `<U1-U0>` with Simpson/trapezoid; `integrate_thermodynamic_derivative()` is the generic quadrature helper. Simpson requires an odd number of uniformly spaced nodes. Standard errors are propagated from per-window errors as if windows were independent.
- `applications/mbar.py` implements MBAR with reduced potentials `U_k/(RT)` and explicit sample-state labels. `mbar_harmonic()` samples each harmonic state and runs the solver. The self-consistent equations are iterated in log space with damping and gauge `f[0]=0`; no MBAR uncertainty/overlap matrix is implemented yet.
- `applications/lj_alchemy.py` adds a trapped soft-core Lennard-Jones pair. It samples the exact radial Boltzmann distribution by inverse CDF, provides FEP/TI/MBAR estimators, and validates against high-accuracy radial quadrature. The soft core keeps `dU/dlambda` finite at lambda=0; do not replace it with direct linear LJ coupling without addressing endpoint singularities.
- `applications/mc_alchemy.py` runs Metropolis Monte Carlo on an actual two-particle `System` with a harmonic bond and a soft-core LJ pair potential. It preserves the center of mass, samples relative configurations, and feeds FEP/TI/MBAR; the reference is independent radial quadrature.
- `applications/solvated_alchemy.py` holds a neutral ligand fixed in an Ewald water box and decouples ligand-water interactions with soft-core LJ potentials applied through `System.indexed_pair_potentials`. Water-water energy is common across lambda and cancels; FEP/TI/MBAR cross-validate the ligand decoupling. The water model's exclusions include all ligand-water pairs so the general nonbonded model does not double count them.
- `applications/flexible_docking.py` adds `SideChain` groups that are translated as rigid blocks. The score is `interaction + 0.5*k*|displacement|^2`; this is a simple strain calibration, not a trained affinity model. `dock_flexible_receptor()` samples side-chain displacements and ligand rigid-body poses.
- `System.indexed_pair_potentials` stores `(i, j, potential)` entries for pair-specific interactions. Preserve it in every `System` rebuild.
- `applications/qm_mm.py` solves a 1D finite-difference Schrödinger equation with `ParticleInBox`, optional harmonic potentials, and fixed `MMPointCharge` electrostatics. Mass is in Da and energies in kJ/mol; `ELECTRON_MASS_DA` is exported for electron tests/examples. The Hamiltonian is dense-diagonalized, so keep grids modest.
- `applications/water.py` builds a three-site water box (`SPCE_WATER` defaults). It relies on `Particle.atom_type`/`charge`, harmonic O-H bonds and H-O-H angles, `NonbondedModel`, reaction-field electrostatics by default or Ewald when passed an `EwaldElectrostatics`, and automatic exclusions. It is rigid by default: each molecule gets three SHAKE/RATTLE distance constraints (O-H, O-H, H-H). Pass `rigid=False` for flexible water. Use a conservative timestep because the bonded terms are stiff.
- `composition/constraints.py` implements SHAKE positions and RATTLE velocities. `System.constraints` is applied automatically by `velocity_verlet_step_system` (SHAKE after the position update, RATTLE after the velocity update). Minimization and barostat scaling also re-apply SHAKE. Preserve `constraints` in every `System` rebuild.
- `applications/docking.py` scores only receptor-ligand pairs with the receptor's `nonbonded_model`; ligand internal energy is constant for a rigid ligand. Sampling is uniform in a translation ball and over SO(3); results are sorted by score. It has no desolvation, entropy, flexible receptor, or local refinement, so scores are interaction energies rather than binding free energies.
- `models/ewald.py` implements periodic Ewald electrostatics. `NonbondedModel` evaluates the real-space `erfc` part pair by pair; `System` adds the reciprocal, self, and net-charge terms once. Excluded 1-2/1-3 pairs are skipped in real space and corrected by subtracting their reciprocal-space contribution; self-image terms remain. Ewald requires a `Box`. `EwaldElectrostatics(alpha, reciprocal_cutoff)` controls accuracy; the model cutoff is the real-space cutoff.
- `protocols/pressure.py` computes kinetic pressure and a finite-difference virial pressure by isotropically scaling positions and the box, then applies the Berendsen barostat. Pressure is exposed in kJ/(mol nm^3); `pressure_bar()` converts to bar (`1 kJ/(mol nm^3) ≈ 16.6054 bar`). The barostat scales positions and box but not velocities. It requires a `Box`, and systems near `cutoff = L/2` use a one-sided virial difference when shrinking would violate the cutoff.
- Preserve `pair_potentials`, `topology`, `box`, `nonbonded_model`, and `constraints` whenever rebuilding a `System` in integrators, minimization, or thermostats. Also preserve each particle's `atom_type` and `charge`.
- Keep layers separate: nonbonded/bonded models, `System`, Velocity Verlet, minimization, thermostats, and PBC each have their own module.

## Testing and examples

- Every physical equation needs an implementation, a unit test, and a numerical experiment.
- Test **physics relationships**, not just that functions run: equilibrium values, symmetry, Newton's third law, energy conservation, infinite-dilution limits, and analytical force vs. central finite difference (`h=1e-6`, `atol=1e-6`). The validation levels are documented in `docs/validation.md`.
- Fundamental constants need independent golden values, closed-form checks, sign invariants, and order-of-magnitude bounds in `tests/test_physical_benchmarks.py`; do not rely only on self-consistency tests.
- Use `np.testing.assert_allclose` for arrays and `pytest.approx` for scalars. Preserve documented exact examples (e.g. kinetic energy `0.06`, total energy `0.56`).
- Tests are fast (<1 s); run the full suite before finishing. Add a focused test file for each new module and a runnable example under `examples/`.

## Docs caveat

- Code blocks in `docs/` are AI-exported and often have newlines/indentation collapsed. Treat them as intent, not copy-paste source; write clean, PEP 8 code matching the described behavior.
- Tutorial math must render on GitHub and in VS Code. Keep display math on one `$$ ... $$` line, use `\ ` rather than `\,` or `~` for spaces inside math, avoid math inside italics, and use `\lvert ... \rvert` instead of literal `|` in table math.
- `docs/all.md` is internal background material and is not part of the user-facing documentation. Do not link it from READMEs or tutorials; use `docs/tutorials/capstone_imatinib.md` for the user-facing imatinib capstone.
- `docs/lib.md` is authoritative for long-term architecture, units, and roadmap; `docs/energies.md` is authoritative for the first milestone's equations and expected values; `docs/architecture.md` describes the architecture as implemented. The tutorial entry point is `docs/tutorials/README.md`; the CADD-focused overview is `docs/tutorials/introduction.md`; the curriculum and prerequisites are in `docs/tutorials/tracks.md`. Track A is `docs/tutorials/tutorial.md` plus `docs/tutorials/physics/`; Track B starts at `docs/tutorials/chemistry/01_chemical_structure.md` and continues through `02_molecular_representation.md`, `03_molecular_properties.md`, and `04_protonation_and_tautomerism.md`; Track C starts at `docs/tutorials/drug_discovery/README.md`. The tracks meet in `docs/tutorials/capstone_imatinib.md`.
