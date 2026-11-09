# Moleculab

An educational molecular physics and simulation library, built from first
principles.  The goal is not to replace GROMACS, OpenMM, or ORCA.  The goal is
to make the physics behind computational drug discovery visible: every
important equation has an implementation, a unit test, and a numerical
experiment.

Moleculab is aimed at a new computational drug discovery (CADD) user who wants
to understand what a molecular simulation engine actually does before treating
it as a black box.

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
- Tests that check physical relationships, not just that code runs.

This is still an educational engine, not a production molecular-dynamics
package.  Solvent models, free-energy methods, docking, and quantum chemistry
remain on the roadmap.

The longer roadmap lives in [`docs/lib.md`](docs/lib.md).  The first
implementation notes are in [`docs/energies.md`](docs/energies.md), and
[`docs/tutorial.md`](docs/tutorial.md) is a guided path from one particle to a
small molecular simulation.

## Install with uv

Moleculab uses a `uv`-managed virtual environment in `.venv`.

```bash
uv sync
```

This installs the package in editable mode together with the development
dependencies (`pytest`).  To run the examples, install the optional plotting
dependency:

```bash
uv sync --extra examples
```

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
- `minimization.py` relaxes a distorted molecule to its equilibrium geometry.
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
- `ewald.py` computes Ewald energies for a periodic ionic system and shows
  convergence with the reciprocal-space cutoff.
- `barostat.py` runs a short NPT-like argon simulation and reports pressure,
  temperature, and box volume before and after coupling.
- `thermodynamic_integration.py` integrates a harmonic alchemical path and
  shows convergence with the number of λ windows.
- `mbar.py` combines samples from several harmonic alchemical states and
  compares MBAR free energies with the exact result.
- `lj_alchemy.py` compares FEP, TI, and MBAR for a soft-core Lennard-Jones pair
  against a high-accuracy radial-quadrature reference.
- `mc_alchemy.py` runs actual Metropolis MC on a bonded two-particle System and
  compares FEP, TI, and MBAR with the same reference.
- `solvated_alchemy.py` decouples a neutral ligand from an Ewald water box and
  cross-validates FEP, TI, and MBAR.
- `flexible_docking.py` relaxes a clashing receptor side chain and compares the
  calibrated flexible score with rigid-receptor docking.
- `qm_mm.py` solves a finite-difference electron-in-a-box model and perturbs it
  with fixed MM point charges.

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
- [ ] Molecular FEP/TI/MBAR for complex systems
- [ ] Trained/empirical binding-affinity scoring functions
- [ ] Real QM/MM, Hartree-Fock, and DFT

The roadmap is a guide, not a race.  Each milestone should have a working
implementation, tests, and a numerical experiment before the next one begins.
