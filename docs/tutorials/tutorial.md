# Moleculab tutorial: from one particle to a small simulation

**New here?** Read the [introduction](introduction.md) first. It explains why
the sequence starts with units and one particle and how it connects to force
fields, MD, docking, free energy, and QM/MM.

This tutorial is a guided path through the implemented library.  It is meant
to be read with the code open.  Each section answers three questions:

1. What is the physical idea?
2. Where is it implemented?
3. How can I check that it is correct?

Run everything from the repository root with the uv-managed environment:

```bash
uv sync
uv run pytest -q
```

For a slower, more visual path through sections 1–10, see
[`individual/README.md`](individual/README.md) and the ten deep-dive tutorials
inside `individual/`. The corresponding software layers, module dependencies,
and design decisions are in [`../architecture.md`](../architecture.md).

Companion guides:

- [Troubleshooting and error diagnosis](troubleshooting.md)
- [Notebook and visualization prompts](visualization_prompts.md)
- [Answer keys and expected magnitudes](answer_keys.md)

## 1. Units first

Moleculab uses one internal unit system:

| Quantity | Unit |
| --- | --- |
| distance | nm |
| time | ps |
| energy | kJ/mol |
| mass | Da |
| force | kJ/(mol nm) |
| temperature | K |
| charge | e |

The key convenience is that Newton's second law becomes numerically

$$
a = F / m
$$

without an extra conversion constant.  The constants live in
`src/moleculab/units.py`.

Check: `uv run pytest tests/test_units.py -q`.

Deep dive: [01 — Units and scales](individual/01_units.md).

## 2. A particle and its energy

A `Particle` stores only `(mass, position, velocity)`.  Acceleration is
calculated from force and mass when needed.

```python
from moleculab import Particle

particle = Particle(
    mass=12.0,
    position=[0.1, 0.0, 0.0],
    velocity=[0.1, 0.0, 0.0],
)
```

Kinetic energy is

$$
K = \tfrac{1}{2} m \mathbf v \cdot \mathbf v.
$$

For this particle, `K = 0.06 kJ/mol`.

Check: `uv run pytest tests/test_particles.py tests/test_physics.py -q`.

Deep dive: [02 — Particles and kinetic energy](individual/02_particle_energy.md).

## 3. A potential must provide energy and force

The harmonic potential is an external spring-like model:

$$
U(\mathbf r) = \tfrac{1}{2} k |\mathbf r - \mathbf r_0|^2,
\qquad
\mathbf F = -k(\mathbf r - \mathbf r_0).
$$

```python
from moleculab import HarmonicPotential

potential = HarmonicPotential(stiffness=100.0, center=[0.0, 0.0, 0.0])
```

At `r = (0.1, 0, 0)`, `U = 0.5 kJ/mol` and `F = (-10, 0, 0) kJ/(mol nm)`.
The force is not invented separately: it is the negative gradient of the
energy.  The test suite verifies this with a finite-difference derivative.

Run: `uv run python examples/first_particle.py`.
Check: `uv run pytest tests/test_potentials.py -q`.

Deep dive: [03 — Potential energy and force](individual/03_potential_force.md).

## 4. Time integration and energy conservation

Velocity Verlet updates position with the current acceleration and velocity
with the average of the old and new accelerations:

$$
\mathbf r_{n+1} = \mathbf r_n + \mathbf v_n \Delta t
+ \tfrac{1}{2}\mathbf a_n \Delta t^2,
$$

$$
\mathbf v_{n+1} = \mathbf v_n
+ \tfrac{1}{2}(\mathbf a_n + \mathbf a_{n+1})\Delta t.
$$

For a harmonic oscillator with `m = 12 Da` and `k = 12*pi^2`, the period is
2 ps.  The example starts at `x = 0.1 nm` and integrates 2000 steps of
`0.001 ps`.  The position returns to its starting value and the total energy
stays nearly constant.

Run: `uv run python examples/harmonic_oscillator.py`.

For an independent check, run
`uv run python examples/harmonic_validation.py`. It compares the numerical
trajectory with the analytic solution and reports second-order convergence as
the timestep is halved.

Check: `uv run pytest tests/test_integrators.py tests/test_harmonic_validation.py -q`.

Deep dive: [04 — Velocity Verlet and energy conservation](individual/04_velocity_verlet.md).

Checkpoint: [one oscillator, end to end](individual/checkpoint_04.md).

## 5. Nonbonded interactions

Two pair potentials are implemented.

### Lennard-Jones

$$
U(r) = 4\epsilon
\left[\left(\frac{\sigma}{r}\right)^{12}
- \left(\frac{\sigma}{r}\right)^6\right].
$$

- `U = 0` at `r = sigma`.
- The minimum is at `r = 2^(1/6) sigma`.
- The well depth there is `-epsilon`.
- The force is repulsive at short range and attractive around the minimum.

### Coulomb

$$
U(r) = k_e \frac{q_i q_j}{r}.
$$

In Moleculab units, `k_e = COULOMB_CONSTANT = 138.935456 kJ nm/(mol e²)`.
Like charges repel; opposite charges attract.

Run: `uv run python examples/nonbonded_interactions.py`.
Check: `uv run pytest tests/test_pair_potentials.py -q`.

Deep dive: [05 — Nonbonded interactions](individual/05_nonbonded.md).

## 6. Many particles: `System`

A `System` sums nonbonded interactions over every unique pair:

$$
U_{\text{nonbonded}} = \sum_{i<j} U_{\text{pair}}(r_{ij}).
$$

For each pair, the force on particle `j` is the negative of the force on
particle `i`, so total momentum is conserved when there are no external
forces.

```python
from moleculab import LennardJonesPotential, Particle, System

system = System(
    particles=[
        Particle(39.948, [0.0, 0.0, 0.0]),
        Particle(39.948, [0.38, 0.0, 0.0]),
    ],
    pair_potentials=[LennardJonesPotential(epsilon=0.997, sigma=0.34)],
)
print(system.potential_energy())
print(system.forces())
```

Check: `uv run pytest tests/test_system.py tests/test_system_integrator.py -q`.

Deep dive: [06 — Many particles and the System](individual/06_system.md).

## 7. Bonds, angles, and torsions

Nonbonded potentials act on all pairs.  Bonded potentials act only on atoms
specified by the molecular graph:

```python
import numpy as np

from moleculab import (
    Angle,
    Bond,
    HarmonicAnglePotential,
    HarmonicBondPotential,
    Topology,
)

topology = Topology(
    bonds=[Bond(0, 1, HarmonicBondPotential(stiffness=3000.0, length0=0.10))],
    angles=[
        Angle(
            1, 0, 2,
            HarmonicAnglePotential(stiffness=500.0, angle0=np.deg2rad(104.5)),
        )
    ],
)
```

Important conventions:

- Bond displacement is `r_j - r_i`.
- Angle vectors are `r_ji = r_i - r_j` and `r_jk = r_k - r_j`.
- Torsion bond vectors are `r_ij = r_j - r_i`, `r_jk = r_k - r_j`,
  `r_kl = r_l - r_k`.
- `angle0` is in radians.

Run: `uv run python examples/bonded_molecule.py` and
`uv run python examples/torsion_scan.py`.
Check: `uv run pytest tests/test_bonded.py tests/test_topology.py tests/test_torsion.py -q`.

Deep dive: [07 — Bonds, angles, and torsions](individual/07_bonded.md).

## 8. Energy minimization

At a minimum, every net force is zero.  Steepest descent moves atoms along the
force:

$$
\mathbf r_{\text{new}} = \mathbf r + \text{step} \cdot \mathbf F.
$$

`step_size` has units of nm²/(kJ/mol), and backtracking reduces it whenever a
trial move raises the energy.

Run: `uv run python examples/minimization.py`.
Check: `uv run pytest tests/test_minimization.py -q`.

Deep dive: [08 — Energy minimization](individual/08_minimization.md).

## 9. Temperature control

The instantaneous kinetic temperature is

$$
T = \frac{2K}{N_{\text{dof}} R},
$$

with `R = 0.008314462618 kJ/(mol K)`.  For an isolated system the three
center-of-mass translation modes are normally removed, so `N_dof = 3N - 3`.

The Berendsen thermostat weakly rescales velocities:

$$
\lambda = \sqrt{1 + \frac{\Delta t}{\tau}
\left(\frac{T_0}{T} - 1\right)}.
$$

Run: `uv run python examples/thermostat.py`.
Check: `uv run pytest tests/test_thermostats.py -q`.

Deep dive: [09 — Temperature and thermostats](individual/09_thermostats.md).

## 10. Periodic boundaries

An orthorhombic `Box` uses the minimum-image convention:

$$
\mathbf d_{\min} = \mathbf d
- \operatorname{round}(\mathbf d / \mathbf L) \odot \mathbf L.
$$

When `System.box` is set, nonbonded interactions use the minimum image.
Bonded terms use the stored coordinates directly.

Run: `uv run python examples/periodic_box.py`.
Check: `uv run pytest tests/test_box.py -q`.

Deep dive: [10 — Periodic boundaries](individual/10_periodic.md).

## 11. Trajectory analysis

A simulation produces more than a final structure.  Two standard analyses are:

### Mean-squared displacement

$$
\mathrm{MSD}(\tau) = \left\langle
|\mathbf r_i(t+\tau) - \mathbf r_i(t)|^2
\right\rangle_{i,t}.
$$

For diffusive motion in three dimensions, `MSD(tau) ~ 6 D tau`, so the slope
gives the diffusion coefficient.  When a periodic box is supplied, Moleculab
uses minimum-image displacements.

### Radial distribution function

$$
g(r_k) = \frac{N_k}
{N_{\mathrm{frames}} N_{\mathrm{pairs}}
(4\pi r_k^2 \Delta r / V)}.
$$

`g(r) = 1` for an ideal gas.  Peaks show solvation shells; the first peak of a
Lennard-Jones fluid appears near `r ~ sigma`.

```python
from moleculab import Box, mean_squared_displacement, radial_distribution_function

box = Box([1.1, 1.1, 1.1])
msd = mean_squared_displacement(trajectory, box=box, times=times)
rdf = radial_distribution_function(trajectory, box=box, r_max=0.55, n_bins=40)
```

Run: `uv run python examples/trajectory_analysis.py`.
Check: `uv run pytest tests/test_analysis.py -q`.

## 12. Explicit solvent

`build_solvent_box` creates an explicit coarse-grained solvent in a periodic
box.  Each solvent molecule is one Lennard-Jones particle:

$$
U(r) = 4\epsilon\left[\left(\frac{\sigma}{r}\right)^{12}
- \left(\frac{\sigma}{r}\right)^6\right].
$$

The builder can size the box from a target number density,

$$
L = \left(\frac{N}{\rho}\right)^{1/3},
$$

place particles on a lattice or randomly, and initialize velocities at a
target temperature.  The result is a normal `System` with a `Box`, so the same
thermostat, integrator, and analysis tools apply.

```python
from moleculab import COARSE_GRAINED_WATER, build_solvent_box

system = build_solvent_box(
    COARSE_GRAINED_WATER,
    64,
    density=30.0,
    temperature=300.0,
    seed=3,
)
```

This is intentionally a **coarse-grained** water model: it captures excluded
volume and dispersion, but not directional hydrogen bonding.  Atomistic water
would require per-atom types, charges, intramolecular exclusions, and
long-range electrostatics.

Run: `uv run python examples/solvent_box.py`.
Check: `uv run pytest tests/test_solvent.py -q`.

## 13. Free-energy perturbation

FEP estimates the free-energy difference between two states by sampling state
0 and reweighting each sample by the work needed to switch to state 1:

$$
\Delta G = -RT \ln
\left\langle \exp\left(-\frac{\Delta U}{RT}\right) \right\rangle_0,
\qquad
\Delta U = U_1 - U_0.
$$

Moleculab uses molar units, so `RT = 0.008314462618 * T` kJ/mol.  The module's
test system is a three-dimensional harmonic oscillator with a stiffness change,
whose exact result is

$$
\Delta G_{\mathrm{exact}} = \frac{3}{2} RT \ln\left(\frac{k_1}{k_0}\right).
$$

```python
from moleculab import (
    HarmonicPotential,
    fep_harmonic_oscillator,
    harmonic_free_energy_difference,
)

potential0 = HarmonicPotential(stiffness=100.0, center=[0.0, 0.0, 0.0])
potential1 = HarmonicPotential(stiffness=200.0, center=[0.0, 0.0, 0.0])

exact = harmonic_free_energy_difference(potential0, potential1, 300.0)
result = fep_harmonic_oscillator(
    potential0, potential1, 300.0, n_samples=100_000, seed=2024
)
print(exact, result.delta_g, result.standard_error)
```

FEP variance grows when the two states have poor overlap.  That is why
practical alchemical calculations use many intermediate states, bidirectional
estimators, or MBAR instead of a single large perturbation.

Run: `uv run python examples/free_energy_perturbation.py`.
Check: `uv run pytest tests/test_fep.py -q`.

## 14. Atomistic water, atom types, and exclusions

The coarse-grained solvent treats each molecule as one particle.  A three-site
water model instead has one oxygen and two hydrogens:

```python
from moleculab import SPCE_WATER, build_water_box

system = build_water_box(
    8,
    density=5.0,
    temperature=300.0,
    seed=7,
)
```

Each `Particle` now carries an `atom_type` and a `charge`:

```python
print(system.particles[0].atom_type, system.particles[0].charge)  # O -0.8476
```

The topology supplies two O-H bonds and one H-O-H angle per molecule.  The
nonbonded model uses:

- per-type Lennard-Jones parameters with Lorentz-Berthelot mixing;
- 1-2 exclusions for bonded pairs and 1-3 exclusions for the two outer atoms of
  each angle;
- reaction-field electrostatics with cutoff `R_c`:

$$
U(r) = k_e q_i q_j \left[
\frac{1}{r}
+ \frac{\epsilon_{rf}-1}{2\epsilon_{rf}+1}\frac{r^2}{R_c^3}
- \frac{3\epsilon_{rf}}{2\epsilon_{rf}+1}\frac{1}{R_c}
\right].
$$

The potential is shifted to zero at the cutoff.  This is a minimal
reaction-field treatment.  For long-range electrostatics, pass an
`EwaldElectrostatics` object to `build_water_box`; the Ewald path supports the
1-2/1-3 exclusions through a reciprocal-space correction.  By default the water
molecule is made rigid with three SHAKE/RATTLE distance constraints per
molecule: two O-H distances and the H-H distance.  The constrained Velocity
Verlet integrator applies SHAKE after the position update and RATTLE after the
velocity update.  Use a conservative timestep (the example uses `1e-4 ps`).

SHAKE corrects a stretched constraint with

$$
\Delta \mathbf r_i = -\frac{\mu}{m_i}\frac{r-d}{r}\mathbf r,
\qquad
\Delta \mathbf r_j = +\frac{\mu}{m_j}\frac{r-d}{r}\mathbf r,
$$

and RATTLE removes the radial relative velocity.  Pass `rigid=False` to build
flexible water for comparison.

Run: `uv run python examples/water_box.py` and
`uv run python examples/rigid_water.py`.
Check: `uv run pytest tests/test_nonbonded.py tests/test_water.py tests/test_constraints.py -q`.

## 15. Rigid-ligand docking

Docking searches for a ligand placement that interacts favorably with a
receptor.  Moleculab keeps the ligand rigid and scores only receptor-ligand
pairs with the existing typed nonbonded model:

$$
E_{\text{pose}} = \sum_{r \in \text{receptor}}
\sum_{l \in \text{ligand}}
\left(U_{\text{LJ}}^{rl} + U_{\text{Coulomb}}^{rl}\right).
$$

A rigid ligand pose is a center-of-mass translation plus a rotation.  The
sampler draws translations uniformly in a ball and rotations uniformly over
the rotation group SO(3):

```python
from moleculab import RigidLigand, dock_rigid_ligand

ligand = RigidLigand([...])
result = dock_rigid_ligand(
    receptor,
    ligand,
    center=[0.0, 0.0, 0.0],
    translation_radius=0.6,
    n_samples=1000,
    seed=42,
)
print(result.best.score, result.best.translation)
```

The score is a physical interaction energy, not a calibrated binding free
energy.  A real docking workflow also treats receptor flexibility, solvation,
entropy, and local pose refinement.

Run: `uv run python examples/rigid_docking.py`.
Check: `uv run pytest tests/test_docking.py -q`.

## 16. Ewald summation for periodic charges

Truncating Coulomb interactions in a periodic box is physically wrong because
every charge interacts with every periodic image.  Ewald summation splits the
lattice sum into a short-range real-space part and a long-range reciprocal-space
part:

$$
U = U_{\mathrm{real}} + U_{\mathrm{recip}} + U_{\mathrm{self}}
+ U_{\mathrm{charge}}.
$$

The real-space part uses

$$
U_{\mathrm{real}}(r) = k_e q_i q_j
\frac{\operatorname{erfc}(\alpha r)}{r},
$$

and the reciprocal part sums over wave vectors
$\mathbf k = 2\pi(n_x/L_x, n_y/L_y, n_z/L_z)$:

$$
U_{\mathrm{recip}} =
\frac{k_e}{2V}\sum_{\mathbf k \ne 0}
\frac{4\pi}{k^2}e^{-k^2/(4\alpha^2)}
\left|\sum_i q_i e^{i\mathbf k\cdot\mathbf r_i}\right|^2.
$$

Moleculab evaluates the real-space terms in the normal typed pair loop and adds
the reciprocal, self, and net-charge terms once in `System`.

```python
from moleculab import EwaldElectrostatics, NonbondedModel

ewald = EwaldElectrostatics(alpha=4.0, reciprocal_cutoff=30.0)
model = NonbondedModel(
    parameters={...},
    cutoff=2.0,
    electrostatics=ewald,
)
```

Accuracy is controlled by `alpha` and `reciprocal_cutoff`; the Ewald energy
should plateau as the reciprocal cutoff grows.  Ewald requires a periodic box.
Intramolecular 1-2/1-3 exclusions are supported: the real-space pair is
skipped, and its reciprocal-space contribution is subtracted by the exclusion
correction.  Self-image interactions remain part of the periodic energy.

For bonded water, supply the Ewald model to the water builder:

```python
from moleculab import EwaldElectrostatics, build_water_box

ewald = EwaldElectrostatics(alpha=6.0, reciprocal_cutoff=40.0)
system = build_water_box(
    2, box_length=3.0, temperature=300.0,
    electrostatics=ewald, seed=11,
)
```

Run: `uv run python examples/ewald.py` and
`uv run python examples/water_ewald.py`.
Check: `uv run pytest tests/test_ewald.py tests/test_water.py -q`.

## 17. Pressure coupling and NPT

Pressure has a kinetic and a virial part:

$$
P = \frac{2K}{3V} + P_{\mathrm{virial}}.
$$

Moleculab estimates the virial by isotropically scaling both the box and all
particle positions and measuring the potential-energy change:

$$
P_{\mathrm{virial}} = -\frac{\partial U}{\partial V}
\approx -\frac{U(V+\Delta V)-U(V-\Delta V)}{2\Delta V}.
$$

This automatically includes bonds, angles, nonbonded interactions, and Ewald
terms.  Pressure is reported in `kJ/(mol nm^3)` in the core and converted to
bar with `pressure_bar()`.

The Berendsen barostat applies

$$
\mu = \left[1-\frac{\kappa \Delta t}{\tau_P}(P_0-P)\right]^{1/3}
$$

and scales the box and positions by `mu`.  Higher target pressure compresses
the system; lower target pressure expands it.

```python
from moleculab import berendsen_barostat, pressure_bar

system = berendsen_barostat(
    system,
    target_pressure=30.0,   # bar
    coupling_time=0.1,      # ps
    dt=0.001,               # ps
    compressibility=1e-2,   # bar^-1 (exaggerated for a short demo)
)
print(pressure_bar(system))
```

Real water has `kappa ~ 4.5e-5 bar^-1` and needs much longer coupling times.
The barostat scales positions and the box but leaves velocities unchanged.

Run: `uv run python examples/barostat.py`.
Check: `uv run pytest tests/test_pressure.py -q`.

## 18. Thermodynamic integration

Thermodynamic integration computes a free-energy difference by integrating the
ensemble average of the Hamiltonian derivative along an alchemical path:

$$
\Delta G = \int_0^1
\left\langle \frac{\partial U_\lambda}{\partial \lambda}
\right\rangle_\lambda d\lambda.
$$

For the linear path

$$
U_\lambda = (1-\lambda)U_0 + \lambda U_1,
$$

the derivative is simply $U_1-U_0$.  The module's test system is a harmonic
oscillator with interpolated stiffness and center.  At every $\lambda$, the
potential is quadratic, so samples are drawn exactly from its Gaussian
Boltzmann distribution.

```python
from moleculab import (
    HarmonicPotential,
    thermodynamic_integration_harmonic,
)

potential0 = HarmonicPotential(stiffness=100.0, center=[0.0, 0.0, 0.0])
potential1 = HarmonicPotential(stiffness=200.0, center=[0.0, 0.0, 0.0])

result = thermodynamic_integration_harmonic(
    potential0, potential1, 300.0,
    n_windows=21, n_samples=20_000, seed=2024,
)
print(result.delta_g, result.standard_error)
```

The exact result is `1.5 * R * T * ln(k1/k0)`.  Simpson's rule converges
quickly for this smooth path; molecular alchemical transformations may need
many more windows and better overlap between neighboring states.

Run: `uv run python examples/thermodynamic_integration.py`.
Check: `uv run pytest tests/test_ti.py -q`.

## 19. MBAR over alchemical windows

MBAR combines samples from several states instead of analyzing one pair at a
time.  For reduced potentials $u_i(x)=U_i(x)/(RT)$ and $N_k$ samples from
state k, the self-consistent equations are

$$
e^{-f_i} = \sum_{n=1}^{N}
\frac{e^{-u_i(x_n)}}{\sum_k N_k e^{f_k-u_k(x_n)}}.
$$

Moleculab fixes $f_0=0$, so relative free energies are $\Delta G_i=RT f_i$.

```python
from moleculab import HarmonicPotential, mbar_harmonic

potentials = [
    HarmonicPotential(stiffness=k, center=[0.0, 0.0, 0.0])
    for k in [100.0, 150.0, 200.0]
]
result = mbar_harmonic(potentials, 300.0, n_samples=20_000, seed=2024)
print(result.delta_g)
```

The exact harmonic result for state i is `1.5 * R * T * ln(k_i/k_0)`.
The current implementation reports free energies and convergence, but not MBAR
uncertainties or the overlap matrix.  Molecular applications need many more
windows and careful overlap diagnostics.

Run: `uv run python examples/mbar.py`.
Check: `uv run pytest tests/test_mbar.py -q`.

## 20. Molecular FEP/TI/MBAR with a soft-core Lennard-Jones pair

The harmonic examples are analytically clear, but molecular interactions have
repulsive cores.  `SoftCoreLennardJones` couples a mobile particle to a fixed
particle through a harmonic trap and a soft-core LJ interaction:

$$
s(r,\lambda)=\alpha(1-\lambda)+(r/\sigma)^6,
$$

$$
U_{\mathrm{LJ}}(r,\lambda)
=4\epsilon\lambda\left(\frac{1}{s^2}-\frac{1}{s}\right).
$$

At `lambda=0` the LJ term vanishes; at `lambda=1` it is the ordinary 12-6
potential.  The soft core makes `dU/dlambda` finite at the endpoint, which
keeps TI and MBAR stable.

Samples are drawn from the exact radial distribution

$$
p_\lambda(r)\propto r^2 e^{-U_\lambda(r)/(RT)}
$$

on a high-resolution grid.  FEP, TI, and MBAR are all validated against a
high-accuracy radial-quadrature reference.

```python
from moleculab import SoftCoreLennardJones

pair = SoftCoreLennardJones()
reference = pair.reference_delta_g()
print(reference)
print(pair.fep(n_windows=11, n_samples=50_000, seed=101).delta_g)
print(pair.ti(n_windows=21, n_samples=50_000, seed=202).delta_g)
print(pair.mbar(n_windows=11, n_samples=20_000, seed=303).delta_g[-1])
```

Run: `uv run python examples/lj_alchemy.py`.
Check: `uv run pytest tests/test_lj_alchemy.py -q`.

## 21. Minimal QM/MM

The QM/MM module solves a one-dimensional time-independent Schrödinger
equation with finite differences:

$$
-\frac{\hbar^2}{2m}\frac{d^2\psi}{dx^2}+V(x)\psi=E\psi.
$$

The box has hard walls at `x=0` and `x=L`.  The Hamiltonian is tridiagonal and
is diagonalized with `numpy.linalg.eigh`.  Without an external potential, the
finite-difference levels converge to the analytic particle-in-a-box values

$$
E_n=\frac{\hbar^2\pi^2 n^2}{2mL^2}.
$$

MM point charges add

$$
V_{\mathrm{MM}}(x)=k_e q_{\mathrm{QM}}
\sum_i\frac{q_i}{|x-x_i|}.
$$

```python
from moleculab import ELECTRON_MASS_DA, MMPointCharge, ParticleInBox

box = ParticleInBox(length=1.0, mass=ELECTRON_MASS_DA, n_grid=400)
result = box.solve(n_states=3)
print(result.energies)

result = box.solve(
    mm_charges=[MMPointCharge(-0.2, 0.5), MMPointCharge(1.2, -0.5)],
    qm_charge=-1.0,
    n_states=3,
)
print(result.energies)
```

The model also accepts a harmonic potential array.  This is a teaching model:
it is one-dimensional, uses fixed MM charges, and does not solve
Hartree-Fock or DFT.

Run: `uv run python examples/qm_mm.py`.
Check: `uv run pytest tests/test_qm_mm.py -q`.

## 22. Multi-particle Monte Carlo alchemy

The trapped-LJ toy samples a radial distribution directly.  `mc_alchemy`
instead runs Metropolis Monte Carlo on a real two-particle `System`:

- two equal-mass particles connected by a harmonic bond;
- a soft-core Lennard-Jones pair interaction;
- center-of-mass-preserving trial moves;
- acceptance based on `System.potential_energy()`.

The relative-coordinate distribution is still one-dimensional, so the
reference free energy can be computed independently by radial quadrature.
FEP, TI, and MBAR use the actual MC samples.

```python
from moleculab import BondedLennardJonesMC

pair = BondedLennardJonesMC()
reference = pair.reference_delta_g()
print(pair.fep(n_windows=11, n_samples=1000, seed=101).delta_g)
print(pair.ti(n_windows=21, n_samples=1000, seed=202).delta_g)
print(pair.mbar(n_windows=11, n_samples=1000, seed=303).delta_g[-1])
```

The reference is

$$
\Delta G = -RT \ln \frac{Z_{\lambda=1}}{Z_{\lambda=0}},
\qquad
Z_\lambda = \int_0^\infty r^2
e^{-[U_{\mathrm{bond}}(r)+U_{\mathrm{LJ}}(r,\lambda)]/(RT)}dr.
$$

Run: `uv run python examples/mc_alchemy.py`.
Check: `uv run pytest tests/test_mc_alchemy.py -q`.

## 23. Solvated-ligand free energy with Ewald water

`SolvatedLigandAlchemy` combines the water model, Ewald electrostatics,
pair-indexed potentials, rigid constraints, and the molecular free-energy
estimators:

- water-water interactions use the water box's typed `NonbondedModel`, including
  Ewald and intramolecular exclusions;
- a neutral ligand is held fixed at a chosen position;
- ligand-water pairs are excluded from the general nonbonded model and handled
  by explicit soft-core Lennard-Jones potentials in
  `System.indexed_pair_potentials`;
- the water is sampled with constrained, thermostatted Velocity Verlet.

The water-water energy is identical in every alchemical state, so it cancels
from the ligand decoupling.  FEP, TI, and MBAR therefore use only the
ligand-water soft-core energies:

$$
\Delta G = \Delta G_{\mathrm{ligand-water}}.
$$

```python
from moleculab import EwaldElectrostatics, SolvatedLigandAlchemy, build_water_box

ewald = EwaldElectrostatics(alpha=4.0, reciprocal_cutoff=10.0)
water = build_water_box(
    4, box_length=2.5, temperature=300.0, seed=3,
    electrostatics=ewald,
)
alchemy = SolvatedLigandAlchemy(water, ligand_epsilon=1.0, ligand_sigma=0.3)
print(alchemy.fep(n_windows=5, n_samples=10, seed=11).delta_g)
print(alchemy.ti(n_windows=5, n_samples=10, seed=12).delta_g)
print(alchemy.mbar(n_windows=5, n_samples=10, seed=13).delta_g[-1])
```

There is no simple analytic reference for the solvated system, so the three
estimators are cross-validated against one another.

Run: `uv run python examples/solvated_alchemy.py`.
Check: `uv run pytest tests/test_solvated_alchemy.py -q`.

## 24. Flexible-receptor docking

Rigid docking assumes the receptor is fixed.  `SideChain` adds explicit
side-chain flexibility: a group of receptor atoms is translated as a rigid
block during the search.  The calibrated score adds a harmonic strain penalty:

$$
S = E_{\mathrm{interaction}}
+ \sum_s \frac{1}{2}k_s|\Delta \mathbf r_s|^2.
$$

The interaction energy is still the typed Lennard-Jones/Coulomb energy
between the ligand and receptor.  The strain term penalizes large side-chain
displacements and is a simple physical calibration, not a trained affinity
model.

```python
from moleculab import SideChain, dock_flexible_receptor

result = dock_flexible_receptor(
    receptor,
    ligand,
    [SideChain(atom_indices=(1,), max_displacement=0.3, stiffness=1.0)],
    center=[0.3, 0.0, 0.0],
    translation_radius=0.05,
    n_ligand_samples=50,
    n_receptor_samples=100,
    seed=1,
)
print(result.best.score)
print(result.best.receptor_displacements)
```

The example starts with a side-chain atom clashing with the ligand.  Rigid
docking scores a huge positive energy; flexible docking moves the side chain
and recovers a favorable pose.

Run: `uv run python examples/flexible_docking.py`.
Check: `uv run pytest tests/test_flexible_docking.py -q`.

## 25. CPU and optional acceleration

Moleculab defaults to NumPy on the CPU.  The optional PyTorch backend is
lazy-loaded and currently accelerates only the Ewald reciprocal energy/force
path.

```python
from moleculab import backend_info, set_backend

set_backend("auto")   # CUDA/ROCm or MPS if available, otherwise CPU
print(backend_info())
```

You can also set `MOLECULAB_BACKEND=auto` in the environment.  Explicit
`set_backend("torch")` runs PyTorch even on CPU for benchmarking; `auto` stays
on CPU when no supported accelerator is visible.

Install the optional ROCm framework with:

```bash
uv sync --extra rocm
```

For NVIDIA/CUDA, install a matching PyTorch wheel manually instead.  The
`rocm` extra resolves `torch` and `triton-rocm` from the PyTorch ROCm 7.2 index
declared in `pyproject.toml`; a plain `uv sync` stays CPU-only.

For AMD ROCm, verify the wheel/device:

```bash
uv pip install --index-url https://download.pytorch.org/whl/rocm7.2 \
  "torch==2.14.1+rocm7.2"
```

Use the wheel matching your ROCm version (`rocm6.4` for ROCm 6.x, for example).
On a Radeon 890M / gfx1151 system, verify:

```bash
rocminfo
hipconfig
python -c "import torch; print(torch.cuda.is_available(), torch.cuda.get_device_name(0))"
```

For APU throughput, `set_precision("float32")` selects the faster ROCm path,
and `to_torch`/`to_numpy` avoid an intermediate NumPy copy before the
non-blocking device move.  See `examples/amd_gpu_tuning.py` for zero-copy,
vectorized GPU pairwise distances, bfloat16 autocast, and optional
`torch.compile`.

Run the benchmark:

```bash
uv run python examples/backend_benchmark.py --backend cpu
uv run python examples/backend_benchmark.py --backend auto
```

`auto` selects the ROCm GPU when `torch.cuda.is_available()` is true.  Vendor
NPUs are not supported by this backend layer yet, and all other modules remain
CPU/NumPy.

Check: `uv run pytest tests/test_backend.py -q`.

## Suggested exercises

1. Change the harmonic stiffness and predict the oscillator period
   `T = 2*pi*sqrt(m/k)` before running the simulation.
2. Find the Lennard-Jones distance where the force changes from repulsive to
   attractive.  Compare with `2^(1/6) sigma`.
3. Add a second torsion term and observe how the energy landscape changes.
4. Minimize a molecule with one stretched bond and one distorted angle.
5. Run the thermostat at several target temperatures and compare the mean
   temperature after warm-up.
6. Place two particles on opposite sides of a periodic box and verify that
   the minimum-image distance is shorter than the raw distance.
7. Run `trajectory_analysis.py`, then vary the thermostat temperature and see
   how the first RDF peak and the MSD slope change.
8. Build solvent boxes at several densities, run each for a short time, and
   compare the resulting RDF first peak with the model's `sigma`.
9. Run the FEP example with a larger stiffness change and observe how the
   standard error and convergence rate change.
10. Inspect `system.topology.nonbonded_exclusions()` for a water box and verify
    that every bonded O-H and 1-3 H-H pair is excluded.
11. Increase the docking `n_samples` or translation radius and observe whether
    the best score changes; then compare the result with the ligand centered
    at the receptor origin.
12. Run the Ewald example with a larger or smaller `reciprocal_cutoff` and
    identify where the energy has converged.
13. Run the barostat example with a different target pressure or
    compressibility and compare the final box volume and pressure.
14. Build rigid and flexible water boxes, run both for the same time, and
    compare the O-H bond-length fluctuations.
15. Run the TI example with more λ windows and compare the error with the
    exact harmonic free-energy difference.
16. Run the MBAR example with fewer samples per state and observe how the
    endpoint free-energy error grows.
17. Run `water_ewald.py` with a different Ewald cutoff and compare the
    potential energy and runtime.
18. Compare FEP, TI, and MBAR in `lj_alchemy.py` using different numbers of
    windows and samples.
19. Vary the QM grid spacing in `qm_mm.py` and confirm that the box energies
    converge quadratically toward the analytic result.
20. Run `mc_alchemy.py` with different Monte Carlo step sizes and sample
    counts, and compare the acceptance/variance behavior.
21. Change the solvated ligand's size or charge and observe how the
    ligand-water decoupling free energy and estimator spread change.
22. Vary the side-chain stiffness and maximum displacement in
    `flexible_docking.py` and compare the gain over rigid docking.
23. Run `backend_benchmark.py` with CPU and the available backend and compare
    the Ewald timing and checksum.

## Where to go next

Before extending the library, work through the
[imatinib–BCR-ABL1 capstone](capstone_imatinib.md). It connects the disease
biology, kinase inhibition, force fields, water, docking, MD, free energy, and
resistance mutations to the components you have built. Use the tutorials as the foundations and the case
study as the integration exercise.

Water deserves special attention throughout. It screens electrostatics, forms
hydrogen bonds, fills pockets, and contributes entropy when displaced. The
water, solvent, Ewald, and solvated-alchemy modules exist because binding is an
aqueous process; a gas-phase interaction energy is not a binding free energy.

The remaining roadmap items are molecular FEP/TI/MBAR for complex systems,
trained/empirical binding-affinity scoring, and real QM/MM/DFT.
Each should follow the same pattern: derive the equation, implement energy and
force consistently, add a unit test, and add a numerical experiment.
