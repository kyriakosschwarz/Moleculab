# Individual tutorials 1–10

These ten tutorials expand the first ten sections of
[`../tutorial.md`](../tutorial.md) into a slower, more visual path from units
to interacting molecular systems.

![The tutorial sequence](assets/course_map.svg)

The corresponding software layers, module dependency rules, and decision
records are described in [`../../architecture.md`](../../architecture.md).

## A consistent course structure

Every tutorial follows the same pattern:

| Stage | Question |
| --- | --- |
| Physical intuition | What is the phenomenon, without equations? |
| Mathematical description | What equation describes it? |
| Units and dimensions | Is the equation physically consistent? |
| Worked example | What does it predict for a concrete case? |
| Python implementation | How is it computed? |
| Tests and verification | How do we know it is right? |
| Experiments | What can the reader explore? |
| Misconceptions | What subtle errors should be avoided? |
| Connections | What comes next? |
| Exercises and summary | Can the reader use it independently? |

![The physics to software pattern](assets/stage_flow.svg)

Material is labelled by importance:

- **Core** — required for the current tutorial.
- **Deeper insight** — useful context, safe to revisit later.
- **Looking ahead** — a preview developed in a later tutorial.

## The running examples

Tutorials 01–04 reuse one particle:

| Property | Value | Unit |
| --- | --- | --- |
| mass | 12 | Da |
| position | (0.1, 0, 0) | nm |
| velocity | (0.1, 0, 0) | nm/ps |
| speed | 0.10 | nm/ps |
| kinetic energy | 0.06 | kJ/mol |

It is an arithmetic example, not a typical room-temperature particle.
Tutorial 02 compares its speed with the thermal speed of a 12 Da particle.

From tutorial 05 onward, the recurring cast grows:

| Object | Purpose |
| --- | --- |
| a Lennard-Jones pair | nonbonded energy, force, cutoff |
| a three-particle chain | pair sums and Newton's third law |
| a distorted triatomic | bonds, angles, minimisation |
| an argon-like pair | temperature and thermostat |
| a two-particle periodic box | minimum image and cutoffs |

## The ten tutorials

| Tutorial | Physics | Code | Visuals |
| --- | --- | --- | --- |
| [01 — Units and scales](01_units.md) | dimensions, SI-to-simulation conversions, molar energy | `units.py` | scale chart, unit bridge, course map |
| [02 — Particles and kinetic energy](02_particle_energy.md) | position, velocity, speed, $K=\tfrac12mv^2$ | `Particle`, `kinetic_energy` | vector diagram, KE curves, stored/derived state, Maxwell distribution, motion GIF |
| [03 — Potential energy and force](03_potential_force.md) | energy landscapes, gradient, harmonic model | `HarmonicPotential`, `potential_energy`, `acceleration` | energy/force curves, 2D quiver, finite-difference figures |
| [04 — Velocity Verlet](04_velocity_verlet.md) | numerical integration, energy conservation, timesteps | `velocity_verlet_step`, `simulate_velocity_verlet` | algorithm flow, position/velocity, energy, phase space, oscillator GIF, drift plots |
| [05 — Nonbonded interactions](05_nonbonded.md) | Lennard-Jones, Coulomb, cutoffs, atom types | `LennardJonesPotential`, `CoulombPotential`, `NonbondedModel` | pair energy and force curves, Coulomb comparison |
| [06 — Many particles](06_system.md) | pair sums, Newton's third law, `O(N²)` cost | `System`, `nonbonded_energy`, `forces` | unique-pair and force-cancellation diagrams |
| [07 — Bonds, angles, torsions](07_bonded.md) | molecular graph, harmonic bond/angle, periodic torsion | `Topology`, bonded potentials | bond, angle, and torsion energy curves |
| [08 — Energy minimization](08_minimization.md) | steepest descent, backtracking, convergence | `minimize_energy`, `steepest_descent_step` | energy and force histories |
| [09 — Temperature and thermostats](09_thermostats.md) | degrees of freedom, kinetic temperature, Berendsen coupling | `kinetic_temperature`, `berendsen_thermostat` | temperature relaxation at several coupling times |
| [10 — Periodic boundaries](10_periodic.md) | minimum image, wrapping, cutoff rule | `Box` | periodic box and minimum-image diagram |

For larger systems, [`src/moleculab/neighbour.py`](../../../src/moleculab/neighbour.py)
provides an optional cell-list candidate screen. Use
`system.forces(cutoff=..., use_neighbour_list=True)`; the default `System` path
remains the explicit O(N²) teaching loop. See
[`examples/neighbour_list.py`](../../../examples/neighbour_list.py).

## Checkpoint after part 4

Before moving from one particle to interacting particles, complete the
[checkpoint](checkpoint_04.md). It combines units, particle state, the harmonic
potential, Velocity Verlet, energy exchange, and the analytic solution in one
short investigation.

## Physical and numerical checks

A correct implementation must preserve physical truths, not just run without
exceptions. Some of these are exact identities; others are numerical or
algorithmic expectations that hold in the appropriate regime.

| Part | Check | Type |
| --- | --- | --- |
| 01 | Every equation is dimensionally consistent; $R=N_Ak_B$. | identity |
| 02 | Mass is positive; kinetic energy is non-negative and independent of direction. | identity |
| 03 | $\mathbf F=-\nabla U$; adding a constant to $U$ changes neither force nor motion. | model identity |
| 04 | The analytic trajectory is reproduced; halving $\Delta t$ reduces the error by about four; energy stays bounded. | numerical expectation |
| 05 | Each pair force is the negative gradient of its energy; the two pair forces are equal and opposite. | model identity |
| 06 | The total force from internal pair interactions is zero; each pair is counted once. | identity |
| 07 | Translation and rotation do not change internal energy; each bonded term's forces sum to zero. | identity |
| 08 | Accepted minimisation steps never raise the energy; convergence is measured by force norm. | algorithm expectation |
| 09 | Heating raises kinetic energy, cooling lowers it, and the mean approaches the target under suitable coupling. | protocol expectation |
| 10 | Wrapping does not change interactions under the minimum image; the cutoff is at most half the shortest box edge. | convention identity |

An **identity** should hold to numerical precision. A **numerical expectation**
is a property of a method in a particular regime, tested with tolerances and
reported convergence. An **algorithm or protocol expectation** describes how a
procedure behaves, not a law of nature.

The recurring question is:

> What physical truth should still hold after we write the code?

## Models, algorithms, and protocols

The tutorials distinguish three kinds of statement:

| Concept | Category |
| --- | --- |
| Lennard-Jones potential | physical model |
| Coulomb potential | physical model |
| Harmonic bond/angle/torsion | physical model |
| Velocity Verlet | numerical algorithm |
| Steepest descent | numerical algorithm |
| Berendsen thermostat | simulation method/protocol |
| Periodic boundary conditions | simulation convention |

Changing the algorithm should not change the physical model. Changing the
protocol should not be confused with changing the energy function.

## Recommended pace

Each tutorial assumes the previous one:

1. **Units** give every later number a physical meaning.
2. **Particle state and kinetic energy** turn motion into a measurable quantity.
3. **Potential energy and force** explain where motion comes from.
4. **Velocity Verlet** combines state, energy, and force into a trajectory.
5. **Checkpoint** verifies the one-particle story end to end.
6. **Nonbonded interactions** add Lennard-Jones and Coulomb physics.
7. **`System`** sums those interactions over many particles.
8. **Bonded terms** add the molecular graph.
9. **Minimisation** finds a nearby low-energy structure.
10. **Thermostats** control the statistical temperature.
11. **Periodic boundaries** remove surfaces and prepare for bulk simulation.

Do not skip ahead to the thermostat or periodic boundaries before the
energy–force relationship and pair sums feel natural. The later molecular
mechanics and MD layers reuse exactly these ideas.

1. Read one tutorial at a time.
2. Run the short code snippet before moving on.
3. Run the listed tests and explain what each one protects against.
4. Do the self-check questions without looking at the answers.
5. Do at least one calculation, one reasoning, and one implementation exercise
   before starting the next tutorial.

## Numerical validation

Tutorial 04 is backed by an independent analytic check, not only by energy
conservation. Run:

```bash
uv run python examples/harmonic_validation.py
```

The script compares the numerical trajectory with
$x(t)=x_0\cos(\omega t)+(v_0/\omega)\sin(\omega t)$ and reports the observed
convergence order. Halving the timestep reduces the position, velocity, and
energy errors by about four, as expected for a second-order method.

The permanent checks are in
[`tests/test_harmonic_validation.py`](../../../tests/test_harmonic_validation.py),
together with the state-preservation tests in
[`tests/test_state_preservation.py`](../../../tests/test_state_preservation.py).

Every core potential is also covered by a consolidated force-gradient audit in
[`tests/test_force_gradients.py`](../../../tests/test_force_gradients.py). It
numerically differentiates the energy and compares the result with the
analytical force.

Full-system symmetries and combined-system behaviour are checked in
[`tests/test_symmetry.py`](../../../tests/test_symmetry.py). The complete
validation plan is in [`docs/validation.md`](../../validation.md).

## Visual conventions

Where a colour carries physical meaning, the figures use a consistent language:

| Colour | Typical meaning |
| --- | --- |
| blue | position, potential energy, or the current particle |
| orange | velocity or kinetic energy |
| green | total energy or a derived quantity |
| red | force, error, or danger |
| gray | reference lines, labels, and neutral annotations |

Categorical figures, such as the course map and the energy-scale chart, may use
additional colours to separate categories.

Static figures are SVG; animations are GIF. They are generated from the actual
library using:

```bash
uv run python docs/tutorials/individual/generate_visuals.py
```

The assets are in [`assets/`](assets/).
