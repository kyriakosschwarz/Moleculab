# Validation milestone

This document defines how the current Moleculab foundation is validated. The
goal is to establish that the library behaves correctly as a connected physical
system, not merely that individual functions return finite numbers.

The tutorials teach the physics; this milestone checks that the implementation
keeps its promises. The architectural context is in
[`architecture.md`](architecture.md).

## Guiding distinction

Not every check is the same kind of claim.

| Type | Meaning | Example |
| --- | --- | --- |
| identity | should hold to numerical precision | $\mathbf F=-\nabla U$ |
| model property | follows from the interaction model | equal-and-opposite pair forces |
| numerical expectation | holds for a method in a regime | second-order Verlet convergence |
| algorithm expectation | describes a procedure | minimisation lowers energy |
| protocol expectation | describes a sampling method | Berendsen approaches a target |

Exact identities are tested with tight tolerances. Numerical expectations are
tested with tolerances appropriate to a finite timestep and reported
quantitatively.

## Level 1 — Individual components

**Question:** does each model produce a consistent energy and force, and does
input validation reject invalid states?

```bash
uv run pytest tests/test_units.py tests/test_particles.py -q
uv run pytest tests/test_potentials.py tests/test_pair_potentials.py -q
uv run pytest tests/test_force_gradients.py -q
uv run pytest tests/test_nonbonded.py tests/test_bonded.py tests/test_torsion.py -q
```

Covered:

- unit labels and conversion factors;
- positive mass and three-component vectors;
- harmonic, Lennard-Jones, Coulomb, bond, angle, torsion, reaction-field, and
  typed nonbonded gradients;
- Lorentz-Berthelot combining and exclusions;
- parameter validation and shape checks.

The force-gradient audit uses only the energy function in its numerical
derivative. It never calls the analytical force internally.

## Level 2 — Combined systems

**Question:** do many-particle sums and symmetries behave correctly?

```bash
uv run pytest tests/test_system.py tests/test_system_integrator.py -q
uv run pytest tests/test_symmetry.py tests/test_box.py -q
uv run pytest tests/test_neighbour.py -q
```

Covered:

- each unique pair is counted once;
- two-particle forces are equal and opposite;
- total internal force is zero;
- translation leaves internal energy and forces unchanged;
- rotation leaves internal energy unchanged and rotates forces;
- wrapping preserves periodic pair interactions;
- lattice translations preserve periodic energy;
- a quarter-turn of a cubic box preserves energy and rotates forces;
- the cutoff does not exceed half the shortest box edge.

## Level 3 — Numerical behaviour

**Question:** does the dynamics behave as the numerical analysis predicts?

```bash
uv run pytest tests/test_integrators.py tests/test_harmonic_validation.py -q
uv run pytest tests/test_minimization.py tests/test_thermostats.py -q
```

Covered:

- Velocity Verlet evaluates the force at the updated position;
- the numerical trajectory matches the analytic harmonic solution;
- position, velocity, and energy errors are second order in $\Delta t$;
- total energy stays bounded for a conservative oscillator;
- the measured period matches the analytic period;
- minimisation lowers energy and converges on force norm;
- heating raises kinetic energy and cooling lowers it;
- state and chemical identity survive every rebuild path.

## Level 4 — Reproducibility

**Question:** can a reader reproduce the documented results from a clean
checkout?

```bash
uv sync
uv run pytest -q

uv run python examples/first_particle.py
uv run python examples/harmonic_oscillator.py
uv run python examples/harmonic_validation.py
uv run python examples/nonbonded_interactions.py
uv run python examples/bonded_molecule.py
uv run python examples/minimization.py
uv run python examples/thermostat.py
uv run python examples/periodic_box.py

uv run python docs/tutorials/individual/generate_visuals.py
uv run python docs/generate_architecture_visuals.py
```

The checkpoint in
[`tutorials/individual/checkpoint_04.md`](tutorials/individual/checkpoint_04.md)
reproduces the analytic oscillator table independently.

## Independent physical benchmarks

Self-consistency tests are not enough. A test that recomputes a constant from
the same formula, or compares a potential to a constant imported from the same
module, can pass even when the constant is wrong. The suite therefore includes
four independent layers:

| Layer | Question | Example |
| --- | --- | --- |
| golden values | does the constant match an independent literature value? | $k_e\approx138.935456$ |
| closed-form smoke tests | does a hand-calculable geometry give the expected number? | two unit charges at 1 nm |
| sign and direction invariants | does the physics have the right sign? | like charges repel |
| magnitude bounds | is the value in the right order of magnitude? | $100<k_e<200$ |

These live in
[`tests/test_physical_benchmarks.py`](../tests/test_physical_benchmarks.py).
They are deliberately redundant with the derivation tests: redundancy is the
point when a single operator typo could corrupt every downstream result.

## Where the invariants live

| Invariant | Test |
| --- | --- |
| dimensional consistency | `tests/test_units.py` |
| golden constants and magnitude bounds | `tests/test_physical_benchmarks.py` |
| closed-form Coulomb and Lennard-Jones references | `tests/test_physical_benchmarks.py` |
| Coulomb sign and direction | `tests/test_physical_benchmarks.py`, `tests/test_pair_potentials.py` |
| non-negative kinetic energy | `tests/test_physics.py` |
| force is negative energy gradient | `tests/test_force_gradients.py` |
| double-well minima, barrier, and gradient | `tests/test_potentials.py`, `tests/test_force_gradients.py` |
| second-order convergence | `tests/test_harmonic_validation.py` |
| bounded oscillator energy | `tests/test_integrators.py` |
| equal-and-opposite pair forces | `tests/test_pair_potentials.py`, `tests/test_symmetry.py` |
| zero total internal force | `tests/test_symmetry.py` |
| translation and rotation invariance | `tests/test_symmetry.py` |
| periodic wrapping invariance | `tests/test_symmetry.py`, `tests/test_box.py` |
| force-based minimisation convergence | `tests/test_minimization.py` |
| thermostat heating and cooling | `tests/test_thermostats.py` |
| state and identity preservation | `tests/test_state_preservation.py` |
| constraint satisfaction after state updates | `tests/test_constraints.py` |
| cell-list candidates contain every true pair | `tests/test_neighbour.py` |
| opt-in neighbour-list energies and forces match the default path | `tests/test_neighbour.py` |

## Current status

```text
$ uv run pytest -q
382 passed, 1 skipped
```

The skip is the PyTorch-unavailable test; PyTorch is installed in this
environment, so the backend path is exercised.

All core examples listed in Level 4 run successfully, and the visual generator
reproduces the tutorial assets.

## Math rendering compatibility

The tutorials must render on both GitHub and the VS Code Markdown preview.
GitHub's Markdown parser and KaTeX have different quirks, so the following
conventions are enforced:

- **Display math stays on one line**: `$$ ... $$`. Multi-line display blocks can
  be parsed as Markdown headings or lists when a line starts with `-`, `=`, or
  similar.
- **Use `\ ` for spaces inside math**, not `\,` or `~`. GitHub turns `\,` into a
  comma and drops `~`; `\ ` survives.
- **Do not put inline math inside italic spans.** GitHub renders it as literal
  dollar-delimited text. Image captions containing math are therefore not
  italicised.
- **Do not use literal `|` inside math in table cells.** It breaks the table
  parser. Use `\lvert ... \rvert`.
- **Use `^\circ` for degrees inside math**, not the Unicode degree sign.
- **Avoid uncommon macros** such as `\boxed` when a plain display equation is
  sufficient.

All 947 math expressions in the individual tutorials are checked against both
KaTeX and MathJax, and the GitHub Markdown API is used to confirm that every
expression is recognised as math.

## Known limitations

- The Berendsen thermostat is an equilibration method, not an exact canonical
  sampler.
- Lennard-Jones uses a plain cutoff by default; shifted and switched forms are
  discussed but not implemented in the simple path.
- Periodic rotation invariance is tested for a cubic box and a symmetry-compatible
  rotation. General rotation invariance is tested for non-periodic systems.
- Boxes are orthorhombic; triclinic cells are not implemented.
- Long-range electrostatics is validated through the existing reaction-field and
  Ewald tests, but a full Ewald accuracy study is outside this milestone.
