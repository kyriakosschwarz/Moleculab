# Validation milestone

This document defines how the current Moleculab foundation is validated. The
goal is to establish that the library behaves correctly as a connected physical
system, not merely that individual functions return finite numbers.

The tutorials teach the physics; this milestone checks that the implementation
keeps its promises. The architectural context is in
[`architecture.md`](../development/architecture.md).

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
uv run pytest tests/test_acid_base.py -q
uv run pytest tests/test_binding_site.py -q
uv run pytest tests/test_interactions.py -q
uv run pytest tests/test_fingerprints.py -q
uv run pytest tests/test_docking_poses.py -q
uv run pytest tests/test_pose_comparison.py -q
uv run pytest tests/test_docking_benchmark.py -q
uv run pytest tests/test_trajectory_contacts.py -q
uv run pytest tests/test_sar.py -q
uv run pytest tests/test_structure_manifest.py -q
uv run pytest tests/test_ligand_identity.py -q
uv run pytest tests/test_ligand_coordinates.py -q
uv run pytest tests/test_binding_site_reconciliation.py -q
uv run pytest tests/test_interaction_interpretation.py -q
uv run pytest tests/test_evidence_fingerprint.py -q
uv run pytest tests/test_chemical_preparation.py -q
uv run pytest tests/test_chemical_state_enumeration.py -q
uv run pytest tests/test_state_interactions.py -q
uv run pytest tests/test_protein_preparation.py -q
uv run pytest tests/test_electrostatics.py -q
uv run pytest tests/test_adme.py -q
```

Covered:

- unit labels and conversion factors;
- positive mass and three-component vectors;
- harmonic, Lennard-Jones, Coulomb, bond, angle, torsion, reaction-field, and
  typed nonbonded gradients;
- Lorentz-Berthelot combining and exclusions;
- Henderson-Hasselbalch limits, microstate normalisation, and Boltzmann
  populations;
- binding-site distances, cutoff selection, residue aggregation, and
  rigid-motion invariance;
- hydrogen-bond distance/angle criteria, donor direction, and negative cases;
- fingerprint schema determinism, duplicate handling, and occupancy;
- pose RMSD, explicit atom mappings, symmetry-aware mapping, and recovery
  thresholds;
- pose clustering, centroid comparison, and interaction-pattern Jaccard;
- benchmark statuses, valid-case denominators, and protocol retention;
- contact matrices, occupancy, missing frames, and lifetime run lengths;
- pIC50 conversion, ligand efficiency, and matched-pair category separation;
- dipole moments, origin dependence, and charge/position validation;
- pH-dependent LogD corrections for monoprotic acids and bases;
- limiting behavior and numerical stability at extreme pH;
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

uv run python examples/physics/first_particle.py
uv run python examples/physics/harmonic_oscillator.py
uv run python examples/physics/harmonic_validation.py
uv run python examples/physics/nonbonded_interactions.py
uv run python examples/physics/bonded_molecule.py
uv run python examples/physics/minimization.py
uv run python examples/physics/thermostat.py
uv run python examples/physics/periodic_box.py
uv run python examples/chemistry/protonation.py
uv run --with rdkit python examples/chemistry/stereochemistry.py
uv run --with rdkit python examples/chemistry/molecular_polarity.py
uv run --with rdkit python examples/chemistry/hydrogen_bonding.py
uv run --with rdkit python examples/chemistry/adme_tradeoffs.py
uv run --with rdkit python examples/drug_discovery/capstone_imatinib.py
uv run --with gemmi python examples/drug_discovery/protein_structure.py
uv run python examples/drug_discovery/binding_site.py
uv run python examples/drug_discovery/interactions.py
uv run python examples/drug_discovery/interaction_fingerprints.py
uv run python examples/drug_discovery/docking_pose_analysis.py
uv run python examples/drug_discovery/compare_binding_modes.py
uv run python examples/drug_discovery/docking_benchmark.py
uv run python examples/drug_discovery/trajectory_contacts.py
uv run python examples/drug_discovery/matched_molecular_pairs.py
uv run python examples/drug_discovery/structural_observations_thermodynamics.py

uv run python docs/tutorials/physics/generate_visuals.py
uv run python docs/generate_architecture_visuals.py
```

The checkpoint in
[`tutorials/physics/checkpoint_04.md`](../tutorials/physics/checkpoint_04.md)
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
[`tests/test_physical_benchmarks.py`](../../tests/test_physical_benchmarks.py).
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
| manifest checksums, selections, decisions, actors, and machine-readable checks are validated | `tests/test_structure_manifest.py` |
| ligand formula, connectivity, bond orders, mapping ambiguity, and unresolved discrepancies are tested | `tests/test_ligand_identity.py` |
| coordinate mapping, occupancy, alternate locations, geometry, and missing atoms are tested | `tests/test_ligand_coordinates.py` |
| pocket selection, residue namespaces, water policy, contacts, and proximity boundaries are tested | `tests/test_binding_site_reconciliation.py` |
| interaction classification, missing evidence, and the no-affinity boundary are tested | `tests/test_interaction_interpretation.py` |
| three-state evidence, assessable denominators, schema compatibility, and unknown-not-zero comparisons are tested | `tests/test_evidence_fingerprint.py` |
| prepared states, generated hydrogens, explicit assumptions, and unselected alternatives are tested | `tests/test_chemical_preparation.py` |
| protonation/tautomer enumeration, pKa provenance, and no-population assumptions are tested | `tests/test_chemical_state_enumeration.py` |
| state-aware contact interpretation preserves state differences without ranking or affinity inference | `tests/test_state_interactions.py` |
| protein-pocket preparation preserves templates, multiple states, and unassigned water/hetero roles | `tests/test_protein_preparation.py` |
| acid-base fractions match closed-form Henderson-Hasselbalch values | `tests/test_acid_base.py` |
| protonation microstates sum to one and match site marginals | `tests/test_acid_base.py` |
| Boltzmann ratios match known free-energy differences | `tests/test_acid_base.py` |
| binding-site distances match hand calculations | `tests/test_binding_site.py` |
| pocket selection is translation- and rotation-invariant | `tests/test_binding_site.py` |
| residue aggregation preserves residue identity without duplicates | `tests/test_binding_site.py` |
| hydrogen-bond distance and angle boundaries are explicit | `tests/test_interactions.py` |
| hydrogen-bond direction and invalid chemistry are tested | `tests/test_interactions.py` |
| interaction detection is rigid-motion and ordering invariant | `tests/test_interactions.py` |
| fingerprint features are deterministic and duplicate-safe | `tests/test_fingerprints.py` |
| occupancy distinguishes missing frames from absent features | `tests/test_fingerprints.py` |
| pose RMSD is zero for identical mapped coordinates and rigid-motion invariant | `tests/test_docking_poses.py` |
| symmetry-aware pose comparison uses explicit mappings | `tests/test_docking_poses.py` |
| pose clustering is deterministic for a fixed distance matrix | `tests/test_pose_comparison.py` |
| geometry and interaction similarity are reported separately | `tests/test_pose_comparison.py` |
| benchmark success rate excludes invalid cases from the denominator | `tests/test_docking_benchmark.py` |
| benchmark statuses agree with pose-recovery results | `tests/test_docking_benchmark.py` |
| contact occupancy excludes missing frames and reports them | `tests/test_trajectory_contacts.py` |
| equal occupancy can yield different contact lifetimes | `tests/test_trajectory_contacts.py` |
| pIC50 and fold changes agree with the concentration ratio | `tests/test_sar.py` |
| matched-pair activity, descriptors, and docking scores stay separate | `tests/test_sar.py` |
| neutral dipoles are origin-independent; charged dipoles require an origin | `tests/test_electrostatics.py` |
| dipole vectors rotate with coordinates and preserve magnitude | `tests/test_electrostatics.py` |
| dipole magnitudes convert between e nm and Debye | `tests/test_electrostatics.py` |
| LogD differs from LogP by the ionisation correction | `tests/test_adme.py` |
| acid and base LogD corrections have the correct pH direction | `tests/test_adme.py` |
| LogD limiting behavior and extreme-pH stability are tested | `tests/test_adme.py` |

## Current status

```text
$ uv run pytest -q
536 passed, 1 skipped
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

All 947 math expressions in the physics deep dives are checked against both
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
