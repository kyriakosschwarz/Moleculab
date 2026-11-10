# Troubleshooting and error diagnosis

> **Use this with the main tutorial and the individual tutorials.** Most
> molecular-simulation failures are not mysterious. They fall into a small
> number of categories: wrong units, wrong forces, an unstable timestep,
> inconsistent boundaries, or invalid input geometry.

## 1. First diagnostic questions

Before changing parameters at random, ask:

1. Are all quantities in the same unit system?
2. Is the force the negative gradient of the energy?
3. Is the timestep small compared with the fastest motion?
4. Are particles overlapping or extremely close?
5. Are periodic boundaries and cutoffs consistent?
6. Are constraints applied after every state update?
7. Is the system minimised before dynamics?
8. Is the trajectory actually changing, or is it frozen?

The relevant tests are:

```bash
uv run pytest tests/test_units.py tests/test_force_gradients.py -q
uv run pytest tests/test_integrators.py tests/test_constraints.py -q
uv run pytest tests/test_box.py tests/test_neighbour.py -q
```

## 2. Symptom table

| Symptom | Likely causes | Quick check | Fix |
| --- | --- | --- | --- |
| coordinates become `NaN` or `inf` | overlapping particles, huge force, too-large timestep | inspect minimum pair distance and force magnitude | minimise, increase separation, reduce `dt` |
| energy grows without bound | unstable integrator, timestep too large, force sign error | halve `dt`; compare force with finite difference | reduce `dt`; fix force; minimise first |
| energy jumps suddenly | cutoff discontinuity, bad geometry, constraint failure | locate the step and pair distance | use a valid cutoff, separate atoms, check constraints |
| particle flies out of the box | no periodic box, bad minimum image, huge force | inspect displacement and force | add `Box` or fix displacement |
| temperature drifts or explodes | bad initial velocities, thermostat misuse, unstable dynamics | track kinetic temperature and energy | reduce `dt`, adjust `tau`, check degrees of freedom |
| minimisation energy increases | accepted step too large, broken line search | inspect trial energies | reduce step size, use backtracking |
| SHAKE/RATTLE does not converge | impossible geometry, overlapping atoms, too many coupled constraints | inspect constraint error and geometry | fix initial geometry; increase separation |
| forces are zero unexpectedly | only bonded terms and no nonbonded model, exclusions, cutoff | print pair distances and model cutoff | add the intended interaction |
| force sign is wrong | displacement convention, missing minus sign | central-difference test | implement `F = -grad U` |
| RDF never approaches 1 | bad normalization, box too small, cutoff beyond `L/2` | inspect `r_max`, density, pair counts | use a valid `r_max` and box |
| MSD saturates | wrapped coordinates at long lag times | compare with unwrapped coordinates | unwrap the trajectory |
| docking score is implausibly low | atomic clash, missing exclusions, wrong charges | inspect closest receptor-ligand distance | minimise, check parameters |
| Ewald energy does not converge | reciprocal cutoff too small, non-neutral system | increase `reciprocal_cutoff`; check net charge | increase accuracy; neutralise |
| barostat box collapses or explodes | target pressure mismatch, compressibility too large, `dt` too large | track volume and pressure | reduce `dt`, adjust coupling, check pressure |
| angle or torsion force is huge | near-linear angle or collinear torsion | inspect angle/dihedral geometry | avoid singular geometry; adjust structure |
| bonded molecule jumps across boundary | minimum image applied to bonded terms | inspect bond vector | use raw coordinates for bonded terms |
| free-energy estimate has huge error | poor overlap, too few samples/windows, endpoint singularity | inspect window overlap and variance | add windows/samples; use soft core |

## 3. Recipes

### 3.1 Timestep too large

**Symptom.** The trajectory looks plausible for a few steps, then energy grows,
the particle oscillates wildly, or coordinates explode.

**Diagnosis.** Halve `dt`. If the energy drift falls by roughly four, the
method is behaving as a second-order integrator and the original timestep was
too large. For a harmonic oscillator, check the period
$T=2\pi\sqrt{m/k}$ and resolve it with many steps.

**Fix.** Reduce `dt`. For atomistic systems, the fastest bond vibrations set the
limit. Constraints can remove the fastest modes and allow a larger step.

### 3.2 Overlapping particles

**Symptom.** Lennard-Jones energy is enormous, forces are huge, and the first
step produces `NaN`.

**Diagnosis.** Print the minimum pair distance. If it is much smaller than
`sigma`, the repulsive wall dominates.

**Fix.** Separate the atoms or minimise before dynamics. A soft-core potential
is appropriate for alchemical paths, not for preparing a physical structure.

### 3.3 Wrong force sign

**Symptom.** A particle accelerates away from the minimum, energy grows, or a
bond stretches instead of restoring.

**Diagnosis.** Run the central finite-difference check:

```python
numerical = -(U(x + h) - U(x - h)) / (2 * h)
```

Compare with the analytical force. The test suite does this for every core
potential.

**Fix.** Use $\mathbf F=-\nabla U$ and check the displacement convention.

### 3.4 Missing minimum image

**Symptom.** Two particles near opposite box faces have a huge interaction even
though they are close through the boundary.

**Diagnosis.** Compare raw separation with `box.minimum_image(...)`.

**Fix.** Apply the minimum image to nonbonded displacements. Do not apply it to
bonded terms.

### 3.5 Cutoff larger than half the box

**Symptom.** A pair interacts with more than one periodic image, or the
neighbour list is ambiguous.

**Diagnosis.** Compare `cutoff` with `min(box.lengths) / 2`.

**Fix.** Reduce the cutoff or enlarge the box. `System` rejects an invalid
typed-model cutoff.

### 3.6 Constraints are not satisfied

**Symptom.** O-H bonds drift, rigid water deforms, or SHAKE/RATTLE raises.

**Diagnosis.** Call `system.max_constraint_error()` or
`assert_constraints_satisfied(system)`.

**Fix.** Apply SHAKE after position updates and RATTLE after velocity updates.
The integrator, minimiser, and barostat already do this; the guard catches a
new algorithm that forgets.

### 3.7 Thermostat behaviour looks wrong

**Symptom.** Temperature overshoots, oscillates, or does not reach the target.

**Diagnosis.** Check the coupling time, the target, the degrees of freedom, and
whether centre-of-mass motion is removed.

**Fix.** Use a sensible `tau`, verify the degree-of-freedom count, and remember
that Berendsen is an equilibration method, not an exact canonical sampler.

### 3.8 Minimisation does not converge

**Symptom.** Energy decreases slowly, oscillates, or the force stays above
tolerance.

**Diagnosis.** Inspect the energy and force history. A small energy change is
not proof of convergence; check the force norm.

**Fix.** Reduce the step size, improve the starting geometry, or use a better
minimiser for a difficult landscape.

### 3.9 Ewald or charge problems

**Symptom.** Electrostatic energy depends strongly on `reciprocal_cutoff`, or a
charged system behaves strangely.

**Diagnosis.** Increase the reciprocal cutoff and check convergence. Check the
net charge; a non-neutral periodic system needs a charge correction.

**Fix.** Increase Ewald accuracy, neutralise the system, and verify exclusions.

### 3.10 Free-energy estimates have huge error

**Symptom.** FEP/TI/MBAR estimates vary wildly between runs or have large
standard errors.

**Diagnosis.** Inspect the overlap between neighbouring windows and the
variance of the energy derivative. Look for endpoint singularities.

**Fix.** Add windows, add samples, use a soft-core path, and check convergence.

## 4. A short diagnostic checklist

- [ ] Units are consistent.
- [ ] Masses are positive, shapes are `(3,)`.
- [ ] The force matches the energy gradient.
- [ ] The timestep resolves the fastest motion.
- [ ] The system was minimised before dynamics.
- [ ] Constraints are satisfied after every step.
- [ ] The cutoff is at most half the shortest box edge.
- [ ] Nonbonded pairs use the minimum image; bonded terms do not.
- [ ] The thermostat and barostat are appropriate for the question.
- [ ] The observable is averaged over enough samples.

If a result looks physically surprising, first check the model, then the units,
then the numerical parameters. Do not tune the thermostat or barostat to hide a
force error.
