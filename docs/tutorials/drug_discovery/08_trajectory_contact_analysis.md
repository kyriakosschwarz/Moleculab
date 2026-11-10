# Track C, module 08 — Trajectory contact analysis

> **Where this fits.** [Module 04](04_interaction_fingerprints.md) measured
> interaction persistence across frames. This module works one level lower:
> geometric atom and residue contacts from coordinate frames.
>
> **What you should already know.** Track C modules 02–04 and Track A's
> trajectory and periodic-boundary material.
>
> **What you will be able to do afterwards.** Build contact matrices, aggregate
> them by residue, compute sampled-frame occupancy, and distinguish occupancy
> from contact lifetime.

## 1. Learning objectives

**Core**

- Define a geometric protein–ligand contact.
- Build atom-level and residue-level contact matrices.
- Compute sampled-frame contact occupancy.
- Compute contact run lengths and sampled lifetimes.
- Handle missing frames separately from frames with no contacts.

**Deeper insight**

- A contact is a geometric observation, not a chemical interaction.
- Occupancy and lifetime answer different questions.
- Two contacts with the same occupancy can have very different temporal
  behaviour.
- Lifetime estimates depend on frame spacing, trajectory length, and the
  definition of contact breaking and reforming.
- Persistence is not binding affinity or a binding free energy.

**Looking ahead**

- The capstone will connect persistent contacts to structural interpretation.
- A future thermodynamics module will explain why persistence alone does not
  establish a favourable free-energy contribution.

## 2. Contact definition

A geometric contact between protein atom $i$ and ligand atom $j$ at frame $t$
is

$$ c_{ij}(t)=1\quad\text{when}\quad \lVert\mathbf r_i(t)-\mathbf r_j(t)\rVert < r_c. $$

The default comparison is strict. An inclusive variant uses $\le r_c$.

Important cautions:

- a contact is not a hydrogen bond;
- a contact is not a hydrophobic interaction;
- a contact is not evidence of a favourable energy;
- the cutoff is a geometric convention, not a physical boundary.

Use [interaction geometry](03_protein_ligand_interactions.md) when the question
is chemical interaction type rather than proximity.

### 2.1 Coordinates and boundaries

Positions are in nm. The contact layer does not apply periodic minimum images.

If the coordinate frames use periodic boundaries, the caller must prepare
coordinates consistently. Wrapping a trajectory after the fact can create
apparent jumps and must be handled explicitly.

## 3. Contact matrices

For one frame, the atom-level contact matrix is

$$ C(t)\in\{0,1\}^{N_{\mathrm{protein}}\times N_{\mathrm{ligand}}}. $$

Rows are protein atoms; columns are ligand atoms. An entry is 1 when the pair
is within the cutoff.

The implementation computes the full pairwise distance matrix explicitly. This
is appropriate for the small teaching systems. A spatial index can be added
later if profiling demonstrates a need.

A contact matrix is useful even before any chemistry:

- it shows which ligand atoms touch which protein atoms;
- it supports residue aggregation;
- it can be computed for any coordinate frame;
- it is independent of a trajectory file format.

## 4. Residue-level contacts

A residue contacts a ligand atom when at least one of its protein atoms
contacts that ligand atom.

This produces a residue-by-ligand contact matrix:

$$ R(t)\in\{0,1\}^{N_{\mathrm{residues}}\times N_{\mathrm{ligand}}}. $$

The residue identity must come from the structure adapter. Use chain and
insertion-code information where needed; a residue array index is not a stable
identity.

Residue-level contacts are easier to interpret than atom-level contacts, but
they hide which protein atom made the contact. Keep the atom-level matrix if
that detail matters.

## 5. Sampled-frame occupancy

For a contact matrix over valid frames,

$$ O_{ij}=\frac{1}{N_{\mathrm{valid}}}\sum_{t\in\mathrm{valid}}c_{ij}(t). $$

Values lie between 0 and 1. An occupancy of 0.5 means the contact was present
in half of the valid analysed frames.

A missing or invalid frame is not a frame with no contacts. The denominator
uses valid frames only, and the result reports total, valid, and missing frame
counts.

This is a descriptive statistic for the sampled frames. It is not automatically
an equilibrium probability.

## 6. Contact lifetimes and transitions

Occupancy alone does not show temporal structure.

A contact present in 50% of frames could be:

- one continuous contact lasting half the trajectory;
- many short contacts appearing intermittently.

For a binary presence sequence, a run is a maximal sequence of consecutive
frames in which the contact is present. The implementation reports:

- run lengths in frames;
- number of separate runs;
- observed fraction;
- mean and maximum run length;
- optional sampled durations when frame spacing is supplied.

### 6.1 Sampled duration

If the frame spacing is $\Delta t$, the implementation uses

$$ \text{sampled lifetime}=K\Delta t $$

for a run of $K$ consecutive frames.

This is a sampled-duration convention. If frames are instantaneous snapshots,
it overestimates the span between the first and last snapshots by one interval.
State the convention explicitly.

A one-frame contact can have a finite sampled lifetime under this convention.
That is useful for consistency, but it should not be confused with a physical
residence time.

### 6.2 The same occupancy can hide different lifetimes

The example constructs two contacts with occupancy 0.5:

```text
continuous contact: runs (6,), mean 0.600 ns
intermittent contact: runs (1, 1, 1, 1, 1, 1), mean 0.100 ns
```

Both are present in half the frames. Their temporal behaviour is completely
different.

### 6.3 Contact transitions

Each run is one appearance of the contact. The number of runs measures how
often the contact appears, but not whether the underlying motion is a genuine
break-and-reform event or an artefact of the frame spacing.

A contact can flicker across the cutoff because of thermal motion, coordinate
noise, or a cutoff that sits near a typical distance. Do not overinterpret short
runs without checking the trajectory and the cutoff sensitivity.

## 7. Missing frames and sampling

Record separately:

- valid frames in which the contact is absent;
- frames that could not be analysed;
- frames where the ligand or protein atoms were missing or incompatible.

If the trajectory is short relative to the contact lifetime, the observed
occupancy and lifetime may not represent the accessible behaviour.

Adjacent MD frames are temporally correlated. A long run of identical contact
states is not a large number of independent observations.

Use block averaging, autocorrelation analysis, or effective-sample-size
estimates when uncertainty matters. This module reports descriptive run
statistics, not calibrated uncertainties.

### 7.1 Occupancy and lifetime have different missing-data requirements

Consider the sequence:

```text
contact:  1  1  ?  1  0
```

The question mark is a missing or invalid frame.

Occupancy uses the four valid frames and gives

$$ O = \frac{3}{4} = 0.75. $$

Lifetime is more delicate. The implementation does not infer continuity across
the missing frame. It reports two observed runs, `(2, 1)`, and records one
missing frame. It does not claim that the first three observations form one
continuous run.

This is a conservative policy. Bridging a missing frame requires an explicit
model and should not happen silently.

### 7.2 Censoring

A contact that is already present in the first frame, or still present in the
last frame, may be left- or right-censored. Its true lifetime could extend
beyond the observed trajectory.

The reported run length is therefore an observed duration, not necessarily the
full physical lifetime.

## 8. Interpretation

Contact persistence can support statements such as:

- this contact is present in most sampled frames;
- this contact appears transiently;
- this contact switches between two ligand atoms;
- this residue remains in contact while the ligand changes orientation.

It does not establish:

- that the contact is favourable;
- that the contact contributes a specific binding free energy;
- that the contact is the cause of binding;
- that the sampled trajectory is equilibrated;
- that the contact would persist at another temperature or protonation state.

Persistence is evidence about the sampled trajectory. Thermodynamic
interpretation requires a free-energy calculation with a defined model.

## 9. Implementation and validation

The implementation is in
`src/moleculab/analysis/trajectory_contacts.py`. It is parser-independent and
does not load trajectory files.

### 9.1 API

| Object or function | Purpose |
| --- | --- |
| `calculate_contact_matrix(...)` | atom-level boolean contact matrix for one frame |
| `ResidueContactMatrix` | residue IDs and residue-level contact matrix |
| `aggregate_contacts_by_residue(...)` | aggregate atom contacts into residues |
| `ContactOccupancyResult` | occupancy, counts, and frame accounting |
| `calculate_contact_occupancy(...)` | occupancy over valid frames |
| `ContactLifetimeResult` | run lengths, observed fraction, missing count, sampled durations |
| `calculate_contact_lifetimes(...)` | run statistics for one contact |

### 9.2 Example

```bash
uv run python examples/drug_discovery/trajectory_contacts.py
```

Observed output:

```text
Atom contact matrix (rows = protein atoms, columns = ligand atoms)
[[1 0]
 [1 0]
 [1 1]
 [0 0]]

Residue contact matrix
  ASP A 123: [1, 0]
  PHE A 87: [1, 1]
  LYS A 45: [0, 0]

Sampled-frame occupancy
  continuous contact: 0.500
  intermittent contact: 0.500
  valid frames: 12/12

Contact lifetimes (frame spacing = 0.1 ns)
  continuous: runs (6,), mean 0.600 ns
  intermittent: runs (1, 1, 1, 1, 1, 1), mean 0.100 ns
```

### 9.3 Tests

```bash
uv run pytest tests/analysis/test_trajectory_contacts.py -v
```

The tests cover:

- contact matrix geometry;
- strict and inclusive cutoffs;
- translation and rotation invariance;
- residue aggregation;
- occupancy over valid frames;
- missing-frame handling;
- contact lifetimes and run lengths;
- same occupancy with different lifetimes;
- missing frames breaking observed continuity;
- all-present and all-absent contacts;
- cutoff and frame-spacing validation.

### 9.4 Acceptance criteria

A trajectory-contact module is correct when:

- contacts are computed in nm with a documented cutoff;
- atom-level and residue-level matrices are distinct;
- residue aggregation uses a stable supplied residue identity;
- missing frames are excluded from the denominator and reported;
- missing frames break observed lifetime continuity rather than being silently
  bridged;
- occupancy is 0 for never-present contacts and 1 for always-present contacts;
- lifetime run lengths are computed from consecutive frames;
- the frame-spacing convention is explicit;
- equal occupancy does not imply equal lifetime;
- contact persistence is not reported as binding affinity.

## 10. Exercises

### Calculation

1. For a 4-by-2 contact matrix, compute the residue-level matrix for three
   residue labels.
2. A contact is present in 6 of 10 valid frames. What is its occupancy?
3. Two frames are missing. Recompute the occupancy over valid frames.
4. A contact has run lengths 2, 5, and 1 with 0.2 ns spacing. What is the mean
   sampled lifetime?
5. Construct two 12-frame patterns with occupancy 0.5 but different mean
   lifetimes.

### Reasoning

1. Why is a contact not the same as a chemical interaction?
2. Why does occupancy not distinguish a long contact from many short ones?
3. Why does the lifetime estimate depend on frame spacing?
4. Why are adjacent MD frames not independent observations?
5. Why is a persistent contact not automatically favourable?
6. Why should missing frames be separated from frames with no contacts?
7. Why is a wrapped trajectory a potential source of false contact transitions?

### Implementation

1. Run `examples/drug_discovery/trajectory_contacts.py` and reproduce the output.
2. Change the cutoff and observe how contacts and lifetimes change.
3. Add a missing frame and confirm the denominator changes.
4. Add a residue with two contacting atoms and verify aggregation.
5. Compare a continuous contact and an intermittent contact with equal
   occupancy.

## 11. Summary

- A geometric contact is a proximity observation, not a chemical interaction.
- Atom-level contact matrices can be aggregated into residue-level matrices.
- Sampled-frame occupancy is the fraction of valid frames in which a contact is
  present.
- Contact lifetimes are run lengths; sampled durations depend on frame spacing.
- Equal occupancy can hide very different temporal behaviour.
- Missing frames are excluded from the denominator and reported separately.
- Contact persistence is not a binding free energy or an affinity measure.
- Trajectory file loading and chemical interpretation remain outside this
  parser-independent layer.

**Previous:** [07 — Docking validation and benchmarking](07_docking_validation.md).

**Next:** [09 — Lead optimisation and structure–activity relationships](09_lead_optimisation_sar.md).
