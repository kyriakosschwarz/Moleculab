# Track C, module 04 — Interaction fingerprints and sampled-frame occupancy

> **Where this fits.** [Module 03](03_protein_ligand_interactions.md) detected
> individual geometric interactions. This module turns those records into a
> stable representation and measures how often each interaction appears across
> a collection of frames.
>
> **What you should already know.** Track C modules 01–03; the difference
> between a structure, a pose, and a trajectory.
>
> **What you will be able to do afterwards.** Build a feature schema, create
> binary and count fingerprints, compute sampled-frame occupancy over valid
> frames, and explain why it is not a binding affinity or an equilibrium
> probability.

## 1. Learning objectives

**Core**

- Distinguish an individual interaction, a fingerprint, and a sampled-frame
  occupancy.
- Define a stable feature identity for a residue-level interaction.
- Build binary and count-based fingerprints.
- Compute sampled-frame occupancy over valid frames.
- Handle missing frames separately from absent interactions.

**Deeper insight**

- Feature indexing must not depend on array position or input order.
- A shared schema is required to compare fingerprints across frames.
- Adjacent trajectory frames are not independent observations.
- Sampled-frame occupancy depends on the detection criteria and the sampling.
- A sampled-frame occupancy is a descriptive statistic, not an equilibrium
  probability.
- A geometric fingerprint is not an energy decomposition.

**Looking ahead**

- Docking can produce a collection of poses, each with a fingerprint.
- MD can produce a trajectory of frames, each with a fingerprint.
- The source of the frames belongs in metadata, not in the occupancy
  algorithm.

## 2. Three levels of representation

Keep these levels distinct:

| Level | Question | Example |
| --- | --- | --- |
| individual interaction | what was detected in one frame? | ligand N7 hydrogen-bonds to ASP A 123 at 0.29 nm |
| interaction fingerprint | which features are present? | ASP A 123 has a hydrogen-bond feature |
| sampled-frame occupancy | how persistent is a feature? | the ASP A 123 feature appears in 72% of valid frames |

The third level is more precisely called **sampled-frame occupancy**. It
describes the analysed frames, not an equilibrium probability.

They are different data products.

- An individual interaction is an observation.
- A fingerprint is a structured summary of observations in one frame.
- Sampled-frame occupancy is an aggregate over frames.

Do not collapse them into a single score.

## 3. Why interaction fingerprints are useful

A raw list of interactions is hard to compare:

- it may contain different numbers of contacts;
- the order of records may change;
- atom indices may change between structures;
- different frames may contain different subsets.

A fingerprint gives the list a fixed coordinate system. Once the feature set is
fixed, each frame becomes a binary vector or a count vector. That enables:

- frame-to-frame comparison;
- residue-level summaries;
- ligand comparison across poses;
- contact persistence over a trajectory;
- machine-readable feature tables.

The fingerprint is only as meaningful as its feature identity. If the feature
identity is unstable, the comparison is unstable.

## 4. Binary and count-based fingerprints

Suppose there are $M$ features and $N$ frames. For frame $t$, define a binary
indicator for feature $j$:

$$ f_j(t)\in\{0,1\}. $$

The frame fingerprint matrix is

$$ F_{tj}=f_j(t),\qquad F\in\{0,1\}^{N\times M}. $$

A count-based fingerprint keeps the number of observations that mapped to each
feature:

$$ c_j(t)=\text{number of observations of feature }j\text{ in frame }t. $$

The binary fingerprint is the presence indicator:

$$ f_j(t)=1\quad\text{if and only if}\quad c_j(t)>0. $$

Duplicates should not inflate the binary feature. If the same hydrogen bond is
reported twice in one frame, the bit remains 1 and the count becomes 2. Counts
are useful for diagnostics; bits are the stable presence representation.

## 5. Defining feature identity

A feature identity should be stable and explicit. A practical residue-level
identity is:

```text
(interaction_type, protein_residue, ligand_residue, ligand_atom_id)
```

For example:

```text
("hydrogen_bond", "ASP A 123", "LIG B 1", "N7")
```

The protein residue identifier must include enough information to distinguish:

- chain A from chain B;
- residue 123 from residue 123 with an insertion code;
- different residue names when the numbering scheme is reused.

The ligand atom identifier must be stable across frames. Do not use the array
position as molecular identity. Array order can change when atoms are reordered,
when a structure is rebuilt, or when a different parser is used.

### 5.1 Residue-level versus atom-pair fingerprints

Different resolutions answer different questions:

| Resolution | Feature identity | Question |
| --- | --- | --- |
| residue-level | interaction type + protein residue | which residues interact? |
| atom-pair | interaction type + protein atom + ligand atom | which specific atom pair interacts? |
| interaction-type | interaction type + residue | which categories occur? |

This module starts with residue-level features while retaining atom-level
interaction records underneath. That makes it possible to change resolution
later without losing the underlying observations.

## 6. Frame-wise detection and a shared schema

Interaction detection happens in module 03. Fingerprint construction happens
separately:

```text
frame structure
    → interaction detector
    → interaction records
    → fingerprint builder
    → frame fingerprint
```

The fingerprint builder does not load trajectories and does not run MD.

### 6.1 A shared schema

For a collection of frames, use one schema:

```python
schema = FingerprintSchema.from_interactions(all_interactions)
frame_fingerprint = build_fingerprint(frame_interactions, schema)
```

If each frame builds its own schema, a feature that is absent from that frame
will disappear from its schema. Cross-frame comparison then becomes invalid.

The schema stores an ordered tuple of `InteractionFeature` objects. Ordering is
deterministic, so feature $j$ means the same thing in every frame.

### 6.2 Detection criteria are metadata

A fingerprint should retain the criteria used by the detector. Changing a
hydrogen-bond distance or angle threshold can change sampled-frame occupancy.
The result is
only reproducible if the criteria travel with it.

The implementation stores `HydrogenBondCriteria` on each frame fingerprint when
interactions are present and carries the common criteria into the sampled-frame
occupancy result. A frame with no interactions has no criteria of its own; it
inherits the common criteria from the other valid frames.

## 7. Sampled-frame occupancy and persistence

For feature $j$, sampled-frame occupancy over valid frames is

$$ O_j=\frac{1}{N_{\mathrm{valid}}}\sum_{t\in\mathrm{valid}}f_j(t). $$

Therefore

$$ 0\le O_j\le 1. $$

A sampled-frame occupancy of 0.72 means the interaction was detected in 72% of
the valid analysed frames.

This is a descriptive statistic for the sampled frames. It is not automatically:

- an equilibrium probability;
- a residence time;
- a binding free energy;
- a measure of interaction strength.

A feature that is absent from the trajectory may still be physically
accessible. Absence is evidence about the sampled frames, not proof of
impossibility.

### 7.1 Missing frames

A missing or invalid frame is not the same as a frame in which every
interaction is absent.

- absent: the frame was analysed and the feature was not detected;
- missing: the frame could not be analysed, for example because atoms were
  missing or coordinates were invalid.

The denominator uses valid frames only. The result records the total, valid,
and missing frame counts.

## 8. Comparing ligands and binding modes

Fingerprints can compare:

- two poses of the same ligand;
- two ligands in the same binding site;
- the same ligand across replicates;
- an experimental structure and a predicted pose.

Before comparing, ensure that:

- the schemas are aligned;
- the residue identifiers are comparable;
- the ligand atom identifiers are comparable;
- the detection criteria are the same.

A shared schema is essential. Comparing fingerprints with different feature
indexing is like comparing vectors with different coordinate axes.

Useful comparisons include:

- shared versus unique features;
- sampled-frame occupancy profiles over the same feature set;
- changes in interaction type for the same residue;
- loss or gain of a contact between two ligands.

A fingerprint difference is a hypothesis about a binding-mode difference. It is
not proof of an affinity difference.

## 9. Trajectory sampling and temporal correlation

Adjacent MD frames are not independent observations. If an interaction persists
for 100 consecutive frames, that is not 100 independent observations.

Consequences:

- sampled-frame occupancy is still a valid descriptive fraction of frames;
- a naive standard error that assumes independence can be too small;
- block averaging, autocorrelation analysis, or effective-sample-size estimates
  are needed for uncertainty;
- short trajectories may not sample the relevant states at all.

The module implements sampled-frame occupancy, not uncertainty estimation.
That limitation should be stated with the result.

### 9.1 Source metadata

Sampled-frame occupancy should work for any collection of frames:

- MD snapshots;
- NMR models;
- docking poses;
- experimental structures from several conditions.

The algorithm should not assume that the frames came from MD. Record the source
and sampling interval as metadata outside the sampled-frame occupancy analysis.

## 10. Limitations and interpretation

- Feature identity is only as stable as the identifiers supplied by the
  adapter.
- Detection criteria affect sampled-frame occupancy.
- Missing frames affect the denominator.
- Temporal correlation affects uncertainty, not the raw fraction.
- A fingerprint is not an energy decomposition.
- A residue-level feature hides atom-level directionality.
- A contact may be geometrically present but energetically unimportant.
- Different protonation or tautomer states can change the feature set.

## 11. Implementation and validation

The implementation is in `src/moleculab/analysis/fingerprints.py`.

### 11.1 API

| Object or function | Purpose |
| --- | --- |
| `InteractionFeature` | stable residue-level feature identity |
| `FingerprintSchema` | ordered, shared feature set |
| `InteractionFingerprint` | binary and count fingerprint for one frame |
| `build_fingerprint(interactions, schema=None)` | build one frame fingerprint |
| `calculate_occupancy(fingerprints)` | sampled-frame occupancy over valid frames |
| `OccupancyResult` | sampled-frame occupancy, counts, frame counts, schema, criteria |

### 11.2 Synthetic example

The example uses ten frames: eight valid and two missing. It has hand-checkable
occupancies.

```bash
uv run python examples/drug_discovery/interaction_fingerprints.py
```

Output:

```text
Feature schema
  0: hydrogen_bond | ASP A 123 | LIG B 1 | N7
  1: hydrogen_bond | LYS A 45 | LIG B 1 | N1
  2: hydrogen_bond | PHE A 87 | LIG B 1 | O2

Frame fingerprints (1 = present)
  frame  1: 1 0 1
  frame  2: 1 0 0
  frame  3: 1 0 0
  frame  4: 1 0 0
  frame  5: 1 0 0
  frame  6: 1 0 0
  frame  7: 0 0 0
  frame  8: 0 0 1
  frame  9: missing
  frame 10: missing

Sampled-frame occupancy
  ASP A 123 / N7: 0.750 (6/8 valid frames)
  LYS A 45 / N1: 0.000 (0/8 valid frames)
  PHE A 87 / O2: 0.250 (2/8 valid frames)
  total frames: 10; valid: 8; missing: 2
  criteria: distance <= 0.35 nm, angle >= 120 degrees
```

The ASP feature is present in 6 of 8 valid frames, the PHE feature in 2 of 8,
and the LYS feature in none. The two missing frames are excluded from the
denominator.

### 11.3 Tests

```bash
uv run pytest tests/test_fingerprints.py -v
```

The tests cover:

- empty interaction lists;
- duplicate observations and binary presence;
- deterministic schema construction;
- atom-order independence of feature identity;
- sampled-frame occupancy 0, 1, and intermediate values;
- missing frames versus absent interactions;
- chain and insertion-code distinction;
- schema and criteria retention;
- mixed-criteria rejection;
- schema mismatch rejection;
- invalid fingerprint shapes and counts.

### 11.4 Acceptance criteria

An interaction-fingerprint module is correct when:

- empty frames produce valid empty fingerprints;
- duplicate observations do not change binary presence;
- feature indexing is deterministic and independent of atom order;
- the same feature means the same thing in every frame;
- sampled-frame occupancy is 0 for never-observed features and 1 for
  always-observed features;
- missing frames are excluded from the denominator and reported separately;
- chain identifiers and insertion codes remain distinct;
- schema and detection criteria are retained with the result;
- sampled-frame occupancy is reported as a frame statistic, not a binding
  affinity.

### 11.5 Current scope

The fingerprint layer currently consumes hydrogen-bond records, the first
detector from module 03. Adding salt bridges, hydrophobic contacts, aromatic
interactions, or cation–π features should extend the feature constructor while
keeping the same schema, fingerprint, and sampled-frame occupancy interfaces.

## 12. Exercises

### Calculation

1. For ten frames with a feature present in seven, what is its sampled-frame
   occupancy?
2. Two frames are invalid. Repeat the calculation using only valid frames.
3. A feature appears in 4 of 6 valid frames. What are its sampled-frame
   occupancy and count?
4. Two identical observations occur in one frame. What are the bit and count?
5. Build a schema with three features and write the binary vectors for three
   synthetic frames.

### Reasoning

1. Why is a residue's array position not a stable feature identity?
2. Why does a shared schema matter for comparing frames?
3. Why is sampled-frame occupancy not an equilibrium probability?
4. Why are adjacent MD frames not independent observations?
5. Why are missing frames different from frames with no interactions?
6. Why should detection criteria travel with the sampled-frame occupancy
   result?
7. Why is a fingerprint not an energy decomposition?

### Implementation

1. Run `examples/drug_discovery/interaction_fingerprints.py` and reproduce every sampled-frame
   occupancy by hand.
2. Add a feature that is present in every valid frame and verify sampled-frame
   occupancy 1.
3. Add a feature that is never present and verify sampled-frame occupancy 0.
4. Add two missing frames and confirm the denominator changes.
5. Reorder the interaction records and confirm the schema is unchanged.
6. Compare two synthetic binding modes using a shared schema.

## 13. Summary

- An individual interaction, a fingerprint, and a sampled-frame occupancy are
  different data products.
- Binary fingerprints record presence; count fingerprints record repeated
  observations.
- Feature identity must be stable and based on structural identifiers, not
  array positions.
- A shared schema is required for cross-frame comparison.
- Sampled-frame occupancy is the fraction of valid frames in which a feature is
  present.
- Missing frames are excluded from the denominator and reported separately.
- Adjacent MD frames are temporally correlated; sampled-frame occupancy is
  descriptive unless uncertainty is estimated properly.
- The source of the frames belongs in metadata, not in the sampled-frame
  occupancy algorithm.
- The next module is docking fundamentals and pose interpretation.

**Previous:** [03 — Protein–ligand interaction geometry](03_protein_ligand_interactions.md).

**Next:** [05 — Docking fundamentals and pose interpretation](05_docking_fundamentals.md).
