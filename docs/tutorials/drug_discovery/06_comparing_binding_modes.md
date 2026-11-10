# Track C, module 06 — Comparing binding modes

> **Where this fits.** [Module 05](05_docking_fundamentals.md) compared one
> predicted pose with a reference. This module compares poses and interaction
> patterns across ligands or across multiple poses.
>
> **What you should already know.** Track C modules 02–05, especially pose RMSD
> and interaction fingerprints.
>
> **What you will be able to do afterwards.** Cluster poses by explicit
> geometric distance, compare pose placement, compare interaction patterns at
> exact and residue levels, and explain when geometric similarity and
> interaction similarity disagree.

## 1. Learning objectives

**Core**

- Explain what it means for two ligands to share a binding mode.
- Cluster poses using pairwise RMSD and a threshold.
- Compare pose geometry using RMSD and centroid displacement.
- Compare interaction fingerprints at exact-feature and residue-pattern levels.
- Distinguish geometric similarity from interaction similarity.

**Deeper insight**

- A shared pose does not imply a shared interaction pattern.
- A shared residue interaction pattern does not imply the same geometry.
- Single-linkage clustering can chain through intermediate poses.
- Clustering and comparison thresholds are protocol choices.
- A composite "binding-mode score" would hide the distinction the analysis is
  trying to teach.

**Looking ahead**

- Module 07 will cover docking validation and benchmarking.
- The capstone will compare two ligands using chemistry, geometry, interactions,
  and pose recovery.

## 2. What is a binding mode?

A binding mode describes how a ligand is placed in a binding site and which
interactions it makes.

It has at least two parts:

| Part | Example evidence |
| --- | --- |
| geometry | pose RMSD, centroid displacement, occupied subpockets |
| interaction pattern | residues contacted, interaction types, directionality |

Two ligands can:

- occupy the same region and make the same interactions;
- occupy the same region and make different interactions;
- occupy different regions and contact the same residues;
- differ in both geometry and interactions.

Those are different scientific claims. Report them separately.

## 3. Pose clustering

Pose clustering groups structures that are geometrically similar under an
explicit atom mapping.

For poses $i$ and $j$, define a distance

$$ d_{ij}=\mathrm{RMSD}_{ij}. $$

For a threshold $r_c$, connect poses when

$$ d_{ij}\le r_c. $$

A simple single-linkage clustering treats each connected component as a
cluster.

### 3.1 What clustering does and does not say

Clustering says:

> under this atom mapping, this distance measure, and this threshold, these
> poses form a connected group.

It does not say:

- that the poses share the same interactions;
- that the cluster represents one thermodynamic state;
- that the cluster is the correct binding mode;
- that a different threshold would produce the same grouping.

### 3.2 Single-linkage chaining

Single-linkage clustering can chain through intermediate poses. If A is close
to B and B is close to C, then A and C may end up in the same cluster even when
their direct distance is larger than the threshold.

This is a property of the method, not a bug. Document it, and inspect the
distance matrix when the cluster count matters.

### 3.3 Symmetry

Symmetry-aware clustering requires an explicit distance function that considers
the relevant equivalent mappings. The default `pose_distance` maps atoms by
stable identifier and is not symmetry-aware.

A symmetric ligand can appear in two clusters under the default distance and in
one cluster under an explicit symmetry mapping.

## 4. Comparing pose geometry

The primary geometric measure is RMSD under an explicit mapping:

$$ \mathrm{RMSD}=\sqrt{\frac{1}{N}\sum_{i=1}^{N}\lVert\mathbf r_i-\mathbf r_i^{\mathrm{ref}}\rVert^2}. $$

RMSD captures placement, orientation, and internal shape together. It does not
separate them.

A centroid displacement is a simpler translation measure:

$$ \mathbf d_{\mathrm{centroid}}=\bar{\mathbf r}_A-\bar{\mathbf r}_B. $$

The centroid distance is the norm of that vector. It ignores orientation and
internal conformation, so it should be reported alongside RMSD, not instead of
it.

### 4.1 Common frame

Both poses must be in the same coordinate frame. Receptor alignment is a
separate operation. Applying the same rigid transformation to both poses leaves
the RMSD unchanged; transforming only one pose generally changes it.

## 5. Comparing interaction patterns

An interaction fingerprint represents which features are present in a pose.
Two comparison levels are useful.

### 5.1 Exact feature comparison

Compare full `InteractionFeature` identities:

```text
(interaction_type, protein_residue, ligand_residue, ligand_atom_id)
```

This level is appropriate when comparing poses of the same ligand or when the
ligand atom identifiers are directly comparable.

### 5.2 Residue-pattern comparison

Project each feature to:

```text
(interaction_type, protein_residue)
```

This level is appropriate for comparing different ligands that may use
different atoms to contact the same residue.

It answers a different question:

> do these ligands make the same kind of interaction with the same protein
> residue?

It does not say that the ligand atoms or the local geometry are the same.

### 5.3 Jaccard similarity

For binary feature sets $A$ and $B$,

$$ J(A,B)=\frac{\lvert A\cap B\rvert}{\lvert A\cup B\rvert}. $$

The value ranges from 0 to 1. If both sets are empty, the implementation
returns 1.0 because the observed patterns are identical.

Jaccard similarity is a set comparison. It does not weight interaction strength,
occupancy, or geometry. Two patterns can have identical binary fingerprints
despite different contact distances, hydrogen-bond angles, or interaction
energies. That is a limitation of the representation, not a defect in the
implementation.

## 6. The interpretation matrix

The two comparison axes produce four cases:

| Geometry | Interaction pattern | Interpretation |
| --- | --- | --- |
| similar | similar | same binding-mode hypothesis |
| similar | different | same placement, different chemistry |
| different | similar | different placement, same residue contacts |
| different | different | different binding modes |

This table is a guide, not a classification algorithm.

The most interesting cases are the off-diagonal ones:

- similar geometry with different interactions means the ligand is in the same
  region but presents different chemistry;
- similar interactions with different geometry means the same residues are
  engaged from a different arrangement.

Neither case is adequately described by a single "binding-mode score".

## 7. Comparing ligands across a shared binding site

When comparing different ligands:

1. Use the same prepared protein structure or explicitly aligned structures.
2. Define the binding site consistently.
3. Choose the same interaction criteria.
4. Compare residue-level patterns, not ligand atom IDs.
5. Report pose geometry separately from interaction patterns.
6. Check protonation, tautomer, and stereochemical assumptions.
7. Do not interpret a similarity value as an affinity difference.

If two ligands use different atom identifiers, exact-feature comparison will
show little or no overlap even when the residue-level pattern is identical.
That is why both levels matter.

### 7.1 Within-ligand versus between-ligand comparison

Within one ligand, atom correspondence and symmetry mappings can be defined
consistently across poses. Atom-mapped RMSD is meaningful.

Between different ligands, atom-by-atom RMSD may be meaningless when the
molecules have different connectivity, atom counts, or atom identities. In that
case, compare spatial placement, pharmacophoric features, or residue-level
interaction patterns rather than assuming an atom-to-atom correspondence.

## 8. Clustering and comparison thresholds

Thresholds are protocol choices, not universal constants.

They depend on:

- ligand size;
- atom mapping and symmetry;
- the reference structure;
- the interaction criteria;
- the purpose of the comparison;
- whether the question is pose recovery, mode clustering, or interaction
  comparison.

A common mistake is to reuse one threshold for every question. A pose-recovery
threshold, a clustering threshold, and an interaction-distance cutoff answer
different questions even when they have the same units.

## 9. Limitations

- RMSD is mapping-dependent.
- Symmetry-aware comparison requires explicit mappings.
- Centroid distance ignores orientation and shape.
- Fingerprint Jaccard ignores interaction strength and geometry.
- Residue-level patterns hide ligand-atom differences.
- Single-linkage clustering can chain.
- Different protonation or tautomer states produce different patterns.
- A shared interaction pattern is not a shared binding free energy.
- Clustering is not a thermodynamic decomposition.

## 10. Implementation and validation

The implementation is in `src/moleculab/analysis/docking/`.

### 10.1 API

| Object or function | Purpose |
| --- | --- |
| `PoseCluster` | member indices and representative |
| `PoseClusteringResult` | clusters, distance matrix, threshold, labels |
| `pose_distance(...)` | RMSD under identifier or explicit mappings |
| `cluster_poses(...)` | single-linkage connected components |
| `PoseGeometryComparison` | RMSD, centroid displacement, best mapping |
| `compare_pose_geometry(...)` | geometry comparison |
| `InteractionPatternComparison` | shared and unique patterns, Jaccard |
| `compare_interaction_patterns(...)` | exact or residue-level comparison |
| `BindingModeComparison` | geometry and interaction results separately |
| `compare_binding_modes(...)` | combined call without a composite score |

### 10.2 Example

```bash
uv run python examples/drug_discovery/compare_binding_modes.py
```

Observed output:

```text
Same geometry, different interactions
  RMSD: 0.000 nm
  exact feature Jaccard: 0.000

Different geometry, same residue-level pattern
  RMSD: 0.500 nm
  residue-pattern Jaccard: 1.000

Residue-level shared features
  shared: (('hydrogen_bond', 'ASP A 123'),)
  only in A: ()
  only in B: ()

Pose clustering (single linkage)
  threshold: 0.20 nm
  cluster 0: ['A', 'B']
  cluster 1: ['C']

Pose geometry
  RMSD: 0.500 nm
  centroid distance: 0.500 nm
```

The first case has identical coordinates but different interaction patterns.
The second has different coordinates but the same residue-level pattern. The
example is designed to make the distinction impossible to miss.

### 10.3 Tests

```bash
uv run pytest tests/test_pose_comparison.py -v
```

The tests cover:

- geometry comparison and centroid displacement;
- symmetry-aware geometry comparison;
- exact feature comparison;
- residue-pattern comparison across different ligand atom IDs;
- low RMSD with different interaction patterns;
- high RMSD with the same residue pattern;
- clustering by distance;
- input-order-independent partitions;
- symmetry-aware clustering;
- atom-ID mismatch rejection;
- single-linkage chaining;
- threshold and empty-input validation.

### 10.4 Acceptance criteria

A binding-mode comparison module is correct when:

- pose geometry and interaction patterns are reported separately;
- RMSD uses an explicit atom mapping;
- symmetry-equivalent mappings are explicit;
- centroid displacement is documented as a translation-only measure;
- exact and residue-level interaction comparisons are distinguished;
- Jaccard similarity is computed from binary feature presence;
- clusters are deterministic for a fixed distance matrix and threshold;
- single-linkage chaining is documented;
- no composite score hides the difference between geometry and interactions.

## 11. Exercises

### Reasoning

1. Why can two ligands have low RMSD but different interaction patterns?
2. Why can two ligands have high RMSD but the same residue-level pattern?
3. Why is a single binding-mode score potentially misleading?
4. Why does single-linkage clustering chain?
5. Why does symmetry-aware clustering require explicit mappings?
6. Why is centroid distance insufficient as a pose comparison?
7. Why should interaction criteria match when comparing fingerprints?

### Calculation

1. For feature sets $\{A,B,C\}$ and $\{B,C,D\}$, compute the Jaccard similarity.
2. For a cluster with poses A, B, and C, and distances 0.1, 0.2, and 0.35 nm,
   which poses join at a 0.25 nm threshold under single linkage?
3. Compute the centroid displacement for two poses separated by 0.3 nm in x.
4. For two identical poses with different fingerprints, what are the RMSD and
   exact-feature Jaccard?
5. For two shifted poses with the same residue pattern, what are the RMSD and
   residue-pattern Jaccard?

### Implementation

1. Run `examples/drug_discovery/compare_binding_modes.py` and reproduce the output.
2. Add a third interaction feature and recompute the exact Jaccard value.
3. Change the clustering threshold and identify when the chain breaks.
4. Add a symmetry-aware distance function and cluster two equivalent poses.
5. Compare two poses of the same ligand with exact features and then with
   residue-level patterns.

## 12. Summary

- A binding mode has geometric and interaction components.
- Pose clustering groups similar geometries under an explicit distance and
  threshold.
- Single-linkage clustering is deterministic but can chain.
- Symmetry-aware comparison requires explicit mappings.
- Exact feature comparison and residue-pattern comparison answer different
  questions.
- Jaccard similarity compares binary interaction patterns, not energies.
- Similar geometry with different interactions, and different geometry with
  similar interactions, are both scientifically meaningful.
- Report geometry and interactions separately; do not hide them in one score.

**Previous:** [05 — Docking fundamentals and pose interpretation](05_docking_fundamentals.md).

**Next:** [07 — Docking validation and benchmarking](07_docking_validation.md).
