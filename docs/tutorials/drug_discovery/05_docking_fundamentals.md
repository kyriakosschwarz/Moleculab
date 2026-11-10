# Track C, module 05 — Docking fundamentals and pose interpretation

> **Where this fits.** [Module 04](04_interaction_fingerprints.md) analysed
> interactions across frames. This module asks how docking generates candidate
> poses, how those poses are scored and ranked, and how to interpret them
> without confusing a docking score with a binding free energy.
>
> **What you should already know.** Track C modules 01–04; Track A's docking
> examples; Track B protonation and tautomerism.
>
> **What you will be able to do afterwards.** Explain the docking search
> problem, distinguish search from scoring, evaluate a pose with geometry and
> chemical evidence, and compute a symmetry-aware ligand RMSD.

## 1. Learning objectives

**Core**

- Explain the purpose of molecular docking.
- Distinguish pose generation from pose scoring and ranking.
- Describe rigid and flexible docking.
- Explain how ligand and receptor flexibility expand the search space.
- Interpret docking scores without treating them as binding free energies.
- Compare a predicted pose with a reference structure using RMSD.

**Deeper insight**

- Search quality and scoring quality are separate failure modes.
- A docking score is a model-dependent ranking quantity.
- Symmetry-aware atom mapping is essential for meaningful ligand RMSD.
- Pose interpretation needs geometric, chemical, and conformational evidence.
- Receptor alignment and ligand alignment are different operations.

**Looking ahead**

- Module 06 will compare poses and binding modes.
- Module 07 will cover docking validation and benchmarking.
- The capstone will compare ligands using descriptors, protonation states,
  binding-site geometry, interactions, fingerprints, and pose recovery.

## 2. What docking attempts to solve

Docking addresses two related but different questions:

1. **Pose prediction:** where and how might a ligand bind?
2. **Affinity prediction:** how strongly does it bind?

Moleculab's existing docking examples address the first question with a simple
search and an interaction-energy score. They do not provide a validated
affinity predictor.

A docking run typically produces:

- one or more candidate poses;
- a score for each pose;
- a ranking;
- sometimes a predicted binding mode and a confidence estimate.

The score is not automatically a binding free energy. It is a ranking
quantity produced by a particular scoring model.

## 3. Input representation and preparation

Docking output quality is bounded by input quality.

Before a search, decide:

- which protein structure to use;
- which chains form the biological assembly;
- whether missing atoms or loops matter;
- which waters, ions, and cofactors to retain;
- the protonation state of protein residues;
- the ligand's protonation and tautomer state;
- the ligand's stereochemistry;
- whether the ligand is rigid or flexible;
- whether the receptor is rigid or flexible.

These are not minor technical details. A docking run with the wrong tautomer or
a missing metal ion can produce a plausible-looking pose for the wrong chemical
system.

The structure adapter and chemical-feature layers from Track C modules 01–03
provide the information needed to make these choices explicit.

## 4. The docking search problem

A ligand pose can be described by:

$$ q=(\mathbf t,\boldsymbol\omega,\boldsymbol\phi), $$

where:

- $\mathbf t$ is translation in three dimensions;
- $\boldsymbol\omega$ represents orientation, with three rotational degrees of
  freedom;
- $\boldsymbol\phi$ represents internal ligand torsion angles.

For a ligand with $n$ rotatable torsions, this gives approximately

$$ 3+3+n $$

search variables before considering receptor flexibility.

The search space is therefore large and highly non-convex. A docking algorithm
must explore translations, rotations, and internal conformations while avoiding
steric clashes and finding chemically plausible contacts.

### 4.1 Search versus scoring

Keep three operations separate:

| Operation | Question |
| --- | --- |
| search | which candidate poses should be evaluated? |
| scoring | how favourable is each candidate under the model? |
| ranking | which candidates are reported first? |

A better search algorithm cannot automatically compensate for a systematically
inaccurate scoring function. Conversely, a good scoring function cannot recover
a pose that the search never generated.

## 5. Ligand flexibility

A rigid ligand has only six external degrees of freedom. A flexible ligand adds
rotatable bonds.

More torsions mean:

- a larger search space;
- more conformations to sample;
- a greater chance of missing the relevant conformation;
- a greater need for conformational strain penalties or filtering.

Ligand conformational preferences come from the same force-field and
conformational ideas as Track A. A pose that is geometrically possible but
requires an implausible torsion may be a poor hypothesis.

## 6. Receptor flexibility

Rigid-receptor docking keeps the protein fixed. Flexible-receptor docking allows
some protein atoms or side chains to move.

Approximations include:

- rigid receptor;
- selected side-chain flexibility;
- multiple receptor conformations;
- induced-fit or ensemble docking;
- fully flexible receptor search in specialised methods.

Moleculab's Track A examples include a simple rigid-ligand search and a
side-chain displacement model with a harmonic strain penalty. Those are teaching
models, not production docking engines.

## 7. Scoring functions

A conceptual docking score can be written as

$$ S_{\mathrm{dock}} = E_{\mathrm{inter}} + E_{\mathrm{intra}} + E_{\mathrm{strain}} + E_{\mathrm{other}}. $$

This is a conceptual decomposition, not a universal docking equation. Different
programs use different terms and functional forms.

Common contributions include:

- van der Waals complementarity;
- electrostatics;
- hydrogen-bond geometry;
- hydrophobic contacts;
- ligand conformational strain;
- torsional penalties;
- desolvation approximations;
- sometimes empirical or knowledge-based terms.

The central distinction is

$$ S_{\mathrm{dock}} \ne \Delta G_{\mathrm{bind}}. $$

Even when a score is reported in kJ/mol or kcal/mol, it remains a
model-dependent ranking quantity unless the method explicitly constructs and
validates a thermodynamic free-energy estimate.

A docking score may correlate with affinity within a congeneric series. That is
useful. It is not the same as a measured or rigorously computed binding free
energy.

## 8. Pose generation and ranking

A docking program may use:

- stochastic search;
- genetic algorithms;
- Monte Carlo simulated annealing;
- incremental construction;
- shape or pharmacophore matching;
- local optimisation;
- combinations of the above.

The search proposes poses. The scoring function evaluates them. The ranking
orders the candidates.

Failure can occur at any stage:

- the search misses the correct pose;
- the scoring function prefers a wrong pose;
- the receptor conformation is wrong;
- the ligand state is wrong;
- the binding site includes or excludes the wrong waters or cofactors;
- the reported pose is a local minimum of the scoring model.

A "best" pose is the best pose under the model, not a guarantee of the correct
binding mode.

## 9. Pose interpretation

Do not rely on the docking score alone. Use independent evidence.

### 9.1 Geometric plausibility

- Does the ligand fit the pocket?
- Are there severe steric clashes?
- Are hydrogen-bond distances and angles plausible?
- Are hydrophobic contacts consistent with the pocket?
- Does the ligand occupy the expected subpockets?

Reuse [binding-site geometry](02_binding_sites.md) and
[interaction geometry](03_protein_ligand_interactions.md).

### 9.2 Chemical plausibility

- Is the protonation state consistent with the environment?
- Is the tautomer consistent with the hydrogen-bond pattern?
- Are donor and acceptor assignments plausible?
- Is the stereochemistry correct?
- Are charged groups desolvated in a reasonable way?

Reuse Track B modules 02 and 04.

### 9.3 Conformational plausibility

- Is the ligand conformation strained?
- Are torsions near preferred values?
- Does the pose require an unrealistic intramolecular contact?
- Is the conformation consistent with the ligand's accessible ensemble?

Reuse Track A force-field and conformer ideas.

### 9.4 Validation evidence

- Does the pose reproduce a known experimental binding mode?
- Is the result stable across reasonable preparation choices?
- Does the same ligand produce the same pose in repeated runs?
- Does the pose survive a different scoring function?

A pose supported by several independent observations is a stronger hypothesis
than a pose supported only by one score.

## 10. Ligand RMSD and symmetry

For predicted and reference ligand coordinates in the same frame, a basic RMSD
is

$$ \mathrm{RMSD} = \sqrt{\frac{1}{N}\sum_{i=1}^{N}\lVert \mathbf r_i - \mathbf r_i^{\mathrm{ref}}\rVert^2}. $$

Several details matter:

- Are the protein structures aligned first?
- Are the pose and reference in the same coordinate frame?
- Are ligand atoms mapped correctly?
- Are chemically equivalent atoms interchangeable?
- Is the mapping complete?
- Is the comparison symmetry-aware?

### 10.1 Receptor alignment versus ligand RMSD

If two receptor structures are not aligned, the ligand RMSD is meaningless.
Receptor alignment is a separate operation. The pose-analysis functions assume
that both coordinate sets are already in a common frame.

Rigid-body invariance has a precise meaning: applying the **same** rigid
transformation to both the pose and the reference leaves their RMSD unchanged.
Applying a transformation to only one coordinate set generally changes the
RMSD, because the relative displacement between the two sets has changed.
Alignment is the operation that removes an arbitrary frame difference.

### 10.2 Symmetry-aware mapping

A symmetric ligand can have an equivalent pose that produces a large RMSD under
an arbitrary atom mapping. The mapping must be supplied explicitly.

This module supports symmetry-aware comparison by accepting several valid
one-to-one mappings and returning the minimum RMSD. It does not perform graph
automorphism perception itself.

### 10.3 Missing or incompatible atoms

Do not silently drop atoms to make an RMSD smaller. The mapping must be
complete by default. A partial comparison is allowed only when the caller
explicitly requests it.

## 11. Pose analysis implementation

This module does not implement a docking search engine. It analyses externally
generated poses.

### 11.1 Pose record

`LigandPose` stores:

- positions in nm;
- stable, unique atom identifiers;
- optional docking score;
- score unit and score source when a score is present;
- source and metadata.

The score is metadata, not a binding free energy.

### 11.2 API

| Function or object | Purpose |
| --- | --- |
| `LigandPose` | coordinates, atom IDs, score provenance |
| `atom_id_mapping(pose, reference)` | map atoms by stable identifiers |
| `pose_rmsd(...)` | RMSD for an explicit mapping |
| `symmetry_aware_rmsd(...)` | minimum RMSD over explicit symmetry mappings |
| `evaluate_pose_recovery(...)` | RMSD threshold result with score provenance |
| `PoseRecovery` | RMSD, threshold, pass flag, best mapping, score metadata |

### 11.3 Example

```bash
uv run python examples/drug_discovery/docking_pose_analysis.py
```

Observed output:

```text
Pose comparison
  atom mapping by ID: ((0, 0), (1, 1), (2, 2))
  RMSD with that mapping: 0.3268 nm
  symmetry-aware RMSD: 0.0115 nm (mapping 2 of 2)

Pose recovery
  RMSD: 0.0115 nm
  threshold: 0.050 nm
  passed: True

Score provenance
  score: -8.1
  unit: arbitrary
  source: example ranker
  This is a docking score, not a binding free energy.

Pocket around the docked pose
  selected receptor atoms: [0, 1, 2, 3]
  centre: (0.000, 0.050, 0.000) nm
```

The identity mapping is large because the two carboxylate oxygens are swapped.
The explicit symmetry mapping recovers the equivalent pose.

### 11.4 Validation and acceptance criteria

```bash
uv run pytest tests/analysis/test_docking_poses.py -v
```

The tests cover:

- RMSD zero for identical mapped coordinates;
- a hand-calculated RMSD;
- invariance to a consistent rigid-body transformation;
- independence from array order with an explicit mapping;
- atom-ID mapping;
- missing atom rejection;
- invalid mapping rejection;
- explicit partial-mapping opt-in;
- symmetry-aware best mapping;
- pose recovery and threshold behaviour;
- score-unit and provenance validation.

A pose-analysis module is correct when:

- RMSD is zero for identical mapped coordinates;
- a consistent rigid-body transformation applied to both coordinate sets
  leaves RMSD unchanged;
- atom mappings are validated before calculation;
- symmetry-equivalent mappings are explicit;
- missing or incompatible atoms are not silently ignored;
- score units and provenance are preserved;
- receptor alignment is treated separately;
- a docking score is never labelled a binding free energy without a justified
  thermodynamic method.

A pose-recovery threshold is a protocol choice, not a universal definition of
success. The example threshold of 0.05 nm is illustrative. Valid thresholds
depend on ligand size, atom mapping, symmetry handling, the validation
protocol, and the purpose of the benchmark.

Pose recovery, chemical plausibility, and affinity prediction are different
questions:

| Question | Measure |
| --- | --- |
| Did the method reproduce a reference geometry? | pose RMSD and recovery threshold |
| Are the interactions chemically reasonable? | interaction geometry and fingerprints |
| Does the method estimate experimental binding strength? | affinity benchmarks |

A method can do well on one and poorly on another. A low RMSD does not by
itself establish correct interactions, and a correct pose does not guarantee a
good affinity prediction.

## 12. Limitations

Docking remains an approximation. Important limitations include:

- rigid or partially flexible receptor;
- incomplete ligand conformational sampling;
- simplified scoring functions;
- implicit or absent solvation;
- ignored or incompletely modelled waters;
- protonation and tautomer uncertainty;
- induced fit and allosteric effects;
- metal coordination and covalent binding;
- entropy and desolvation;
- score calibration and domain dependence;
- no guarantee that the top-ranked pose is correct.

A docking score is a hypothesis generator. It is not a measured affinity.

## 13. Exercises

### Reasoning

1. Why is pose prediction different from affinity prediction?
2. Why can a better search algorithm not fix a biased scoring function?
3. Why does ligand flexibility increase the search space?
4. Why is a docking score not a binding free energy?
5. Why must receptor alignment precede ligand RMSD comparison?
6. Why is symmetry-aware atom mapping necessary?
7. Why is a top-ranked pose not proof of the correct binding mode?

### Calculation

1. For a ligand with five rotatable torsions, how many search variables are
   there before receptor flexibility?
2. Compute the RMSD for two atoms displaced by 0.1 nm and 0.2 nm.
3. Compute the RMSD for an identical mapped pose.
4. Given a pose with RMSD 0.18 nm and a recovery threshold of 0.20 nm, does it
   pass?
5. Construct two symmetry mappings for a ligand with two equivalent oxygens.

### Implementation

1. Run `examples/drug_discovery/docking_pose_analysis.py` and reproduce the output.
2. Swap the two equivalent ligand atoms in the reference and confirm the
   symmetry-aware RMSD remains small.
3. Apply the same rotation and translation to both pose and reference and
   confirm RMSD is unchanged.
4. Remove an atom ID and confirm that mapping fails rather than silently
   dropping it.
5. Add score metadata and confirm that unit and source are retained.
6. Reuse `binding_site_geometry` to describe the pocket around the docked pose.

## 14. Summary

- Docking generates and ranks candidate ligand poses.
- Search, scoring, and ranking are separate operations.
- Ligand and receptor flexibility expand a large, non-convex search space.
- A docking score is a model-dependent ranking quantity, not automatically a
  binding free energy.
- Pose interpretation should combine pocket geometry, interaction geometry,
  chemical state, conformational plausibility, and validation evidence.
- Ligand RMSD requires explicit atom mapping and a common coordinate frame.
- Symmetry-aware RMSD uses explicit equivalent mappings.
- The implementation analyses external poses; it does not contain a docking
  search engine.
- The next step is comparing poses and binding modes across ligands.

**Previous:** [04 — Interaction fingerprints and sampled-frame occupancy](04_interaction_fingerprints.md).

**Next:** [06 — Comparing binding modes](06_comparing_binding_modes.md).
