# Track C — Structure-based drug discovery

> **Question.** What does the physics mean for a drug?
>
> **Entry point.** This page.
>
> **Prerequisites.** Track A sections 1–14, especially nonbonded interactions,
> water, and electrostatics; Track B modules 01–04 for chemical identity and
> protonation state.
>
> **Status.** All ten structural and thermodynamic modules are available. The
> [imatinib–BCR-ABL1 capstone](../capstone_imatinib.md) is the integration
> exercise, and its [visual learning layer](visualizations/README.md) now
> includes molecular, binding-site, evidence, sensitivity, and thermodynamic
> figures.

Track C is where the engine becomes a molecular-recognition argument. A
protein structure is not just a large collection of particles, and a docking
score is not automatically a binding free energy. The track therefore starts
with geometry and representation before it introduces scoring.

## 1. What this track adds

| Question | Track C material |
| --- | --- |
| What does a protein structure contain? | chains, residues, coordinates, missing atoms, alternate locations |
| Where is the binding site? | pockets, cavities, waters, shape and volume |
| Which residues contact the ligand? | hydrogen bonds, salt bridges, hydrophobic and aromatic contacts |
| Which contacts persist? | residue-level interaction fingerprints and sampled-frame occupancy |
| How are poses generated? | rigid and flexible docking |
| What does a score mean? | interaction energy versus binding free energy |

## 2. Modules

| Module | Focus | Status |
| --- | --- | --- |
| [01 — Protein structure fundamentals](01_protein_structure.md) | chains, residues, coordinates, structure quality | available |
| [02 — Binding sites and cavities](02_binding_sites.md) | pocket geometry, waters, volume, complementarity | available |
| [03 — Protein–ligand interaction geometry](03_protein_ligand_interactions.md) | hydrogen bonds, salt bridges, hydrophobic and aromatic contacts | available |
| [04 — Interaction fingerprints](04_interaction_fingerprints.md) | residue-ligand contact features and sampled-frame occupancy | available |
| [05 — Docking foundations](05_docking_fundamentals.md) | search, scoring, poses, and limitations | available |
| [06 — Comparing binding modes](06_comparing_binding_modes.md) | pose clustering, geometry, interaction-pattern comparison | available |
| [07 — Docking validation and benchmarking](07_docking_validation.md) | redocking, cross-docking, pose recovery, score interpretation | available |
| [08 — Trajectory contact analysis](08_trajectory_contact_analysis.md) | contact matrices, occupancy, and lifetimes | available |
| [09 — Lead optimisation and SAR](09_lead_optimisation_sar.md) | matched pairs, activity scales, ligand efficiency, and evidence categories | available |
| [10 — Structural observations to thermodynamics](10_structural_observations_to_thermodynamics.md) | occupancy, ensembles, free energy, and why contacts are not affinities | available |
| [Visual learning layer](visualizations/README.md) | reproducible molecular, evidence, sensitivity, and thermodynamic figures | Batches 1–3 available |

Start with [module 01](01_protein_structure.md). It establishes the structural
representation and import boundary that the later modules assume.

### 2.1 Available now

The docking foundations live in the Track A main tutorial:

- [Rigid-ligand docking](../physics/tutorial.md#15-rigid-ligand-docking) introduces
  receptor-ligand interaction scoring, translation/rotation sampling, and the
  distinction between an interaction energy and a binding free energy.
- [Flexible-receptor docking](../physics/tutorial.md#24-flexible-receptor-docking) adds
  explicit side-chain displacement and a harmonic strain penalty.

The supporting physics is in:

- [Atomistic water, atom types, and exclusions](../physics/tutorial.md#14-atomistic-water-atom-types-and-exclusions);
- [Ewald summation for periodic charges](../physics/tutorial.md#16-ewald-summation-for-periodic-charges);
- [Trajectory analysis](../physics/tutorial.md#11-trajectory-analysis).

The applied module code is in `src/moleculab/applications/docking.py` and
`src/moleculab/applications/flexible_docking.py`. Examples are
`examples/drug_discovery/rigid_docking.py` and `examples/drug_discovery/flexible_docking.py`.

### 2.2 Module 01: protein structure fundamentals

[Module 01](01_protein_structure.md) teaches how to inspect a protein structure
before scoring anything:

- PDB and mmCIF files, chains, residues, and atom names;
- experimental versus predicted structures;
- missing atoms, alternate locations, and incomplete side chains;
- hydrogen atoms and protonation states;
- water, ions, cofactors, and ligand records;
- symmetry mates and biological assemblies;
- what a structure does not contain: missing loops, disorder, and uncertainty.

The acceptance criterion is representational: a reader should be able to
identify the chains, residues, ligands, waters, and missing regions in a real
structure, and explain which of them are safe to use in a calculation.

Runnable inspection: `uv run --with gemmi python examples/drug_discovery/protein_structure.py`.

### 2.3 Module 02: binding sites and cavities

[Module 02](02_binding_sites.md) treats the binding site as a geometric object:

- pocket detection and cavity volume;
- buried versus exposed surface;
- shape complementarity;
- conserved and structural waters;
- the difference between an ATP site, an allosteric site, and a surface groove.

The first implementation is deliberately modest: minimum protein–ligand
distances, cutoff-based atom selection, residue aggregation, and a pocket centre
and bounding box. It is structural analysis, not an energy calculation.

### 2.4 Module 03: protein–ligand interaction geometry

[Module 03](03_protein_ligand_interactions.md) identifies interactions from
coordinates before assigning energies:

- donor–acceptor distances and hydrogen-bond angles;
- salt-bridge distances and charge complementarity;
- hydrophobic contacts and aromatic stacking geometry;
- steric clashes;
- metal coordination.

A geometric contact says that an interaction is possible in that structure. It
does not say how much that contact contributes to binding.

This milestone implements the chemical-feature record and hydrogen-bond
detector. The remaining interaction types are specified as next detectors so
that each can have its own geometry and tests.

### 2.5 Module 04: interaction fingerprints

[Module 04](04_interaction_fingerprints.md) records which ligand atoms contact
which protein residues and with which interaction type. Over a collection of
frames, the same representation supports:

- sampled-frame occupancy;
- hydrogen-bond persistence;
- residue interaction timelines;
- pose comparison;
- ligand-induced changes in the binding site.

The essential distinction is:

> A geometric interaction fingerprint is not an energy decomposition.

### 2.6 Module 05: docking foundations

[Module 05](05_docking_fundamentals.md) treats docking as pose generation and
ranking, not as a direct binding-free-energy method:

- rigid-body search;
- scoring with the typed nonbonded model;
- translation and rotation sampling;
- flexible side chains and strain;
- pose interpretation and symmetry-aware RMSD;
- why the best score is a hypothesis, not a measured affinity.

The existing Track A docking examples remain the simple search/scoring
demonstrations. Module 05 adds the conceptual and pose-analysis layer.

### 2.7 Module 06: comparing binding modes

[Module 06](06_comparing_binding_modes.md) compares poses and interaction
patterns:

- single-linkage pose clustering with explicit distance functions;
- symmetry-aware clustering;
- RMSD and centroid displacement;
- exact-feature and residue-pattern Jaccard comparison;
- cases where geometric similarity and interaction similarity disagree.

It reports geometry and interaction patterns separately and does not create a
single binding-mode score.

### 2.8 Module 07: docking validation and benchmarking

[Module 07](07_docking_validation.md) covers redocking, cross-docking,
pose-recovery thresholds, benchmark composition, selection bias, protocol
metadata, and the distinction between pose recovery, chemical plausibility, and
affinity prediction.

The benchmark layer records valid failures and invalid cases separately. The
success rate uses only valid cases in the denominator; invalid cases remain
visible in the report.

### 2.9 Module 08: trajectory contact analysis

[Module 08](08_trajectory_contact_analysis.md) analyses geometric contacts over
coordinate frames:

- atom-level and residue-level contact matrices;
- sampled-frame contact occupancy;
- contact run lengths and lifetimes;
- the difference between a continuous contact and many short contacts;
- missing-frame handling and sampling limitations.

It is parser-independent and does not load trajectory files.

### 2.10 Module 09: lead optimisation and SAR

[Module 09](09_lead_optimisation_sar.md) connects structural evidence to
medicinal-chemistry decisions:

- matched molecular pairs and activity changes;
- pIC50, IC50, Ki, and Kd;
- physicochemical descriptors;
- ligand efficiency and its limitations;
- docking and interaction evidence as mechanistic hypotheses;
- selectivity, liabilities, and data-category separation.

The comparison layer is parser-independent; RDKit or another chemistry toolkit
provides descriptors and matched-pair detection at the boundary.

### 2.11 Module 10: structural observations to thermodynamics

[Module 10](10_structural_observations_to_thermodynamics.md) is the conceptual
bridge to Track A's free-energy material:

- structural observations versus thermodynamic quantities;
- ensemble populations and the partition function;
- energy, degeneracy, and entropy;
- desolvation, water reorganisation, and conformational entropy;
- binding free energy and the dissociation constant;
- why occupancy and contact lifetime are not binding free energy;
- a two-ligand demonstration where contact count and free-energy components
  disagree.

It does not introduce a new free-energy method.

## 3. What you can do now

A useful Track C exercise with the current library is:

1. Read a receptor and a ligand from the docking examples.
2. Identify the receptor residues within a chosen distance of the ligand.
3. Classify the contacts by element and functional group.
4. Compare the contact set for two poses with different scores.
5. Explain why a better score does not prove a better binding free energy.

This exercise uses only geometry and the existing nonbonded model. It does not
require a new scoring function.

## 4. Connection to Track B

Track C needs the chemical state from Track B:

- a ligand's protonation changes its charge and hydrogen-bonding pattern;
- a tautomer changes which atoms can donate or accept hydrogen bonds;
- a protein residue's pKa can shift in a binding site;
- a charged group may be desolvated when it is buried.

Without a chosen chemical state, an interaction fingerprint is ambiguous. The
Track B modules [protonation and tautomerism](../chemistry/04_protonation_and_tautomerism.md)
and [molecular representation](../chemistry/02_molecular_representation.md)
are the relevant preparation.

## 5. Connection to Track A

Track C reuses:

- `System` and the typed nonbonded model for interaction energies;
- `Box` and Ewald electrostatics for periodic systems;
- minimisation and MD for relaxation and sampling;
- trajectory analysis for contact persistence;
- free-energy estimators for the thermodynamics behind a score.

The same warning applies throughout: the engine computes what its model says.
Track C interprets that output in a drug-discovery context, and must state the
model's limitations.

## 6. Reference

- [Units and conventions](../../reference/units.md) — Angstrom, nm, degrees, kJ/mol,
  kcal/mol, molar concentrations, pIC50, pH, pressure, and structural units.

## 7. Acceptance criteria

A Track C module is complete when:

- the protein representation is explicit about missing and modelled atoms;
- binding-site geometry is described without pretending it is an energy;
- interaction types are assigned from clear geometric criteria;
- a contact fingerprint can be reproduced from a known structure;
- sampled-frame occupancy is checked against a manually constructed trajectory;
- docking results are reported as poses and scores, not as binding affinities;
- the link to the relevant Track A and Track B foundations is explicit.

**Previous:** [Track B — Chemical structure and properties](../chemistry/README.md).

**Next:** [Capstone — imatinib and BCR-ABL1](../capstone_imatinib.md).
