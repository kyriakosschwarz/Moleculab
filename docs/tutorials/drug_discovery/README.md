# Track C — Structure-based drug discovery

> **Question.** What does the physics mean for a drug?
>
> **Entry point.** This page.
>
> **Prerequisites.** Track A sections 1–14, especially nonbonded interactions,
> water, and electrostatics; Track B modules 01–04 for chemical identity and
> protonation state.
>
> **Status.** The docking foundations are available. Protein-structure and
> interaction-fingerprint modules are planned. The
> [imatinib–BCR-ABL1 capstone](../capstone_imatinib.md) is the current
> integration exercise.

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
| Which contacts persist? | residue-level interaction fingerprints and trajectory occupancy |
| How are poses generated? | rigid and flexible docking |
| What does a score mean? | interaction energy versus binding free energy |

## 2. Modules

| Module | Focus | Status |
| --- | --- | --- |
| 01 — Protein structure fundamentals (planned) | chains, residues, coordinates, structure quality | planned |
| 02 — Binding sites and cavities | pocket geometry, waters, volume, complementarity | planned |
| 03 — Protein–ligand interaction geometry | hydrogen bonds, salt bridges, hydrophobic and aromatic contacts | planned |
| 04 — Interaction fingerprints | residue-ligand contact features and trajectory occupancy | planned |
| 05 — Docking foundations | search, scoring, poses, and limitations | available in Track A |
| 06 — Trajectory contact analysis | contact timelines and persistence | planned |

Module 01 is the recommended starting point once it is available. Its planned
scope is described below so the rest of the track has a stable target.

### 2.1 Available now

The docking foundations live in the Track A main tutorial:

- [Rigid-ligand docking](../tutorial.md#15-rigid-ligand-docking) introduces
  receptor-ligand interaction scoring, translation/rotation sampling, and the
  distinction between an interaction energy and a binding free energy.
- [Flexible-receptor docking](../tutorial.md#24-flexible-receptor-docking) adds
  explicit side-chain displacement and a harmonic strain penalty.

The supporting physics is in:

- [Atomistic water, atom types, and exclusions](../tutorial.md#14-atomistic-water-atom-types-and-exclusions);
- [Ewald summation for periodic charges](../tutorial.md#16-ewald-summation-for-periodic-charges);
- [Trajectory analysis](../tutorial.md#11-trajectory-analysis).

The applied module code is in `src/moleculab/applications/docking.py` and
`src/moleculab/applications/flexible_docking.py`. Examples are
`examples/rigid_docking.py` and `examples/flexible_docking.py`.

### 2.2 Module 01: protein structure fundamentals

The first module should teach how to inspect a protein structure before scoring
anything:

- PDB files, chains, residues, and atom names;
- experimental versus predicted structures;
- missing atoms, alternate locations, and incomplete side chains;
- hydrogen atoms and protonation states;
- water, ions, cofactors, and ligand records;
- symmetry mates and biological assemblies;
- what a structure does not contain: missing loops, disorder, and uncertainty.

The acceptance criterion is geometric and representational: a reader should be
able to identify the chains, residues, ligands, waters, and missing regions in a
real structure, and explain which of them are safe to use in a calculation.

### 2.3 Module 02: binding sites and cavities

This module should treat the binding site as a geometric object:

- pocket detection and cavity volume;
- buried versus exposed surface;
- shape complementarity;
- conserved and structural waters;
- the difference between an ATP site, an allosteric site, and a surface groove.

### 2.4 Module 03: protein–ligand interaction geometry

This module should identify interactions from coordinates before assigning
energies:

- donor–acceptor distances and hydrogen-bond angles;
- salt-bridge distances and charge complementarity;
- hydrophobic contacts and aromatic stacking geometry;
- steric clashes;
- metal coordination.

A geometric contact says that an interaction is possible in that structure. It
does not say how much that contact contributes to binding.

### 2.5 Module 04: interaction fingerprints

A residue-level fingerprint should record which ligand atoms contact which
protein residues and with which interaction type. Over a trajectory, the same
representation supports:

- contact occupancy;
- hydrogen-bond persistence;
- residue interaction timelines;
- pose comparison;
- ligand-induced changes in the binding site.

The essential distinction is:

> A geometric interaction fingerprint is not an energy decomposition.

### 2.6 Module 05: docking foundations

This module is already covered by the Track A docking sections. It should be
read as a Track C topic:

- rigid-body search;
- scoring with the typed nonbonded model;
- translation and rotation sampling;
- flexible side chains and strain;
- why the best score is a hypothesis, not a measured affinity.

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

## 6. Acceptance criteria

A Track C module is complete when:

- the protein representation is explicit about missing and modelled atoms;
- binding-site geometry is described without pretending it is an energy;
- interaction types are assigned from clear geometric criteria;
- a contact fingerprint can be reproduced from a known structure;
- contact occupancy is checked against a manually constructed trajectory;
- docking results are reported as poses and scores, not as binding affinities;
- the link to the relevant Track A and Track B foundations is explicit.

**Previous:** [Track B — Chemical structure and properties](../chemistry/README.md).  
**Next:** [Capstone — imatinib and BCR-ABL1](../capstone_imatinib.md).
