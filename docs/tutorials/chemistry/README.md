# Track B — Chemical structure and properties

> **Question.** What is the molecule, and what state is it in?
>
> **Entry point.** This page.
>
> **Prerequisites.** Track A sections 3, 4b, and 7 for the first module;
> Track A sections 14 and 16 for protonation and electrostatics.
>
> **Companion tracks.** [Track A — molecular physics](../physics/tutorial.md) ·
> [Track C — structure-based drug discovery](../drug_discovery/README.md)

This track connects molecular structure and chemical properties to the physics
engine. It does not replace the physics sequence. Track A supplies the energy
model; Track B supplies the chemical object, its identity, and its state.

For the three-track overview, see the [tutorial map](../README.md). For
prerequisites and the full curriculum, see [the learning tracks](../curriculum.md).

## Modules

| Module | Focus | Prerequisite | Status |
| --- | --- | --- | --- |
| [01 — Chemical structure](01_chemical_structure.md) | functional groups, rings, stereochemistry, conformations | Track A 3, 4b, 7 | available |
| [02 — Molecular representation](02_molecular_representation.md) | graphs, SMILES, aromaticity, stereochemistry, 3D conformers | module 01 | available |
| [03 — Molecular properties](03_molecular_properties.md) | descriptors, fingerprints, Tanimoto, chemical space | module 02 | available |
| [04 — Protonation and tautomerism](04_protonation_and_tautomerism.md) | pKa, pH populations, microstates, tautomers | modules 02–03; Track A 14, 16 | available |
| [05 — Stereochemistry and molecular recognition](05_stereochemistry.md) | R/S, E/Z, enantiomers, diastereomers, meso compounds | modules 01–02 | available |
| [06 — Molecular polarity and electronic properties](06_molecular_polarity.md) | electronegativity, formal and partial charge, dipoles, Coulomb connection | modules 02 and 04 | available |
| [07 — Hydrogen bonding and molecular recognition](07_hydrogen_bonding.md) | donors/acceptors, geometry, solvent competition, force-field connection | modules 04 and 06 | available |
| [08 — Medicinal chemistry and ADME fundamentals](08_medicinal_chemistry_adme.md) | solubility, permeability, LogP/LogD, metabolism, toxicity, and trade-offs | modules 03–07 | available |
| Solvation and LogP/LogD | transfer free energy and lipophilicity | module 04; Track A free energy | proposed |
| ADME-Tox | properties, metabolism, and empirical prediction | modules 03–04 | proposed |

## Reading order

- Start with [chemical structure](01_chemical_structure.md) after Track A
  section 4b or 7.
- Continue with [molecular representation](02_molecular_representation.md) and
  [molecular properties](03_molecular_properties.md).
- Read [protonation and tautomerism](04_protonation_and_tautomerism.md) after
  Track A electrostatics and water, or alongside them.
- Continue with [stereochemistry and molecular recognition](05_stereochemistry.md)
  after the representation module.
- Read [molecular polarity and electronic properties](06_molecular_polarity.md)
  alongside the electrostatics material in Track A.
- Continue with [hydrogen bonding and molecular recognition](07_hydrogen_bonding.md)
  before the interaction-analysis modules in Track C.
- Read [medicinal chemistry and ADME fundamentals](08_medicinal_chemistry_adme.md)
  after the protonation, polarity, and hydrogen-bonding modules.
- Then continue to [Track C](../drug_discovery/README.md) and the
  [capstone](../capstone_imatinib.md).
- Solvation, LogP/LogD, and ADME-Tox are proposed modules, not implemented
  tutorials.

## What this track adds to the engine

Track B answers questions the force field cannot answer by itself:

- Which atoms and bonds are present?
- Which ring is aromatic?
- Which stereoisomer is being modelled?
- Which protonation and tautomer states are populated at a given pH?
- Which descriptors and similarity measures summarize the molecule?

The force field then receives atom types, partial charges, bonded parameters,
and topology. The chemistry toolkit does the chemical perception; the physics
engine does the energy and dynamics.

## Design rule

Keep chemical perception, physical modelling, numerical simulation, and
drug-discovery interpretation separate. A toolkit may identify an amide or an
aromatic ring; the force field computes its energy; MD samples its motion; an
interaction analyser interprets its contacts. Those are related operations, not
one monolithic class.

RDKit is an optional external toolkit for chemical perception. Moleculab's
physics layer must not depend on RDKit types.

**Previous:** [Track A — molecular physics and simulation](../physics/tutorial.md).

**Next:** [Track C — structure-based drug discovery](../drug_discovery/README.md).
