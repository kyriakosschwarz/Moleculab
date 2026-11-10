# Chemistry track

This companion track connects molecular structure and chemical properties to
the physics engine. It does not replace the physics sequence.

## Modules

| Module | Focus | Status |
| --- | --- | --- |
| [01 — Chemical structure](01_chemical_structure.md) | functional groups, rings, stereochemistry, conformations | available |
| Cheminformatics | representations, descriptors, fingerprints, similarity | proposed |
| pKa, protonation, tautomerism | discrete chemical states and pH populations | proposed |
| Solvation and LogP/LogD | transfer free energy and lipophilicity | proposed |
| ADME-Tox | properties, metabolism, and empirical prediction | proposed |

## Reading order

The chemistry track can be read alongside the physics track:

- after part 04b, read [chemical structure](01_chemical_structure.md);
- after the electrostatics and water material, continue to pKa and solvation;
- after docking, continue to interaction fingerprints and pharmacophores.

See [the three learning tracks](../tracks.md) for prerequisites and the full
roadmap.

## Design rule

Keep chemical perception, physical modelling, numerical simulation, and
drug-discovery interpretation separate. A toolkit may identify an amide or an
aromatic ring; the force field computes its energy; MD samples its motion; an
interaction analyser interprets its contacts. Those are related operations, not
one monolithic class.
