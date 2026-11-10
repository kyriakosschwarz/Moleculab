# Learning tracks

Moleculab's material is easier to follow when it is split into three connected
tracks rather than one ever-growing tutorial.

| Track | Question | Core material |
| --- | --- | --- |
| **A — Molecular physics and simulation** | How does the physical model work? | units, particles, potentials, integrators, molecular mechanics, MD, thermodynamics, free energy, QM/MM |
| **B — Chemical structure and properties** | What is the molecule? | functional groups, rings, stereochemistry, conformations, protonation, tautomers, descriptors, similarity, solvation, ADME |
| **C — Structure-based drug discovery** | What does the physics mean for a drug? | protein structures, binding sites, docking, pharmacophores, interaction fingerprints, trajectory analysis, lead optimization |

The tracks are not isolated. They share the same molecule and the same physical
engine. The [imatinib–BCR-ABL1 capstone](capstone_imatinib.md) is the point
where all three meet.

## Prerequisites and entry points

| Material | Prerequisite | Outcome |
| --- | --- | --- |
| Physics parts 01–10 | basic vectors and Python | build and run a molecular-mechanics system |
| [Chemical structure](chemistry/01_chemical_structure.md) | parts 03, 04b, 07 | reason about functional groups, rings, stereochemistry, and conformations |
| Cheminformatics module | chemical structure | descriptors, fingerprints, similarity |
| pKa and protonation module | electrostatics and water | choose chemically realistic protonation states |
| Structural biology module | docking and proteins | judge experimental versus predicted structures |
| Interaction fingerprints | docking | residue-level geometric contacts and occupancy |
| Solvation, LogP/LogD | free energy | connect solvation thermodynamics to lipophilicity |
| ADME-Tox | properties and similarity | separate physical simulation from empirical prediction |

## Track A — molecular physics and simulation

This is the existing sequence:

1. units and scales;
2. particle state and kinetic energy;
3. potential energy and force;
4. Velocity Verlet;
5. checkpoint;
6. chemical bridge: Lennard-Jones, rotamers, double wells;
7. nonbonded interactions;
8. many-particle systems;
9. bonds, angles, torsions;
10. minimisation;
11. temperature and thermostats;
12. periodic boundaries.

The main tutorial then continues through trajectory analysis, solvent, water,
docking, Ewald, pressure, free energy, FEP/TI/MBAR, QM/MM, and flexible docking.

## Track B — chemical structure and properties

The first companion module is
[Chemical structure and conformational behaviour](chemistry/01_chemical_structure.md).
It covers:

- functional groups and their structural consequences;
- rings, aromaticity, and restricted rotation;
- stereochemistry;
- ethane and butane rotamers;
- ring strain;
- the difference between a force-field torsion term and the electronic effects
  it represents.

Later modules in this track are proposed, not yet implemented:

| Module | Focus | Status |
| --- | --- | --- |
| Cheminformatics | representations, descriptors, fingerprints, Tanimoto similarity | proposed |
| pKa, protonation, tautomerism | discrete chemical states and pH populations | proposed |
| Solvation and LogP/LogD | transfer free energy, neutral versus ionized distribution | proposed |
| ADME-Tox | lipophilicity, solubility, permeability, metabolic soft spots | proposed |

These modules should use established chemistry toolkits such as RDKit where
appropriate rather than reimplementing cheminformatics algorithms.

## Track C — structure-based drug discovery

This track starts from the docking material and adds:

- protein structure quality and preparation;
- binding-site inspection;
- geometric interaction fingerprints;
- hydrogen-bond, salt-bridge, hydrophobic, and aromatic contacts;
- contact occupancy over an MD trajectory;
- pharmacophore hypotheses and pose validation.

Two distinctions are essential:

1. **Geometric fingerprints are not energy decomposition.** A contact table says
   which residues are near which ligand groups and for how long. It does not
   prove that those residues dominate the binding free energy.
2. **Per-residue energy contributions are model-dependent.** Solvent
   reorganisation, entropy, polarisation, and cooperativity mean the total
   binding free energy is not a simple sum of independent residue terms.

## Cross-cutting corrections

Several common simplifications are deliberately avoided:

- **pH does not make partial charges vary continuously.** It changes the
  population of discrete protonation states, each with its own topology and
  charge distribution.
- **LogP is not LogD.** LogP usually refers to the neutral compound; LogD is the
  distribution coefficient at a specified pH and includes ionized states.
- **PAINS are not toxicity filters.** They flag assay interference and
  potentially misleading screening results.
- **AlphaFold is not an experimental ligand-bound structure.** A confident
  prediction does not guarantee the correct binding-site conformation or
  induced fit.
- **A hydrogen bond is not always a separate force-field term.** In a
  fixed-charge force field it usually emerges from electrostatics and van der
  Waals interactions.
- **Water is part of the thermodynamics.** Desolvation and water networks can
  decide whether an apparently favourable contact actually helps binding.

## Recommended implementation order

1. Chemical structure and conformational behaviour (this track's first module).
2. Cheminformatics representations and similarity.
3. pKa, protonation, and tautomerism.
4. Protein–ligand geometric interaction fingerprints.
5. Solvation, desolvation, LogP/LogD, and free energy.

The chemistry layer should remain separate from the physics engine: chemical
perception identifies an amide or an aromatic ring, the force field computes its
energy, MD samples its motion, and an interaction analyser interprets its
contacts. Those are related operations, not one monolithic class.
