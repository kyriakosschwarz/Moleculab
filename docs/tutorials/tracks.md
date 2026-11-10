# Curriculum and learning tracks

Moleculab is organized as three connected tracks rather than one ever-growing
tutorial. The [tutorial map](README.md) is the quick entry point; this page is
the detailed curriculum, prerequisite map, and status table.

```text
Track A — Molecular physics and simulation
    units → particles → forces → integration → molecular mechanics
          → environment → electrostatics → free energy

Track B — Chemical structure and properties
    structure → representation → properties → protonation states

Track C — Structure-based drug discovery
    protein structure → binding site → interactions → fingerprints → docking

                    \ | /
                     capstone
             imatinib and BCR-ABL1
```

The tracks share the same molecule and the same physical engine. They are
separated by the question being asked, not by a software boundary between
unrelated subjects.

## 1. Choosing a path

| Background or goal | Recommended path |
| --- | --- |
| new to molecular simulation | Track A sections 1–10, then Track B 01–04, then Track A 11–24, then Track C |
| chemistry-first | Track B 01–04, then Track A 1–10 and 14–17, then Track C |
| drug-discovery focus | Track A 1–10 and 14–17, Track B 01–04, Track C 01–05, then the capstone |
| software and architecture | Track A 1–10, then [the architecture](../architecture.md), then selected application sections |
| just want the integration story | [Introduction](introduction.md), then the [capstone](capstone_imatinib.md) |

The default order is not arbitrary. Track A sections 1–10 establish the energy,
force, integration, and molecular-mechanics vocabulary. Track B supplies the
chemical identity and state of the molecule. Track A sections 11–24 add the
environment and thermodynamics. Track C then asks recognition questions with
those tools.

## 2. Track A — Molecular physics and simulation

**Question:** how does the physical model work?

**Entry point:** [From one particle to a small simulation](tutorial.md)

**Slower visual path:** [Physics deep dives 1–10](physics/README.md)

**Prerequisites:** basic vectors, calculus, and Python. No prior MD experience
is assumed.

| Stage | Material | Covers | Status |
| --- | --- | --- | --- |
| Foundations | [tutorial sections 1–4](tutorial.md), [physics 01–04](physics/README.md), [checkpoint](physics/checkpoint_04.md) | units, particle state, energy and force, Velocity Verlet, analytic validation | available |
| Chemical bridge | [tutorial section 4b](physics/04b_chemical_bridge.md) | Lennard-Jones dissociation, torsional rotamers, double wells | available |
| Molecular mechanics | [tutorial sections 5–10](tutorial.md), [physics 05–10](physics/README.md) | nonbonded interactions, `System`, bonds/angles/torsions, minimisation, thermostats, periodic boundaries | available |
| Environment and analysis | [tutorial sections 11–14 and 16–17](tutorial.md) | MSD, RDF, solvent, atomistic water, Ewald, pressure/NPT | available |
| Free energy and alchemy | [tutorial sections 18–23](tutorial.md) | FEP, TI, MBAR, soft-core LJ, Monte Carlo alchemy, solvated-ligand decoupling | available |
| Track C bridge | [tutorial sections 15 and 24](tutorial.md) | rigid-ligand and flexible-receptor docking | available |
| Performance | [tutorial section 25](tutorial.md) | CPU default and optional PyTorch Ewald backend | available |

Track A is the dependency chain for the other tracks. The physics deep dives
expand the first ten concepts into a slower, figure-rich path. The main
tutorial keeps the same section numbering for stable links.

## 3. Track B — Chemical structure and properties

**Question:** what is the molecule, and what state is it in?

**Entry point:** [Chemical structure and properties](chemistry/README.md)

**Prerequisites:** Track A sections 3, 4b, and 7 are useful for the first
module; Track A electrostatics and water are useful for protonation.

| Module | Focus | Prerequisite | Status |
| --- | --- | --- | --- |
| [01 — Chemical structure](chemistry/01_chemical_structure.md) | functional groups, rings, stereochemistry, conformations | Track A 3, 4b, 7 | available |
| [02 — Molecular representation](chemistry/02_molecular_representation.md) | graphs, SMILES, aromaticity, stereochemistry, 3D conformers | module 01 | available |
| [03 — Molecular properties](chemistry/03_molecular_properties.md) | descriptors, fingerprints, Tanimoto, chemical space | module 02 | available |
| [04 — Protonation and tautomerism](chemistry/04_protonation_and_tautomerism.md) | pKa, pH populations, microstates, tautomers | modules 02–03; Track A 14, 16 | available |
| [05 — Stereochemistry and molecular recognition](chemistry/05_stereochemistry.md) | R/S, E/Z, enantiomers, diastereomers, meso compounds | modules 01–02 | available |
| Solvation and LogP/LogD | transfer free energy, neutral versus ionized distribution | module 04; Track A free energy | proposed |
| ADME-Tox | lipophilicity, solubility, permeability, metabolic soft spots | modules 03–04 | proposed |

Track B uses RDKit as an optional external toolkit for chemical perception.
Moleculab's force field receives atom types, charges, and topology; it does not
import RDKit types.

## 4. Track C — Structure-based drug discovery

**Question:** what does the physics mean for a drug?

**Entry point:** [Structure-based drug discovery](drug_discovery/README.md)

**Prerequisites:** Track A sections 1–14 and Track B modules 01–04. Docking
also uses Track A sections 15–17 and the free-energy material.

| Module | Focus | Foundation | Status |
| --- | --- | --- | --- |
| [01 — Protein structure fundamentals](drug_discovery/01_protein_structure.md) | chains, residues, coordinates, missing atoms, structure quality | Track A particles and topology | available |
| [02 — Binding sites and cavities](drug_discovery/02_binding_sites.md) | pocket geometry, waters, volume, complementarity | Track A periodic box and nonbonded model | available |
| [03 — Protein–ligand interaction geometry](drug_discovery/03_protein_ligand_interactions.md) | hydrogen bonds, salt bridges, hydrophobic and aromatic contacts | Track A electrostatics and water | available |
| [04 — Interaction fingerprints](drug_discovery/04_interaction_fingerprints.md) | residue-ligand contacts and sampled-frame occupancy | Track A trajectory analysis | available |
| [05 — Docking foundations](drug_discovery/05_docking_fundamentals.md) | search, scoring, poses, limitations | Track A sections 15 and 24 | available |
| [06 — Comparing binding modes](drug_discovery/06_comparing_binding_modes.md) | pose clustering, geometry, and interaction patterns | Track C modules 04–05 | available |
| [07 — Docking validation and benchmarking](drug_discovery/07_docking_validation.md) | redocking, cross-docking, pose recovery, score interpretation | Track C modules 04–06 | available |
| [08 — Trajectory contact analysis](drug_discovery/08_trajectory_contact_analysis.md) | contact matrices, occupancy, and lifetimes | Track A trajectory analysis | available |
| [09 — Lead optimisation and SAR](drug_discovery/09_lead_optimisation_sar.md) | matched pairs, activity scales, ligand efficiency, and evidence categories | Track B and Track C | available |
| [10 — Structural observations to thermodynamics](drug_discovery/10_structural_observations_to_thermodynamics.md) | occupancy, ensembles, free energy, and contact interpretation | Track A free energy | available |

The current library supports the geometric and scoring pieces of Track C. It
does not provide protein preparation, a trained affinity model, or a validated
prospective docking pipeline.

## 5. Where the tracks meet

The [imatinib–BCR-ABL1 capstone](capstone_imatinib.md) is the integration
exercise. It uses all three tracks:

| Capstone question | Track A | Track B | Track C |
| --- | --- | --- | --- |
| What is imatinib? | molecular mechanics | identity, tautomer, protonation | ligand representation |
| Where does it bind? | periodic box, water, electrostatics | chemical state | protein structure, binding site |
| Which contacts matter? | nonbonded interactions | hydrogen-bond donors and acceptors | interaction geometry, fingerprints |
| How stable is the complex? | MD, free energy | protonation-coupled states | pose comparison, occupancy |
| Why do mutations resist? | force-field changes | chemical-state changes | binding-site geometry changes |

The capstone is not a turnkey pipeline. It is a guided map from the tutorial
components to a realistic problem.

## 6. Cross-cutting corrections

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
- **A geometric interaction is not an energy contribution.** A contact table
  says what is near what. It does not prove what dominates the binding free
  energy.

## 7. Recommended study and development order

For a reader:

1. Track A foundations and molecular mechanics.
2. Track B chemical structure, representation, properties, and protonation.
3. Track A environment, electrostatics, and free energy.
4. Track C protein structure, binding sites, interaction geometry, and
   fingerprints.
5. The capstone.

For a module author:

1. Define the question and the track it belongs to.
2. State prerequisites and the model boundary.
3. Derive the equation or the geometric criterion.
4. Implement only the part that belongs to that layer.
5. Add an independent test and a runnable example.
6. Add the module to the track table and the tutorial map.
7. Check links and math rendering on GitHub and in VS Code.

The chemistry layer should remain separate from the physics engine: chemical
perception identifies an amide or an aromatic ring, the force field computes its
energy, MD samples its motion, and an interaction analyser interprets its
contacts. Those are related operations, not one monolithic class.
