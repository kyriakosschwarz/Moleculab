# Capstone: imatinib and BCR-ABL1

> **Where this fits.** This is the integration exercise at the end of all three
> tracks. It does not introduce a new module. It asks you to connect the
> physics, chemical state, protein structure, docking, dynamics, and
> free-energy ideas to one real drug-discovery story.
>
> **What you need.** The core Track A tutorials, especially 05–10 and the
> sections on solvent, water, docking, Ewald, and free energy; Track B modules
> 01–04 for chemical identity and protonation; and Track C modules 01–10 for
> protein structure, binding sites, interactions, docking, validation,
> trajectory contacts, SAR, and the structure–thermodynamics bridge.

The capstone is the point where the tracks meet:

| Track | Contribution to the case study |
| --- | --- |
| [A — molecular physics](tutorial.md) | force field, water, MD, docking score, free energy |
| [B — chemical structure and properties](chemistry/README.md) | molecular identity, tautomer and protonation state, descriptors |
| [C — structure-based drug discovery](drug_discovery/README.md) | protein structure, binding site, interactions, pose interpretation |

For the three-track overview, see the [tutorial map](README.md); for the full
curriculum, see [the learning tracks](tracks.md).

## 1. Why imatinib?

Imatinib (Gleevec/Glivec) is a small-molecule tyrosine kinase inhibitor used to
treat chronic myeloid leukemia (CML). It is an excellent teaching case because:

- the disease has a clear molecular cause;
- the drug target and binding site are structurally characterized;
- the mechanism is understood at the molecular level;
- experimental protein-ligand structures are public;
- resistance mutations connect structure, binding, and clinical outcome;
- the example touches nearly every layer of a CADD workflow.

The capstone is not a claim that Moleculab is a validated imatinib design tool.
It is a guided map from the tutorial foundations to a realistic problem.

## 2. Disease and target

CML is associated with the Philadelphia chromosome, produced by a reciprocal
translocation between chromosomes 9 and 22. The translocation creates the
BCR-ABL1 fusion gene.

The ABL1 side contributes a tyrosine kinase domain. In the fusion protein, the
kinase activity is abnormally persistent, driving uncontrolled proliferation.

The relevant drug-binding region is the ATP-binding pocket of the ABL1 kinase
domain. The full BCR-ABL1 fusion protein is larger and more complicated than
the kinase domain alone, so structural work usually starts from the kinase
domain and treats it as the relevant model system.

## 3. What imatinib does

Imatinib is an ATP-competitive inhibitor. It occupies the ATP-binding region and
stabilizes an inactive kinase conformation. That reduces the kinase activity of
the fusion protein.

Two distinctions matter:

- imatinib does not simply destroy the protein; it changes its catalytic
  behaviour by binding a particular conformation;
- inhibiting the fusion protein is not the same as being absolutely selective
  for the fusion over normal ABL1.

The structural question is therefore not only “does it fit?” but “which
conformation does it stabilize, and how does that change the free energy?”

## 4. Structural starting points

Useful public structures include:

| Structure | Contents | Caveat |
| --- | --- | --- |
| PDB 2HYY | human ABL1 kinase domain with imatinib | kinase domain, not the full BCR-ABL1 fusion |
| PDB 1IEP | mouse c-Abl with imatinib | mouse proxy; used in the classic Vina tutorial |

The PDB chemical component identifier for imatinib is `STI`.

Before simulating anything:

1. inspect the chains, ligands, waters, ions, and missing atoms;
2. decide which chain and which ligand copy match each other;
3. check protonation and hydrogen placement;
4. keep the experimental pose separate from an ideal ligand conformer.

Production preparation is a project in itself. Moleculab provides the physical
concepts and some of the scoring machinery, not a complete protein-preparation
pipeline.

## 5. Map the case study onto the tutorials

| Case-study concept | Tutorial foundation |
| --- | --- |
| atoms, positions, charges | `Particle`, units |
| bond, angle, torsion terms | `Topology`, bonded potentials |
| Lennard-Jones and electrostatics | nonbonded interactions, `NonbondedModel` |
| protein-ligand complex | `System` composition |
| energy minimisation | steepest descent |
| MD equilibration | Velocity Verlet, thermostats, barostats |
| periodic environment | `Box`, minimum image |
| solvent and ions | water and solvent models |
| long-range electrostatics | reaction field, Ewald |
| binding pose | rigid and flexible docking |
| binding affinity | FEP, TI, MBAR, alchemical decoupling |
| electronic effects | QM/MM concepts |

## 6. A realistic workflow

The following is a conceptual workflow. External tools are needed for
protein-scale preparation and production MD; the Moleculab tutorials explain the
underlying models.

### Step 1 — inspect the structure

Use a molecular viewer to identify:

- the kinase domain;
- imatinib in the ATP pocket;
- the hinge region and gatekeeper residue;
- the DFG motif and activation loop;
- crystallographic waters and ions;
- missing loops or atoms.

Do not skip this step. A simulation of a badly prepared structure is a
simulation of the wrong system.

### Step 2 — prepare the model

Protein and ligand parameters usually come from different force fields. A
typical production workflow combines:

- a protein force field such as AMBER14SB;
- a small-molecule force field such as OpenFF Sage;
- a water model such as TIP3P or SPC/E;
- counterions and a periodic box.

Moleculab's `NonbondedModel`, `Topology`, `applications/water.py`, and `Box` show how those
pieces enter the energy function.

### Step 3 — minimise and equilibrate

Start from a prepared structure and remove severe strain:

1. minimise the solvated complex;
2. heat gradually;
3. equilibrate density and temperature;
4. check energy, temperature, pressure, and structural stability.

The minimisation, thermostat, and barostat tutorials explain the algorithms.
Production MD adds constraints, multiple stages, and much longer timescales.

### Step 4 — redock and compare with the crystal pose

A useful validation exercise is **redocking**:

1. remove imatinib from the prepared receptor;
2. generate poses with a docking program;
3. compare the best pose with the experimental pose using RMSD;
4. ask whether the scoring function ranks the crystal-like pose well.

A pose RMSD around 2 Å is often used as a practical redocking benchmark, but
symmetry, atom mapping, and alignment choices matter.

Moleculab's `applications/docking.py` and `applications/flexible_docking.py` show how interaction-energy
scoring is constructed. They do not replace a production docking engine.

### Step 5 — run MD and analyse

A protein-ligand MD trajectory can be used to study:

- ligand stability in the pocket;
- hydrogen-bond persistence;
- protein flexibility;
- water occupancy;
- conformational changes;
- RMSD and RMSF;
- interactions with resistance mutations.

Moleculab's `analysis/trajectory.py` provides MSD and RDF; RMSD/RMSF and hydrogen-bond
analysis are natural extensions.

### Step 6 — estimate a binding free energy

Docking gives a pose and an interaction score. A binding free energy requires
statistical mechanics:

- FEP for an alchemical transformation;
- TI for a thermodynamic path;
- MBAR for combining multiple sampled states;
- careful treatment of ligand-water and protein-water interactions.

The free-energy tutorials show the estimators on small, analytically controlled
systems. A production calculation needs many windows, replicates, and a
convergence analysis.

### Step 7 — include electronic effects where needed

Classical force fields use fixed charges and bonded parameters. For chemistry
that depends on electronic structure, a QM/MM treatment may be needed.

Moleculab's `applications/qm_mm.py` is a one-dimensional finite-difference model, not a
production QM/MM engine. It teaches the coupling idea: a quantum region
embedded in a classical environment.

## 7. Water is part of the mechanism

Water is not background. In the imatinib–ABL1 story, water matters because:

- the binding pocket is solvated before the ligand arrives;
- some waters are displaced on binding;
- water can mediate hydrogen bonds between ligand and protein;
- water screens electrostatic interactions;
- released water changes the entropy and enthalpy of binding;
- water occupancy can change the preferred pose or the apparent affinity.

This is why a gas-phase interaction energy is not a binding free energy. The
same lesson appears in the docking tutorial: a score is not a measured
affinity.

Moleculab represents water in several ways:

| Model | Role |
| --- | --- |
| `applications/solvent.py` | coarse-grained solvent box |
| `applications/water.py` | three-site atomistic water with constraints |
| `models/ewald.py` | long-range periodic electrostatics |
| `applications/solvated_alchemy.py` | ligand-water decoupling |
| `applications/fep.py`, `applications/ti.py`, `applications/mbar.py` | free-energy estimators |

The library's default water model is SPC/E. TIP3P is another common choice.
Both are models; neither is “the true water.”

## 8. Resistance mutations

The T315I gatekeeper mutation is a famous example of resistance to imatinib.
It illustrates the connection between structure, binding, and free energy:

- a mutation can change the shape or chemistry of the pocket;
- that can weaken ligand binding;
- a weaker binding free energy can reduce clinical effectiveness;
- a different inhibitor may be needed.

The educational point is not to memorize one mutation. It is to see why a
binding free energy, not a single docking score, is the relevant quantity.

## 9. Exercises

### Structural reasoning

1. Identify the ATP-binding pocket and the hinge region in PDB 2HYY.
2. List the crystallographic waters near imatinib. Which could mediate
   hydrogen bonds?
3. Compare the experimental imatinib pose with an ideal ligand conformer.

### Docking

1. Redock imatinib and compute pose RMSD to the crystal structure.
2. Compare rigid and flexible-receptor docking.
3. Explain why a good docking score does not prove a strong binding affinity.

### Dynamics

1. Run a small solvated system and track a hydrogen bond over time.
2. Compute RDF for water around a charged group.
3. Explain how water displacement could change the binding entropy.

### Free energy

1. Explain the difference between an interaction energy and a binding free
   energy.
2. Describe an alchemical path for decoupling imatinib from the pocket.
3. Explain why multiple windows and replicates are needed.

### Water

1. Why can a binding pocket water be favourable before binding and unfavourable
   after ligand displacement?
2. How does electrostatic screening change the apparent interaction strength?
3. Why does a vacuum docking score overstate electrostatic interactions?

## 10. What Moleculab can and cannot do here

Moleculab can:

- explain and demonstrate the physical models;
- build small solvated systems;
- run MD, minimisation, docking, and free-energy estimators on controlled
  examples;
- show how water, Ewald, and alchemical decoupling enter the calculation.

Moleculab cannot yet:

- prepare a full protein-ligand system automatically;
- assign production protein and ligand force-field parameters;
- run a publication-scale imatinib MD campaign;
- provide a validated imatinib binding free energy;
- replace OpenMM, GROMACS, AutoDock Vina, OpenFF, or ORCA.

The capstone is successful when you can explain **why** each production step
exists and **which physical model** it is using.

## 11. Further tools

To go beyond the educational engine:

- **PDB / RCSB** for experimental structures;
- **OpenMM** and **GROMACS** for production MD;
- **OpenFF** and **AMBER** for force-field parameters;
- **AutoDock Vina** for docking;
- **ChimeraX** or **PyMOL** for structure inspection;
- **ORCA**, **Psi4**, or **PySCF** for quantum chemistry.

The goal is not to replace these tools. The goal is to understand what they are
doing.

**Related:** [tutorial map](README.md), [introduction](introduction.md),
[Track A main path](tutorial.md), [Track A physics deep dives](physics/README.md),
[Track B](chemistry/README.md), [Track C](drug_discovery/README.md),
[learning tracks](tracks.md), [architecture](../architecture.md), and
[validation](../validation.md).
