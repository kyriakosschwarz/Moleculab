# Tutorials: three connected tracks

Moleculab is learned along three connected tracks. They share one molecule and
one physics engine, but they answer different questions.

```text
Track A — Molecular physics and simulation
    How does the physical model work?

Track B — Chemical structure and properties
    What is the molecule, and what state is it in?

Track C — Structure-based drug discovery
    What does the physics mean for a drug?

                    \ | /
                     capstone
             imatinib and BCR-ABL1
```

The tracks are not three separate courses with no overlap. Track A supplies the
engine. Track B supplies the chemical object and its state. Track C turns the
engine and the chemical object into a molecular-recognition question. The
[imatinib–BCR-ABL1 capstone](capstone_imatinib.md) is where they meet.

For the detailed prerequisites and module-by-module curriculum, see
[the learning tracks](tracks.md). This page is the quick orientation.

## Start here

| If you want to... | Start with | Then |
| --- | --- | --- |
| understand the simulation engine | [Track A](tutorial.md) | [physics deep dives](physics/README.md), [Track C](drug_discovery/README.md) |
| understand chemical identity and states | [Track B](chemistry/README.md) | [Track A electrostatics and water](tutorial.md), [Track C](drug_discovery/README.md) |
| study molecular recognition | [Track A](tutorial.md) foundations, then [Track C](drug_discovery/README.md) | [Track B](chemistry/README.md), [capstone](capstone_imatinib.md) |
| see the whole picture | [introduction](introduction.md) | [tracks](tracks.md), [capstone](capstone_imatinib.md) |

A practical default is:

1. Track A, sections 1–10, with the [physics deep dives](physics/README.md)
   as needed.
2. Track B, modules 01–04.
3. Track A, sections 11–24, for environment, electrostatics, free energy,
   and the docking bridge.
4. Track C, modules 01–05.
5. The [capstone](capstone_imatinib.md).

## Track A — Molecular physics and simulation

**Question:** how does the physical model work?

**Entry point:** [From one particle to a small simulation](tutorial.md)

**Slower path:** [Physics deep dives 1–10](physics/README.md)

Track A builds units, particle state, energy and force, integration, molecular
mechanics, periodic boundaries, solvent, electrostatics, and free energy. It is
the dependency chain that the other tracks assume.

The main tutorial also contains two Track C bridge sections:

- [rigid-ligand docking](tutorial.md#15-rigid-ligand-docking);
- [flexible-receptor docking](tutorial.md#24-flexible-receptor-docking).

They use Track A physics, but they ask a Track C question.

## Track B — Chemical structure and properties

**Question:** what is the molecule, and what state is it in?

**Entry point:** [Chemical structure and properties](chemistry/README.md)

Track B covers functional groups, stereochemistry, conformations, molecular
representation, descriptors, similarity, protonation, and tautomerism. It uses
RDKit as an optional external toolkit for chemical perception; the physics
engine does not depend on it.

## Track C — Structure-based drug discovery

**Question:** what does the physics mean for a drug?

**Entry point:** [Structure-based drug discovery](drug_discovery/README.md)

Track C starts from protein structure and binding-site geometry, then builds
towards interaction fingerprints, docking, and trajectory contact analysis.
Modules [01](drug_discovery/01_protein_structure.md),
[02](drug_discovery/02_binding_sites.md),
[03](drug_discovery/03_protein_ligand_interactions.md),
[04](drug_discovery/04_interaction_fingerprints.md),
[05](drug_discovery/05_docking_fundamentals.md),
[06](drug_discovery/06_comparing_binding_modes.md),
[07](drug_discovery/07_docking_validation.md),
[08](drug_discovery/08_trajectory_contact_analysis.md),
[09](drug_discovery/09_lead_optimisation_sar.md), and
[10](drug_discovery/10_structural_observations_to_thermodynamics.md) are
available. The docking foundations are already available as Track A
sections 15 and 24. [Units and conventions](drug_discovery/units.md) collects
the common structural and affinity units.

## Integration and reference

| Material | Use it for |
| --- | --- |
| [Introduction](introduction.md) | why the tracks exist and how CADD methods connect |
| [Learning tracks](tracks.md) | prerequisites, module tables, and study orders |
| [Capstone: imatinib and BCR-ABL1](capstone_imatinib.md) | the cross-track integration exercise |
| [Track C units](drug_discovery/units.md) | Angstrom, nm, kcal/mol, affinity, and structural units |
| [Troubleshooting](troubleshooting.md) | symptom table and diagnostic recipes |
| [Visualization prompts](visualization_prompts.md) | notebook-sized plotting prompts |
| [Answer keys](answer_keys.md) | answer sketches and expected magnitudes |
| [Architecture](../architecture.md) | how the software is structured |
| [Validation](../validation.md) | how correctness is checked |

## A note on naming

The tutorial files predate the three-track organisation. Their titles and
section numbers are kept stable so existing links and citations continue to
work:

- `tutorial.md` is the Track A main path.
- `physics/` contains Track A deep dives for the first ten concepts; it was
  formerly named `individual/`.
- `chemistry/` is Track B.
- `drug_discovery/` is Track C.
- `capstone_imatinib.md` integrates all three tracks.
