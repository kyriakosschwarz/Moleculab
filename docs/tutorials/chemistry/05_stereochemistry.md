# Track B, module 05 — Stereochemistry and molecular recognition

> **Where this fits.** Track B modules 01–04 established structure,
> representation, properties, and protonation. This module adds the spatial
> configuration that formula, graph, and many descriptors do not capture.
>
> **What you should already know.** Molecular graphs, SMILES, conformers,
> formal charge, and the difference between a graph and a geometry.
>
> **What you will be able to do afterwards.** Assign R/S and E/Z
> stereochemistry, distinguish enantiomers from diastereomers, explain why
> stereochemistry affects recognition, and validate stereochemical
> representation with RDKit.

## 1. Learning objectives

**Core**

- Distinguish constitutional isomers, configurational stereoisomers, and
  conformers.
- Identify tetrahedral stereocentres and assign R/S configuration.
- Distinguish enantiomers from diastereomers.
- Assign E/Z configuration at double bonds.
- Explain why stereochemistry affects protein binding.

**Deeper insight**

- A SMILES `@`/`@@` marker is relative to atom order, not an absolute R/S
  label.
- E/Z is more general than cis/trans.
- Meso compounds contain stereocentres but are achiral overall.
- Unspecified stereochemistry is not the same as a racemate.
- Formula, connectivity, and many descriptors do not distinguish
  stereoisomers.

**Looking ahead**

- The capstone requires the correct stereoisomer, not just the correct graph.
- Chiral ligands can have different binding modes in a chiral pocket.
- Docking and free-energy calculations inherit the stereochemical assumptions
  of their input.

## 2. Why stereochemistry matters

Two molecules can have the same:

- molecular formula;
- connectivity;
- functional groups;
- molecular weight;
- many calculated descriptors;

and still differ in three-dimensional configuration.

That difference can change:

- whether a protein recognises the molecule;
- which subpocket is occupied;
- hydrogen-bond geometry;
- the accessible conformations;
- metabolic handling;
- toxicity.

Stereochemistry is therefore part of chemical identity, not a cosmetic detail
of the drawing.

## 3. Configuration versus conformation

A **conformer** is a different spatial arrangement reached by rotation about
single bonds or other low-barrier motions.

A **configurational stereoisomer** requires bond breaking or a high-barrier
process to interconvert.

| Change | Example | Interconversion |
| --- | --- | --- |
| conformation | butane anti ↔ gauche | rotation about a single bond |
| configuration | R ↔ S at a stereocentre | bond breaking or inversion |
| configuration | E ↔ Z at a double bond | rotation about a double bond |

At ordinary temperatures, butane conformers interconvert rapidly. Lactic-acid
enantiomers do not interconvert without breaking bonds.

The example demonstrates butane conformers with the same graph:

```text
seed 7:  dihedral C-C-C-C = 180.00 degrees
seed 21: dihedral C-C-C-C =  60.00 degrees
```

Both are conformers of the same molecule, not separate stereoisomers.

## 4. Chirality and R/S configuration

A tetrahedral carbon with four different substituents is a stereocentre.

The R/S label follows the Cahn–Ingold–Prelog priority rules:

1. Rank the four substituents by atomic number, with higher atomic number
   first.
2. If necessary, compare the next atoms along each substituent.
3. Treat double bonds as duplicated atoms for priority.
4. Orient the lowest-priority substituent away from the viewer.
5. Clockwise is R; counterclockwise is S.

The label is a property of the configuration, not of the drawing orientation.

### 4.1 Lactic acid

Lactic acid has one stereocentre:

```text
R: C[C@@H](O)C(=O)O    centers [(1, 'R')]
S: C[C@H](O)C(=O)O     centers [(1, 'S')]
```

Both have formula $\mathrm{C_3H_6O_3}$ and the same connectivity. They are
enantiomers.

The `@`/`@@` markers in SMILES are not absolute R/S labels. They describe the
local order around the atom and can change if the atom order in the SMILES is
changed. Use a stereochemistry assignment step, not the raw `@`/`@@` text, to
obtain CIP labels.

## 5. Enantiomers, diastereomers, and constitutional isomers

| Relationship | Definition | Example |
| --- | --- | --- |
| constitutional isomers | different connectivity | ethanol and dimethyl ether |
| enantiomers | non-superimposable mirror images | R- and S-lactic acid |
| diastereomers | stereoisomers that are not mirror images | R,R- and meso-tartaric acid |

Enantiomers have identical properties in an achiral environment. They differ in
optical rotation and in interactions with chiral environments such as proteins.

Diastereomers generally differ in physical properties and can be separated by
ordinary achiral methods.

### 5.1 Tartaric acid

Tartaric acid has two stereocentres. The example shows:

```text
R,R:   centers [(1, 'R'), (5, 'R')]
S,S:   centers [(1, 'S'), (5, 'S')]
meso:  centers [(1, 'R'), (5, 'S')]
```

R,R and S,S are enantiomers. The meso R,S form is a diastereomer and has an
internal mirror plane, so the molecule is achiral overall despite having
stereocentres.

## 6. E/Z double-bond stereochemistry

Double bonds can also be stereogenic when each carbon bears two different
substituents.

Use CIP priorities on each carbon:

- if the higher-priority substituents are on the same side, the configuration
  is Z;
- if they are on opposite sides, the configuration is E.

For simple two-substituent alkenes, E is often called trans and Z cis, but the
E/Z notation is more general.

The example shows:

```text
E: C/C=C/C    double-bond stereo [(1, 'Bond_Trans')]
Z: C/C=C\C    double-bond stereo [(1, 'Bond_Cis')]
```

E and Z are configurational isomers. They do not interconvert by ordinary bond
rotation without breaking the double bond.

## 7. Meso compounds and symmetry

A molecule can contain stereocentres and still be achiral if it has an internal
plane or centre of symmetry.

Tartaric acid's R,S form is the standard example. The two stereocentres have
opposite configurations, and the molecule is superimposable on its mirror
image.

This is why counting stereocentres is not the same as determining chirality.

## 8. Stereochemistry and protein binding

Proteins are chiral environments. A binding pocket is built from amino acids
with defined stereochemistry, so enantiomers can interact differently.

Possible consequences:

- one enantiomer binds productively; the other does not;
- both bind, but with different affinities;
- one is a substrate and the other an inhibitor;
- one has the desired activity and the other has an off-target effect;
- the stereoisomers have different metabolic or toxicity profiles.

Thalidomide is a well-known cautionary example: the R and S enantiomers differ
in biological effects, and the drug can racemise in vivo. The exact historical
and pharmacological details are more complex than a single label, so use the
example to illustrate the principle rather than as a complete case study.

A docking or MD workflow must therefore use the intended stereoisomer. A
correct molecular graph with the wrong configuration is still the wrong
molecule.

## 9. RDKit representation and validation

RDKit can assign and validate stereochemistry. The runnable example is:

```bash
uv run --with rdkit python examples/stereochemistry.py
```

### 9.1 Reading stereochemistry

```python
from rdkit import Chem
from rdkit.Chem import FindMolChiralCenters

molecule = Chem.MolFromSmiles("C[C@H](O)C(=O)O")
Chem.AssignStereochemistry(molecule, cleanIt=True, force=True)
print(FindMolChiralCenters(molecule, includeUnassigned=True))
print(Chem.MolToSmiles(molecule))
print(Chem.MolToInchiKey(molecule))
```

Observed for the S enantiomer:

```text
[(1, 'S')]
C[C@H](O)C(=O)O
JVTAAEKCZFNVCJ-REOHCLBHSA-N
```

### 9.2 Reading E/Z stereochemistry

```python
for info in Chem.FindPotentialStereo(molecule):
    if info.type.name == "Bond_Double":
        print(info.centeredOn, info.descriptor.name)
```

For E-2-butene this reports a `Bond_Trans` descriptor. For Z-2-butene it
reports `Bond_Cis`.

### 9.3 Canonical SMILES and InChIKey

Canonical SMILES preserves stereochemistry, and InChIKey distinguishes
stereoisomers. For the lactic-acid enantiomers:

```text
R: JVTAAEKCZFNVCJ-UWTATZPHSA-N
S: JVTAAEKCZFNVCJ-REOHCLBHSA-N
```

The first block is shared; the stereochemical layer differs.

Do not infer identity from the molecular formula or a 2D descriptor table
alone.

### 9.4 Unspecified stereochemistry

A molecule with no stereochemical information is not necessarily a racemate. It
may be:

- unknown;
- deliberately unspecified;
- a mixture;
- a single enantiomer whose configuration was not recorded.

RDKit reports an unassigned centre as `?` when `includeUnassigned=True` is used:

```text
unspecified stereocentre: [(1, '?')]
```

Unspecified is not the same as racemic.

### 9.5 Validation checks

The example also checks the following:

1. CIP labels are meaningful only after stereochemistry has been specified or
   assigned.
2. A stereochemical SMILES round trip preserves the intended configuration.
3. Enantiomers can share ordinary achiral descriptors such as molecular weight
   and TPSA but still behave differently in a chiral environment.
4. Changing the conformer-generation seed changes the torsion without changing
   configurational identity.

For lactic acid:

```text
R=[(1, 'R')], S=[(1, 'S')]
MW=90.0780, TPSA=57.53, HBD=2, HBA=2
enantiomer InChIKeys differ: True
R SMILES round trip preserves configuration: True
S SMILES round trip preserves configuration: True
```

Canonical SMILES is a toolkit-dependent canonical representation. InChIKey is a
compact identifier derived from InChI. Neither replaces the original input,
standardisation record, and provenance.

## 10. Common misconceptions

- **R/S is not the same as d/l.** d/l is an older optical-rotation convention
  and does not map directly to R/S.
- **R/S is not the sign of optical rotation.** The label describes
  configuration, not the direction of rotation.
- **E/Z is not always cis/trans.** The terms coincide only for simple
  two-substituent alkenes.
- **A molecule with stereocentres is not necessarily chiral.** Meso compounds
  are achiral overall.
- **The `@`/`@@` text is not an absolute label.** It depends on atom order.
- **Formula and descriptors do not distinguish stereoisomers.**
- **Unspecified stereochemistry is not a racemate by default.**

## 11. Acceptance criteria

A stereochemistry module is correct when:

- R/S labels are assigned from CIP priorities, not inferred from `@`/`@@`;
- enantiomers are recognized as non-superimposable mirror images;
- diastereomers are distinguished from enantiomers;
- E/Z is assigned from substituent priorities;
- meso compounds are recognized as achiral overall;
- conformers are distinguished from configurational stereoisomers;
- stereochemical identity is preserved in canonical SMILES and InChIKey;
- CIP labels are assigned only after stereochemistry is specified or
  determined;
- a stereochemical SMILES round trip preserves the intended configuration;
- unspecified stereochemistry is not silently treated as a single enantiomer
  or a racemate;
- the connection to chiral protein recognition is explicit.

## 12. Exercises

### Inspection

1. Assign R/S to the two lactic-acid enantiomers.
2. Identify the stereocentres in tartaric acid and classify the R,R, S,S, and
   meso forms.
3. Assign E/Z to 2-butene and to a tri-substituted alkene.
4. Find a molecule with a stereocentre but no net chirality.
5. Compare the InChIKeys of an enantiomer pair and a diastereomer pair.

### Reasoning

1. Why is a conformer not a configurational stereoisomer?
2. Why does the `@`/`@@` marker depend on atom order?
3. Why can two molecules with identical formulas have different activities?
4. Why is unspecified stereochemistry ambiguous?
5. Why is E/Z more general than cis/trans?
6. Why is a meso compound achiral despite having stereocentres?
7. Why must docking use the intended stereoisomer?

### Implementation

1. Run `examples/stereochemistry.py` and reproduce the output.
2. Change the atom order in a chiral SMILES and confirm that the R/S label is
   unchanged after assignment.
3. Compare canonical SMILES and InChIKey for R- and S-lactic acid.
4. Generate a molecule with unspecified stereochemistry and explain what the
   record does and does not say.
5. Enumerate the stereoisomers of a molecule with two stereocentres and identify
   any meso forms.
6. Round-trip a chiral SMILES through canonicalisation and confirm that the
   configuration and InChIKey are preserved.
7. Confirm that enantiomers share common achiral descriptors but differ in
   chiral recognition.

## 13. Summary

- Stereochemistry is part of chemical identity.
- Conformers interconvert by rotation; configurational stereoisomers require
  bond breaking or a high-barrier process.
- R/S describes tetrahedral configuration; E/Z describes double-bond
  configuration.
- Enantiomers are mirror images; diastereomers are not.
- Meso compounds contain stereocentres but are achiral overall.
- Formula and many descriptors cannot distinguish stereoisomers.
- SMILES `@`/`@@` is order-dependent; use a stereochemistry assignment step for
  CIP labels.
- Proteins are chiral environments, so stereochemistry affects recognition.
- Unspecified stereochemistry is not a racemate by default.

**Previous:** [04 — Protonation and tautomerism](04_protonation_and_tautomerism.md).

**Next:** [06 — Molecular polarity and electronic properties](06_molecular_polarity.md).
