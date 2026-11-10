# Molecular representation: from atoms to chemical identity

> **Where this fits.** The previous module used the physics engine to study a
> torsional model. This module asks what a computer actually needs in order to
> treat a molecule as a chemical object rather than a bag of particles.
>
> **What you should already know.** Atoms, bonds, topology, torsions, and the
> difference between a graph and a geometry.
>
> **What you will be able to do afterwards.** Distinguish formulas, graphs,
> SMILES, stereochemistry, and conformers; inspect a molecule with RDKit; and
> explain why chemical perception belongs outside the force-field layer.

## 1. Learning objectives

**Core**

- Distinguish a molecular formula, graph, SMILES string, and 3D conformer.
- Read formal charges, bond orders, aromaticity, and stereochemistry.
- Generate and inspect a 3D conformer with RDKit.

**Deeper insight**

- Explain why a graph is not a geometry.
- Explain why a SMILES string does not specify a unique 3D structure.
- Explain why Moleculab should use a cheminformatics toolkit for perception
  rather than reimplementing SMARTS or aromaticity rules.

**Looking ahead**

- Descriptors and fingerprints use the graph representation.
- Protonation and tautomerism change the chemical graph and charges.
- Docking needs a chemically correct ligand graph and a chosen 3D conformer.

## 2. The representations

A molecule can be represented at several levels.

| Representation | What it specifies | What it does not specify |
| --- | --- | --- |
| molecular formula | element counts | connectivity or geometry |
| 2D graph | atoms, bonds, bond orders | a unique 3D conformation |
| formal charges | integer charge per atom | partial charges or polarisation |
| stereochemistry | configuration at stereocentres and double bonds | a particular conformer |
| SMILES | connectivity, bond orders, aromaticity, stereochemistry | a unique 3D structure |
| 3D coordinates | atom positions | chemical perception by itself |
| conformer | one 3D arrangement | the full conformational ensemble |

The central distinction is:

> A molecular graph is not a molecular geometry.

Two conformers have the same graph and different coordinates. A SMILES string
usually describes connectivity and stereochemistry, but not one experimentally
meaningful geometry.

### 2.1 Chemical identity versus representation

A SMILES string is a representation, not a unique identifier. The same molecule
can be written in several valid ways. Ethanol, for example, can be written as
`CCO`, `OCC`, or `C(O)C`. All three parse to the same molecular graph:

```python
from rdkit import Chem

for smiles in ["CCO", "OCC", "C(O)C"]:
    molecule = Chem.MolFromSmiles(smiles)
    print(smiles, "->", Chem.MolToSmiles(molecule))
```

RDKit prints the same canonical SMILES for all three:

```text
CCO -> CCO
OCC -> CCO
C(O)C -> CCO
```

Canonical SMILES is a normalised representation produced by a particular
canonicalisation algorithm. It is not a universal chemical identity standard:
different toolkits, and even different versions of the same toolkit, can choose
different canonical forms. A canonical SMILES is convenient for deduplication
inside one pipeline, but it should not be treated as an absolute identifier
across databases.

This matters as soon as you:

- deduplicate a compound library;
- join two datasets by structure;
- compare stereoisomers or protonation states;
- store a molecule and later reconstruct it.

### 2.2 What a molecular record must retain

A single SMILES string is not enough to describe everything a simulation
workflow needs. A molecular record may need to keep several layers separately:

1. **Original input representation**, with its source and provenance.
2. **Standardised representation**, after an explicit normalisation step.
3. **Stereochemical information**, including whether it was assigned or left
   unspecified.
4. **Protonation and tautomeric state**, since these change the graph and the
   charges.
5. **Generated conformer coordinates**, with the method and settings used.
6. **Toolkit and settings**, so the result can be reproduced.

Not every workflow needs all six fields. The point is that they are different
pieces of information that can otherwise be silently conflated. Moleculab does
not implement a full chemical standardiser yet; this module only establishes why
the distinction exists.

## 3. Why this matters for Moleculab

Moleculab's `Topology` stores bonds, angles, and torsions by atom index. It does
not know whether an atom is carbon or oxygen, whether a ring is aromatic, or
whether a bond is a single or double bond. Those are chemical-perception
questions.

The force field sees:

- atom types;
- partial charges;
- bonded parameters;
- Lennard-Jones parameters;
- exclusions.

The chemistry toolkit sees:

- elements;
- bond orders;
- aromatic rings;
- formal charges;
- stereochemistry;
- functional groups.

Keeping these layers separate prevents a force-field class from becoming a
half-written cheminformatics library.

## 4. Using RDKit as an external chemistry toolkit

RDKit is an established open-source cheminformatics toolkit. It should be an
optional external dependency, not a requirement of the physics engine.

For a one-off session:

```bash
uv run --with rdkit python
```

Or use any environment where RDKit is installed. The code in this module is
intentionally small; the point is the representation, not a new API.

### 4.1 Build a molecule from SMILES

```python
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

mol = Chem.MolFromSmiles("CCO")  # ethanol
print("atoms:", mol.GetNumAtoms())
print("bonds:", mol.GetNumBonds())
print("formula:", rdMolDescriptors.CalcMolFormula(mol))
```

### 4.2 Inspect atoms

```python
for atom in mol.GetAtoms():
    print(
        atom.GetIdx(),
        atom.GetSymbol(),
        "formal charge =", atom.GetFormalCharge(),
        "aromatic =", atom.GetIsAromatic(),
    )
```

### 4.3 Inspect bonds

```python
for bond in mol.GetBonds():
    print(
        bond.GetBeginAtomIdx(),
        bond.GetEndAtomIdx(),
        bond.GetBondType(),
        "aromatic =", bond.GetIsAromatic(),
    )
```

### 4.4 Rings and aromaticity

```python
rings = mol.GetRingInfo().AtomRings()
print("rings:", rings)
aromatic_atoms = [
    atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIsAromatic()
]
print("aromatic atoms:", aromatic_atoms)
```

Aromaticity is a chemical perception result. Different toolkits can use
different aromaticity models. The force field receives the resulting atom types
and parameters, not a quantum-mechanical aromaticity calculation.

### 4.5 Formal charge

```python
print("net formal charge:", Chem.GetFormalCharge(mol))
```

Formal charge is not the same as a partial charge. A carboxylate has a formal
charge of -1, while its partial charges are distributed over the atoms. The
force field needs the partial charges; the chemistry layer needs the formal
charge and protonation state.

### 4.6 Stereochemistry

```python
from rdkit.Chem import FindMolChiralCenters

mol = Chem.MolFromSmiles("C[C@H](O)C(=O)O")  # one lactic-acid enantiomer
print(FindMolChiralCenters(mol, includeUnassigned=True))
```

Stereochemistry changes the 3D shape and can change biological activity. A
docking or simulation workflow must use the intended stereoisomer.

## 5. From graph to 3D conformer

RDKit can generate a 3D conformer:

```python
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

mol = Chem.AddHs(Chem.MolFromSmiles("CCCC"))  # butane with hydrogens
AllChem.EmbedMolecule(mol, randomSeed=7)
conformer = mol.GetConformer()
coordinates = np.array(
    [
        list(conformer.GetAtomPosition(index))
        for index in range(mol.GetNumAtoms())
    ]
)
print(coordinates.shape)
```

The coordinates are one conformer, not the conformational ensemble. A molecule
with several rotatable bonds has many possible 3D structures.

### 5.1 Embedding, optimisation, and conformer sampling

Three operations are easy to confuse:

| Operation | What it does | What it does not guarantee |
| --- | --- | --- |
| embedding | generates initial 3D coordinates consistent with the graph | a low-energy structure |
| geometry optimisation | relaxes coordinates on a force-field energy surface | the global minimum |
| conformer sampling | generates several conformers and selects a representative set | a Boltzmann-weighted ensemble |

`AllChem.EmbedMolecule` embeds one conformer. `AllChem.MMFFOptimizeMolecule`
then relaxes it. `AllChem.EmbedMultipleConfs` generates several starting
conformers.

```python
from rdkit import Chem
from rdkit.Chem import AllChem

molecule = Chem.AddHs(Chem.MolFromSmiles("CCCC"))
AllChem.EmbedMolecule(molecule, randomSeed=7)
AllChem.MMFFOptimizeMolecule(molecule)
print("conformers after optimisation:", molecule.GetNumConformers())

sampled = Chem.AddHs(Chem.MolFromSmiles("CCCC"))
conformer_ids = AllChem.EmbedMultipleConfs(sampled, numConfs=5, randomSeed=7)
print("sampled conformers:", len(conformer_ids))
```

A conformer generated by a molecular-mechanics force field is a model result,
not an experimentally determined structure. Embedding several conformers does
not prove that the global minimum has been found, and it does not by itself
produce equilibrium populations. Those require energy evaluation and, in
general, a free-energy calculation.

### 5.2 Compare two conformers

```python
for seed in [7, 21]:
    mol = Chem.AddHs(Chem.MolFromSmiles("CCCC"))
    AllChem.EmbedMolecule(mol, randomSeed=seed)
    conf = mol.GetConformer()
    positions = np.array(
        [
            list(conf.GetAtomPosition(index))
            for index in range(mol.GetNumAtoms())
        ]
    )
    print(seed, positions[0], positions[-1])
```

The two conformers have the same graph and different coordinates. That is the
graph-versus-geometry distinction in one example. Do not require independently
generated conformers to have identical coordinates; they should share
connectivity and stereochemistry, not positions.

## 6. A small molecule collection

Use these to practise:

| Molecule | SMILES | Chemistry to notice |
| --- | --- | --- |
| ethanol | `CCO` | hydroxyl, hydrogen bonding |
| acetic acid | `CC(=O)O` | carboxylic acid, ionisation |
| benzene | `c1ccccc1` | aromatic ring, planarity |
| aspirin | `CC(=O)Oc1ccccc1C(=O)O` | ester, carboxylic acid, aromatic ring |
| imatinib | `CC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCN(CC3)C)NC4=NC=CC(=N4)C5=CN=CC=C5` | amide, piperazine, pyridine/pyrimidine, aromatic rings |

For each molecule, inspect:

1. atom count and formula;
2. bond orders;
3. aromatic atoms and rings;
4. formal charge;
5. stereocentres;
6. a generated 3D conformer.

For imatinib, check the formula `C29H31N7O`. RDKit reports two different masses:

```text
average molar mass      493.615 g/mol
monoisotopic mass       493.2590086120001 g/mol
```

`Descriptors.MolWt` uses standard atomic weights and gives the average molar
mass of a macroscopic sample. `Descriptors.ExactMolWt` uses the most abundant
isotope of each element and gives the monoisotopic mass used in mass
spectrometry and exact-mass work. They are not interchangeable.

The formula alone also does not specify stereochemistry, protonation state,
tautomeric state, or a unique 3D structure. Matching a formula is a useful
sanity check, not a proof that the chemical record is complete.

If the formula is wrong, the SMILES is wrong. If the formula is right, the
structure can still be the wrong stereoisomer, protonation state, or tautomer.

## 7. The adapter boundary

The intended software boundary is:

```text
RDKit / chemistry toolkit
    chemical perception, graphs, aromaticity, SMILES, conformers
            ↓  adapter
Moleculab chemistry layer
    molecule records, chemical states, descriptors
            ↓
Moleculab physics
    topology, atom types, charges, force field, simulation
```

The adapter should translate between representations without making the force
field depend on RDKit types. It should not create a second competing molecular
graph unless there is a specific reason.

A possible future module layout is illustrative, not a request to create all
files now:

```text
moleculab/chemistry/
├── representation.py
├── descriptors.py
└── conformers.py
```

## 8. What not to implement now

Do not build:

- a custom SMARTS engine;
- a custom aromaticity model;
- a custom pKa predictor;
- a custom protein structure predictor;
- a full ADME model.

Those are substantial projects. Use established toolkits for chemical
perception and focus Moleculab's code on the physics and the analysis that the
curriculum is teaching.

## 9. Acceptance criteria

A molecular-representation module is correct when:

- a molecule built from SMILES has the expected formula and formal charge;
- atom and bond counts match the graph;
- aromatic atoms and rings are identified consistently;
- different valid SMILES for the same graph normalise to the same representation;
- a 3D conformer has one coordinate per atom;
- two conformers of the same molecule share a graph but differ in coordinates;
- an embedded conformer can be optimised without changing connectivity;
- average molar mass and monoisotopic mass are clearly distinguished;
- the physics layer receives atom types, charges, and topology without importing
  RDKit types.

## 10. Exercises

### Identity

1. Parse `CCO`, `OCC`, and `C(O)C`, print the canonical SMILES for each, and
   compare them.
2. Compute the InChIKey for the same three strings and compare.
3. Explain why a canonical SMILES from one toolkit is not a universal chemical
   identifier.
4. List the six fields of a molecular record from section 2.2 and say which ones
   your current workflow would lose.

### Inspection

1. Build ethanol, acetic acid, benzene, aspirin, and imatinib from SMILES.
2. For each, print the formula, atom count, bond count, rings, and formal charge.
3. Identify the functional groups by eye from the graph.
4. Find the stereocentres in lactic acid and compare the two enantiomers.

### Conformers

1. Generate three butane conformers with different random seeds.
2. Compare their end-to-end distances.
3. Generate two imatinib conformers and compare their shapes.
4. Embed one butane conformer, record a coordinate, then optimise it and record
   the same coordinate again. Explain the change.
5. Generate five conformers, compare their energies, and explain why the lowest
   sampled energy is not proof of the global minimum.
6. Explain why the same SMILES can produce different coordinates.

### Reflection

1. Why is a molecular graph not a geometry?
2. Why does a force field need atom types and partial charges rather than only
   a SMILES string?
3. Why should aromaticity perception live outside the force-field module?
4. Why is `Chem.GetFormalCharge` not a partial-charge calculation?
5. Why is an embedded conformer not the same as an equilibrium ensemble?
6. Why is average molar mass not the same as monoisotopic mass?

## 11. Summary

- A formula, graph, SMILES, stereochemical label, and 3D conformer are different
  representations.
- A graph describes connectivity; a conformer describes one geometry.
- A SMILES string is a representation, not a unique identifier; canonical SMILES
  is toolkit-specific.
- A molecular record may need to retain provenance, standardisation,
  stereochemistry, protonation/tautomer state, conformers, and toolkit settings.
- Average molar mass and monoisotopic mass are different quantities.
- Embedding, geometry optimisation, and conformer sampling are distinct
  operations, and generated conformers are model results.
- Formal charge is not a partial charge.
- Aromaticity and stereochemistry are chemical perception results.
- RDKit is an appropriate external toolkit for perception; Moleculab's physics
  should not depend on it directly.
- The adapter boundary keeps chemical interpretation, physical modelling, and
  numerical simulation separate.

**Previous:** [01 — Chemical structure](01_chemical_structure.md).  
**Next:** [03 — Molecular properties](03_molecular_properties.md).
