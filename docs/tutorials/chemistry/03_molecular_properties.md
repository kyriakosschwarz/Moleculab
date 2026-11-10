# Molecular properties and chemical space

> **Where this fits.** Module 02 showed how to build and inspect a chemically
> meaningful molecule. This module asks what can be calculated from that
> molecule and how compounds can be compared.
>
> **What you should already know.** Molecular graphs, SMILES, functional
> groups, and conformers.
>
> **What you will be able to do afterwards.** Interpret common molecular
> descriptors, explain Morgan fingerprints, calculate Tanimoto similarity, and
> discuss why structural similarity does not guarantee biological similarity.

## 1. Learning objectives

**Core**

- Interpret molecular weight, heavy-atom count, rotatable bonds, hydrogen-bond
  donors and acceptors, TPSA, cLogP, and formal charge.
- Explain what a molecular fingerprint encodes.
- Calculate Tanimoto similarity for binary vectors.
- Compare molecules with RDKit fingerprints.

**Deeper insight**

- Distinguish average molecular weight from monoisotopic mass.
- Distinguish graph-based descriptors from conformer-dependent properties.
- Explain why descriptors are models, not measurements.
- Explain why high structural similarity does not imply similar activity.
- Explain the idea of chemical space and scaffold hopping.

**Looking ahead**

- Protonation changes formal charge and descriptors.
- Solvation and LogP/LogD require thermodynamics, not only graph descriptors.
- Pharmacophores describe spatial feature arrangements that fingerprints do
  not capture.

## 2. Basic descriptors

Descriptors are numbers computed from a molecular representation. They are
useful for filtering, comparison, and hypothesis generation. They are not
measurements of biological activity.

| Descriptor | Interpretation | Limitation |
| --- | --- | --- |
| molecular formula | element counts | no connectivity or geometry |
| molecular weight | average molar mass | average, not a unique mass |
| exact/monoisotopic mass | mass of the most abundant isotopes | not the average sample mass |
| heavy atom count | non-hydrogen atoms | ignores hydrogen bonding and shape |
| rotatable bonds | approximate flexibility | not a direct entropy calculation |
| H-bond donors | donor groups | environment and geometry matter |
| H-bond acceptors | acceptor groups | environment and geometry matter |
| TPSA | approximate polar surface area | graph-based, not a full solvent model |
| cLogP | predicted lipophilicity | empirical model, pH-independent form |
| formal charge | net integer charge | not a partial-charge distribution |

### 2.1 Average versus monoisotopic mass

RDKit gives both:

- `Descriptors.MolWt` is the average molar mass, using standard atomic weights.
- `Descriptors.ExactMolWt` is the monoisotopic mass, using the most abundant
  isotope of each element.

For aspirin, the values are about 180.16 and 180.0423 g/mol. They are both
useful, but they are not interchangeable. Mass spectrometry usually reports
monoisotopic or isotopic masses; a macroscopic sample has an average molar
mass.

### 2.2 Rotatable bonds are not conformational entropy

A rotatable bond count is a quick structural indicator. It does not compute
the number of accessible conformations, their energies, or the entropy. An
amide bond is usually counted as restricted; a ring bond is not freely
rotatable. The descriptor is a heuristic, not a partition function.

### 2.3 cLogP is a prediction

cLogP is an empirical prediction of the partition coefficient of the neutral
compound. It is not a measured LogP, and it does not include the pH-dependent
distribution of ionized states. LogD is a different quantity.

### 2.4 Formal charge is not partial charge

A carboxylate has formal charge -1. Its partial charges are distributed over
the atoms. The force field uses partial charges; the chemistry layer uses formal
charge and protonation state.

### 2.5 Conformational properties require 3D sampling

The descriptors above are mostly graph-based: molecular weight, TPSA, cLogP,
and rotatable-bond count can be computed without coordinates. Other properties
depend on the three-dimensional structure and on which conformers are populated:

- end-to-end distance and radius of gyration;
- solvent-accessible surface area and molecular shape;
- intramolecular hydrogen-bond patterns;
- strain energy;
- conformational entropy;
- the spatial arrangement of pharmacophore features.

One embedded or minimised conformer is one sample, not a Boltzmann-weighted
average. A property computed from that sample is a model result conditional on
the chosen conformer and force field.

Two molecules with the same number of rotatable bonds can still have very
different conformational landscapes. Ring constraints, amide planarity, steric
clashes, and intramolecular hydrogen bonds all change which conformations are
accessible and how much configurational volume they occupy.

## 3. Molecular fingerprints

A fingerprint encodes structural features as a numerical vector. A simple
binary fingerprint is

$$ A=(1,0,1,1,0,0,1,\ldots). $$

Each bit says whether a feature is present. Morgan fingerprints are built from
atom-centred environments: for each atom, the toolkit records the atom and its
neighbours, then expands the neighbourhood radius.

- radius 0: the atom itself;
- radius 1: immediate neighbours;
- radius 2: neighbours of neighbours.

The feature identifiers are hashed into a fixed-length bit vector. Hashing can
cause collisions, so a fingerprint is a representation, not a unique chemical
identity.

## 4. Tanimoto similarity

For two binary fingerprints $A$ and $B$, the Tanimoto coefficient is

$$ T(A,B)=\frac{|A\cap B|}{|A\cup B|}. $$

If $a$ is the number of active bits in $A$, $b$ in $B$, and $c$ the number
shared, then

$$ T(A,B)=\frac{c}{a+b-c}. $$

The value ranges from 0 (no shared bits) to 1 (identical active sets).

### 4.1 A hand-checkable example

Take

$$ A=(1,0,1,1,0,0,1), \qquad B=(1,1,0,1,0,1,1). $$

Then $a=4$, $b=5$, and $c=3$ (bits 0, 3, and 6). Therefore

$$ T(A,B)=\frac{3}{4+5-3}=0.5. $$

Do this by hand before trusting a toolkit result.

## 5. A runnable comparison

The example compares aspirin, ibuprofen, naproxen, acetaminophen, and caffeine.

```bash
uv run --with rdkit python examples/molecular_properties.py
```

Descriptor output:

```text
      molecule    formula       MW    exact  HBD  HBA    TPSA  rot charge
       aspirin     C9H8O4   180.16 180.0423    1    3    63.6    2      0
     ibuprofen   C13H18O2   206.28 206.1307    1    1    37.3    4      0
      naproxen   C14H14O3   230.26 230.0943    1    2    46.5    3      0
 acetaminophen   C8H9NO2   151.16 151.0633    2    2    49.3    1      0
      caffeine  C8H10N4O2   194.19 194.0804    0    3    61.8    0      0
```

Morgan fingerprint similarity matrix:

```text
                     aspirin     ibuprofen      naproxen acetaminophen      caffeine
       aspirin         1.000         0.195         0.205         0.222         0.089
     ibuprofen         0.195         1.000         0.421         0.184         0.087
      naproxen         0.205         0.421         1.000         0.167         0.080
 acetaminophen         0.222         0.184         0.167         1.000         0.098
      caffeine         0.089         0.087         0.080         0.098         1.000
```

Ibuprofen and naproxen are the most similar pair here (0.421), which is
chemically plausible: both are arylpropionic-acid NSAIDs. Aspirin and
acetaminophen are much less similar under this fingerprint even though both are
analgesic/anti-inflammatory drugs. That is the key lesson:

> Structural similarity is a hypothesis generator, not a biological equivalence
> test.

## 6. Chemical space and scaffold hopping

Chemical space is the abstract space of possible molecules and their
properties. A compound library occupies a small, uneven part of it.

- **Similarity searching** asks which known compounds resemble a query.
- **Scaffold hopping** asks whether a different molecular scaffold can present
  a similar arrangement of chemical features.
- **Diversity selection** asks which compounds cover different regions of
  chemical space.

Fingerprint similarity is useful for all three, but it is not the same as
three-dimensional pharmacophore matching. Two molecules can have different
graphs and similar spatial features, or similar graphs and different binding
modes.

## 7. What descriptors cannot tell you

Descriptors do not include:

- the protein binding site;
- conformational strain on binding;
- desolvation and water networks;
- entropy;
- protonation at the relevant pH;
- metabolism, permeability, or toxicity;
- experimental uncertainty.

They are filters and hypotheses, not answers.

## 8. Acceptance criteria

A molecular-properties module is correct when:

- descriptor values agree with independent toolkit calculations;
- average and monoisotopic masses are clearly distinguished;
- graph-based descriptors and conformer-dependent properties are distinguished;
- a hand-calculated Tanimoto value matches the implementation;
- fingerprint similarity is symmetric and lies in [0, 1];
- a molecule is most similar to itself;
- the tutorial states that similarity is not biological equivalence.

## 9. Exercises

### Calculation

1. Compute Tanimoto similarity for the two binary vectors in section 4.1 by
   hand.
2. Build two vectors with $a=5$, $b=4$, and $c=2$; compute the similarity.
3. For aspirin, compute the difference between average and monoisotopic mass.
4. Compare the rotatable-bond counts of ibuprofen and caffeine and explain the
   difference qualitatively.

### Reasoning

1. Why is a rotatable-bond count not a direct measure of conformational
   entropy?
2. Why can two molecules with similar fingerprints have different activities?
3. Why can two molecules with different scaffolds share a pharmacophore?
4. Why is cLogP not a measured LogP?
5. Why is formal charge not a partial charge?

### Implementation

1. Run `examples/molecular_properties.py`.
2. Add another molecule and extend the similarity matrix.
3. Vary the Morgan radius and observe how the matrix changes.
4. Compare ECFP-like Morgan fingerprints with a simple atom-count descriptor.
5. Find a pair with moderate similarity but different functional groups and
   explain the result.

## 10. Summary

- Descriptors encode selected structural and physicochemical properties.
- Graph-based descriptors are not the same as conformer-dependent properties.
- Average mass and monoisotopic mass are different quantities.
- Morgan fingerprints encode atom-centred environments and are hashed into a
  bit vector.
- Tanimoto similarity compares binary fingerprints.
- Structural similarity is not biological equivalence.
- Chemical space, similarity searching, and scaffold hopping are hypothesis
  tools, not binding predictions.
- The next step is protonation and tautomerism, where chemical state changes
  the graph, charge, and interactions.

**Previous:** [02 — Molecular representation](02_molecular_representation.md).  
**Next:** [04 — Protonation and tautomerism](04_protonation_and_tautomerism.md).
