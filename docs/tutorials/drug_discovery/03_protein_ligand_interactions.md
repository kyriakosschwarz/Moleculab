# Track C, module 03 — Protein–ligand interaction geometry

> **Where this fits.** [Module 02](02_binding_sites.md) selected the protein
> atoms around a ligand. This module asks which noncovalent interactions are
> geometrically plausible inside that region.
>
> **What you should already know.** Track C modules 01–02; Track B module 04
> for protonation and tautomerism.
>
> **What you will be able to do afterwards.** Describe hydrogen-bond geometry,
> explain why distance alone is insufficient, represent chemical features
> independently of a toolkit, and interpret interaction lists without treating
> them as binding free energies.

## 1. Learning objectives

**Core**

- Distinguish geometric proximity from a chemical interaction.
- Describe hydrogen-bond donor, hydrogen, and acceptor geometry.
- Compute and interpret the D-H...A angle.
- Represent donors, acceptors, charges, aromaticity, and hydrophobicity as
  prepared features.
- Detect hydrogen bonds with explicit distance and angle criteria.

**Deeper insight**

- Donor/acceptor identity depends on chemical environment and protonation.
- Salt bridges, aromatic interactions, cation–π contacts, hydrophobic contacts,
  and steric clashes each need a different geometric definition.
- An interaction list is not an energy decomposition.
- Interaction fingerprints, sampled-frame occupancy, and energy contributions
  answer
  different questions.

**Looking ahead**

- Module 04 will turn interaction lists into residue-level fingerprints and
  sampled-frame occupancy.
- Docking will eventually score poses, but a pose score is still not a binding
  free energy.

## 2. Interaction fundamentals

Molecular recognition is not one force. A protein–ligand complex contains a
mixture of:

- hydrogen bonds;
- salt bridges and other electrostatic interactions;
- aromatic stacking;
- cation–π interactions;
- hydrophobic contacts;
- steric repulsion;
- solvent-mediated effects.

Each has a geometric signature, but geometry alone does not determine its
energetic contribution. The same hydrogen-bond distance can correspond to a
strong, weak, or strain-compensated interaction depending on the environment.

This module uses a three-stage pipeline:

```text
structure adapter
    Gemmi, RDKit, PDB/mmCIF
        ↓
chemical feature representation
    element, donor/acceptor, charge, aromaticity, hydrophobicity
        ↓
interaction geometry
    distances, angles, planes, overlaps
        ↓
typed interaction records
```

The interaction layer does not parse structures and does not import a chemistry
toolkit. It consumes a prepared feature representation.

### 2.1 Why a feature representation is necessary

The binding-site module needed only positions. Interaction detection needs
chemical roles:

- an oxygen can be a donor, an acceptor, both, or neither;
- a nitrogen can be neutral, protonated, or part of an aromatic ring;
- a carbon can be hydrophobic or attached to a polar group;
- a carboxylate is charged only in its deprotonated state;
- a histidine ring can change protonation and tautomer state.

Those are not properties of the element alone. They come from the chemical
graph, protonation state, tautomer, and sometimes an explicitly placed hydrogen.

## 3. Hydrogen bonds

A hydrogen bond is usually written D-H...A:

- D is the donor heavy atom;
- H is a hydrogen attached to D;
- A is the acceptor.

The donor–acceptor distance is

$$ d_{DA}=\lVert \mathbf r_D-\mathbf r_A\rVert. $$

The D-H...A angle at the hydrogen is

$$ \theta_{DHA}=\cos^{-1}\left(\frac{(\mathbf r_D-\mathbf r_H)\cdot(\mathbf r_A-\mathbf r_H)}{\lVert \mathbf r_D-\mathbf r_H\rVert\ \lVert \mathbf r_A-\mathbf r_H\rVert}\right). $$

A plausible hydrogen bond generally has:

- a donor and an acceptor with compatible chemistry;
- a donor–acceptor distance within the chosen cutoff;
- a reasonably linear D-H...A arrangement.

The thresholds are operational definitions, not laws of nature.

### 3.1 Why distance alone is insufficient

A short donor–acceptor distance with a poor angle is not a good hydrogen bond.
The example includes a short contact with an unsuitable angle and reports zero
interactions.

Conversely, a well-aligned donor and acceptor can still fail to make a
meaningful contribution if the local environment is hydrophobic or if the
geometry costs strain elsewhere.

### 3.2 Donor and acceptor perception

Do not infer donor/acceptor status from the element alone.

Examples:

- an alcohol oxygen is a donor and an acceptor;
- an ether oxygen is an acceptor but not a donor;
- a protonated amine is a donor;
- a neutral amine can be a donor and an acceptor;
- an amide nitrogen is usually not a good acceptor;
- a carboxylate oxygen is an acceptor and is negatively charged when
  deprotonated.

This is where Track B module 04 becomes essential: protonation and tautomerism
change which atoms can donate or accept.

### 3.3 Hydrogen placement

The angular criterion requires a hydrogen position. X-ray structures often do
not contain hydrogens. The structure adapter must therefore:

- use explicitly modelled hydrogens;
- place hydrogens according to a defined rule;
- or skip the angular criterion and report a distance-only contact with that
  limitation recorded.

The implementation in this module requires hydrogen positions for an angular
hydrogen bond. A donor without hydrogen positions is skipped.

## 4. Salt bridges

A salt bridge is an electrostatic interaction between oppositely charged
groups, such as a carboxylate and an ammonium group.

Geometrically, a first approximation uses the distance between the charged
atoms or charged-group centroids:

$$ d_{+-}=\lVert \mathbf r_{+}-\mathbf r_{-}\rVert. $$

Important details:

- the charges must actually be present at the chosen pH;
- a carboxylate–ammonium pair is not the same as a neutral acid–amine pair;
- a close pair of opposite formal charges may still be weakly contributing if
  the charges are desolvated or geometrically strained;
- a salt bridge is not simply any close pair of atoms with opposite formal
  charges in a force field.

Salt-bridge detection is a natural next detector after hydrogen bonds because
it reuses the same feature representation but does not require hydrogen
positions.

## 5. Aromatic interactions

Aromatic interactions are described with ring geometry:

- ring centroid;
- ring normal;
- angle between two ring normals;
- lateral displacement between centroids.

Common limiting arrangements are:

- **parallel stacking**: normals are nearly parallel, with the rings offset
  laterally;
- **edge-to-face**: one ring normal is roughly perpendicular to the other.

A centroid distance alone is not enough. Two rings can be close but poorly
oriented, and two aromatic rings can interact through different electronic
arrangements.

This module specifies the geometry but does not yet implement an aromatic-ring
detector.

## 6. Cation–π interactions

A cation–π interaction involves a positive charge near the face of an aromatic
ring.

Useful geometric quantities:

- cation-to-ring-centroid distance;
- angle between the cation–centroid vector and the ring normal;
- distance of the cation from the ring plane.

A cation near the edge of a ring is not the same as a cation over the ring face.
As with salt bridges, the charge state must be explicit.

## 7. Hydrophobic contacts

Hydrophobic contacts are usually described as proximity between nonpolar atoms
or groups.

A simple geometric criterion is an atom–atom distance below a cutoff, but:

- not every carbon–carbon contact is a meaningful hydrophobic contribution;
- buried surface area and desolvation matter;
- a contact can be part of a larger packing arrangement;
- polar atoms nearby can change the interpretation.

Hydrophobic-contact detection is a natural extension of the feature
representation once hydrophobic atom flags are available.

## 8. Steric clashes

A steric clash is an unfavourable overlap between atoms. A simple geometric
measure compares the interatomic distance with the sum of van der Waals radii:

$$ \mathrm{overlap}=\frac{r_i+r_j-d_{ij}}{r_i+r_j}. $$

A positive overlap indicates interpenetration of the van der Waals spheres. It
is a geometric flag, not a calculated repulsive energy.

A force-field repulsion depends on the specific potential, atom types, and
parameters. Do not report a clash flag as an energy.

## 9. Interaction fingerprints

An interaction fingerprint records which interactions are present in a pose or
structure. A useful record contains:

- interaction type;
- protein atom or residue identifier;
- ligand atom identifier;
- relevant distance;
- relevant angle or plane geometry;
- the detection criteria;
- chemical-feature provenance.

Keep three concepts distinct:

| Concept | Question |
| --- | --- |
| interaction fingerprint | which interactions are observed? |
| sampled-frame occupancy | how often does each interaction appear across frames? |
| energy decomposition | how much does each term contribute under a physical model? |

They are related but not interchangeable. A contact table does not prove an
energy contribution.

## 10. Units and conventions

Track C uses nm for internal coordinates and degrees for reported angles. The
wider drug-discovery literature uses Angstrom, kcal/mol, micromolar and
nanomolar concentrations, and pIC50 values.

The reference page [Units and conventions](units.md) covers:

- Angstrom to nm;
- degrees and radians;
- kJ/mol and kcal/mol;
- M, mM, µM, nM, and pIC50;
- pH and pKa;
- bar and atm;
- Da and kDa;
- Å² and Å³;
- occupancy and B-factors.

Convert at the boundary and keep one unit system inside the calculation.

## 11. Limitations

Geometric interaction detection has clear limits:

- a detected contact is not proof of a favourable interaction;
- thresholds are conventions;
- donor/acceptor assignment depends on protonation and tautomerism;
- hydrogens may be absent or placed by a model;
- waters and cofactors can mediate or block interactions;
- alternate locations create mutually exclusive geometries;
- a single static structure does not show dynamics or occupancy;
- interaction counts are not binding affinities.

The strongest statement a geometric detector can make is:

> under these structural and chemical assumptions, these atoms satisfy this
> operational interaction definition.

## 12. Implementation and validation

The first implementation is in `src/moleculab/analysis/interactions.py`. It
provides the chemical-feature record and hydrogen-bond detection.

### 12.1 Feature record

`InteractionAtom` contains:

- position in nm;
- element, atom type, residue identifier, and atom name;
- donor and acceptor flags;
- formal charge;
- aromatic and hydrophobic flags;
- hydrogen positions for donors.

The flags are supplied by the caller. The class does not infer chemistry.

`HydrogenBond` records the chosen hydrogen position and its index in the
donor's hydrogen list, so the selection of the most linear hydrogen can be
reproduced.

### 12.2 Hydrogen-bond API

```python
from moleculab import (
    HydrogenBondCriteria,
    InteractionAtom,
    detect_hydrogen_bonds,
)

protein_donor = InteractionAtom(
    position=[0.0, 0.0, 0.0],
    element="N",
    name="NZ",
    residue_id="LYS A 1",
    is_donor=True,
    hydrogen_positions=([0.1, 0.0, 0.0],),
)
ligand_acceptor = InteractionAtom(
    position=[0.30, 0.0, 0.0],
    element="O",
    name="O1",
    residue_id="LIG B 1",
    is_acceptor=True,
)

bonds = detect_hydrogen_bonds(
    [protein_donor],
    [ligand_acceptor],
    HydrogenBondCriteria(max_distance=0.35, min_angle_degrees=120.0),
)
```

The runnable example is:

```bash
uv run python examples/interactions.py
```

Observed output:

```text
Hydrogen-bond detection
  type: hydrogen_bond
  protein: LYS A 1 NZ
  ligand:  LIG B 1 O1
  donor is protein: True
  donor-acceptor distance: 0.300 nm
                        = 3.000 Angstrom
  D-H...A angle: 180.0 degrees

Short distance, unsuitable angle
  interactions: 0

Donor without hydrogen positions
  interactions: 0

Ligand donor to protein acceptor
  donor is protein: False
  distance: 0.300 nm
  angle: 180.0 degrees
```

### 12.3 Criteria

The defaults are:

| Criterion | Default | Comparison |
| --- | --- | --- |
| donor–acceptor distance | 0.35 nm | inclusive |
| D-H...A angle | 120° | inclusive |

These are teaching defaults, not universal constants. Report the criteria with
the result.

### 12.4 Tests

```bash
uv run pytest tests/test_interactions.py -v
```

The tests cover:

- a linear hydrogen bond;
- ligand-donor to protein-acceptor direction;
- invalid donor/acceptor combinations;
- a short distance with an unsuitable angle;
- distance and angle boundaries;
- selection of the most linear hydrogen;
- donors without hydrogen positions;
- translation and rotation invariance;
- reproducibility under atom reordering;
- feature and criteria validation.

### 12.5 What this milestone does not implement

The module specifies the geometry of salt bridges, aromatic stacking,
cation–π interactions, hydrophobic contacts, and steric clashes. Their
detectors are not implemented yet. Keeping them separate avoids pretending that
an element lookup is chemical perception.

### 12.6 Acceptance criteria

A hydrogen-bond detector is correct when:

- it requires a donor, an acceptor, and a hydrogen position;
- it checks both protein-donor and ligand-donor directions;
- it uses a documented distance and angle criterion;
- it rejects a short contact with an unsuitable angle;
- it chooses a hydrogen deterministically when several are present;
- it records the chosen hydrogen position and index;
- it reports distance, angle, atom identifiers, and criteria;
- translation and rotation do not change the detected geometry;
- atom ordering does not change the set of detected interactions;
- it returns geometric records, not a binding score.

## 13. Exercises

### Geometry

1. Compute the D-H...A angle for a linear arrangement.
2. Construct a short donor–acceptor contact with an angle below the cutoff and
   confirm that it is rejected.
3. Move the acceptor until the distance is exactly at the cutoff. Is it
   selected under the inclusive rule?
4. Rotate the entire complex and confirm that distance and angle are unchanged.
5. Add a second hydrogen to a donor and predict which one will be selected.

### Chemistry

1. Explain why an ether oxygen is an acceptor but not a donor.
2. Explain how protonation changes a carboxylate–amine pair.
3. Explain why an amide nitrogen is usually a poor acceptor.
4. Explain why a histidine tautomer can change hydrogen-bond direction.
5. Explain why a close carbon–carbon contact is not automatically a
   hydrophobic interaction.

### Implementation

1. Run `examples/interactions.py` and reproduce the output.
2. Add a negative example with two donors and no acceptor.
3. Add a second ligand acceptor and verify the sorted result order.
4. Replace the hydrogen position and observe the angle change.
5. Write a small adapter that converts your own structure record into
   `InteractionAtom` objects.

## 14. Summary

- Noncovalent interactions have geometric signatures but are not determined by
  distance alone.
- Hydrogen bonds require donor/acceptor chemistry, a donor–acceptor distance,
  and a D-H...A angle.
- Donor/acceptor perception depends on protonation, tautomerism, and hydrogen
  placement.
- Salt bridges, aromatic interactions, cation–π contacts, hydrophobic contacts,
  and steric clashes each need their own geometric definition.
- Interaction fingerprints, sampled-frame occupancy, and energy decomposition
  are different
  concepts.
- The first implementation provides a toolkit-independent feature record and a
  tested hydrogen-bond detector.
- The next module will turn interaction lists into residue-level fingerprints
  and sampled-frame occupancy.

**Previous:** [02 — Binding-site geometry](02_binding_sites.md).

**Next:** [04 — Interaction fingerprints and sampled-frame occupancy](04_interaction_fingerprints.md).
