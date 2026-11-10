# Track C, module 02 — Binding-site geometry

> **Where this fits.** [Module 01](01_protein_structure.md) established what a
> protein structure file contains. This module asks which part of that
> structure surrounds a reference ligand, and how to describe that region
> geometrically without confusing proximity with binding.
>
> **What you should already know.** Track A sections 1–10 for positions and
> distance; Track C module 01 for structure representation.
>
> **What you will be able to do afterwards.** Select protein atoms around a
> ligand, aggregate them into residues, compute a pocket centre and bounding
> box, and explain the assumptions behind each geometric choice.

## 1. Learning objectives

**Core**

- Define a binding site relative to a reference ligand.
- Compute atom-to-atom distances in nanometres.
- Select protein atoms using a distance cutoff.
- Aggregate selected atoms into residues without duplicates.
- Compute a pocket centre, bounding box, and dimensions.

**Deeper insight**

- A cutoff is a geometric convention, not a chemical interaction.
- A pocket selection is specific to the reference ligand and cutoff.
- Alternate locations, hydrogens, waters, ions, and cofactors need explicit
  policies.
- Ordinary crystallographic structures do not imply periodic minimum-image
  distances.
- Solvent exposure and structural uncertainty affect the interpretation of a
  pocket.

**Looking ahead**

- Module 03 will classify hydrogen bonds, salt bridges, aromatic contacts, and
  hydrophobic contacts from geometry.
- A geometric contact is not an energy contribution.
- Docking will eventually search poses inside or near this region.

## 2. What is a binding site?

A binding site is a region of a protein where a ligand, cofactor, ion, or
another protein binds. For a structure that already contains a ligand, the
simplest operational definition is:

> the protein atoms near the reference ligand.

That definition is deliberately geometric. It does not tell you whether the
ligand is a drug, whether the site is functional, or whether the interaction is
favourable. It only identifies a neighbourhood.

A binding site can also be defined without a ligand:

- from a known pocket-detection algorithm;
- from a known residue set;
- from a cavity in the protein surface;
- from a conserved structural motif.

This module uses a reference ligand because it gives a reproducible,
inspectable starting point.

### 2.1 Ligand-centred versus residue-based definitions

A **ligand-centred** definition starts from ligand atoms and selects nearby
protein atoms. It is local and reproducible, but it depends on the ligand pose.

A **residue-based** definition lists known pocket residues. It is useful for
comparing several structures, but it depends on prior annotation.

Both can be useful. The important thing is to state which one is being used.

## 3. Defining a pocket relative to a reference ligand

Let protein atom $i$ have position $\mathbf r_i$ and ligand atom $j$ have
position $\mathbf r_j$, both in nm. The Euclidean distance is

$$ d_{ij} = \lVert \mathbf r_i - \mathbf r_j \rVert = \sqrt{(x_i-x_j)^2+(y_i-y_j)^2+(z_i-z_j)^2}. $$

For each protein atom, the minimum distance to the ligand is

$$ d_i^{\min} = \min_j d_{ij}. $$

For a cutoff $r_c$, a simple pocket selection is

$$ \mathcal P(r_c) = \{i \in \text{protein atoms} : d_i^{\min} < r_c\}. $$

The strict inequality is a convention. The implementation exposes
`inclusive=True` for $d_i^{\min} \le r_c$.

This is an $O(NM)$ calculation: every protein atom is compared with every
ligand atom. For the small structures in this course that is appropriate.
A neighbour list or spatial index can come later if profiling shows a need.

### 3.1 Coordinates and periodic boundaries

The geometry layer assumes Cartesian coordinates in nm. It does **not** apply
minimum-image wrapping.

Crystallographic symmetry and simulation periodicity are different issues:

- a crystal structure may have symmetry mates that are not in the deposited
  coordinate file;
- a simulation box may use periodic minimum images;
- a binding-site calculation on a deposited structure usually should not
  silently wrap coordinates.

If periodic distances are required, make that a separate, explicit decision.

## 4. Atom-to-atom and residue-to-ligand distances

The first output is one minimum distance per protein atom:

$$ d_i^{\min} = \min_j d_{ij}. $$

The runnable example gives:

```text
Atom minimum distances (nm)
   ALA A 1: 0.100
   ALA A 1: 0.200
   GLY A 2: 0.350
   SER A 3: 0.400
   VAL A 4: 1.000
```

Atom-level selection is the foundation. Residue-level selection comes next.

### 4.1 Residue aggregation

A residue belongs to the selected pocket if at least one of its eligible atoms
satisfies the cutoff. This rule avoids requiring every atom of a side chain to
be close to the ligand.

The residue identity must include the chain and insertion code:

```text
chain identifier + residue number + insertion code + residue name
```

The geometry module does not invent residue identities. The structure adapter
supplies residue labels, and `group_indices_by_residue` groups the selected
atom indices without duplicating a residue.

## 5. Distance cutoffs and their physical meaning

A cutoff is a geometric convention. It is not a hydrogen-bond criterion, a
force-field cutoff, or an energy threshold.

Different questions need different cutoffs:

- a tight contact list might use a smaller cutoff;
- a pocket definition might use a larger cutoff;
- a hydrogen-bond analysis also needs donor/acceptor identity and an angle;
- a hydrophobic-contact analysis needs carbon or hydrophobic atom types;
- an energy calculation needs a force field, not a distance alone.

The example uses 0.5 nm for a compact teaching pocket. That number is not a
universal definition of a binding site.

### 5.1 Strict versus inclusive

The default selection is strict:

$$ d_i^{\min} < r_c. $$

With `inclusive=True`, the boundary is included:

$$ d_i^{\min} \le r_c. $$

The distinction matters at the boundary:

```text
Cutoff boundary
  strict <  0.500: []
  inclusive <= 0.500: [0]
```

Report the convention with the result. Two calculations that differ only at
the boundary should not be compared as if they used the same rule.

## 6. Pocket centre, bounding box, and dimensions

Once atoms are selected, a simple geometric summary is:

- **centre**: the arithmetic mean of the selected atom positions;
- **bounding box**: the axis-aligned minimum and maximum coordinates;
- **dimensions**: the side lengths of that box.

The centre is not the ligand centre and not the centre of the cavity. It is the
mean of the selected protein atoms, so it changes when the selection or atom
set changes.

The example gives:

```text
Pocket selection (strict, cutoff 0.500 nm)
  selected atom indices: [0, 1, 2, 3]
  selected atoms: 4
  residues: ALA A 1, GLY A 2, SER A 3
  centre: (0.175, 0.050, 0.087) nm
  box dimensions: (0.600, 0.200, 0.350) nm
```

The bounding box is axis-aligned. Rotating the coordinate frame changes the box
dimensions even though the internal distances do not change. For rotation-
invariant shape descriptors, use a different representation.

## 7. Residue-level pocket composition

Residue composition answers:

- which residues line the pocket?
- which chains contribute?
- is the site mostly hydrophobic, polar, or charged?
- which residues are conserved across structures?

The last two questions require chemistry and sequence information. This module
only supplies the residue list.

In the example, the selected atoms belong to:

```text
ALA A 1, GLY A 2, SER A 3
```

The repeated ALA atoms are grouped into one residue. The residue list is not an
interaction fingerprint: it says which residues are near the ligand, not which
ones hydrogen-bond to it.

## 8. Waters, ions, cofactors, and alternate locations

A structure contains more than protein atoms. The geometry layer does not decide
which of them belong in a pocket calculation.

### 8.1 Waters

Crystallographic waters may be:

- structural waters that mediate contacts;
- bulk-like waters trapped in the crystal;
- partially occupied or uncertain;
- absent from the model even when present in solution.

Including or excluding waters is a modelling decision. If they are included,
keep them labelled separately from protein atoms.

### 8.2 Ions and cofactors

Metal ions, chloride, phosphate, and cofactors can be essential for a binding
site. Do not discard them merely because they are hetero records.

### 8.3 Alternate locations

Alternate locations represent mutually exclusive modelled positions. Do not
count both conformations as simultaneously occupied atoms unless the analysis
explicitly intends to treat them as an ensemble.

The geometry module receives whatever positions the adapter supplies. The
adapter must record the policy:

- preserve all alternate locations;
- select one by occupancy;
- select one by an explicit rule;
- reject ambiguous structures.

### 8.4 Hydrogens

Hydrogen atoms are often absent from X-ray structures. If they are present, the
selection should state whether they are included. A heavy-atom pocket
definition is not the same as an all-atom definition.

## 9. Solvent accessibility and pocket exposure

A distance cutoff can select a deep cavity or a shallow surface groove. The same
number of atoms does not imply the same kind of site.

Useful geometric context includes:

- solvent-accessible surface area (SASA);
- buried versus exposed atoms;
- pocket volume;
- distance from the protein surface;
- shape complementarity between ligand and pocket.

This module does not yet implement SASA or volume. It provides the selection and
bounding box that those analyses would build on.

A geometrically plausible pocket is not automatically a biologically relevant
binding site. It could be:

- a crystal contact;
- a surface groove;
- a site that binds only a fragment;
- a site whose importance depends on a conformational change.

## 10. Structural uncertainty and limitations

A binding-site geometry inherits all the uncertainty of the structure:

- missing atoms and residues;
- alternate locations;
- low resolution or poor model quality;
- a predicted rather than experimental structure;
- a ligand pose that may be uncertain;
- a cutoff chosen for convenience.

It also inherits the limitations of the definition:

- a reference ligand selects only its own neighbourhood;
- a cutoff is not a thermodynamic boundary;
- proximity is not an interaction;
- a residue list is not an energy decomposition;
- no periodic minimum image is applied unless explicitly requested.

## 11. Python implementation and validation

The implementation is in `src/moleculab/analysis/binding_site.py`. It depends
only on NumPy and Cartesian positions.

### 11.1 API

| Function | Purpose |
| --- | --- |
| `minimum_distances(protein_positions, ligand_positions)` | one minimum distance per protein atom |
| `select_pocket_atoms(protein_positions, ligand_positions, cutoff, inclusive=False)` | protein atom indices inside the cutoff |
| `pocket_geometry(protein_positions, atom_indices, cutoff=None)` | centre, bounding box, and dimensions |
| `binding_site_geometry(protein_positions, ligand_positions, cutoff, inclusive=False)` | selection and geometry in one call |
| `group_indices_by_residue(atom_indices, residue_ids)` | residue-level grouping without duplicates |

### 11.2 Runnable example

```bash
uv run python examples/drug_discovery/binding_site.py
```

The output above comes from this example.

### 11.3 Tests

```bash
uv run pytest tests/test_binding_site.py -v
```

The tests cover:

- exact distances for simple configurations;
- empty protein and ligand arrays;
- output shape and validation;
- strict and inclusive cutoff behaviour;
- translation and rotation invariance;
- pocket centre and bounding box;
- residue aggregation without duplicate residues;
- no implicit periodic minimum image;
- coordinate-scale consistency.

### 11.4 Acceptance criteria

A binding-site geometry module is correct when:

- distances are computed in nm and match hand calculations;
- selection uses a documented strict or inclusive cutoff;
- residue aggregation includes a residue if at least one eligible atom is
  selected;
- residue groups contain no duplicates;
- translation and rotation leave internal distances unchanged;
- alternate-location and hydrogen policies are explicit at the adapter
  boundary;
- no periodic minimum image is applied silently;
- pocket geometry is reported as geometry, not as affinity or interaction
  energy.

## 12. Exercises

### Calculation

1. Compute $d_i^{\min}$ by hand for the five protein atoms in the example.
2. Which atoms enter the pocket for $r_c = 0.3$ nm? Which residues?
3. Change the cutoff to 0.45 nm and predict the selection before running it.
4. Add a protein atom exactly at 0.5 nm and compare strict and inclusive
   selection.
5. Compute the centre and bounding box for a different atom selection.

### Reasoning

1. Why is a binding site defined relative to a reference ligand not a complete
   list of all possible sites?
2. Why is proximity not proof of a hydrogen bond?
3. Why should alternate locations not be counted as simultaneously occupied?
4. Why is a crystal structure not automatically a periodic simulation box?
5. Why does a pocket bounding box change under rotation while distances do not?
6. Why is a residue list not an interaction fingerprint?

### Implementation

1. Run `examples/drug_discovery/binding_site.py` and reproduce the output.
2. Add a water position and decide whether it belongs in the selection.
3. Add residue labels with insertion codes and verify that grouping remains
   unambiguous.
4. Apply a translation and rotation to both protein and ligand and confirm that
   distances and selected atom indices are unchanged.
5. Compare the pocket geometry for two different ligand poses.

## 13. Summary

- A binding site can be defined geometrically as the protein atoms near a
  reference ligand.
- Atom-level selection uses the minimum protein–ligand distance and a documented
  cutoff.
- Residues join the pocket if at least one eligible atom is selected.
- The centre, bounding box, and dimensions are geometric summaries, not
  binding affinities.
- Waters, ions, cofactors, alternate locations, and hydrogens require explicit
  policies.
- Ordinary crystallographic coordinates do not imply periodic minimum images.
- Solvent exposure and structural uncertainty affect how a pocket should be
  interpreted.
- The next step is interaction geometry: hydrogen bonds, salt bridges, aromatic
  contacts, and hydrophobic contacts.

**Previous:** [01 — Protein structure fundamentals](01_protein_structure.md).

**Next:** [03 — Protein–ligand interaction geometry](03_protein_ligand_interactions.md).
