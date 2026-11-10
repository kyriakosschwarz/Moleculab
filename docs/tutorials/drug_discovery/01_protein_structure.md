# Track C, module 01 — Protein structure fundamentals

> **Where this fits.** This is the first module in
> [Track C](README.md). Track A supplied the physics engine; Track B supplied
> chemical identity and state. This module asks what a protein structure file
> actually contains before any docking score or free-energy calculation.
>
> **What you should already know.** Track A sections 1–10, especially
> particles, `System`, and topology; Track B modules 01–02 for molecular
> representation.
>
> **What you will be able to do afterwards.** Read the hierarchy of a protein
> structure, distinguish experimental from predicted structures, identify
> missing and alternate data, preserve provenance, and convert coordinates from
> Angstrom to nanometre at the import boundary.

## 1. Learning objectives

**Core**

- Distinguish atoms, residues, chains, secondary structure, domains, and
  biological assemblies.
- Read chain identifiers, residue numbers, insertion codes, and alternate
  locations.
- Distinguish experimental structures from predicted models.
- Recognize missing atoms, missing residues, waters, ions, cofactors, and
  multiple models.
- Convert imported coordinates from Angstrom to nanometre.

**Deeper insight**

- A residue number is a label, not a guaranteed sequence index.
- A structure file is a model with uncertainty, not a complete molecule.
- Missing information must be represented as missing, not silently filled with
  zero or rebuilt without a record.
- Provenance is part of the data, not an optional comment.

**Looking ahead**

- Binding-site geometry depends on which chains, waters, and alternate
  locations are retained.
- Protein–ligand interactions require a chosen protonation state.
- Predicted structures need confidence and provenance, not just coordinates.

## 2. The structural hierarchy

Protein structure is described at several levels:

| Level | What it describes |
| --- | --- |
| primary structure | the amino-acid sequence |
| secondary structure | local hydrogen-bonded patterns: helices, sheets, turns, loops |
| tertiary structure | the three-dimensional fold of one chain |
| quaternary structure | the arrangement of several chains |
| domain | a compact, often independently folding region |
| biological assembly | the functional oligomeric form, which may differ from the asymmetric unit |

A structure file represents this hierarchy with atoms grouped into residues,
residues grouped into chains, and chains grouped into models or assemblies.

### 2.1 Atoms, residues, and chains

An atom record contains at least:

- an atom name, such as `CA` or `CB`;
- an element;
- a residue name, such as `ALA` or `HOH`;
- a chain identifier, such as `A`;
- a residue number and insertion code;
- coordinates;
- occupancy and B-factor, when available.

A residue is not the same as an amino acid in a sequence alignment. The same
residue number can appear in different chains, and insertion codes distinguish
residues that share a number. Chain identifiers are local to the file; `A` in
one structure is not the same chain as `A` in another.

A stable residue identity therefore needs more than the residue number:

```text
chain identifier + residue number + insertion code + residue name
```

### 2.2 Sequence versus structure

The sequence tells you which residues should be present. The structure tells
you which residues were observed, modelled, or built. They are not the same
thing.

A chain can have:

- residues that are present in the sequence but absent from the model;
- residues with only a backbone or only some side-chain atoms;
- residues with alternate locations;
- residues modelled from a template rather than directly observed.

A gap in residue numbering is a prompt to inspect the sequence and metadata. It
is not proof by itself that a residue is missing: numbering schemes can skip
numbers for historical or construct reasons.

## 3. Experimental and predicted structures

Different methods establish different things.

| Source | What it provides | Important limitations |
| --- | --- | --- |
| X-ray crystallography | electron density from a crystal; atomic coordinates and B-factors | crystal packing, missing flexible regions, ligand and water modelling, radiation damage |
| cryo-EM | density maps and fitted models, often for large assemblies | local resolution varies, conformational heterogeneity, model fitting choices |
| NMR | an ensemble of models in solution | lower resolution for large systems, ensemble interpretation, fewer restraints |
| homology model | a structure built from related templates | accuracy depends on template identity and alignment |
| predicted model | coordinates inferred from sequence and learned patterns | confidence varies by region; no experimental density; ligands and waters are not automatically present |

A predicted structure is not an experimental structure. It can be extremely
useful, but its provenance and confidence must remain explicit.

### 3.1 Provenance fields

At minimum, keep:

- source format;
- source identifier or file path;
- experimental or predicted status;
- method, if known;
- resolution or confidence measure, if known;
- the date or version of the source;
- any processing steps applied after import.

A coordinate file without provenance is difficult to reproduce or compare.

## 4. Structural data and uncertainty

A structure file is not a complete molecular truth. It is a model with gaps and
uncertainty.

| Feature | What it means | What not to do |
| --- | --- | --- |
| missing atom | the atom was not observed or not modelled | treat it as present at the origin |
| missing residue | a region is absent or disordered | silently rebuild it |
| alternate location | two or more modelled positions for an atom | average them without recording the choice |
| occupancy | relative occupancy of an alternate location | assume occupancies always sum to one |
| B-factor | displacement or uncertainty parameter | interpret it as a simulation temperature |
| water | a modelled solvent molecule | assume the list is complete or biologically unique |
| ion or cofactor | a non-polymer component | discard it before checking its role |
| multiple models | an ensemble, often from NMR | treat model 1 as the only truth |
| missing hydrogens | common in X-ray structures | assume the protonation state is known |

### 4.1 Alternate locations and occupancy

An atom can have more than one modelled position. The file may label these with
alternate-location identifiers such as `A` and `B`, and give each an occupancy.

Do not silently discard alternate locations during import. Either:

- preserve them;
- select one according to an explicit rule; or
- reject the structure and ask the user to decide.

Whichever choice is made, record it in the structure provenance.

### 4.2 Missing information

Missing atoms and residues are not zero coordinates. They are absent data. A
force field or docking calculation may need a complete structure, but completing
it is a modelling step with its own assumptions.

The import boundary should preserve:

- which atoms were present in the source;
- which residues were incomplete;
- which regions were absent;
- whether hydrogens were present;
- whether waters, ions, and cofactors were retained.

## 5. File formats and the import boundary

### 5.1 PDB and mmCIF

The PDB format is widely supported but uses fixed columns and has limitations
for large structures, multi-model data, and rich metadata. mmCIF is the modern
format used by the PDB archive and supports larger structures, assemblies, and
more explicit metadata.

Prefer mmCIF for new work. Support PDB where a tool or dataset requires it, but
do not let the limitations of the legacy format define the internal data model.

### 5.2 Angstrom to nanometre

Structural databases report coordinates in Angstrom. Moleculab uses nanometre.
The conversion is exact:

$$ 1\ \mathrm{Angstrom} = 0.1\ \mathrm{nm}. $$

Convert once, at the import boundary. Do not mix Angstrom coordinates with
Moleculab positions later in the pipeline.

### 5.3 The adapter boundary

The intended direction is:

```text
mmCIF / PDB / structure toolkit
    parsing, hierarchy, metadata, alternate locations
            ↓  adapter
Moleculab structure record
    atoms, residues, chains, positions in nm, provenance
            ↓
Track C analysis
    binding sites, contacts, fingerprints
            ↓
Track A physics
    topology, force field, simulation
```

The physics engine should not import a structure parser. The structure layer
should not silently rewrite chemistry.

A useful atom record contains:

```text
atom name, element, residue name, chain id, residue number,
insertion code, alternate location, position in nm,
occupancy, B-factor, model, hetero flag, water flag
```

A useful structure record also contains the source format, source identifier,
provenance, and any import decisions.

## 6. Inspecting a structure

The runnable example uses `gemmi`, an optional external structure toolkit:

```bash
uv run --with gemmi python examples/protein_structure.py
```

It reads a small synthetic PDB excerpt and reports the hierarchy, alternate
locations, missing information, and coordinate conversion:

```text
Structure summary
  source format: Pdb
  models: 1
  chains: A, B
  atom records: 9
  polymer residues: 2
  waters: 1
  other hetero residues: 1

Polymer residues
  ALA A 1: 6 atom records
  GLY A 3: 1 atom record
  residue-number gap in chain A: 1 -> 3

Alternate locations
  ALA A 1 CB: altloc '' occupancy 0.60
  ALA A 1 CB: altloc 'A' occupancy 0.40

Coordinate boundary
  first atom N: (11.104, 6.134, -6.504) Angstrom
  converted: (1.1104, 0.6134, -0.6504) nm

Provenance: experimental teaching excerpt (declared by the caller)
Missing atoms are left missing; the example does not rebuild them.
```

Several details are deliberate:

- ALA 1 has six atom records but only five distinct atom names, because CB has
  two alternate locations. Counting atom records is not the same as counting
  unique atoms.
- GLY 3 has only one atom record. The example reports the incomplete residue
  rather than filling in the missing atoms.
- The numbering gap from residue 1 to residue 3 is reported as a gap, not
  interpreted as a definitive missing residue without further evidence.
- The blank alternate-location marker from the toolkit is normalized to `""`
  for display.
- The first coordinate is converted from Angstrom to nanometre exactly once.

### 6.1 Reading a structure with a toolkit

A typical toolkit loop looks like this:

```python
import gemmi

structure = gemmi.read_structure("structure.cif")
model = structure[0]

for chain in model:
    for residue in chain:
        for atom in residue:
            position_nm = (
                atom.pos.x * 0.1,
                atom.pos.y * 0.1,
                atom.pos.z * 0.1,
            )
```

The toolkit supplies parsed hierarchy and metadata. The adapter supplies the
conversion to Moleculab units and the provenance record.

## 7. What a structure does not tell you

Even a high-quality structure does not directly tell you:

- the protonation state of every residue;
- the positions of missing loops or disordered regions;
- the equilibrium distribution of conformations;
- which waters are functionally important;
- whether a predicted pocket is biologically relevant;
- the binding affinity of a ligand;
- whether the crystal form matches the biological assembly.

Those are separate modelling questions. The structure is the starting point,
not the answer.

## 8. Acceptance criteria

A protein-structure import is correct when:

- Angstrom-to-nanometre conversion gives the expected coordinates;
- atom and residue counts are preserved;
- chain identifiers, residue numbers, and insertion codes remain
  distinguishable;
- alternate locations and occupancies are represented explicitly;
- missing atoms and residues remain missing rather than being silently
  interpreted as zero or as a complete structure;
- experimental, predicted, and unknown provenance remain distinguishable;
- the source format and identifier are retained;
- structure parsing does not leak toolkit objects into the physics layer.

## 9. Exercises

### Inspection

1. Run `examples/protein_structure.py` and reproduce the output.
2. Change the PDB excerpt to add a second model. How should the structure record
   represent both models?
3. Add an insertion code to a residue and explain why `chain + number` is no
   longer a unique residue identity.
4. Find a real mmCIF entry and inspect its chain identifiers, residue numbering,
   alternate locations, and missing residues.
5. Compare the same structure in PDB and mmCIF form. Which metadata is easier to
   preserve in each format?

### Reasoning

1. Why is a residue number not always a sequence index?
2. Why should alternate locations not be averaged silently?
3. Why is a predicted structure not interchangeable with an experimental one?
4. Why should geometry code avoid Angstrom coordinates?
5. Why is a crystallographic water not automatically a functional water?
6. Why is a structure file not a complete molecule?

### Implementation

1. Write a structure record with atoms, residues, chains, positions in nm,
   provenance, and source format.
2. Preserve the original alternate-location labels and occupancies.
3. Record which import decisions were made.
4. Do not fill missing atoms in the import layer.
5. Keep the structure adapter independent of the Track A force-field classes.

## 10. Summary

- Protein structure has a hierarchy: atoms, residues, chains, secondary
  structure, domains, and assemblies.
- Residue numbers and chain identifiers are labels, not global identities.
- Experimental and predicted structures have different provenance and
  uncertainty.
- Alternate locations, occupancies, missing atoms, waters, ions, and multiple
  models are data, not noise to discard silently.
- Convert Angstrom to nanometre once at the import boundary.
- Prefer mmCIF for new structure work and preserve metadata.
- A structure file is a model with gaps, not a complete molecular truth.

**Previous:** [Track C index](README.md).

**Next:** [02 — Binding-site geometry](02_binding_sites.md).
