# Track B, module 06 — Molecular polarity and electronic properties

> **Where this fits.** Track B modules 01–05 established structure,
> representation, properties, protonation, and stereochemistry. This module
> connects chemical electron distribution to the electrostatics used in Track A.
>
> **What you should already know.** Formal charge, protonation, molecular
> graphs, and the Coulomb potential.
>
> **What you will be able to do afterwards.** Distinguish formal charge, partial
> charge, and protonation state; compute a dipole moment; explain what the
> $q_i$ values in a force field do and do not represent.

## 1. Learning objectives

**Core**

- Explain electronegativity and bond polarity.
- Distinguish formal charge from partial atomic charge.
- Distinguish partial charge from protonation state.
- Compute a molecular dipole moment from charges and positions.
- Explain what a force-field charge represents.

**Deeper insight**

- A dipole moment depends on the coordinate origin for a charged system.
- Partial charges are model-dependent, not experimental observables.
- Molecular electrostatic potential is a spatial property, not a single number.
- Hydrogen-bond donor/acceptor status depends on chemistry, protonation, and
  geometry.
- Coulomb interactions use the chosen charge model; changing the model changes
  the electrostatics.

**Looking ahead**

- Track A's Coulomb, reaction-field, and Ewald models use these charges.
- Track C's interaction analysis depends on protonation and charge assignment.
- The capstone must state which charge model and protonation state were used.

## 2. Electronegativity and bond polarity

Electronegativity describes how strongly an atom attracts electrons in a bond.

A bond between atoms of different electronegativity is polar: electron density
is unevenly distributed. This is a qualitative chemical picture; the actual
electron distribution is a quantum-mechanical property.

Bond polarity affects:

- dipole moments;
- hydrogen-bond strength and direction;
- electrostatic interactions;
- molecular recognition;
- the choice of force-field parameters.

Do not confuse a polar bond with a formal charge. A neutral molecule can have
strong bond dipoles and a large molecular dipole moment.

## 3. Formal charge

Formal charge is an integer bookkeeping quantity from a Lewis structure.

For an atom,

$$ \text{formal charge} = \text{valence electrons} - \text{nonbonding electrons} - \tfrac{1}{2}\text{bonding electrons}. $$

Formal charges:

- are integers;
- depend on the chosen Lewis structure;
- are not the same as partial charges;
- change when protonation changes.

Example: acetate has a formal charge of −1, localized in the Lewis picture on
one oxygen, while the actual negative charge is delocalised over the carboxylate
group.

## 4. Partial charge

A partial atomic charge is a real-valued model quantity representing how much
electron density is assigned to an atom.

Partial charges:

- are model-dependent;
- are not directly measurable;
- depend on the molecular geometry and the charge model;
- are used by force fields to compute Coulomb interactions.

A common educational model is Gasteiger charge assignment. It is a fast
empirical scheme, not an experimental measurement and not a unique answer.

Different charge models can give different partial charges for the same
molecule. The force field's electrostatics inherit that choice.

## 5. Protonation state

Protonation state describes whether a proton is present at a particular site.

It is different from both formal and partial charge:

| Concept | Example | Nature |
| --- | --- | --- |
| protonation state | acetic acid versus acetate | discrete chemical state |
| formal charge | acetate has −1 | integer bookkeeping |
| partial charge | O has a model-dependent negative value | continuous model parameter |

Changing protonation changes:

- the molecular graph;
- the formal charge;
- the number of hydrogen-bond donors and acceptors;
- the partial-charge model input;
- the force-field topology.

Do not model a protonation change by merely scaling the partial charges of a
fixed graph.

## 6. Dipole moment

For charges $q_i$ at positions $\mathbf r_i$, the dipole moment about an
origin $\mathbf r_0$ is

$$ \boldsymbol\mu = \sum_i q_i(\mathbf r_i-\mathbf r_0). $$

The magnitude is $|\boldsymbol\mu|$.

For a neutral system, the result is independent of the chosen origin. For a
charged system, it is origin-dependent, so an explicit origin is required.

The implementation stores the dipole in e nm and converts to Debye:

$$ 1\ \mathrm{e\ nm}\approx48.032\ \mathrm{Debye}. $$

A dipole moment is a vector. Two molecules can have similar magnitudes but
different orientations.

### 6.1 What a dipole moment is not

A dipole moment:

- is not a formal charge;
- is not a partial-charge model;
- is not a full electrostatic potential;
- depends on geometry and charge assignment;
- does not by itself determine binding affinity.

## 7. Molecular electrostatic potential

The molecular electrostatic potential (MEP) describes the electrostatic
potential at points around a molecule.

It is a spatial field, not a single number. It depends on:

- the charge distribution;
- the molecular conformation;
- the protonation state;
- the surrounding environment.

A colour-coded MEP surface is a visualization of a model, not a direct
experimental measurement.

MEP is useful for:

- identifying likely hydrogen-bond donor and acceptor regions;
- comparing electrostatic complementarity;
- generating pharmacophore hypotheses;
- explaining why a charged group prefers a particular environment.

## 8. Hydrogen-bond donors and acceptors

Hydrogen-bonding capacity is chemical, not merely electrostatic.

A donor needs:

- a hydrogen attached to an electronegative atom;
- a suitable acceptor;
- a plausible geometry.

An acceptor needs:

- a lone pair or suitable electron density;
- a compatible geometry.

Protonation changes donor/acceptor status:

- a neutral carboxylic acid is a donor;
- a carboxylate is an acceptor and a salt-bridge partner;
- an ammonium group is a donor;
- a neutral amine can donate and accept.

This is why Track B module 04 is a prerequisite for interaction analysis.

## 9. Force-field charges and Coulomb interactions

Track A's Coulomb term is

$$ U_{ij}(r)=k_e\frac{q_iq_j}{r}, $$

where $q_i$ and $q_j$ are the charges in elementary-charge units and $r$
is in nm.

The $q_i$ values are force-field parameters. They are chosen to reproduce
selected properties, not measured directly.

Consequences:

- changing the charge model changes the electrostatics;
- changing protonation changes the topology and charges;
- a partial charge is not a formal charge;
- a force-field charge is not an experimental observable.

A force field can be physically consistent while still using an approximate
charge model.

## 10. RDKit example and validation

The runnable example is:

```bash
uv run --with rdkit python examples/chemistry/molecular_polarity.py
```

It demonstrates:

- formal charges for acetic acid, acetate, methylamine, and methylammonium;
- per-atom formal charges;
- Gasteiger partial charges for water;
- a dipole moment computed from those partial charges and coordinates;
- the connection to the Coulomb expression.

Observed output:

```text
Formal charge and protonation state
   acetic acid: net=+0 atom formal charges=[0, 0, 0, 0]
       acetate: net=-1 atom formal charges=[0, 0, 0, -1]
   methylamine: net=+0 atom formal charges=[0, 0]
  methylammonium: net=+1 atom formal charges=[0, 1]

Water: Gasteiger partial charges
   O: partial charge = -0.412 e
   H: partial charge = +0.206 e
   H: partial charge = +0.206 e
  net partial charge: +0.000000 e
  dipole magnitude: 0.0246 e nm
  dipole magnitude: 1.179 Debye
```

The Gasteiger dipole is a model result. It is not an experimental dipole and
will change with the charge model and conformer.

The dipole calculation is implemented in
`src/moleculab/analysis/electrostatics.py` and tested in
`tests/models/test_electrostatics.py`.

### 10.1 Acceptance criteria

A molecular-polarity module is correct when:

- formal charge is treated as integer bookkeeping;
- partial charge is treated as a model-dependent quantity;
- protonation state is not conflated with either;
- dipole moments use explicit charges and positions;
- charged-system dipole calculations require an explicit origin;
- units are consistent and Debye conversion is documented;
- the connection to the Coulomb term is explicit;
- no partial-charge model is presented as experimentally unique.

## 11. Common misconceptions

- **Formal charge is not partial charge.**
- **Partial charge is not an experimental observable.**
- **Protonation state is not a continuous charge scaling.**
- **A neutral molecule can have a large dipole moment.**
- **A dipole moment is a vector, not just a magnitude.**
- **A charged molecule's dipole moment depends on the origin.**
- **Gasteiger charges are not the unique correct charges.**
- **An MEP surface is a model visualization, not a measurement.**

## 12. Exercises

### Reasoning

1. Why is formal charge an integer while partial charge is not?
2. Why does a charged molecule's dipole moment depend on the origin?
3. Why can a neutral molecule have a large dipole moment?
4. Why is protonation state not the same as formal charge?
5. Why is Gasteiger charge not an experimental result?
6. Why does the choice of partial charges affect Coulomb energy?
7. Why does a carboxylate behave differently from a carboxylic acid in a
   hydrogen-bond analysis?

### Calculation

1. Two charges $+1e$ and $-1e$ are separated by 0.2 nm. Compute the dipole
   magnitude in e nm and Debye.
2. Repeat with charges $+0.5e$ and $-0.5e$.
3. For a neutral system, shift all coordinates by 1 nm and confirm the dipole is
   unchanged.
4. For a charged system, compute the dipole about two different origins and
   explain why the results differ.
5. Compute the Coulomb energy for two unit charges at 0.5 nm.

### Implementation

1. Run `examples/chemistry/molecular_polarity.py` and reproduce the output.
2. Change the protonation state of acetic acid and compare formal charges.
3. Compare the Gasteiger charges of water and ammonia.
4. Embed and optimise water with two seeds and compare the dipole magnitudes.
5. Explain why a different charge model would change the Coulomb energy.

## 13. Summary

- Electronegativity and bond polarity shape the electron distribution.
- Formal charge is integer bookkeeping from a Lewis structure.
- Partial charge is a model-dependent real-valued assignment.
- Protonation state is a discrete chemical state.
- A dipole moment is a charge-weighted displacement vector.
- Neutral dipoles are origin-independent; charged dipoles are not.
- Molecular electrostatic potential is a spatial field, not a single number.
- Hydrogen-bond donor/acceptor status depends on chemistry and protonation.
- Force-field charges are parameters, and Coulomb interactions inherit them.
- No charge model is experimentally unique; state the model and provenance.

**Previous:** [05 — Stereochemistry and molecular recognition](05_stereochemistry.md).

**Next:** [07 — Hydrogen bonding and molecular recognition](07_hydrogen_bonding.md).
