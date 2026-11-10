# Chemical structure and conformational behaviour

> **Where this fits.** This is the first module in
> [Track B — chemical structure and properties](README.md). It follows the
> physics of bonds, angles, and torsions (Track A parts 04b and 07) and asks a
> chemical question: why do real molecules prefer particular shapes?
>
> **What you should already know.** Atoms, bonds, angles, torsions, and the
> periodic torsion potential.
>
> **What you will be able to do afterwards.** Recognise common functional
> groups, explain how rings and stereochemistry restrict conformations, read a
> torsional energy curve, and distinguish a force-field torsion term from the
> electronic effects it represents.

## 1. Learning objectives

**Core**

- Recognise the common functional groups that shape molecular interactions.
- Explain why rings restrict rotation.
- Read a torsional energy curve with minima and barriers.
- Relate a force-field torsion term to a conformational preference.

**Deeper insight**

- Explain why aromaticity and conjugation are electronic effects, not torsion
  terms.
- Explain how conformational restriction connects to entropy and binding.
- Explain why a single threefold term is a teaching model, not a quantitative
  butane force field.

**Looking ahead**

- pKa and protonation change the charge state and hydrogen-bonding pattern.
- Interaction fingerprints describe which residues contact which ligand groups.
- Desolvation and water networks connect conformation to binding free energy.

## 2. Functional groups

A functional group is a small pattern of atoms whose chemistry is recognisable
across many molecules. The goal is not a catalogue to memorise. The goal is to
connect structure to interactions and energy.

| Group | Structural feature | Consequence |
| --- | --- | --- |
| hydroxyl | O-H | hydrogen-bond donor and acceptor, polar |
| carbonyl | C=O | strong dipole, hydrogen-bond acceptor |
| carboxylic acid | COOH | ionisable, hydrogen bonding, pH-sensitive |
| amide | N-C=O | planar, restricted rotation, hydrogen bonding |
| amine | N | basic, can be protonated |
| aromatic ring | conjugated ring | planar, hydrophobic and polarisable |
| thiol | S-H | sulfur chemistry, oxidation-sensitive |
| ether | C-O-C | polar but no O-H donor |
| ester | C(=O)-O-C | polar, hydrolytically labile |

The force field does not “know” that a group is called an amide. It sees atom
types, bonds, angles, torsions, partial charges, and Lennard-Jones parameters.
The chemical name is a human interpretation of that graph.

## 3. Rings and restricted rotation

A ring is a closed path of bonded atoms. It changes the available
conformations.

- **Benzene** is planar and aromatic. The planarity and delocalised electron
  density are electronic properties; a force field represents their effects
  through atom types, partial charges, bonded parameters, and planarity
  restraints.
- **Cyclohexane** has chair, boat, and other conformations. The chair is
  normally the lowest-energy form; the boat is higher because of eclipsing and
  steric strain.
- **Piperidine** and **piperazine** are saturated six-membered rings containing
  nitrogen. They appear often in drug-like molecules and have their own
  conformational preferences.
- **Cyclopropane** and **cyclobutane** have large ring strain because their
  bond angles deviate strongly from the preferred tetrahedral geometry.

Typical ring-strain magnitudes are on the order of:

| Ring | Approximate strain (kJ/mol) |
| --- | --- |
| cyclopropane | about 115 |
| cyclobutane | about 110 |
| cyclopentane | about 25 |
| cyclohexane chair | about 0 |

These are order-of-magnitude references. The exact value depends on the force
field and on substituents.

A ring restricts the number of rotatable bonds. That reduces the conformational
space and therefore the conformational entropy. It can also pre-organise a
ligand for binding, which may help or hurt depending on the target.

## 4. Stereochemistry

Stereoisomers have the same connectivity but different three-dimensional
arrangement.

- **Enantiomers** are mirror images.
- **Diastereomers** are not mirror images.
- A **chiral centre** has four different substituents.
- **cis/trans** isomerism describes relative arrangement across a ring or a
  double bond.

Different stereoisomers can have very different biological activity. A docking
or simulation workflow must use the correct stereochemistry; it is not a
cosmetic detail.

## 5. Conformations: ethane and butane

### 5.1 Ethane

Ethane has a threefold rotational barrier. A common teaching form is

$$ U(\phi)=\frac{V_3}{2}\left[1+\cos(3\phi)\right]. $$

The barrier is about 12 kJ/mol for ethane. The three staggered conformations
are equivalent, and the three eclipsed conformations are the maxima.

### 5.2 Butane

Butane has a more interesting torsional profile:

- **anti** conformation: the two methyl groups are opposite; usually the
  lowest energy;
- **gauche** conformations: the methyl groups are about 60° apart; local
  minima;
- **eclipsed** conformations: higher energy.

A single threefold term gives three equivalent minima. A real butane force
field adds more terms to reproduce the anti/gauche energy difference and the
different eclipsed barriers. This is an important lesson:

> The periodic torsion term is the right **shape** for rotational periodicity,
> but a quantitative force field may need several terms.

## 6. A computational torsion model

Moleculab's `PeriodicTorsionPotential` implements

$$ U(\phi)=k\left[1+\cos(n\phi-\delta)\right]. $$

The energy ranges from 0 to $2k$. For $n=3$ there are three minima per full
rotation.

The example builds a four-atom chain with harmonic bonds and angles and one
periodic torsion. It starts near an eclipsed geometry, minimises the energy,
and relaxes to a staggered dihedral.

```bash
uv run python examples/chemical_structure.py
```

Observed output:

```text
initial dihedral:   90.00 deg
initial energy:     53.959 kJ/mol
final dihedral:     60.00 deg
final energy:        0.000 kJ/mol
accepted steps:      1433

phi =   0.0 deg: U = 4.00 kJ/mol
phi =  60.0 deg: U = 0.00 kJ/mol
phi = 120.0 deg: U = 4.00 kJ/mol
phi = 180.0 deg: U = 0.00 kJ/mol

Full rotational scan
  minimum: 0.00 kJ/mol
  maximum: 4.00 kJ/mol
  barrier: 4.00 kJ/mol
  U(0) = U(360): 4.00 = 4.00 kJ/mol

Conformational entropy: anti/gauche toy model at 300 K
  state 0: E = 0.0 kJ/mol, p = 0.418
  state 1: E = 0.9 kJ/mol, p = 0.291
  state 2: E = 0.9 kJ/mol, p = 0.291
  conformational entropy = 9.01 J/(mol K)
```

The phase convention matters. The implementation uses

$$ U(\phi)=k\left[1+\cos(n\phi-\delta)\right]. $$

For the example, $k=2$ kJ/mol, $n=3$, and $\delta=0$. Here $k$ is the
**energy amplitude** of the torsion term, not a harmonic force constant. Its
unit is kJ/mol. The minima are at 60°, 180°, and 300°, where
$\cos(3\phi)=-1$. The maxima are at 0°, 120°, and 240°. The energy ranges from
0 to $2k=4$ kJ/mol. The minimum being zero is an energy reference choice, not a
claim that the conformation has zero absolute energy.

The full scan checks the physical periodic property

$$ U(\phi)=U(\phi+360^\circ). $$

For this threefold term, the stronger statement is

$$ U(\phi)=U(\phi+120^\circ). $$

The tests in `tests/test_torsion.py` check periodicity, the three minima, the
three barriers, and the barrier height.

The generalized force conjugate to the dihedral coordinate is

$$ F_\phi=-\frac{dU}{d\phi}=kn\sin(n\phi-\delta). $$

This is a force in kJ/(mol rad), not a Cartesian force in kJ/(mol nm). The
tests compare it with a central finite difference at asymmetric angles, where
the derivative is not trivially zero.

### 6.1 Why the starting geometry matters

A minimiser finds a nearby local minimum. Start in a different well and it may
relax to a different rotamer. This is the same local-minimum issue as the
double-well tutorial, now with a chemical interpretation.

### 6.2 Conformational entropy

A molecule with several low-energy rotamers has more conformational entropy
than a rigid molecule. For a simplified set of states with energies $E_i$, the
Boltzmann probability is

$$ p_i=\frac{e^{-E_i/(RT)}}{\sum_j e^{-E_j/(RT)}}. $$

The molar conformational entropy is

$$ S=-R\sum_i p_i\ln p_i. $$

For a toy butane model with one anti state at 0 kJ/mol and two gauche states at
0.9 kJ/mol, at 300 K:

```text
state 0: E = 0.0 kJ/mol, p = 0.418
state 1: E = 0.9 kJ/mol, p = 0.291
state 2: E = 0.9 kJ/mol, p = 0.291
conformational entropy = 9.01 J/(mol K)
```

For comparison, three equally populated states would give $R\ln 3\approx9.13$
J/(mol K). The lower value reflects the anti preference.

This is a **discrete three-state toy model**, not a quantitative prediction of
butane's conformational entropy. Real conformational basins occupy continuous
regions of configuration space, and a basin's statistical weight depends on its
energy and the volume of configurations it contains. The continuous
configurational partition function is

$$ Z_{\mathrm{conf}}=\int e^{-U(\mathbf R)/(RT)}d\mathbf R. $$

The toy model is a bridge to that idea, not a replacement for it.

The entropy result is also the bridge to binding thermodynamics: a flexible
ligand can lose conformational entropy when it binds. The lowest-energy
conformer alone does not determine the population or the free energy.

### 6.3 The ideal model versus real butane

The threefold term is an idealised model. Real butane has:

- one anti minimum;
- two gauche minima at higher energy than anti;
- different eclipsed barriers;
- steric and electronic contributions beyond a single cosine.

A real force field combines several torsional and nonbonded terms. The tutorial
uses the simple term because it makes the periodicity, minima, and barriers
visible. It should not be read as a quantitative butane potential.

## 7. What a force field does and does not represent

A classical force field represents:

- bond stretching;
- angle bending;
- torsional preferences;
- van der Waals interactions;
- electrostatics;
- often planarity and improper terms.

It does **not** directly represent:

- aromaticity as a quantum phenomenon;
- conjugation;
- charge transfer;
- polarisation;
- bond breaking;
- electronic excited states.

Those effects enter only through parameters fitted to experiment or quantum
chemistry, or through a QM/MM model. A torsion term can reproduce the shape of
a rotational barrier without explaining why the barrier exists.

## 8. Water and the environment

Conformation is not only an internal property. Water changes the relative
stability of polar and nonpolar groups. A folded conformation that buries a
charged group may be favourable in vacuum but costly in water. Desolvation is
therefore part of the conformational and binding story.

This is why Track B connects to the water, solvent, and free-energy modules of
Track A. A gas-phase torsion scan is a starting point, not the final answer.

## 9. Common misconceptions

1. **“A force field knows what an aromatic ring is.”** It sees atom types and
   parameters. Aromaticity is an electronic interpretation.
2. **“A single torsion term is a quantitative butane model.”** It captures
   periodicity; real force fields combine terms.
3. **“The lowest-energy conformer is the bioactive conformer.”** Binding can
   select a higher-energy conformation if the interaction energy compensates.
4. **“Rings are just cosmetic.”** They restrict conformations and change
   entropy, strain, and presentation.
5. **“Stereochemistry does not matter for docking.”** It can change the entire
   shape and interaction pattern.
6. **“Minimisation finds the global conformer.”** It finds a nearby local
   minimum.

## 10. Exercises

### Calculation

1. For $k=2$ kJ/mol and $n=3$, compute the torsion energy at 0°, 60°, 120°, and
   180°.
2. For $k=6$ kJ/mol, what is the full barrier height $2k$?
3. If a molecule has four independent rotamers, how many conformations are
   possible in a simple discrete model with three states per rotamer?

### Reasoning

1. Why is benzene planar while cyclohexane is not?
2. Why can a ring reduce conformational entropy?
3. Why does a force field need parameters for aromaticity rather than an
   explicit aromaticity term?
4. Why might the bioactive conformation be higher in energy than the global
   minimum?
5. Why is a single torsion term insufficient for quantitative butane?

### Implementation

1. Run `examples/chemical_structure.py` and change the initial geometry.
2. Scan the torsion energy from -180° to 180° and plot it.
3. Add a second torsion term and observe the new minima and barriers.
4. Minimise the chain from several starting dihedrals and compare the final
   rotamers.
5. Build a four-atom chain with a ring-like constraint and discuss which
   rotations are removed.

## 11. Summary

- Functional groups connect molecular structure to interactions.
- Rings restrict rotation and introduce strain and conformational entropy.
- Stereochemistry changes shape and biological activity.
- Torsional potentials explain why rotamers and barriers exist.
- A force field represents electronic effects through parameters, not by
  simulating electrons.
- Water and desolvation change the relative stability of conformations.
- The first chemistry module is a bridge from molecular mechanics to medicinal
  chemistry, not a replacement for the physics sequence.

**Previous:** [04b — From springs to chemistry](../physics/04b_chemical_bridge.md).

**Next:** [02 — Molecular representation](02_molecular_representation.md).
