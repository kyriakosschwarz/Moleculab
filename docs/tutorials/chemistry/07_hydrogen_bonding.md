# Track B, module 07 — Hydrogen bonding and molecular recognition

> **Where this fits.** [Module 06](06_molecular_polarity.md) explained charge
> distributions and dipoles. This module turns that picture into a directional
> chemical interaction: the hydrogen bond.
>
> **What you should already know.** Formal and partial charge, protonation,
> electronegativity, and the Coulomb interaction.
>
> **What you will be able to do afterwards.** Identify donors and acceptors,
> interpret hydrogen-bond geometry, explain solvent competition and
> desolvation, and connect chemical hydrogen bonding to Track C's geometric
> detector.

## 1. Learning objectives

**Core**

- Identify hydrogen-bond donors and acceptors from chemical structure.
- Explain why electronegativity and lone-pair availability matter.
- Interpret donor–acceptor distance and D-H...A angle.
- Distinguish intramolecular from intermolecular hydrogen bonds.
- Explain why a hydrogen bond is not automatically favourable to binding.

**Deeper insight**

- Donor/acceptor status depends on protonation and tautomerism.
- Hydrogen-bond criteria are operational definitions, not physical laws.
- Solvent competes for the same donor and acceptor groups.
- Fixed-charge force fields usually represent hydrogen bonding through
  electrostatics and van der Waals terms rather than a separate bond term.
- A persistent hydrogen bond can coexist with an unfavourable binding free
  energy.

**Looking ahead**

- Track C's interaction detector uses the same geometric criteria.
- The capstone must distinguish geometric hydrogen bonds from thermodynamic
  contributions.

## 2. What a hydrogen bond is

A hydrogen bond is usually written D-H...A:

- **D** is the donor heavy atom;
- **H** is a hydrogen attached to D;
- **A** is the acceptor.

The donor–acceptor distance is

$$ d_{DA}=\lVert\mathbf r_D-\mathbf r_A\rVert. $$

The D-H...A angle at the hydrogen is

$$ \theta_{DHA}=\cos^{-1}\left(\frac{(\mathbf r_D-\mathbf r_H)\cdot(\mathbf r_A-\mathbf r_H)}{\lVert\mathbf r_D-\mathbf r_H\rVert\ \lVert\mathbf r_A-\mathbf r_H\rVert}\right). $$

A plausible hydrogen bond generally has:

- a suitable donor and acceptor;
- a short donor–acceptor distance;
- a reasonably linear D-H...A arrangement.

The thresholds are conventions. They are not universal physical laws.

## 3. Donors and acceptors are chemical features

Do not infer donor/acceptor status from the element alone.

| Group | Donor? | Acceptor? | Note |
| --- | --- | --- | --- |
| alcohol O-H | yes | yes | oxygen has lone pairs |
| ether O | no | yes | no O-H hydrogen |
| carboxylic acid O-H | yes | yes | neutral acid |
| carboxylate O | no | yes | deprotonated, often delocalised |
| amine N-H | yes | yes | neutral amine |
| ammonium N-H | yes | no | positive charge, no lone pair |
| amide N-H | yes | weak/no | lone pair delocalised into carbonyl |
| amide carbonyl O | no | yes | strong acceptor |

Protonation and tautomerism can change the table. This is why Track B module 04
comes before hydrogen-bond analysis.

RDKit can count donors and acceptors:

```text
ethanol:          HBD=1 HBA=1
acetic acid:      HBD=1 HBA=1
acetate:          HBD=0 HBA=2
methylamine:      HBD=1 HBA=1
methylammonium:   HBD=1 HBA=0
acetamide:        HBD=1 HBA=1
```

These counts are cheminformatics features. They do not by themselves tell you
whether a hydrogen bond is geometrically possible or favourable.

## 4. Geometry and directionality

Hydrogen bonds are directional. A short donor–acceptor distance with a poor
angle is not a good hydrogen bond.

The example constructs:

- a linear D-H...A contact at 0.30 nm and 180°;
- a short contact with an unsuitable angle;
- a donor with no hydrogen positions.

Only the linear case is detected:

```text
linear donor-acceptor contact:
  distance=0.300 nm angle=180.0 degrees
short contact with unsuitable angle: 0 H-bond(s)
donor without hydrogen positions: 0 H-bond(s)
```

The detector uses `HydrogenBondCriteria` from Track C:

| Criterion | Default | Comparison |
| --- | --- | --- |
| donor–acceptor distance | 0.35 nm | inclusive |
| D-H...A angle | 120° | inclusive |

Changing these thresholds changes the detected interaction set. Report the
criteria with the result.

## 5. Directionality and competing interactions

A donor can point toward more than one acceptor. An acceptor can receive more
than one donor. The preferred arrangement depends on:

- lone-pair directionality;
- steric accessibility;
- nearby charges;
- solvent competition;
- cooperative hydrogen-bond networks.

Hydrogen bonds are not isolated pairwise objects in a condensed phase. They
are part of a network that can reorganise when a ligand binds.

## 6. Intramolecular versus intermolecular hydrogen bonds

An **intramolecular** hydrogen bond forms between two groups in the same
molecule. It can stabilise a conformation, but it may also:

- reduce the number of groups available to the protein;
- be solvent-exposed;
- create a strained ring;
- change the conformational ensemble.

An **intermolecular** hydrogen bond forms between different molecules, such as
a ligand and a protein, or a solute and water.

A hydrogen-bond count that mixes both types can be misleading. State whether
the count is intra- or intermolecular.

## 7. Hydrogen bonds in water and solvent

Water is not a passive background. It donates and accepts hydrogen bonds
competitively.

When a ligand binds, it may need to:

- break ligand–water hydrogen bonds;
- break protein–water hydrogen bonds;
- displace or reorganise waters;
- form new ligand–protein hydrogen bonds.

The net free-energy change is a balance:

$$ \Delta G_{\mathrm{bind}} = \Delta H - T\Delta S, $$

with contributions from interactions, desolvation, conformational entropy, and
water reorganisation.

A ligand–protein hydrogen bond that replaces a ligand–water hydrogen bond may
provide little net stabilisation. A buried charged group may pay a large
desolvation penalty.

## 8. Why hydrogen-bond count is not binding affinity

A persistent hydrogen bond can be:

- favourable;
- roughly neutral after desolvation;
- part of a strained binding mode;
- less important than a water-mediated interaction;
- offset by loss of conformational entropy.

Counting hydrogen bonds is a hypothesis-generating exercise. It is not a free
energy calculation.

This connects directly to [Track C module 10](../drug_discovery/10_structural_observations_to_thermodynamics.md).

## 9. Protein–ligand hydrogen bonds

In a protein–ligand complex, hydrogen bonds can be:

- direct;
- water-mediated;
- part of a charged salt bridge;
- weak or bifurcated;
- satisfied by a different donor/acceptor after a protonation change.

A useful interaction analysis records:

- donor and acceptor identity;
- distance and angle;
- whether the interaction is direct or water-mediated;
- the protonation and tautomer assumptions;
- the detection criteria.

Track C's `HydrogenBond` records provide exactly this kind of structured
observation.

## 10. How force fields represent hydrogen bonding

In a fixed-charge force field, a hydrogen bond is usually not a separate
bonded term. It emerges from:

- Coulomb interactions between partial charges;
- Lennard-Jones interactions;
- the molecular geometry that places donor and acceptor correctly.

Some force fields add explicit hydrogen-bond or polarisation terms, but a
simple fixed-charge model does not.

Consequences:

- partial charges affect hydrogen-bond strength;
- protonation changes donor/acceptor status and charges;
- the force field cannot create a hydrogen bond where the geometry is poor;
- a geometric detector and a force-field energy answer different questions.

## 11. RDKit and Track C connection

The runnable example is:

```bash
uv run --with rdkit python examples/hydrogen_bonding.py
```

It uses RDKit for donor/acceptor counts and Track C's
`detect_hydrogen_bonds` for the geometric criteria.

Observed output:

```text
Donor and acceptor counts (RDKit)
          ethanol: HBD=1 HBA=1
      acetic acid: HBD=1 HBA=1
          acetate: HBD=0 HBA=2
      methylamine: HBD=1 HBA=1
   methylammonium: HBD=1 HBA=0
        acetamide: HBD=1 HBA=1

Geometric criteria
  linear donor-acceptor contact:
    distance=0.300 nm angle=180.0 degrees
  short contact with unsuitable angle: 0 H-bond(s)
  donor without hydrogen positions: 0 H-bond(s)
```

The chemistry layer explains what the donor/acceptor flags mean. Track C's
geometry layer evaluates whether the arrangement satisfies the criteria. They
are separate responsibilities.

## 12. Common misconceptions

- **Every N or O is an acceptor.** Amide nitrogen is usually a poor acceptor;
  protonation can remove acceptor capacity.
- **A short distance is enough.** Directionality matters.
- **A hydrogen bond is always favourable.** Desolvation and entropy can offset
  it.
- **Hydrogen-bond count predicts affinity.** It does not.
- **Intramolecular and intermolecular hydrogen bonds are interchangeable in a
  count.** They have different consequences.
- **A force field has an explicit hydrogen-bond term.** In a fixed-charge model,
  it usually emerges from electrostatics and van der Waals terms.
- **A geometric H-bond is a thermodynamic contribution.** It is a structural
  observation.

## 13. Acceptance criteria

A hydrogen-bond chemistry module is correct when:

- donors and acceptors are treated as chemical features;
- protonation-dependent changes are explicit;
- distance and angle criteria are documented;
- short-distance/poor-angle contacts are rejected;
- intramolecular and intermolecular hydrogen bonds are distinguished;
- solvent competition and desolvation are discussed;
- the connection to force-field electrostatics is explicit;
- geometric detection is not presented as a free-energy calculation.

## 14. Exercises

### Inspection

1. Classify the donor/acceptor status of ethanol, acetic acid, acetate,
   methylamine, methylammonium, and acetamide.
2. Draw a linear and a bent D-H...A arrangement and predict whether each passes
   the default criteria.
3. Identify an intramolecular hydrogen bond in a simple molecule.
4. Find a protein–ligand hydrogen bond that could be water-mediated.
5. Compare the hydrogen-bond donor/acceptor counts of a neutral and protonated
   amine.

### Reasoning

1. Why is an amide nitrogen a poor acceptor?
2. Why is a short donor–acceptor distance insufficient?
3. Why can a hydrogen bond be neutral in net free-energy terms?
4. Why does solvent competition matter?
5. Why are intramolecular hydrogen bonds not equivalent to intermolecular ones?
6. Why is a fixed-charge force field's hydrogen bonding implicit?
7. Why is a geometric hydrogen-bond count not an affinity prediction?

### Implementation

1. Run `examples/hydrogen_bonding.py` and reproduce the output.
2. Change the distance cutoff and observe which contacts survive.
3. Change the angle cutoff and observe which contacts survive.
4. Add a donor without hydrogen positions and confirm that it is skipped.
5. Connect the example to a Track C interaction fingerprint and explain what
   additional information the fingerprint preserves.

## 15. Summary

- Hydrogen bonds are directional donor–hydrogen...acceptor interactions.
- Donor/acceptor status is chemical and protonation-dependent, not an element
  lookup.
- Distance and angle criteria are operational definitions.
- Intramolecular and intermolecular hydrogen bonds have different consequences.
- Water competes for donors and acceptors, so desolvation matters.
- A hydrogen-bond count is not a binding free energy.
- Fixed-charge force fields usually represent hydrogen bonding through
  electrostatics and van der Waals interactions.
- Chemistry explains the features; Track C geometry evaluates the arrangement.

**Previous:** [06 — Molecular polarity and electronic properties](06_molecular_polarity.md).

**Next:** Medicinal chemistry and ADME fundamentals (module 08, planned).
