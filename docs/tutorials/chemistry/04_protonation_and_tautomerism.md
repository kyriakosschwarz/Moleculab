# Protonation and tautomerism: chemical states at a given pH

> **Where this fits.** Modules 02 and 03 of
> [Track B](README.md) established molecular identity, representation,
> descriptors, and similarity. This module adds the chemical state: which
> protons are attached, which charges are present, and how those populations
> change with pH.
>
> **What you should already know.** Formal charge, molecular graphs, SMILES,
> and the Boltzmann distribution from the conformational-entropy discussion.
>
> **What you will be able to do afterwards.** Derive the Henderson-Hasselbalch
> equation, compute protonation fractions and multi-site microstate
> populations, distinguish tautomers from resonance structures, and explain
> why docking and molecular dynamics need an explicit chemical-state choice.

## 1. Learning objectives

**Core**

- Write acid dissociation equilibria and define $K_a$ and $\mathrm{p}K_a$.
- Derive the Henderson-Hasselbalch equation.
- Compute protonated and deprotonated fractions as functions of pH.
- Compute populations of several independent ionizable sites.
- Distinguish tautomers from resonance structures.

**Deeper insight**

- pH changes the populations of discrete protonation states; it does not
  continuously rescale the partial charges of one fixed graph.
- A pKa is environment-dependent, not a universal molecular constant.
- Enumerating states, estimating populations, and predicting pKa are three
  different computational problems.
- Tautomerisation changes connectivity and hydrogen-bonding patterns, not just
  the position of a proton in a drawing.

**Looking ahead**

- Protonation changes formal charge, partial charges, atom types, and force-field
  parameters.
- Docking and standard MD require a chosen protonation and tautomer state.
- Protein-ligand interactions depend on the charge state of both partners.

## 2. Acids, bases, and proton transfer

A Brønsted acid donates a proton and a Brønsted base accepts one. For a
monoprotic acid HA in water,

$$ \mathrm{HA} + \mathrm{H_2O} \rightleftharpoons \mathrm{A^-} + \mathrm{H_3O^+}. $$

The acid dissociation constant is

$$ K_a = \frac{[\mathrm{H^+}][\mathrm{A^-}]}{[\mathrm{HA}]}, $$

where the concentrations are approximate activities in dilute solution. pH is

$$ \mathrm{pH} = -\log_{10} a_{\mathrm{H^+}} \approx -\log_{10}[\mathrm{H^+}]. $$

A strong acid has a large $K_a$; a weak acid has a small one. The conjugate base
of HA is $\mathrm{A^-}$. For a base B, the conjugate acid is $\mathrm{BH^+}$:

$$ \mathrm{BH^+} \rightleftharpoons \mathrm{B} + \mathrm{H^+}. $$

The same equilibrium constant describes both cases. What changes is the label
we attach to the protonated and deprotonated forms.

Physiological pH is about 7.4, but it is not the only relevant condition.
Lysosomes are more acidic, the stomach is strongly acidic, and a protein
binding site can have a local electrostatic environment very different from
bulk water. The chemical state must be chosen for the condition being modelled.

## 3. pKa and the Henderson-Hasselbalch equation

By definition,

$$ \mathrm{p}K_a = -\log_{10} K_a. $$

Rearranging the equilibrium expression gives the ratio

$$ \frac{[\mathrm{A^-}]}{[\mathrm{HA}]} = \frac{K_a}{[\mathrm{H^+}]} = 10^{\mathrm{pH}-\mathrm{p}K_a}. $$

Taking the base-10 logarithm gives the Henderson-Hasselbalch equation:

$$ \mathrm{pH} = \mathrm{p}K_a + \log_{10}\left(\frac{[\mathrm{A^-}]}{[\mathrm{HA}]}\right). $$

The fractions follow from normalisation. For an acid,

$$ f_{\mathrm{A^-}} = \frac{1}{1+10^{\mathrm{p}K_a-\mathrm{pH}}}, \qquad f_{\mathrm{HA}} = \frac{1}{1+10^{\mathrm{pH}-\mathrm{p}K_a}}. $$

For a base, the fraction of conjugate acid is

$$ f_{\mathrm{BH^+}} = \frac{1}{1+10^{\mathrm{pH}-\mathrm{p}K_a}}, \qquad f_{\mathrm{B}} = \frac{1}{1+10^{\mathrm{p}K_a-\mathrm{pH}}}. $$

The labels differ, but the mathematics is the same. At $\mathrm{pH}=\mathrm{p}K_a$
the two forms are equally populated. One pH unit below the pKa gives a 10:1
ratio in favour of the protonated form; one unit above gives 10:1 in favour of
the deprotonated form.

### 3.1 Worked example: acetic acid

Acetic acid has $\mathrm{p}K_a \approx 4.76$. At pH 7.4,

$$ \frac{[\mathrm{A^-}]}{[\mathrm{HA}]} = 10^{7.4-4.76} = 10^{2.64} \approx 437. $$

Therefore

$$ f_{\mathrm{A^-}} = \frac{437}{438} \approx 0.9977. $$

Almost all of the acetic acid is deprotonated at physiological pH. At pH 4.76
the fraction is exactly 0.5.

## 4. Protonation-state populations as a function of pH

The important conceptual point is that pH selects among **discrete chemical
states**. It does not turn a fixed molecule into a continuously charged object.
For acetic acid the two states are neutral HA and anionic $\mathrm{A^-}$. For an amine
they are neutral B and cationic $\mathrm{BH^+}$.

At pH 7.4 the example in this module gives:

```text
       species   pKa  deprotonated  protonated
   acetic acid  4.76        0.9977      0.0023
      ammonium  9.25        0.0139      0.9861
     imidazole  6.00        0.9617      0.0383
```

For ammonium, "deprotonated" means neutral $\mathrm{NH_3}$ and "protonated" means
$\mathrm{NH_4^+}$. For imidazole, the neutral form is the deprotonated form and the
cation is the protonated form.

A titration table makes the steepness visible:

```text
Acetic acid titration (pKa 4.76)
   pH    f(A-)    f(HA)
 2.00   0.0017   0.9983
 3.76   0.0909   0.9091
 4.76   0.5000   0.5000
 5.76   0.9091   0.0909
 7.40   0.9977   0.0023
```

Near the pKa the populations change rapidly with pH. This is the buffering
region. Far from the pKa one form dominates.

A standard molecular-mechanics force field cannot vary protonation continuously
during a normal simulation. The topology contains a fixed number of atoms,
bonds, angles, and charges. Changing the protonation state therefore means
choosing a different chemical graph and rebuilding the force-field parameters.

## 5. Tautomers and proton rearrangement

Tautomers are constitutional isomers that differ in the position of a proton
and the arrangement of double bonds. Common examples are keto-enol,
lactam-lactim, and imine-enamine pairs.

Tautomerism is not resonance. Resonance structures have the same nuclear
positions and differ only in the distribution of electrons. Tautomers have
different connectivities and require the movement of at least one proton.

Tautomers of the same molecule have the same molecular formula and the same
number of protons. Their equilibrium is therefore a composition-conserving
isomerisation, and the Boltzmann expression below is appropriate. A protonation
equilibrium is different: it changes the number of protons and therefore the
composition of the species. Protonation populations need the pH-dependent
proton chemical potential introduced in section 6.1, not a bare tautomer
Boltzmann ratio.

For two tautomers A and B,

$$ K_T = \frac{[\mathrm{B}]}{[\mathrm{A}]} = e^{-\Delta G^\circ/(RT)}, $$

where $\Delta G^\circ = G_B^\circ - G_A^\circ$. The equilibrium populations are

$$ p_A = \frac{1}{1+e^{-\Delta G^\circ/(RT)}}, \qquad p_B = \frac{1}{1+e^{+\Delta G^\circ/(RT)}}. $$

At 298.15 K the example gives:

```text
delta G     p(A)     p(B)
    0.0   0.5000   0.5000
   -2.0   0.3086   0.6914
    2.0   0.6914   0.3086
    5.0   0.8826   0.1174
```

A negative $\Delta G^\circ$ favours B. A difference of only a few kJ/mol is
enough to change the populations substantially.

Tautomer populations are an equilibrium statement. Tautomerisation can be slow
or catalysed, so the equilibrium may not be reached on the timescale of an
experiment or a simulation. The calculation also assumes that the relevant
tautomers have been identified in the first place.

Tautomers matter because they can have different hydrogen-bond donors and
acceptors, different dipole moments, different ring electronics, and different
metal-binding geometries. A force-field topology built for one tautomer is not
automatically valid for another.

## 6. Molecules with multiple ionizable groups

A molecule can have several ionizable sites. If the sites are independent, the
probability of a microstate is the product of the site fractions:

$$ P(\mathbf{s}) = \prod_i p_i(s_i), \qquad s_i \in \{0,1\}. $$

Here $s_i=1$ means site $i$ is protonated and $s_i=0$ means it is deprotonated.
A molecule with $n$ sites has $2^n$ microstates. The mean net charge is

$$ \langle Q \rangle = \sum_{\mathbf{s}} P(\mathbf{s}) Q(\mathbf{s}). $$

The independent-site approximation ignores electrostatic interactions between
sites. It is a useful starting point, not a general truth.

For two sites with pKa values 4.0 and 9.0 at pH 7.0, using an acid site with
charges 0 and -1 and a base site with charges +1 and 0:

```text
  state  probability  charge
     DD     0.009891    -1.0
     HD     0.000010     0.0
     DH     0.989110     0.0
     HH     0.000989     1.0
  mean protonation: 0.0010, 0.9901
  mean net charge:  -0.0089
```

The labels list site 0 first and site 1 second; H means protonated and D means
deprotonated. The dominant state is DH: the acid is deprotonated and the base is
protonated, giving a near-neutral molecule.

### 6.1 Coupled sites

When sites interact, their pKa values are not independent. A nearby negative
charge stabilises a protonated base; a nearby positive charge makes protonation
less favourable. The correct bookkeeping uses microstate free energies:

$$ P_{\mathbf{s}} = \frac{e^{-G_{\mathbf{s}}/(RT)}}{\sum_{\mathbf{t}} e^{-G_{\mathbf{t}}/(RT)}}. $$

For a protonated site at a given pH, the pH-dependent contribution is

$$ G_i^{\mathrm{prot}} - G_i^{\mathrm{deprot}} = RT\ln 10\left(\mathrm{pH}-\mathrm{p}K_{a,i}\right). $$

This term is the bridge between the Henderson-Hasselbalch picture and the
Boltzmann picture. It includes the chemical potential of the proton at the
chosen pH. With zero coupling between sites, the microstate model reduces to
the product of independent Henderson-Hasselbalch fractions. The coupling terms
describe interactions between sites; the pH term describes the proton
chemical potential.

Adding a positive interaction energy when both sites are protonated penalises
the HH state. With a +2.0 kJ/mol HH coupling, the example gives:

```text
  state   energy  probability
     DD     0.00     0.009897
     HD    17.12     0.000010
     DH   -11.42     0.989652
     HH     7.71     0.000442
```

The independent model would put HH at 5.71 kJ/mol. The coupling raises it to
7.71 kJ/mol and reduces its population. In a protein, such couplings can shift
a residue's effective pKa by several units.

The number of microstates grows exponentially. For a molecule with ten
ionizable sites there are 1024 states, and for twenty sites there are more than
one million. Many states may be negligible, but enumerating them is still a
combinatorial problem.

## 7. Protein residues

Amino-acid side chains are weak acids and bases. The following values are
typical model-peptide values, not universal constants:

| Group | Typical pKa | Protonated charge | Deprotonated charge | Note |
| --- | --- | --- | --- | --- |
| Asp side chain | 3.9 | 0 | -1 | carboxylic acid |
| Glu side chain | 4.1 | 0 | -1 | carboxylic acid |
| His side chain | 6.0 | +1 | 0 | imidazole |
| Cys side chain | 8.3 | 0 | -1 | thiol |
| Tyr side chain | 10.1 | 0 | -1 | phenol |
| Lys side chain | 10.5 | +1 | 0 | amine |
| Arg side chain | 12.5 | +1 | 0 | guanidinium |
| N-terminus | 8.0 | +1 | 0 | amine |
| C-terminus | 3.1 | 0 | -1 | carboxylic acid |

At pH 7.4, Asp and Glu are mostly deprotonated, Lys and Arg are mostly
protonated, and His is close to its pKa. This is why histidine is often
described as a pH-sensitive switch in enzymes and why its protonation state can
change dramatically in a binding site.

The table is deliberately approximate. In a protein, the local dielectric
environment, nearby charges, hydrogen bonds, and desolvation can shift a pKa by
several units. A buried aspartate may remain protonated; a histidine in a
cationic pocket may have a much lower pKa than 6.0.

Protein net charge therefore depends on pH. The isoelectric point is the pH at
which the net charge is zero. Below it the protein tends to be positively
charged; above it, negatively charged. This is a qualitative statement; the
actual charge follows from all ionizable groups and their environment.

## 8. Consequences for formal charge and interactions

Protonation changes formal charge, but formal charge is not the whole story.
Formal charge is an integer bookkeeping quantity; partial charges are a
model-dependent continuous distribution used by the force field. Protonation
also changes:

- the number of hydrogen-bond donors and acceptors;
- the direction and strength of salt bridges;
- the electrostatic potential around the molecule;
- the desolvation penalty for burying a charge;
- the preferred tautomer and the available hydrogen-bond geometry;
- the atom types, bonded parameters, partial charges, and exclusions used by a
  force field.

A carboxylate is a hydrogen-bond acceptor and a salt-bridge partner. The neutral
carboxylic acid is a hydrogen-bond donor. A neutral amine can accept a hydrogen
bond; an ammonium group donates hydrogen bonds and carries a positive charge.
The difference is not a small perturbation.

For a force field, changing protonation means rebuilding the molecular graph and
its parameters. Simply scaling the partial charges of a fixed topology would not
represent a different chemical state.

In a binding site, the protein and ligand protonation states are coupled. A
charged ligand may prefer a complementary charged residue, but burying that
charge costs desolvation free energy. A protonation change can therefore alter
both the interaction energy and the conformational preference. Those effects
belong to the binding free energy, not to a single geometric contact.

## 9. Consequences for docking and molecular dynamics

Standard molecular dynamics uses a fixed topology. The protonation and tautomer
states must be chosen before the simulation starts. Constant-pH methods exist,
but they are specialised algorithms with their own approximations and are not
the default behaviour of a conventional force field.

Docking has a similar problem. A ligand with several ionizable groups can be
prepared in several states. The choice affects formal charge, hydrogen bonding,
and the electrostatic score. Many docking workflows enumerate a small set of
plausible states and dock each one. That is state enumeration, not pKa
prediction.

A docking score is also not a pH-dependent binding free energy. It may include
electrostatics and desolvation terms, but it usually does not solve the coupled
protonation equilibrium of the protein-ligand complex.

Practical consequences:

1. Record the intended pH and the protonation/tautomer state with the structure.
2. Do not assume that the neutral form is the relevant form at physiological pH.
3. Check whether the force field has parameters for the chosen state.
4. Treat a pKa shift in a binding site as a model result that needs validation.
5. Keep geometric interaction analysis separate from a claim about free-energy
   contribution.

## 10. RDKit representation and protonation-state enumeration

RDKit can represent charged and neutral forms and can enumerate tautomers. It
should not be treated as a universal pKa predictor.

### 10.1 Representing protonation states

```python
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

states = {
    "acetic acid": "CC(=O)O",
    "acetate": "CC(=O)[O-]",
    "methylamine": "CN",
    "methylammonium": "C[NH3+]",
    "imidazole": "c1c[nH]cn1",
    "imidazolium": "c1c[nH+]c[nH]1",
}

for name, smiles in states.items():
    molecule = Chem.MolFromSmiles(smiles)
    print(
        f"{name:14s} charge={Chem.GetFormalCharge(molecule):+d} "
        f"formula={rdMolDescriptors.CalcMolFormula(molecule)}"
    )
```

Observed output:

```text
acetic acid    charge=+0 formula=C2H4O2
acetate        charge=-1 formula=C2H3O2-
methylamine    charge=+0 formula=CH5N
methylammonium charge=+1 formula=CH6N+
imidazole      charge=+0 formula=C3H4N2
imidazolium    charge=+1 formula=C3H5N2+
```

The formula and formal charge change when the proton is added or removed. These
are different chemical graphs, not two drawings of the same graph.

### 10.2 Neutralisation is not a pKa model

RDKit's `Uncharger` can convert a charged representation to a neutral one where
a reasonable neutral form exists:

```python
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

uncharger = rdMolStandardize.Uncharger()
for smiles in ["CC(=O)[O-]", "C[NH3+]", "c1c[nH+]c[nH]1"]:
    molecule = Chem.MolFromSmiles(smiles)
    print(smiles, "->", Chem.MolToSmiles(uncharger.uncharge(molecule)))
```

Observed output:

```text
CC(=O)[O-] -> CC(=O)O
C[NH3+] -> CN
c1c[nH+]c[nH]1 -> c1c[nH]cn1
```

This is a standardisation operation. It does not tell you which form dominates
at a given pH.

### 10.3 Tautomer enumeration

```python
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

molecule = Chem.MolFromSmiles("Cc1c[nH]cn1")
enumerator = rdMolStandardize.TautomerEnumerator()
for tautomer in enumerator.Enumerate(molecule):
    print(Chem.MolToSmiles(tautomer))
```

Observed output:

```text
Cc1c[nH]cn1
Cc1cnc[nH]1
```

Enumeration gives possible structures. It does not give their energies,
populations, or interconversion rates. The toolkit is doing graph chemistry,
not thermodynamics.

### 10.4 Three different computational problems

| Problem | Question | Input | Output |
| --- | --- | --- | --- |
| enumeration | which protonation or tautomer graphs exist? | molecular graph, toolkit rules | a list of valid states |
| population | what fraction of each state is present? | pKa values or microstate energies, pH, temperature | state probabilities |
| prediction | what pKa does this group have? | structure and environment | an estimated pKa |

RDKit helps mainly with the first problem. The `acid_base` helpers in this
module solve the second once pKa values or microstate energies are supplied.
The third is a separate research problem: pKa prediction in solution is
difficult, and prediction inside a protein binding site is harder still because
the local environment is part of the answer.

Conflating these problems is a common source of overconfidence. A toolkit that
can draw a charged form has not thereby predicted which form dominates at pH 7.4.

## 11. Validation and acceptance criteria

The implementation in `src/moleculab/applications/acid_base.py` is checked
against closed-form values, not against itself:

```bash
uv run pytest tests/test_acid_base.py -v
uv run python examples/protonation.py
```

Key checks:

| Check | Expected result |
| --- | --- |
| acid at pH = pKa | deprotonated fraction 0.5 |
| one pH unit above pKa | deprotonated fraction 10/11 |
| one pH unit below pKa | deprotonated fraction 1/11 |
| acetic acid at pH 7.4 | deprotonated fraction 0.997714 |
| ammonium at pH 7.4 | protonated fraction 0.986071 |
| two-site microstates | probabilities sum to 1 |
| two-site mean charge | -0.008902 |
| Boltzmann ratio for $\Delta G=RT\ln 9$ | high-energy population 0.1 |
| tautomer with $\Delta G=0$ | equal populations |

Acceptance criteria for this module:

- the Henderson-Hasselbalch equation is derived, not merely stated;
- acid and base conventions are distinguished explicitly;
- pH is described as changing populations of discrete states;
- independent-site and coupled-site models are distinguished;
- tautomers are distinguished from resonance structures;
- same-composition tautomer weighting is distinguished from pH-dependent
  protonation weighting;
- enumeration, population estimation, and pKa prediction are distinguished;
- protein pKa values are marked as approximate and environment-dependent;
- the force-field consequences of changing protonation are stated;
- the example and tests reproduce the documented numbers.

## 12. Exercises

### Calculation

1. Derive the deprotonated fraction from the Henderson-Hasselbalch equation.
2. For a pKa of 4.76, compute the fraction deprotonated at pH 3.76, 4.76, and
   5.76.
3. At what pH is a carboxylic acid with pKa 4.76 90% deprotonated?
4. For two independent sites with pKa 4.0 and 9.0 at pH 7.0, compute all four
   microstate probabilities by hand.
5. What $\Delta G^\circ$ gives a 99:1 tautomer ratio at 298.15 K?
6. For a protein group with pKa 6.0, what fraction is protonated at pH 7.4?

### Reasoning

1. Why does pH change the populations of states rather than continuously
   changing the partial charges of one graph?
2. Why is a tautomer not a resonance structure?
3. Why can a protein pKa differ substantially from a model-peptide pKa?
4. Why is tautomer enumeration not a prediction of tautomer populations?
5. Why is a docking score not automatically a pH-dependent binding free energy?
6. Why does adding a proton require rebuilding a force-field topology?

### Implementation

1. Run `examples/protonation.py` and reproduce the tables.
2. Change the acetic-acid pKa and observe how the titration curve moves.
3. Add a third ionizable site to the microstate calculation and explain the
   growth from four to eight states.
4. Use RDKit to build the neutral and charged forms of glycine and compare their
   formal charges and formulas.
5. Enumerate the tautomers of 4-methylimidazole and explain why there are two.
6. Compute the mean net charge of a small peptide from a table of residue pKa
   values, stating clearly that the result uses an independent-site model.

## 13. Summary

- Acids and bases exchange protons; $K_a$ and pKa quantify the equilibrium.
- The Henderson-Hasselbalch equation relates pH, pKa, and the ratio of
  protonated to deprotonated forms.
- pH selects populations of discrete chemical states; it does not continuously
  charge one fixed graph.
- Multiple ionizable sites require microstate bookkeeping. Independent sites
  factorise; coupled sites require microstate free energies.
- Tautomers are constitutional isomers, not resonance structures. Their
  populations follow a Boltzmann distribution when they have the same
  composition, but the equilibrium may not be fast.
- Protonation changes the number of protons, so its populations require the
  pH-dependent proton chemical potential, not a bare tautomer ratio.
- Protein residue pKa values are approximate and environment-dependent.
- Protonation changes formal charge, hydrogen bonding, electrostatics,
  desolvation, and force-field parameters.
- Docking and MD require an explicit chemical-state choice; enumeration is not
  prediction.
- RDKit can represent and enumerate states, but it is not a universal pKa
  engine.

**Previous:** [03 — Molecular properties](03_molecular_properties.md).  
**Next:** [05 — Stereochemistry and molecular recognition](05_stereochemistry.md).
