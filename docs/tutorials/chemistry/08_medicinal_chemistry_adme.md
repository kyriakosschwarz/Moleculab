# Track B, module 08 — Medicinal chemistry and ADME fundamentals

> **Where this fits.** Track B modules 04–07 established chemical state,
> stereochemistry, polarity, and hydrogen bonding. This module connects those
> properties to drug-like behaviour: absorption, distribution, metabolism,
> excretion, and toxicity.
>
> **What you should already know.** Protonation, pKa, polarity, hydrogen
> bonding, descriptors, and matched-pair SAR.
>
> **What you will be able to do afterwards.** Interpret common physicochemical
> trade-offs, distinguish LogP from LogD, explain the limits of drug-likeness
> rules, and separate measured ADME outcomes from calculated descriptors.

## 1. Learning objectives

**Core**

- Explain the four ADME processes.
- Relate solubility, permeability, ionisation, and lipophilicity to absorption.
- Distinguish LogP from LogD.
- Explain common metabolic and toxicity concerns.
- Interpret physicochemical trade-offs without reducing them to one score.

**Deeper insight**

- A calculated descriptor is not an experimental ADME measurement.
- Rule-of-Five compliance does not establish oral bioavailability.
- Structural alerts and PAINS flags are hypotheses, not toxicity predictions.
- Improving potency can worsen solubility, permeability, or clearance.
- Imatinib shows that high molecular weight and TPSA do not automatically make
  a compound non-drug-like.

**Looking ahead**

- The capstone must distinguish measured, calculated, predicted, and
  hypothesised properties.
- A production ADME model requires dedicated data, validation, and domain
  expertise.

## 2. What ADME means

| Process | Question |
| --- | --- |
| absorption | does the compound enter the bloodstream? |
| distribution | where does it go in the body? |
| metabolism | how is it chemically transformed? |
| excretion | how is it removed? |
| toxicity | does it cause harm or off-target effects? |

These are biological and pharmacokinetic questions. Molecular descriptors can
suggest hypotheses, but they do not measure these outcomes.

## 3. Absorption

Absorption depends on:

- solubility in the gut;
- permeability across membranes;
- ionisation at the relevant pH;
- lipophilicity;
- efflux transporters;
- formulation.

A useful qualitative balance is:

- very polar compounds may have low passive permeability;
- very lipophilic compounds may have low aqueous solubility;
- ionised compounds may partition poorly into membranes;
- intermediate lipophilicity often balances the two.

This is a trade-off, not a single threshold.

### 3.1 Solubility

Solubility depends on crystal packing, ionisation, particle form, and solvent.
A calculated LogP is not a solubility measurement.

### 3.2 Permeability

Passive permeability depends on membrane partitioning and molecular size and
polarity. TPSA is a useful descriptor, but it is not a permeability assay.

## 4. Distribution

After absorption, a compound can:

- bind plasma proteins;
- partition into tissues;
- enter cells;
- be transported actively;
- remain in the plasma.

Plasma protein binding reduces the free fraction available for target binding.
A high total plasma concentration does not necessarily mean a high free
concentration.

## 5. Metabolism

Metabolism usually makes a compound more polar and easier to excrete.

Common phase I transformations:

- oxidation;
- reduction;
- hydrolysis.

Common phase II transformations:

- glucuronidation;
- sulfation;
- acetylation;
- glutathione conjugation.

Metabolic soft spots are structural positions that are readily oxidised or
conjugated. They can limit half-life or produce reactive metabolites.

A metabolic soft spot is a structural hypothesis. It is not a measured
metabolic rate.

## 6. Excretion

Excretion removes the compound and its metabolites.

Routes include:

- renal clearance into urine;
- biliary clearance into faeces;
- other minor routes.

Renal clearance depends on glomerular filtration, secretion, and reabsorption.
These are physiological processes, not simple descriptor functions.

## 7. Toxicity and structural alerts

Toxicity can arise from:

- reactive metabolites;
- off-target binding;
- mitochondrial effects;
- immune reactions;
- accumulation;
- formulation effects.

Structural alerts are substructures associated with known liabilities. They
are hypotheses, not predictions of clinical toxicity.

PAINS are not toxicity filters. They flag assay interference or promiscuous
behaviour that can mislead screening.

A structural alert should trigger investigation, not automatic rejection.

## 8. Physicochemical trade-offs

Common descriptors include:

| Descriptor | What it suggests | What it does not measure |
| --- | --- | --- |
| molecular weight | size | permeability or solubility |
| cLogP | neutral-form lipophilicity | LogD at a given pH |
| TPSA | polar surface area | membrane permeability |
| HBD/HBA | hydrogen-bonding capacity | actual interaction strength |
| rotatable bonds | flexibility | conformational entropy |
| formal charge | net charge at the chosen state | distribution in the body |

A change that improves one property can worsen another.

Examples:

- adding lipophilic groups can increase potency but reduce solubility;
- adding polar groups can improve solubility but reduce permeability;
- removing a metabolic soft spot can improve stability but change target
  binding;
- increasing molecular weight can improve interactions but worsen oral
  absorption.

## 9. LogP versus LogD

### 9.1 LogP

For a neutral compound,

$$ \log P = \log_{10}\left(\frac{[\mathrm{solute}]_{\mathrm{octanol}}}{[\mathrm{solute}]_{\mathrm{water}}}\right). $$

LogP describes the neutral form. It does not include ionisation.

### 9.2 LogD

LogD is the distribution coefficient at a specified pH:

$$ \log D = \log_{10}\left(\frac{[\mathrm{solute}]_{\mathrm{octanol,total}}}{[\mathrm{solute}]_{\mathrm{water,total}}}\right). $$

It includes all relevant ionisation states.

For a monoprotic acid,

$$ \log D = \log P - \log_{10}\left(1 + 10^{\mathrm{pH}-\mathrm{p}K_a}\right). $$

For a monoprotic base,

$$ \log D = \log P - \log_{10}\left(1 + 10^{\mathrm{p}K_a-\mathrm{pH}}\right). $$

At pH = pKa, LogD is lower than LogP by $\log_{10}2 \approx 0.301$.

These equations assume that the ionised species has negligible octanol
partitioning, that the neutral species has the supplied LogP, and that ion
pairing, aggregation, and other non-ideal effects are absent. Real systems can
deviate from those assumptions.

The pKa values are inputs. The implementation does not predict pKa.

## 10. Rule of Five and drug-likeness rules

Lipinski's Rule of Five is a set of approximate oral-drug guidelines:

| Property | Common threshold |
| --- | --- |
| molecular weight | ≤ 500 |
| cLogP | ≤ 5 |
| hydrogen-bond donors | ≤ 5 |
| hydrogen-bond acceptors | ≤ 10 |

The rules are heuristics. They:

- were derived from a particular compound set;
- do not prove oral bioavailability;
- do not cover transporters, metabolism, or safety;
- are not a complete drug-design objective;
- can be intentionally violated by successful drugs.

Imatinib is a useful reminder: its molecular weight is about 494 and its TPSA
is about 86, yet it is an oral drug. Descriptor thresholds are not destiny.

## 11. Measured, calculated, predicted, and hypothesised

Keep the categories distinct:

| Category | Example | Status |
| --- | --- | --- |
| measured | solubility, permeability, clearance, IC50 | experimental |
| calculated | MW, cLogP, TPSA, HBD/HBA | descriptor |
| predicted | ADME model output | model estimate |
| hypothesised | metabolic soft spot, structural alert | mechanistic hypothesis |

A calculated TPSA is not a permeability measurement. A predicted clearance is
not a measured clearance. A structural alert is not a toxicity result.

## 12. RDKit example and validation

The example computes descriptors for several compounds and compares LogP with
LogD at pH 7.4:

```bash
uv run --with rdkit python examples/adme_tradeoffs.py
```

Observed output:

```text
  molecule       MW   cLogP    TPSA  HBD  HBA  rot charge
   aspirin   180.16    1.31    63.6    1    3    2      0
 ibuprofen   206.28    3.07    37.3    1    1    4      0
  atenolol   266.34    0.45    84.6    3    4    8      0
  caffeine   194.19   -1.03    61.8    0    3    0      0
  imatinib   493.62    4.59    86.3    2    7    7      0

LogP versus LogD at pH 7.4
  ibuprofen (acid, pKa 4.9): LogP=3.073, LogD=0.572
  atenolol (base, pKa 9.6): LogP=0.452, LogD=-1.751
```

The example interprets:

- ibuprofen is more lipophilic and less polar than atenolol;
- ibuprofen's acid is ionised at pH 7.4, lowering LogD;
- atenolol's base is protonated at pH 7.4, making LogD negative;
- imatinib has high molecular weight and TPSA but is still an oral drug.

The values are calculated descriptors and illustrative pKa inputs. They are
not measured ADME outcomes.

### 12.1 Acceptance criteria

An ADME fundamentals module is correct when:

- ADME outcomes are distinguished from molecular descriptors;
- LogP and LogD are clearly different;
- LogD includes pH and ionisation;
- pKa is treated as an input, not predicted;
- solubility, permeability, metabolism, and toxicity are not reduced to one
  score;
- Rule-of-Five limitations are explicit;
- structural alerts and PAINS are not presented as toxicity predictions;
- measured, calculated, predicted, and hypothesised categories remain
  separate.

## 13. Common misconceptions

- **LogP is LogD.** LogP is the neutral-form partition coefficient; LogD is
  pH-dependent.
- **A high LogP always means good absorption.** It can reduce solubility.
- **High TPSA always means poor permeability.** It is a descriptor, not an
  assay.
- **Rule-of-Five compliance means oral bioavailability.** It is a guideline.
- **A structural alert means toxicity.** It is a hypothesis.
- **PAINS are toxicity filters.** They flag assay interference.
- **A calculated ADME value is a measurement.**
- **Improving potency improves the candidate.** It can worsen other properties.

## 14. Exercises

### Reasoning

1. Why is LogD more relevant than LogP at physiological pH?
2. Why can a high-LogP compound have poor oral absorption?
3. Why is TPSA not a permeability measurement?
4. Why is Rule-of-Five compliance not proof of bioavailability?
5. Why are structural alerts not toxicity predictions?
6. Why can a metabolic soft spot be a hypothesis rather than a measurement?
7. Why is imatinib a useful counterexample to simple drug-likeness rules?

### Calculation

1. For LogP = 3.0, pKa = 4.9, and pH = 7.4, compute LogD for an acid.
2. For LogP = 1.0, pKa = 9.6, and pH = 7.4, compute LogD for a base.
3. At pH = pKa, by how much does LogD differ from LogP?
4. A compound has MW 520, LogP 5.4, HBD 4, HBA 8. How many Rule-of-Five flags?
5. Explain why a negative LogD does not prove a compound is inactive.

### Implementation

1. Run `examples/adme_tradeoffs.py` and reproduce the output.
2. Add another compound and interpret its physicochemical trade-offs.
3. Change the pH and observe how LogD changes for an acid and a base.
4. Compare a compound with high TPSA and one with low TPSA.
5. Explain why the example does not predict oral bioavailability.

## 15. Summary

- ADME describes absorption, distribution, metabolism, excretion, and toxicity.
- Physicochemical descriptors suggest hypotheses but are not measurements.
- LogP describes the neutral form; LogD includes ionisation at a specified pH.
- pKa is an input to LogD, not a prediction.
- Rule-of-Five is a guideline, not a bioavailability test.
- Structural alerts and PAINS are not toxicity predictions.
- Potency improvements can worsen solubility, permeability, or clearance.
- Measured, calculated, predicted, and hypothesised data must remain separate.
- A successful drug can violate simple descriptor rules.

**Previous:** [07 — Hydrogen bonding and molecular recognition](07_hydrogen_bonding.md).

**Next:** [Learning tracks](../tracks.md) and the
[imatinib capstone](../capstone_imatinib.md).
