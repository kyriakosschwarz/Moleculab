# Track C, module 09 — Lead optimisation and structure–activity relationships

> **Where this fits.** Track C modules 01–08 built structural and trajectory
> analysis tools. This module turns to medicinal-chemistry decisions: what
> molecular change should be made, and what evidence supports it?
>
> **What you should already know.** Track B modules 02–04 and Track C modules
> 03–08.
>
> **What you will be able to do afterwards.** Compare a matched molecular pair,
> interpret pIC50 changes, calculate a ligand-efficiency proxy, and separate
> measured activity, calculated descriptors, docking scores, and mechanistic
> hypotheses.

## 1. Learning objectives

**Core**

- Explain what a structure–activity relationship is.
- Compare a matched molecular pair.
- Convert IC50 to pIC50 and interpret fold changes.
- Interpret common physicochemical descriptors.
- Calculate a ligand-efficiency proxy.
- Distinguish measured, calculated, computational, and hypothetical evidence.

**Deeper insight**

- IC50, Ki, and Kd are not interchangeable without assumptions.
- Ligand efficiency is a size-normalised potency proxy, not a free energy.
- Docking and interaction fingerprints generate mechanistic hypotheses, not
  proof of affinity.
- Potency, selectivity, solubility, permeability, and metabolic stability are
  different objectives.
- A matched pair is evidence about a local transformation, not a universal SAR
  rule.

**Looking ahead**

- The capstone will connect structural analysis to real medicinal-chemistry
  reasoning.
- Future modules could add explicit selectivity and ADME datasets, but not a
  black-box affinity predictor.

## 2. From structural observation to decision

Structural analysis answers questions such as:

- where does the ligand bind?
- which residues contact it?
- which interactions are present?
- how persistent are the contacts?

Lead optimisation asks a different question:

> which chemical change is likely to improve the candidate?

That requires evidence from several categories:

| Category | Example |
| --- | --- |
| measured activity | IC50, Ki, or Kd from an assay |
| calculated descriptors | molecular weight, cLogP, TPSA, HBD/HBA, rotatable bonds |
| computational score | docking score or interaction fingerprint |
| mechanistic hypothesis | “the fluorine may improve shape complementarity” |

These categories must remain separate. A docking score is not a measured
potency, and a descriptor is not evidence of binding.

## 3. Structure–activity relationships and matched pairs

A structure–activity relationship (SAR) relates a chemical change to a change
in measured activity.

A **matched molecular pair** is two compounds that differ by one well-defined
transformation, such as:

- hydrogen to fluorine;
- methyl to ethyl;
- chlorine to methyl;
- one ring atom substitution.

For compounds A and B with the same activity measurement type,

$$ \Delta \mathrm{pActivity} = \mathrm{pActivity}_B - \mathrm{pActivity}_A. $$

A positive value means B is more potent on the measured scale.

The fold change is

$$ \text{fold change} = 10^{\Delta \mathrm{pActivity}}. $$

Equivalently, for IC50,

$$ \Delta \mathrm{p}IC_{50}=\mathrm{p}IC_{50,B}-\mathrm{p}IC_{50,A}, \qquad \frac{IC_{50,A}}{IC_{50,B}}=10^{\Delta \mathrm{p}IC_{50}}. $$

A positive delta means compound B has the lower IC50. A one-unit increase in
pIC50 corresponds to a ten-fold lower IC50.

Matched pairs are useful because they reduce confounding differences between
compounds. They do not remove them: a small transformation can change
conformation, solubility, protonation, or off-target binding.

## 4. Activity scales

For an IC50 in molar units,

$$ \mathrm{p}IC_{50} = -\log_{10}\left(\frac{IC_{50}}{1\ \mathrm{M}}\right). $$

Examples:

| IC50 | pIC50 |
| --- | --- |
| 10 µM | 5 |
| 1 µM | 6 |
| 100 nM | 7 |
| 10 nM | 8 |

A change from 1 µM to 100 nM is a one-unit increase in pIC50 and a ten-fold
increase in potency.

### 4.1 IC50, Ki, and Kd

These quantities are related but not interchangeable:

- **IC50** is the concentration producing half-maximal inhibition in a
  particular assay.
- **Ki** is an inhibition constant derived under a binding model.
- **Kd** is an equilibrium dissociation constant.

Converting between them requires assumptions about the assay, substrate
concentration, and binding mechanism. Do not compare pIC50, pKi, and pKd as if
they were the same measurement.

### 4.2 Same measurement type is not the same as comparable measurements

Two compounds can both report an IC50 and still have measurements that are not
directly comparable.

Possible differences include:

- assay format;
- substrate concentration;
- target construct;
- temperature and buffer;
- laboratory or protocol;
- data-processing method.

Same measurement type is a software validation condition. Comparable
measurements are an experimental interpretation condition.

The implementation stores `source` and `protocol_id` on `AssayMeasurement`.
A matched-pair result reports both values and sets `assay_context_matches` when
the sources and protocol identifiers are present and identical.

This flag is conservative. Matching metadata does not prove biological
comparability, and different metadata does not automatically make a comparison
invalid. It does make the limitation visible.

## 5. Physicochemical properties

Common calculated descriptors include:

| Descriptor | Interpretation | Caution |
| --- | --- | --- |
| molecular weight | size | salts and counterions matter |
| cLogP | predicted lipophilicity of the neutral form | not a measured LogP |
| TPSA | approximate polar surface area | graph-based, not a solvent model |
| H-bond donors | donor count | environment matters |
| H-bond acceptors | acceptor count | environment matters |
| rotatable bonds | flexibility indicator | not a direct entropy calculation |
| formal charge | net charge at the chosen state | protonation-dependent |

Descriptors are calculated from a molecular representation. They are useful for
filtering and hypothesis generation, not for proving biological activity.

## 6. Ligand efficiency

A size-normalised potency proxy is

$$ \mathrm{LE} = \frac{RT\ln 10\ \mathrm{p}IC_{50}}{N_{\mathrm{heavy}}}, $$

where $N_{\mathrm{heavy}}$ is the number of heavy atoms. The result is in
kJ/mol per heavy atom.

At 298.15 K, $RT\ln 10 \approx 5.708$ kJ/mol. If pIC50 increases by one unit,
LE increases by about $5.708/N_{\mathrm{heavy}}$ kJ/mol per heavy atom.

### 6.1 Limitations

Ligand efficiency:

- is based on an IC50 proxy, not a measured binding free energy;
- depends on the assay and the choice of activity measurement;
- does not account for selectivity, solubility, permeability, or metabolism;
- can favour small, weakly potent fragments;
- should not be used as a universal ranking function.

It is a useful conversation tool, not a decision rule.

## 7. Structural interpretation

Docking poses, interaction fingerprints, and trajectory contacts can suggest why
an analogue behaves differently.

Use them to generate hypotheses such as:

- the new substituent fills a hydrophobic subpocket;
- the transformation removes a hydrogen bond;
- a clash forces a different ligand conformation;
- a contact is present in most frames but not in others.

Do not claim:

- that a docking score explains the measured potency;
- that one contact causes the affinity difference;
- that a shared interaction pattern proves a shared binding mode;
- that a persistent contact is a favourable free-energy contribution.

The structural evidence and the assay evidence are independent observations
that may or may not agree.

## 8. Selectivity and liabilities

Potency is only one objective.

A useful analogue may improve:

- selectivity against an anti-target;
- solubility;
- permeability;
- metabolic stability;
- safety or toxicity profile;
- synthetic accessibility.

A change that improves potency can worsen another property. Report the
trade-off rather than a single “better” label.

## 9. Data categories and uncertainty

Keep at least four categories separate:

| Category | Example | Uncertainty |
| --- | --- | --- |
| measured | IC50, Ki, Kd | assay variability, model assumptions |
| calculated | molecular weight, cLogP, TPSA | descriptor definition and version |
| computational | docking score, interaction fingerprint | scoring model and sampling |
| hypothesis | “fluorine improves packing” | mechanistic conjecture |

A matched-pair comparison should report:

- the transformation;
- the measured activity change and measurement type;
- descriptor changes;
- docking-score change and provenance;
- the hypothesis, clearly labelled as a hypothesis.

Do not combine them into one “compound quality” score.

## 10. Implementation and validation

The implementation is in `src/moleculab/analysis/sar.py`. It is
parser-independent: matched-pair detection, descriptor calculation, and
structure preparation belong to RDKit or another chemistry toolkit.

### 10.1 API

| Object or function | Purpose |
| --- | --- |
| `AssayMeasurement` | measured value, measurement type, source, protocol, notes |
| `CompoundRecord` | compound identity, assay, descriptors, docking provenance |
| `MatchedMolecularPair` | two compounds plus a transformation and hypothesis |
| `MatchedPairComparison` | measured, descriptor, docking, assay-context, and hypothesis fields |
| `pIC50(...)` | convert molar IC50 to pIC50 |
| `ligand_efficiency(...)` | size-normalised potency proxy |
| `compare_matched_pair(...)` | compare a prepared pair without a composite score |

### 10.2 Chemistry boundary

RDKit can provide descriptors and matched-pair detection. The core comparison
layer consumes the prepared records:

```text
SMILES and assay data
    → RDKit or another toolkit
    → compound record: descriptors, transformation, activity
    → Moleculab matched-pair comparison
    → measured / calculated / computational / hypothesis fields
```

Do not build a custom matched-pair or ADME engine inside the comparison layer.

### 10.3 Example

```bash
uv run python examples/matched_molecular_pairs.py
```

Observed output:

```text
Measured activity
  IC50 A: 1000.0 nM
  IC50 B: 100.0 nM
  pIC50 A: 6.000
  pIC50 B: 7.000
  delta pIC50: +1.000
  fold change in IC50: 10.000x
  assay context matches: True (assay-1/protocol-1 vs assay-1/protocol-1)

Calculated descriptors
  delta MW: +20.00
  delta cLogP: +0.50

Computational docking
  delta score: -0.600 (arbitrary units)

Ligand efficiency (IC50 proxy)
  A: 1.712 kJ/mol per heavy atom
  B: 1.816 kJ/mol per heavy atom
```

The example keeps measured activity, calculated descriptors, docking scores,
and the mechanistic hypothesis separate.

### 10.4 Tests

```bash
uv run pytest tests/test_sar.py -v
```

The tests cover:

- pIC50 conversion;
- assay measurement types and validation;
- ligand-efficiency calculation and validation;
- compound score provenance;
- matched-pair identity and transformation validation;
- activity deltas and fold changes;
- missing activity and missing docking scores;
- measurement-type mismatch;
- assay source and protocol preservation;
- assay-context flagging;
- descriptor-set mismatch;
- docking unit and source mismatch.

### 10.5 Acceptance criteria

A matched-pair comparison module is correct when:

- pIC50 uses molar concentration;
- IC50, Ki, and Kd remain distinct;
- same measurement type is not treated as proof of assay comparability;
- assay source and protocol metadata are preserved and flagged;
- fold change is consistent with the pActivity difference;
- descriptor sets are compared explicitly;
- docking scores retain unit and source;
- missing values produce missing deltas, not zero;
- measured, calculated, computational, and hypothesis fields remain separate;
- no single composite score is returned;
- ligand efficiency is labelled as a proxy with stated limitations.

## 11. Exercises

### Calculation

1. Convert 10 µM, 1 µM, and 100 nM to pIC50.
2. A compound changes from 2 µM to 50 nM. What is the pIC50 change and fold
   change?
3. Compute the ligand efficiency for pIC50 7 and 30 heavy atoms.
4. A transformation adds 20 Da and increases pIC50 by 0.5. Report the descriptor
   and activity changes separately.
5. A docking score changes by −0.6 arbitrary units. Explain why this is not a
   measured potency change.

### Reasoning

1. Why are IC50, Ki, and Kd not interchangeable?
2. Why is a matched pair useful but not conclusive?
3. Why is ligand efficiency not a free energy?
4. Why can a potency improvement worsen a development candidate?
5. Why is a docking score not an assay result?
6. Why should a mechanistic hypothesis be labelled separately?
7. Why is a composite “compound quality” score risky?

### Implementation

1. Run `examples/matched_molecular_pairs.py` and reproduce the output.
2. Add a second matched pair with a different transformation.
3. Add a missing assay value and confirm the activity delta is missing.
4. Add a descriptor and propagate its delta.
5. Change the docking-score source and confirm the comparison is rejected.

## 12. Summary

- SAR relates a chemical change to a change in measured activity.
- Matched molecular pairs reduce confounding differences but do not eliminate
  them.
- pIC50 is a logarithmic potency scale; a one-unit change is a ten-fold change.
- IC50, Ki, and Kd are different measurements.
- Descriptors are calculated properties, not assays.
- Ligand efficiency is a size-normalised potency proxy with clear limitations.
- Docking and interaction analysis generate mechanistic hypotheses.
- Potency, selectivity, solubility, permeability, and metabolism are separate
  objectives.
- Measured, calculated, computational, and hypothetical evidence must remain
  distinguishable.
- A single composite score hides the evidence the analysis is meant to teach.

**Previous:** [08 — Trajectory contact analysis](08_trajectory_contact_analysis.md).

**Next:** [10 — From structural observations to thermodynamics](10_structural_observations_to_thermodynamics.md).
