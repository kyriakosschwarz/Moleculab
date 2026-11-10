# Track C, module 10 — From structural observations to thermodynamics

> **Where this fits.** Track C modules 01–09 built structural, interaction,
> trajectory, docking, and SAR analysis. This module slows down and asks what
> those observations can and cannot say about binding free energy.
>
> **What you should already know.** Track A free-energy and statistical
> mechanics material; Track C modules 03–04 and 08.
>
> **What you will be able to do afterwards.** Explain why occupancy and contact
> counts are not binding free energies, connect ensemble populations to
> thermodynamic quantities, and interpret a structural observation without
> overclaiming.

## 1. Learning objectives

**Core**

- Distinguish structural observations from thermodynamic quantities.
- Explain ensemble populations and the configurational partition function.
- Explain why the lowest-energy configuration need not dominate.
- Identify conformational-entropy and desolvation contributions.
- Relate binding free energy to the dissociation constant.
- Explain why contact occupancy is not binding free energy.

**Deeper insight**

- A persistent interaction can be offset by strain, desolvation, or entropy.
- Enthalpy–entropy compensation is common in binding.
- A trajectory describes sampled configurations, not automatically an
  equilibrium ensemble.
- A free-energy calculation requires a defined model and adequate sampling.
- Structural evidence and thermodynamic evidence complement each other but are
  not interchangeable.

**Looking ahead**

- The capstone integrates all three tracks.
- Track A's FEP, TI, and MBAR methods are the foundation for quantitative
  free-energy calculations.

## 2. What structural analysis measures

Track C has produced descriptive observations:

- binding-site atom and residue contacts;
- hydrogen bonds and other geometric interactions;
- interaction fingerprints;
- sampled-frame occupancy;
- contact lifetimes and transitions;
- docking poses and pose-recovery statistics.

These are statements about structures, poses, or sampled frames. They are not
automatically statements about thermodynamics.

A contact table can say:

> this residue is near this ligand atom in 80% of the analysed frames.

It cannot, by itself, say:

> this contact contributes −8 kJ/mol to binding.

The second statement needs a thermodynamic model and a free-energy calculation.

## 3. From configurations to ensembles

A molecular configuration is one set of coordinates. A trajectory samples many
configurations.

The configurational partition function is

$$ Z_{\mathrm{conf}}=\int e^{-\beta U(\mathbf R)}\,d\mathbf R, $$

where $\beta=1/(k_BT)$ and $U(\mathbf R)$ is the potential energy.

The Helmholtz free energy is

$$ A=-k_BT\ln Z_{\mathrm{conf}}. $$

For a discrete set of states with energies $E_i$ and degeneracies $g_i$,

$$ p_i=\frac{g_i e^{-\beta E_i}}{Z}, \qquad Z=\sum_i g_i e^{-\beta E_i}. $$

Two ideas follow:

1. The population of a state depends on its energy.
2. It also depends on how many configurations or states have that energy.

A single low-energy structure is not the whole ensemble.

## 4. Energy versus free energy

The lowest-energy configuration is not necessarily the most populated state.

A higher-energy state can dominate when it has many more accessible
configurations. The example demonstrates this with two states:

```text
state A: E = 0 kJ/mol, degeneracy 1
state B: E = 5 kJ/mol, degeneracy 100
population A = 0.070
population B = 0.930
```

This is not a chemical binding model. It is a demonstration that energy alone
does not determine population.

Entropy is the bridge between the number of accessible states and the free
energy.

## 5. Solvent and conformational entropy

Binding changes the accessible configurations of the whole system, not just the
ligand.

Important contributions include:

- **ligand conformational entropy**: a flexible ligand may lose torsional
  freedom when it binds;
- **protein conformational entropy**: side chains and loops may become more or
  less mobile;
- **desolvation**: charged and polar groups lose favourable interactions with
  water;
- **water reorganisation**: binding can displace, trap, or reorganise waters;
- **residual motion**: the bound complex still has vibrational, rotational, and
  conformational entropy.

A persistent hydrogen bond may be locally favourable but can be accompanied by:

- a large desolvation penalty;
- a loss of ligand conformational entropy;
- a strained ligand conformation;
- a protein conformational change;
- water networks that must be disrupted.

These terms can offset or reverse the apparent benefit of the contact.

## 6. Binding thermodynamics

For the binding equilibrium

$$ \mathrm{P}+\mathrm{L}\rightleftharpoons\mathrm{PL}, $$

the dissociation constant is

$$ K_d=\frac{[\mathrm{P}][\mathrm{L}]}{[\mathrm{PL}]}. $$

With a standard-state concentration $C^\circ=1\ \mathrm{M}$,

$$ \Delta G^\circ=RT\ln\left(\frac{K_d}{C^\circ}\right). $$

A smaller $K_d$ gives a more negative $\Delta G^\circ$.

The relation is

$$ \frac{K_d}{C^\circ}=e^{\Delta G^\circ/(RT)}. $$

This is a thermodynamic statement under a defined standard state. It is not
obtained by counting contacts.

$K_d$ is also not the same as IC50 or Ki. Converting between assay measurements
and thermodynamic constants requires assumptions about the assay and binding
mechanism.

## 7. Why occupancy is not binding free energy

A sampled-frame occupancy is an indicator average:

$$ O=\frac{1}{N_{\mathrm{valid}}}\sum_t \mathbf{1}_{\mathrm{contact}}(t). $$

It describes how often a geometric criterion is satisfied in the analysed
frames.

Binding free energy is a free-energy difference between macroscopic states. It
depends on:

- the interaction energy;
- the configurational entropy of ligand, protein, and complex;
- desolvation and solvent entropy;
- protonation and tautomer states;
- the standard state and temperature;
- the force field or energy model;
- the quality of phase-space sampling.

Occupancy and free energy are related through the ensemble, but one is not a
substitute for the other.

A high-occupancy contact can be favourable, neutral, or even part of a
strained, high-energy binding mode. The occupancy alone does not tell you which.

## 8. A two-ligand demonstration

The example constructs two hypothetical ligands:

| Ligand | Persistent contacts | Interaction | Desolvation | Ligand entropy | Protein entropy | Water | Total |
| --- | --- | --- | --- | --- | --- | --- | --- |
| A | 4 | −45 | +20 | +12 | +8 | +5 | 0 |
| B | 2 | −30 | +8 | +4 | +3 | +2 | −13 |

All values are kJ/mol and are purely illustrative. The important result is:

- A has more persistent contacts;
- B has the more favourable hypothetical total.

The numbers are not a force field or a free-energy calculation. They show the
bookkeeping that a contact count ignores.

The example also shows the degeneracy effect:

```text
population A = 0.070
population B = 0.930
```

A higher-energy but much more numerous set of states can dominate.

## 9. What MD can and cannot establish

A well-sampled MD trajectory can provide:

- structural fluctuations;
- contact occupancies and lifetimes;
- conformational populations;
- hydrogen-bond persistence;
- water residence behaviour;
- some entropy estimates through specialised methods.

It does not automatically provide:

- a rigorous binding free energy;
- an experimental affinity;
- a unique protonation or tautomer state;
- a complete picture of slow conformational changes;
- a guarantee that the force field is accurate.

A free-energy difference requires a defined alchemical or thermodynamic path
and adequate sampling. Track A's FEP, TI, and MBAR modules are the starting
point for that work.

## 10. What structural analysis is still good for

This module is not arguing that structural analysis is unimportant. It is
arguing that the claims must match the method.

Structural analysis is excellent for:

- generating mechanistic hypotheses;
- identifying plausible interactions;
- detecting clashes or missing interactions;
- comparing binding modes;
- checking whether a proposed mechanism is geometrically possible;
- deciding which interactions deserve further thermodynamic study.

Thermodynamic calculations are needed when the question is:

- how strongly does this ligand bind?
- what is the free-energy cost of this transformation?
- how much does this contact contribute under a defined model?
- which of two binding modes is more stable?

## 11. Connecting to Track A

Track A provides the foundations for quantitative free-energy work:

- FEP and alchemical transformations;
- thermodynamic integration;
- MBAR;
- soft-core potentials for endpoint stability;
- Ewald electrostatics and solvent models;
- thermostats, barostats, and periodic boundaries.

Track C provides the structural interpretation and the recognition hypotheses
that motivate those calculations.

The capstone connects the tracks:

| Track | Contribution |
| --- | --- |
| A — physics | energy models, MD, ensembles, free-energy estimators |
| B — chemistry | molecular identity, protonation, descriptors, transformations |
| C — drug discovery | structures, binding sites, interactions, docking, SAR |

## 12. Implementation and example

This module intentionally adds no large analysis API.

The illustrative example is:

```bash
uv run python examples/drug_discovery/structural_observations_thermodynamics.py
```

It demonstrates:

- contact-count ranking versus a hypothetical component-based ranking;
- the difference between interaction, desolvation, and entropy terms;
- the effect of degeneracy on state populations;
- the statement that occupancy is not a free energy.

The numbers are teaching values, not a binding-affinity predictor.

## 13. Exercises

### Reasoning

1. Why does a persistent contact not prove a favourable binding contribution?
2. Why can a higher-energy state be more populated than a lower-energy one?
3. What contributions are missing from a contact count?
4. Why is occupancy not an equilibrium probability automatically?
5. Why is a free-energy calculation more demanding than a contact analysis?
6. Why is IC50 not the same as Kd?
7. Why can a ligand with fewer contacts bind more strongly?

### Calculation

1. For $K_d=1\ \mathrm{nM}$ and $T=298.15$ K, compute $\Delta G^\circ$.
2. For two states with equal degeneracy and a 5 kJ/mol gap, compute their
   Boltzmann populations at 298.15 K.
3. Repeat with a degeneracy ratio of 100:1 in favour of the higher-energy state.
4. A contact occupancy is 0.8 over 100 valid frames. What does that number
   establish, and what does it not establish?
5. A ligand loses 12 kJ/mol of conformational entropy on binding but gains
   20 kJ/mol of interaction energy. What is the net contribution before other
   terms?

### Implementation

1. Run the example and reproduce the component table.
2. Change the desolvation penalty of ligand A until it becomes more favourable
   than B.
3. Change the degeneracy of the higher-energy state and observe the population.
4. Add a missing frame and confirm occupancy and lifetime use different
   denominators.
5. Explain why the example does not call any FEP or TI function.

## 14. Summary

- Structural observations describe poses, contacts, and sampled frames.
- Thermodynamic quantities describe the statistical distribution of states.
- The partition function and free energy connect configurations to populations.
- Energy alone does not determine population; degeneracy and entropy matter.
- Desolvation, ligand entropy, protein entropy, and water reorganisation can
  offset favourable contacts.
- Binding free energy is related to $K_d$ through a thermodynamic standard
  state.
- Occupancy and contact lifetime are not binding free energy.
- Structural analysis generates hypotheses; free-energy methods test them.
- The capstone integrates structure, chemistry, physics, and thermodynamics.

**Previous:** [09 — Lead optimisation and structure–activity relationships](09_lead_optimisation_sar.md).

**Next:** [Imatinib–BCR-ABL1 capstone](../capstone_imatinib.md).
