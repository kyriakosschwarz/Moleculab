# Units and conventions in structure-based drug discovery

> **Use this as a reference.** Structural biology, cheminformatics, simulation,
> and medicinal chemistry use different conventional units. This page collects
> the ones used in Track C and shows how they relate to Moleculab's internal
> units.

## 1. Length

| Quantity | Common unit | Moleculab unit | Conversion |
| --- | --- | --- | --- |
| atomic coordinates | Angstrom | nm | $1\ \mathrm{Angstrom}=0.1\ \mathrm{nm}$ |
| bond length | Angstrom | nm | same |
| structure resolution | Angstrom | not used internally | resolution is a quality measure, not a coordinate |
| coarse-grained bead size | nm | nm | — |

Angstrom is written as Å in prose. One Angstrom is $10^{-10}$ m and one
nanometre is $10^{-9}$ m, so one Angstrom is exactly 0.1 nm.

Typical scale references:

- a C-C single bond is about 1.5 Å = 0.15 nm;
- a hydrogen bond donor-acceptor distance is often 2.5–3.5 Å = 0.25–0.35 nm;
- a protein domain is often several nm across.

Convert once, at the structure-import boundary:

```python
position_nm = position_angstrom * 0.1
```

Do not mix Angstrom and nanometre coordinates in the same calculation.

## 2. Angles

Code uses radians; human-facing reports often use degrees.

$$ 180^\circ = \pi\ \mathrm{rad}. $$

A hydrogen-bond angle is commonly reported as the D-H...A angle in degrees.
Keep the convention explicit: the angle at the hydrogen is not the same as the
angle at the donor or acceptor.

## 3. Energy

Moleculab uses kJ/mol internally. Medicinal-chemistry tables often use
kcal/mol.

$$ 1\ \mathrm{kcal/mol}=4.184\ \mathrm{kJ/mol}. $$

Useful reference values at 298.15 K:

| Quantity | kJ/mol | kcal/mol |
| --- | --- | --- |
| $RT$ | 2.479 | 0.592 |
| $RT\ln 10$ | 5.708 | 1.364 |

A ten-fold change in an equilibrium constant corresponds to about
$5.7$ kJ/mol at room temperature.

## 4. Temperature

Moleculab uses kelvin.

$$ T(\mathrm{K}) = T(^\circ\mathrm{C}) + 273.15. $$

Common values:

- 25 °C = 298.15 K;
- 37 °C = 310.15 K.

## 5. Mass

Moleculab uses daltons (Da).

$$ 1\ \mathrm{kDa}=1000\ \mathrm{Da}. $$

Average molar mass and monoisotopic mass are different quantities; see
[Track B module 02](../tutorials/chemistry/02_molecular_representation.md).

## 6. Concentration and affinity

Common concentration units:

| Unit | In molar |
| --- | --- |
| M | $1$ |
| mM | $10^{-3}$ |
| µM | $10^{-6}$ |
| nM | $10^{-9}$ |
| pM | $10^{-12}$ |

Affinity values such as $K_d$, $K_i$, and $IC_{50}$ are concentrations, not
free energies. They are also not interchangeable: $K_d$ describes an
equilibrium dissociation constant, while $IC_{50}$ depends on the assay
conditions.

A common summary is

$$ \mathrm{p}IC_{50} = -\log_{10}\left(\frac{IC_{50}}{1\ \mathrm{M}}\right). $$

A ten-fold increase in potency raises pIC50 by one unit.

For a dissociation constant, a useful approximate relation is

$$ \Delta G^\circ = RT\ln\left(\frac{K_d}{1\ \mathrm{M}}\right). $$

The sign convention matters: a tighter binder has a more negative
$\Delta G^\circ$. The equation assumes a standard-state concentration of 1 M
and does not by itself include protein or ligand conformational entropy.

## 7. pH and pKa

pH and pKa are dimensionless logarithmic quantities.

$$ \mathrm{pH} = -\log_{10} a_{\mathrm{H^+}} \approx -\log_{10}[\mathrm{H^+}]. $$

Physiological pH is about 7.4; it is not the only relevant condition. See
[Track B module 04](../tutorials/chemistry/04_protonation_and_tautomerism.md).

## 8. Pressure

Moleculab uses kJ/(mol nm³) internally and provides `pressure_bar()` for
conversion.

$$ 1\ \mathrm{bar}=10^5\ \mathrm{Pa}, \qquad 1\ \mathrm{atm}=1.01325\ \mathrm{bar}. $$

## 9. Time

Molecular simulation commonly uses picoseconds and nanoseconds.

$$ 1\ \mathrm{ns}=1000\ \mathrm{ps}, \qquad 1\ \mathrm{\mu s}=1000\ \mathrm{ns}. $$

## 10. Volume and surface area

Structural biology often reports volumes and surfaces in Å³ and Å².

$$ 1\ \mathrm{nm}^3 = 1000\ \mathrm{Angstrom}^3, \qquad 1\ \mathrm{nm}^2 = 100\ \mathrm{Angstrom}^2. $$

Solvent-accessible surface area (SASA) is often reported in Å². A pocket volume
may be reported in Å³. Do not compare these directly with Moleculab lengths
without converting.

## 11. Charge, occupancy, and B-factors

- Formal charge is in units of the elementary charge $e$.
- Occupancy is a fraction between 0 and 1.
- B-factors are commonly reported in Å² and describe displacement or
  uncertainty. They are not simulation temperatures.
- Partial charges are force-field parameters, not formal charges.

## 12. Quick conversion examples

```python
position_nm = position_angstrom * 0.1
distance_nm = distance_angstrom * 0.1
energy_kj = energy_kcal * 4.184
temperature_k = temperature_celsius + 273.15
concentration_m = concentration_nm * 1e-9
```

These conversions belong at the boundaries of the analysis. Keep the internal
calculation in one consistent unit system.

**Related:** [Track C index](../tutorials/drug_discovery/README.md), [protein structure](../tutorials/drug_discovery/01_protein_structure.md),
[binding-site geometry](../tutorials/drug_discovery/02_binding_sites.md), and
[interaction geometry](../tutorials/drug_discovery/03_protein_ligand_interactions.md).
