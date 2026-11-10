# Answer keys and expected magnitudes

> **Use this after trying the exercises.** These are answer sketches and
> order-of-magnitude checks, not a substitute for running the code. Stochastic
> exercises depend on seeds, sample counts, and system size.

## 1. Magnitude benchmarks

These are reference scales, not universal constants.

| Quantity | Reasonable scale | Notes |
| --- | --- | --- |
| thermal energy $RT$ at 300 K | 2.5 kJ/mol | exact for the stated temperature |
| typical hydrogen bond | 5–40 kJ/mol | strongly environment-dependent |
| Lennard-Jones well | 0.1–2 kJ/mol | depends on atom pair |
| covalent bond | 300–500 kJ/mol | depends on bond |
| bare $\pm e$ Coulomb at 0.3 nm | about 463 kJ/mol | unscreened vacuum value |
| O-H bond length | about 0.1 nm | water model dependent |
| H-O-H angle | about 104.5° | water model dependent |
| liquid water density | about 1 g/cm³ | about 33.4 molecules/nm³ |
| water self-diffusion at 300 K | about $2\times10^{-9}$ to $3\times10^{-9}$ m²/s | TIP3P/SPC/E order of magnitude |
| water MSD slope $6D$ | about $0.012$ to $0.016$ nm²/ps | from the diffusion coefficient |
| stable protein backbone RMSD | about 0.1–0.3 nm | relative to starting structure |
| stable ligand RMSD in a pocket | about 0.1–0.2 nm | depends on pose and fit |
| redocking success threshold | about 0.2 nm | symmetry and atom mapping matter |
| Ewald relative energy convergence | about $10^{-4}$ to $10^{-5}$ | depends on cutoff and system |

The diffusion conversion is

$$ 1\ \mathrm{nm^2/ps}=10^{-6}\ \mathrm{m^2/s}. $$

So $2.3\times10^{-9}$ m²/s is $0.0023$ nm²/ps.

## 2. Main tutorial exercise sketches

### 1. Harmonic stiffness and period

$$ T=2\pi\sqrt{\frac{m}{k}}. $$

Doubling $k$ divides the period by $\sqrt{2}$. Doubling $m$ multiplies it by
$\sqrt{2}$.

### 2. Lennard-Jones force sign change

The force is zero at

$$ r_{\min}=2^{1/6}\sigma. $$

Inside $r_{\min}$ the force is repulsive; outside it is attractive until the
tail becomes negligible.

### 3. Second torsion term

The energy landscape becomes a sum of periodic terms. It can gain more minima,
higher barriers, or a shifted preferred dihedral, depending on periodicity and
phase.

### 4. Minimising a stretched bond and distorted angle

The energy decreases monotonically for accepted steps. The bond length moves
toward $r_0$, the angle toward $\theta_0$, and the final maximum force falls
below the tolerance.

### 5. Thermostat at several targets

The mean temperature after warm-up should approach each target. The spread and
relaxation time depend on the coupling time and the degrees of freedom.

### 6. Minimum image

For a 1 nm box, a raw separation of 0.9 nm has a minimum-image separation of
0.1 nm. The interaction energy and force must use the minimum image.

### 7. Temperature, RDF, and MSD

Higher temperature generally broadens the RDF and increases the MSD slope. The
first RDF peak stays near the characteristic nearest-neighbour distance.

### 8. Solvent density and RDF

The first RDF peak should remain near the model's `sigma` scale; its height and
the long-range approach to 1 depend on density and temperature.

### 9. Larger FEP stiffness change

A larger perturbation increases the free-energy difference and usually the
variance. More samples or intermediate states are needed to keep the estimate
stable.

### 10. Water exclusions

For each rigid water molecule, the O-H bonds exclude the two 1-2 pairs and the
H-O-H angle excludes the 1-3 H-H pair. The exact count is three exclusions per
molecule.

### 11. Docking samples and radius

More samples can improve the best score, but the best score may not be at the
receptor centre. A docking score is an interaction energy, not a binding free
energy.

### 12. Ewald reciprocal cutoff

Increasing `reciprocal_cutoff` improves accuracy until the energy changes by
less than the chosen tolerance. Decreasing it loses accuracy and may produce
visible energy shifts.

### 13. Barostat target pressure

A higher target pressure compresses the box; a lower target pressure expands
it. The response magnitude depends on compressibility and coupling time.

### 14. Rigid versus flexible water

Rigid water keeps the O-H and H-H distances fixed. Flexible water allows bond
and angle fluctuations, so the O-H distribution has a finite width.

### 15. TI windows

More lambda windows usually reduce the integration error. For the harmonic
reference, compare with

$$ \Delta G=\frac32RT\ln\frac{k_1}{k_0}. $$

### 16. MBAR samples per state

Fewer samples per state increase the endpoint free-energy error. MBAR can still
combine states, but poor overlap cannot be fixed by reweighting alone.

### 17. Ewald cutoff in water

A larger Ewald cutoff generally improves accuracy and costs more time. Compare
energy and runtime, not just one of them.

### 18. FEP, TI, and MBAR

The three estimators should agree within their statistical errors. Their
variances and convergence behaviour differ, especially with few windows or poor
overlap.

### 19. QM grid spacing

For a second-order finite-difference Hamiltonian, halving the grid spacing
reduces the energy error by roughly four, until other errors dominate.

### 20. Monte Carlo step size

Too large a step gives low acceptance; too small a step explores slowly. The
best step size balances acceptance and decorrelation.

### 21. Solvated ligand size or charge

Changing ligand size or charge changes the ligand-water decoupling free energy.
Charged ligands usually require more windows and careful treatment of
long-range electrostatics.

### 22. Flexible docking

Flexible side chains can improve the interaction score at the cost of strain.
The best result balances interaction energy and the harmonic strain penalty.

### 23. Backend benchmark

The CPU and optional accelerated Ewald paths should agree within numerical
tolerance. Compare checksum, energy, force, and wall-clock time.

## 3. Individual tutorial self-checks

Tutorials 01–10 already include detailed self-check answers at the end of each
file. Use those for the one-particle and molecular-mechanics foundations.

## 4. What a good answer looks like

A strong answer:

- states the physical principle;
- gives an order-of-magnitude estimate;
- identifies the relevant test or invariant;
- distinguishes an exact identity from a numerical expectation;
- states the model limitations.

For example, “the energy went down” is weaker than “the accepted steps lowered
the energy and the maximum force fell below $10^{-3}$ kJ/(mol nm), which is the
stated convergence criterion.”

**Related:** [tutorial](tutorial.md),
[troubleshooting](troubleshooting.md),
[visualization prompts](visualization_prompts.md), and
[validation](../validation.md).
