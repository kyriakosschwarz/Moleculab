# Introduction to Drug Discovery Physics

**User:** Anonymous  
**Created:** 11/7/2026 8:36:20  
**Updated:** 11/9/2026 16:45:28  
**Exported:** 11/9/2026 16:52:30  
**Link:** [https://chatgpt.com/c/6aeed4f4-267c-83eb-890c-285bcea45d16](https://chatgpt.com/c/6aeed4f4-267c-83eb-890c-285bcea45d16)  

## Prompt:
11/7/2026, 8:37:07 AM

introduction to the physics and calculations of drug discovery applications like FF, FE, MD, DFT, etc.

## Response:
11/7/2026, 8:37:10 AM · gpt-6

> **Building a molecular-modeling map**
>
> I’m organizing the introduction around what each method represents, what it calculates, and how the methods connect. I’ll include core equations, practical workflows, assumptions, and limitations, then suggest a small worked example for practice.
>
> **Built a molecular-modeling map**
>
> 
>
> Worked for 40s

# Introduction to the Physics and Calculations of Computational Drug Discovery

Computational drug discovery combines **physics, chemistry, statistical mechanics, quantum mechanics, and numerical methods** to understand how molecules interact with biological targets, predict binding affinity, and design better drug candidates.

The most important computational approaches you'll encounter are\:

- **FF** - Force Fields
- **MD** - Molecular Dynamics
- **FE** - Free Energy calculations
- **DFT** - Density Functional Theory
- **QM/MM** - Quantum Mechanics / Molecular Mechanics
- **Docking** - Prediction of ligand binding poses
- **Enhanced Sampling** - Methods for exploring molecular conformations and rare events

These methods are closely related, but they describe different aspects of molecular physics.

Let's build the foundation from the physics upward.

# 1\. The physical foundations of drug discovery

## 1.1 What are we trying to calculate?

Imagine a small molecule (ligand) binding to a protein.

The fundamental question is\:

**How favorable is the interaction between the ligand and the protein?**

Several physical contributions determine the answer\:

1. Electrostatic interactions.
2. Van der Waals interactions.
3. Hydrogen bonding.
4. Hydrophobic effects.
5. Solvent interactions.
6. Entropic contributions.
7. Protein conformational changes.

The binding process can be represented as\:

$$
P + L \rightleftharpoons PL
$$

where\:

- $P$\: protein
- $L$\: ligand
- $PL$\: protein-ligand complex

The thermodynamic quantity of interest is the binding free energy\:

$$
\boxed{\Delta G_{\mathrm{bind}} =
G_{PL}-G_P-G_L}
$$

A more negative binding free energy generally indicates more favorable binding.

Importantly, binding is not determined simply by how strongly two molecules attract each other. Solvent rearrangement, molecular flexibility, and entropy also matter.

## 1.2 Thermodynamics: Energy versus Free Energy

One of the most important concepts in computational chemistry is the distinction between energy and free energy.

The Gibbs free energy is\:

$$
\boxed{G=H-TS}
$$

where\:

- $G$\: Gibbs free energy
- $H$\: enthalpy
- $T$\: absolute temperature
- $S$\: entropy

For a binding process\:

$$
\Delta G=\Delta H-T\Delta S
$$

This means binding can be favorable because of\:

- Favorable interactions (enthalpy).
- Favorable changes in disorder (entropy).
- A combination of both.

For example, a ligand may form strong hydrogen bonds with a protein but lose substantial conformational entropy upon binding.

Therefore, strong interactions do not automatically imply strong binding.

## 1.3 Statistical mechanics: Connecting physics to probability

Molecules are constantly moving. They do not have one fixed structure or one fixed energy.

Statistical mechanics describes the probability of observing different molecular configurations.

The Boltzmann distribution is\:

$$
\boxed{P_i=\frac{e^{-E_i/k_BT}}{Z}}
$$

where\:

- $P_i$\: probability of configuration $i$
- $E_i$\: energy of configuration $i$
- $k_B$\: Boltzmann constant
- $T$\: temperature
- $Z$\: partition function

The partition function is\:

$$
Z=\sum_i e^{-E_i/k_BT}
$$

For continuous molecular coordinates\:

$$
Z=\int e^{-\beta U(\mathbf{x})}\,d\mathbf{x}
$$

with\:

$$
\beta=\frac{1}{k_BT}
$$

This equation is central to molecular dynamics and free energy calculations.

**The key insight\:** Drug binding is fundamentally a statistical-mechanical problem. We are interested in the probability of molecular states, not just the energy of a single structure.

---

# 2\. FF - Molecular Mechanics and Force Fields

(Foundation)

A force field (FF) is a mathematical model describing the potential energy of a molecular system.

It approximates how atoms interact without explicitly calculating their electronic structure.

The potential energy is commonly written as\:

$$
U_{\mathrm{total}}=
U_{\mathrm{bonded}}+
U_{\mathrm{nonbonded}}
$$

## 2.1 Bonded interactions

These describe interactions between atoms connected through chemical bonds.

**Bond stretching**

$$
U_{\mathrm{bond}}=
\frac{1}{2}k_b(r-r_0)^2
$$

- $r$\: current bond length
- $r_0$\: equilibrium bond length
- $k_b$\: force constant

**Angle bending**

$$
U_{\mathrm{angle}}=
\frac{1}{2}k_\theta(\theta-\theta_0)^2
$$

**Torsional interactions**

$$
U_{\mathrm{torsion}}=
\sum_n \frac{V_n}{2}
[1+\cos(n\phi-\delta_n)]
$$

These terms describe molecular flexibility and rotation around chemical bonds.

## 2.2 Nonbonded interactions

These describe interactions between atoms that are not directly connected by chemical bonds.

Two particularly important contributions are\:

**Lennard-Jones potential (van der Waals)**

$$
U_{\mathrm{LJ}}(r)=
4\epsilon
\left[
\left(\frac{\sigma}{r}\right)^{12}
-
\left(\frac{\sigma}{r}\right)^6
\right]
$$

The two terms describe\:

- Short-range repulsion.
- Longer-range dispersion attraction.

**Coulomb electrostatic potential**

$$
U_{\mathrm{Coulomb}}(r)=
\frac{1}{4\pi\epsilon_0}
\frac{q_iq_j}{r}
$$

where $q_i$ and $q_j$ are atomic partial charges.

Electrostatic interactions are particularly important for\:

- Hydrogen bonding.
- Salt bridges.
- Charged amino acids.
- Ligand-protein recognition.

## 2.3 Common force fields

| Force field | Typical applications |
|--------|--------|
| AMBER | Proteins, nucleic acids |
| CHARMM | Proteins, lipids, membranes |
| OPLS | Proteins and organic molecules |
| GROMOS | Biomolecular simulations |
| GAFF | Small organic molecules |
| OpenFF | Small-molecule parameterization |

Some important software packages include\:

- GROMACS - Molecular simulations.
- AMBER - Biomolecular simulation.
- OpenMM - Molecular simulation framework.
- CHARMM - Molecular simulation and force fields.

**Important limitation\:** Force fields are approximations. Their quality depends on parameterization, and conventional force fields generally do not describe chemical bond breaking and formation.

---

# 3\. MD - Molecular Dynamics

(Molecular motion)

Molecular Dynamics simulates the movement of atoms over time by solving Newton's equations of motion.

The fundamental equation is\:

$$
\boxed{F_i=m_i a_i}
$$

The force on each atom is obtained from the potential energy\:

$$
\boxed{F_i=-\nabla_i U}
$$

Combining these equations\:

$$
m_i\frac{d^2\mathbf{r}_i}{dt^2}
=
-\nabla_i U
$$

This is the central equation of classical molecular dynamics.

## 3.1 How an MD simulation works

The basic algorithm is\:

1. Start with an initial molecular structure.
2. Calculate the potential energy.
3. Calculate forces on all atoms.
4. Integrate Newton's equations.
5. Update atomic positions and velocities.
6. Repeat for millions of time steps.

Numerical integration methods include\:

- Verlet.
- Velocity Verlet.
- Leapfrog.

A typical MD timestep is around 1-2 femtoseconds.

For example\:

$$
1\ \mathrm{fs}=10^{-15}\ \mathrm{s}
$$

A simulation lasting 100 ns therefore requires approximately 50-100 million integration steps, depending on the timestep.

## 3.2 What do we learn from MD?

MD allows us to study\:

- Protein flexibility.
- Ligand stability in binding pockets.
- Hydrogen bond persistence.
- Water molecules in binding sites.
- Conformational changes.
- Protein-ligand interactions.
- Structural fluctuations.

One common analysis is RMSD\:

$$
RMSD=
\sqrt{
\frac{1}{N}
\sum_{i=1}^{N}
|\mathbf r_i-\mathbf r_i^{ref}|^2
}
$$

RMSD measures structural deviation from a reference structure.

Other important analyses include\:

| Analysis | What it measures |
|--------|--------|
| RMSD | Structural deviation |
| RMSF | Atomic or residue flexibility |
| Radius of gyration | Molecular compactness |
| Hydrogen bonds | Interaction persistence |
| SASA | Solvent-accessible surface area |
| PCA | Major collective motions |

**Important\:** A stable-looking MD trajectory does not prove that a ligand binds strongly. MD explores molecular motion, but binding affinity requires appropriate thermodynamic analysis.

---

# 4\. FE - Free Energy Calculations

(Binding affinity prediction)

Free Energy calculations are among the most important methods for quantitative computational drug discovery.

They aim to estimate how favorable a molecular process is thermodynamically.

The fundamental relationship is\:

$$
\boxed{\Delta G=-k_BT\ln K}
$$

Here, $K$ is a dimensionless equilibrium constant appropriate to the process.

For ligand binding, the standard binding free energy is\:

$$
\boxed{\Delta G^\circ_{\mathrm{bind}}=
RT\ln\left(\frac{K_d}{C^\circ}\right)}
$$

where\:

- $K_d$\: dissociation constant.
- $C^\circ$\: standard concentration, usually 1 M.
- $R$\: gas constant.

Lower $K_d$ means stronger binding and a more negative binding free energy.

## 4.1 Molecular dynamics versus free energy

This is a critical distinction\:

- MD generates molecular configurations.
- FE calculations use statistical mechanics to estimate differences in thermodynamic stability.

MD is often part of an FE calculation, but they are not the same thing.

## 4.2 Main free energy methods

**A. Free Energy Perturbation (FEP)**

FEP estimates the free energy difference between two molecular states.

$$
\boxed{
\Delta G_{A\rightarrow B}
=
-k_BT
\ln
\left\langle
e^{-\beta(U_B-U_A)}
\right\rangle_A
}
$$

This is the Zwanzig equation.

It calculates the free energy change when a system is transformed from state A into state B.

For example\:

- Ligand A → Ligand B.
- Changing a substituent.
- Replacing one chemical group with another.

**B. Thermodynamic Integration (TI)**

TI introduces a coupling parameter $\lambda$ to gradually transform one molecular state into another.

$$
\boxed{
\Delta G=
\int_0^1
\left\langle
\frac{\partial U(\lambda)}
{\partial \lambda}
\right\rangle_\lambda d\lambda
}
$$

The system is simulated at several values of $\lambda$, and the resulting ensemble averages are integrated numerically.

**C. MM/PBSA and MM/GBSA**

These methods estimate binding free energies using molecular mechanics and implicit-solvent models.

A simplified expression is\:

$$
\Delta G_{\mathrm{bind}}
\approx
\Delta E_{\mathrm{MM}}
+\Delta G_{\mathrm{solv}}
-T\Delta S
$$

They are often less computationally demanding than rigorous alchemical methods, but their accuracy depends strongly on the system and approximations.

**D. Umbrella Sampling**

This method uses biasing potentials to sample regions of conformational space that would otherwise be difficult to visit.

It can be used to calculate a potential of mean force (PMF), for example along a ligand dissociation coordinate.

## 4.3 Relative versus absolute free energy

Two common approaches are\:

| Method | Objective |
|--------|--------|
| Relative binding free energy (RBFE) | Compare the binding affinity of related ligands |
| Absolute binding free energy (ABFE) | Estimate the binding affinity of an individual ligand |

RBFE is particularly useful in lead optimization, where medicinal chemists modify an existing compound and want to know whether the modification improves binding.

---

---

# 5\. DFT - Density Functional Theory

(Quantum mechanics)

DFT is fundamentally different from classical molecular mechanics.

While force fields approximate molecular interactions using parameterized equations, DFT uses quantum mechanics to describe the electronic structure of molecules.

It is based on the idea that the electronic energy of a system can be expressed as a functional of its electron density.

## 5.1 Why do we need quantum mechanics?

Atoms consist of nuclei and electrons.

When molecules interact, electrons determine\:

- Chemical bonding.
- Charge distribution.
- Electron transfer.
- Polarization.
- Reaction mechanisms.
- Bond breaking and formation.

Classical force fields cannot accurately describe all these phenomena.

The electronic structure problem begins with the Schrödinger equation\:

$$
\boxed{\hat H\Psi=E\Psi}
$$

where\:

- $\hat H$\: Hamiltonian operator.
- $\Psi$\: electronic wavefunction.
- $E$\: energy eigenvalue.

For a many-electron system, solving the full Schrödinger equation directly is extremely expensive.

DFT provides a more computationally practical framework.

## 5.2 The central concept of DFT

Instead of explicitly solving for the many-electron wavefunction, DFT uses the electron density\:

$$
\rho(\mathbf r)
$$

The ground-state energy is expressed as\:

$$
E[\rho]=
T[\rho]+
V_{\mathrm{ext}}[\rho]+
E_{\mathrm{H}}[\rho]+
E_{\mathrm{xc}}[\rho]
$$

where\:

- $T[\rho]$\: electronic kinetic energy.
- $V_{\mathrm{ext}}[\rho]$\: electron-nuclear interaction.
- $E_{\mathrm H}[\rho]$\: electron-electron electrostatic energy.
- $E_{\mathrm{xc}}[\rho]$\: exchange-correlation energy.

The exchange-correlation functional is particularly important because it accounts for complicated many-electron effects.

In practical calculations, the Kohn-Sham formulation converts the problem into a set of self-consistent equations\:

$$
\left[
-\frac{\hbar^2}{2m_e}\nabla^2+
V_{\mathrm{eff}}(\mathbf r)
\right]\psi_i
=
\epsilon_i\psi_i
$$

The electron density is reconstructed from the occupied orbitals\:

$$
\rho(\mathbf r)=
\sum_i f_i|\psi_i(\mathbf r)|^2
$$

where $f_i$ represents orbital occupation.

## 5.3 How is DFT used in drug discovery?

DFT is particularly useful when electronic effects are important.

Examples include\:

1. Calculating molecular geometries.
2. Optimizing ligand structures.
3. Calculating partial atomic charges.
4. Studying protonation and tautomerism.
5. Investigating hydrogen bonding.
6. Calculating molecular electrostatic potentials.
7. Studying reaction mechanisms.
8. Understanding metal coordination in metalloproteins.

For example, consider a ligand that binds to a zinc-containing enzyme.

The interaction between the ligand and zinc may involve electronic effects that are poorly represented by a simple classical force field.

DFT can provide a more accurate description of the coordination chemistry.

## 5.4 Common DFT functionals

| Functional | General characteristics |
|--------|--------|
| PBE | General-purpose GGA functional |
| B3LYP | Hybrid functional, widely used in molecular chemistry |
| PBE0 | Hybrid functional |
| M06-2X | Often used for molecular thermochemistry |
| ωB97X-D | Range-separated hybrid with dispersion treatment |

Common software\:

- ORCA - Quantum chemistry.
- Gaussian - Quantum chemistry.
- CP2K - Electronic structure and molecular simulations.
- Quantum ESPRESSO - Plane-wave electronic structure calculations.

**Important limitation\:** DFT is not exact in practical applications. Results depend on the functional, basis set, treatment of dispersion, solvation model, and convergence settings.

---

# 6\. QM/MM - Combining Quantum and Classical Mechanics

(Multiscale modelling)

QM/MM combines the advantages of DFT and molecular mechanics.

The central idea is simple\:

- Treat the chemically important region with quantum mechanics.
- Treat the surrounding environment with classical mechanics.

For example, imagine an enzyme containing a ligand in its active site.

You might treat\:

- Ligand: QM.
- Catalytic amino acids: QM.
- Metal ion: QM.
- Remaining protein: MM.
- Solvent: MM.

The total energy can be approximated as\:

$$
E_{\mathrm{QM/MM}}=
E_{\mathrm{QM}}+
E_{\mathrm{MM}}+
E_{\mathrm{QM-MM}}
$$

The interaction term accounts for coupling between the quantum and classical regions.

### Why is this useful?

QM/MM is especially valuable for\:

- Enzyme reaction mechanisms.
- Proton transfer.
- Covalent inhibitor mechanisms.
- Metal-containing enzymes.
- Electronic polarization and charge transfer.
- Chemical transformations in active sites.

It is computationally more expensive than conventional MD but often substantially cheaper than treating the entire protein using quantum mechanics.

---

# 7\. Docking - Predicting How a Ligand Binds

Molecular docking attempts to predict the position and orientation of a ligand inside a protein binding pocket.

The main components are\:

1. Search for possible ligand conformations.
2. Place the ligand in the binding pocket.
3. Evaluate the resulting poses.
4. Rank poses using a scoring function.

A simplified scoring function might be\:

$$
S=
w_1E_{\mathrm{vdW}}+
w_2E_{\mathrm{elec}}+
w_3E_{\mathrm{HB}}+
w_4E_{\mathrm{hydrophobic}}
$$

where $w_i$ are weighting factors.

Common docking programs include\:

- AutoDock Vina\.
- Glide\.
- GOLD\.
- DOCK\.

Docking is computationally efficient and useful for virtual screening.

However, a docking score is not equivalent to a rigorous binding free energy.

It is primarily a ranking approximation.

---

# 8\. How all these methods fit together

This is perhaps the most important part of understanding computational drug discovery.

Each technique answers a different question.

| Technique | Main question | Physical foundation |
|--------|--------|--------|
| FF | What is the potential energy of this molecular configuration? | Classical mechanics |
| MD | How do atoms move over time? | Newtonian mechanics |
| FE | How favorable is a molecular transformation or binding process? | Statistical mechanics |
| DFT | How are electrons distributed and how do they interact? | Quantum mechanics |
| QM/MM | How do electronic effects influence a chemical process in a molecular environment? | Quantum + classical mechanics |
| Docking | How might a ligand fit into a binding pocket? | Search and approximate scoring |
| Enhanced sampling | How can we explore rarely visited molecular configurations? | Statistical mechanics |

## 8.1 The computational drug discovery workflow

A typical workflow looks like this\:

(Step 1)

**Protein and ligand preparation**

Protein structures, protonation, ligand structures and parameterization

(Step 2)

**Molecular docking**

Generate and rank possible binding poses

(Step 3)

**Force field parameterization**

Define molecular interactions and potential energy functions

(Step 4)

**Molecular Dynamics**

Explore protein-ligand motion and conformational flexibility

(Step 5)

**Free Energy Calculations**

Estimate binding affinity and compare ligand modifications

(Optional advanced calculations)

**DFT / QM-MM**

Investigate electronic interactions, catalytic mechanisms and chemical reactivity

**Drug candidate prioritization**

Combine computational predictions with experimental validation

This is a conceptual workflow rather than a mandatory sequence. For example, DFT calculations may be performed during ligand parameterization, before MD, or independently to investigate a particular interaction.

# 9\. A numerical example: What does binding affinity mean?

Let's connect the physics to a real drug discovery question.

Suppose you have two ligands\:

- Ligand A: $K_d=1\ \mu M$
- Ligand B: $K_d=10\ nM$

Which one binds more strongly?

We calculate\:

$$
\Delta G^\circ=RT\ln(K_d/C^\circ)
$$

At 298 K\:

$$
RT\approx 0.592\ \mathrm{kcal/mol}
$$

For ligand A\:

$$
\Delta G_A\approx -8.18\ \mathrm{kcal/mol}
$$

For ligand B\:

$$
\Delta G_B\approx -10.91\ \mathrm{kcal/mol}
$$

Therefore\:

$$
\Delta\Delta G=\Delta G_B-\Delta G_A
$$

$$
\boxed{\Delta\Delta G=-2.73\ \mathrm{kcal/mol}}
$$

Ligand B binds approximately 100 times more strongly than ligand A.

**The important physical insight**

A difference of only approximately 1.36 kcal/mol corresponds to a tenfold difference in binding affinity at room temperature.

This explains why relatively small changes in molecular interactions can have significant pharmacological consequences.

# 10\. Important computational concepts to understand

Beyond the individual methods, there are several mathematical concepts that connect almost everything in computational drug discovery.

| Mathematical concept | Application |
|--------|--------|
| Numerical integration | Molecular dynamics trajectories |
| Differential equations | Newton's equations of motion |
| Linear algebra | Quantum chemistry and eigenvalue problems |
| Probability distributions | Molecular configurations |
| Monte Carlo methods | Statistical sampling |
| Optimization | Geometry optimization and energy minimization |
| Numerical differentiation | Force calculations |
| Thermodynamic integration | Free energy calculations |
| Fourier transforms | Molecular analysis and computational methods |
| Bayesian statistics | Uncertainty and predictive modelling |

You should also become familiar with\:

- Potential energy surfaces (PES).
- Free energy landscapes.
- Conformational sampling.
- Periodic boundary conditions.
- Temperature and pressure coupling.
- Equilibration and convergence.
- Statistical uncertainty and error estimation.

# 11\. Recommended software ecosystem

Here is a practical selection of tools to learn.

**GROMACS**

Molecular dynamics, energy minimization, trajectory analysis and simulation.

Website: [gromacs.org](https://www.gromacs.org/)

---

**OpenMM**

Programmable molecular simulations, force fields and GPU acceleration.

[openmm.org](https://openmm.org/)

---

**RDKit**

Molecular structures, cheminformatics, molecular descriptors and ligand preparation.

[rdkit.org](https://www.rdkit.org/)

---

**ORCA**

DFT, quantum chemistry, molecular properties and electronic structure.

[faccts.de/orca](https://www.faccts.de/orca/)

---

**AutoDock Vina**

Molecular docking and virtual screening.

[vina.scripps.edu](https://vina.scripps.edu/)

---

**PLUMED**

Enhanced sampling, collective variables and free energy analysis.

[plumed.org](https://www.plumed.org/)

# 12\. Suggested learning sequence

If your objective is to understand not just how to run these programs, but also the physics and mathematics behind their calculations, I recommend this sequence\:

| Stage | Subject | What you should be able to do |
|--------|--------|--------|
| 1 | Classical mechanics | Understand forces and potential energy |
| 2 | Thermodynamics | Understand energy, entropy and free energy |
| 3 | Statistical mechanics | Understand molecular ensembles and Boltzmann distributions |
| 4 | Force fields | Understand molecular energy functions |
| 5 | Molecular dynamics | Run and analyze MD simulations |
| 6 | Quantum mechanics | Understand Schrödinger and Kohn-Sham equations |
| 7 | DFT | Calculate molecular electronic properties |
| 8 | Free energy methods | Understand FEP, TI, RBFE and ABFE |
| 9 | QM/MM | Study chemical processes in proteins |

My recommendation is to learn **statistical mechanics before advanced free energy calculations**, because it provides the theoretical foundation for understanding why MD and FE work.

---

# 13\. Check your understanding

Here are five questions to test whether the fundamental concepts are clear.

**Knowledge check** (/5 answered)

[Check answers]

**A particularly useful next step\:** Work through one complete protein-ligand example, starting from the molecular potential energy function, deriving the MD equations, and then calculating binding free energy using statistical mechanics. That connects FF, MD, FE and DFT into one coherent understanding.

## Prompt:
11/7/2026, 8:44:54 AM

Work through one complete protein-ligand example, with a well known commercial drug, its target and the associated diseases (explaining the biology as well), using accessible data and free software end to end

## Response:
11/7/2026, 8:45:04 AM · gpt-6

> **Refining the workflow**
>
> I’m continuing the structural research and checking an important caveat: 1IEP is a mouse c-Abl kinase complex, so I’ll present it as a structural proxy rather than human BCR-ABL1.
>
> **Refined the workflow**
>
> 
>
> **Building a reproducible pipeline**
>
> I’m extending the plan with free software, downloadable structures, and command examples. I’ll distinguish docking, molecular dynamics, quantum calculations, and binding-energy estimates, while noting that rigorous free-energy comparisons require matched ligand series and replicated simulations.
>
> **Choosing structures carefully**
>
> I’m using 2HYY for a human ABL structural reference, while noting its multiple copies require matching each protein chain with its crystallographic ligand. I’ll reserve 1IEP for Vina’s established tutorial workflow.
>
> **Integrating the analysis**
>
> I’m connecting docking, dynamics, quantum calculations, and affinity estimation into one reproducible workflow. I’ll clearly separate quick exploratory methods from rigorous free-energy calculations, which need carefully matched analogs and repeated simulations.
>
> **Built a reproducible, carefully structured analysis pipeline**
>
> 
>
> **Searched 13 websites**
>
> 
>
> Worked for 58s

# Complete computational drug discovery case study: Imatinib (Gleevec) and the BCR-ABL1 protein

Let's work through a real drug discovery example, from human disease biology to molecular interactions, force fields, molecular dynamics and free energy calculations.

I suggest using **imatinib**, marketed as Gleevec in the US and Glivec in some other countries.

It is an excellent teaching example because\:

- It is a well-known commercial anticancer drug.
- Its molecular target is well characterized.
- Its mechanism of action is understood at the molecular level.
- Experimental protein-ligand crystal structures are publicly available.
- We can use free software and publicly accessible molecular data.
- It provides a natural connection between molecular recognition, protein conformational changes and drug resistance.

## 1\. First understand the disease: Chronic Myeloid Leukemia

### 1.1 What is CML?

Chronic Myeloid Leukemia \(CML) is a cancer of blood-forming cells in the bone marrow.

Normally, bone marrow stem and progenitor cells produce blood cells under tightly regulated signaling mechanisms.

In CML, a genetic abnormality causes a signaling protein to become abnormally active, promoting uncontrolled proliferation of myeloid cells.

The key molecular abnormality is the Philadelphia chromosome.

### 1.2 The Philadelphia chromosome

The Philadelphia chromosome results from a reciprocal chromosomal translocation\:

$$
t(9;22)(q34;q11)
$$

This brings together two genes\:

- BCR: located on chromosome 22.
- ABL1: located on chromosome 9.

The result is a fusion gene\:

$$
\boxed{BCR-ABL1}
$$

This gene produces an abnormal fusion protein with constitutively active tyrosine kinase activity.

### 1.3 What does the BCR-ABL1 protein do?

The normal ABL1 protein is a non-receptor tyrosine kinase.

Tyrosine kinases transfer phosphate groups from ATP to tyrosine residues on proteins.

The general reaction is\:

$$
\mathrm{Protein-Tyr+ATP}
\rightarrow
\mathrm{Protein-Tyr-P+ADP}
$$

Phosphorylation can activate downstream signaling proteins.

The BCR-ABL1 fusion protein has abnormally persistent kinase activity.

This activates several signaling pathways, including\:

| Signaling pathway | Biological consequence |
|--------|--------|
| RAS/MAPK | Cell proliferation |
| PI3K/AKT | Cell survival |
| JAK/STAT | Growth and survival signaling |

The result is increased proliferation and survival of leukemic cells.

### 1.4 Why is BCR-ABL1 such a good drug target?

The cancer cells depend heavily on the abnormal kinase signaling.

This creates a therapeutic opportunity: inhibit the kinase and suppress the abnormal signaling.

This is an example of targeted cancer therapy.

---

## 2\. The drug: Imatinib

Imatinib is a small-molecule tyrosine kinase inhibitor.

Its commercial formulation is known as Gleevec or Glivec.

Its important molecular properties include\:

| Property | Description |
|--------|--------|
| Drug | Imatinib |
| Drug class | Tyrosine kinase inhibitor |
| Main target | BCR-ABL1 |
| Other targets | KIT, PDGFRA and PDGFRB |
| Main disease application | Chronic myeloid leukemia |
| Other applications | Certain gastrointestinal stromal tumors and other kinase-driven cancers |

It is important to distinguish the BCR-ABL1 fusion protein from normal ABL1: imatinib inhibits the kinase activity of the fusion protein, but it is not exclusively selective for the fusion protein.

### 2.1 How does imatinib work?

Imatinib binds in the ATP-binding region of the ABL kinase domain.

It preferentially recognizes an inactive, DFG-out conformation of the kinase.

Normally, ATP enters the kinase binding pocket and enables phosphorylation.

Imatinib occupies the binding region and stabilizes an inactive kinase conformation.

Consequently\:

1. ATP-dependent kinase activity is inhibited.
2. BCR-ABL1-mediated phosphorylation decreases.
3. Downstream growth and survival signaling is suppressed.
4. Leukemic cell proliferation decreases.

This is an important distinction: imatinib does not simply destroy the protein. It inhibits its catalytic activity by binding to a particular structural conformation.

### 2.2 What happens in drug resistance?

Some BCR-ABL1 mutations reduce the effectiveness of imatinib.

One famous example is\:

$$
\boxed{\mathrm{T315I}}
$$

This mutation is located in the ATP-binding region and can interfere with imatinib binding.

It is useful to remember this example because it connects molecular structure, binding free energy and clinical resistance.

---

# 3\. Obtain the actual experimental protein-ligand structure

Now we move from biology to computational chemistry.

The first major task is to obtain an experimental structure of the drug bound to its target.

## 3.1 Our starting structure: PDB 2HYY

**PDB 2HYY**

(Experimental X-ray structure)

Human ABL1 kinase domain + imatinib · 2.40 Å resolution

[Open the experimental structure](https://www.rcsb.org/structure/2HYY)

This is particularly useful because it contains the human ABL1 kinase domain bound to imatinib. It is not the entire BCR-ABL1 fusion protein, but it contains the relevant ABL1 catalytic domain. ([www2.rcsb.org](https://www2.rcsb.org/structure/2HYY), [rcsb.org](https://www.rcsb.org/structure/2HYY))

An additional useful structure is [PDB 1IEP](https://www.rcsb.org/structure/1IEP), which is the original 2.10 Å imatinib-bound c-Abl structure used in the official AutoDock Vina example. Note that it is from mouse ABL. ([rcsb.org](https://rcsb.org/structure/1IEP), [autodock-vina.readthedocs.io](https://autodock-vina.readthedocs.io/en/latest/docking_basic.html))

For our main example, use 2HYY.

## 3.2 Download the molecular data

Create a project directory\:

```bash
mkdir imatinib_project
cd imatinib_project

mkdir data structures ligands docking md dft analysis
```

Download the protein structure\:

```bash
cd data

wget https://files.rcsb.org/download/2HYY.cif
```

Download the imatinib chemical component\:

```bash
wget https://files.rcsb.org/ligands/download/STI_ideal.sdf
```

You now have\:

| File | Contents |
|--------|--------|
| `2HYY.cif` | Experimental protein-ligand crystal structure |
| `STI_ideal.sdf` | Chemical structure of imatinib |

The ligand identifier `STI` is the PDB chemical component identifier for imatinib in this structure.

**An important distinction\:** The crystal ligand coordinates show the experimentally observed binding pose. The ideal ligand coordinates are a standardized chemical representation. Keep both available because they serve different purposes.

## 3.3 Visualize the structure

Install free visualization software\:

- UCSF ChimeraX\: [Download ChimeraX](https://www.cgl.ucsf.edu/chimerax/)\.
- PyMOL\: [Open-source PyMOL](https://github.com/schrodinger/pymol-open-source)\.

Your first objective is to visualize\:

1. The ABL1 kinase domain.
2. Imatinib inside the ATP-binding pocket.
3. Amino acid residues surrounding the ligand.
4. The kinase activation loop and DFG motif.

At this stage, don't run any simulations. First understand the actual experimental structure.

---

# 4\. Molecular chemistry: Understand the ligand

Before performing docking or MD, we need to understand the chemistry of imatinib.

Its molecular formula is\:

$$
\mathrm{C}_{29}\mathrm{H}_{31}\mathrm{N}_7\mathrm O
$$

It has\:

- 29 carbon atoms.
- 31 hydrogen atoms.
- 7 nitrogen atoms.
- 1 oxygen atom.

The molecule contains several chemically distinct regions\:

- A pyridine-containing heteroaromatic region.
- A pyrimidine-containing heteroaromatic region.
- An amide group.
- An aromatic benzene ring.
- A methylpiperazine-containing region.

These groups contribute differently to ligand recognition.

For example\:

- Nitrogen atoms can act as hydrogen-bond acceptors.
- The amide group provides polar interaction possibilities.
- Aromatic groups contribute to hydrophobic and dispersion interactions.
- The piperazine group influences protonation, solubility and electrostatics.

## 4.1 Calculate basic molecular properties with RDKit

Install RDKit\:

```bash
conda install -c conda-forge rdkit
```

Run\:

```python
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

mol = Chem.SDMolSupplier(
    "data/STI_ideal.sdf",
    removeHs=False
)[0]

print("Molecular formula:",
      rdMolDescriptors.CalcMolFormula(mol))

print("Molecular weight:",
      Descriptors.MolWt(mol))

print("LogP:",
      Descriptors.MolLogP(mol))

print("H-bond donors:",
      Descriptors.NumHDonors(mol))

print("H-bond acceptors:",
      Descriptors.NumHAcceptors(mol))

print("Rotatable bonds:",
      Descriptors.NumRotatableBonds(mol))
```

These are useful molecular descriptors.

They help characterize\:

- Molecular size.
- Lipophilicity.
- Hydrogen-bonding capacity.
- Flexibility.

For drug discovery, these descriptors help interpret why some compounds behave differently in biological environments.

**Important\:** Imatinib has ionizable groups, particularly its piperazine-containing region. Protonation state matters. For the simulations, you should establish the intended protonation state at approximately physiological pH rather than assume that the neutral ideal SDF is automatically appropriate.

---

# 5\. Docking: Can we reproduce the experimentally observed binding pose?

Now we perform our first genuine computational drug discovery calculation.

We know that imatinib binds to ABL1 from the crystal structure.

Can docking independently recover a similar binding pose?

This is called *redocking*\.

## 5.1 Install AutoDock Vina

AutoDock Vina is free and open-source.

**AutoDock Vina**

[Official website](https://vina.scripps.edu/) · [Documentation](https://autodock-vina.readthedocs.io/)

Install\:

```bash
conda install -c conda-forge vina meeko
```

Meeko prepares molecular structures for docking.

The official Vina documentation includes a worked imatinib/c-Abl example, which is particularly useful for checking your setup. ([autodock-vina.readthedocs.io](https://autodock-vina.readthedocs.io/en/latest/docking_basic.html), [autodock-vina.readthedocs.io](https://autodock-vina.readthedocs.io/en/stable/))

## 5.2 Prepare the receptor and ligand

The basic workflow is\:

1. Select one ABL1 protein chain.
2. Extract its corresponding crystal ligand.
3. Remove the ligand from the receptor.
4. Add missing hydrogen atoms.
5. Assign receptor protonation states.
6. Convert the receptor to PDBQT.
7. Prepare imatinib as a ligand PDBQT.

One important point: don't simply convert the original CIF file to PDBQT without checking its chains, ligand copies, waters and other components.

For the initial exercise, use the official Vina 1IEP tutorial as the reference preparation workflow. It provides the necessary receptor and ligand preparation instructions. ([autodock-vina.readthedocs.io](https://autodock-vina.readthedocs.io/en/latest/docking_basic.html))

## 5.3 Define the docking box

The docking box must cover the original ligand-binding region.

If the original ligand has coordinates\:

$$
(x_i,y_i,z_i)
$$

calculate the ligand centroid\:

$$
x_c=\frac{1}{N}\sum_i x_i
$$

$$
y_c=\frac{1}{N}\sum_i y_i
$$

$$
z_c=\frac{1}{N}\sum_i z_i
$$

These coordinates provide the center of your docking box.

For a first experiment, choose a box that contains the complete ligand-binding pocket with sufficient surrounding space.

## 5.4 Run docking

Once the receptor and ligand are prepared\:

```bash
vina \
  --receptor receptor.pdbqt \
  --ligand imatinib.pdbqt \
  --center_x X \
  --center_y Y \
  --center_z Z \
  --size_x 25 \
  --size_y 25 \
  --size_z 25 \
  --exhaustiveness 16 \
  --out docking/imatinib_poses.pdbqt
```

Replace `X`, `Y` and `Z` with the actual pocket coordinates.

The result contains predicted ligand poses and approximate docking scores.

## 5.5 How do we evaluate the result?

We compare the predicted ligand pose against the experimental crystal pose.

The key metric is ligand RMSD\:

$$
RMSD=
\sqrt{
\frac{1}{N}
\sum_{i=1}^{N}
\left|
\mathbf r_i^{pred}
-
\mathbf r_i^{crystal}
\right|^2
}
$$

The atom correspondence and molecular alignment must be handled correctly when calculating this value.

A lower RMSD indicates better reproduction of the experimental binding pose.

A value around 2 Å or less is often used as a practical redocking benchmark, although symmetry, atom mapping and alignment choices matter.

**What you should learn\:** Docking is primarily a pose-prediction and ranking method. A good docking score does not establish a measured binding affinity.

---

# 6\. Force fields: Build the physical model

Now we move from docking to molecular mechanics.

The docking calculation gives us a possible pose. We want to understand the molecular interactions that stabilize that pose.

Our model will contain\:

- ABL1 protein.
- Imatinib.
- Water.
- Counterions where necessary.

The total potential energy is\:

$$
U=
U_{\mathrm{bonded}}
+
U_{\mathrm{electrostatic}}
+
U_{\mathrm{vdW}}
$$

We need suitable parameters for both the protein and ligand.

## 6.1 Protein force field

For the protein, use\:

**AMBER14SB**

It provides parameters for proteins, including backbone and side-chain interactions.

## 6.2 Ligand force field

For imatinib, use\:

**OpenFF Sage**

OpenFF provides small-molecule force fields designed for organic molecules.

**Open Force Field**

[Official website](https://openforcefield.org/) · [Documentation](https://docs.openforcefield.org/)

The OpenFF documentation contains a protein-ligand simulation example combining a Sage ligand model, an AMBER protein model, water and OpenMM. This is an excellent starting point for our system. ([docs.openforcefield.org](https://docs.openforcefield.org/en/latest/examples/openforcefield/openff-toolkit/toolkit_showcase/toolkit_showcase.html), [docs.openforcefield.org](https://docs.openforcefield.org/en/latest/examples.html))

## 6.3 Why do we need different force fields?

Proteins and small organic molecules have different parameterization requirements.

The protein force field describes\:

- Peptide bonds.
- Backbone angles.
- Side-chain torsions.
- Protein electrostatics.

The ligand force field describes\:

- Aromatic systems.
- Heterocycles.
- Amide geometry.
- Piperazine conformations.
- Ligand electrostatics.

Combining them requires a consistent parameterization and treatment of interactions between the two molecular components.

---

# 7\. Molecular Dynamics: Simulate the imatinib-ABL1 complex

Now we can perform the central simulation.

## 7.1 Install OpenMM

OpenMM is an open-source molecular simulation engine.

Website: [OpenMM](https://openmm.org/)

Install the main tools\:

```bash
conda install -c conda-forge openmm pdbfixer mdtraj
```

For ligand parameterization, also install the OpenFF toolkit and its simulation integration packages, following the current installation instructions\:

[OpenFF installation and tutorials](https://docs.openforcefield.org/)

## 7.2 Prepare the protein-ligand complex

The preparation process is scientifically important.

You need to\:

1. Select a single protein chain and its matching ligand copy.
2. Resolve missing residues and atoms where appropriate.
3. Determine protonation states.
4. Add hydrogens.
5. Parameterize imatinib.
6. Add explicit water.
7. Add counterions.
8. Check the system for clashes and abnormal geometry.

For a first simulation, preserve the experimentally observed ligand pose rather than using a freshly docked pose.

This gives us a known starting configuration.

## 7.3 Solvate the system

We place the complex inside a periodic water box.

The water molecules are important because protein-ligand binding occurs in an aqueous environment.

Use a TIP3P water model compatible with the selected protein force field.

Add ions to neutralize the system, and consider physiological ionic strength for more realistic conditions.

The simulation box should have sufficient padding around the protein to avoid interactions with its periodic images.

## 7.4 Minimize the energy

Before starting MD, perform energy minimization.

We are trying to remove severe atomic overlaps and unfavorable initial geometries.

Mathematically\:

$$
\boxed{
\nabla U(\mathbf r)=0
}
$$

At a local minimum, the forces are approximately zero.

In practice, minimization algorithms such as L-BFGS or conjugate gradient are used.

## 7.5 Equilibration

Next, gradually bring the system to the desired thermodynamic conditions.

A typical educational setup might be\:

| Parameter | Example |
|--------|--------|
| Temperature | 300 K |
| Pressure | 1 bar |
| Water model | TIP3P |
| Protein force field | AMBER14SB |
| Ligand force field | OpenFF Sage |
| MD timestep | 2 fs, with appropriate constraints |
| Initial production simulation | 10-20 ns |

These are illustrative starting parameters, not a validated protocol for this particular kinase system.

Typical equilibration stages\:

1. Restrained minimization.
2. NVT equilibration (constant particle number, volume and temperature).
3. NPT equilibration (constant particle number, pressure and temperature).
4. Production MD.

## 7.6 Production MD

Once the system is equilibrated, we let the protein and ligand move.

The fundamental equation is\:

$$
m_i\frac{d^2\mathbf r_i}{dt^2}
=
-\nabla_i U
$$

OpenMM calculates the forces and propagates the molecular system over time.

The OpenMM tutorial collection includes examples of system construction, equilibration, simulation and advanced sampling. ([openmm.github.io](https://openmm.github.io/openmm-cookbook/latest/tutorials.html), [openmm.github.io](https://openmm.github.io/openmm-cookbook/dev/index.html))

## 7.7 Analyze the MD trajectory

After the simulation, we want to answer a biological question\:

**Does imatinib remain in its experimentally observed binding pocket, and how does the surrounding protein behave?**

Use MDTraj\:

```bash
conda install -c conda-forge mdtraj
```

### A. RMSD

```python
import mdtraj as md
import matplotlib.pyplot as plt

traj = md.load(
    "analysis/production.dcd",
    top="analysis/system.pdb"
)

protein = traj.topology.select("protein and backbone")

rmsd = md.rmsd(
    traj,
    traj,
    atom_indices=protein
)

plt.plot(traj.time, rmsd)
plt.xlabel("Time (ps)")
plt.ylabel("RMSD (nm)")
plt.title("ABL1 protein RMSD")
plt.show()
```

For the actual analysis, align frames to an appropriate reference structure first and calculate ligand RMSD separately after aligning the protein binding-site region.

### B. RMSF

RMSF measures how much individual atoms or residues fluctuate during the simulation.

$$
RMSF_i=
\sqrt{
\left\langle
|\mathbf r_i-\langle\mathbf r_i\rangle|^2
\right\rangle
}
$$

It helps identify flexible regions of ABL1.

Pay particular attention to\:

- The activation loop.
- The ATP-binding pocket.
- The regions surrounding imatinib.

### C. Hydrogen bond analysis

We want to understand which interactions persist throughout the simulation.

A hydrogen bond is commonly identified using\:

- Donor-acceptor distance.
- Donor-hydrogen-acceptor angle.

For example, you can calculate the fraction of simulation frames in which a particular interaction exists.

$$
f_{\mathrm{HB}}=
\frac{N_{\mathrm{frames\ with\ HB}}}
{N_{\mathrm{total\ frames}}}
$$

This provides a measure of interaction persistence.

It is not itself a binding free energy.

### D. Ligand-protein distance

Track distances between important ligand atoms and nearby protein residues.

For example\:

$$
d(t)=
|\mathbf r_{\mathrm{ligand}}(t)
-\mathbf r_{\mathrm{residue}}(t)|
$$

This allows you to observe whether important interactions remain geometrically stable.

**Expected output of this stage\:**

Protein RMSD

Protein RMSF

Hydrogen-bond occupancy

Interaction distances

These analyses are complementary. A low RMSD does not prove strong binding, and a persistent hydrogen bond does not necessarily indicate favorable overall binding free energy.

---

# 8\. Free energy: How strongly does imatinib bind?

This is where the project becomes particularly interesting.

We have\:

- An experimentally determined protein-ligand structure.
- A parameterized molecular system.
- An MD trajectory.
- Information about molecular interactions.

But we still have not calculated the binding free energy.

Remember\:

$$
\Delta G_{\mathrm{bind}}=
G_{PL}-G_P-G_L
$$

## 8.1 First approach: MM/GBSA

For an introductory project, MM/GBSA is a useful next step.

It approximates the binding free energy using molecular mechanics and an implicit-solvent model.

$$
\Delta G_{\mathrm{bind}}
\approx
\Delta E_{\mathrm{MM}}
+
\Delta G_{\mathrm{GB}}
+
\Delta G_{\mathrm{SA}}
-
T\Delta S
$$

The terms represent\:

- $\Delta E_{\mathrm{MM}}$\: molecular mechanics energy difference.
- $\Delta G_{\mathrm{GB}}$\: generalized Born solvation contribution.
- $\Delta G_{\mathrm{SA}}$\: nonpolar solvent-accessible surface contribution.
- $-T\Delta S$\: entropy contribution.

For an educational calculation, one can first omit the entropy term and clearly report the result as an approximate MM/GBSA energy estimate, not a rigorous binding free energy.

A free tool for this approach is\:

gmx_MMPBSA - [Documentation](https://valdes-tresanco-ms.github.io/gmx_MMPBSA/)\.

It is designed to perform MM/PBSA and MM/GBSA calculations using molecular simulation trajectories.

It is primarily designed around GROMACS and AMBER workflows, so it is easiest to use if you generate the corresponding compatible simulation files.

## 8.2 Second approach: Free Energy Perturbation

For a more rigorous drug discovery calculation, I would eventually move to FEP.

Imagine we have two related imatinib derivatives\:

$$
L_A \longrightarrow L_B
$$

We want to calculate whether the modification improves binding.

The relative binding free energy is\:

$$
\boxed{
\Delta\Delta G_{\mathrm{bind}}=
\Delta G_{\mathrm{bind},B}
-
\Delta G_{\mathrm{bind},A}
}
$$

A thermodynamic cycle is used\:

The transformation is calculated twice\:

1. Transform ligand A into ligand B in the protein binding site.
2. Transform ligand A into ligand B in water.

Then\:

$$
\boxed{
\Delta\Delta G_{\mathrm{bind}}=
\Delta G_{\mathrm{complex}}^{A\rightarrow B}
-
\Delta G_{\mathrm{solvent}}^{A\rightarrow B}
}
$$

This is especially valuable for medicinal chemistry because it can help predict whether a chemical modification will improve binding affinity.

The OpenMM Cookbook provides tutorials on multistate sampling and alchemical free energy calculations\:

[OpenMM free energy tutorials](https://openmm.github.io/openmm-cookbook/latest/tutorials.html)\. ([openmm.github.io](https://openmm.github.io/openmm-cookbook/latest/tutorials.html), [docs.openforcefield.org](https://docs.openforcefield.org/en/latest/examples.html))

**An important scientific limitation\:** We cannot calculate a meaningful relative binding free energy from one ligand alone. We need a second chemically defined ligand, a carefully designed transformation, adequate sampling and convergence analysis.

For a proper project, I would compare imatinib against one or more experimentally characterized close analogues, rather than arbitrarily modifying the drug and treating the resulting calculation as experimental validation.

---

# 9\. DFT: Investigating the electronic properties of imatinib

Now let's introduce quantum chemistry into our example.

MD and force fields describe atomic interactions approximately.

DFT can investigate the electronic structure of imatinib itself and selected molecular fragments.

For example, we can investigate\:

- Electron density.
- Molecular electrostatic potential.
- Charge distribution.
- Hydrogen-bonding tendencies.
- Molecular orbital energies.
- The effects of protonation.

## 9.1 Free quantum chemistry software

I recommend\:

CP2K - [Official website](https://www.cp2k.org/)\.

It is an open-source electronic-structure and molecular simulation package.

Another option is Psi4 - [Official website](https://psicode.org/)\.

## 9.2 What would we calculate?

Start with the isolated imatinib molecule.

First, obtain a chemically valid 3D geometry and optimize it.

The geometry optimization seeks\:

$$
\boxed{
\frac{\partial E}{\partial R_i}=0
}
$$

for the nuclear coordinates $R_i$\.

Then perform a DFT single-point energy calculation.

For example, with CP2K, a PBE calculation is a reasonable introductory exercise, with appropriate basis sets, pseudopotentials and convergence settings.

Next, calculate and visualize the molecular electrostatic potential\:

$$
V(\mathbf r)=
\sum_A
\frac{Z_A}{|\mathbf r-\mathbf R_A|}
-
\int
\frac{\rho(\mathbf r')}
{|\mathbf r-\mathbf r'|}
d\mathbf r'
$$

This helps identify electron-rich and electron-poor regions of the molecule.

You can then connect those regions to potential hydrogen-bonding and electrostatic interactions in the ABL1 binding pocket.

**A useful advanced experiment\:** Compare the electronic properties of neutral imatinib and a protonated imatinib state. This demonstrates why protonation-state selection matters in ligand preparation and molecular simulation.

Do not interpret isolated-ligand DFT calculations as a prediction of binding affinity. To understand the full binding interaction quantum mechanically, you would eventually need an appropriate protein-ligand QM/MM model.

---

# 10\. Connect the molecular calculations to disease biology

Let's return to the original biological question.

Why does this entire computational workflow matter?

Because molecular interactions determine whether the drug can inhibit the disease-associated target.

For imatinib, there are several relevant disease applications.

| Disease | Relevant molecular target | Biological connection |
|--------|--------|--------|
| Chronic myeloid leukemia (CML) | BCR-ABL1 | Constitutively active fusion kinase |
| Philadelphia chromosome-positive ALL | BCR-ABL1 | Abnormal kinase signaling in leukemic cells |
| Gastrointestinal stromal tumors (GIST) | KIT | Abnormal KIT kinase signaling in some tumors |
| Certain myeloproliferative disorders | PDGFRA / PDGFRB | Abnormal receptor tyrosine kinase signaling |

The National Cancer Institute provides a useful overview of CML biology and targeted treatment: [NCI: Chronic Myeloid Leukemia Treatment](https://www.cancer.gov/types/leukemia/hp/cml-treatment-pdq)\. ([cancer.gov](https://www.cancer.gov/types/leukemia/hp/cml-treatment-pdq), [cancer.gov](https://www.cancer.gov/publications/dictionaries/cancer-terms/def/philadelphia-chromosome))

The important point is that imatinib is not just a molecule that fits into a protein pocket.

It inhibits particular signaling mechanisms that contribute to disease.

This is the connection between\:

**Molecular structure → Molecular interactions → Protein function → Cell signaling → Disease biology.**

---

# 11\. Complete software and data toolkit

Here is the software ecosystem I would use for this project.

All the core packages below are freely available, although individual package licenses and commercial-use conditions should be checked.

**Your computational drug discovery laboratory**

| Task | Software | Cost |
|--------|--------|--------|
| Protein structure | RCSB PDB | Free |
| Ligand structures | PubChem / PDB | Free |
| Molecular visualization | ChimeraX / PyMOL | Free versions |
| Cheminformatics | RDKit | Free |
| Docking | AutoDock Vina | Open source |
| Docking preparation | Meeko | Open source |
| Protein preparation | PDBFixer | Open source |
| Force fields | AMBER / OpenFF | Free |
| Molecular dynamics | OpenMM | Open source |
| Trajectory analysis | MDTraj | Open source |
| Free energy | OpenMM tools / gmx\_MMPBSA | Open source |
| Quantum chemistry | CP2K / Psi4 | Open source |
| Plotting | Python / Matplotlib | Open source |

---

**Programming environment**

- Python
- Jupyter Notebook
- Conda or Mamba
- Linux or WSL2

A GPU is helpful for longer MD simulations, but initial molecular preparation, docking and small educational calculations can be performed on a CPU.

# 12\. What should the final project deliver?

I would organize the work into the following deliverables.

| Deliverable | Scientific question |
|--------|--------|
| 1\. Disease biology report | Why is BCR-ABL1 important in CML? |
| 2\. Molecular structure | How does imatinib bind ABL1? |
| 3\. Ligand characterization | What are the chemical properties of imatinib? |
| 4\. Docking results | Can we reproduce the experimental binding pose? |
| 5\. Force field model | How are protein and ligand interactions represented? |
| 6\. MD trajectory | How does the complex behave over time? |
| 7\. MD analysis | Which interactions persist and which regions fluctuate? |
| 8\. MM/GBSA or FEP | What is the estimated thermodynamic contribution to binding? |
| 9\. DFT analysis | What electronic properties contribute to molecular recognition? |
| 10\. Biological interpretation | How do these findings relate to kinase inhibition? |

## 13\. A realistic four-week learning project

I would divide this into four weeks, with each week producing something tangible.

(Week 1)

**Biology, structures and docking**

- Understand CML and BCR-ABL1.
- Download PDB 2HYY.
- Visualize imatinib in the binding pocket.
- Calculate molecular descriptors with RDKit.
- Perform redocking with AutoDock Vina.

Output: A protein-ligand structural analysis and docking report.

---

(Week 2)

**Force fields and molecular dynamics**

- Understand AMBER and OpenFF.
- Parameterize the protein and ligand.
- Solvate the system.
- Perform minimization and equilibration.
- Run a short MD simulation.

Output: An MD trajectory of the imatinib-ABL1 complex.

---

(Week 3)

**Trajectory analysis and thermodynamics**

- Calculate RMSD and RMSF.
- Analyze hydrogen bonds.
- Investigate protein-ligand distances.
- Understand MM/GBSA.
- Learn the thermodynamic cycle behind FEP.

Output: Molecular dynamics analysis and an introductory binding-energy report.

---

(Week 4)

**Quantum chemistry and biological interpretation**

- Optimize imatinib geometry using DFT.
- Calculate electronic properties.
- Visualize the electrostatic potential.
- Connect electronic effects to ligand-protein recognition.
- Interpret the results in terms of kinase inhibition.

Output: A complete computational case study.

# 14\. Important scientific limitations

I want to emphasize five points so that this project teaches you real computational drug discovery rather than simply how to operate software.

1\. **A crystal structure is not a complete biological system.** The 2HYY structure contains the ABL1 kinase domain, not the entire BCR-ABL1 fusion protein or a living cell.

2\. **Docking scores are not experimental binding affinities.** Redocking validates aspects of the structural workflow, not the predictive accuracy of the scoring function.

3\. **Force field quality matters.** Incorrect protonation, ligand charges or torsional parameters can compromise an entire MD simulation.

4\. **A short MD trajectory is not proof of binding stability.** Sampling, replicas and convergence must be considered.

5\. **Free energy calculations require careful statistical sampling.** A numerical answer is not necessarily a reliable answer. Uncertainty estimates and convergence tests are essential.

---

# 15\. The learning resources I would prioritize

Rather than reading dozens of unrelated papers, I would start with these resources.

**1\. Molecular Modelling: Principles and Applications - Andrew Leach**

Excellent foundation for molecular mechanics, quantum chemistry, statistical mechanics and molecular simulation.

---

**2\. Understanding Molecular Simulation - Frenkel and Smit**

Particularly valuable for statistical mechanics, sampling, Monte Carlo and free energy methods.

---

**3\. Molecular Driving Forces - Ken Dill and Sarina Bromberg**

A particularly accessible introduction to statistical thermodynamics and molecular interactions.

---

**4\. Practical software tutorials**

- [AutoDock Vina: Imatinib docking tutorial](https://autodock-vina.readthedocs.io/en/latest/docking_basic.html)
- [OpenMM: Molecular dynamics and free energy](https://openmm.github.io/openmm-cookbook/latest/tutorials.html)
- [OpenFF: Protein-ligand simulations](https://docs.openforcefield.org/en/latest/examples.html)
- [CP2K: Quantum chemistry tutorials](https://www.cp2k.org/howto)

---

# 16\. Test yourself: choose the next computational experiment

To make this a genuinely hands-on learning experience, I would suggest choosing which part of the project to implement first.

**1\. How detailed should the practical implementation be?**

Beginner - explain every command and equation

Complete - executable code, commands and explanations

Advanced - rigorous physics, algorithms and mathematical derivations

---

**2\. Which calculation should we implement first?**

Molecular docking and binding interactions

Force fields and molecular dynamics

Free energy calculations (FEP / MMGBSA)

DFT and electronic structure

Complete integrated workflow

---

**3\. Where will you run the calculations?**

My own computer (Linux / Windows / Mac)

Google Colab

University or HPC cluster

Not sure yet

[Build my practical project]

**My recommendation\:** Start with the complete MD workflow, using the experimentally observed imatinib-ABL1 complex as your starting structure. It provides the best bridge between classical mechanics, statistical mechanics, molecular interactions and the later free energy calculations.

## Prompt:
11/9/2026, 12:15:39 PM

gradually build a small python or rust library from scratch to really understand the physics in depth, in an educational manner

## Response:
11/9/2026, 12:15:51 PM · gpt-6

> **Building the first lesson**
>
> I’m starting with Python for clarity, while keeping Rust as a later performance path. I’ll introduce units, project structure, Lennard-Jones energy and forces, then add small numerical tests before expanding toward dynamics.
>
> **Formalizing Lennard-Jones forces**
>
> I’m defining the pair potential and vector force carefully, including units and sign conventions. I’ll test the equilibrium distance, minimum energy, and force using finite differences before adding particles or integration.
>
> **Built the first lesson and formalized Lennard-Jones forces**
>
> 
>
> Worked for 38s

# Building a molecular simulation library from scratch: From physics to drug discovery

I think this is an excellent way to learn computational drug discovery. Rather than treating GROMACS, OpenMM, ORCA or AutoDock as black boxes, we'll gradually build our own small molecular simulation engine and understand exactly what happens inside it.

We'll use **Python initially**, with NumPy for numerical calculations and Matplotlib for visualization. Later, we can rewrite performance-critical components in Rust.

The goal is not to reproduce a commercial molecular simulation package. It is to build a small, scientifically meaningful library whose components implement the same fundamental physical principles.

Our eventual application will be the imatinib-ABL1 protein-ligand system.

But we'll start with two atoms, derive the physics, write the code, test it against analytical results, and progressively build up to a protein-ligand simulation.

## 1\. What are we going to build?

Let's call our educational library\:

# Moleculab

(Educational project)

**A molecular physics and simulation engine**

A Python library that gradually develops capabilities for\:

- Molecular structures and chemical topology.
- Molecular mechanics and force fields.
- Energy and force calculations.
- Molecular dynamics.
- Thermodynamics and statistical mechanics.
- Molecular interactions and docking.
- Free energy calculations.
- Quantum chemistry and DFT concepts.
- Protein-ligand simulation.
  
  Working project name; not an existing validated scientific package.

Our design principle will be\:

**Every important physical equation should have a corresponding implementation, unit test and numerical experiment.**

For example, when implementing the Lennard-Jones potential, we will not just write a function. We will\:

1. Derive the potential mathematically.
2. Understand its physical meaning.
3. Derive its force analytically.
4. Implement both energy and force.
5. Compare the analytical force with a numerical derivative.
6. Plot the potential and force.
7. Investigate how changing its parameters changes molecular behavior.

This is how we will approach every major component.

## 2\. The architecture of our library

Eventually, our library will have the following architecture\:

**Moleculab - Public API**

Create molecules, simulate, analyze and visualize

**Structure**

Atoms, bonds, residues, molecules

**Force fields**

Bonds, angles, torsions, nonbonded

**Physics**

Energy, forces, integration

**Simulation**

MD, thermostats, barostats

**Thermodynamics**

Ensembles, free energy, sampling

**Quantum chemistry**

Orbitals, electron density, DFT

**Numerical foundation**

Linear algebra · Numerical integration · Optimization · Statistics

A critical architectural decision: we will separate the physical model from the numerical algorithms.

For example\:

- Lennard-Jones describes a physical interaction model.
- Velocity Verlet is a numerical integration algorithm.
- NVT describes a thermodynamic ensemble.
- A thermostat is an algorithm for sampling that ensemble.

Keeping these concepts separate makes the library easier to understand and extend.

## 3\. Development roadmap

We'll build the library in small, testable milestones.

| Milestone | Implementation | Physics learned |
|--------|--------|--------|
| 0 | Units and numerical conventions | Dimensional analysis |
| 1 | Atoms and coordinates | Classical mechanics |
| 2 | Lennard-Jones potential | Interatomic interactions |
| 3 | Coulomb potential | Electrostatics |
| 4 | Molecular topology | Chemical connectivity |
| 5 | Bond and angle potentials | Molecular mechanics |
| 6 | Torsional potentials | Molecular conformations |
| 7 | Energy minimization | Optimization |
| 8 | Velocity Verlet | Numerical integration |
| 9 | Molecular dynamics | Newtonian mechanics |
| 10 | Temperature and thermostats | Statistical mechanics |
| 11 | Solvent and periodic boundaries | Condensed matter physics |
| 12 | Trajectory analysis | Statistical analysis |
| 13 | Protein-ligand interactions | Molecular recognition |
| 14 | Free energy | Statistical thermodynamics |
| 15 | FEP and thermodynamic integration | Alchemical transformations |
| 16 | Basic quantum chemistry | Schrödinger equation |
| 17 | Hartree-Fock and DFT foundations | Electronic structure |
| 18 | QM/MM concepts | Multiscale physics |

We won't rush through the milestones. Each one should have a working implementation before we move on.

### Python first, Rust later

I recommend this division\:

| Python | Rust |
|--------|--------|
| Derivations and experimentation | High-performance calculations |
| NumPy numerical arrays | Efficient memory management |
| Matplotlib plots | Parallel computation |
| Jupyter notebooks | Production-oriented engine |
| Rapid prototyping | Optimized numerical kernels |

Eventually, we can expose Rust functions to Python using PyO3.

For now, Python makes it much easier to inspect intermediate calculations and understand the mathematics.

## 4\. Establish our physical units

This is one of the most important decisions in computational chemistry.

We will use the following unit system\:

| Quantity | Unit |
|--------|--------|
| Distance | nm |
| Time | ps |
| Energy | kJ/mol |
| Mass | Da |
| Force | kJ/(mol nm) |
| Temperature | K |
| Charge | Elementary charge $e$ |

Why these units?

They are widely used in biomolecular simulation and are convenient for molecular dynamics.

For example\:

$$
1\ \mathrm{nm}=10^{-9}\ \mathrm m
$$

$$
1\ \mathrm{ps}=10^{-12}\ \mathrm s
$$

$$
1\ \mathrm{Da}\approx1.66054\times10^{-27}\ \mathrm{kg}
$$

An interesting property of this unit system is that the numerical Newtonian acceleration equation becomes particularly convenient\:

$$
\boxed{a_i=\frac{F_i}{m_i}}
$$

when force, mass and acceleration are expressed in our selected simulation units.

We will still implement explicit unit checks and conversion utilities later, because mixing units is a common source of serious simulation errors.

---

# 5\. Milestone 1: The first molecular interaction

Let's start with two atoms.

Imagine two neutral carbon-like atoms separated by a distance $r$\.

They interact through two principal physical effects\:

1. Attraction due to dispersion interactions.
2. Repulsion when their electron clouds overlap strongly.

The Lennard-Jones potential models these effects.

## 5.1 Deriving the Lennard-Jones potential

The potential is\:

$$
\boxed{
U(r)=4\epsilon
\left[
\left(\frac{\sigma}{r}\right)^{12}
-
\left(\frac{\sigma}{r}\right)^6
\right]
}
$$

There are two parameters\:

- $\epsilon$\: depth of the potential well.
- $\sigma$\: distance at which the potential crosses zero.

Let's understand each term.

**Attractive term\:**

$$
U_{\mathrm{attraction}}(r)=
-4\epsilon\left(\frac{\sigma}{r}\right)^6
$$

This approximates London dispersion interactions.

**Repulsive term\:**

$$
U_{\mathrm{repulsion}}(r)=
4\epsilon\left(\frac{\sigma}{r}\right)^{12}
$$

This approximates short-range repulsion.

The twelfth power is a computationally convenient approximation to the much steeper repulsion associated with electron-cloud overlap.

### What does the potential look like?

The minimum occurs at\:

$$
\boxed{r_{\min}=2^{1/6}\sigma}
$$

At that distance\:

$$
\boxed{U(r_{\min})=-\epsilon}
$$

This gives us two excellent analytical tests for our implementation.

## 5.2 Derive the force

Newton's equations require forces, not just energies.

The force is the negative derivative of potential energy\:

$$
F(r)=-\frac{dU(r)}{dr}
$$

Starting from\:

$$
U(r)=4\epsilon
\left[
\sigma^{12}r^{-12}
-
\sigma^6r^{-6}
\right]
$$

Differentiate\:

$$
\frac{dU}{dr}=
4\epsilon
\left[
-12\sigma^{12}r^{-13}
+
6\sigma^6r^{-7}
\right]
$$

Therefore\:

$$
\boxed{
F(r)=
\frac{24\epsilon}{r}
\left[
2\left(\frac{\sigma}{r}\right)^{12}
-
\left(\frac{\sigma}{r}\right)^6
\right]
}
$$

This is the radial force.

To obtain the three-dimensional vector force, define\:

$$
\mathbf r_{ij}=\mathbf r_i-\mathbf r_j
$$

Then\:

$$
\boxed{
\mathbf F_i=
\frac{24\epsilon}{r^2}
\left[
2\left(\frac{\sigma}{r}\right)^{12}
-
\left(\frac{\sigma}{r}\right)^6
\right]
\mathbf r_{ij}
}
$$

And Newton's third law gives\:

$$
\mathbf F_j=-\mathbf F_i
$$

Notice something important\:

- When atoms are too close, they repel.
- When atoms are sufficiently far apart, they attract.
- At the equilibrium separation, the force is zero.

This is the first physical model we will implement.

---

## 5.3 Implement the potential in Python

Let's create the project structure\:

```text
moleculab/
│
├── pyproject.toml
├── README.md
│
├── src/
│   └── moleculab/
│       ├── __init__.py
│       ├── units.py
│       ├── potentials.py
│       ├── particles.py
│       └── system.py
│
├── tests/
│   ├── test_potentials.py
│   └── test_forces.py
│
└── notebooks/
    └── 01_lennard_jones.ipynb
```

For the first milestone, we only need `models/potentials.py`\.

Install our numerical dependencies\:

```bash
pip install numpy matplotlib pytest
```

Create `src/moleculab/models/potentials.py`\:

```python
import numpy as np

def lennard_jones_energy(r, epsilon, sigma):
    """
    Lennard-Jones potential energy.

    Parameters
    ----------
    r : float
        Distance between atoms (nm).

    epsilon : float
        Potential well depth (kJ/mol).

    sigma : float
        Zero-crossing distance (nm).

    Returns
    -------
    float
        Potential energy (kJ/mol).
    """

    if r <= 0:
        raise ValueError("Distance must be positive.")

    if epsilon <= 0 or sigma <= 0:
        raise ValueError("epsilon and sigma must be positive.")

    sr6 = (sigma / r) ** 6

    energy = 4.0 * epsilon * (sr6**2 - sr6)

    return energy

def lennard_jones_force(displacement, epsilon, sigma):
    """
    Vector force on atom i due to atom j.

    displacement = position_i - position_j

    Units:
        kJ/(mol nm)
    """

    displacement = np.asarray(displacement, dtype=float)

    r = np.linalg.norm(displacement)

    if r <= 0:
        raise ValueError("Atoms cannot occupy identical positions.")

    if epsilon <= 0 or sigma <= 0:
        raise ValueError("epsilon and sigma must be positive.")

    sr6 = (sigma / r) ** 6

    coefficient = (
        24.0 * epsilon / r**2
        * (2.0 * sr6**2 - sr6)
    )

    return coefficient * displacement
```

A few important implementation details\:

- `np.linalg.norm()` calculates the Euclidean distance.
- `sr6` calculates $(\sigma/r)^6$\.
- We reuse this term rather than calculating powers repeatedly.
- The force is a vector, whereas the potential energy is a scalar.

This distinction will become essential when we implement many-atom systems.

## 5.4 Test the physics, not just the code

Scientific software needs tests based on physical laws.

Create `tests/test_potentials.py`\:

```python
import numpy as np

from moleculab.potentials import (
    lennard_jones_energy,
    lennard_jones_force,
)

def test_energy_at_minimum():

    epsilon = 0.997
    sigma = 0.34

    r_min = 2 ** (1 / 6) * sigma

    energy = lennard_jones_energy(
        r_min, epsilon, sigma
    )

    assert np.isclose(energy, -epsilon)

def test_force_at_minimum():

    epsilon = 0.997
    sigma = 0.34

    r_min = 2 ** (1 / 6) * sigma

    force = lennard_jones_force(
        [r_min, 0, 0],
        epsilon,
        sigma
    )

    assert np.allclose(force, [0, 0, 0], atol=1e-10)

def test_newton_third_law():

    epsilon = 0.997
    sigma = 0.34

    displacement = np.array([0.4, 0.1, 0.0])

    force_i = lennard_jones_force(
        displacement, epsilon, sigma
    )

    force_j = lennard_jones_force(
        -displacement, epsilon, sigma
    )

    assert np.allclose(force_i, -force_j)
```

We have now tested three fundamental physical properties\:

1. The minimum energy is $-\epsilon$\.
2. The force is zero at equilibrium.
3. The forces on the two atoms are equal and opposite.

Run\:

```bash
pytest -v
```

## 5.5 A more powerful test: Compare analytical and numerical forces

This is a particularly important exercise.

We derived the force analytically\:

$$
F(r)=-\frac{dU}{dr}
$$

But how do we know our implementation is correct?

We can approximate the derivative numerically using finite differences\:

$$
\frac{dU}{dr}\approx
\frac{U(r+h)-U(r-h)}{2h}
$$

Therefore\:

$$
\boxed{
F_{\mathrm{numerical}}(r)=
-\frac{U(r+h)-U(r-h)}{2h}
}
$$

Implement\:

```python
def numerical_force(r, epsilon, sigma, h=1e-6):

    u_plus = lennard_jones_energy(
        r + h, epsilon, sigma
    )

    u_minus = lennard_jones_energy(
        r - h, epsilon, sigma
    )

    return -(u_plus - u_minus) / (2 * h)
```

Compare\:

```python
epsilon = 0.997
sigma = 0.34
r = 0.40

f_analytical = lennard_jones_force(
    [r, 0, 0], epsilon, sigma
)[0]

f_numerical = numerical_force(
    r, epsilon, sigma
)

print("Analytical:", f_analytical)
print("Numerical:", f_numerical)
print("Difference:", abs(f_analytical - f_numerical))
```

The two values should agree closely, with small differences caused by numerical approximation and floating-point precision.

This technique is called *finite-difference verification of analytical gradients*\. It will be one of our most useful tools throughout the project, especially when implementing complicated torsional and electrostatic forces.

## 5.6 Visualize our first molecular potential

Create `notebooks/01_lennard_jones.ipynb` or run the following Python script\:

```python
import numpy as np
import matplotlib.pyplot as plt

from moleculab.potentials import lennard_jones_energy

epsilon = 0.997
sigma = 0.34

distances = np.linspace(
    0.25, 0.9, 500
)

energies = [
    lennard_jones_energy(r, epsilon, sigma)
    for r in distances
]

r_min = 2 ** (1 / 6) * sigma

plt.figure(figsize=(8, 5))

plt.plot(distances, energies)

plt.axhline(0, linestyle="--")
plt.axvline(r_min, linestyle="--")

plt.xlabel("Distance (nm)")
plt.ylabel("Potential energy (kJ/mol)")

plt.title("Lennard-Jones potential")

plt.grid(True)
plt.show()
```

You should see the characteristic Lennard-Jones potential well.

The values above are illustrative parameters for a carbon-like interaction, not validated parameters for imatinib or ABL1.

---

# 6\. Understand what we have actually achieved

We have not yet simulated a molecule.

We have implemented a physical interaction between two idealized particles.

But this is a significant foundation.

Our first implementation already contains\:

| Concept | Implemented? |
|--------|--------|
| Potential energy | (Yes) |
| Analytical force | (Yes) |
| Vector mathematics | (Yes) |
| Newton's third law | (Tested) |
| Numerical differentiation | (Tested) |
| Multiple interacting atoms | (Next) |
| Time integration | (Later) |
| Molecular dynamics | (Later) |

The next logical step is to generalize this two-particle model to an arbitrary number of particles.

For $N$ atoms\:

$$
\boxed{
U_{\mathrm{total}}=
\sum_{i<j} U_{ij}(r_{ij})
}
$$

And the force on atom $i$ is\:

$$
\boxed{
\mathbf F_i=
\sum_{j\neq i}\mathbf F_{ij}
}
$$

This introduces an important computational problem: how to calculate all pairwise interactions efficiently.

For $N$ atoms, the simplest algorithm evaluates\:

$$
\frac{N(N-1)}{2}
$$

unique atom pairs.

That is an $O(N^2)$ algorithm. Later we will implement cutoffs, neighbor lists and eventually spatial cell lists to improve performance.

## 7\. Make the project installable

Add a minimal `pyproject.toml`\:

```toml
[build-system]
requires = ["setuptools>=68"]
build-backend = "setuptools.build_meta"

[project]
name = "moleculab"
version = "0.1.0"
description = "An educational molecular physics engine"
requires-python = ">=3.10"
dependencies = ["numpy", "matplotlib"]

[tool.setuptools.packages.find]
where = ["src"]
```

Add `src/moleculab/__init__.py`\:

```python
from .potentials import (
    lennard_jones_energy,
    lennard_jones_force,
)
```

Install your library\:

```bash
pip install -e .
```

Now you can import your own library from any Python script in the environment\:

```python
import moleculab

print(moleculab.lennard_jones_energy(
    0.40, 0.997, 0.34
))
```

---

## 8\. Your first physics exercises

Before moving forward, I would like you to understand these four questions.

**Milestone 1 knowledge check** (/4)

[Check my answers]

## 9\. What we will implement next

(Next milestone)

**Milestone 2 - An N-particle molecular system**

We will introduce\:

1. An `Atom` class with mass, position, velocity and charge.
2. A `System` class containing multiple atoms.
3. Pairwise distance calculations.
4. A general potential-energy calculator.
5. A general force calculator.
6. Conservation of total momentum.
7. A simulation of 10-100 interacting particles.

Then we will introduce velocity Verlet\:

$$
\mathbf r(t+\Delta t)=
\mathbf r(t)+
\mathbf v(t)\Delta t+
\frac12\mathbf a(t)\Delta t^2
$$

$$
\mathbf v(t+\Delta t)=
\mathbf v(t)+
\frac12[\mathbf a(t)+\mathbf a(t+\Delta t)]
\Delta t
$$

This will be our first actual molecular dynamics engine.

**The longer-term destination\:** We will progressively replace our simple particle system with chemically connected atoms, amino acids, water molecules and eventually the imatinib-ABL1 complex. We will understand each additional physical term before implementing it.

One important rule for this project: we will never add an equation to the library without understanding what it represents, its assumptions, its units and how to test it.

Continue Moleculab lessons tomorrow

## Prompt:
11/9/2026, 12:22:50 PM

the&#x20;

## Lennard-Jones potential

example is already too complicated we want a much milder and gradual introduction to physics to gain a solid understanding from first principles

## Response:
11/9/2026, 12:22:52 PM · gpt-6

You're right. We started building a molecular simulation engine before establishing the physics needed to understand it. That's backwards.

Let's reset **Moleculab**\. We will slow down considerably and build your understanding from first principles. No Lennard-Jones equations, no molecular dynamics code, and no complicated terminology until the underlying ideas make sense.

The goal is not just to learn how to run computational drug discovery software. It's to understand what the software is calculating, why those calculations make physical sense, and where their limitations come from.

# A new starting point: physics from zero

We'll follow a simple rule\:

**One new physical idea at a time. Understand it, experiment with it, and only then move on.**

Here's the revised learning path\:

1. **What is the physical world made of?**
   
   Matter, atoms, molecules, mass, size, and the difference between a physical object and a mathematical model.
2. **Position and distance**
   
   How do we describe where things are? What is a coordinate system? How do we measure distances?
3. **Motion**
   
   Position, time, speed, velocity, and acceleration.
4. **Forces**
   
   What makes things move? Pushes, pulls, gravity, and Newton's laws.
5. **Energy**
   
   What is energy? Kinetic energy, potential energy, and conservation of energy.
6. **Interactions between atoms**
   
   Why atoms attract or repel one another. Only now do we introduce simple models of atomic interactions.
7. **From atoms to molecules**
   
   Chemical bonds, molecular geometry, and molecular flexibility.
8. **From physics to computational chemistry**
   
   Force fields, energy minimization, molecular dynamics, quantum mechanics, and free-energy calculations.

# Lesson 1 - What are we actually studying?

Let's begin with something familiar.

Water is made of molecules. Each water molecule contains\:

- Two hydrogen atoms.
- One oxygen atom.

An atom is a tiny unit of matter. Different kinds of atoms have different physical properties.

For now, let's focus on three ideas\:

**1\. Mass**

Mass describes how much matter an object contains. An oxygen atom has more mass than a hydrogen atom.

**2\. Position**

An atom occupies a location in space. We can describe that location using coordinates.

**3\. Interaction**

Atoms can influence one another. They may attract or repel one another depending on their physical and chemical properties.

These three ideas will eventually be the foundation of molecular simulation.

## A simple thought experiment

Imagine two balls floating in empty space.

We want to understand\:

- Where are the balls?
- How far apart are they?
- Are they moving?
- Do they influence one another?

These are exactly the kinds of questions a molecular simulation asks about atoms.

## Let's check our understanding

No mathematics required yet.

**Question 1**

Imagine an oxygen atom and a hydrogen atom. Which statement about their mass is correct?

They have exactly the same mass.

The oxygen atom has more mass.

The hydrogen atom has more mass.

[Check answer]

## How we will use programming

We won't write a physics engine yet.

When the concepts are clear, we'll use very small Python experiments to make them tangible. For example, we'll start by representing a single object and describing its position\:

```python
position = 5  # position along a line
```

Then we'll explore what happens when that object moves. Eventually, we'll extend these simple ideas to atoms, molecules, proteins, and drug binding.

**The destination remains computational drug discovery, but the journey starts with basic physics.**

One important question before we continue: **Would you prefer that I assume no prior knowledge of physics or mathematics at all, and explain even basic concepts such as units, equations, and graphs from scratch?**

## Prompt:
11/9/2026, 12:28:30 PM

assume solid math and programming prior knowledge, and basic physics and chemistry knowledge

## Response:
11/9/2026, 12:28:34 PM · gpt-6-mini

Perfect. That gives us a much better starting point.

We'll assume you are comfortable with\:

- **Mathematics\:** algebra, calculus, derivatives, integrals, vectors, matrices, and basic differential equations.
- **Programming\:** Python, functions, classes, NumPy, testing, and numerical methods.
- **Physics\:** basic mechanics, forces, energy, and thermodynamics.
- **Chemistry\:** atoms, molecules, chemical bonds, charges, polarity, and molecular structure.

We won't spend time teaching those subjects from scratch. Instead, we'll develop the deeper physical intuition behind them and connect it to computational chemistry.

# Our revised teaching philosophy

I suggest three layers for every topic\:

1. **Physical intuition\:** What is happening in nature, and why?
2. **Mathematical formulation\:** How do we describe it quantitatively? What assumptions are involved?
3. **Computational implementation\:** How do we translate the mathematics into code, and how do we verify the result?

We will resist introducing a mathematical model just because it is commonly used in molecular simulation. First, we need to understand what physical phenomenon it represents.

# Revised curriculum

| Stage | Topic | Central question |
|--------|--------|--------|
| 1 | Statistical mechanics | How do microscopic particles produce macroscopic properties? |
| 2 | Classical mechanics at molecular scales | What determines atomic motion? |
| 3 | Potential energy surfaces | What does molecular energy actually mean? |
| 4 | Interatomic interactions | Why do atoms attract and repel? |
| 5 | Molecular mechanics | How can we model molecular structure and flexibility? |
| 6 | Numerical optimization | How do we find stable molecular configurations? |
| 7 | Molecular dynamics | How do we simulate molecular motion? |
| 8 | Statistical mechanics of binding | How do we connect simulations to binding affinity? |
| 9 | Free-energy calculations | How can we calculate relative binding affinities? |
| 10 | Quantum mechanics | Where do classical models break down? |
| 11 | DFT and QM/MM | How do we model electronic structure and chemical reactions? |
| 12 | Computational drug discovery | How do we combine these methods to study real drug-target systems? |

# Lesson 1: The central idea - microscopic versus macroscopic physics

Let's start with a fundamental question\:

**Why do we need statistical mechanics to understand molecular systems?**

Consider a glass of water.

At the macroscopic level, we can measure\:

- Temperature
- Pressure
- Density
- Heat capacity

But water consists of an enormous number of interacting molecules, each with its own position, velocity, and energy.

The macroscopic properties emerge from the collective behavior of those microscopic particles.

This is the key conceptual bridge\:

$$
\text{Microscopic physics}
\quad\longrightarrow\quad
\text{Statistical mechanics}
\quad\longrightarrow\quad
\text{Macroscopic properties}
$$

For computational drug discovery, this matters enormously.

A protein and a ligand do not have one single configuration. They continuously fluctuate between many configurations. Binding affinity is therefore not simply the energy of one protein-ligand structure.

It is fundamentally connected to the statistical behavior of molecular configurations.

## A first important distinction

Consider two quantities\:

**Potential energy\:** The energy associated with a particular molecular configuration.

**Free energy\:** A thermodynamic quantity that accounts for the statistical contributions of accessible configurations, including both energetic and entropic effects.

This distinction will become central when we study docking, molecular dynamics, MM/GBSA, and free-energy perturbation.

## A small thought experiment

Imagine two possible molecular configurations\:

- Configuration A has lower potential energy but is accessible in relatively few ways.
- Configuration B has slightly higher potential energy but can occur in many more configurations.

Which one is thermodynamically favored?

The answer is not necessarily A. The number of accessible configurations matters too. That's the role of entropy.

---

## What we'll build toward

Eventually, we will derive the canonical partition function\:

$$
Z = \sum_i e^{-\beta E_i}
$$

where

$$
\beta = \frac{1}{k_B T}
$$

This equation is one of the foundations of molecular simulation. It connects molecular energy, temperature, and statistical populations.

But we won't rush into it. First, we'll understand why it has this form and what it means physically.

**Next lesson\:** Why molecular systems need statistical mechanics, and how probability distributions emerge from microscopic energy and temperature. We'll derive the Boltzmann distribution from physical reasoning before writing any simulation code.

## Prompt:
11/9/2026, 2:30:51 PM

elaborate in detail on all the terms and their connections from "A protein and a ligand..." onwards, and up to "any simulation code".

## Response:
11/9/2026, 2:30:55 PM · gpt-6-mini

Let's slow down and unpack this carefully. This is one of the most important conceptual foundations in computational drug discovery because it explains why we need statistical mechanics, why a molecular structure is not enough to predict drug binding, and why methods such as molecular dynamics (MD) and free-energy perturbation (FEP) exist.

We'll build the explanation in a sequence, connecting each concept to the next.

# 1\. A protein and a ligand are dynamic molecular systems

First, let's define what we mean by a protein and a ligand.

## 1.1 What is a protein?

A protein is a biological macromolecule composed of a chain of amino acids.

There are 20 common proteinogenic amino acids. Each has a shared backbone structure and a chemically distinct side chain.

A protein's amino-acid sequence determines its chemical composition and strongly influences how it folds into a three-dimensional structure.

For example, consider a protein kinase such as BCR-ABL1, the molecular target of imatinib.

A kinase is an enzyme that transfers phosphate groups, typically from ATP, to other molecules. BCR-ABL1 is an abnormal tyrosine kinase associated with chronic myeloid leukemia (CML).

The protein's three-dimensional structure contains regions with different biological functions. One of these is the ATP-binding pocket.

## 1.2 What is a ligand?

In molecular biology and drug discovery, a ligand is a molecule that binds to another molecular entity, such as a protein.

A ligand might be\:

- A naturally occurring metabolite.
- A substrate of an enzyme.
- A signaling molecule.
- A drug molecule.

Imatinib, for example, is a small-molecule ligand of the ABL1 kinase domain.

The ligand has its own chemical structure, including\:

- Atoms.
- Covalent bonds.
- Bond angles.
- Torsional degrees of freedom.
- Partial charges.
- Electronic structure.

The protein and ligand interact through a combination of physical and chemical interactions, including electrostatics, dispersion, hydrogen bonding, and steric repulsion.

## 1.3 Why isn't the protein-ligand structure fixed?

A crystal structure or cryo-EM structure provides a structural model derived from experimental measurements.

It does not mean that the protein and ligand are completely stationary in reality.

At nonzero temperature, atoms have thermal motion. Bonds vibrate, angles fluctuate, side chains rotate, and larger parts of the protein may undergo conformational changes.

The ligand also moves. Its flexible groups rotate around certain bonds, and it can adopt different conformations.

Consequently, a protein-ligand complex is not one rigid arrangement of atoms. It is a dynamic system that explores a collection of molecular configurations.

This is our first major connection\:

**A molecular structure is a snapshot; molecular behavior is a process.**

---

# 2\. What exactly is a molecular configuration?

We need a precise mathematical description of a molecular configuration before we can discuss its energy.

Suppose we have a system containing $N$ atoms.

Each atom has a position in three-dimensional space\:

$$
\mathbf r_i =
\begin{pmatrix}
x_i\\
y_i\\
z_i
\end{pmatrix}
$$

where $i$ identifies the atom.

The complete set of atomic positions is\:

$$
\mathbf R =
(\mathbf r_1,\mathbf r_2,\ldots,\mathbf r_N)
$$

We call $\mathbf R$ the *configuration* of the molecular system.

For example, if a system has 1000 atoms, its configuration contains 3000 Cartesian coordinates\:

$$
\mathbf R \in \mathbb R^{3000}
$$

This is a point in a 3000-dimensional configuration space.

## 2.1 Configuration space versus physical space

This distinction is important.

Physical space has three dimensions: x, y, and z.

But the configuration space of a molecular system has three dimensions for every atom.

A single point in configuration space represents one complete arrangement of all atoms.

A different point represents a different arrangement.

For a protein-ligand complex, changing one side-chain torsion angle changes the point in configuration space, even if the rest of the protein remains nearly unchanged.

## 2.2 What about molecular motion?

A molecular configuration describes positions, but not velocities.

To describe the complete classical mechanical state of a molecular system, we also need the momenta\:

$$
\mathbf p_i = m_i\mathbf v_i
$$

where\:

- $m_i$ is the mass of atom $i$\.
- $\mathbf v_i$ is its velocity.
- $\mathbf p_i$ is its momentum.

The complete state is therefore\:

$$
(\mathbf R,\mathbf P)
$$

where $\mathbf P$ represents all atomic momenta.

This is a point in *phase space*\.

For $N$ atoms, classical phase space has $6N$ dimensions: $3N$ position coordinates and $3N$ momentum coordinates.

The distinction\:

| Concept | Contains | Meaning |
|--------|--------|--------|
| Configuration | Positions | Where are the atoms? |
| Phase-space state | Positions and momenta | Where are the atoms and how are they moving? |
| Trajectory | A sequence of states over time | How does the system evolve? |

This distinction will become useful when we discuss molecular dynamics.

---

# 3\. Why does a molecular configuration have energy?

Now we can introduce potential energy.

Potential energy is a scalar quantity associated with the arrangement of a system.

For a molecular system, we can write\:

$$
U = U(\mathbf R)
$$

This means that the potential energy depends on the atomic positions.

Change the atomic positions, and the potential energy generally changes.

For example\:

- Stretching a covalent bond changes its energy.
- Bending a bond angle changes its energy.
- Rotating a flexible molecular group can change its energy.
- Bringing two atoms too close can produce a strong repulsive interaction.
- Moving charged groups relative to one another changes electrostatic energy.

The function $U(\mathbf R)$ is called the *potential energy function*\.

In molecular mechanics, it is commonly called the *force field energy function* when it is described using a classical force field.

## 3.1 A simple physical analogy

Imagine a ball on a hilly landscape.

The height of the landscape represents potential energy.

A ball at the bottom of a valley has lower potential energy than a ball at the top of a hill.

For a molecular system, however, the landscape is not a two-dimensional geographical surface.

It exists in a high-dimensional configuration space.

We call it the *potential energy surface*, or PES.

$$
U = U(\mathbf R)
$$

Each point on the surface corresponds to one molecular configuration and has an associated potential energy.

## 3.2 What does a minimum on this surface mean?

A local minimum is a configuration where small changes in atomic positions do not decrease the potential energy.

Mathematically, at a stationary point\:

$$
\nabla U(\mathbf R)=0
$$

For a local minimum, the curvature must also satisfy the appropriate stability conditions, usually expressed through the Hessian matrix.

The Hessian is\:

$$
H_{ij} =
\frac{\partial^2 U}{\partial q_i\partial q_j}
$$

where $q_i$ and $q_j$ are coordinates in configuration space.

A minimum on the potential energy surface represents a locally stable arrangement under the chosen energy model.

But here's a crucial point\:

**The lowest-energy configuration is not automatically the most populated configuration at a given temperature.**

To understand why, we need to introduce statistical mechanics.

---

# 4\. Why does temperature make molecular systems dynamic?

At absolute zero, a classical system at its global potential energy minimum would have no thermal kinetic energy in the idealized picture.

At ordinary temperatures, molecules have thermal energy and continually move.

The total classical energy is typically written\:

$$
E(\mathbf R,\mathbf P)
=
U(\mathbf R)+K(\mathbf P)
$$

Here\:

- $U(\mathbf R)$ is potential energy.
- $K(\mathbf P)$ is kinetic energy.

For a system of classical particles\:

$$
K(\mathbf P)
=
\sum_{i=1}^{N}
\frac{\mathbf p_i^2}{2m_i}
$$

This equation tells us that the total energy depends on both atomic positions and momenta.

The atoms move because forces act on them, and their momenta change over time.

At ordinary biological temperatures, proteins and ligands constantly fluctuate.

For example\:

- A side chain may rotate between different orientations.
- A loop may move between open and closed conformations.
- A ligand may change its torsion angles.
- Water molecules may enter and leave the binding pocket.

These fluctuations are not defects or numerical noise. They are fundamental properties of molecular systems.

## 4.1 What does temperature mean at the molecular level?

Temperature is a macroscopic thermodynamic property, but it is connected to microscopic motion.

For a classical system in thermal equilibrium, the equipartition theorem tells us that each independent quadratic degree of freedom contributes, on average\:

$$
\frac{1}{2}k_BT
$$

to the energy.

Here\:

- $k_B$ is the Boltzmann constant.
- $T$ is absolute temperature in kelvin.

For a simple three-dimensional particle, the three translational kinetic degrees of freedom contribute\:

$$
\langle K\rangle =
\frac{3}{2}k_BT
$$

This is an average over an ensemble or sufficiently long equilibrium observation, not a statement that every particle always has exactly this kinetic energy.

At room temperature, $k_BT$ is approximately\:

$$
k_BT \approx 2.48\ \text{kJ/mol}
$$

This is a useful energy scale in molecular chemistry.

It helps us judge whether an energy difference is small or large relative to thermal fluctuations.

---

# 5\. Why one molecular structure cannot tell us everything

Let's return to our protein-ligand complex.

Imagine that imatinib is bound in the ABL1 binding pocket.

There may be many accessible configurations of the complex\:

- The ligand can rotate around flexible bonds.
- Protein side chains can change orientation.
- Hydrogen bonds can fluctuate.
- Water molecules can rearrange.
- The protein backbone can undergo small movements.

Each configuration has a different potential energy.

Let's imagine three simplified configurations\:

**Configuration A**

Low potential energy, relatively uncommon arrangement.

---

**Configuration B**

Slightly higher potential energy, but many similar arrangements are accessible.

---

**Configuration C**

Higher potential energy, relatively unfavorable arrangement.

Which configuration will we observe most frequently?

We cannot answer that by looking only at the potential energy of one configuration.

We must consider two things\:

1. How energetically favorable is the configuration?
2. How many configurations of similar energy are accessible?

The second point introduces entropy.

---

# 6\. Entropy: why the number of accessible configurations matters

Entropy is often introduced as a measure of disorder. That intuition can be useful, but for molecular statistical mechanics, a more precise starting point is that entropy is related to the number and probability distribution of accessible microscopic states.

For a system with discrete states, the Gibbs entropy is\:

$$
S=-k_B\sum_i p_i\ln p_i
$$

where\:

- $S$ is entropy.
- $p_i$ is the probability of microscopic state $i$\.
- $k_B$ is the Boltzmann constant.

If there are $W$ equally probable states, then\:

$$
p_i=\frac{1}{W}
$$

and therefore\:

$$
S=k_B\ln W
$$

This is the familiar Boltzmann entropy formula.

The important intuition is that a system with many accessible microscopic states can have higher entropy than one restricted to only a few states.

## 6.1 A simple example

Suppose a ligand can adopt two broad classes of configurations.

- Class A: 1 accessible configuration, with energy 0.
- Class B: 100 accessible configurations, each with energy 2 kJ/mol.

At room temperature\:

$$
k_BT\approx2.48\ \text{kJ/mol}
$$

The energetic penalty of class B is only about $0.81k_BT$\.

Although each individual B configuration has higher energy, there are many more ways for the system to occupy class B.

This creates competition between energy and entropy.

We will quantify this competition next.

---

# 7\. From energy and entropy to free energy

The Gibbs free energy is\:

$$
G=H-TS
$$

where\:

- $G$ is Gibbs free energy.
- $H$ is enthalpy.
- $T$ is temperature.
- $S$ is entropy.

At constant temperature and pressure, Gibbs free energy is the relevant thermodynamic potential for determining spontaneous changes and equilibrium.

For many molecular simulation discussions, we instead first introduce the Helmholtz free energy\:

$$
A=U-TS
$$

where $U$ is internal energy.

Why introduce two kinds of free energy?

Because they apply naturally to different thermodynamic constraints\:

| Free energy | Natural conditions |
|--------|--------|
| Helmholtz $A$ | Constant temperature and volume |
| Gibbs $G$ | Constant temperature and pressure |

Molecular simulations commonly use constant-volume or constant-pressure ensembles, so both concepts are relevant.

## 7.1 What does free energy mean physically?

Free energy combines energetic and entropic contributions.

A state can be favorable because\:

- It has low energy.
- It has high entropy.
- Or a balance of both.

Return to our ligand example.

Binding may create favorable interactions between the ligand and protein. This can lower the interaction energy.

But binding may also restrict the ligand's movement and reduce its conformational entropy.

The protein may lose some conformational freedom too.

Meanwhile, water molecules previously occupying the binding pocket may be released into bulk solvent, potentially increasing their entropy.

Therefore, binding is a competition between multiple energetic and entropic contributions.

This is why a ligand with many strong-looking interactions does not necessarily have the strongest binding affinity.

**Free energy, not the interaction energy of one structure, governs equilibrium binding.**

---

# 8\. Probability: how often do molecular configurations occur?

We now have all the ingredients needed to connect energy, temperature, and probability.

For a system in thermal equilibrium with a heat bath at temperature $T$, the probability of observing a microscopic state $i$ with energy $E_i$ is\:

$$
p_i =
\frac{e^{-\beta E_i}}{Z}
$$

where\:

$$
\beta=\frac{1}{k_BT}
$$

and\:

$$
Z=\sum_i e^{-\beta E_i}
$$

This is the **Boltzmann distribution**\.

Let's unpack every term.

| Term | Meaning | Physical role |
|--------|--------|--------|
| $p_i$ | Probability of state $i$ | How often we expect to observe that state |
| $E_i$ | Energy of state $i$ | Energetic favorability |
| $k_B$ | Boltzmann constant | Connects microscopic energy to temperature |
| $T$ | Absolute temperature | Sets the scale of thermal fluctuations |
| $\beta$ | $1/(k_BT)$ | Converts energy into a dimensionless quantity |
| $Z$ | Partition function | Normalizes probabilities and summarizes the ensemble |

The exponential factor is important\:

$$
e^{-\beta E_i}
$$

Lower-energy states receive greater statistical weight.

However, the probability of a broad collection of configurations depends on the total weight of all those configurations, not just the weight of one configuration.

## 8.1 Why do we need the partition function?

Without $Z$, the Boltzmann weights are not normalized probabilities.

The partition function ensures\:

$$
\sum_i p_i=1
$$

It also contains information about the thermodynamics of the system.

For example, in the canonical ensemble\:

$$
A=-k_BT\ln Z
$$

This equation connects the microscopic energy landscape to the Helmholtz free energy.

The partition function is therefore much more than a normalization constant. It is a central mathematical object in statistical mechanics.

---

# 9\. From individual states to molecular ensembles

So far, we've discussed discrete states for simplicity.

Real molecular systems have continuous coordinates and momenta.

Instead of summing over a finite set of states, we integrate over phase space.

For a classical system, the canonical probability density is proportional to\:

$$
\rho(\mathbf R,\mathbf P)
\propto
e^{-\beta H(\mathbf R,\mathbf P)}
$$

where the Hamiltonian is\:

$$
H(\mathbf R,\mathbf P)
=
K(\mathbf P)+U(\mathbf R)
$$

The partition function becomes a phase-space integral, with appropriate normalization factors\:

$$
Z =
\frac{1}{h^{3N}N!}
\int
e^{-\beta H(\mathbf R,\mathbf P)}
\,d\mathbf R\,d\mathbf P
$$

For classical indistinguishable particles, this expression includes Planck's constant $h$ and the factor $N!$\. More complete treatments account for particle types, constraints, and other details.

For many molecular simulation calculations, we can integrate out the momenta analytically. The remaining configurational partition function is proportional to\:

$$
Z_\mathrm{conf}
=
\int e^{-\beta U(\mathbf R)}\,d\mathbf R
$$

This equation is particularly important.

It tells us that the equilibrium distribution of molecular configurations is determined by the potential energy function and temperature.

It is the conceptual foundation for classical equilibrium molecular simulations.

---

# 10\. How does all this connect to computational drug discovery?

Let's connect the pieces directly.

Suppose we want to understand how strongly imatinib binds to ABL1.

We need to compare the thermodynamic properties of different states\:

- ABL1 with imatinib bound.
- ABL1 without imatinib.
- Imatinib free in solution.

The binding reaction can be represented as\:

$$
P+L\rightleftharpoons PL
$$

where\:

- $P$ is the protein.
- $L$ is the ligand.
- $PL$ is the protein-ligand complex.

The standard binding free energy is related to the equilibrium association constant\:

$$
\Delta G^\circ_\mathrm{bind}
=
-RT\ln(K_a C^\circ)
$$

where\:

- $R$ is the gas constant.
- $K_a$ is the association constant.
- $C^\circ$ is the standard concentration, conventionally 1 mol/L, making the logarithm's argument dimensionless.

Equivalently, using the dissociation constant $K_d=1/K_a$\:

$$
\Delta G^\circ_\mathrm{bind}
=
RT\ln\left(\frac{K_d}{C^\circ}\right)
$$

These relationships connect molecular binding to measurable experimental quantities.

A lower $K_d$ corresponds to tighter binding and a more favorable standard binding free energy.

## 10.1 Why a molecular dynamics simulation is useful

A molecular dynamics simulation attempts to generate a time-dependent trajectory of molecular configurations.

It allows us to examine\:

- Protein flexibility.
- Ligand flexibility.
- Hydrogen-bond persistence.
- Solvent rearrangement.
- Conformational changes.
- The range of configurations sampled by the complex.

However, an MD trajectory is not automatically a free-energy calculation.

It provides samples of molecular behavior under a specified model and ensemble. Calculating free energies requires additional statistical-mechanical analysis.

## 10.2 Why free-energy perturbation is useful

Suppose we want to compare two related drug candidates, ligand A and ligand B.

Rather than calculating their absolute binding free energies independently, we can construct a thermodynamic cycle that transforms A into B in solution and in the protein binding site.

The difference between these transformations gives the relative binding free energy\:

$$
\Delta\Delta G_\mathrm{bind}
=
\Delta G_{A\rightarrow B}^{\mathrm{bound}}
-
\Delta G_{A\rightarrow B}^{\mathrm{solution}}
$$

This is the central idea behind relative free-energy perturbation (FEP).

It is rooted directly in the statistical mechanics we have just developed.

---

# 11\. The complete conceptual connection

Here is the logical chain we have built\:

**Protein + ligand**

A collection of interacting atoms

**Molecular configuration**

The positions of all atoms

**Potential energy surface**

Energy associated with each configuration

**Thermal fluctuations**

Molecules explore different configurations at finite temperature

**Statistical mechanics**

Connects energies, temperature, and probabilities

**Free energy**

Combines energetic and entropic contributions

**Binding affinity**

A thermodynamic description of protein-ligand binding

# 12\. What we have deliberately not covered yet

There are several important questions that naturally follow from this discussion\:

1. **Where does the potential energy function come from?** How do we calculate the energy of a protein with thousands of atoms?
2. **How do we calculate forces?** How does the potential energy surface determine molecular motion?
3. **How do we generate molecular trajectories?** How do Newton's equations become a numerical simulation?
4. **How do we know that a simulation samples equilibrium correctly?** What is the difference between a trajectory, an ensemble, and an equilibrium distribution?
5. **How do we calculate free energy from simulation data?** What are sampling, convergence, and statistical uncertainty?

These questions lead to molecular mechanics, force fields, MD, and free-energy methods.

But we should not jump there yet.

## Check your understanding

I'd like to check three conceptual distinctions before introducing the next layer.

[Check all answers]

My recommendation is that our next lesson should focus on **the difference between a microscopic state, a molecular configuration, a trajectory, and a statistical ensemble**\. We'll develop that distinction with concrete examples and simple mathematics, because it is essential to understanding everything from molecular dynamics to FEP.

## Prompt:
11/9/2026, 2:41:43 PM

explain 3 and 4 in much greater detail so that it makes intuitive sense

## Response:
11/9/2026, 2:41:45 PM · gpt-6-mini

Let's go much deeper into points **3\. Potential energy surfaces** and **4\. Thermal fluctuations**\.

These two concepts are closely connected. Understanding their relationship is essential because almost every computational method in molecular drug discovery relies on it.

The central question we want to answer is\:

**How does a molecule's three-dimensional structure determine its energy, and why does that molecule continually change its structure at a given temperature?**

We'll approach this physically first, then mathematically, and finally connect everything to protein-ligand binding.

# 3\. Potential energy surfaces: how molecular structure determines energy

## 3.1 Start with a simple physical example

Imagine a ball moving over a landscape.

There are hills and valleys.

- When the ball is at the bottom of a valley, it has relatively low gravitational potential energy.
- When it is at the top of a hill, it has relatively high gravitational potential energy.
- If the ball moves, its potential energy changes depending on where it is.

The landscape defines a relationship between position and potential energy.

For a ball at height $h$, we can write\:

$$
U(h)=mgh
$$

where\:

- $m$ is mass.
- $g$ is gravitational acceleration.
- $h$ is height.

The important idea is that **the potential energy is a function of the position of the object**\.

Molecular systems work according to the same general principle, although their energy functions are much more complicated.

Instead of a ball's height, we consider the positions of atoms.

## 3.2 From one particle to a molecule

Consider a very simple molecule consisting of two atoms connected by a chemical bond.

The atoms cannot be separated by an arbitrary distance without changing the molecule's energy.

If the atoms are too close, their electron clouds overlap unfavorably, creating strong repulsion.

If they are too far apart, the bond becomes stretched and eventually breaks.

There is usually a preferred equilibrium bond length.

Let's call this distance $r_0$\.

We can describe the energy as a function of the interatomic distance\:

$$
U=U(r)
$$

A simplified qualitative picture looks like this\:

The curve has three important regions.

**Region A: At very short distances**

The atoms are pushed strongly apart. The energy increases dramatically.

**Region B: Around the equilibrium distance**

The energy reaches a minimum. This is the energetically preferred bond length.

**Region C: At large distances**

The atoms become increasingly separated. The bond becomes stretched, and the energy approaches the dissociation limit.

The minimum is important because it represents a stable equilibrium arrangement for the isolated diatomic molecule under this model.

Notice something interesting\:

The molecule is not simply trying to minimize its energy by bringing its atoms as close together as possible.

At short distances, repulsion dominates. At longer distances, attraction and bond formation matter.

The energy function captures this competition.

## 3.3 Why is the energy a function of atomic coordinates?

Now let's generalize from two atoms to a molecule containing many atoms.

Suppose a molecule has $N$ atoms.

Each atom has three Cartesian coordinates\:

$$
\mathbf r_i=(x_i,y_i,z_i)
$$

We can collect all coordinates into one vector\:

$$
\mathbf R =
(x_1,y_1,z_1,x_2,y_2,z_2,\ldots,x_N,y_N,z_N)
$$

The molecular potential energy is then\:

$$
\boxed{U=U(\mathbf R)}
$$

This is a compact way of saying\:

> The energy of a molecule depends on the positions of all its atoms.

For a molecule with 100 atoms, there are nominally 300 Cartesian coordinates.

Changing even one coordinate can change the energy.

For example\:

- Move one hydrogen atom: a bond angle might change.
- Rotate one carbon-carbon bond: the molecular conformation changes.
- Move two charged groups closer: their electrostatic interaction changes.
- Move a ligand deeper into a protein pocket: multiple interactions change simultaneously.

The energy function takes all these changes into account.

## 3.4 What exactly is a potential energy surface?

The term *surface* can be misleading.

For a diatomic molecule, we only need one coordinate, the bond distance $r$\. Its energy function is a one-dimensional curve.

For a molecule with three independent coordinates, the energy can be visualized as a three-dimensional landscape.

For a real protein-ligand complex with thousands of atoms, the energy function exists in a space with thousands of dimensions.

This is the potential energy surface, or PES.

$$
U:\mathbb R^{3N}\rightarrow\mathbb R
$$

It maps every possible molecular configuration to one energy value.

Imagine the following\:

- One point on the PES represents one complete molecular arrangement.
- Another point represents the same molecule with a slightly rotated side chain.
- A third point represents a different ligand conformation.
- Each point has an associated energy.

The PES is therefore a mathematical representation of the energetic consequences of changing molecular geometry.

## 3.5 What do valleys, hills, and passes mean for molecules?

Let's use the landscape analogy again.

In molecular systems\:

| Landscape feature | Molecular interpretation |
|--------|--------|
| Valley | A relatively stable molecular configuration or conformational basin |
| Bottom of valley | A local energy minimum |
| Hill | Energetically unfavorable configurations |
| Pass between valleys | A transition region between conformational basins |
| Height of a pass | Energy barrier associated with a transition |

For example, imagine a protein side chain that can point in two different directions.

Both orientations might be relatively stable.

Between them, the side chain may need to rotate through a less favorable orientation.

The two stable orientations correspond to different regions of the potential energy surface. The intermediate orientation may have higher energy.

This is why molecular conformations can be separated by energy barriers.

**Important distinction\:** A minimum is not necessarily the global minimum.

A molecule may have many local minima, corresponding to different conformations.

A local minimum is stable against sufficiently small perturbations, but another conformation may have lower energy.

## 3.6 How does energy produce force?

This is the connection that makes the potential energy surface physically useful.

In classical mechanics, force is related to the negative spatial derivative of potential energy\:

$$
\boxed{\mathbf F=-\nabla U}
$$

For one dimension\:

$$
F(x)=-\frac{dU}{dx}
$$

Why the negative sign?

Because force points in the direction in which potential energy decreases most rapidly.

Imagine a ball on a slope.

- If the ball is on a downward slope, the force points downhill.
- If the ball is on an upward slope, the force points in the opposite direction.
- At the bottom of a valley, the slope is zero, so the net force is zero.

The same principle applies to atoms.

For atom $i$\:

$$
\mathbf F_i=-\frac{\partial U}{\partial\mathbf r_i}
$$

This means that if we know the potential energy function, we can calculate the force acting on each atom.

That force can then be used to calculate atomic acceleration through Newton's second law\:

$$
\mathbf F_i=m_i\mathbf a_i
$$

Combining the two equations\:

$$
\boxed{
m_i\mathbf a_i=
-\frac{\partial U}{\partial\mathbf r_i}
}
$$

This equation is one of the foundations of classical molecular dynamics.

We will return to it later.

## 3.7 What does a minimum really mean?

At a local energy minimum\:

$$
\nabla U(\mathbf R_0)=0
$$

This means that the forces on the atoms are zero at that configuration, subject to any constraints and boundary conditions.

But this alone is not enough to establish that we have a minimum.

A stationary point could also be a maximum or a saddle point.

We examine the second derivatives using the Hessian matrix\:

$$
H_{ij}=
\frac{\partial^2 U}{\partial q_i\partial q_j}
$$

At a strict local minimum, the Hessian is positive definite in the unconstrained directions.

Physically, this means that small displacements increase the potential energy.

Think of a ball at the bottom of a bowl. Move it slightly in any direction, and it moves to a higher position on the bowl.

There is also a subtle distinction worth remembering\:

- A minimum on a potential energy surface is a statement about potential energy.
- A thermodynamically favored state is a statement about free energy.

They are not necessarily the same thing.

We'll connect them in section 4.

---

# 4\. Thermal fluctuations: why molecules never stay still

Now we have an energy landscape.

But a protein does not simply settle into the lowest-energy configuration and remain there.

Why?

Because biological molecules exist at finite temperatures.

Let's understand what that means.

## 4.1 What does temperature do to a molecule?

Consider a simple molecule with a chemical bond.

At zero temperature, in an idealized classical picture, the system would occupy its lowest-energy configuration with zero kinetic energy.

At finite temperature, the atoms have kinetic energy.

They move, collide, vibrate, and change their positions.

A chemical bond therefore does not have exactly one fixed length.

Instead, its length fluctuates around an equilibrium value.

For example, suppose the equilibrium bond length is\:

$$
r_0=0.15\ \mathrm{nm}
$$

The instantaneous bond length might fluctuate around this value\:

$$
0.149,\quad 0.151,\quad 0.150,\quad 0.152\ \mathrm{nm}
$$

These are illustrative values, not a particular measured bond trajectory.

The key idea is that the bond length is a dynamic variable.

The same is true for\:

- Bond angles.
- Torsion angles.
- Side-chain orientations.
- Protein loops.
- Ligand conformations.

A protein is therefore better understood as a dynamic molecular system than as a static structure.

## 4.2 What causes the fluctuations?

At finite temperature, atoms possess kinetic energy.

For an atom\:

$$
K=\frac{1}{2}mv^2
$$

In three dimensions\:

$$
K=
\frac{1}{2}m(v_x^2+v_y^2+v_z^2)
$$

At thermal equilibrium, the kinetic energy is distributed statistically.

For a classical system, the equipartition theorem gives\:

$$
\langle K\rangle=\frac{3}{2}k_BT
$$

for three independent translational degrees of freedom.

Here\:

- $k_B$ is the Boltzmann constant.
- $T$ is the absolute temperature.
- The angle brackets indicate an ensemble or equilibrium average.

At 300 K\:

$$
k_BT\approx4.14\times10^{-21}\ \mathrm J
$$

On a molar scale\:

$$
RT\approx2.49\ \mathrm{kJ/mol}
$$

This is a useful reference value for thinking about molecular interactions.

The higher the temperature, the larger the typical thermal energy scale.

That generally means greater thermal fluctuations, although the exact structural response depends on the shape of the energy landscape and the molecular constraints.

## 4.3 An important distinction: thermal energy does not mean every atom has the same energy

The equipartition theorem describes average energy.

Individual atoms have different instantaneous velocities and kinetic energies.

At thermal equilibrium, the velocity distribution of classical particles follows the Maxwell-Boltzmann distribution.

For one Cartesian velocity component\:

$$
p(v_x)\propto
\exp\left(-\frac{mv_x^2}{2k_BT}\right)
$$

This is a Gaussian distribution centered around zero.

Some atoms move slowly, some move quickly, and a small fraction move much faster than average.

There is no single velocity that characterizes every atom.

This statistical behavior is essential to understanding molecular motion.

## 4.4 Why don't atoms simply fly apart?

This is where section 3 and section 4 connect.

Thermal motion tries to move atoms around. The potential energy surface determines which movements are energetically favorable or unfavorable.

Consider a ball in a bowl again.

The ball has kinetic energy and can move around.

However, the shape of the bowl constrains its movement.

If the ball has little kinetic energy, it stays near the bottom.

If it has more kinetic energy, it explores a larger region of the bowl.

For a molecule\:

- Thermal motion drives atomic displacement.
- The potential energy surface creates forces.
- Those forces influence the motion.
- The resulting motion explores different molecular configurations.

This is the fundamental interplay\:

$$
\boxed{
\text{Thermal motion}
+
\text{Interatomic forces}
=
\text{Molecular dynamics}
}
$$

Neither thermal motion nor the potential energy surface alone describes the complete behavior.

## 4.5 The role of energy barriers

Let's imagine a molecule with two possible conformations.

Conformation A is in one energy minimum.

Conformation B is in another.

To move from A to B, the molecule must pass through a higher-energy region.

The difference between the energy of the barrier and the initial minimum is called the activation energy or barrier height, depending on the context.

Let's call it\:

$$
\Delta U^\ddagger
$$

If the barrier is small compared with the thermal energy scale, transitions may occur relatively frequently.

If the barrier is large compared with the thermal energy scale, transitions may be rare.

A useful dimensionless comparison is\:

$$
\frac{\Delta U^\ddagger}{k_BT}
$$

For example\:

| Barrier | Relative to thermal energy |
|--------|--------|
| 1 kJ/mol | $0.4\,k_BT$ |
| 5 kJ/mol | $2.0\,k_BT$ |
| 10 kJ/mol | $4.0\,k_BT$ |
| 25 kJ/mol | $10.0\,k_BT$ |
| 50 kJ/mol | $20.1\,k_BT$ |

These are illustrative barrier values at approximately 300 K.

The larger this ratio, generally the less frequently thermal fluctuations can overcome the barrier.

But we should be careful: barrier height alone does not determine the transition rate. The shape of the energy landscape, friction, solvent interactions, and dynamical pathways also matter.

This is particularly important in proteins.

A protein may have two conformations that are both energetically accessible, but the transition between them may be slow because of a substantial barrier.

That is why a short molecular dynamics simulation may fail to observe certain conformational changes even if those conformations are physically relevant.

## 4.6 Thermal fluctuations are not random structural damage

This is an important conceptual point.

When we say that proteins fluctuate, we do not mean that their structures are randomly destroyed.

Their motion is constrained by\:

- Covalent bonds.
- Noncovalent interactions.
- Steric constraints.
- Solvent interactions.
- The overall molecular architecture.

A protein typically fluctuates around a set of structurally accessible conformations.

Some regions are relatively rigid. Others are flexible.

For example\:

- A well-packed protein core is often relatively constrained.
- Surface loops may be more flexible.
- Binding-site side chains may sample multiple orientations.
- Ligands may have flexible substituents that rotate.

These fluctuations are often functionally important.

Proteins are not merely rigid containers for drug molecules. Their flexibility can influence binding, selectivity, and biological activity.

---

# 5\. The crucial connection: the potential energy surface determines the fluctuations

Now let's combine both concepts.

Imagine we have a molecular system at temperature $T$\.

It has\:

- A potential energy function $U(\mathbf R)$\.
- A collection of possible configurations $\mathbf R$\.
- Thermal energy associated with its temperature.

The system does not explore all configurations equally.

At equilibrium, the probability density of observing a particular configuration in the canonical ensemble is\:

$$
\boxed{
p(\mathbf R)=
\frac{e^{-\beta U(\mathbf R)}}{Z_\mathrm{conf}}
}
$$

where\:

$$
\beta=\frac{1}{k_BT}
$$

and\:

$$
Z_\mathrm{conf}
=
\int e^{-\beta U(\mathbf R)}d\mathbf R
$$

This is the configurational Boltzmann distribution.

Let's interpret it.

### Low-energy configurations

If $U(\mathbf R)$ is low, then\:

$$
e^{-\beta U(\mathbf R)}
$$

is relatively large.

Such configurations receive greater statistical weight.

### High-energy configurations

If $U(\mathbf R)$ is high, the exponential factor is smaller.

These configurations receive less statistical weight.

### Temperature changes the distribution

At higher temperature, $\beta$ becomes smaller.

Energy differences have a weaker influence on the relative Boltzmann weights.

At lower temperature, energy differences have a stronger influence.

This does not mean that a molecule simply occupies its minimum-energy structure at low temperature. It means that higher-energy configurations become progressively less populated relative to lower-energy ones, under the specified equilibrium model.

## 5.1 A numerical example

Suppose two molecular configurations have energies\:

$$
U_A=0
$$

$$
U_B=5\ \mathrm{kJ/mol}
$$

At 300 K\:

$$
RT\approx2.49\ \mathrm{kJ/mol}
$$

The relative Boltzmann weight of B compared with A is\:

$$
\frac{w_B}{w_A}
=
e^{-\Delta U/(RT)}
$$

Therefore\:

$$
\frac{w_B}{w_A}
=
e^{-5/2.49}
\approx0.134
$$

So, for two individual states with equal statistical measure, B has about 13.4% of the weight of A.

If these were the only two possible states, their normalized probabilities would be approximately\:

$$
p_A=0.882
$$

$$
p_B=0.118
$$

But real molecules have continuous configurations, and the volume of configuration space associated with each region also matters.

That is the subtle connection to entropy.

---

# 6\. Why the lowest-energy configuration is not necessarily the most populated region

This is perhaps the most important point in the entire discussion.

Imagine two valleys on a potential energy landscape.

- Valley A is deep but very narrow.
- Valley B is slightly higher but much broader.

Valley A contains a very favorable minimum-energy configuration.

Valley B contains many configurations with similar, slightly higher energies.

Even though A has the lowest minimum, B may contain more total statistical probability.

Why?

Because the probability of a region is obtained by integrating its Boltzmann weight over the region\:

$$
P(\text{region})
=
\frac{
\int_{\text{region}}e^{-\beta U(\mathbf R)}d\mathbf R
}{
Z_\mathrm{conf}
}
$$

A broad region can accumulate a large probability because it contains many configurations.

This is how entropy enters naturally.

We do not need to think of entropy as a mysterious extra force that pushes molecules toward disorder.

Instead, we can understand it as a consequence of how much configuration space is available to the system and how that space is populated.

**Energy determines the statistical weight of each configuration. The number and distribution of configurations determine the total weight of a region.**

That is the essential connection between potential energy and entropy.

---

# 7\. What this means for a protein-ligand complex

Let's return to imatinib and ABL1.

Imagine two different ligands bind to the same protein.

Ligand A forms several favorable interactions\:

- Hydrogen bonds.
- Electrostatic interactions.
- Hydrophobic contacts.
- Dispersion interactions.

Ligand B forms fewer favorable interactions.

It might seem obvious that A should bind more strongly.

But consider what happens during binding.

## 7.1 The ligand loses conformational freedom

A free ligand in solution may rotate around several flexible bonds.

When it binds, some of those rotations become restricted.

It loses accessible conformational states.

This often corresponds to an unfavorable entropy contribution.

## 7.2 The protein may lose conformational freedom

The binding pocket may become more constrained.

Some protein side chains and loops may lose accessible conformations.

That can also be entropically unfavorable.

## 7.3 Water molecules reorganize

Before binding, the binding pocket may contain water molecules.

When a ligand enters, some of these waters may be displaced into bulk solvent.

Depending on their initial environment and the states they can access afterward, their release may be favorable or unfavorable to binding.

This is why water is not simply a background medium in molecular drug discovery.

Water is an active participant in the thermodynamics of binding.

## 7.4 The protein and ligand gain favorable interactions

The bound complex may gain favorable electrostatic, dispersion, hydrogen-bonding, and other interactions.

These contributions can make binding energetically favorable.

The final binding free energy reflects the combined result of all these changes.

We can express the idea schematically as\:

$$
\Delta G_\mathrm{bind}
=
\Delta H_\mathrm{bind}
-
T\Delta S_\mathrm{bind}
$$

This is a simplified thermodynamic decomposition, with the understanding that the enthalpy and entropy refer to the complete binding process, including solvent and conformational effects.

It explains why simply counting hydrogen bonds or measuring the interaction energy of a crystal structure cannot reliably predict binding affinity.

---

# 8\. One important distinction before we proceed: potential energy versus free energy

These two concepts are often confused.

| Potential energy | Free energy |
|--------|--------|
| Assigned to a molecular configuration | Describes a thermodynamic state or collection of configurations |
| Depends on atomic coordinates in a chosen model | Depends on temperature and the statistical distribution |
| Determines forces through its gradient | Determines equilibrium thermodynamic preferences |
| Can be calculated for one structure | Requires statistical information over accessible states |

For example\:

- A molecular mechanics force field calculates potential energy for a given structure.
- Molecular dynamics uses forces derived from that potential energy to generate a trajectory.
- Statistical mechanics connects the sampled configurations to thermodynamic properties.
- Free-energy methods estimate the thermodynamic differences between states.

This is the conceptual chain behind computational drug discovery.

---

# 9\. A small mathematical experiment: see the difference between energy and probability

Before writing molecular simulation code, we can use a very simple mathematical example.

Imagine a hypothetical one-dimensional molecular coordinate $x$, with potential energy\:

$$
U(x)=\frac{1}{2}kx^2
$$

This is a harmonic potential, similar to the approximation of a molecular bond vibrating near its equilibrium length.

The energy is lowest at $x=0$\.

As we move away from zero, the energy increases quadratically.

The equilibrium probability density is\:

$$
p(x)\propto
e^{-U(x)/(k_BT)}
$$

This produces a Gaussian distribution.

The interesting result is that the probability distribution has a width that depends on temperature and the stiffness of the potential\:

$$
\sigma_x^2=\frac{k_BT}{k}
$$

Here\:

- Higher temperature produces larger fluctuations.
- Higher stiffness $k$ produces smaller fluctuations.

Let's calculate and plot this conceptually.

**Explore a harmonic molecular coordinate**

Change temperature and stiffness to see how the probability distribution changes.

**Temperature**

## K

> **Slider**
>
> Current value unavailable.
>
> **Range:** 100 to 600
>
> **Step:** 10

**Potential stiffness**

## kJ mol⁻¹ nm⁻²

> **Slider**
>
> Current value unavailable.
>
> **Range:** 2 to 40
>
> **Step:** 1

---

### Potential energy

### Equilibrium probability shape

**Predicted fluctuation width**

# nm

Standard deviation, $\sigma_x=\sqrt{k_BT/k}$\. The coordinate is illustrative, not a fitted molecular bond.

The plots show relative shapes, not normalized probability densities. The potential and probability curves are calculated from the same harmonic model.

Notice the physical relationship\:

- Increasing temperature broadens the probability distribution.
- Increasing stiffness narrows the distribution.
- The minimum of the potential remains at the same position.
- The most probable position remains at the minimum in this symmetric harmonic example.

This is a very simple model, but it captures the essential relationship between potential energy and thermal fluctuations.

In real molecules, the energy surfaces are multidimensional and generally nonlinear, but the underlying principles remain the same.

---

# 10\. What I want you to take away

There are five fundamental ideas here\:

1\. **A potential energy surface is a function of molecular geometry.** It tells us the energetic cost of different arrangements of atoms.

2\. **Forces are gradients of potential energy.** The shape of the energy surface determines how atoms accelerate.

3\. **Temperature produces thermal fluctuations.** Atoms do not sit motionless at energy minima; they continually move and explore configurations.

4\. **The potential energy surface constrains those fluctuations.** Low-energy regions are statistically favored, while energy barriers can make transitions rare.

5\. **The equilibrium behavior depends on both energy and accessible configuration space.** This is the foundation of entropy and free energy.

The next important conceptual step is understanding why a molecular system doesn't just explore configurations according to energy alone, but according to a probability distribution over its entire phase space. That will take us naturally into microscopic states, ensembles, the Boltzmann distribution, and eventually the statistical mechanics underlying MD and FEP.

## Prompt:
11/9/2026, 3:16:23 PM

2. The architecture of our library
Eventually, our library will have the following architecture:
Moleculab – Public API

Create molecules, simulate, analyze and visualize

Structure

Atoms, bonds, residues, molecules

Force fields

Bonds, angles, torsions, nonbonded

Physics

Energy, forces, integration

Simulation

MD, thermostats, barostats

Thermodynamics

Ensembles, free energy, sampling

Quantum chemistry

Orbitals, electron density, DFT

Numerical foundation

Linear algebra · Numerical integration · Optimization · Statistics

A critical architectural decision: we will separate the physical model from the numerical algorithms.

3. Development roadmap
We'll build the library in small, testable milestones.
Milestone	Implementation	Physics learned
0	Units and numerical conventions	Dimensional analysis
1	Atoms and coordinates	Classical mechanics
2	Lennard-Jones potential	Interatomic interactions
3	Coulomb potential	Electrostatics
4	Molecular topology	Chemical connectivity
5	Bond and angle potentials	Molecular mechanics
6	Torsional potentials	Molecular conformations
7	Energy minimization	Optimization
8	Velocity Verlet	Numerical integration
9	Molecular dynamics	Newtonian mechanics
10	Temperature and thermostats	Statistical mechanics
11	Solvent and periodic boundaries	Condensed matter physics
12	Trajectory analysis	Statistical analysis
13	Protein-ligand interactions	Molecular recognition
14	Free energy	Statistical thermodynamics
15	FEP and thermodynamic integration	Alchemical transformations
16	Basic quantum chemistry	Schrödinger equation
17	Hartree-Fock and DFT foundations	Electronic structure
18	QM/MM concepts	Multiscale physics

We won't rush through the milestones. Each one should have a working implementation before we move on.
Python first, Rust later
I recommend this division:
Python	Rust
Derivations and experimentation	High-performance calculations
NumPy numerical arrays	Efficient memory management
Matplotlib plots	Parallel computation
Jupyter notebooks	Production-oriented engine
Rapid prototyping	Optimized numerical kernels

Eventually, we can expose Rust functions to Python using PyO3.

For now, Python makes it much easier to inspect intermediate calculations and understand the mathematics.
4. Establish our physical units
This is one of the most important decisions in computational chemistry.

We will use the following unit system:
Quantity	Unit
Distance	nm
Time	ps
Energy	kJ/mol
Mass	Da
Force	kJ/(mol nm)
Temperature	K
Charge	Elementary charge $e$

Why these units?

They are widely used in biomolecular simulation and are convenient for molecular dynamics.

For example:
$$
1\ \mathrm{nm}=10^{-9}\ \mathrm m
$$

$$
1\ \mathrm{ps}=10^{-12}\ \mathrm s
$$

$$
1\ \mathrm{Da}\approx1.66054\times10^{-27}\ \mathrm{kg}
$$

An interesting property of this unit system is that the numerical Newtonian acceleration equation becomes particularly convenient:
$$
\boxed{a_i=\frac{F_i}{m_i}}
$$

when force, mass and acceleration are expressed in our selected simulation units.

We will still implement explicit unit checks and conversion utilities later, because mixing units is a common source of serious simulation errors.

we can now start building the library based on what we have learned so far with potential, kinetic eneries, etc.

the code should be as simple and understandable as possible for educational purposes

## Response:
11/9/2026, 3:16:28 PM · gpt-6-mini

> **Building the first model**
>
> I’m starting with a tiny, testable 1D harmonic oscillator: particles store mass, position, and velocity; a potential supplies energy and force; the system computes kinetic, potential, and total energy.
>
> **Built the first model**
>
> 
>
> ****
>
> 
>
> Worked for 4m 24s

Excellent. We can now start building Moleculab, but I want to make one important adjustment to our implementation strategy\:

**We should not begin with a molecular force field. We should first build a tiny physics engine that can represent a particle, calculate its kinetic and potential energy, and calculate the force acting on it.**

That gives us the foundations of molecular mechanics and molecular dynamics without introducing unnecessary chemical complexity.

We'll make the code deliberately simple. Every function should correspond to a physical concept we understand.

# Moleculab - first implementation

(Milestones 0-1) (Educational prototype)

Our first implementation will support\:

1. Representing a particle.
2. Storing its mass, position, and velocity.
3. Calculating kinetic energy.
4. Defining a potential energy function.
5. Calculating potential energy.
6. Calculating force from potential energy.
7. Calculating total mechanical energy.

We'll use a simple harmonic potential as our first physical model.

This is not yet a molecular force field. It is a deliberately simple model that lets us understand the relationship between position, energy, and force.

---

# 1\. First, define what our library represents

Consider a single particle moving in three-dimensional space.

We need three fundamental properties.

**Mass**

$$
m
$$

Mass determines how strongly a particle responds to a given force.

**Position**

$$
\mathbf r =
\begin{pmatrix}
x\\y\\z
\end{pmatrix}
$$

Position describes where the particle is.

**Velocity**

$$
\mathbf v =
\begin{pmatrix}
v_x\\v_y\\v_z
\end{pmatrix}
$$

Velocity describes how its position changes over time.

Together, these properties describe the mechanical state of our particle\:

$$
\text{Particle state}=(m,\mathbf r,\mathbf v)
$$

For now, we will not introduce acceleration as a stored property. We can calculate it from force and mass whenever we need it.

## Our first object: `Particle`

In Python\:

```python
import numpy as np

class Particle:

    def __init__(self, mass, position, velocity):

        self.mass = float(mass)

        self.position = np.array(position, dtype=float)

        self.velocity = np.array(velocity, dtype=float)
```

Let's instantiate a particle\:

```python
particle = Particle(
    mass=12.0,
    position=[1.0, 0.0, 0.0],
    velocity=[0.1, 0.0, 0.0]
)
```

This describes a particle with\:

| Property | Value |
|--------|--------|
| Mass | 12 Da |
| X position | 1.0 nm |
| Y position | 0.0 nm |
| Z position | 0.0 nm |
| X velocity | 0.1 nm/ps |
| Y velocity | 0.0 nm/ps |
| Z velocity | 0.0 nm/ps |

We are using a three-dimensional representation from the beginning because real molecules exist in three-dimensional space.

---

# 2\. Establish our physical units

Before calculating energy, we need consistent units.

We will use your proposed unit system\:

| Quantity | Unit | Symbol |
|--------|--------|--------|
| Distance | Nanometre | nm |
| Time | Picosecond | ps |
| Energy | Kilojoule per mole | kJ/mol |
| Mass | Dalton | Da |
| Force | Kilojoule per mole per nanometre | kJ/(mol nm) |
| Temperature | Kelvin | K |
| Charge | Elementary charge | e |

Let's put these conventions in one place.

**File: `src/moleculab/foundation/units.py`**

```python
# Distance
DISTANCE_UNIT = "nm"

# Time
TIME_UNIT = "ps"

# Energy
ENERGY_UNIT = "kJ/mol"

# Mass
MASS_UNIT = "Da"

# Force
FORCE_UNIT = "kJ/(mol nm)"

# Temperature
TEMPERATURE_UNIT = "K"

# Charge
CHARGE_UNIT = "e"

# Physical constants

BOLTZMANN_CONSTANT = 0.008314462618
# kJ / (mol K)
```

Why store the Boltzmann constant?

Because later we will need the relationship\:

$$
k_BT
$$

to understand thermal energy, temperature, and statistical mechanics.

Notice that we express the Boltzmann constant per mole\:

$$
R = 0.008314462618\ \mathrm{kJ/(mol\,K)}
$$

This is numerically the molar gas constant, $R=N_Ak_B$, rather than the single-particle Boltzmann constant.

For molecular simulation, this molar form is convenient because our energy unit is kJ/mol.

---

# 3\. Kinetic energy

Now let's implement our first physical equation.

Kinetic energy is the energy associated with motion.

For a particle\:

$$
K=\frac12 mv^2
$$

In three dimensions, speed squared is\:

$$
v^2=v_x^2+v_y^2+v_z^2
$$

Therefore\:

$$
\boxed{
K=\frac12m(v_x^2+v_y^2+v_z^2)
}
$$

Or, using vector notation\:

$$
\boxed{
K=\frac12m\mathbf v\cdot\mathbf v
}
$$

The dot product is simply the sum of the squared velocity components.

## A useful property of our units

With mass in Da and velocity in nm/ps, the numerical value of this equation gives energy in kJ/mol.

This is particularly convenient.

For example, a particle with\:

$$
m=12\ \mathrm{Da}
$$

and\:

$$
\mathbf v=(0.1,0,0)\ \mathrm{nm/ps}
$$

has kinetic energy\:

$$
K=\frac12(12)(0.1)^2
$$

$$
\boxed{K=0.06\ \mathrm{kJ/mol}}
$$

No additional numerical conversion factor is needed.

## Implementing kinetic energy

**File: `src/moleculab/models/physics.py`**

```python
def kinetic_energy(particle):

    velocity_squared = np.dot(
        particle.velocity,
        particle.velocity
    )

    energy = 0.5 * particle.mass * velocity_squared

    return energy
```

We also need to import NumPy\:

```python
import numpy as np
```

Let's test it\:

```python
particle = Particle(
    mass=12.0,
    position=[1.0, 0.0, 0.0],
    velocity=[0.1, 0.0, 0.0]
)

print(kinetic_energy(particle))
```

Expected output\:

```text
0.06
```

We have now implemented our first physical calculation.

Notice something important: kinetic energy does not depend on position.

It depends only on mass and velocity.

That distinction will become important when we combine kinetic and potential energy.

---

# 4\. Potential energy: introduce a simple physical model

We now need a model that assigns potential energy to different particle positions.

For this purpose, we'll use the harmonic potential.

Imagine a particle attached to a spring.

If we displace the particle from its equilibrium position, the spring exerts a restoring force.

The farther we move it, the stronger the restoring force becomes.

For a spring\:

$$
F=-kx
$$

where\:

- $F$ is the force.
- $k$ is the spring stiffness.
- $x$ is the displacement from equilibrium.

The corresponding potential energy is\:

$$
\boxed{U(x)=\frac12kx^2}
$$

This equation is extremely useful because it connects three concepts\:

- Position determines potential energy.
- The derivative of potential energy determines force.
- Force determines acceleration.

## 4.1 Generalize to three dimensions

Instead of one coordinate, we use a position vector.

Let the equilibrium position be\:

$$
\mathbf r_0
$$

The displacement is\:

$$
\Delta\mathbf r=\mathbf r-\mathbf r_0
$$

The squared distance from equilibrium is\:

$$
|\Delta\mathbf r|^2
=
\Delta\mathbf r\cdot\Delta\mathbf r
$$

Our potential energy becomes\:

$$
\boxed{
U(\mathbf r)=
\frac12 k|\mathbf r-\mathbf r_0|^2
}
$$

This is the three-dimensional harmonic potential with equal stiffness in all directions.

For now, it is enough to think of a particle attached to an imaginary spring in three dimensions.

## 4.2 Why is this a good educational model?

It has a number of useful properties\:

- There is one obvious equilibrium position.
- The potential energy is zero at equilibrium.
- Energy increases as the particle moves away.
- The force always points toward equilibrium.
- The mathematics is simple.
- We can later use it to study oscillations and molecular vibrations.

It is also a useful approximation for the behavior of real molecular systems near a stable energy minimum.

However, real chemical bonds are not perfectly harmonic over arbitrary distances.

---

# 5\. Implement the potential as a separate physical model

Here we reach our first important architectural decision.

We want to separate\:

- The mathematical description of a physical interaction.
- The general physics calculations performed using that description.

For example, a harmonic potential defines an energy function.

A separate physics function can calculate the potential energy of a particle using that model.

Let's create a new file.

**File: `src/moleculab/models/potentials.py`**

```python
import numpy as np

class HarmonicPotential:

    def __init__(self, stiffness, center):

        self.stiffness = float(stiffness)

        self.center = np.array(center, dtype=float)

    def energy(self, position):

        displacement = position - self.center

        distance_squared = np.dot(
            displacement,
            displacement
        )

        return 0.5 * self.stiffness * distance_squared

    def force(self, position):

        displacement = position - self.center

        return -self.stiffness * displacement
```

Let's examine the implementation.

### The energy method

```python
displacement = position - self.center
```

This calculates how far the particle is from equilibrium.

```python
distance_squared = np.dot(
    displacement,
    displacement
)
```

This calculates the squared distance.

```python
return 0.5 * self.stiffness * distance_squared
```

This implements\:

$$
U=\frac12k|\mathbf r-\mathbf r_0|^2
$$

### The force method

```python
return -self.stiffness * displacement
```

This implements\:

$$
\boxed{\mathbf F=-k(\mathbf r-\mathbf r_0)}
$$

The negative sign ensures that the force points toward the equilibrium position.

The important thing is that the force and energy describe the same physical model.

We will later test that consistency mathematically.

---

# 6\. Connect the potential energy model to physics

We now have a potential that can calculate energy and force.

Let's connect this to the general equations we learned.

Potential energy\:

$$
U(\mathbf r)=\frac12k|\mathbf r-\mathbf r_0|^2
$$

Force\:

$$
\mathbf F=-\nabla U
$$

For our harmonic potential, taking the gradient gives\:

$$
\nabla U=k(\mathbf r-\mathbf r_0)
$$

Therefore\:

$$
\mathbf F=-k(\mathbf r-\mathbf r_0)
$$

This is exactly what our `force()` method implements.

This relationship is extremely important\:

**We do not independently invent an energy equation and a force equation. The force must be mathematically consistent with the energy function.**

In computational molecular physics, this is one of the most important correctness requirements.

## 6.1 Calculate acceleration

Newton's second law tells us\:

$$
\mathbf F=m\mathbf a
$$

Rearranging\:

$$
\boxed{\mathbf a=\frac{\mathbf F}{m}}
$$

We can implement this as a general physics function.

Add this to `models/physics.py`\:

```python
def acceleration(particle, force):

    return force / particle.mass
```

Now we can calculate the acceleration of a particle in our harmonic potential\:

```python
force = potential.force(particle.position)

a = acceleration(particle, force)
```

Notice our separation of responsibilities\:

- `HarmonicPotential` knows how to calculate its force.
- `acceleration()` knows Newton's second law.
- `Particle` stores the particle's physical properties.

This is the beginning of the separation between physical models and general physics algorithms.

---

# 7\. Potential energy, kinetic energy, and total energy

We now have two forms of mechanical energy.

**Kinetic energy\:**

$$
K=\frac12m|\mathbf v|^2
$$

**Potential energy\:**

$$
U=U(\mathbf r)
$$

The total mechanical energy is\:

$$
\boxed{E=K+U}
$$

This is the total energy of our particle in the harmonic potential, assuming no other energy contributions.

Let's implement this relationship.

Update `models/physics.py`\:

```python
import numpy as np

def kinetic_energy(particle):

    velocity_squared = np.dot(
        particle.velocity,
        particle.velocity
    )

    return 0.5 * particle.mass * velocity_squared

def potential_energy(particle, potential):

    return potential.energy(particle.position)

def total_energy(particle, potential):

    kinetic = kinetic_energy(particle)

    potential_e = potential_energy(
        particle,
        potential
    )

    return kinetic + potential_e

def acceleration(particle, force):

    return force / particle.mass
```

We now have a small collection of functions corresponding directly to physical equations.

| Function | Physical concept | Equation |
|--------|--------|--------|
| `kinetic_energy()` | Energy of motion | $K=\frac12mv^2$ |
| `potential_energy()` | Energy of configuration | $U=U(\mathbf r)$ |
| `total_energy()` | Total mechanical energy | $E=K+U$ |
| `acceleration()` | Newton's second law | $\mathbf a=\mathbf F/m$ |

---

# 8\. A subtle but important point: our units

Let's examine the units of the harmonic potential.

We defined\:

$$
U=\frac12kr^2
$$

Our energy is in kJ/mol, and distance is in nm.

Therefore, stiffness must have units\:

$$
[k]=\frac{\mathrm{kJ}}{\mathrm{mol\,nm^2}}
$$

The force is\:

$$
F=kr
$$

Its units are therefore\:

$$
[F]=\frac{\mathrm{kJ}}{\mathrm{mol\,nm}}
$$

This matches our chosen force unit.

Now consider acceleration\:

$$
a=\frac Fm
$$

Using our simulation units\:

$$
\frac{
\mathrm{kJ/(mol\,nm)}
}{
\mathrm{Da}
}
$$

This may initially look like an unusual acceleration unit.

But with the conversions between joules, moles, kilograms, nanometres, and picoseconds, it simplifies numerically to\:

$$
\boxed{\mathrm{nm/ps^2}}
$$

So, for our chosen units, the numerical equation is simply\:

$$
\boxed{a_{\mathrm{numeric}}=
\frac{F_{\mathrm{numeric}}}{m_{\mathrm{numeric}}}}
$$

This is why we can implement acceleration without an additional conversion constant.

It is a particularly convenient property of this unit system.

---

# 9\. Put the pieces together

Let's create a complete example.

We will place a particle at a position away from the equilibrium point of a harmonic potential.

The particle also has an initial velocity.

**File: `examples/first_particle.py`**

```python
import numpy as np

from moleculab.particles import Particle
from moleculab.potentials import HarmonicPotential
from moleculab.physics import (
    kinetic_energy,
    potential_energy,
    total_energy,
    acceleration,
)

# 1. Create a particle

particle = Particle(
    mass=12.0,
    position=[0.1, 0.0, 0.0],
    velocity=[0.1, 0.0, 0.0]
)

# 2. Define a potential energy model

potential = HarmonicPotential(
    stiffness=100.0,
    center=[0.0, 0.0, 0.0]
)

# 3. Calculate potential energy

U = potential_energy(particle, potential)

# 4. Calculate kinetic energy

K = kinetic_energy(particle)

# 5. Calculate total energy

E = total_energy(particle, potential)

# 6. Calculate force

F = potential.force(particle.position)

# 7. Calculate acceleration

a = acceleration(particle, F)

# 8. Display results

print("Kinetic energy:", K, "kJ/mol")

print("Potential energy:", U, "kJ/mol")

print("Total energy:", E, "kJ/mol")

print("Force:", F, "kJ/(mol nm)")

print("Acceleration:", a, "nm/ps^2")
```

Let's calculate the expected values by hand first.

Our particle has\:

$$
m=12\ \mathrm{Da}
$$

$$
\mathbf r=(0.1,0,0)\ \mathrm{nm}
$$

$$
\mathbf v=(0.1,0,0)\ \mathrm{nm/ps}
$$

Our harmonic potential has\:

$$
k=100\ \mathrm{kJ/(mol\,nm^2)}
$$

and\:

$$
\mathbf r_0=(0,0,0)
$$

**Kinetic energy\:**

$$
K=\frac12(12)(0.1)^2
$$

$$
\boxed{K=0.06\ \mathrm{kJ/mol}}
$$

**Potential energy\:**

$$
U=\frac12(100)(0.1)^2
$$

$$
\boxed{U=0.5\ \mathrm{kJ/mol}}
$$

**Total energy\:**

$$
E=0.06+0.5
$$

$$
\boxed{E=0.56\ \mathrm{kJ/mol}}
$$

**Force\:**

$$
F=-100(0.1)
$$

$$
\boxed{\mathbf F=(-10,0,0)}
$$

**Acceleration\:**

$$
a=\frac{-10}{12}
$$

$$
\boxed{\mathbf a=(-0.8333,0,0)\ \mathrm{nm/ps^2}}
$$

The negative sign tells us that acceleration points toward the equilibrium position.

This is exactly the physical behavior we expect from a spring.

---

# 10\. Our initial library structure

We can now organize the implementation into a small package.

```text
moleculab/
│
├── pyproject.toml
├── README.md
│
├── src/
│   └── moleculab/
│       │
│       ├── __init__.py
│       ├── units.py
│       ├── particles.py
│       ├── potentials.py
│       └── physics.py
│
├── examples/
│   └── first_particle.py
│
└── tests/
    ├── test_particles.py
    ├── test_potentials.py
    └── test_physics.py
```

Let's complete the two remaining basic files.

**File: `src/moleculab/foundation/particles.py`**

```python
import numpy as np

class Particle:

    def __init__(self, mass, position, velocity):

        if mass <= 0:
            raise ValueError("Mass must be positive.")

        self.mass = float(mass)

        self.position = np.asarray(
            position,
            dtype=float
        )

        self.velocity = np.asarray(
            velocity,
            dtype=float
        )

        if self.position.shape != (3,):
            raise ValueError("Position must have 3 components.")

        if self.velocity.shape != (3,):
            raise ValueError("Velocity must have 3 components.")
```

We added basic validation because a particle must have positive mass and three-dimensional position and velocity vectors.

**File: `src/moleculab/__init__.py`**

```python
from .particles import Particle

from .potentials import HarmonicPotential

from .physics import (
    kinetic_energy,
    potential_energy,
    total_energy,
    acceleration,
)
```

Now users can write\:

```python
from moleculab import Particle, HarmonicPotential
```

rather than importing from each individual module.

# 11\. Test the physics, not just the Python code

This is particularly important for a scientific library.

A function can execute without errors and still implement the wrong physics.

We therefore need tests that verify physical relationships.

We'll start with three tests.

## Test 1: Kinetic energy

**File: `tests/test_physics.py`**

```python
import numpy as np

from moleculab import Particle
from moleculab import kinetic_energy

def test_kinetic_energy():

    particle = Particle(
        mass=12.0,
        position=[0.0, 0.0, 0.0],
        velocity=[0.1, 0.0, 0.0]
    )

    energy = kinetic_energy(particle)

    assert np.isclose(energy, 0.06)
```

This verifies that our implementation follows\:

$$
K=\frac12mv^2
$$

## Test 2: Potential energy at equilibrium

**File: `tests/test_potentials.py`**

```python
import numpy as np

from moleculab import HarmonicPotential

def test_energy_at_equilibrium():

    potential = HarmonicPotential(
        stiffness=100.0,
        center=[0.0, 0.0, 0.0]
    )

    energy = potential.energy(
        np.array([0.0, 0.0, 0.0])
    )

    assert np.isclose(energy, 0.0)
```

At equilibrium\:

$$
\mathbf r=\mathbf r_0
$$

Therefore\:

$$
U=\frac12k(0)^2=0
$$

## Test 3: Force is the negative derivative of energy

This is our most interesting test.

We know\:

$$
F(x)=-\frac{dU}{dx}
$$

We can independently approximate the derivative numerically\:

$$
F_\mathrm{numerical}(x)
\approx
-\frac{U(x+h)-U(x-h)}{2h}
$$

Here $h$ is a very small displacement.

Let's test whether the force returned by our potential agrees with this numerical derivative.

```python
import numpy as np

from moleculab import HarmonicPotential

def test_force_matches_energy_derivative():

    potential = HarmonicPotential(
        stiffness=100.0,
        center=[0.0, 0.0, 0.0]
    )

    position = np.array([0.1, 0.0, 0.0])

    h = 1e-6

    position_plus = position + np.array([h, 0.0, 0.0])

    position_minus = position - np.array([h, 0.0, 0.0])

    energy_plus = potential.energy(position_plus)

    energy_minus = potential.energy(position_minus)

    numerical_force = -(energy_plus - energy_minus) / (2 * h)

    analytical_force = potential.force(position)[0]

    assert np.isclose(
        numerical_force,
        analytical_force,
        atol=1e-6
    )
```

This test is important because it checks that our energy and force functions describe the same physical system.

Run the tests with\:

```bash
pytest -v
```

---

# 12\. Install and run Moleculab

Our `pyproject.toml` can initially be very simple\:

```toml
[build-system]
requires = ["setuptools>=68"]
build-backend = "setuptools.build_meta"

[project]
name = "moleculab"
version = "0.1.0"
description = "An educational molecular physics library"
requires-python = ">=3.10"
dependencies = ["numpy"]

[tool.setuptools.packages.find]
where = ["src"]
```

From the root of the project\:

```bash
pip install -e .
pip install pytest
```

Then\:

```bash
pytest -v
```

And run our example\:

```bash
python examples/first_particle.py
```

# 13\. What have we actually built?

Let's take stock.

**Current capabilities** (Implemented)

Represent a particle in three-dimensional space.

Store mass, position, and velocity.

Calculate kinetic energy.

Define a harmonic potential.

Calculate potential energy.

Calculate forces from the potential.

Calculate acceleration.

Calculate total mechanical energy.

Test mathematical consistency of the physics.

---

**Not implemented yet**

- Multiple interacting particles.
- Molecular bonds.
- Molecular force fields.
- Time integration.
- Molecular dynamics.
- Temperature control.
- Statistical ensembles.

There is one important limitation to understand: we have not yet simulated anything. We can calculate the energy, force, and acceleration at a given instant, but we have not calculated what happens next.

That will be our next step.

# 14\. A crucial architectural distinction

Before moving on, let's revisit your original architectural principle.

We have already started separating physical models from general physics calculations.

| Component | Responsibility |
|--------|--------|
| `Particle` | Represents physical properties of a particle |
| `HarmonicPotential` | Defines one particular physical energy model |
| `kinetic_energy()` | Calculates kinetic energy |
| `potential_energy()` | Evaluates a potential model |
| `acceleration()` | Applies Newton's second law |
| `total_energy()` | Combines energy contributions |

Later, we can replace the harmonic potential with\:

- Lennard-Jones interactions.
- Coulomb interactions.
- Bond potentials.
- Angle potentials.
- Torsional potentials.

The particle representation and much of the physics infrastructure will remain unchanged.

That is the benefit of separating the physical model from the numerical algorithms.

One more distinction: our harmonic potential represents a particle in an external potential, not a chemical bond between two atoms. When we introduce molecular topology, we'll create a separate model for interactions between particles.

---

# 15\. What I propose we do next

I would not introduce Lennard-Jones yet.

Instead, I suggest the next milestone be\:

(Recommended next milestone)

**A particle moving in a harmonic potential.**

We already know its position, velocity, acceleration, kinetic energy, and potential energy. Now we can investigate how these quantities change over time.

We'll introduce the simplest numerical integration method needed to solve Newton's equation\:

$$
m\frac{d^2x}{dt^2}=-kx
$$

This is the harmonic oscillator.

We'll use it to understand\:

1. How position changes over time.
2. How velocity changes over time.
3. How kinetic and potential energy convert into each other.
4. Why total mechanical energy should remain constant in an ideal system.
5. Why numerical integration introduces errors.
6. How we can measure and control those errors.

This is an excellent bridge between elementary mechanics and molecular dynamics.

**My recommendation\:** finish and test this small foundation before adding more features. We want Moleculab to grow from a correct, understandable physical model-not from a collection of equations and code that we haven't fully understood.



---
Powered by [ChatGPT Exporter](https://www.chatgptexporter.com)