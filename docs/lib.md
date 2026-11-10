> **As-built architecture:** see [`architecture.md`](architecture.md) for the
> implemented layers, module dependencies, and decision records. This file is
> the long-term roadmap. For the three-track tutorial entry, see
> [`tutorials/README.md`](tutorials/README.md); for the CADD-focused course
> overview, see [`tutorials/introduction.md`](tutorials/introduction.md).

Building a molecular simulation library from scratch: From physics to drug discovery
I think this is an excellent way to learn computational drug discovery. Rather than treating GROMACS, OpenMM, ORCA or AutoDock as black boxes, we'll gradually build our own small molecular simulation engine and understand exactly what happens inside it.

We'll use Python initially, with NumPy for numerical calculations and Matplotlib for visualization. Later, we can rewrite performance-critical components in Rust.

The goal is not to reproduce a commercial molecular simulation package. It is to build a small, scientifically meaningful library whose components implement the same fundamental physical principles.

Our eventual application will be the imatinib–ABL1 protein-ligand system.

But we'll start with two atoms, derive the physics, write the code, test it against analytical results, and progressively build up to a protein-ligand simulation.
1. What are we going to build?
Let's call our educational library:
Moleculab
Educational project

A molecular physics and simulation engine

A Python library that gradually develops capabilities for:
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

Our design principle will be:

Every important physical equation should have a corresponding implementation, unit test and numerical experiment.

For example, when implementing the Lennard-Jones potential, we will not just write a function. We will:
1. Derive the potential mathematically.
2. Understand its physical meaning.
3. Derive its force analytically.
4. Implement both energy and force.
5. Compare the analytical force with a numerical derivative.
6. Plot the potential and force.
7. Investigate how changing its parameters changes molecular behavior.
This is how we will approach every major component.
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

For example:
- Lennard-Jones describes a physical interaction model.
- Velocity Verlet is a numerical integration algorithm.
- NVT describes a thermodynamic ensemble.
- A thermostat is an algorithm for sampling that ensemble.
Keeping these concepts separate makes the library easier to understand and extend.
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

$$ 1\ \mathrm{nm}=10^{-9}\ \mathrm m $$

$$ 1\ \mathrm{ps}=10^{-12}\ \mathrm s $$

$$ 1\ \mathrm{Da}\approx1.66054\times10^{-27}\ \mathrm{kg} $$

An interesting property of this unit system is that the numerical Newtonian acceleration equation becomes particularly convenient:

$$ a_i=\frac{F_i}{m_i} $$

when force, mass and acceleration are expressed in our selected simulation units.

We will still implement explicit unit checks and conversion utilities later, because mixing units is a common source of serious simulation errors.

