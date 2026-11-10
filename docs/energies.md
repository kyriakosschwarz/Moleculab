Excellent. We can now start building Moleculab, but I want to make one important adjustment to our implementation strategy:

We should not begin with a molecular force field. We should first build a tiny physics engine that can represent a particle, calculate its kinetic and potential energy, and calculate the force acting on it.

That gives us the foundations of molecular mechanics and molecular dynamics without introducing unnecessary chemical complexity.

We'll make the code deliberately simple. Every function should correspond to a physical concept we understand.

# Moleculab — first implementation



Milestones 0–1

Educational prototype



Our first implementation will support:

1. Representing a particle.
2. Storing its mass, position, and velocity.
3. Calculating kinetic energy.
4. Defining a potential energy function.
5. Calculating potential energy.
6. Calculating force from potential energy.
7. Calculating total mechanical energy.

We'll use a simple harmonic potential as our first physical model.

This is not yet a molecular force field. It is a deliberately simple model that lets us understand the relationship between position, energy, and force.

# 1. First, define what our library represents

Consider a single particle moving in three-dimensional space.

[4.1 Displacement and Velocity Vectors - University Physics Volume 1 | OpenStax](https://images.openai.com/static-rsc-4/dm1phUkG03qq97FDCfRiUjfQVIZTMatxNXfOwDeQsu3HCGMvJCNnkKz-FjesQjLJ1CRTWY-_kr6YfJiamdD5akrNCcZW_Xeki_HNDjijOCo7u-yld70GAWijCJx5lRCJWmodJZr4boCSrZYK5xKFOVAyJ1SvRiRjU3EyGmpFLa8?purpose=inline)

[openstax.org](https://openstax.org/books/university-physics-volume-1/pages/4-1-displacement-and-velocity-vectors?utm_source=chatgpt.com)

We need three fundamental properties.

Mass

\$$ m \ $$

Mass determines how strongly a particle responds to a given force.

Position

\$$ \mathbf r = \begin{pmatrix} x\\\y\\\z \end{pmatrix} \ $$

Position describes where the particle is.

Velocity

\$$ \mathbf v = \begin{pmatrix} v_x\\\v_y\\\v_z \end{pmatrix} \ $$

Velocity describes how its position changes over time.

Together, these properties describe the mechanical state of our particle:

\$$ \text{Particle state}=(m,\mathbf r,\mathbf v) \ $$

For now, we will not introduce acceleration as a stored property. We can calculate it from force and mass whenever we need it.

## Our first object: `Particle`

In Python:

```
import numpy as npclass Particle:    def __init__(self, mass, position, velocity):        self.mass = float(mass)        self.position = np.array(position, dtype=float)        self.velocity = np.array(velocity, dtype=float)
```

Let's instantiate a particle:

```
particle = Particle(    mass=12.0,    position=[1.0, 0.0, 0.0],    velocity=[0.1, 0.0, 0.0])
```

This describes a particle with:

| Property   | Value     |
| ---------- | --------- |
| Mass       | 12 Da     |
| X position | 1.0 nm    |
| Y position | 0.0 nm    |
| Z position | 0.0 nm    |
| X velocity | 0.1 nm/ps |
| Y velocity | 0.0 nm/ps |
| Z velocity | 0.0 nm/ps |

We are using a three-dimensional representation from the beginning because real molecules exist in three-dimensional space.

# 2. Establish our physical units

Before calculating energy, we need consistent units.

We will use your proposed unit system:

| Quantity    | Unit                             | Symbol      |
| ----------- | -------------------------------- | ----------- |
| Distance    | Nanometre                        | nm          |
| Time        | Picosecond                       | ps          |
| Energy      | Kilojoule per mole               | kJ/mol      |
| Mass        | Dalton                           | Da          |
| Force       | Kilojoule per mole per nanometre | kJ/(mol nm) |
| Temperature | Kelvin                           | K           |
| Charge      | Elementary charge                | e           |

Let's put these conventions in one place.

File: `src/moleculab/foundation/units.py`

```
# DistanceDISTANCE_UNIT = "nm"# TimeTIME_UNIT = "ps"# EnergyENERGY_UNIT = "kJ/mol"# MassMASS_UNIT = "Da"# ForceFORCE_UNIT = "kJ/(mol nm)"# TemperatureTEMPERATURE_UNIT = "K"# ChargeCHARGE_UNIT = "e"# Physical constantsBOLTZMANN_CONSTANT = 0.008314462618# kJ / (mol K)
```

Why store the Boltzmann constant?

Because later we will need the relationship:

\$$ k_BT \ $$

to understand thermal energy, temperature, and statistical mechanics.

Notice that we express the Boltzmann constant per mole:

\$$ R = 0.008314462618\ \mathrm{kJ/(mol\\ K)} \ $$

This is numerically the molar gas constant, \$R=N_Ak_B\$, rather than the single-particle Boltzmann constant.

For molecular simulation, this molar form is convenient because our energy unit is kJ/mol.

# 3. Kinetic energy

Now let's implement our first physical equation.

Kinetic energy is the energy associated with motion.

For a particle:

\$$ K=\frac12 mv^2 \ $$

In three dimensions, speed squared is:

\$$ v^2=v_x^2+v_y^2+v_z^2 \ $$

Therefore:

\$$  K=\frac12m(v_x^2+v_y^2+v_z^2)  \ $$

Or, using vector notation:

\$$  K=\frac12m\mathbf v\cdot\mathbf v  \ $$

The dot product is simply the sum of the squared velocity components.

## A useful property of our units

With mass in Da and velocity in nm/ps, the numerical value of this equation gives energy in kJ/mol.

This is particularly convenient.

For example, a particle with:

\$$ m=12\ \mathrm{Da} \ $$

and:

\$$ \mathbf v=(0.1,0,0)\ \mathrm{nm/ps} \ $$

has kinetic energy:

\$$ K=\frac12(12)(0.1)^2 \ $$

\$$ K=0.06\ \mathrm{kJ/mol} \ $$

No additional numerical conversion factor is needed.

## Implementing kinetic energy

File: `src/moleculab/models/physics.py`

```
def kinetic_energy(particle):    velocity_squared = np.dot(        particle.velocity,        particle.velocity    )    energy = 0.5 * particle.mass * velocity_squared    return energy
```

We also need to import NumPy:

```
import numpy as np
```

Let's test it:

```
particle = Particle(    mass=12.0,    position=[1.0, 0.0, 0.0],    velocity=[0.1, 0.0, 0.0])print(kinetic_energy(particle))
```

Expected output:

```
0.06
```

We have now implemented our first physical calculation.

Notice something important: kinetic energy does not depend on position.

It depends only on mass and velocity.

That distinction will become important when we combine kinetic and potential energy.

# 4. Potential energy: introduce a simple physical model

We now need a model that assigns potential energy to different particle positions.

For this purpose, we'll use the harmonic potential.

Imagine a particle attached to a spring.

[Oscillation Meaning Science at Thomas Lujan blog](https://images.openai.com/static-rsc-4/JZYMlTugFQcCAdx6boWD2PHWyb9NpIexJGARp4Ag3DBiGpuZM592XjHbQrm7w30hE_ElU8GCeAIjYLPVKL6zhtNPpUQqaNkM900a8VDOM1pl96WBn3v3bgmgmqN1QaZgb5XtNsOwIxTUCLgWT1FdG0MEZjlDzMYo50-a1D0pAxI?purpose=inline)

[storage.googleapis.com](https://storage.googleapis.com/dwmuaoorrdmlie/oscillation-meaning-science.html?utm_source=chatgpt.com)

If we displace the particle from its equilibrium position, the spring exerts a restoring force.

The farther we move it, the stronger the restoring force becomes.

For a spring:

\$$ F=-kx \ $$

where:

- \$F\$ is the force.
- \$k\$ is the spring stiffness.
- \$x\$ is the displacement from equilibrium.

The corresponding potential energy is:

\$$ U(x)=\frac12kx^2 \ $$

This equation is extremely useful because it connects three concepts:

- Position determines potential energy.
- The derivative of potential energy determines force.
- Force determines acceleration.

## 4.1 Generalize to three dimensions

Instead of one coordinate, we use a position vector.

Let the equilibrium position be:

\$$ \mathbf r_0 \ $$

The displacement is:

\$$ \Delta\mathbf r=\mathbf r-\mathbf r_0 \ $$

The squared distance from equilibrium is:

\$$ |\Delta\mathbf r|^2 = \Delta\mathbf r\cdot\Delta\mathbf r \ $$

Our potential energy becomes:

\$$  U(\mathbf r)= \frac12 k|\mathbf r-\mathbf r_0|^2  \ $$

This is the three-dimensional harmonic potential with equal stiffness in all directions.

For now, it is enough to think of a particle attached to an imaginary spring in three dimensions.

## 4.2 Why is this a good educational model?

It has a number of useful properties:

- There is one obvious equilibrium position.
- The potential energy is zero at equilibrium.
- Energy increases as the particle moves away.
- The force always points toward equilibrium.
- The mathematics is simple.
- We can later use it to study oscillations and molecular vibrations.

It is also a useful approximation for the behavior of real molecular systems near a stable energy minimum.

However, real chemical bonds are not perfectly harmonic over arbitrary distances.

# 5. Implement the potential as a separate physical model

Here we reach our first important architectural decision.

We want to separate:

- The mathematical description of a physical interaction.
- The general physics calculations performed using that description.

For example, a harmonic potential defines an energy function.

A separate physics function can calculate the potential energy of a particle using that model.

Let's create a new file.

File: `src/moleculab/models/potentials.py`

```
import numpy as npclass HarmonicPotential:    def __init__(self, stiffness, center):        self.stiffness = float(stiffness)        self.center = np.array(center, dtype=float)    def energy(self, position):        displacement = position - self.center        distance_squared = np.dot(            displacement,            displacement        )        return 0.5 * self.stiffness * distance_squared    def force(self, position):        displacement = position - self.center        return -self.stiffness * displacement
```

Let's examine the implementation.

### The energy method

```
displacement = position - self.center
```

This calculates how far the particle is from equilibrium.

```
distance_squared = np.dot(    displacement,    displacement)
```

This calculates the squared distance.

```
return 0.5 * self.stiffness * distance_squared
```

This implements:

\$$ U=\frac12k|\mathbf r-\mathbf r_0|^2 \ $$

### The force method

```
return -self.stiffness * displacement
```

This implements:

\$$ \mathbf F=-k(\mathbf r-\mathbf r_0) \ $$

The negative sign ensures that the force points toward the equilibrium position.

The important thing is that the force and energy describe the same physical model.

We will later test that consistency mathematically.

# 6. Connect the potential energy model to physics

We now have a potential that can calculate energy and force.

Let's connect this to the general equations we learned.

Potential energy:

\$$ U(\mathbf r)=\frac12k|\mathbf r-\mathbf r_0|^2 \ $$

Force:

\$$ \mathbf F=-\nabla U \ $$

For our harmonic potential, taking the gradient gives:

\$$ \nabla U=k(\mathbf r-\mathbf r_0) \ $$

Therefore:

\$$ \mathbf F=-k(\mathbf r-\mathbf r_0) \ $$

This is exactly what our `force()` method implements.

This relationship is extremely important:

We do not independently invent an energy equation and a force equation. The force must be mathematically consistent with the energy function.

In computational molecular physics, this is one of the most important correctness requirements.

## 6.1 Calculate acceleration

Newton's second law tells us:

\$$ \mathbf F=m\mathbf a \ $$

Rearranging:

\$$ \mathbf a=\frac{\mathbf F}{m} \ $$

We can implement this as a general physics function.

Add this to `models/physics.py`:

```
def acceleration(particle, force):    return force / particle.mass
```

Now we can calculate the acceleration of a particle in our harmonic potential:

```
force = potential.force(particle.position)a = acceleration(particle, force)
```

Notice our separation of responsibilities:

- `HarmonicPotential` knows how to calculate its force.
- `acceleration()` knows Newton's second law.
- `Particle` stores the particle's physical properties.

This is the beginning of the separation between physical models and general physics algorithms.

# 7. Potential energy, kinetic energy, and total energy

We now have two forms of mechanical energy.

Kinetic energy:

\$$ K=\frac12m|\mathbf v|^2 \ $$

Potential energy:

\$$ U=U(\mathbf r) \ $$

The total mechanical energy is:

\$$ E=K+U \ $$

This is the total energy of our particle in the harmonic potential, assuming no other energy contributions.

Let's implement this relationship.

Update `models/physics.py`:

```
import numpy as npdef kinetic_energy(particle):    velocity_squared = np.dot(        particle.velocity,        particle.velocity    )    return 0.5 * particle.mass * velocity_squareddef potential_energy(particle, potential):    return potential.energy(particle.position)def total_energy(particle, potential):    kinetic = kinetic_energy(particle)    potential_e = potential_energy(        particle,        potential    )    return kinetic + potential_edef acceleration(particle, force):    return force / particle.mass
```

We now have a small collection of functions corresponding directly to physical equations.

| Function             | Physical concept        | Equation                    |
| -------------------- | ----------------------- | --------------------------- |
| `kinetic_energy()`   | Energy of motion        | \$K=\frac12mv^2\$         |
| `potential_energy()` | Energy of configuration | \$U=U(\mathbf r)\$        |
| `total_energy()`     | Total mechanical energy | \$E=K+U\$                 |
| `acceleration()`     | Newton's second law     | \$\mathbf a=\mathbf F/m\$ |

# 8. A subtle but important point: our units

Let's examine the units of the harmonic potential.

We defined:

\$$ U=\frac12kr^2 \ $$

Our energy is in kJ/mol, and distance is in nm.

Therefore, stiffness must have units:

\$$ [k]=\frac{\mathrm{kJ}}{\mathrm{mol\\ nm^2}} \ $$

The force is:

\$$ F=kr \ $$

Its units are therefore:

\$$ [F]=\frac{\mathrm{kJ}}{\mathrm{mol\\ nm}} \ $$

This matches our chosen force unit.

Now consider acceleration:

\$$ a=\frac Fm \ $$

Using our simulation units:

\$$ \frac{ \mathrm{kJ/(mol\\ nm)} }{ \mathrm{Da} } \ $$

This may initially look like an unusual acceleration unit.

But with the conversions between joules, moles, kilograms, nanometres, and picoseconds, it simplifies numerically to:

\$$ \mathrm{nm/ps^2} \ $$

So, for our chosen units, the numerical equation is simply:

\$$ a\_{\mathrm{numeric}}= \frac{F\_{\mathrm{numeric}}}{m\_{\mathrm{numeric}}} \ $$

This is why we can implement acceleration without an additional conversion constant.

It is a particularly convenient property of this unit system.

# 9. Put the pieces together

Let's create a complete example.

We will place a particle at a position away from the equilibrium point of a harmonic potential.

The particle also has an initial velocity.

File: `examples/first_particle.py`

```
import numpy as npfrom moleculab.foundation.particles import Particlefrom moleculab.models.potentials import HarmonicPotentialfrom moleculab.models.physics import (    kinetic_energy,    potential_energy,    total_energy,    acceleration,)# 1. Create a particleparticle = Particle(    mass=12.0,    position=[0.1, 0.0, 0.0],    velocity=[0.1, 0.0, 0.0])# 2. Define a potential energy modelpotential = HarmonicPotential(    stiffness=100.0,    center=[0.0, 0.0, 0.0])# 3. Calculate potential energyU = potential_energy(particle, potential)# 4. Calculate kinetic energy
```

Let's calculate the expected values by hand first.

Our particle has:

\$$ m=12\ \mathrm{Da} \ $$

\$$ \mathbf r=(0.1,0,0)\ \mathrm{nm} \ $$

\$$ \mathbf v=(0.1,0,0)\ \mathrm{nm/ps} \ $$

Our harmonic potential has:

\$$ k=100\ \mathrm{kJ/(mol\\ nm^2)} \ $$

and:

\$$ \mathbf r_0=(0,0,0) \ $$

Kinetic energy:

\$$ K=\frac12(12)(0.1)^2 \ $$

\$$ K=0.06\ \mathrm{kJ/mol} \ $$

Potential energy:

\$$ U=\frac12(100)(0.1)^2 \ $$

\$$ U=0.5\ \mathrm{kJ/mol} \ $$

Total energy:

\$$ E=0.06+0.5 \ $$

\$$ E=0.56\ \mathrm{kJ/mol} \ $$

Force:

\$$ F=-100(0.1) \ $$

\$$ \mathbf F=(-10,0,0) \ $$

Acceleration:

\$$ a=\frac{-10}{12} \ $$

\$$ \mathbf a=(-0.8333,0,0)\ \mathrm{nm/ps^2} \ $$

The negative sign tells us that acceleration points toward the equilibrium position.

This is exactly the physical behavior we expect from a spring.

# 10. Our initial library structure

We can now organize the implementation into a small package.

```
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

File: `src/moleculab/foundation/particles.py`

```
import numpy as npclass Particle:    def __init__(self, mass, position, velocity):        if mass <= 0:            raise ValueError("Mass must be positive.")        self.mass = float(mass)        self.position = np.asarray(            position,            dtype=float        )        self.velocity = np.asarray(            velocity,            dtype=float        )        if self.position.shape != (3,):            raise ValueError("Position must have 3 components.")        if self.velocity.shape != (3,):            raise ValueError("Velocity must have 3 components.")
```

We added basic validation because a particle must have positive mass and three-dimensional position and velocity vectors.

File: `src/moleculab/__init__.py`

```
from .particles import Particlefrom .potentials import HarmonicPotentialfrom .physics import (    kinetic_energy,    potential_energy,    total_energy,    acceleration,)
```

Now users can write:

```
from moleculab import Particle, HarmonicPotential
```

rather than importing from each individual module.

# 11. Test the physics, not just the Python code

This is particularly important for a scientific library.

A function can execute without errors and still implement the wrong physics.

We therefore need tests that verify physical relationships.

We'll start with three tests.

## Test 1: Kinetic energy

File: `tests/test_physics.py`

```
import numpy as npfrom moleculab import Particlefrom moleculab import kinetic_energydef test_kinetic_energy():    particle = Particle(        mass=12.0,        position=[0.0, 0.0, 0.0],        velocity=[0.1, 0.0, 0.0]    )    energy = kinetic_energy(particle)    assert np.isclose(energy, 0.06)
```

This verifies that our implementation follows:

\$$ K=\frac12mv^2 \ $$

## Test 2: Potential energy at equilibrium

File: `tests/test_potentials.py`

```
import numpy as npfrom moleculab import HarmonicPotentialdef test_energy_at_equilibrium():    potential = HarmonicPotential(        stiffness=100.0,        center=[0.0, 0.0, 0.0]    )    energy = potential.energy(        np.array([0.0, 0.0, 0.0])    )    assert np.isclose(energy, 0.0)
```

At equilibrium:

\$$ \mathbf r=\mathbf r_0 \ $$

Therefore:

\$$ U=\frac12k(0)^2=0 \ $$

## Test 3: Force is the negative derivative of energy

This is our most interesting test.

We know:

\$$ F(x)=-\frac{dU}{dx} \ $$

We can independently approximate the derivative numerically:

\$$ F\_\mathrm{numerical}(x) \approx -\frac{U(x+h)-U(x-h)}{2h} \ $$

Here \$h\$ is a very small displacement.

Let's test whether the force returned by our potential agrees with this numerical derivative.

```
import numpy as npfrom moleculab import HarmonicPotentialdef test_force_matches_energy_derivative():    potential = HarmonicPotential(        stiffness=100.0,        center=[0.0, 0.0, 0.0]    )    position = np.array([0.1, 0.0, 0.0])    h = 1e-6    position_plus = position + np.array([h, 0.0, 0.0])    position_minus = position - np.array([h, 0.0, 0.0])    energy_plus = potential.energy(position_plus)    energy_minus = potential.energy(position_minus)    numerical_force = -(energy_plus - energy_minus) / (2 * h)    analytical_force = potential.force(position)[0]    assert np.isclose(        numerical_force,        analytical_force,        atol=1e-6    )
```

This test is important because it checks that our energy and force functions describe the same physical system.

Run the tests with:

```
pytest -v
```

# 12. Install and run Moleculab

Our `pyproject.toml` can initially be very simple:

```
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

From the root of the project:

```
pip install -e .
pip install pytest
```

Then:

```
pytest -v
```

And run our example:

```
python examples/first_particle.py
```

# 13. What have we actually built?

Let's take stock.

Current capabilities

Implemented

- Represent a particle in three-dimensional space.
- Store mass, position, and velocity.
- Calculate kinetic energy.
- Define a harmonic potential.
- Calculate potential energy.
- Calculate forces from the potential.
- Calculate acceleration.
- Calculate total mechanical energy.
- Test mathematical consistency of the physics.

Not implemented yet

- Multiple interacting particles.
- Molecular bonds.
- Molecular force fields.
- Time integration.
- Molecular dynamics.
- Temperature control.
- Statistical ensembles.

There is one important limitation to understand: we have not yet simulated anything. We can calculate the energy, force, and acceleration at a given instant, but we have not calculated what happens next.

That will be our next step.

# 14. A crucial architectural distinction

Before moving on, let's revisit your original architectural principle.

We have already started separating physical models from general physics calculations.

| Component            | Responsibility                               |
| -------------------- | -------------------------------------------- |
| `Particle`           | Represents physical properties of a particle |
| `HarmonicPotential`  | Defines one particular physical energy model |
| `kinetic_energy()`   | Calculates kinetic energy                    |
| `potential_energy()` | Evaluates a potential model                  |
| `acceleration()`     | Applies Newton's second law                  |
| `total_energy()`     | Combines energy contributions                |

Later, we can replace the harmonic potential with:

- Lennard-Jones interactions.
- Coulomb interactions.
- Bond potentials.
- Angle potentials.
- Torsional potentials.

The particle representation and much of the physics infrastructure will remain unchanged.

That is the benefit of separating the physical model from the numerical algorithms.

One more distinction: our harmonic potential represents a particle in an external potential, not a chemical bond between two atoms. When we introduce molecular topology, we'll create a separate model for interactions between particles.

# 15. What I propose we do next

I would not introduce Lennard-Jones yet.

Instead, I suggest the next milestone be:

Recommended next milestone

A particle moving in a harmonic potential.

We already know its position, velocity, acceleration, kinetic energy, and potential energy. Now we can investigate how these quantities change over time.

We'll introduce the simplest numerical integration method needed to solve Newton's equation:

\$$ m\frac{d^2x}{dt^2}=-kx \ $$

This is the harmonic oscillator.

We'll use it to understand:

1. How position changes over time.
2. How velocity changes over time.
3. How kinetic and potential energy convert into each other.
4. Why total mechanical energy should remain constant in an ideal system.
5. Why numerical integration introduces errors.
6. How we can measure and control those errors.

This is an excellent bridge between elementary mechanics and molecular dynamics.

My recommendation: finish and test this small foundation before adding more features. We want Moleculab to grow from a correct, understandable physical model—not from a collection of equations and code that we haven't fully understood.