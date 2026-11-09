# 03 — Potential energy, force, and the energy landscape

> **Where we are.** Tutorial 02 described motion and its kinetic energy. Now we
> ask where motion comes from: a potential energy function and the force that
> points downhill on it.
>
> **What you should already know.** Vectors, kinetic energy, and the unit
> conventions from tutorials 01–02.
>
> **What you will be able to do afterwards.** Read a potential energy curve,
> derive the force from it, explain why force points downhill, verify energy and
> force consistency numerically, and recognize the harmonic model as the local
> approximation behind many molecular interactions.

## 1. The central idea: a landscape of energy

A potential energy function assigns an energy to every possible configuration
of the system:

$$
U=U(\mathbf r_1,\mathbf r_2,\dots).
$$

For one particle in one dimension, we can draw $U(x)$ as a curve. For one
particle in two dimensions, it is a surface. For many particles, it is a
high-dimensional surface that we cannot draw, but the same intuition applies.

A ball on a hilly landscape is the standard analogy:

- at a valley bottom, the ball is stable;
- on a slope, it accelerates downhill;
- the steeper the slope, the stronger the acceleration;
- a hilltop is an unstable equilibrium.

The force is the downhill direction of the landscape. This is the single most
important idea in molecular mechanics.

## 2. Conservative forces and work

A force is **conservative** if the work it does depends only on the endpoints,
not on the path. Gravity and spring forces are conservative. Friction is not.

For a conservative force, we can define a potential energy $U$ such that

$$
W=-\Delta U.
$$

If the particle moves from $A$ to $B$,

$$
W_{A\to B}=U(A)-U(B).
$$

The force does positive work when the potential energy decreases. A ball
falling downward loses potential energy and gains kinetic energy.

In one dimension,

$$
W=\int_{x_A}^{x_B}F(x)\,dx.
$$

Comparing with $W=U(A)-U(B)$ gives

$$
F(x)=-\frac{dU}{dx}.
$$

In three dimensions,

$$
\mathbf F=-\nabla U.
$$

The gradient points uphill. The negative gradient points downhill. Force points
toward lower potential energy.

![Work and potential energy](assets/work_energy.svg)

*The shaded area under the force curve is the work done. The same area appears
as the change in potential energy. Force and energy are two views of one
landscape.*

## 3. The harmonic potential

The simplest landscape is a parabola:

$$
U(\mathbf r)=\frac12k|\mathbf r-\mathbf r_0|^2.
$$

Here:

- $k$ is the stiffness, in kJ/(mol nm²);
- $\mathbf r_0$ is the equilibrium position;
- $|\mathbf r-\mathbf r_0|$ is the distance from equilibrium.

The force is

$$
\mathbf F
=
-\nabla U
=
-k(\mathbf r-\mathbf r_0).
$$

This is Hooke's law in three dimensions. The force is:

- zero at equilibrium;
- proportional to displacement;
- directed back toward equilibrium.

The potential energy is a bowl. The force always points toward the bottom of the
bowl.

![Harmonic energy and force](assets/harmonic_energy_force.svg)

*Left: the energy bowl. Right: the force, which is the negative slope of the
bowl. At $x=0.1$ nm, $U=0.5$ kJ/mol and $F=-10$ kJ/(mol nm).*

### 3.1 Why the stiffness controls the shape

![Harmonic potentials for several stiffnesses](assets/harmonic_curvature.svg)

*A larger $k$ gives a narrower, steeper bowl. For the same displacement, the
force and the curvature are larger.*

The second derivative of $U$ is

$$
\frac{d^2U}{dx^2}=k.
$$

So the stiffness is literally the curvature of the energy landscape.

## 4. Why a parabola is more than a toy

Real molecular potentials are not exactly harmonic. But near a stable minimum,
every smooth potential looks approximately parabolic.

Expand $U(x)$ around a minimum $x_0$:

$$
U(x)=U(x_0)
+U'(x_0)(x-x_0)
+\frac12U''(x_0)(x-x_0)^2
+\cdots.
$$

At a minimum, $U'(x_0)=0$, so

$$
U(x)\approx U(x_0)
+\frac12k(x-x_0)^2,
\qquad
k=U''(x_0).
$$

The constant $U(x_0)$ does not affect force. The quadratic term does. This is
why harmonic potentials describe:

- bond vibrations near equilibrium;
- angle bending near the preferred angle;
- the bottom of a ligand binding pocket;
- small oscillations of a protein around a stable conformation.

The harmonic model breaks down far from the minimum. Real bonds dissociate, and
real nonbonded interactions have a finite well depth and a long-range tail.
But locally, the parabola is the right first approximation.

## 5. Worked example

Take

$$
k=100\ \mathrm{kJ/(mol\,nm^2)},
\qquad
\mathbf r=(0.1,0,0)\ \mathrm{nm},
\qquad
\mathbf r_0=(0,0,0).
$$

Potential energy:

$$
U=\frac12(100)(0.1)^2=0.5\ \mathrm{kJ/mol}.
$$

Force:

$$
\mathbf F=-100(0.1,0,0)
=(-10,0,0)\ \mathrm{kJ/(mol\,nm)}.
$$

For a particle of mass 12 Da,

$$
\mathbf a=\frac{\mathbf F}{m}
=\frac{(-10,0,0)}{12}
=(-0.8333,0,0)\ \mathrm{nm/ps^2}.
$$

The acceleration points back toward the origin. If the particle is displaced to
the right, it accelerates left; if displaced left, it accelerates right.

## 6. Energy landscapes in two dimensions

![Harmonic energy landscape and force arrows](assets/harmonic_quiver.svg)

*Contours of constant energy in two dimensions. The force arrows point
perpendicular to the contours, toward the center. The center is the minimum.*

The two-dimensional picture makes the gradient intuitive. The force is
perpendicular to the level curves and points from high energy to low energy.
At the center, the gradient is zero and the force vanishes.

This is the geometry behind energy minimization: to find a stable structure,
walk downhill until the gradient is zero.

## 7. The zero of energy is arbitrary

You can add a constant to $U$:

$$
U'(x)=U(x)+C.
$$

Then

$$
\frac{dU'}{dx}=\frac{dU}{dx},
$$

so the force is unchanged. Only energy differences matter for force and
dynamics.

This is an early glimpse of a deep idea in thermodynamics: absolute energies
are often meaningless, while energy differences and free-energy differences are
physical. Binding free energy is a difference, not an absolute energy.

## 8. Force as a numerical derivative

We can test whether a force implementation agrees with its energy
implementation by finite differences:

$$
F(x)\approx-\frac{U(x+h)-U(x-h)}{2h}.
$$

The central difference cancels the first-order error, leaving an error that
scales as $h^2$ for smooth functions.

![Finite-difference check of force](assets/force_derivative.svg)

*The tangent at a point is the derivative of the energy. A central difference
uses two nearby points to estimate that tangent.*

There is a tradeoff in choosing $h$:

- too large: truncation error;
- too small: floating-point cancellation.

![Finite-difference error versus step size](assets/finite_difference_error.svg)

*For this harmonic potential, the error is minimized around $h\sim10^{-5}$ nm.
Making $h$ smaller than that makes the result worse because of roundoff.*

The test in [`tests/test_potentials.py`](../../../tests/test_potentials.py)
uses a small step such as $h=10^{-6}$ nm and compares forces with a tolerance.
It is not a formality: a sign error or a missing factor in the force would be
caught immediately.

## 9. Code walkthrough

The harmonic model is in
[`src/moleculab/potentials.py`](../../../src/moleculab/potentials.py):

```python
from dataclasses import dataclass, field
import numpy as np

@dataclass
class HarmonicPotential:
    stiffness: float
    center: np.ndarray = field(default_factory=lambda: np.zeros(3))

    def energy(self, position):
        displacement = np.asarray(position) - self.center
        return 0.5 * self.stiffness * np.dot(
            displacement, displacement
        )

    def force(self, position):
        displacement = np.asarray(position) - self.center
        return -self.stiffness * displacement
```

The general physics functions are in
[`src/moleculab/physics.py`](../../../src/moleculab/physics.py):

```python
def potential_energy(particle, potential):
    return potential.energy(particle.position)

def acceleration(particle, force):
    return force / particle.mass
```

Why this separation?

| Component | Responsibility |
| --- | --- |
| `Particle` | stores state |
| `HarmonicPotential` | defines the landscape |
| `potential_energy()` | evaluates the landscape |
| `acceleration()` | applies Newton's second law |

A potential does not need to know about time integration. An integrator does
not need to know whether the landscape is harmonic or Lennard-Jones. This makes
it possible to add new physical models without rewriting the dynamics.

Try the code:

```bash
uv run python - <<'PY'
from moleculab import (
    HarmonicPotential, Particle,
    potential_energy, acceleration,
)

particle = Particle(12.0, [0.1, 0.0, 0.0], [0.0, 0.0, 0.0])
potential = HarmonicPotential(100.0, [0.0, 0.0, 0.0])
force = potential.force(particle.position)
print("U =", potential_energy(particle, potential))
print("F =", force)
print("a =", acceleration(particle, force))
PY
```

Expected:

```text
U = 0.5
F = [-10.   0.   0.]
a = [-0.83333333  0.          0.        ]
```

Tests:

```bash
uv run pytest tests/test_potentials.py tests/test_physics.py -q
```

## 10. External, bonded, and nonbonded potentials

The harmonic potential in this tutorial is an **external** potential: it
describes one particle in a specified energy landscape.

Later, Moleculab will use several different kinds of interaction:

- **bonded:** a spring between two specific atoms;
- **angle:** a harmonic penalty for bending an angle;
- **torsion:** a periodic energy for rotation around a bond;
- **nonbonded:** Lennard-Jones and Coulomb interactions between many pairs.

They all follow the same principle: energy defines force, and force is the
negative gradient. But their inputs and scopes differ. Keeping them separate is
an architectural choice, not a mathematical one.

## 11. Why energy–force consistency is non-negotiable

A simulation can run with an inconsistent force. It will simply produce wrong
physics. Common failure modes:

- **Sign error:** the force points uphill, so the system explodes.
- **Missing factor:** the oscillator frequency is wrong by a square-root factor.
- **Wrong center:** equilibrium is displaced.
- **Force from a different model:** energy conservation fails even if the
  trajectory looks plausible.

The finite-difference test is a direct check that the force is the derivative of
the energy. It should be part of the definition of every new potential.

## 12. Numerical experiments

### Experiment 1: slope and force

Pick $x=0.1$ nm and compute $U(x-h)$, $U(x)$, $U(x+h)$ for $h=10^{-4}$ nm.
Estimate the slope and compare with the analytical force.

### Experiment 2: stiffness

For $k=50$, $100$, and $200$, compute the force at the same displacement. What
is the ratio? What happens to the curvature?

### Experiment 3: energy zero

Add 100 kJ/mol to the potential energy at every position. Does the force
change? Does the trajectory change? Why?

### Experiment 4: finite-difference step

Vary $h$ from $10^{-8}$ to $10^{-1}$ nm and plot the error in the numerical
force. Can you reproduce the U-shaped curve in the figure?

## 13. Common mistakes

1. **Using $+\nabla U$ instead of $-\nabla U$.** Force points downhill.
2. **Forgetting the equilibrium center.** Use $\mathbf r-\mathbf r_0$.
3. **Mixing stiffness units.** $k$ must be kJ/(mol nm²) for this equation.
4. **Assuming the harmonic model is valid everywhere.** It is a local
   approximation.
5. **Confusing an external potential with a chemical bond.** A bond involves
   two atoms and a relative coordinate.
6. **Testing only the energy.** Always test the force derivative too.
7. **Thinking the absolute energy matters.** Only differences affect force.

## 14. Connections to CADD

Potential energy is the heart of molecular mechanics:

- a force field is a sum of potential energy terms;
- docking scores are interaction energies;
- energy minimization finds low-energy structures;
- MD integrates the forces from these potentials;
- free-energy methods compare the thermodynamics of different potentials.

When a docking score or a binding energy seems implausible, ask:

- what potential terms were included?
- what units were used?
- is the force consistent with the energy?
- is the model valid in this configuration?

## 15. Self-check with answers

1. **For $k=200$ and displacement 0.05 nm, what are $U$ and $|F|$?**
   $U=\tfrac12(200)(0.05)^2=0.25$ kJ/mol; $|F|=200(0.05)=10$ kJ/(mol nm).

2. **In which direction does the force point when $x>0$?**
   Toward equilibrium, so in the $-x$ direction.

3. **Why does adding a constant to $U$ not change the force?**
   The derivative of a constant is zero.

4. **What happens to curvature and force if $k$ doubles?**
   Both double for the same displacement.

5. **Why is a parabola a reasonable local bond model?**
   Near a smooth minimum, the first derivative vanishes and the leading
   non-constant term is quadratic.

6. **What does a nonzero force at equilibrium indicate?**
   Either the chosen point is not the true minimum, or the force implementation
   is inconsistent with the energy.

7. **Why is the finite-difference test important?**
   It independently checks the derivative relationship between energy and
   force.

## 16. Exercises

1. Derive the force for $U(x)=\tfrac12k(x-x_0)^2$ by hand.
2. Write a small Python script that verifies the harmonic force with a central
   difference for several positions.
3. Plot $U(x)$ and $F(x)$ for two different values of $k$.
4. Explain why a harmonic bond cannot describe bond breaking.
5. Read `HarmonicPotential` and list every place where a unit or a sign
   convention is encoded.

## 17. Summary

- Potential energy is a landscape; force is its negative gradient.
- The harmonic potential is the universal local model near a stable minimum.
- Stiffness is curvature.
- Energy and force must be derived from the same model and tested together.
- The zero of energy is arbitrary; differences matter.
- These ideas will become force fields, minimization, dynamics, and docking.

**Previous:** [02 — Particles and kinetic energy](02_particle_energy.md).  
**Next:** [04 — Velocity Verlet and energy conservation](04_velocity_verlet.md).
