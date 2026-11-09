# 02 — Particles, velocity, and kinetic energy

> **Where we are.** Tutorial 01 established the unit system. We now use those
> units to describe the simplest possible physical object: one classical point
> particle.
>
> **What you should already know.** Vectors, components, and the unit
> conventions from tutorial 01.
>
> **What you will be able to do afterwards.** Describe a particle's state,
> convert between velocity and speed, derive and compute kinetic energy, and
> explain why a simulation must store velocity but not acceleration.

## 1. From a molecule to a particle model

A real molecule is complicated. It has nuclei, electrons, bonds, vibrations,
and quantum behavior. A classical molecular simulation replaces each atom with
a **point particle** carrying only:

- mass;
- position;
- velocity.

This is an approximation. It ignores electronic structure, bond making and
breaking, and quantum effects. It is nevertheless extremely useful for
questions about molecular motion, structure, and nonbonded interactions at
ordinary temperatures.

The modelling choice is deliberate: start with the smallest object whose
motion we can describe exactly, then add chemistry later.

## 2. Position is a vector

A position tells us where a particle is relative to an origin. In three
dimensions:

$$
\mathbf r=
\begin{pmatrix}
x\\y\\z
\end{pmatrix}.
$$

In Moleculab, `position` is a NumPy array with shape `(3,)`, and every
component is in nm.

The origin is a choice. Moving the origin changes the coordinate numbers but
not the physics. What matters for interactions is usually a **displacement**:

$$
\Delta\mathbf r=\mathbf r_i-\mathbf r_j.
$$

The distance between the particles is the norm:

$$
r_{ij}=|\mathbf r_i-\mathbf r_j|.
$$

This distinction between a position and a displacement will matter repeatedly:
a bond or a force depends on the relative position, not on where the pair is in
the box.

## 3. Velocity is the rate of change of position

Velocity is

$$
\mathbf v=\frac{d\mathbf r}{dt}.
$$

In a simulation with a small time step $\Delta t$, a finite-difference
approximation is

$$
\mathbf v\approx\frac{\Delta\mathbf r}{\Delta t}.
$$

Velocity has direction. If a particle moves 0.1 nm in 1 ps to the right, its
velocity is $(0.1,0,0)$ nm/ps. If it moves to the left, the velocity is
$(-0.1,0,0)$ nm/ps.

Speed is the magnitude:

$$
|\mathbf v|=\sqrt{v_x^2+v_y^2+v_z^2}.
$$

Speed has no direction and is never negative.

![Velocity components and speed](assets/particle_velocity.svg)

*Velocity is built from components. The dashed lines show the component
vectors; the arrow is the full velocity. Here $v_x=0.08$ nm/ps,
$v_y=0.06$ nm/ps, and $|\mathbf v|=0.10$ nm/ps.*

### Average versus instantaneous velocity

If a particle moves from A to B over a finite time, the average velocity is

$$
\mathbf v_{\text{avg}}
=
\frac{\mathbf r_B-\mathbf r_A}{t_B-t_A}.
$$

Instantaneous velocity is the limit as the time interval becomes very small.
MD works with instantaneous velocities at each step.

## 4. The state of a particle

Moleculab stores the state

$$
(\mathbf r,\mathbf v;m).
$$

Acceleration is **not** stored. It is computed from the force and mass:

$$
\mathbf a=\frac{\mathbf F}{m}.
$$

Why not store acceleration as well? Because then we would have two sources of
truth. If the force changes, a stored acceleration can become inconsistent
with the current position. The state should contain only the information needed
to continue the dynamics; force and acceleration are derived quantities.

This is analogous to a game engine storing position and velocity, not “current
speed and current acceleration” separately.

## 5. Kinetic energy

Kinetic energy is the energy associated with motion:

$$
K=\frac12m|\mathbf v|^2.
$$

In three dimensions,

$$
K=\frac12m(v_x^2+v_y^2+v_z^2)
=\frac12m\mathbf v\cdot\mathbf v.
$$

### 5.1 Where does the factor 1/2 come from?

Imagine a particle starting from rest and accelerating under a constant force
$F$. Its acceleration is $a=F/m$. After time $t$, its speed is

$$
v=at.
$$

The distance travelled is

$$
d=\frac12at^2.
$$

The work done by the force is

$$
W=Fd=ma\cdot\frac12at^2=\frac12m(at)^2=\frac12mv^2.
$$

So the work done becomes exactly the kinetic energy. The factor $1/2$ is not an
arbitrary convention; it follows from the geometry of constant acceleration.

### 5.2 A second derivation: work and integration

For a force that changes with position,

$$
W=\int F\,dx.
$$

Using $F=m\,dv/dt$ and $dx=v\,dt$,

$$
W=\int m\frac{dv}{dt}v\,dt
=\int mv\,dv
=\frac12mv^2\bigg|_{\text{start}}^{\text{end}}.
$$

If the particle starts from rest, $W=\tfrac12mv^2$.

### 5.3 Kinetic energy is a scalar

Even though velocity is a vector, kinetic energy is a scalar. Direction matters
through the dot product, but the result has no direction. A particle moving
north at 100 m/s has the same kinetic energy as one moving east at 100 m/s.

### 5.4 Kinetic energy is frame-dependent

If you observe the particle from a moving train, its velocity changes, so its
kinetic energy changes. This is not a paradox. Energy is conserved within a
chosen inertial frame. Later, when we discuss temperature, we will often remove
the center-of-mass motion to avoid counting the same bulk motion twice.

## 6. Worked examples

### Example 1: the standard Moleculab particle

For

$$
m=12\ \mathrm{Da},
\qquad
\mathbf v=(0.1,0,0)\ \mathrm{nm/ps},
$$

we have

$$
K=\frac12(12)(0.1)^2=0.06\ \mathrm{kJ/mol}.
$$

Notice that no conversion constant appears. This is the payoff of the unit
system.

### Example 2: convert to SI for intuition

A speed of $0.1$ nm/ps is

$$
0.1\frac{10^{-9}\ \mathrm{m}}{10^{-12}\ \mathrm{s}}
=100\ \mathrm{m/s}.
$$

That is fast by everyday standards, but normal for a light molecule at a few
hundred kelvin.

### Example 3: doubling speed

If speed doubles from 0.1 to 0.2 nm/ps, kinetic energy increases by a factor of
four:

$$
K_{\text{new}}
=\frac12(12)(0.2)^2
=0.24\ \mathrm{kJ/mol}.
$$

This quadratic dependence is important in thermostats and in collision
analysis.

## 7. Kinetic energy curves

![Kinetic energy versus speed for several masses](assets/kinetic_energy_curves.svg)

*Kinetic energy grows as the square of speed and linearly with mass. The blue
point marks the worked example: 12 Da at 0.1 nm/ps gives 0.06 kJ/mol.*

For a fixed speed, a heavier particle carries more kinetic energy. For a fixed
kinetic energy, a heavier particle moves more slowly:

$$
v=\sqrt{\frac{2K}{m}}.
$$

This is why hydrogen atoms, being light, move much faster than oxygen or
carbon atoms at the same temperature.

## 8. Thermal motion and the Maxwell–Boltzmann distribution

At a finite temperature, particles do not all have the same speed. In a
classical gas, the probability density of speed is

$$
f(v)\propto v^2e^{-mv^2/(2RT)}.
$$

The factor $v^2$ counts the number of velocity vectors with a given magnitude
in three dimensions. The exponential is the Boltzmann factor: high speeds cost
energy.

![Maxwell–Boltzmann speed distributions](assets/maxwell_speed.svg)

*At higher temperature, the distribution shifts to higher speeds and becomes
broader. The dotted lines mark the mean speed.*

Useful intuition:

- at low temperature, most particles are slow;
- at high temperature, the distribution has a long high-speed tail;
- light particles have broader distributions than heavy particles at the same
  temperature.

This distribution is a preview of statistical mechanics. Later, temperature,
thermostats, and free energy will all depend on averages over distributions
like this one.

## 9. Kinetic energy in a many-particle system

For $N$ particles, total kinetic energy is the sum:

$$
K_{\text{total}}=\sum_{i=1}^{N}\frac12m_i|\mathbf v_i|^2.
$$

It is often useful to split motion into:

- **center-of-mass motion** of the whole system;
- **internal motion** relative to the center of mass.

The center-of-mass velocity is

$$
\mathbf V_{\text{CM}}
=
\frac{\sum_i m_i\mathbf v_i}{\sum_i m_i}.
$$

If the system is isolated, $\mathbf V_{\text{CM}}$ is constant. The internal
kinetic energy is the part that corresponds to vibrations, rotations, and
relative motion. Temperature is usually defined from the internal kinetic
energy, not from bulk translation.

This distinction will reappear in the thermostat tutorial.

## 10. Code walkthrough

The particle representation is in
[`src/moleculab/particles.py`](../../../src/moleculab/particles.py).

A simplified version:

```python
from dataclasses import dataclass, field
import numpy as np

@dataclass
class Particle:
    mass: float
    position: np.ndarray = field(default_factory=lambda: np.zeros(3))
    velocity: np.ndarray = field(default_factory=lambda: np.zeros(3))

    def __post_init__(self):
        self.mass = float(self.mass)
        self.position = np.asarray(self.position, dtype=float)
        self.velocity = np.asarray(self.velocity, dtype=float)
```

The real implementation adds validation and copy semantics. Why?

### 10.1 Positive mass

A zero or negative mass would make acceleration undefined or unphysical.
`Particle` raises an error immediately.

### 10.2 Exactly three components

A position with two components would silently broadcast against a 3D vector in
later calculations, producing a confusing result. Shape validation catches the
error at construction.

### 10.3 Copies, not references

If `Particle` stored the input list directly, then changing that list later
would change the particle without any explicit assignment. The constructor
copies the arrays so the particle owns its state.

### 10.4 Speed as a property

```python
@property
def speed(self):
    return float(np.linalg.norm(self.velocity))
```

The speed is computed from the velocity. It is not stored, for the same reason
acceleration is not stored: derived quantities should not get out of sync.

### 10.5 Kinetic energy function

```python
def kinetic_energy(particle):
    velocity_squared = float(
        np.dot(particle.velocity, particle.velocity)
    )
    return 0.5 * particle.mass * velocity_squared
```

`np.dot` is the dot product. It is equivalent to summing the squares of the
components, but it expresses the vector equation directly.

Try the code:

```bash
uv run python - <<'PY'
from moleculab import Particle, kinetic_energy

particle = Particle(
    mass=12.0,
    position=[0.0, 0.0, 0.0],
    velocity=[0.1, 0.0, 0.0],
)
print("speed:", particle.speed)
print("kinetic energy:", kinetic_energy(particle))
PY
```

Expected:

```text
speed: 0.1
kinetic energy: 0.06
```

Tests:

```bash
uv run pytest tests/test_particles.py tests/test_physics.py -q
```

## 11. Animated example: constant velocity

![Animated particle moving at constant velocity](assets/particle_motion.gif)

*A particle with constant velocity moves in a straight line. No force acts, so
the velocity arrow is constant and kinetic energy stays constant.*

This animation separates two ideas that are easy to conflate:

- position changes with time;
- velocity can stay constant while position changes.

A force is what changes velocity. In this animation there is no force, so the
velocity arrow never changes.

## 12. Numerical experiments

### Experiment 1: mass and energy

For a fixed speed of 0.1 nm/ps, compute kinetic energy for masses 1, 12, and
39.948 Da. Confirm the ratio is the same as the mass ratio.

### Experiment 2: speed and energy

For a fixed mass of 12 Da, compute kinetic energy for speeds 0.05, 0.1, 0.2,
and 0.4 nm/ps. Confirm that doubling speed quadruples energy.

### Experiment 3: thermal speed scale

Compute the mean speed from the Maxwell–Boltzmann distribution for mass 12 Da
at 300 K. Compare it with 0.1 nm/ps. Does the thermal speed match the example
particle?

### Experiment 4: reference frames

Take a particle with velocity (0.1, 0, 0) nm/ps. Add a bulk velocity
(0.1, 0, 0) to every particle in the system. How does the kinetic energy change?
Why does this matter for temperature?

## 13. Common mistakes

1. **Using speed where a vector is needed.** A force calculation needs
   direction.
2. **Adding speeds like scalars.** Use vector addition.
3. **Using mass in the wrong unit.** The shortcut only works with Da.
4. **Forgetting the square.** Kinetic energy is quadratic in speed.
5. **Treating kinetic energy as conserved by itself.** It can convert to
   potential energy.
6. **Confusing one particle's kinetic energy with temperature.** Temperature is
   a statistical property of many degrees of freedom.
7. **Ignoring the center-of-mass motion.** Bulk translation can contaminate a
   temperature estimate.

## 14. Connections to later tutorials

- **Potential energy and force:** kinetic energy will exchange with potential
  energy.
- **Velocity Verlet:** velocity is the quantity the integrator updates.
- **Temperature:** average kinetic energy per degree of freedom defines
  temperature.
- **Thermostats:** thermostats rescale velocities, not positions.
- **Free energy:** kinetic contributions are often integrated out, but they
  matter for the partition function.
- **Docking:** scoring functions usually ignore explicit velocities, but the
  same particle representation underlies the molecular model.

## 15. Self-check with answers

1. **A 1 Da particle moves at 0.2 nm/ps. What is $K$?**
   $0.5\times1\times0.2^2=0.02$ kJ/mol.

2. **A 12 Da particle has $K=0.54$ kJ/mol. What is its speed?**
   $v=\sqrt{2K/m}=\sqrt{1.08/12}=0.3$ nm/ps.

3. **Does changing position change kinetic energy?**
   No. Kinetic energy depends on mass and velocity only.

4. **Two particles have equal speeds but opposite velocities. Do they have the
   same kinetic energy?**
   Yes. Kinetic energy depends on speed, not direction.

5. **Why is acceleration not stored?**
   It is determined by force and mass; storing it would risk inconsistency.

6. **What is the difference between velocity and speed?**
   Velocity is a vector; speed is its magnitude.

7. **Why is kinetic energy frame-dependent?**
   Because velocity depends on the reference frame.

8. **What does the $v^2$ factor in the Maxwell–Boltzmann distribution mean?**
   It counts the number of velocity vectors with a given magnitude in three
   dimensions.

## 16. Exercises

1. Write a function `kinetic_energy_from_speed(mass, speed)` and compare it
   with the vector version.
2. Compute the kinetic energy of a 1 Da particle at the mean thermal speed at
   300 K.
3. Explain why a heavier particle moves more slowly at the same temperature.
4. Modify the animation idea: what would the velocity arrow look like if a
   constant force acted? Sketch it before coding.
5. Read `tests/test_particles.py` and list the invalid states it rejects.

## 17. Summary

- A classical particle is described by mass, position, and velocity.
- Position and velocity are 3D vectors; speed is the magnitude of velocity.
- Kinetic energy is $K=\tfrac12m\mathbf v\cdot\mathbf v$.
- It is a scalar, non-negative, and frame-dependent.
- The code stores only the state needed to continue the dynamics and derives
  everything else.
- Thermal motion is a distribution of speeds, not a single speed.

**Previous:** [01 — Units and scales](01_units.md).  
**Next:** [03 — Potential energy and force](03_potential_force.md).
