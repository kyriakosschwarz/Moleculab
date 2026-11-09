# 04 — Velocity Verlet and the conservation of energy

> **Where we are.** Tutorial 03 gave us a potential landscape and a force. Now
> we let the particle move through time and watch kinetic and potential energy
> trade places.
>
> **What you should already know.** Kinetic energy, potential energy, force as
> the negative gradient, and the unit system.
>
> **What you will be able to do afterwards.** Explain why a numerical
> integrator is needed, derive Velocity Verlet, implement one step, analyze a
> trajectory, choose a timestep, and diagnose energy drift.

## 1. From force to motion

Newton's second law is

$$
\mathbf F=m\mathbf a
=m\frac{d^2\mathbf r}{dt^2}.
$$

For a potential,

$$
m\frac{d^2\mathbf r}{dt^2}=-\nabla U.
$$

This is a second-order differential equation. Given an initial position and
velocity, the equation determines the future trajectory.

For a real molecular system, the right-hand side is a complicated function of
many coordinates. There is no useful closed-form solution. We must integrate
numerically.

The computational loop is:

1. Start with positions and velocities.
2. Compute forces from the potential.
3. Update positions and velocities by a small time step.
4. Repeat.

The quality of step 3 determines whether the simulation is physically
trustworthy.

## 2. The simplest attempt: explicit Euler

The most direct update is

$$
\mathbf r_{n+1}
=
\mathbf r_n+\mathbf v_n\Delta t,
$$

$$
\mathbf v_{n+1}
=
\mathbf v_n+\mathbf a_n\Delta t.
$$

This is explicit Euler. It uses the current force for the entire step.

It is easy to implement and easy to understand. It is also a poor choice for
oscillatory systems. The reason is that it treats the acceleration as constant
over the step, while in reality the force changes as the particle moves.

For a harmonic oscillator, explicit Euler systematically gains energy. The
particle spirals outward instead of staying on a closed orbit.

![Euler versus Velocity Verlet](assets/euler_vs_verlet.svg)

*For the same timestep, explicit Euler gains energy while Velocity Verlet keeps
the total energy bounded. A method can look stable for a short time and still
be physically wrong.*

## 3. A better update: Taylor expansion

Expand position in a Taylor series:

$$
\mathbf r(t+\Delta t)
=
\mathbf r(t)
+\mathbf v(t)\Delta t
+\frac12\mathbf a(t)\Delta t^2
+\frac16\mathbf j(t)\Delta t^3
+\cdots.
$$

Truncating after the quadratic term gives

$$
\mathbf r_{n+1}
=
\mathbf r_n
+\mathbf v_n\Delta t
+\frac12\mathbf a_n\Delta t^2.
$$

The position update is now second-order accurate. The neglected term is
proportional to $\Delta t^3$, so the local error is smaller than Euler's.

For velocity, we can use the average acceleration over the step:

$$
\mathbf v_{n+1}
=
\mathbf v_n
+\frac12(\mathbf a_n+\mathbf a_{n+1})\Delta t.
$$

This is the **Velocity Verlet** algorithm.

## 4. The algorithm in detail

One step:

1. Compute the force at the current position.
2. Compute the current acceleration:

$$
\mathbf a_n=\frac{\mathbf F_n}{m}.
$$

3. Update the position:

$$
\mathbf r_{n+1}
=
\mathbf r_n+\mathbf v_n\Delta t
+\frac12\mathbf a_n\Delta t^2.
$$

4. Compute the force at the **new** position.
5. Compute the new acceleration $\mathbf a_{n+1}$.
6. Update the velocity:

$$
\mathbf v_{n+1}
=
\mathbf v_n
+\frac12(\mathbf a_n+\mathbf a_{n+1})\Delta t.
$$

The extra force evaluation at the new position is the key difference from
Euler. It is also why Velocity Verlet is stable for oscillatory systems.

## 5. Why Velocity Verlet is special

Velocity Verlet has several properties that matter for MD:

- It is **second-order accurate** in the timestep.
- It is **time-reversible** in exact arithmetic.
- It is **symplectic**: it preserves the geometric structure of Hamiltonian
  mechanics.
- It has excellent long-term energy behavior for conservative systems.
- It does not require storing acceleration between steps.

“Symplectic” is a technical term, but the intuition is simpler: the method does
not systematically pump energy into or out of an oscillator. The numerical
trajectory stays near the true energy surface instead of spiraling away.

This is why Velocity Verlet, or closely related Verlet-family methods, are
standard in molecular dynamics.

## 6. The harmonic oscillator as an exact benchmark

For a one-dimensional harmonic oscillator,

$$
m\frac{d^2x}{dt^2}=-kx.
$$

The exact solution is

$$
x(t)=A\cos(\omega t+\phi),
\qquad
\omega=\sqrt{\frac{k}{m}}.
$$

The period is

$$
T=\frac{2\pi}{\omega}
=
2\pi\sqrt{\frac{m}{k}}.
$$

We choose

$$
m=12\ \mathrm{Da},
\qquad
k=12\pi^2\ \mathrm{kJ/(mol\,nm^2)}.
$$

Then

$$
\omega=\sqrt{\frac{12\pi^2}{12}}=\pi\ \mathrm{rad/ps},
\qquad
T=2\ \mathrm{ps}.
$$

Starting at $x=0.1$ nm with zero velocity:

$$
x(t)=0.1\cos(\pi t)\ \mathrm{nm}.
$$

After 2 ps, the particle should be back at $x=0.1$ nm. This gives us an exact
test for both the period and the trajectory.

## 7. Position and velocity

![Oscillator position and velocity](assets/oscillator_position_velocity.svg)

*Position is a cosine; velocity is a sine. At the turning points, velocity is
zero. At equilibrium, speed is maximal.*

The quarter-period phase difference is important:

- maximum displacement → maximum potential energy, zero kinetic energy;
- equilibrium → zero potential energy, maximum kinetic energy;
- the particle does not stop at equilibrium; inertia carries it through.

## 8. Energy exchange

![Kinetic, potential, and total energy](assets/oscillator_energy.svg)

*Kinetic and potential energy oscillate out of phase. Their sum remains nearly
constant.*

For the exact oscillator,

$$
K=\frac12mv^2,
\qquad
U=\frac12kx^2,
\qquad
E=K+U=\text{constant}.
$$

The simulation does not enforce constant energy explicitly. Velocity Verlet
approximately conserves it because of its geometric properties. Small bounded
oscillations in the total energy are normal. A systematic drift is a warning
sign.

## 9. Phase space

![Phase-space orbit of the oscillator](assets/oscillator_phase.svg)

*Plotting velocity against position gives a closed ellipse. The shape encodes
the conserved energy and the oscillation frequency.*

Phase space is a useful diagnostic:

- a closed orbit indicates stable conservative motion;
- an outward spiral means energy gain;
- an inward spiral means numerical damping;
- a sudden jump means a bad step or a force discontinuity.

In many dimensions, phase space is high-dimensional, but projections onto
important coordinates can still reveal instabilities.

## 10. Animated oscillator

![Animated mass on a spring](assets/oscillator.gif)

*A mass on a spring. The force points toward equilibrium; position, velocity,
and energy evolve together.*

Watch the feedback loop:

1. Position determines the spring force.
2. Force determines acceleration.
3. Acceleration changes velocity.
4. Velocity changes position.

The animation is not just decorative. It is a visual representation of the
integration algorithm.

## 11. Code walkthrough

The integrator is in
[`src/moleculab/integrators.py`](../../../src/moleculab/integrators.py):

```python
def velocity_verlet_step(particle, potential, dt):
    force_now = potential.force(particle.position)
    acceleration_now = acceleration(particle, force_now)

    new_position = (
        particle.position
        + particle.velocity * dt
        + 0.5 * acceleration_now * dt**2
    )

    force_next = potential.force(new_position)
    acceleration_next = acceleration(particle, force_next)

    new_velocity = (
        particle.velocity
        + 0.5 * (acceleration_now + acceleration_next) * dt
    )

    return Particle(
        mass=particle.mass,
        position=new_position,
        velocity=new_velocity,
    )
```

The trajectory loop is:

```python
def simulate_velocity_verlet(particle, potential, dt, steps):
    current = particle.copy()
    positions = np.empty((steps + 1, 3))
    velocities = np.empty((steps + 1, 3))
    energies = np.empty(steps + 1)

    for index in range(steps + 1):
        positions[index] = current.position
        velocities[index] = current.velocity
        energies[index] = total_energy(current, potential)
        if index < steps:
            current = velocity_verlet_step(current, potential, dt)

    return Trajectory(...)
```

Three design choices are worth noticing:

1. **The input particle is not mutated.** The function returns a new particle.
   This makes reasoning and testing easier.
2. **The force is evaluated at the new position.** That is the defining feature
   of Velocity Verlet.
3. **The trajectory stores every state.** You can analyze the step where
   something went wrong instead of seeing only the final result.

Run the example:

```bash
uv run python examples/harmonic_oscillator.py
```

Tests:

```bash
uv run pytest tests/test_integrators.py -q
```

## 12. Timestep and stability

![The same oscillator at several timesteps](assets/dt_stability.svg)

*Small timesteps track the exact solution closely. Large timesteps distort the
period and amplitude; very large timesteps can destabilize the trajectory.*

The safe timestep depends on the fastest motion in the system. For a harmonic
oscillator, the period is

$$
T=2\pi\sqrt{\frac{m}{k}}.
$$

A useful rule of thumb is to resolve the fastest period with many steps, not
just a few. In atomistic MD, the fastest motions are often bond vibrations
involving hydrogen, which is why typical timesteps are around 1–2 fs
($0.001$–$0.002$ ps).

## 13. Energy drift versus timestep

![Energy drift versus timestep](assets/energy_drift.svg)

*The maximum relative energy drift decreases roughly as $\Delta t^2$ for
Velocity Verlet. Smaller timesteps improve conservation but cost more CPU time.*

This is the central practical tradeoff of MD:

- larger timestep → fewer force evaluations → faster simulation;
- larger timestep → larger integration error → distorted physics;
- too large a timestep → instability or a blown-up trajectory.

A production MD code often uses constraints (SHAKE/RATTLE) to remove the
fastest bond vibrations and allow a larger timestep. Moleculab introduces rigid
water constraints later for the same reason.

## 14. Why Euler is not enough

Explicit Euler is first-order and does not preserve the symplectic structure.
For an oscillator, it tends to gain energy. A simple intuitive picture:

- the particle moves in a straight line for the whole step;
- then the force is applied as if it had been constant;
- the force is strongest at the turning points, so the correction is applied
  too late.

Velocity Verlet uses the force at both ends of the step. The averaged
acceleration is a much better approximation to the true motion.

## 15. What the trajectory contains

A `Trajectory` stores:

- time;
- position at every step;
- velocity at every step;
- kinetic energy;
- potential energy;
- total energy.

This is enough for:

- plotting motion;
- checking energy conservation;
- computing period and frequency;
- computing phase-space orbits;
- diagnosing instabilities;
- later, analyzing diffusion or structural fluctuations.

## 16. Connections to molecular dynamics and CADD

This small oscillator is the core of MD:

- a force field supplies the potential;
- the integrator advances positions and velocities;
- a thermostat controls temperature;
- a barostat controls pressure;
- analysis extracts structural and thermodynamic information.

In a protein-ligand simulation, the same loop runs for millions of steps with
thousands of atoms. The physics is more complex, but the numerical principle is
identical: compute forces, update positions, update velocities, repeat.

Understanding this tiny system makes it much easier to interpret:

- why timesteps are measured in femtoseconds;
- why energy conservation is a diagnostic;
- why constraints are used for water and hydrogen bonds;
- why a trajectory can look “stable” but still be numerically wrong.

## 17. Common mistakes

1. **Using the old acceleration for the entire velocity update.** Velocity
   Verlet averages old and new acceleration.
2. **Forgetting to recompute force at the new position.** This reduces the
   method toward Euler.
3. **Choosing a timestep that is too large.** The trajectory may look plausible
   for a few steps and then explode.
4. **Judging a simulation only by its final snapshot.** Check energy drift and
   conservation over time.
5. **Confusing period with timestep.** The period is 2 ps here; the timestep is
   0.001 ps.
6. **Expecting exact energy conservation.** Floating-point roundoff and finite
   timesteps produce small bounded oscillations.
7. **Using the wrong mass or stiffness units.** The period formula assumes
   consistent units.

## 18. Self-check with answers

1. **Why does Velocity Verlet evaluate the force twice per step?**
   The new force gives the new acceleration used in the averaged velocity
   update.

2. **For $m=12$ and $k=12\pi^2$, what is the period?**
   $T=2\pi\sqrt{12/(12\pi^2)}=2$ ps.

3. **What happens to the energy drift if the timestep is halved?**
   The error generally drops by about a factor of four for second-order
   methods.

4. **Why does an Euler oscillator gain energy?**
   It uses only the old acceleration, so it consistently overshoots.

5. **At the turning point, what are $K$ and $U$?**
   Kinetic energy is zero; potential energy is maximal.

6. **Why is 1–2 fs a typical MD timestep?**
   It resolves fast bond vibrations while remaining computationally practical.

7. **What does a closed phase-space orbit indicate?**
   Stable, approximately conservative motion.

8. **Why is the trajectory stored at every step?**
   To support analysis and debugging of the integration.

## 19. Exercises

1. Derive the period of a harmonic oscillator from $m d²x/dt²=-kx$.
2. Modify the example to start with zero position and nonzero velocity. How do
   the cosine and sine roles swap?
3. Compute the maximum relative energy drift for timesteps 0.01, 0.001, and
   0.0001 ps. Plot the result on a log-log axis.
4. Implement explicit Euler for the same oscillator and compare its energy
   drift with Velocity Verlet.
5. Explain why a smaller timestep is not always better in a real simulation
   (think about total wall-clock time and accumulated roundoff).

## 20. Summary

- Velocity Verlet updates position with old velocity and old acceleration, then
  updates velocity with the average of old and new acceleration.
- It is second-order, time-reversible, and symplectic.
- It keeps the total energy of a conservative oscillator bounded.
- The timestep must resolve the fastest motion in the system.
- Energy drift is a diagnostic, not a decoration.
- This small integrator is the engine that will later drive many-particle
  molecular dynamics.

**Previous:** [03 — Potential energy and force](03_potential_force.md).  
**Return to the tutorial index:** [Moleculab tutorial](../tutorial.md).
