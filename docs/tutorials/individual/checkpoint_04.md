# Checkpoint after part 4 — one oscillator, end to end

> **Where we are.** You have finished the one-particle sequence: units,
> particle state, energy and force, and Velocity Verlet. Before adding
> interactions, prove to yourself that the whole chain works.
>
> **What you need.** Tutorials 01–04 and a Python session with Moleculab.
>
> **What this is.** A short investigation, not a new tutorial. Do the tasks
> before reading the reference solution.

## 1. The goal

Build one particle in a harmonic well, simulate it, compare it with the exact
solution, and explain every curve.

This checkpoint combines:

- units and the acceleration shortcut;
- the particle state and kinetic energy;
- the harmonic potential and its force;
- Velocity Verlet and energy exchange;
- analytic validation.

## 2. Setup

Use

$$ m=12\ \mathrm{Da}, \qquad k=12\pi^2\ \mathrm{kJ/(mol\ nm^2)}, \qquad x_0=0.1\ \mathrm{nm}, \qquad v_0=0. $$

Then

$$ \omega=\sqrt{\frac{k}{m}}=\pi\ \mathrm{rad/ps}, \qquad T=\frac{2\pi}{\omega}=2\ \mathrm{ps}. $$

The exact solution is

$$ x(t)=0.1\cos(\pi t)\ \mathrm{nm}. $$

## 3. Tasks

### Task 1: predict before running

Without running code, answer:

1. Where is the force at $t=0$?
2. Which direction is the acceleration?
3. What are $K$, $U$, and $E$ at $t=0$?
4. At $t=0.5$ ps, what are $x$, $v$, $K$, and $U$?
5. At $t=1$ ps, what are $x$, $v$, $K$, and $U$?

Write your answers down.

### Task 2: simulate

Run the oscillator for two periods with $\Delta t=0.001$ ps.

### Task 3: plot

Plot position, velocity, kinetic energy, potential energy, and total energy on
two panels. Mark the turning points and the equilibrium crossings.

### Task 4: compare with the analytic solution

At every step, compute

$$ x_{\text{exact}}(t)=0.1\cos(\pi t). $$

Find the maximum absolute position error over the two periods.

### Task 5: timestep study

Repeat with $\Delta t=0.02$, $0.005$, and $0.00125$ ps. Compute the maximum
position error for each. Halve the timestep and check the error ratio.

### Task 6: explain

Write two or three sentences for each:

1. Why does the particle accelerate toward the origin?
2. Why does kinetic energy become maximal at the origin?
3. Why does the total energy stay nearly constant?
4. Why is the numerical trajectory not exactly the analytic one?
5. Why does the error fall by about four when the timestep is halved?

## 4. Reference solution

```python
import numpy as np

from moleculab import HarmonicPotential, Particle, simulate_velocity_verlet

mass = 12.0
stiffness = 12.0 * np.pi**2
omega = np.sqrt(stiffness / mass)

particle = Particle(mass, [0.1, 0.0, 0.0], [0.0, 0.0, 0.0])
potential = HarmonicPotential(stiffness, [0.0, 0.0, 0.0])

for dt in [0.02, 0.005, 0.00125]:
    steps = int(round(4.0 / dt))
    trajectory = simulate_velocity_verlet(particle, potential, dt=dt, steps=steps)
    exact = 0.1 * np.cos(omega * trajectory.time)
    error = np.max(np.abs(trajectory.position[:, 0] - exact))
    print(f"dt = {dt:g} ps  steps/period = {round(2.0 / dt):4d}  "
          f"max error = {error:.3e} nm")
```

The first trajectory also gives the energy arrays:

```python
trajectory = simulate_velocity_verlet(
    particle, potential, dt=0.001, steps=4000
)
initial_energy = trajectory.total_energy[0]
drift = np.max(np.abs(trajectory.total_energy - initial_energy))
print("initial energy:", initial_energy)
print("maximum energy drift:", drift)
```

## 5. Expected results

For the exact benchmark:

| $\Delta t$ (ps) | steps/period | max position error (nm) |
| --- | --- | --- |
| 0.02 | 100 | $1.816\times10^{-4}$ |
| 0.005 | 400 | $1.135\times10^{-5}$ |
| 0.00125 | 1600 | $7.094\times10^{-7}$ |

At $\Delta t=0.001$ ps, the maximum energy drift over two periods is about
$2\times10^{-6}$ kJ/mol.

The ratios are approximately four:

$$ \frac{1.816\times10^{-4}}{1.135\times10^{-5}}\approx16 \quad\text{for a factor-of-four change in }\Delta t, $$

which is consistent with second order.

## 6. Answers to the predictions

| Time (ps) | $x$ (nm) | $v$ (nm/ps) | $K$ (kJ/mol) | $U$ (kJ/mol) | $E$ (kJ/mol) |
| --- | --- | --- | --- | --- | --- |
| 0 | 0.1 | 0 | 0 | 0.5922 | 0.5922 |
| 0.5 | 0 | $-0.3142$ | 0.5922 | 0 | 0.5922 |
| 1.0 | $-0.1$ | 0 | 0 | 0.5922 | 0.5922 |

The total energy is the same at all three times, up to numerical error.

## 7. Interpretation

1. **Acceleration toward the origin.** The force is
   $\mathbf F=-k(\mathbf r-\mathbf r_0)$. At $x>0$, it points toward $-x$.
2. **Kinetic energy at the origin.** The potential is minimal there, so the
   particle has converted potential energy into kinetic energy.
3. **Total energy.** Velocity Verlet is symplectic for this conservative system,
   so the energy error stays bounded instead of drifting systematically.
4. **Difference from the analytic curve.** The integrator uses a finite
   timestep and a truncated Taylor expansion. It approximates the true
   trajectory; it does not solve it exactly.
5. **Fourth-power error reduction.** The local position error is third order,
   and the accumulated error over a fixed time is second order. Halving
   $\Delta t$ therefore reduces the error by about $2^2=4$.

## 8. Are you ready for part 5?

You are ready if you can do all of the following without looking at the notes:

- [ ] State the units of $m$, $k$, $x$, $v$, $a$, $K$, and $U$.
- [ ] Compute $K$ and $U$ from a particle state.
- [ ] Explain why $\mathbf F=-\nabla U$.
- [ ] Write one Velocity Verlet step.
- [ ] Predict the turning points and equilibrium crossings.
- [ ] Explain why the total energy is bounded.
- [ ] Explain why the trajectory error scales as $\Delta t^2$.

If one item is unclear, revisit the corresponding tutorial before starting
[05 — Nonbonded interactions](05_nonbonded.md).

**Previous:** [04 — Velocity Verlet](04_velocity_verlet.md).

**Next:** [05 — Nonbonded interactions](05_nonbonded.md).
