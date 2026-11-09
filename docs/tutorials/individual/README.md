# Individual tutorials 1–4

These four tutorials expand the first four sections of
[`../tutorial.md`](../tutorial.md) into a slower, more visual path from units
to a working Velocity Verlet simulation.

| Tutorial | Physics | Code | Visuals |
| --- | --- | --- | --- |
| [01 — Units and scales](01_units.md) | dimensions, SI-to-simulation conversions, molar energy | `units.py` | scale chart, unit-bridge diagram |
| [02 — Particles and kinetic energy](02_particle_energy.md) | position, velocity, speed, $K=\tfrac12mv^2$ | `Particle`, `kinetic_energy` | vector diagram, KE curves, Maxwell distribution, motion GIF |
| [03 — Potential energy and force](03_potential_force.md) | energy landscapes, gradient, harmonic model | `HarmonicPotential`, `potential_energy`, `acceleration` | energy/force curves, 2D quiver, finite-difference figure |
| [04 — Velocity Verlet](04_velocity_verlet.md) | numerical integration, energy conservation, timesteps | `velocity_verlet_step`, `simulate_velocity_verlet` | position/velocity, energy, phase space, oscillator GIF, drift plots |

## Recommended pace

Each tutorial assumes the previous one. The knowledge ladder is deliberate:

1. **Units** give every later number a physical meaning.
2. **Particle state and kinetic energy** turn motion into a measurable quantity.
3. **Potential energy and force** explain where motion comes from.
4. **Velocity Verlet** combines state, energy, and force into a trajectory.

Do not skip ahead to the integrator before the energy–force relationship feels
natural. The later molecular mechanics and MD layers reuse exactly these four
ideas.

1. Read one tutorial at a time.
2. Run the short code snippet before moving on.
3. Run the listed tests and explain what each one protects against.
4. Do the self-check questions without looking at the answers.
5. Only then start the next tutorial.

## Visual assets

Static figures are SVG; animations are GIF. They are generated from the actual
library using:

```bash
uv run python docs/tutorials/individual/generate_visuals.py
```

The assets are in [`assets/`](assets/).
