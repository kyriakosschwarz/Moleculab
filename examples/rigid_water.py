"""Demonstrate SHAKE/RATTLE rigid-water constraints.

The example first distorts one O-H bond and applies SHAKE directly.  It then
runs a short constrained Velocity Verlet trajectory and reports the maximum
constraint violation and temperature.

Run it with:

    uv run python examples/rigid_water.py
"""

import numpy as np

from moleculab import (
    berendsen_thermostat,
    build_water_box,
    kinetic_temperature,
    shake_positions,
    velocity_verlet_step_system,
)


def maximum_constraint_error(positions, constraints):
    errors = []
    for constraint in constraints:
        distance = np.linalg.norm(
            positions[constraint.i] - positions[constraint.j]
        )
        errors.append(abs(distance - constraint.distance))
    return float(max(errors))


def main() -> None:
    system = build_water_box(
        2,
        box_length=3.0,
        temperature=300.0,
        seed=7,
    )
    masses = np.array([particle.mass for particle in system.particles])

    print("Rigid-water constraints")
    print(f"  water molecules: {system.n_particles // 3}")
    print(f"  distance constraints: {len(system.constraints)}")
    print(f"  constraint targets: {system.constraints[0].distance:.4f} nm (O-H)")
    print()

    distorted = system.positions.copy()
    distorted[1] += np.array([0.03, 0.01, 0.0])
    before = maximum_constraint_error(distorted, system.constraints)
    shaken = shake_positions(distorted, masses, system.constraints)
    after = maximum_constraint_error(shaken.positions, system.constraints)

    print("SHAKE position correction")
    print(f"  maximum O-H/error before: {before:.6f} nm")
    print(f"  maximum O-H/error after:  {after:.3e} nm")
    print(f"  SHAKE iterations: {shaken.iterations}")
    print()

    dt = 1e-4
    steps = 500
    for _ in range(steps):
        system = velocity_verlet_step_system(system, dt)
        system = berendsen_thermostat(
            system,
            temperature=300.0,
            coupling_time=0.05,
            dt=dt,
        )

    final_error = maximum_constraint_error(system.positions, system.constraints)
    print("Constrained Velocity Verlet")
    print(f"  simulated time: {steps * dt:.3f} ps")
    print(f"  final temperature: {kinetic_temperature(system):.2f} K")
    print(f"  final maximum constraint error: {final_error:.3e} nm")
    print()
    print("SHAKE corrects positions; RATTLE removes radial velocity components.")
    print("The integrator applies both automatically when System.constraints")
    print("is non-empty.")


if __name__ == "__main__":
    main()
