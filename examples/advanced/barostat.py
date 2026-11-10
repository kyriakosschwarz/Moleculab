"""Run a short NPT-like simulation with the Berendsen barostat.

The system is a low-density coarse-grained argon box.  The pressure has a
kinetic contribution from particle velocities and a virial contribution from
the interactions.  The barostat scales the periodic box and every particle
position together.

The compressibility is deliberately exaggerated so the box responds visibly
during a short demo.  Real water has kappa ~ 4.5e-5 bar^-1 and needs much
longer pressure-coupling times.

Run it with:

    uv run python examples/advanced/barostat.py
"""

import numpy as np

from moleculab import (
    ARGON,
    berendsen_barostat,
    berendsen_thermostat,
    build_solvent_box,
    kinetic_pressure,
    kinetic_temperature,
    pressure_bar,
    velocity_verlet_step_system,
    virial_pressure,
)


def main() -> None:
    system = build_solvent_box(
        ARGON,
        27,
        density=0.5,
        temperature=300.0,
        seed=5,
    )

    target_pressure = 30.0  # bar
    dt = 0.001
    steps = 1000

    print("Berendsen NPT demo")
    print(f"  particles: {system.n_particles}")
    print(f"  initial box: {system.box.lengths} nm")
    print(f"  initial volume: {system.box.volume:.4f} nm^3")
    print(f"  initial temperature: {kinetic_temperature(system):.2f} K")
    print(f"  initial pressure: {pressure_bar(system):.2f} bar")
    print(f"    kinetic part: {kinetic_pressure(system):.4f} kJ/(mol nm^3)")
    print(f"    virial part:  {virial_pressure(system):.4f} kJ/(mol nm^3)")
    print(f"  target pressure: {target_pressure:.1f} bar")
    print()

    pressures = []
    for _ in range(steps):
        system = velocity_verlet_step_system(system, dt)
        system = berendsen_thermostat(
            system,
            temperature=300.0,
            coupling_time=0.2,
            dt=dt,
        )
        system = berendsen_barostat(
            system,
            target_pressure=target_pressure,
            coupling_time=0.1,
            dt=dt,
            compressibility=1e-2,
        )
        pressures.append(pressure_bar(system))

    print("After NPT-like dynamics")
    print(f"  final box: {system.box.lengths} nm")
    print(f"  final volume: {system.box.volume:.4f} nm^3")
    print(f"  final temperature: {kinetic_temperature(system):.2f} K")
    print(f"  final pressure: {pressure_bar(system):.2f} bar")
    print(f"  mean late pressure: {np.mean(pressures[steps // 2:]):.2f} bar")
    print(f"  final potential energy: {system.potential_energy():.4f} kJ/mol")
    print()
    print("The volume changed because the barostat couples the box to the")
    print("pressure.  Particle positions were scaled with the box.")


if __name__ == "__main__":
    main()
