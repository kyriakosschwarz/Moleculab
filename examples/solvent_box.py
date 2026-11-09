"""Build and simulate an explicit coarse-grained water box.

This is the minimal explicit-solvent path in Moleculab.  Each water molecule is
a single Lennard-Jones particle, so the model captures density, excluded
volume, and dispersion but not hydrogen bonding.  The example shows how to:

* build a periodic solvent box at a target density;
* initialize velocities at a target temperature;
* run a short NVT-like loop;
* analyze the solvent with MSD and g(r).

Run it with:

    uv run python examples/solvent_box.py
"""

import numpy as np

from moleculab import (
    COARSE_GRAINED_WATER,
    berendsen_thermostat,
    build_solvent_box,
    kinetic_temperature,
    mean_squared_displacement,
    radial_distribution_function,
    velocity_verlet_step_system,
)


def main() -> None:
    n_molecules = 64
    density = 30.0
    target_temperature = 300.0
    coupling_time = 0.2
    dt = 0.001
    steps = 500
    sample_every = 10

    system = build_solvent_box(
        COARSE_GRAINED_WATER,
        n_molecules,
        density=density,
        temperature=target_temperature,
        seed=3,
    )

    print("Explicit coarse-grained water box")
    print(f"  model: {COARSE_GRAINED_WATER.name}")
    print(f"  molecules: {n_molecules}")
    print(f"  box: {system.box.lengths} nm")
    print(f"  density: {n_molecules / system.box.volume:.2f} molecules/nm^3")
    print(f"  initial temperature: {kinetic_temperature(system):.2f} K")
    print()

    frames = []
    temperatures = []
    for step in range(steps):
        system = velocity_verlet_step_system(system, dt)
        system = berendsen_thermostat(
            system,
            temperature=target_temperature,
            coupling_time=coupling_time,
            dt=dt,
        )
        temperatures.append(kinetic_temperature(system))
        if (step + 1) % sample_every == 0:
            frames.append(system.positions.copy())

    trajectory = np.asarray(frames)
    times = np.arange(trajectory.shape[0]) * sample_every * dt

    print("NVT-like trajectory")
    print(f"  simulated time: {steps * dt:.2f} ps")
    print(f"  final temperature: {kinetic_temperature(system):.2f} K")
    print(f"  mean after warm-up: {np.mean(temperatures[steps // 2:]):.2f} K")
    print()

    msd = mean_squared_displacement(trajectory, box=system.box, times=times)
    print("Mean-squared displacement")
    print(f"{'lag (ps)':>10} {'MSD (nm^2)':>14}")
    for index in [0, 10, 20, 30, 40, 49]:
        print(f"{msd.lag[index]:10.3f} {msd.msd[index]:14.6f}")
    print()

    rdf = radial_distribution_function(
        trajectory,
        box=system.box,
        r_max=float(np.min(system.box.lengths) / 2.0),
        n_bins=40,
    )
    print("Radial distribution function")
    print(f"  first peak: r = {rdf.r[np.argmax(rdf.g)]:.3f} nm, "
          f"g = {rdf.g.max():.3f}")
    print(f"  g(r) near the box half-width: {rdf.g[-1]:.3f}")
    print()
    print("This united-atom model shows excluded-volume structure but not")
    print("the directional hydrogen-bond network of atomistic water.")


if __name__ == "__main__":
    main()
