"""Simulate a small periodic Lennard-Jones fluid and analyze its trajectory.

The example demonstrates two standard trajectory analyses:

* mean-squared displacement (MSD), which reports how far particles move after
  a lag time;
* radial distribution function g(r), which reports local structure relative to
  an ideal gas.

Run it with:

    uv run python examples/physics/trajectory_analysis.py
"""

import numpy as np

from moleculab import (
    Box,
    LennardJonesPotential,
    Particle,
    System,
    berendsen_thermostat,
    kinetic_temperature,
    mean_squared_displacement,
    radial_distribution_function,
    velocity_verlet_step_system,
)


def lattice_positions(box_length: float, n_side: int) -> np.ndarray:
    positions = []
    spacing = box_length / n_side
    for i in range(n_side):
        for j in range(n_side):
            for k in range(n_side):
                positions.append(
                    [
                        (i + 0.5) * spacing,
                        (j + 0.5) * spacing,
                        (k + 0.5) * spacing,
                    ]
                )
    return np.array(positions)


def main() -> None:
    n_side = 3
    n_particles = n_side**3
    box_length = 1.1
    box = Box([box_length, box_length, box_length])
    sigma = 0.34
    epsilon = 0.997
    mass = 39.948

    rng = np.random.default_rng(7)
    positions = lattice_positions(box_length, n_side)
    velocities = rng.normal(scale=0.05, size=(n_particles, 3))
    velocities -= np.mean(velocities, axis=0)

    system = System(
        particles=[
            Particle(mass, positions[index], velocities[index])
            for index in range(n_particles)
        ],
        pair_potentials=[LennardJonesPotential(epsilon=epsilon, sigma=sigma)],
        box=box,
    )

    dt = 0.001
    steps = 1000
    sample_every = 10
    target_temperature = 100.0
    coupling_time = 0.2

    frames = []
    for step in range(steps):
        system = velocity_verlet_step_system(system, dt)
        system = berendsen_thermostat(
            system,
            temperature=target_temperature,
            coupling_time=coupling_time,
            dt=dt,
        )
        if (step + 1) % sample_every == 0:
            frames.append(system.positions.copy())

    trajectory = np.asarray(frames)
    times = np.arange(trajectory.shape[0]) * sample_every * dt

    print("Periodic Lennard-Jones fluid")
    print(f"  particles: {n_particles}")
    print(f"  box: {box.lengths} nm")
    print(f"  simulated time: {steps * dt:.2f} ps")
    print(f"  final temperature: {kinetic_temperature(system):.2f} K")
    print()

    msd = mean_squared_displacement(trajectory, box=box, times=times)
    print("Mean-squared displacement")
    print(f"{'lag (ps)':>10} {'MSD (nm^2)':>14}")
    for index in [0, 10, 25, 50, 75, 99]:
        print(f"{msd.lag[index]:10.3f} {msd.msd[index]:14.6f}")
    print()

    rdf = radial_distribution_function(
        trajectory,
        box=box,
        r_max=box_length / 2.0,
        n_bins=40,
    )
    print("Radial distribution function")
    print(f"{'r (nm)':>10} {'g(r)':>10}")
    for index in [0, 5, 10, 15, 20, 26, 32, 39]:
        print(f"{rdf.r[index]:10.3f} {rdf.g[index]:10.4f}")
    print()
    print("g(r) -> 1 at long range for a homogeneous fluid.")
    print("The first peak near r ~ sigma is the first solvation shell.")


if __name__ == "__main__":
    main()
