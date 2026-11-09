"""Explore periodic boundary conditions with a small Lennard-Jones pair.

Two particles near opposite faces of a box can be much closer through the
periodic boundary than their raw coordinates suggest.  The minimum-image
convention finds that shorter separation.

Run it with:

    uv run python examples/periodic_box.py
"""

import numpy as np

from moleculab import Box, LennardJonesPotential, Particle, System


def main() -> None:
    box = Box([1.0, 1.0, 1.0])
    print("Box:", box, "volume:", box.volume, "nm^3")
    print()

    for displacement in ([0.9, 0.0, 0.0], [-0.9, 0.0, 0.0], [0.4, 0.0, 0.0]):
        minimum = box.minimum_image(displacement)
        print(
            f"displacement {displacement} -> minimum image {minimum}, "
            f"distance {np.linalg.norm(minimum):.3f} nm"
        )

    print()
    sigma = 0.1
    potential = LennardJonesPotential(epsilon=1.0, sigma=sigma)
    system = System(
        particles=[
            Particle(1.0, [0.05, 0.0, 0.0]),
            Particle(1.0, [0.95, 0.0, 0.0]),
        ],
        pair_potentials=[potential],
        box=box,
    )

    raw_distance = abs(0.95 - 0.05)
    minimum_distance = box.distance([0.05 - 0.95, 0.0, 0.0])
    print(f"Raw separation:     {raw_distance:.3f} nm")
    print(f"Minimum-image:      {minimum_distance:.3f} nm")
    print(f"LJ energy:          {system.potential_energy():.6f} kJ/mol")
    print(f"Forces: {system.forces()}")
    print()
    print("At r = sigma the Lennard-Jones energy is zero and the force is repulsive.")


if __name__ == "__main__":
    main()
