"""Compare a cell-list neighbour screen with brute-force pair enumeration.

The cell list does not replace the force calculation.  It only proposes
candidate pairs; the caller still computes the exact minimum-image distance and
applies the potential.

Run it with:

    uv run python examples/physics/neighbour_list.py
"""

import time

import numpy as np

from moleculab import (
    Box,
    LennardJonesPotential,
    Particle,
    System,
    build_neighbour_list,
)

BOX_LENGTH = 6.0
CUTOFF = 1.0
N_PARTICLES = 300


def brute_force_pairs(positions, cutoff, box):
    pairs = []
    for i in range(len(positions)):
        for j in range(i + 1, len(positions)):
            displacement = box.minimum_image(positions[i] - positions[j])
            if np.linalg.norm(displacement) < cutoff:
                pairs.append((i, j))
    return set(pairs)


def main() -> None:
    rng = np.random.default_rng(21)
    positions = rng.uniform(0.0, BOX_LENGTH, size=(N_PARTICLES, 3))
    box = Box([BOX_LENGTH, BOX_LENGTH, BOX_LENGTH])

    neighbours = build_neighbour_list(positions, CUTOFF, box=box)
    true_pairs = brute_force_pairs(positions, CUTOFF, box)

    candidates = set(neighbours)
    filtered = set()
    for i, j in candidates:
        displacement = box.minimum_image(positions[i] - positions[j])
        if np.linalg.norm(displacement) < CUTOFF:
            filtered.add((i, j))

    all_pairs = N_PARTICLES * (N_PARTICLES - 1) // 2
    print("Cell-list neighbour screen")
    print(f"  particles:        {N_PARTICLES}")
    print(f"  box:              {BOX_LENGTH:g} nm")
    print(f"  cutoff:           {CUTOFF:g} nm")
    print(f"  cells:            {neighbours.n_cells}")
    print(f"  cell size:        {tuple(round(v, 4) for v in neighbours.cell_size)} nm")
    print()
    print(f"  all pairs:        {all_pairs}")
    print(f"  candidate pairs:  {len(candidates)} ({len(candidates) / all_pairs:.1%})")
    print(f"  true pairs:       {len(true_pairs)}")
    print(f"  filtered pairs:   {len(filtered)}")
    print()
    print(f"  all true pairs included: {true_pairs <= candidates}")
    print(f"  filtered == brute force: {filtered == true_pairs}")
    print()
    print("Extra candidates are expected: the cell list is a screen, not a cutoff.")

    system = System(
        particles=[Particle(12.0, position) for position in positions],
        pair_potentials=[LennardJonesPotential(epsilon=1.0, sigma=0.3)],
        box=box,
    )

    start = time.perf_counter()
    forces_brute = system.forces(cutoff=CUTOFF)
    brute_time = time.perf_counter() - start

    start = time.perf_counter()
    forces_neighbour = system.forces(
        cutoff=CUTOFF,
        use_neighbour_list=True,
    )
    neighbour_time = time.perf_counter() - start

    print()
    print("System force parity")
    print(f"  max force difference: {np.max(np.abs(forces_brute - forces_neighbour)):.3e}")
    print(f"  brute-force time:     {brute_time:.4f} s")
    print(f"  neighbour-list time:  {neighbour_time:.4f} s")
    if neighbour_time > 0.0:
        print(f"  speed-up:             {brute_time / neighbour_time:.2f}x")


if __name__ == "__main__":
    main()
