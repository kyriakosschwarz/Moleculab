"""Dock a rigid two-atom ligand into a simple octahedral receptor pocket.

The receptor is six typed particles arranged around the origin.  The ligand is
a rigid two-atom molecule.  Docking samples random translations and rotations
and ranks poses by the receptor-ligand Lennard-Jones plus reaction-field
Coulomb energy.

This is deliberately minimal: only rigid-body search and nonbonded scoring are
included.  There is no flexible receptor, desolvation term, or empirical
scoring function.

Run it with:

    uv run python examples/rigid_docking.py
"""

import numpy as np

from moleculab import (
    Box,
    LennardJonesParameters,
    NonbondedModel,
    Particle,
    ReactionFieldElectrostatics,
    RigidLigand,
    System,
    dock_rigid_ligand,
    score_pose,
)


def make_receptor() -> System:
    radius = 0.7
    positions = []
    for axis in range(3):
        for sign in (-1.0, 1.0):
            position = np.zeros(3)
            position[axis] = sign * radius
            positions.append(position)

    particles = [
        Particle(12.0, position, atom_type="R", charge=-0.15)
        for position in positions
    ]
    model = NonbondedModel(
        parameters={
            "R": LennardJonesParameters(epsilon=0.8, sigma=0.35),
            "L": LennardJonesParameters(epsilon=0.5, sigma=0.3),
        },
        cutoff=1.8,
        electrostatics=ReactionFieldElectrostatics(dielectric=78.3),
    )
    return System(
        particles=particles,
        nonbonded_model=model,
        box=Box([4.0, 4.0, 4.0]),
    )


def make_ligand() -> RigidLigand:
    return RigidLigand(
        [
            Particle(12.0, [-0.075, 0.0, 0.0], atom_type="L", charge=0.15),
            Particle(12.0, [0.075, 0.0, 0.0], atom_type="L", charge=0.15),
        ]
    )


def main() -> None:
    receptor = make_receptor()
    ligand = make_ligand()

    centered = score_pose(
        receptor,
        ligand,
        translation=[0.0, 0.0, 0.0],
        rotation=np.eye(3),
    )
    result = dock_rigid_ligand(
        receptor,
        ligand,
        center=[0.0, 0.0, 0.0],
        translation_radius=0.6,
        n_samples=1000,
        seed=42,
        max_results=5,
    )

    print("Rigid-ligand docking")
    print(f"  receptor particles: {receptor.n_particles}")
    print(f"  ligand atoms: {ligand.n_atoms}")
    print(f"  nonbonded cutoff: {receptor.nonbonded_model.cutoff:.2f} nm")
    print(f"  sampled poses: 1000")
    print()
    print(f"  ligand centered at origin: {centered.score:.4f} kJ/mol")
    print(f"  best docked pose:          {result.best.score:.4f} kJ/mol")
    print(f"  best COM translation:      {result.best.translation} nm")
    print(f"  improvement:               {centered.score - result.best.score:.4f} kJ/mol")
    print()
    print(f"{'rank':>4} {'score (kJ/mol)':>16} {'COM x':>9} {'COM y':>9} {'COM z':>9}")
    for rank, pose in enumerate(result.poses, start=1):
        print(
            f"{rank:4d} {pose.score:16.4f} "
            f"{pose.translation[0]:9.4f} "
            f"{pose.translation[1]:9.4f} "
            f"{pose.translation[2]:9.4f}"
        )

    print()
    print("The score is a physical interaction energy, not a calibrated binding")
    print("free energy.  A real docking workflow adds solvation, entropy, and")
    print("receptor flexibility, and often refines the top poses locally.")


if __name__ == "__main__":
    main()
