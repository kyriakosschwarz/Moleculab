"""Flexible-receptor docking with an explicit side-chain strain penalty.

The receptor has a fixed attractive atom and a flexible side-chain atom that
clashes with the ligand binding site.  Rigid docking cannot move the side chain
and scores the clash as a large positive energy.  Flexible docking samples
side-chain displacements and adds a harmonic strain score,

    S = E_interaction + 1/2 k |delta r|^2.

Run it with:

    uv run python examples/drug_discovery/flexible_docking.py
"""

import numpy as np

from moleculab import (
    LennardJonesParameters,
    NonbondedModel,
    Particle,
    RigidLigand,
    SideChain,
    System,
    dock_flexible_receptor,
    dock_rigid_ligand,
)


def main() -> None:
    model = NonbondedModel(
        parameters={
            "R": LennardJonesParameters(epsilon=1.0, sigma=0.3),
            "L": LennardJonesParameters(epsilon=1.0, sigma=0.3),
        },
        cutoff=1.0,
    )
    receptor = System(
        particles=[
            Particle(1.0, [0.0, 0.0, 0.0], atom_type="R"),
            Particle(1.0, [0.3, 0.0, 0.1], atom_type="R"),
        ],
        nonbonded_model=model,
    )
    ligand = RigidLigand(
        [Particle(1.0, [0.0, 0.0, 0.0], atom_type="L")]
    )
    side_chain = SideChain(atom_indices=(1,), max_displacement=0.3, stiffness=1.0)

    rigid = dock_rigid_ligand(
        receptor,
        ligand,
        center=[0.3, 0.0, 0.0],
        translation_radius=0.05,
        n_samples=200,
        seed=1,
    )
    flexible = dock_flexible_receptor(
        receptor,
        ligand,
        [side_chain],
        center=[0.3, 0.0, 0.0],
        translation_radius=0.05,
        n_ligand_samples=50,
        n_receptor_samples=100,
        seed=1,
        max_results=5,
    )

    print("Flexible-receptor docking")
    print(f"  receptor particles: {receptor.n_particles}")
    print(f"  flexible side chain atoms: {side_chain.atom_indices}")
    print(f"  side-chain max displacement: {side_chain.max_displacement:.2f} nm")
    print(f"  side-chain stiffness: {side_chain.stiffness:.2f} kJ/(mol nm^2)")
    print()
    print(f"  rigid best score:   {rigid.best.score:.4f} kJ/mol")
    print(f"  flexible best score: {flexible.best.score:.4f} kJ/mol")
    print(f"    interaction:      {flexible.best.interaction_score:.4f} kJ/mol")
    print(f"    strain:           {flexible.best.strain_score:.4f} kJ/mol")
    print(f"    side-chain move:  {flexible.best.receptor_displacements[0]}")
    print()
    print(f"{'rank':>4} {'total':>12} {'interaction':>14} {'strain':>12}")
    for rank, pose in enumerate(flexible.poses, start=1):
        print(
            f"{rank:4d} {pose.score:12.4f} "
            f"{pose.interaction_score:14.4f} {pose.strain_score:12.4f}"
        )

    print()
    print("The strain term is a simple physical calibration, not a trained")
    print("binding-affinity scoring function.  It prevents arbitrarily large")
    print("side-chain displacements from dominating the interaction energy.")


if __name__ == "__main__":
    main()
