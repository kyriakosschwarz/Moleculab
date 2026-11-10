"""Decouple a neutral ligand from an explicit Ewald water box.

The water-water interactions use Ewald electrostatics and rigid SHAKE/RATTLE
constraints.  A neutral ligand is held fixed at the box center and decoupled
from the water through soft-core Lennard-Jones potentials.  FEP, TI, and MBAR
are computed from actual thermostatted MD samples.

There is no simple analytic reference for this solvated system, so the three
estimators are cross-validated against one another.

Run it with:

    uv run python examples/advanced/solvated_alchemy.py
"""

import numpy as np

from moleculab import (
    EwaldElectrostatics,
    SolvatedLigandAlchemy,
    build_water_box,
)


def main() -> None:
    ewald = EwaldElectrostatics(alpha=4.0, reciprocal_cutoff=10.0)
    water = build_water_box(
        4,
        box_length=2.5,
        temperature=300.0,
        seed=3,
        electrostatics=ewald,
    )
    alchemy = SolvatedLigandAlchemy(
        water,
        ligand_epsilon=1.0,
        ligand_sigma=0.3,
        temperature=300.0,
        dt=1e-4,
        equilibration_steps=50,
        sample_every=5,
        coupling_time=0.05,
    )

    print("Solvated-ligand soft-core decoupling")
    print(f"  water molecules: {alchemy.n_water // 3}")
    print(f"  water particles: {alchemy.n_water}")
    print(f"  box: {water.box.lengths} nm")
    print(f"  Ewald cutoff: {ewald.reciprocal_cutoff:.1f} nm^-1")
    print(f"  ligand position: {alchemy.ligand_position} nm")
    print(f"  ligand epsilon/sigma: {alchemy.ligand_epsilon:.2f} / "
          f"{alchemy.ligand_sigma:.2f}")
    print()

    n_windows = 5
    n_samples = 10
    fep = alchemy.fep(n_windows=n_windows, n_samples=n_samples, seed=11)
    ti = alchemy.ti(n_windows=n_windows, n_samples=n_samples, seed=12)
    mbar = alchemy.mbar(n_windows=n_windows, n_samples=n_samples, seed=13)
    values = np.array([fep.delta_g, ti.delta_g, mbar.delta_g[-1]])

    print(f"{'method':>8} {'dG (kJ/mol)':>16} {'std error':>14}")
    print(f"{'FEP':>8} {fep.delta_g:16.6f} {fep.standard_error:14.6f}")
    print(f"{'TI':>8} {ti.delta_g:16.6f} {ti.standard_error:14.6f}")
    print(f"{'MBAR':>8} {mbar.delta_g[-1]:16.6f} {'--':>14}")
    print()
    print(f"  cross-estimator spread: {np.ptp(values):.6f} kJ/mol")
    print()
    print("The water-water energy is common to all alchemical states and")
    print("cancels from the ligand decoupling free energy.  The ligand-water")
    print("soft-core path provides the only lambda-dependent energy.")


if __name__ == "__main__":
    main()
