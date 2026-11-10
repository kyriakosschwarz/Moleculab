# Free-energy study records

This directory holds controlled numerical experiments that are not tied to an
experimental structure.

- `harmonic-fep-study.json` is the Batch 15 study. It compares an exact
  three-dimensional harmonic partition function with forward free-energy
  perturbation for two transformations:
  - a mild stiffness change with good overlap;
  - a pure center shift with zero exact free-energy difference but poor
    overlap.
- The record also contains an exact two-state reference with energies,
  degeneracies, populations, partition function, and free-energy difference.
- The record reports configurational partition functions, exact free-energy
  differences, sample-size convergence, effective sample size, maximum weight
  fraction, replicate spread, and confidence-interval coverage.

The harmonic model is a teaching and validation model. It is **not** an
imatinib binding calculation. It contains no protein, ligand, solvent, force
field, or experimental affinity.

Regenerate the record from the repository root:

```bash
uv run python examples/advanced/free_energy_convergence.py
```

The absolute configurational free energy omits the thermal-wavelength term, so
only free-energy differences are interpreted.
