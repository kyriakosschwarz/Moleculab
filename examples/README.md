# Examples

Runnable examples are grouped by learning track and by purpose. Run them from
the repository root with the uv-managed environment.

The examples are teaching demonstrations, not a separate curriculum. Each one is
linked from the relevant tutorial or reference page.

## Physics

| Example | Demonstrates |
| --- | --- |
| [`first_particle.py`](physics/first_particle.py) | one particle in a harmonic potential |
| [`harmonic_oscillator.py`](physics/harmonic_oscillator.py) | Velocity Verlet and energy conservation |
| [`harmonic_validation.py`](physics/harmonic_validation.py) | convergence against the analytic oscillator |
| [`nonbonded_interactions.py`](physics/nonbonded_interactions.py) | Lennard-Jones and Coulomb pair potentials |
| [`bonded_molecule.py`](physics/bonded_molecule.py) | bonds and angles in a small molecule |
| [`torsion_scan.py`](physics/torsion_scan.py) | a torsional energy scan |
| [`minimization.py`](physics/minimization.py) | steepest-descent relaxation |
| [`thermostat.py`](physics/thermostat.py) | Berendsen temperature control |
| [`periodic_box.py`](physics/periodic_box.py) | minimum-image periodic boundaries |
| [`trajectory_analysis.py`](physics/trajectory_analysis.py) | MSD and RDF |
| [`neighbour_list.py`](physics/neighbour_list.py) | cell-list candidate screening |

## Chemistry

| Example | Demonstrates |
| --- | --- |
| [`chemical_structure.py`](chemistry/chemical_structure.py) | torsional conformation and conformational entropy |
| [`chemical_examples.py`](chemistry/chemical_examples.py) | chemical energy landscapes |
| [`molecular_properties.py`](chemistry/molecular_properties.py) | RDKit descriptors and Tanimoto similarity |
| [`protonation.py`](chemistry/protonation.py) | acid-base populations and tautomers |
| [`stereochemistry.py`](chemistry/stereochemistry.py) | R/S, E/Z, enantiomers, diastereomers, meso compounds |
| [`molecular_polarity.py`](chemistry/molecular_polarity.py) | formal/partial charge and dipole moments |
| [`hydrogen_bonding.py`](chemistry/hydrogen_bonding.py) | donor/acceptor features and geometric criteria |
| [`adme_tradeoffs.py`](chemistry/adme_tradeoffs.py) | descriptors, LogP/LogD, and ADME trade-offs |

## Drug discovery

| Example | Demonstrates |
| --- | --- |
| [`protein_structure.py`](drug_discovery/protein_structure.py) | structure hierarchy, uncertainty, and Å→nm conversion |
| [`binding_site.py`](drug_discovery/binding_site.py) | pocket selection and geometry |
| [`interactions.py`](drug_discovery/interactions.py) | hydrogen-bond geometry |
| [`interaction_fingerprints.py`](drug_discovery/interaction_fingerprints.py) | fingerprints and sampled-frame occupancy |
| [`rigid_docking.py`](drug_discovery/rigid_docking.py) | rigid-ligand pose search and scoring |
| [`flexible_docking.py`](drug_discovery/flexible_docking.py) | side-chain flexibility and strain |
| [`docking_pose_analysis.py`](drug_discovery/docking_pose_analysis.py) | explicit-mapping and symmetry-aware RMSD |
| [`compare_binding_modes.py`](drug_discovery/compare_binding_modes.py) | pose clustering and interaction-pattern comparison |
| [`docking_benchmark.py`](drug_discovery/docking_benchmark.py) | valid-case success rates and invalid-case reporting |
| [`trajectory_contacts.py`](drug_discovery/trajectory_contacts.py) | contact occupancy and lifetimes |
| [`matched_molecular_pairs.py`](drug_discovery/matched_molecular_pairs.py) | matched-pair SAR with separated evidence categories |
| [`structural_observations_thermodynamics.py`](drug_discovery/structural_observations_thermodynamics.py) | why contact counts do not determine free energy |
| [`capstone_imatinib.py`](drug_discovery/capstone_imatinib.py) | reproducible imatinib chemical record and provenance ledger |
| [`capstone_structure.py`](drug_discovery/capstone_structure.py) | 2HYY retrieval, immutable source file, and [provenance manifest](../data/structures/2HYY.manifest.json) ([cache policy](../data/structures/README.md)) |
| [`capstone_ligand_identity.py`](drug_discovery/capstone_ligand_identity.py) | deposited-ligand to imatinib identity reconciliation and [atom mapping](../data/structures/2HYY.ligand-identity.json) |
| [`capstone_ligand_coordinates.py`](drug_discovery/capstone_ligand_coordinates.py) | coordinate reconciliation, occupancy, B-factors, and ligand [geometry](../data/structures/2HYY.ligand-coordinates.json) |
| [`capstone_binding_site.py`](drug_discovery/capstone_binding_site.py) | protein pocket, residues, waters, and ligand–protein [contacts](../data/structures/2HYY.binding-site.json) |
| [`capstone_interactions.py`](drug_discovery/capstone_interactions.py) | chemistry-aware interaction [classification](../data/structures/2HYY.interactions.json), preserving unknown evidence |
| [`capstone_evidence_fingerprint.py`](drug_discovery/capstone_evidence_fingerprint.py) | three-state [evidence fingerprint](../data/structures/2HYY.evidence-fingerprint.json); unknown is not zero |
| [`capstone_chemical_preparation.py`](drug_discovery/capstone_chemical_preparation.py) | separate [chemical-state preparation](../data/structures/2HYY.chemical-preparation.json) with explicit assumptions |
| [`capstone_chemical_states.py`](drug_discovery/capstone_chemical_states.py) | protonation and tautomer [state enumeration](../data/structures/2HYY.chemical-states.json); no population estimate |

## Advanced

These extend the physics and alchemical material beyond the introductory path.

| Example | Demonstrates |
| --- | --- |
| [`solvent_box.py`](advanced/solvent_box.py) | coarse-grained solvent |
| [`water_box.py`](advanced/water_box.py) | three-site water |
| [`rigid_water.py`](advanced/rigid_water.py) | SHAKE/RATTLE constraints |
| [`water_ewald.py`](advanced/water_ewald.py) | Ewald electrostatics with water |
| [`ewald.py`](advanced/ewald.py) | periodic Ewald summation |
| [`barostat.py`](advanced/barostat.py) | pressure coupling and NPT |
| [`free_energy_perturbation.py`](advanced/free_energy_perturbation.py) | harmonic FEP |
| [`thermodynamic_integration.py`](advanced/thermodynamic_integration.py) | harmonic TI |
| [`mbar.py`](advanced/mbar.py) | harmonic MBAR |
| [`lj_alchemy.py`](advanced/lj_alchemy.py) | soft-core LJ FEP/TI/MBAR |
| [`mc_alchemy.py`](advanced/mc_alchemy.py) | Monte Carlo alchemy |
| [`solvated_alchemy.py`](advanced/solvated_alchemy.py) | solvated-ligand decoupling |
| [`qm_mm.py`](advanced/qm_mm.py) | one-dimensional QM/MM model |

## Performance

| Example | Demonstrates |
| --- | --- |
| [`backend_benchmark.py`](performance/backend_benchmark.py) | CPU versus optional backend timing |
| [`amd_gpu_tuning.py`](performance/amd_gpu_tuning.py) | ROCm-oriented tuning experiments |

## Conventions

- Examples print their assumptions and expected output.
- Synthetic data are labelled as teaching stand-ins.
- Measured, calculated, predicted, and hypothesised quantities are not mixed.
- Optional dependencies are invoked with `uv run --with ...`.
