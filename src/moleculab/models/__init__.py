"""Physical models: potentials, bonded terms, topology, nonbonded forces, Ewald.

A model defines what the system is and how its energy depends on coordinates.
Models must not run time integration, rescale velocities, or score poses.
"""
