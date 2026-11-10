"""Applications: docking, solvent, water, free energy, acid-base, QM/MM, and provenance.

Applications answer scientific questions using the lower layers. They may
combine models, algorithms, and protocols, but they do not redefine the core
physics.  Structure manifests record experimental inputs and explicit
selections; ligand identity records compare supplied deposited and reference
graphs; ligand coordinate records attach deposited coordinates to that mapping
and record occupancy, alternate locations, and geometry; binding-site records
select the explicit protein environment and atom contacts. None of them parses
coordinates, chooses a structure, or invents an atom mapping on the caller's
behalf.
"""
