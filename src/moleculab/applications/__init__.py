"""Applications: docking, solvent, water, free energy, acid-base, QM/MM, provenance, and structural records.

Applications answer scientific questions using the lower layers. They may
combine models, algorithms, and protocols, but they do not redefine the core
physics.  Structure manifests record experimental inputs and explicit
selections; ligand identity records compare supplied deposited and reference
graphs; ligand coordinate records attach deposited coordinates to that mapping
and record occupancy, alternate locations, and geometry; binding-site records
select the explicit protein environment and atom contacts; interaction records
classify only explicitly supplied chemical evidence; evidence fingerprints keep
present, absent, and unknown distinct; preparation records hold separate
chemical states without mutating the experimental inputs; state enumeration
keeps protonation and tautomer alternatives unselected; state-aware assessment
applies selected states to the same experimental contacts without ranking them;
protein preparation records explicit pocket-side template chemistry; combination
matrices cross ligand and protein states without ranking them. None of them
parses coordinates, chooses a structure, or invents an atom mapping on the
caller's behalf.
"""
