"""Trajectory, binding-site, interaction, fingerprint, and pose analysis.

Analysis reads position arrays and an optional box. It does not import System
and does not modify a simulation state. Binding-site geometry and interaction
detection work with prepared Cartesian and chemical-feature inputs; structure
parsing and chemical perception belong to other layers. Fingerprints and
occupancy consume interaction records and do not assume where the frames came
from. Docking pose analysis consumes externally generated coordinates and
scores; it is not a docking search engine.
"""
