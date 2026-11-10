"""Trajectory, contact, binding-site, interaction, fingerprint, pose, and SAR analysis.

Analysis reads position arrays and an optional box. It does not import System
and does not modify a simulation state. Binding-site geometry and interaction
detection work with prepared Cartesian and chemical-feature inputs; structure
parsing and chemical perception belong to other layers. Fingerprints and
occupancy consume interaction records and do not assume where the frames came
from. Trajectory contacts consume coordinate-derived matrices supplied by the
caller. Docking pose analysis consumes externally generated coordinates and
scores; it is not a docking search engine. SAR analysis consumes prepared assay
values, descriptors, and docking scores without predicting activity.
"""
