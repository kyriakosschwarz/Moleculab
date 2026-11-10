"""Trajectory, binding-site, and interaction analysis.

Analysis reads position arrays and an optional box. It does not import System
and does not modify a simulation state. Binding-site geometry and interaction
detection work with prepared Cartesian and chemical-feature inputs; structure
parsing and chemical perception belong to other layers.
"""
