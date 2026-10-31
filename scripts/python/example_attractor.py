"""Attractor scaling — Grasshopper Python 3 script component (Rhino 8).

Inputs:
    pts:     list[Point3d]  grid points
    attr:    Point3d        attractor point
    min_r:   float          minimum radius
    max_r:   float          maximum radius
Outputs:
    circles: list[Circle]   circles scaled by distance to the attractor
"""
import Rhino.Geometry as rg

dists = [p.DistanceTo(attr) for p in pts]
d_min, d_max = min(dists), max(dists)
span = (d_max - d_min) or 1.0

circles = []
for p, d in zip(pts, dists):
    t = (d - d_min) / span
    r = min_r + t * (max_r - min_r)
    circles.append(rg.Circle(rg.Plane(p, rg.Vector3d.ZAxis), r))
