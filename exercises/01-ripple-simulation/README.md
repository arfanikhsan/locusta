# 01 – Physical Ripple Simulation

**Status:** 🟨 In progress

## Brief
Simulate ripples spreading across a water surface from one or more drop points,
and turn the result into a surface in Rhino.

## Approach
1. Make a grid of points as the water surface
2. Measure each point's distance to the drop point
3. Calculate the height: `z = A · sin(k·d − ω·t) · e^(−decay·d)`
4. Move the points up by `z` and rebuild the surface or mesh
5. Animate the time `t` with a slider
6. Add more drop points and add their heights together to make interference patterns

## Parameters
| Input | Description |
|-------|-------------|
| Amplitude (A) | Wave height |
| Wavelength | Distance between crests (k = 2π / wavelength) |
| Speed (ω) | How fast the ripples travel |
| Decay | How quickly the waves fade with distance |
| Time (t) | Animation frame |

## Files
- `definition.gh` — Grasshopper definition
- `images/` — screenshots

## Notes & what I learned
-
