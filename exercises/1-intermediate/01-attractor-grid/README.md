# 01 – Attractor Grid

**Level:** Intermediate · **Status:** ⬜ Not started

![result](images/result.png)

## Brief
Build a grid of geometry whose size, rotation or height reacts to one or more attractor points or curves.

## Learning goals
- Distance-based remapping (Remap Numbers, Bounds)
- Multiple attractors (closest point, weighted sums)
- Using Graph Mapper to shape the falloff

## Key components
- `Square Grid`
- `Distance`
- `Bounds`
- `Remap Numbers`
- `Graph Mapper`
- `Scale / Extrude`

## Steps
1. [ ] Create a point grid
2. [ ] Measure distance from each point to the attractor
3. [ ] Remap distances to a useful domain
4. [ ] Drive circle radius / box height with the remapped values
5. [ ] Shape the response with a Graph Mapper

## Stretch goals
- [ ] Use a curve as attractor
- [ ] Blend two attractors with weights
- [ ] Bake as a perforated panel

## Plugins
- None

## Files
- `definition.gh` — Grasshopper definition
- `model.3dm` — Rhino model (if needed)
- `images/` — canvas and viewport screenshots

## Parameters
| Input | Type | Range | Description |
|-------|------|-------|-------------|
|  |  |  |  |

## Notes & what I learned
-
