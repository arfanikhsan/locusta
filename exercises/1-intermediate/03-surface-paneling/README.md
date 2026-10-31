# 03 – Surface Paneling

**Level:** Intermediate · **Status:** ⬜ Not started

![result](images/result.png)

## Brief
Panelize a doubly-curved surface with quad, triangle and diamond panels, then vary the panels by a parameter (openness, depth).

## Learning goals
- Isotrim / SubSrf workflow
- Surface normals and frames
- Panel variation driven by analysis (curvature, attractor)

## Key components
- `Isotrim`
- `Divide Domain²`
- `Evaluate Surface`
- `Surface Closest Point`
- `Morph (Box Morph)`

## Steps
1. [ ] Loft a doubly-curved base surface
2. [ ] Subdivide it into a panel grid
3. [ ] Offset panels along their normals
4. [ ] Morph a component into each cell with Box Morph
5. [ ] Vary the opening size by an attractor or curvature

## Stretch goals
- [ ] Planarize panels and check deviation
- [ ] Number panels for fabrication

## Plugins
- LunchBox (optional)

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
