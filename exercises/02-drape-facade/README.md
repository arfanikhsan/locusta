# 02 – Drape Facade

**Status:** 🟨 In progress

## Brief
Simulate draped fabric on a facade, like a satin curtain hanging behind a
sculptural frame. The fabric hangs from a top rail, gathers into vertical
folds, wraps over objects (colliders) and can pool on the floor.

## Files
- `drape_sim.py` — cloth simulation for a Python script component
- `definition.gh` — save your Grasshopper definition here
- `images/` — screenshots

## Setup in Grasshopper (Rhino 8)

1. Place a **Python 3 Script** component (Maths tab → Script, or double-click the canvas and type `python`).
2. Open it and paste in all of `drape_sim.py`.
3. Zoom in on the component and use the ⊕ / ⊖ buttons to create these **inputs**
   (names must match exactly). Right-click each one to set **Type hint**, and set **List Access** for `colliders`:

   | Input | Type hint | Connect | Suggested range |
   |-------|-----------|---------|-----------------|
   | `rail` | Curve | Curve (set from Rhino) | top edge of the fabric |
   | `drop` | float | Number Slider | 1 – 10 |
   | `fullness` | float | Number Slider | 1.0 – 2.5 |
   | `pleats` | int | Number Slider | 0 – 30 |
   | `pleat_depth` | float | Number Slider | 0 (= auto) – 0.5 |
   | `variation` | float | Number Slider | 0 – 1 |
   | `res_u` | int | Number Slider | 20 – 120 |
   | `res_v` | int | Number Slider | 10 – 80 |
   | `iterations` | int | Number Slider | 50 – 500 |
   | `gravity` | float | Number Slider | 0.002 – 0.03 |
   | `stiffness` | float | Number Slider | 0.5 – 1.0 |
   | `bend` | float | Number Slider | 0 – 0.5 |
   | `wind` | Vector3d | Vector XYZ | e.g. (0, 0.003, 0) |
   | `colliders` | GeometryBase, **List Access** | Brep / Mesh | the frame, sculpture… |
   | `thickness` | float | Number Slider | 0.02 – 0.1 |
   | `floor_z` | float | Number Slider | ground height |
   | `seed` | int | Number Slider | 0 – 100 |

   Every input is optional. With nothing connected you get a 6 m curtain
   hanging 4 m from a rail at z = 4.

4. Create these **outputs**: `mesh`, `surface`, `points`, `info`.
5. Connect a Panel to `info` to read the fold count and how much the fabric stretched.

> Rhino 7: use the **GhPython** component instead. The same script works.

## How it works
- The fabric is a grid of points. The top row is pinned along the rail.
- The fabric is `fullness` times wider than the rail, so the extra width
  has to fold. The folds are started as a sine wave (`pleats`, `pleat_depth`),
  and `variation` makes them uneven so they look natural.
- Each step, gravity pulls the points down (Verlet integration), then springs
  between neighbouring points pull them back to their rest lengths:
  - **stretch** springs: up/down and left/right neighbours
  - **shear** springs: diagonal neighbours
  - **bend** springs: points two apart, which make the fabric stiffer (higher `bend` = larger, softer folds)
  - **tethers**: a point can't move farther from its pinned top point than
    the fabric length between them, so the cloth doesn't stretch like rubber
- Points inside a collider are pushed out to the `thickness` distance. Points
  below `floor_z` stop there, so extra fabric pools on the floor.

## Facade like the reference image
1. Draw the facade rail: a line, or a slight arc, along the top of the ground floor.
2. Model the frame as a curve, then **Pipe** it (or **Pipe Variable** for a
   thickness that changes along the curve). Smooth it into a SubD for the chrome look.
3. Connect the frame Brep to `colliders`. Move the rail a little behind the
   frame so the fabric drapes over it.
4. Set `drop` a little longer than the facade height and set `floor_z` to
   ground level, so the hem pools.
5. For a satin look: `fullness` 1.8–2.2, `pleats` 12–20, `bend` 0.15, `variation` 0.4.
6. Smooth the `mesh` output (Weaverbird Catmull-Clark, or the **SubD from Mesh**
   component) and render with a glossy pink material.

## Tips
- Start with low `res_u` and `res_v` (about 30 × 20) while you set things up.
  Increase them for the final result, because simulation time grows with point count.
- If the fabric still looks like it's falling, raise `iterations`.
- If it passes through thin colliders, increase `thickness` or the resolution.
- Use Data Dam or a Boolean Toggle (via Stream Gate) so the simulation
  doesn't re-run every time you move a slider.

## Alternative: Kangaroo2 version
Build the same thing with standard components to compare:
`Mesh Plane` → `Anchor` (top vertices) + `Length (Line)` on mesh edges + `Load` (gravity)
+ `Collider` / `SolidPoint` (frame) + `Angle` (bending) → `Kangaroo2 Solver` → `Show`.

## Notes & what I learned
-
