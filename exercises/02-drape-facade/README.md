# 02 – Drape Facade

**Status:** 🟨 In progress

## Brief
Simulate draped fabric on a facade, like a satin curtain hanging behind a
sculptural frame. The fabric hangs from a top rail, gathers into vertical
folds, wraps over objects (colliders) and can pool on the floor.

## Files
- `drape_sim.py` — live cloth simulation for the Rhino 8 Python 3 script component (uses NumPy)
- `definition.gh` — save your Grasshopper definition here
- `images/` — screenshots

## Setup in Grasshopper (Rhino 8)

1. Place a **Python 3 Script** component (Maths tab → Script, or double-click the canvas and type `python`).
2. Open it and paste in all of `drape_sim.py`. The first line, `# r: numpy`,
   makes Rhino install NumPy the first time it runs (needs internet once).
3. Zoom in on the component and use the ⊕ / ⊖ buttons to create these **inputs**
   (names must match exactly). Right-click each one to set **Type hint**, and set **List Access** for `colliders`:

   | Input | Type hint | Connect | Suggested range |
   |-------|-----------|---------|-----------------|
   | `run` | bool | Boolean Toggle | True = simulate live |
   | `reset` | bool | Button | restart from the flat fabric |
   | `grab` | bool | Boolean Toggle | True = drag the cloth with the mouse |
   | `rail` | Curve | Curve (set from Rhino) | top edge of the fabric |
   | `drop` | float | Number Slider | 1 – 10 |
   | `fullness` | float | Number Slider | 1.0 – 2.5 |
   | `pleats` | int | Number Slider | 0 – 30 |
   | `pleat_depth` | float | Number Slider | 0 (= auto) – 0.5 |
   | `variation` | float | Number Slider | 0 – 1 |
   | `res_u` | int | Number Slider | 20 – 120 |
   | `res_v` | int | Number Slider | 10 – 80 |
   | `steps` | int | Number Slider | 1 – 20 (iterations per frame) |
   | `iterations` | int | Number Slider | 100 – 1000 (max after a change) |
   | `gravity` | float | Number Slider | 0.002 – 0.03 |
   | `stiffness` | float | Number Slider | 0.5 – 1.0 |
   | `bend` | float | Number Slider | 0 – 0.5 |
   | `wind` | Vector3d | Vector XYZ | e.g. (0, 0.003, 0) |
   | `colliders` | GeometryBase, **List Access** | Brep / Mesh / SubD | the frame, sculpture… |
   | `thickness` | float | Number Slider | 0.02 – 0.1 |
   | `floor_z` | float | Number Slider | ground height |
   | `seed` | int | Number Slider | 0 – 100 |

   Every input is optional. With nothing connected you get a 6 m curtain
   hanging 4 m from a rail at z = 4.

4. Create these **outputs**: `mesh`, `surface`, `points`, `info`.
5. Connect a Panel to `info` to see the state, iteration count, time per frame and stretch.

> Rhino 7: this version needs Rhino 8 (Python 3 + NumPy). The older pure-Python
> version is in the git history (commit `dc30d7c`) and works in GhPython.

## Moving the fabric with the cursor
1. Set `grab` = True and `run` = True.
2. In a Rhino viewport, **left-click on the fabric and drag**. The nearest
   point (within 15 pixels) follows the cursor, on a plane facing the camera.
   Let go and the fabric swings back. The label shows `grabbing`.
3. Clicks that miss the fabric work as normal Rhino clicks. Right-drag still
   rotates the view, so rotate first to pull the cloth in another direction.
4. Set `grab` = False when you're done, so left-clicks select objects again.

## Live mode
The component works like the Kangaroo solver: it runs `steps` iterations,
redraws, then schedules itself again, so you watch the cloth settle and can
keep dragging sliders. The component label shows `running`, `settled` or `paused`.

- **Shape sliders** (`rail`, `drop`, `fullness`, `pleats`, `pleat_depth`,
  `variation`, `res_u`, `res_v`, `seed`) rebuild the fabric and start again.
- **Physics sliders** (`gravity`, `stiffness`, `bend`, `wind`, `colliders`,
  `thickness`, `floor_z`) change the running cloth live, without a restart.
  Try dragging `wind` and watch the curtain blow.
- It stops by itself once the cloth stops moving, or after `iterations`.
  Any slider change wakes it up again.
- `run` = False pauses it. The `reset` button drops the fabric again.
- `surface` is only built once the cloth has settled (or is paused), because
  fitting a NURBS surface every frame would slow the preview down.

### Why it's lighter now
- All points are updated together with NumPy array maths instead of a Python
  loop per point. The springs are grouped into sets that never share a point,
  so each set can be solved in one array operation.
- Only points inside the collider's bounding box call Rhino's closest-point
  search, and only every second iteration.
- In a plain-Python test (61 × 40 points, sphere obstacle) this version ran
  150 iterations in about 0.7 s, vs about 4.3 s before, at about 25 ms per frame.

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

## Colliders not working? Check these
- The `info` panel shows `colliders: N (F faces)`. If it says `none` or
  `skipped`, the input type isn't supported or nothing is connected.
- Set `colliders` to **List Access** and the type hint to GeometryBase (or No Type Hint).
- Collisions use the collider's outward normals. Closed solids (Pipe with caps,
  SubD, closed mesh) work best. An open surface pushes the fabric to the
  side its normals face; use `Flip` if it's pushed the wrong way.
- The fabric is created hanging straight down from the rail. If the collider
  sits right in that plane, the fabric starts inside it and is pushed out to the
  nearest side. Move the rail slightly behind or in front of the frame.
- If fabric still passes through a thin pipe, raise `thickness` or `res_v`.

## Tips
- Start with low `res_u` and `res_v` (about 30 × 20) while you set things up.
  Increase them for the final result, because simulation time grows with point count.
- If the preview feels choppy, lower `steps`. If it settles too slowly, raise it.
- If it stops before the fabric has finished falling, raise `iterations`.
- If it passes through thin colliders, increase `thickness` or the resolution.
- Bake the mesh, or connect `surface` to the next step, once the label says `settled`.

## Alternative: Kangaroo2 version
Build the same thing with standard components to compare:
`Mesh Plane` → `Anchor` (top vertices) + `Length (Line)` on mesh edges + `Load` (gravity)
+ `Collider` / `SolidPoint` (frame) + `Angle` (bending) → `Kangaroo2 Solver` → `Show`.

## Notes & what I learned
-
