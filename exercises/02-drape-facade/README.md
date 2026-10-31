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
   | `drop` | float | Number Slider | hanging length, e.g. 400 (cm) |
   | `fullness` | float | Number Slider | 1.0 – 2.5 |
   | `pleats` | int | Number Slider | 3 – 20 (fewer = bigger folds) |
   | `pleat_depth` | float | Number Slider | 0 (= auto), or a depth in model units |
   | `variation` | float | Number Slider | 0 – 1 |
   | `res_u` | int | Number Slider | 20 – 120 |
   | `res_v` | int | Number Slider | 10 – 80 |
   | `steps` | int | Number Slider | 1 – 20 (iterations per frame) |
   | `iterations` | int | Number Slider | 100 – 1000 (max after a change) |
   | `gravity` | float | Number Slider | 0.02 – 0.3 (relative, unit-free) |
   | `stiffness` | float | Number Slider | 0.5 – 1.0 |
   | `bend` | float | Number Slider | 0 – 1 |
   | `wind` | Vector3d | Vector XYZ | relative like gravity, e.g. (0, 0.03, 0) |
   | `colliders` | GeometryBase, **List Access** | Brep / Mesh / SubD | the frame, sculpture… |
   | `thickness` | float | Number Slider | fabric thickness in model units, e.g. 5 – 50 (cm) |
   | `floor_z` | float | Number Slider | ground height |
   | `seed` | int | Number Slider | 0 – 100 |

   Every input is optional. With nothing connected you get a 600-wide curtain
   hanging 400 from a rail at z = 400 (e.g. cm).

4. Create these **outputs**: `mesh`, `thick_mesh`, `surface`, `points`, `info`.
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

## Units and thickness
- **Any unit works.** `gravity` and `wind` are relative to `drop`, so the
  same slider values behave the same in mm, cm or m. Only `drop`, `thickness`,
  `pleat_depth` and `floor_z` are in model units.
- **`thickness` is the fabric's own thickness.** Thick fabric can't bend
  tighter than its thickness, so:
  - bending works over a span of 2 × `thickness`: **thicker = fewer, bigger,
    rounder folds** (cartoon fabric), whatever the resolution
  - neighbouring folds are kept at least `thickness` apart
  - the `mesh` is the fabric's middle surface, kept `thickness / 2` away from
    colliders and the floor, so the thick fabric just touches them
  - `thick_mesh` shows the fabric with its real thickness
- For cartoon fabric (in cm): `thickness` 20–50, `pleats` 4–8, `bend` 0.3–0.6,
  `fullness` 1.5–2. Lower `bend` if the info panel shows a high max stretch.

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
- While it runs, `mesh` is a lighter preview (at most about 50 × 50).
  The full-resolution `mesh`, `thick_mesh`, `surface` and `points` are built
  once the cloth has settled (or is paused).

### Why it's lighter now
- All points are updated together with NumPy array maths. Points are stored as
  `P[xyz, i, j]`, so vector lengths are three plain array additions.
- Springs are grouped into sets that never share a point, so each set is
  solved in one array operation.
- **Colliders become a distance grid.** When the colliders change, Rhino
  measures the distance to them once at about 40,000 grid points (the `info`
  panel shows how long that took). After that, every collision check is
  pure NumPy, with no Rhino calls, even at 100 × 100 points. Moving sliders
  doesn't rebuild the grid.
- Test (plain Python, 100 × 100 points, arch collider, cm units): about 55–65 ms
  per frame of 5 iterations.

## How it works
- The fabric is a grid of points. The top row is pinned along the rail.
- The fabric is `fullness` times wider than the rail, so the extra width
  has to fold. The folds are started as a sine wave (`pleats`, `pleat_depth`),
  and `variation` makes them uneven so they look natural. The folded top edge
  is resampled so the points are evenly spaced along the fabric.
- Each step, gravity pulls the points down (Verlet integration), then springs
  pull them back to their rest lengths:
  - **stretch** springs: up/down and left/right neighbours
  - **shear** springs: diagonal neighbours
  - **bend** springs: points 2 apart and points 2 × `thickness` apart. They only
    push apart, so they stop tight creases but leave big folds alone
  - **tethers**: a point can't move farther from its pinned top point than
    the fabric length between them, so the cloth doesn't stretch like rubber
- **Self-collision:** points in the same row that are close in space but far
  apart along the fabric are pushed to at least `thickness` apart.
- **Colliders:** points closer than `thickness / 2` are pushed out along the
  distance grid's gradient, with some friction. Points below
  `floor_z + thickness / 2` stop there, so extra fabric pools on the floor.

## Facade like the reference image
1. Draw the facade rail: a line, or a slight arc, along the top of the ground floor.
2. Model the frame as a curve, then **Pipe** it (or **Pipe Variable** for a
   thickness that changes along the curve). Smooth it into a SubD for the chrome look.
3. Connect the frame Brep to `colliders`. Move the rail a little behind the
   frame so the fabric drapes over it.
4. Set `drop` a little longer than the facade height and set `floor_z` to
   ground level, so the hem pools.
5. For a satin look: `fullness` 1.8–2.2, `pleats` 12–20, `bend` 0.15, `variation` 0.4,
   thin `thickness`. For a puffy cartoon look, see "Units and thickness" above.
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
  cuts through that plane, the fabric starts inside it and gets pierced
  (the `info` panel shows a WARNING). Move the rail behind the frame, at least
  `thickness / 2` plus the pipe radius away.
- The distance grid is about 40,000 cells around the collider. A very thin
  pipe on a huge frame can fall between cells: make the pipe thicker, or
  split the colliders.

## Tips
- Start with low `res_u` and `res_v` (about 30 × 20) while you set things up.
  Increase them for the final result, because simulation time grows with point count.
- If the preview feels choppy, lower `steps`. If it settles too slowly, raise it.
- If it stops before the fabric has finished falling, raise `iterations`.
- Bake the mesh, or connect `surface` to the next step, once the label says `settled`.

## Alternative: Kangaroo2 version
Build the same thing with standard components to compare:
`Mesh Plane` → `Anchor` (top vertices) + `Length (Line)` on mesh edges + `Load` (gravity)
+ `Collider` / `SolidPoint` (frame) + `Angle` (bending) → `Kangaroo2 Solver` → `Show`.

## Notes & what I learned
-
