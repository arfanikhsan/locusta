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
   (names must match exactly). Right-click each one to set **Type hint**, and set
   **List Access** for `colliders` and `anchors`:

   | Input | Type hint | Connect | Suggested range |
   |-------|-----------|---------|-----------------|
   | `run` | bool | Boolean Toggle | True = simulate live |
   | `reset` | bool | Button | restart the fabric |
   | `grab` | bool | Boolean Toggle | True = drag the cloth with the mouse |
   | `anchors` | Point3d, **List Access** | Point (set from Rhino) | optional: points that hold the fabric |
   | `rail` | GeometryBase | Curve **or** planar strip Surface | Curve = thin cloth, Surface = thick cloth |
   | `start` | int | Number Slider | -1 to 3 (1 behind, -1 front, 0 hanging, 2 from the top, 3 lowered from above) |
   | `drop` | float | Number Slider | hanging length, e.g. 400 (cm) |
   | `fullness` | float | Number Slider | 1.0 – 2.5 |
   | `pleats` | int | Number Slider | 3 – 20 (fewer = bigger folds) |
   | `pleat_depth` | float | Number Slider | 0 (= auto), or a depth in model units |
   | `variation` | float | Number Slider | 0 – 1 |
   | `res_u` | int | Number Slider | 20 – 120 |
   | `res_v` | int | Number Slider | 10 – 120 |
   | `steps` | int | Number Slider | 1 – 20 (iterations per frame) |
   | `iterations` | int | Number Slider | 300 – 3000 (max after a change) |
   | `gravity` | float | Number Slider | 0.02 – 0.3 (relative to the grid spacing, unit-free) |
   | `stiffness` | float | Number Slider | 0.5 – 1.0 |
   | `bend` | float | Number Slider | 0 – 1 (default 0.1) |
   | `compression` | float | Number Slider | 0 – 1, resistance to squashing (low = more wrinkles, default 0.2) |
   | `friction` | float | Number Slider | 0 – 3, grip on colliders and floor (default 0.8) |
   | `drag` | float | Number Slider | 0 – 0.1, air resistance: low = heavy fabric (default 0.02) |
   | `wind` | Vector3d | Vector XYZ | relative like gravity, e.g. (0, 0.03, 0) |
   | `colliders` | GeometryBase, **List Access** | Brep / Mesh / SubD | the frame, sculpture… |
   | `thickness` | float | Number Slider | gap to colliders and floor, model units (e.g. 2 – 10 cm) |
   | `floor_z` | float | Number Slider | ground height |
   | `seed` | int | Number Slider | 0 – 100 |

   Every input is optional. With nothing connected you get a 600-wide curtain
   hanging 400 from a rail at z = 400 (e.g. cm).

4. Create these **outputs**: `mesh`, `thick_mesh`, `surface`, `points`, `info`.
5. Connect a Panel to `info`. It shows the state, iteration, time per frame,
   thin/thick mode, stretch, colliders and grab diagnostics.

> Rhino 7: this version needs Rhino 8 (Python 3 + NumPy). The older pure-Python
> version is in the git history (commit `dc30d7c`) and works in GhPython.

## Thin or thick fabric: curve or surface rail
- **Curve rail → thin cloth.** Works as before: the curve is the top edge.
- **Planar strip surface → thick cloth.** Draw a long, narrow planar surface
  (e.g. `Rectangle` → `Boundary Surfaces`, or `Plane Surface`) along the top of
  the facade. Its long direction is the rail (the strip's centre line) and its
  **short width is the fabric thickness**. A strip 600 × 30 gives a 30 cm thick fabric.
  - Thick fabric can't bend tighter than its thickness, so it makes **fewer,
    bigger, rounder folds** (cartoon fabric). The thick bending fades in over
    about 4 × thickness below the rail, so the fabric follows the pleats at the
    top and relaxes into round folds lower down.
  - Folds are kept at least one thickness apart.
  - The fabric keeps half its thickness (plus `thickness`) away from colliders
    and the floor, and `thick_mesh` shows it with its real thickness.
- **`thickness` keeps its old meaning:** the extra gap between the fabric and
  colliders / floor.
- For cartoon fabric (in cm): strip width 20–50, `pleats` 4–8, `bend` 0.3–0.6,
  `fullness` 1.5–2. Lower `bend` if the info panel shows a high max stretch.

## Start positions
- `start = 1` (default): the fabric starts as a flat, horizontal sheet at
  rail height, **behind** the rail, and falls and swings down under gravity.
  Behind means away from where frames usually stand, so it doesn't swing
  through them.
- `start = -1`: the same, from the other side. If it starts on the wrong
  side, use this (or flip the rail curve's direction).
- `start = 0`: starts already hanging straight down (the old behaviour).
- `start = 2`: **from the top**. A free, flat sheet at rail height, centred
  on the rail, with nothing pinned. It drops straight down onto whatever is
  underneath, like a tablecloth thrown over a table. The rail sets the
  sheet's length and position, `drop` its width. Put the rail above your
  objects; use `anchors` if some points should stay put.
- `start = 3`: **lowered from above**. The fabric hangs from the rail edge,
  raised by its own length, and the edge is lowered to the rail at a steady
  speed. Whatever is longer than the space below lands on the floor or the
  colliders and has to pile up. Make `drop` longer than the rail height (e.g.
  rail at 400, `drop` 650) to get a heap at the base.
- To slow the fall down, lower `gravity` or `steps`. The motion is speed
  limited (half a grid spacing per step) so the fabric can't fold through itself.

## Heavy, thick fabric that piles up in layers
For fabric that compresses into thick rounded folds when it lands, stacks
into layered ridges, bulges where it accumulates and spreads into soft
masses at the base (see `images/sim_heavy_pile.png`):

| Setting | Value | Why |
|---|---|---|
| `rail` | strip surface, 20–30 cm wide | gives the fabric its thickness |
| `start` | 3 | lowered from above, so the extra length piles up |
| `drop` | longer than the rail height | the extra length is what piles up |
| `drag` | 0.002 – 0.005 | heavy: keeps its momentum, slumps instead of floating |
| `compression` | 0.8 – 1.0 | it buckles into ridges instead of squashing flat |
| `bend` | 0.1 – 0.2 | low enough that it can fold into many ridges |
| `friction` | 1.5 – 2.5 | the base grips the floor, so new fabric piles on top |
| `res_u` × `res_v` | about 60 × 130 | enough points for the ridges |

How it works:
- **Thick fabric can't pass through itself**, in any direction: parts of the
  fabric that are far apart along it but closer than the thickness are pushed
  apart. Layers stack instead of flattening. (The search uses a hash grid on
  a sub-grid spaced half a thickness, and runs every 2 steps.)
- **Friction holds resting fabric** with a force set by its weight: small
  slips are undone completely (static friction), so the base doesn't slide
  out from under the falling fabric.
- **Momentum**: low `drag` keeps the falling fabric moving, so it buckles
  when it arrives instead of laying down in one neat fold.

Expect about 60–90 ms per frame for this at 60 × 130 points with thickness
(self-collision is the expensive part). Use fewer points while you set it up.

## Wrinkles and friction
Real cloth stretches very little but hardly resists being squashed: where it
lands on something, it buckles into wrinkles instead of sliding off.
- `compression` (default 0.2): how much the fabric resists squashing. Low
  values let light fabric drape and wrinkle softly. For fabric that should
  buckle into ridges when it piles up, use 0.8 – 1.0: too low, and it just
  squashes flat instead of folding.
- `friction` (default 0.8): how much it grips colliders and the floor. It
  grips with the fabric's weight, like real friction: small sliding stops
  completely (it sticks), larger sliding is slowed down. 0 = ice, 1.5+ = rubbery.
- `bend` (default 0.1): stiffer fabric makes fewer, larger wrinkles.
- Wrinkles need enough points: use at least about 80 × 80 for a drape over an object.

## Units
Any unit works. `gravity` and `wind` are relative to the grid spacing, so the
same slider values behave the same in mm, cm or m (and the solver stays stable
at high resolution). Only `drop`, `thickness`,
`pleat_depth`, `floor_z` and the rail / colliders are in model units.

## Live mode
The component works like the Kangaroo solver: it runs `steps` iterations,
redraws, then schedules itself again, so you watch the cloth settle and can
keep dragging sliders. The component label shows `running`, `settled`,
`paused` or `grabbing`.

- **Shape inputs** (`rail`, `start`, `drop`, `fullness`, `pleats`, `pleat_depth`,
  `variation`, `res_u`, `res_v`, `seed`) rebuild the fabric and start again.
- **Physics inputs** (`gravity`, `stiffness`, `bend`, `wind`, `colliders`,
  `thickness`, `floor_z`, `anchors`) change the running cloth live.
- It stops by itself once the cloth stops moving, or after `iterations`.
  Any change wakes it up again.
- `run` = False pauses it. The `reset` button restarts it.
- While it runs, `mesh` is a lighter preview (at most about 50 × 50). The
  full-resolution `mesh`, `thick_mesh`, `surface` and `points` are built once
  the cloth has settled (or is paused).

## Moving colliders without lag
- Collisions use a **distance grid** of the colliders, so checking every
  point every step is pure NumPy, with no Rhino calls.
- **Moving or rotating a collider** (e.g. with the gumball) reuses the grid
  instantly: the script recognises a rigid move and just moves the grid.
- **Changing a collider's shape** (radius, editing points…): the script
  switches to asking Rhino directly (the old method) straight away, and
  builds the new grid a little each frame in the background. The `info`
  panel shows `building collider grid NN%`. Nothing freezes.

## Moving the fabric
### With the cursor (`grab`)
1. Set `grab` = True and `run` = True.
2. In a Rhino viewport, **left-click on the fabric and drag**. The nearest
   point (within 20 pixels) follows the cursor, on a plane facing the camera.
   Let go and the fabric swings back. The label shows `grabbing`.
3. Set `grab` = False when you're done, so left-clicks select objects again.

If it doesn't react, read the `grab ON:` part of the `info` panel:
- `mouse events 0` → Rhino isn't sending mouse events to the script. Toggle
  `grab` off and on, and if that doesn't help, restart Rhino once, then try again.
- events rise but `clicks 0` → the click isn't detected as a left click.
- `clicks` rise but `hits 0` → the click missed the fabric; click closer to the mesh.
- `error …` → send me that text.

(Earlier versions failed on every left click with a `CoordinateSystem` error:
that class lives in `Rhino.DocObjects`, not `Rhino.Geometry`. Fixed.)

### With anchor points (always works)
Connect points to `anchors` (e.g. a Point parameter with points set from
Rhino). Each point grabs the fabric point nearest to it and holds it there.
Move a point in Rhino (gumball) and the fabric follows; delete it to let go.
Anchors are paired with fabric points when the number of anchors changes.

## How it works
- The fabric is a grid of points. The top row is pinned along the rail.
- The fabric is `fullness` times wider than the rail, so the extra width
  has to fold. The folds are started as a sine wave (`pleats`, `pleat_depth`),
  and `variation` makes them uneven so they look natural. The folded top edge
  is resampled so the points are evenly spaced along the fabric.
- Each step, gravity pulls the points down (Verlet integration, speed
  limited). Contacts with colliders and the floor are handled right away,
  with friction, then springs pull the points back to their rest lengths:
  - **stretch** springs: up/down and left/right neighbours. Stiff when
    stretched, soft (`compression`) when squashed, so the fabric wrinkles
  - **shear** springs: diagonal neighbours
  - **bend** springs: points 2 apart, and for thick fabric points 2 × fabric
    thickness apart. They only push apart, so they stop tight creases but
    leave big folds alone
  - **long-range limits**: points 4, 16 and 64 apart along the fabric can't be
    farther apart in space than the fabric between them, so even a free sheet
    doesn't stretch like rubber
  - **tethers** (pinned fabric): a point can't move farther from its pinned
    top point than the fabric length between them
- **Self-collision** (thick fabric): points in the same row that are close in
  space but far apart along the fabric are pushed at least one thickness apart.
- **Colliders:** points closer than the gap are pushed out along the distance
  grid's gradient, with Coulomb friction. Points below the floor gap stop there,
  so extra fabric pools on the floor. A final push-only pass at the end of
  each step makes sure nothing stays inside.
- The `info` panel shows `stretch` as the 99th percentile, and the single
  worst spring as `max` (often one spot at a sharp corner).

## Facade like the reference image
1. Draw the facade rail: a line or slight arc along the top of the ground
   floor (thin satin), or a narrow strip surface there (thick cartoon fabric).
2. Model the frame as a curve, then **Pipe** it (or **Pipe Variable** for a
   thickness that changes along the curve). Smooth it into a SubD for the chrome look.
3. Connect the frame to `colliders`. Keep the rail behind the frame.
4. Set `drop` a little longer than the facade height and set `floor_z` to
   ground level, so the hem pools.
5. For a satin look: curve rail, `fullness` 1.8–2.2, `pleats` 12–20, `bend` 0.15, `variation` 0.4.
6. Smooth the `mesh` (or `thick_mesh`) output (Weaverbird Catmull-Clark, or
   the **SubD from Mesh** component) and render with a glossy pink material.

## Colliders not working? Check these
- The `info` panel shows `colliders: N (F faces)`. If it says `none` or
  `skipped`, the input type isn't supported or nothing is connected.
- Set `colliders` to **List Access** and the type hint to GeometryBase (or No Type Hint).
- Collisions use the collider's outward normals. Closed solids (Pipe with caps,
  SubD, closed mesh) work best. An open surface pushes the fabric to the
  side its normals face; use `Flip` if it's pushed the wrong way.
- If the fabric starts inside a collider, the `info` panel shows a WARNING
  and the fabric may get pierced. Move the rail away from the frame.
- The distance grid is about 30,000 cells around the collider. A very thin
  pipe on a huge frame can fall between cells: make the pipe thicker, or
  split the colliders.

## Tips
- Start with low `res_u` and `res_v` (about 30 × 20) while you set things up.
  Increase them for the final result, because simulation time grows with point count.
- If the preview feels choppy, lower `steps`. If it settles too slowly, raise it.
- If it stops before the fabric has finished moving, raise `iterations`.
- Bake the mesh, or connect `surface` to the next step, once the label says `settled`.

## Alternative: Kangaroo2 version
Build the same thing with standard components to compare:
`Mesh Plane` → `Anchor` (top vertices) + `Length (Line)` on mesh edges + `Load` (gravity)
+ `Collider` / `SolidPoint` (frame) + `Angle` (bending) → `Kangaroo2 Solver` → `Show`.

## Notes & what I learned
-
