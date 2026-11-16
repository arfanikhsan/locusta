# 04 – Grasshopper UV Texture

**Status:** ⬜ Not started

## Brief
Turn an image into real geometry with **standard Grasshopper components only**:
no script, no plugins. The image is laid onto the surface through the
surface's **UV parameter space**, and each mesh vertex moves along the surface
normal by the brightness of the pixel under it.

This is exercise 03 (Image Relief) built by hand. It's slower to set up and
has fewer options, but every step is visible on the canvas, so it's a good
way to learn how UV, normals and data matching work.

## What is UV?
Every NURBS surface is a rubber sheet stretched over a rectangle of
parameters: **u** runs along one direction, **v** along the other.
- Any point on the surface has a **(u, v)** address, like a map coordinate.
- An image also has a 2D address: **(x, y)** from 0 to 1.
- **UV mapping** = use the surface's (u, v) as the image's (x, y). Whatever
  shape the surface has, the image follows it.

**Reparameterize** (right-click a surface input → Reparameterize) makes u and v
run from 0 to 1, the same range as the image. Do this on every `S` input below.

## Files
- `definition.gh`: save your Grasshopper definition here
- `images/`: screenshots
- Test image: `../03-image-relief/water_sample.png` (seamless, so it tiles cleanly)

## The definition

```
Surface ─┬─> Mesh Surface (U, V) ─> Deconstruct Mesh ─┬─ V (vertices) ─┐
         │                                            └─ F (faces) ────┼──────────────┐
         │                                                             v              │
         ├────────────────────────────> Surface Closest Point (P = V) ─> uvP          │
         │                                                               │            │
         └─> Evaluate Surface (uv = uvP) ─> N (normal)                   │            │
                                             │                           v            │
                                             │   Deconstruct Point ─> × repeats ─> Modulus 1
                                             │                                        │
                                             │                       Construct Point <┘
                                             │                              │
                                             │                       Image Sampler ─> brightness 0–1
                                             │                              │
                                             │                       Graph Mapper (shape)
                                             │                              │
                                             └──────────> Amplitude <── × depth
                                                              │
                                    Move (V by vector) <──────┘
                                              │
                                  Construct Mesh (moved V + F) ─> relief mesh
```

### Step by step

| # | Component | Where (tab › panel) | Connect / set |
|---|-----------|---------------------|---------------|
| 1 | **Surface** param | Params › Geometry | Set one surface from Rhino. Right-click → **Reparameterize** |
| 2 | **Mesh Surface** | Mesh › Primitive | `S` = surface, `U` / `V` = sliders 50–300 (the mesh resolution) |
| 3 | **Deconstruct Mesh** | Mesh › Analysis | `M` = mesh. Gives vertices `V` and faces `F` |
| 4 | **Surface Closest Point** | Surface › Analysis | `P` = vertices, `S` = surface (reparameterized). Output `uvP` = each vertex's (u, v) |
| 5 | **Evaluate Surface** | Surface › Analysis | `S` = surface (reparameterized), `uv` = `uvP`. Output `N` = normal at each vertex |
| 6 | **Deconstruct Point** | Vector › Point | `P` = `uvP` → `X` = u, `Y` = v |
| 7 | **Multiplication** ×2 | Maths › Operators | u × `repeat_u`, v × `repeat_v` (sliders, 1 = image once) |
| 8 | **Modulus** ×2 | Maths › Operators | `A` = result, `B` = 1 → wraps back into 0–1, so the image **tiles** |
| 9 | **Construct Point** | Vector › Point | `X`, `Y` = the wrapped values, `Z` = 0 |
| 10 | **Image Sampler** | Params › Input | Double-click it → choose your image, channel **Brightness**, X and Y domain **0 to 1**. Input = the points from step 9 |
| 11 | **Graph Mapper** | Params › Input | Optional. Right-click → Graph types → Bezier/Parabola to reshape the relief (sharper crests, flatter troughs) |
| 12 | **Multiplication** | Maths › Operators | brightness × `depth` (slider, model units) |
| 13 | **Amplitude** | Vector › Vector | `V` = normal `N`, `A` = the result of step 12 |
| 14 | **Move** | Transform › Euclidean | `G` = vertices `V`, `T` = vectors from step 13 |
| 15 | **Construct Mesh** | Mesh › Primitive | `V` = moved points, `F` = faces from step 3 |

Preview: `Gradient` (Params › Input) on the brightness values →
**Mesh Colours** (Mesh › Primitive) colours the mesh so you can see the image on it.

### Data matching
Everything after step 3 is **one list in the same order as the mesh
vertices**: vertex 0, its uv, its normal, its brightness, its moved position.
That's why `Construct Mesh` can reuse the original faces: the faces point at
vertex numbers, and the numbers didn't change. Keep it that way:
- Don't **Flatten**, **Graft**, sort or cull anything between step 3 and step 15.
- If you connect more than one surface, the lists become a tree (one branch per
  surface). That's fine, as long as every component keeps the same branches.

## Tasks
1. **Flat panel.** A 1000 × 1000 planar surface, `Mesh Surface` 150 × 150,
   `depth` 10, `repeat` 1. Check that the image isn't mirrored or upside down.
   If it is, swap the Image Sampler's X or Y domain to **1 to 0**.
2. **Tiling.** Set `repeat_u` = `repeat_v` = 3. The image repeats 9 times.
   Try a photo that isn't seamless and look at the seams it makes.
3. **Real size.** Make the image always 300 units wide, whatever the
   surface size: **Dimensions** (Surface › Analysis) gives the surface's
   size in U and V. `repeat_u = U size / 300`, `repeat_v = V size / 300`.
4. **Curved surface.** Loft two or three curves into a wavy wall and
   connect it instead of the flat panel. Nothing else changes.
5. **Pillar.** A cylinder side surface (`Extrude` a circle, no caps).
   - Around the pillar, `repeat_u` must be a **whole number**, or the
     pattern won't meet itself at the seam.
   - `Mesh Surface` leaves the seam **open** (two rows of vertices in the same
     place). Add **Weld Mesh** (Mesh › Util) after `Construct Mesh` to close it.
6. **Fade the edges.** Make the relief go to 0 at the open edges, so a panel
   meets its frame cleanly. **Expression** (Maths › Script) with `x` = u and
   `y` = v from step 6:
   `min(1, min(min(x, 1-x), min(y, 1-y)) / 0.05)` → multiply it into step 12.
7. **Solid.** Bake the mesh. In Rhino use `OffsetMesh` with **Solid = Yes** to
   give it a thickness for 3D printing.

## Things you'll notice (and why)
- **The image stretches on some surfaces.** UV is the surface's *parameters*,
  not distance. Where the control points are bunched up, the parameters are
  too, and the image squeezes there. Fix: bake the surface, use Rhino's
  `Rebuild` command so its control points are evenly spaced, and set it back
  into the Surface param. Exercise 03 avoids this by measuring real distance along the
  surface instead.
- **The image is squashed on a long surface.** UV always runs 0–1, whatever the
  size. A 3000 × 1000 wall gets the image 3× wider than tall. Task 3 fixes this.
- **Jagged / noisy relief.** The mesh is coarser than the image, so each vertex
  hits a random pixel. Raise `Mesh Surface` U/V, or blur the image first
  (any image editor).
- **Slow?** `Surface Closest Point` is the heaviest step. 300 × 300 = 90,000
  vertices is fine; much more will lag on every slider move.
- **Trimmed surfaces.** UV covers the *untrimmed* surface. `Mesh Surface`
  ignores trims, so the relief covers the whole untrimmed rectangle. Trim the
  mesh afterwards (`Mesh Split`) or mesh the trimmed Brep with `Mesh Brep`
  (vertices are then uneven, but the same steps still work).

## Compare with 03 – Image Relief
| | 04 (components) | 03 (script) |
|---|---|---|
| Setup | ~15 components | 1 script component |
| Image placement | surface parameters (can stretch) | real distance (no stretch) |
| Tiling with real size | by hand (task 3) | `tile_size` |
| Pillar seam | `Weld Mesh` by hand | welded automatically |
| Solid / thickness | Rhino `OffsetMesh` | `thickness`, `cap` |
| Learning value | see every step | fast to use |

## Notes & what I learned
-
