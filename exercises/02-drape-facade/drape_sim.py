"""
Drape Simulation - Grasshopper Python Script component
Works in the Rhino 8 "Python 3 Script" component and in Rhino 7 GhPython.

A cloth hangs from a top rail curve. The fabric is wider than the rail
(fullness > 1), so it gathers into vertical folds, falls under gravity,
wraps over collider objects and can pool on a floor.

INPUTS  (right-click each input -> set Access and Type hint)
    rail        Curve         item   Top edge the fabric hangs from (default: 6 m line at z=4)
    drop        float         item   Hanging length                          (default 4.0)
    fullness    float         item   Fabric width / rail length, >1 = folds  (default 1.6)
    pleats      int           item   Number of folds along the rail          (default 10)
    pleat_depth float         item   Fold depth, 0 = auto from fullness      (default 0)
    variation   float         item   0-1, randomness of fold depth           (default 0.3)
    res_u       int           item   Points across                           (default 61)
    res_v       int           item   Points down                             (default 40)
    iterations  int           item   Simulation steps                        (default 150)
    gravity     float         item   Pull per step                           (default 0.01)
    stiffness   float         item   0-1, stretch resistance                 (default 1.0)
    bend        float         item   0-1, bending resistance                 (default 0.1)
    wind        Vector3d      item   Optional constant push, e.g. (0, 0.003, 0)
    colliders   GeometryBase  list   Optional Mesh / Brep obstacles
    thickness   float         item   Gap kept between fabric and colliders   (default 0.05)
    floor_z     float         item   Optional floor height for pooling fabric
    seed        int           item   Random seed                             (default 1)

OUTPUTS
    mesh        Quad mesh of the fabric (smooth it with Weaverbird / SubD for renders)
    surface     NURBS surface through the simulated points
    points      Simulated points
    info        Summary text
"""
import math
import random

# ----------------------------------------------------------------------------
# core (pure Python, no Rhino needed)
# ----------------------------------------------------------------------------

def _offset_length(amp, nu, dl, pleats, mult):
    total = 0.0
    prev = 0.0
    for i in range(nu):
        t = float(i) / (nu - 1)
        o = amp * mult[i] * math.sin(2.0 * math.pi * pleats * t)
        if i > 0:
            total += math.sqrt(dl * dl + (o - prev) ** 2)
        prev = o
    return total


def auto_amplitude(rail_len, fullness, nu, pleats, mult):
    """Fold depth so that the folded fabric width equals rail_len * fullness."""
    target = rail_len * fullness
    dl = rail_len / (nu - 1)
    lo, hi = 0.0, rail_len
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if _offset_length(mid, nu, dl, pleats, mult) < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def build_cloth(top, side, rail_len, drop, fullness, pleats, pleat_depth,
                variation, res_v, seed):
    """top: list of (x,y,z) rail points, side: list of horizontal unit normals.
    Returns X, Y, Z, inv_mass, constraints, nu, nv."""
    rnd = random.Random(seed)
    nu, nv = len(top), res_v
    max_pleats = max(1, (nu - 1) // 4)          # at least 4 points per fold
    pleats = max(0, min(int(pleats), max_pleats))

    # smooth random variation of fold depth along the rail
    ph1, ph2 = rnd.uniform(0, 6.283), rnd.uniform(0, 6.283)
    mult = [1.0 + variation * 0.5 * (math.sin(6.283 * 1.3 * i / (nu - 1) + ph1) +
                                     math.sin(6.283 * 2.7 * i / (nu - 1) + ph2))
            for i in range(nu)]

    if pleats == 0:
        amp = 0.0
    elif pleat_depth is None or pleat_depth <= 0:
        amp = auto_amplitude(rail_len, fullness, nu, pleats, mult)
    else:
        amp = pleat_depth

    du = rail_len * fullness / (nu - 1)          # rest spacing across
    dv = drop / (nv - 1)                          # rest spacing down

    X, Y, Z, inv = [], [], [], []
    for i in range(nu):
        t = float(i) / (nu - 1)
        o = amp * mult[i] * math.sin(2.0 * math.pi * pleats * t)
        bx = top[i][0] + side[i][0] * o
        by = top[i][1] + side[i][1] * o
        bz = top[i][2] + side[i][2] * o
        for j in range(nv):
            jitter = 0.0 if j == 0 else rnd.uniform(-0.02, 0.02) * dv
            X.append(bx + side[i][0] * jitter)
            Y.append(by + side[i][1] * jitter)
            Z.append(bz - j * dv)
            inv.append(0.0 if j == 0 else 1.0)   # top row pinned

    def idx(i, j):
        return i * nv + j

    stretch, shear, bending = [], [], []
    diag = math.sqrt(du * du + dv * dv)
    for i in range(nu):
        for j in range(nv):
            if i + 1 < nu:
                stretch.append((idx(i, j), idx(i + 1, j), du))
            if j + 1 < nv:
                stretch.append((idx(i, j), idx(i, j + 1), dv))
            if i + 1 < nu and j + 1 < nv:
                shear.append((idx(i, j), idx(i + 1, j + 1), diag))
                shear.append((idx(i + 1, j), idx(i, j + 1), diag))
            if i + 2 < nu:
                bending.append((idx(i, j), idx(i + 2, j), 2 * du))
            if j + 2 < nv:
                bending.append((idx(i, j), idx(i, j + 2), 2 * dv))
    # tethers: a point can never be further from its pinned top point than
    # the fabric length between them (stops the cloth stretching like rubber)
    tethers = [(idx(i, j), idx(i, 0), j * dv) for i in range(nu) for j in range(1, nv)]
    return X, Y, Z, inv, (stretch, shear, bending, tethers), nu, nv, amp, pleats


def _project(X, Y, Z, inv, cons, k):
    sqrt = math.sqrt
    for a, b, rest in cons:
        wa, wb = inv[a], inv[b]
        ws = wa + wb
        if ws == 0.0:
            continue
        dx = X[b] - X[a]
        dy = Y[b] - Y[a]
        dz = Z[b] - Z[a]
        d = sqrt(dx * dx + dy * dy + dz * dz)
        if d < 1e-12:
            continue
        s = (d - rest) / (d * ws) * k
        X[a] += dx * s * wa; Y[a] += dy * s * wa; Z[a] += dz * s * wa
        X[b] -= dx * s * wb; Y[b] -= dy * s * wb; Z[b] -= dz * s * wb


def _tether(X, Y, Z, tethers):
    sqrt = math.sqrt
    for k, a, maxlen in tethers:
        dx = X[k] - X[a]
        dy = Y[k] - Y[a]
        dz = Z[k] - Z[a]
        d = sqrt(dx * dx + dy * dy + dz * dz)
        if d > maxlen:
            s = maxlen / d
            X[k] = X[a] + dx * s
            Y[k] = Y[a] + dy * s
            Z[k] = Z[a] + dz * s


def simulate(X, Y, Z, inv, constraints, iterations, accel, stiffness, bend,
             floor_z=None, collide=None, damping=0.98, solver_passes=4,
             collide_every=2):
    """Verlet integration + position based constraints."""
    n = len(X)
    PX, PY, PZ = X[:], Y[:], Z[:]
    ax, ay, az = accel
    stretch, shear, bending, tethers = constraints
    free = [k for k in range(n) if inv[k] > 0]
    for it in range(iterations):
        for k in free:
            x, y, z = X[k], Y[k], Z[k]
            X[k] = x + (x - PX[k]) * damping + ax
            Y[k] = y + (y - PY[k]) * damping + ay
            Z[k] = z + (z - PZ[k]) * damping + az
            PX[k], PY[k], PZ[k] = x, y, z
        for _ in range(solver_passes):
            _project(X, Y, Z, inv, stretch, stiffness)
            _project(X, Y, Z, inv, shear, stiffness * 0.5)
            if bend > 0:
                _project(X, Y, Z, inv, bending, bend)
        _tether(X, Y, Z, tethers)
        if floor_z is not None:
            for k in free:
                if Z[k] < floor_z:
                    Z[k] = floor_z
                    PX[k], PY[k], PZ[k] = X[k], Y[k], Z[k]   # friction
        if collide is not None and (it % collide_every == 0 or it == iterations - 1):
            collide(X, Y, Z, PX, PY, PZ, inv)
    return X, Y, Z


def max_strain(X, Y, Z, cons, inv):
    """Largest stretch (tension) of the simulated constraints, as a ratio."""
    worst = 0.0
    for a, b, rest in cons:
        if inv[a] + inv[b] == 0.0:
            continue
        d = math.sqrt((X[b] - X[a]) ** 2 + (Y[b] - Y[a]) ** 2 + (Z[b] - Z[a]) ** 2)
        worst = max(worst, (d - rest) / rest)
    return worst

# ----------------------------------------------------------------------------
# grasshopper
# ----------------------------------------------------------------------------
try:
    import Rhino.Geometry as rg
    IN_RHINO = True
except ImportError:
    IN_RHINO = False


def _val(name, default):
    v = globals().get(name)
    return default if v is None else v


def _collider_mesh(geos):
    joined = rg.Mesh()
    for g in geos or []:
        if g is None:
            continue
        if isinstance(g, rg.Mesh):
            joined.Append(g)
            continue
        if isinstance(g, (rg.Extrusion, rg.Surface)):
            g = g.ToBrep()
        if isinstance(g, rg.Brep):
            parts = rg.Mesh.CreateFromBrep(g, rg.MeshingParameters.Default)
            for m in parts or []:
                joined.Append(m)
    if joined.Vertices.Count == 0:
        return None
    joined.Normals.ComputeNormals()
    joined.FaceNormals.ComputeFaceNormals()
    return joined


def _make_collider(cmesh, thickness):
    bbox = cmesh.GetBoundingBox(True)
    bbox.Inflate(thickness * 2)
    search = thickness * 4

    def collide(X, Y, Z, PX, PY, PZ, inv):
        for k in range(len(X)):
            if inv[k] == 0.0:
                continue
            pt = rg.Point3d(X[k], Y[k], Z[k])
            if not bbox.Contains(pt):
                continue
            mp = cmesh.ClosestMeshPoint(pt, search)
            if mp is None:
                continue
            cp = mp.Point
            nrm = cmesh.NormalAt(mp)
            nrm.Unitize()
            dist = (X[k] - cp.X) * nrm.X + (Y[k] - cp.Y) * nrm.Y + (Z[k] - cp.Z) * nrm.Z
            if dist < thickness:
                X[k] = cp.X + nrm.X * thickness
                Y[k] = cp.Y + nrm.Y * thickness
                Z[k] = cp.Z + nrm.Z * thickness
                PX[k], PY[k], PZ[k] = X[k], Y[k], Z[k]    # friction
    return collide


def run_grasshopper():
    crv = _val("rail", None)
    if crv is None:
        crv = rg.LineCurve(rg.Point3d(0, 0, 4), rg.Point3d(6, 0, 4))
    nu = max(5, int(_val("res_u", 61)))
    nv = max(3, int(_val("res_v", 40)))

    params = crv.DivideByCount(nu - 1, True)
    top, side = [], []
    for t in params:
        p = crv.PointAt(t)
        tan = crv.TangentAt(t)
        nrm = rg.Vector3d.CrossProduct(tan, rg.Vector3d.ZAxis)
        if not nrm.Unitize():
            nrm = rg.Vector3d.YAxis
        top.append((p.X, p.Y, p.Z))
        side.append((nrm.X, nrm.Y, nrm.Z))

    X, Y, Z, inv, cons, nu, nv, amp, used_pleats = build_cloth(
        top, side, crv.GetLength(),
        float(_val("drop", 4.0)), float(_val("fullness", 1.6)),
        int(_val("pleats", 10)), float(_val("pleat_depth", 0.0)),
        float(_val("variation", 0.3)), nv, int(_val("seed", 1)))

    w = _val("wind", None)
    g = float(_val("gravity", 0.01))
    accel = (w.X if w else 0.0, w.Y if w else 0.0, (w.Z if w else 0.0) - g)

    thickness = float(_val("thickness", 0.05))
    cmesh = _collider_mesh(_val("colliders", None))
    collide = _make_collider(cmesh, thickness) if cmesh else None

    simulate(X, Y, Z, inv, cons, int(_val("iterations", 150)), accel,
             float(_val("stiffness", 1.0)), float(_val("bend", 0.1)),
             floor_z=_val("floor_z", None), collide=collide)

    pts = [rg.Point3d(X[k], Y[k], Z[k]) for k in range(len(X))]
    m = rg.Mesh()
    for p in pts:
        m.Vertices.Add(p)
    for i in range(nu - 1):
        for j in range(nv - 1):
            a = i * nv + j
            b = (i + 1) * nv + j
            m.Faces.AddFace(a, b, b + 1, a + 1)
    m.Normals.ComputeNormals()
    m.Compact()

    srf = rg.NurbsSurface.CreateFromPoints(pts, nu, nv, 3, 3)
    txt = "points {}x{} | folds {} | fold depth {:.3f} | max stretch {:.1%}".format(
        nu, nv, used_pleats, amp, max_strain(X, Y, Z, cons[0], inv))
    return m, srf, pts, txt


if IN_RHINO:
    mesh, surface, points, info = run_grasshopper()
    try:
        ghenv.Component.Message = "drape"
    except Exception:
        pass
