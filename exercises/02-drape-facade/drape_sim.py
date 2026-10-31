# r: numpy
"""
Drape Simulation (live) - Rhino 8 "Python 3 Script" component

A thick cloth hangs from a top rail curve. The fabric is wider than the rail
(fullness > 1), so it gathers into vertical folds, falls under gravity,
wraps over collider objects and can pool on a floor.

UNITS: everything scales with your model, so cm, m or mm all work. Gravity
and wind are relative to the drop length; `thickness` is in model units.

THICKNESS: the fabric behaves as if it is `thickness` thick. Thick fabric
can't bend tighter than its thickness, so folds get bigger and rounder
(cartoon fabric), neighbouring folds stay `thickness` apart, and the
fabric's middle surface keeps thickness/2 away from colliders and the floor.
The `thick_mesh` output shows the fabric with that thickness.

LIVE MODE: like Kangaroo, the cloth is stepped a few iterations per frame
and redrawn, so you watch it settle and can drag sliders while it runs.
  - Shape inputs (rail, drop, fullness, pleats, pleat_depth, variation,
    res_u, res_v, seed) rebuild the cloth from scratch.
  - Physics inputs (gravity, stiffness, bend, wind, colliders, thickness,
    floor_z) apply to the running cloth without a restart.

INPUTS  (right-click each input -> set Access and Type hint)
    run         bool          item   Toggle: True = simulate live            (default True)
    reset       bool          item   Button: restart from the flat fabric
    grab        bool          item   Toggle: drag the cloth with the left mouse button
    rail        Curve         item   Top edge the fabric hangs from (default: 600 line at z=400)
    drop        float         item   Hanging length                          (default 400)
    fullness    float         item   Fabric width / rail length, >1 = folds  (default 1.6)
    pleats      int           item   Number of folds along the rail          (default 8)
    pleat_depth float         item   Fold depth, 0 = auto from fullness      (default 0)
    variation   float         item   0-1, randomness of fold depth           (default 0.3)
    res_u       int           item   Points across                           (default 61)
    res_v       int           item   Points down                             (default 40)
    steps       int           item   Iterations per frame                    (default 5)
    iterations  int           item   Max iterations after a change           (default 400)
    gravity     float         item   Pull per step, relative to drop         (default 0.1)
    stiffness   float         item   0-1, stretch resistance                 (default 1.0)
    bend        float         item   0-1, bending resistance                 (default 0.3)
    wind        Vector3d      item   Optional push, relative like gravity, e.g. (0, 0.03, 0)
    colliders   GeometryBase  list   Optional Mesh / Brep / SubD / Surface obstacles
    thickness   float         item   Fabric thickness in model units         (default drop / 40)
    floor_z     float         item   Optional floor height for pooling fabric
    seed        int           item   Random seed                             (default 1)

OUTPUTS
    mesh        Quad mesh of the fabric's middle surface (coarser while running)
    thick_mesh  Closed mesh with the fabric thickness (when settled / paused)
    surface     NURBS surface (when settled / paused)
    points      Simulated points (when settled / paused)
    info        Status text
"""
import math
import time
import numpy as np

# Bump this whenever the Cloth / FieldCollider code changes. Saved simulation
# state from another version is thrown away instead of reused (objects saved
# by an older script don't have the new attributes and would crash).
STATE_VERSION = 4

# ----------------------------------------------------------------------------
# core (NumPy only, no Rhino needed)
# ----------------------------------------------------------------------------

def _fold_offsets(amp, t, pleats, mult):
    return amp * mult * np.sin(2.0 * np.pi * pleats * t)


def auto_amplitude(rail_len, fullness, nu, pleats, mult):
    """Fold depth so that the folded fabric width equals rail_len * fullness."""
    target = rail_len * fullness
    dl = rail_len / (nu - 1)
    t = np.linspace(0.0, 1.0, nu)
    lo, hi = 0.0, rail_len
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        o = _fold_offsets(mid, t, pleats, mult)
        if np.sqrt(dl * dl + np.diff(o) ** 2).sum() < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def _constraint_sets(nu, nv, du, dv):
    """Springs as pairs of array slices. Every set touches each point at most
    once, so a whole set can be solved at once with array maths and still
    behave like a one-spring-at-a-time (Gauss-Seidel) solver."""
    A = slice(None)
    diag = math.hypot(du, dv)
    sets = []
    for r in (0, 1):                                  # stretch across / down
        sets.append(((A, slice(r, nu - 1, 2), A), (A, slice(r + 1, nu, 2), A), du, "stretch"))
        sets.append(((A, A, slice(r, nv - 1, 2)), (A, A, slice(r + 1, nv, 2)), dv, "stretch"))
    for r in (0, 1):                                  # shear diagonals
        sets.append(((A, slice(r, nu - 1, 2), slice(0, nv - 1)),
                     (A, slice(r + 1, nu, 2), slice(1, nv)), diag, "shear"))
        sets.append(((A, slice(r, nu - 1, 2), slice(1, nv)),
                     (A, slice(r + 1, nu, 2), slice(0, nv - 1)), diag, "shear"))
    return sets


def _length(d):
    """Length of vectors stored as d[0], d[1], d[2] (keeps a leading axis)."""
    return np.sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2])[None]


class FieldCollider(object):
    """Collider stored as a signed distance grid: phi < 0 inside, > 0 outside.
    Sampling it is pure NumPy, so every point can be tested every step."""

    def __init__(self, phi, origin, cell, friction=0.5):
        phi = np.asarray(phi, float)
        gx, gy, gz = np.gradient(phi, cell)
        self.F = np.stack([phi, gx, gy, gz], axis=-1)   # value + gradient
        self.origin = np.asarray(origin, float)
        self.cell = float(cell)
        self.shape = np.array(phi.shape)
        self.friction = friction

    def sample(self, X):
        """Trilinear phi and gradient at points X (N, 3). Points outside the
        grid get phi = inf (no collision)."""
        q = (X - self.origin) / self.cell
        i0 = np.floor(q).astype(int)
        inside = np.all((i0 >= 0) & (i0 < self.shape - 1), axis=1)
        i0 = np.clip(i0, 0, self.shape - 2)
        f = np.clip(q - i0, 0.0, 1.0)
        out = 0.0
        for dx in (0, 1):
            wx = f[:, 0] if dx else 1 - f[:, 0]
            for dy in (0, 1):
                wy = f[:, 1] if dy else 1 - f[:, 1]
                for dz in (0, 1):
                    wz = f[:, 2] if dz else 1 - f[:, 2]
                    out = out + (wx * wy * wz)[:, None] * \
                        self.F[i0[:, 0] + dx, i0[:, 1] + dy, i0[:, 2] + dz]
        out[~inside, 0] = np.inf
        return out[:, 0], out[:, 1:]

    def __call__(self, P, prev, free, gap):
        """P, prev: (3, nu, nv) arrays, edited in place."""
        X = P.reshape(3, -1)
        V = prev.reshape(3, -1)
        cand = free.reshape(-1).copy()
        # cheap box test first: only points near the grid are sampled
        lo = self.origin
        hi = self.origin + (self.shape - 1) * self.cell
        for a in range(3):
            cand &= (X[a] > lo[a]) & (X[a] < hi[a])
        idx = np.nonzero(cand)[0]
        if idx.size == 0:
            return
        Xi = X[:, idx].T
        phi, g = self.sample(Xi)
        hit = phi < gap
        if not hit.any():
            return
        idx, Xi, phi, n = idx[hit], Xi[hit], phi[hit], g[hit]
        n /= np.maximum(np.linalg.norm(n, axis=1, keepdims=True), 1e-12)
        vel = Xi - V[:, idx].T
        Xi = Xi + (gap - phi)[:, None] * n                # push out to the gap
        vn = (vel * n).sum(1, keepdims=True)
        vt = vel - vn * n                                 # keep some sliding
        X[:, idx] = Xi.T
        V[:, idx] = (Xi - vt * (1.0 - self.friction)).T


class Cloth(object):
    """Points are stored as P[xyz, i, j]: i across the rail, j down."""

    def __init__(self, top, side, rail_len, drop, fullness, pleats,
                 pleat_depth, variation, res_v, seed):
        """top, side: (nu, 3) arrays of rail points and sideways unit normals."""
        top = np.asarray(top, float)
        side = np.asarray(side, float)
        rnd = np.random.default_rng(seed)
        nu, nv = len(top), int(res_v)
        pleats = max(0, min(int(pleats), max(1, (nu - 1) // 4)))  # >= 4 points per fold

        # the folded top edge is built on a dense version of the rail ...
        nd = 20 * (nu - 1) + 1
        t = np.linspace(0.0, 1.0, nd)
        tr = np.linspace(0.0, 1.0, nu)
        top_d = np.stack([np.interp(t, tr, top[:, a]) for a in range(3)], 1)
        side_d = np.stack([np.interp(t, tr, side[:, a]) for a in range(3)], 1)
        side_d /= np.maximum(np.linalg.norm(side_d, axis=1, keepdims=True), 1e-12)
        # smooth random variation of fold depth along the rail
        ph1, ph2 = rnd.uniform(0, 2 * np.pi, 2)
        mult = 1.0 + variation * 0.5 * (np.sin(2 * np.pi * 1.3 * t + ph1) +
                                        np.sin(2 * np.pi * 2.7 * t + ph2))
        if pleats == 0:
            amp = 0.0
        elif pleat_depth is None or pleat_depth <= 0:
            amp = auto_amplitude(rail_len, fullness, nd, pleats, mult)
        else:
            amp = float(pleat_depth)
        dense = top_d + side_d * _fold_offsets(amp, t, pleats, mult)[:, None]
        # ... then resampled so neighbouring points are equally far apart
        # along the fabric: the cloth starts at exactly its rest size
        seg = np.linalg.norm(np.diff(dense, axis=0), axis=1)
        cum = np.concatenate([[0.0], np.cumsum(seg)])
        at = np.linspace(0.0, cum[-1], nu)
        base = np.stack([np.interp(at, cum, dense[:, a]) for a in range(3)], 1)

        du = cum[-1] / (nu - 1)
        dv = drop / (nv - 1)
        P = np.repeat(base.T[:, :, None], nv, axis=2)
        P[2] -= np.arange(nv)[None, :] * dv

        self.P = np.ascontiguousarray(P)  # (3, nu, nv) positions
        self.prev = self.P.copy()         # previous positions (Verlet)
        self.W = np.ones((1, nu, nv))     # inverse mass, 0 = pinned
        self.W[:, :, 0] = 0.0             # top row hangs from the rail
        self.tether = (np.arange(nv) * dv)[None, None, :] * 1.0001
        self.sets = _constraint_sets(nu, nv, du, dv)
        self.nu, self.nv, self.du, self.dv = nu, nv, du, dv
        self.scale = drop / 40.0          # unit-free gravity / wind reference
        self.amp, self.pleats = amp, pleats
        self.total = 0                    # iterations since the last change
        self.motion = float("inf")        # largest move in the last step

    @property
    def points(self):
        """Positions as (nu, nv, 3) (a view, no copy)."""
        return np.moveaxis(self.P, 0, -1)

    def _project(self, sa, sb, rest, k):
        a, b = self.P[sa], self.P[sb]     # views: edits write into self.P
        wa, wb = self.W[sa], self.W[sb]
        d = b - a
        L = _length(d)
        ws = wa + wb
        s = k * (L - rest) / (np.maximum(L, 1e-12) * np.where(ws > 0, ws, 1.0))
        corr = d * s
        a += wa * corr
        b -= wb * corr

    def _bend(self, span, k):
        """Bending: springs between points `span` apart, across and down.
        They only push apart (never pull), so they resist folds tighter than
        about span x spacing but leave big, round folds alone. The long span
        is 2 x thickness, so thicker fabric makes fewer, rounder folds."""
        P, W = self.P, self.W
        for axis, rest in ((1, span * self.du), (2, span * self.dv)):
            n = P.shape[axis]
            if n <= span:
                continue
            ia = [slice(None)] * 3
            ib = [slice(None)] * 3
            ia[axis] = slice(0, n - span)
            ib[axis] = slice(span, n)
            ia, ib = tuple(ia), tuple(ib)
            d = P[ib] - P[ia]
            L = _length(d)
            wa, wb = W[ia], W[ib]
            ws = wa + wb
            s = 0.5 * k * np.minimum(L - rest, 0.0) / \
                (np.maximum(L, 1e-12) * np.where(ws > 0, ws, 1.0))
            corr = d * s
            P[ia] += wa * corr
            P[ib] -= wb * corr

    def _self_collide(self, thickness, free):
        """Keeps folds `thickness` apart: points in the same row that are far
        apart along the fabric but close in space push each other away."""
        if thickness <= 0:
            return
        P = self.P
        m0 = max(2, int(math.ceil(1.7 * thickness / self.du)))
        m1 = min(self.nu - 1, m0 + min(30, max(4, int(math.ceil(4 * thickness / self.du)))))
        if m0 > m1:
            return
        delta = np.zeros_like(P)
        count = np.zeros_like(self.W)
        f = free.astype(float)
        for m in range(m0, m1 + 1):
            d = P[:, m:] - P[:, :-m]
            L = _length(d)
            close = (L < thickness) & ((f[:, m:] + f[:, :-m]) > 0)
            if not close.any():
                continue
            corr = np.where(close, d * ((thickness - L) / np.maximum(L, 1e-12) * 0.5), 0.0)
            delta[:, m:] += corr * f[:, m:]
            delta[:, :-m] -= corr * f[:, :-m]
            count[:, m:] += close
            count[:, :-m] += close
        P += delta / np.maximum(count, 1.0)

    def step(self, n, accel, stiffness, bend, thickness, floor_z=None,
             collide=None, damping=0.98, passes=4, max_speed=None, grab=None):
        """accel is relative (x scale). grab: ((i, j), target_xyz)."""
        P, W = self.P, self.W
        accel = (np.asarray(accel, float) * self.scale)[:, None, None]
        k_of = {"stretch": stiffness, "shear": stiffness * 0.5}
        span = max(2, min(int(round(2.0 * thickness / self.du)), (self.nu - 1) // 3))
        gap = 0.5 * thickness
        held_w = None
        if grab is not None:
            (gi, gj), target = grab
            held_w = W[0, gi, gj]
            W[0, gi, gj] = 0.0                    # held points act pinned
            self.prev[:, gi, gj] = P[:, gi, gj]   # so it flings when released
            P[:, gi, gj] = target
        free = W > 0
        try:
            for it in range(n):
                before = P.copy()
                vel = (P - self.prev) * damping + accel
                if max_speed:
                    # never move further than max_speed per step, so points
                    # can't jump through a thin collider between checks
                    sp = _length(vel)
                    vel *= np.minimum(1.0, max_speed / np.maximum(sp, 1e-12))
                self.prev[:] = P
                P += vel * W
                for _ in range(passes):
                    for sa, sb, rest, kind in self.sets:
                        self._project(sa, sb, rest, k_of[kind])
                    if bend > 0:
                        self._bend(2, bend)             # smooths small creases
                        if span > 2:
                            self._bend(span, bend)      # big, thick folds
                self._self_collide(thickness, free)
                # tethers: no point gets further from its top point than the
                # fabric length between them (stops rubbery stretching)
                d = P - P[:, :, :1]
                L = _length(d)
                over = (L > self.tether) & free
                if over.any():
                    P[:] = np.where(over, P[:, :, :1] + d * (self.tether / np.maximum(L, 1e-12)), P)
                if floor_z is not None:
                    fz = floor_z + gap
                    low = (P[2:3] < fz) & free
                    if low.any():
                        P[2:3] = np.where(low, fz, P[2:3])
                        self.prev[:] = np.where(low, P, self.prev)   # friction
                if collide is not None:
                    collide(P, self.prev, free, gap)
                self.total += 1
                self.motion = float(np.abs(P - before).max())
        finally:
            if held_w is not None:
                W[0, gi, gj] = held_w

    def max_strain(self):
        """Largest stretch (tension) of the stretch springs, as a ratio."""
        worst = 0.0
        for sa, sb, rest, kind in self.sets:
            if kind != "stretch":
                continue
            L = _length(self.P[sb] - self.P[sa])[0]
            live = (self.W[sa] + self.W[sb])[0] > 0
            if live.any():
                worst = max(worst, float(((L[live] - rest) / rest).max()))
        return worst

    def nearest_on_screen(self, M, x, y, radius):
        """Index of the cloth point closest to screen pixel (x, y), using the
        4x4 world->screen matrix M. None if nothing within radius pixels."""
        pts = self.points.reshape(-1, 3)
        h = pts @ M[:3, :3].T + M[:3, 3]
        w = pts @ M[3, :3] + M[3, 3]
        ok = w > 1e-9                                # in front of the camera
        sx = h[:, 0] / np.where(ok, w, 1.0)
        sy = h[:, 1] / np.where(ok, w, 1.0)
        d2 = np.where(ok, (sx - x) ** 2 + (sy - y) ** 2, np.inf)
        k = int(np.argmin(d2))
        if d2[k] > radius * radius:
            return None
        return divmod(k, self.nv)


def grid_indices(n, stride):
    """Every stride-th index, always including the last one."""
    idx = list(range(0, n, stride))
    if idx[-1] != n - 1:
        idx.append(n - 1)
    return np.array(idx)

# ----------------------------------------------------------------------------
# grasshopper
# ----------------------------------------------------------------------------
try:
    import Rhino
    import Rhino.UI
    import Rhino.Geometry as rg
    import scriptcontext as sc
    import Grasshopper
    IN_RHINO = True
except ImportError:
    IN_RHINO = False


def _val(name, default):
    v = globals().get(name)
    return default if v is None else v


def _to_geometry(g):
    """Rhino object ids (Guid) -> geometry, Grasshopper wrappers -> geometry."""
    if g is None:
        return None
    if hasattr(g, "Value") and not isinstance(g, rg.GeometryBase):
        g = g.Value                                   # GH_Brep, GH_Mesh, ...
    if not isinstance(g, rg.GeometryBase):
        try:
            obj = Rhino.RhinoDoc.ActiveDoc.Objects.FindId(g)
            g = obj.Geometry if obj else None
        except Exception:
            g = None
    return g


def _collider_mesh(geos):
    """Joins every collider into one mesh. Returns (mesh, report text)."""
    joined = rg.Mesh()
    used, skipped = 0, []
    for g in geos or []:
        g = _to_geometry(g)
        if g is None:
            continue
        if isinstance(g, (rg.Extrusion, rg.Surface)):
            g = g.ToBrep()
        part = None
        if isinstance(g, rg.Mesh):
            part = g
        elif isinstance(g, rg.SubD):
            part = rg.Mesh.CreateFromSubD(g, 2)
        elif isinstance(g, rg.Brep):
            part = rg.Mesh()
            for m in rg.Mesh.CreateFromBrep(g, rg.MeshingParameters.QualityRenderMesh) or []:
                part.Append(m)
        if part is not None and part.Vertices.Count > 0:
            joined.Append(part)
            used += 1
        else:
            skipped.append(type(g).__name__)
    if joined.Vertices.Count == 0:
        return None, "colliders: none" + (" (unsupported: {})".format(
            ", ".join(skipped)) if skipped else "")
    if joined.SolidOrientation() == -1:               # closed but inside-out
        joined.Flip(True, True, True)
    joined.Normals.ComputeNormals()
    joined.FaceNormals.ComputeFaceNormals()
    txt = "colliders: {} ({} faces)".format(used, joined.Faces.Count)
    if skipped:
        txt += " skipped: " + ", ".join(skipped)
    return joined, txt


def _build_field(cmesh, pad, max_nodes=40000):
    """Samples the signed distance to the collider on a grid around it, once.
    After that, collisions for every point are pure NumPy (no Rhino calls)."""
    bb = cmesh.GetBoundingBox(True)
    lo = np.array([bb.Min.X, bb.Min.Y, bb.Min.Z]) - pad
    hi = np.array([bb.Max.X, bb.Max.Y, bb.Max.Z]) + pad
    size = hi - lo
    cell = (np.prod(size) / max_nodes) ** (1.0 / 3.0)
    cell = max(cell, size.max() / 200.0)
    shape = np.maximum(2, np.ceil(size / cell).astype(int) + 1)
    band = pad * 2.0
    phi = np.full(shape, band)
    for i in range(shape[0]):
        x = lo[0] + i * cell
        for j in range(shape[1]):
            y = lo[1] + j * cell
            for k in range(shape[2]):
                z = lo[2] + k * cell
                mp = cmesh.ClosestMeshPoint(rg.Point3d(x, y, z), band)
                if mp is None:
                    continue
                cp = mp.Point
                nrm = cmesh.NormalAt(mp)
                nrm.Unitize()
                d = math.sqrt((x - cp.X) ** 2 + (y - cp.Y) ** 2 + (z - cp.Z) ** 2)
                side = (x - cp.X) * nrm.X + (y - cp.Y) * nrm.Y + (z - cp.Z) * nrm.Z
                phi[i, j, k] = d if side >= 0 else -d
    field = FieldCollider(phi, lo, cell)
    field.pad = pad
    return field


def _rail_frames(crv, nu):
    top, side = [], []
    for t in crv.DivideByCount(nu - 1, True):
        p = crv.PointAt(t)
        nrm = rg.Vector3d.CrossProduct(crv.TangentAt(t), rg.Vector3d.ZAxis)
        if not nrm.Unitize():
            nrm = rg.Vector3d.YAxis
        top.append((p.X, p.Y, p.Z))
        side.append((nrm.X, nrm.Y, nrm.Z))
    return top, side


def _geo_key(geos):
    key = []
    for g in geos or []:
        g = _to_geometry(g)
        if g is not None:
            b = g.GetBoundingBox(False)
            key.append((type(g).__name__,) + tuple(round(v, 5) for v in
                       (b.Min.X, b.Min.Y, b.Min.Z, b.Max.X, b.Max.Y, b.Max.Z)))
    return tuple(key)


def _grid_mesh(st, P):
    """Quad mesh for a (a, b, 3) point grid. Faces are built once per grid
    size and reused; each frame only the vertices are moved."""
    a, b = P.shape[:2]
    cache = st.get("mesh_cache")
    if cache is None or cache[0] != (a, b):
        m = rg.Mesh()
        for _ in range(a * b):
            m.Vertices.Add(0.0, 0.0, 0.0)
        for i in range(a - 1):
            for j in range(b - 1):
                p, q = i * b + j, (i + 1) * b + j
                m.Faces.AddFace(p, q, q + 1, p + 1)
        cache = ((a, b), m)
        st["mesh_cache"] = cache
    m = cache[1].Duplicate()
    for k, (x, y, z) in enumerate(P.reshape(-1, 3).tolist()):
        m.Vertices.SetVertex(k, x, y, z)
    m.Normals.ComputeNormals()
    return m


def _schedule_next(comp, ms):
    doc = comp.OnPingDocument()
    if doc is None:
        return
    def again(_doc):
        comp.ExpireSolution(False)
    try:
        doc.ScheduleSolution(ms, Grasshopper.Kernel.GH_Document.GH_ScheduleDelegate(again))
    except Exception:
        doc.ScheduleSolution(ms, again)

# --- cursor grabbing ---------------------------------------------------------
# Rhino sends viewport mouse events to a MouseCallback. The class is created
# only once per Rhino session (a .NET type can't be defined twice), and it just
# forwards events to handler functions that every script run replaces.

def _mouse_callback():
    cb = sc.sticky.get("drape_sim_mouse")
    if cb is not None:
        return cb

    def forward(kind, e):
        h = sc.sticky.get("drape_sim_handlers", {}).get(kind)
        if h is not None:
            try:
                h(e)
            except Exception as ex:
                Rhino.RhinoApp.WriteLine("drape grab: {}".format(ex))

    class DrapeMouse(Rhino.UI.MouseCallback):
        __namespace__ = "LocustaDrape"

        def OnMouseDown(self, e):
            forward("down", e)

        def OnMouseMove(self, e):
            forward("move", e)

        def OnMouseUp(self, e):
            forward("up", e)

    cb = DrapeMouse()
    sc.sticky["drape_sim_mouse"] = cb
    return cb


def _is_left(e):
    try:
        return str(e.MouseButton).endswith("Left")
    except Exception:
        return str(e.Button).endswith("Left")


def _install_grab(comp, st):
    def wake():
        cloth = st.get("cloth")
        if cloth is not None:
            cloth.total, cloth.motion = 0, float("inf")
        if not st.get("running"):          # already looping otherwise
            comp.ExpireSolution(True)

    def down(e):
        cloth = st.get("cloth")
        if cloth is None or not _is_left(e):
            return
        vp = e.View.ActiveViewport
        xf = vp.GetTransform(rg.CoordinateSystem.World, rg.CoordinateSystem.Screen)
        M = np.array([[getattr(xf, "M{}{}".format(r, c)) for c in range(4)]
                      for r in range(4)])
        hit = cloth.nearest_on_screen(M, e.ViewportPoint.X, e.ViewportPoint.Y, 15)
        if hit is None:
            return                         # missed the cloth: normal Rhino click
        x, y, z = cloth.points[hit]
        st["grab"] = {"idx": hit, "target": np.array([x, y, z]),
                      "plane": rg.Plane(rg.Point3d(x, y, z), vp.CameraDirection)}
        e.Cancel = True                    # don't start a selection window
        wake()

    def move(e):
        g = st.get("grab")
        if g is None:
            return
        ok, line = e.View.ActiveViewport.GetFrustumLine(e.ViewportPoint.X,
                                                        e.ViewportPoint.Y)
        if ok:
            hit, t = rg.Intersection.LinePlane(line, g["plane"])
            if hit:
                p = line.PointAt(t)
                g["target"] = np.array([p.X, p.Y, p.Z])
        e.Cancel = True
        wake()

    def up(e):
        if st.pop("grab", None) is not None:
            e.Cancel = True
            wake()

    sc.sticky["drape_sim_handlers"] = {"down": down, "move": move, "up": up}


def run_grasshopper():
    comp = ghenv.Component
    key = "drape_sim_state_v{}".format(STATE_VERSION)
    for k in list(sc.sticky.keys()):
        if (k == "drape_sim" or str(k).startswith("drape_sim_state")) and k != key:
            sc.sticky.pop(k, None)          # drop state saved by older versions
    store = sc.sticky.setdefault(key, {})
    sid = str(comp.InstanceGuid)
    st = store.get(sid) or {}

    drop = float(_val("drop", 400.0))
    crv = _val("rail", None)
    if crv is None:
        crv = rg.LineCurve(rg.Point3d(0, 0, drop), rg.Point3d(1.5 * drop, 0, drop))
    nu = max(5, int(_val("res_u", 61)))
    nv = max(3, int(_val("res_v", 40)))
    top, side = _rail_frames(crv, nu)
    shape = (drop, float(_val("fullness", 1.6)),
             int(_val("pleats", 8)), float(_val("pleat_depth", 0.0)),
             float(_val("variation", 0.3)), nv, int(_val("seed", 1)))
    shape_key = (tuple(round(v, 6) for p in top for v in p), shape)

    cloth = st.get("cloth")
    if cloth is None or bool(_val("reset", False)) or st.get("shape_key") != shape_key:
        cloth = Cloth(top, side, crv.GetLength(), *shape)
        keep = {k: st[k] for k in ("field", "cmesh", "ctext", "ckey") if k in st}
        st.clear()                          # keep the same dict: grab handlers use it
        st.update(keep)                     # the collider field doesn't depend on the cloth
        st.update({"cloth": cloth, "shape_key": shape_key})
        store[sid] = st

    # physics: applied live, any change wakes the cloth up again
    w = _val("wind", None)
    g = float(_val("gravity", 0.1))
    accel = (w.X if w else 0.0, w.Y if w else 0.0, (w.Z if w else 0.0) - g)
    stiffness = float(_val("stiffness", 1.0))
    bend = float(_val("bend", 0.3))
    floor_z = _val("floor_z", None)
    thickness = max(1e-6, float(_val("thickness", drop / 40.0)))
    gap = 0.5 * thickness

    geos = _val("colliders", None)
    if geos is not None and not isinstance(geos, (list, tuple)) and \
            (isinstance(geos, rg.GeometryBase) or not hasattr(geos, "__iter__")):
        geos = [geos]                       # item access instead of list
    ckey = _geo_key(geos)
    if st.get("ckey") != ckey:
        st["cmesh"], st["ctext"] = _collider_mesh(geos)
        st["ckey"] = ckey
        st["field"] = None
    cmesh = st.get("cmesh")
    field = st.get("field")
    build_s = None
    if cmesh is not None and (field is None or field.pad < 2.0 * gap + 2.0 * cloth.dv):
        # (re)build the distance grid only when the collider changes or the
        # gap outgrows the grid's margin; moving sliders doesn't rebuild it
        t0 = time.time()
        diag = cmesh.GetBoundingBox(True).Diagonal.Length
        field = _build_field(cmesh, max(0.1 * diag, 3.0 * gap + 2.0 * cloth.dv))
        st["field"] = field
        build_s = time.time() - t0
        st["ctext_build"] = " | field {}x{}x{} built in {:.1f} s".format(
            field.shape[0], field.shape[1], field.shape[2], build_s)
    max_speed = max(gap, 0.25 * cloth.dv) if field is not None else None
    warn = ""
    if field is not None:
        if cloth.total == 0 or build_s is not None:
            phi0 = field.sample(cloth.points.reshape(-1, 3))[0]
            st["inside"] = int((phi0 < 0).sum())
        if st.get("inside"):
            warn = (" | WARNING: fabric started inside the colliders at {} points, "
                    "it may get pierced. Move the rail behind the frame").format(st["inside"])

    phys_key = (accel, stiffness, bend, floor_z, thickness, ckey)
    if st.get("phys_key") != phys_key:
        st["phys_key"] = phys_key
        cloth.total = 0
        cloth.motion = float("inf")

    # cursor grabbing
    grab_on = bool(_val("grab", False))
    grab_err = ""
    try:
        cb = _mouse_callback()
        if grab_on:
            _install_grab(comp, st)
        cb.Enabled = grab_on
    except Exception as ex:
        grab_on = False
        grab_err = " | grab unavailable: {}".format(ex)
    if not grab_on:
        st.pop("grab", None)
    held = st.get("grab")
    grab = (held["idx"], held["target"]) if held else None

    run = bool(_val("run", True))
    max_it = max(1, int(_val("iterations", 400)))

    def is_settled():
        return held is None and (cloth.total >= max_it or cloth.motion < 1e-4 * cloth.dv)

    settled = is_settled()
    t0 = time.time()
    if run and not settled:
        cloth.step(max(1, int(_val("steps", 5))), accel, stiffness, bend, thickness,
                   floor_z=floor_z, collide=field, max_speed=max_speed, grab=grab)
        settled = is_settled()
    ms = (time.time() - t0) * 1000.0

    still = settled or not run
    if still:
        P = cloth.points
    else:                            # coarser preview while it moves: lighter
        stride = max(1, int(math.ceil(max(cloth.nu, cloth.nv) / 50.0)))
        P = cloth.points[np.ix_(grid_indices(cloth.nu, stride), grid_indices(cloth.nv, stride))]
    m = _grid_mesh(st, P)
    thick = srf = pts = None
    if still:                        # heavier outputs only once it stops moving
        pts = [rg.Point3d(x, y, z) for x, y, z in cloth.points.reshape(-1, 3).tolist()]
        srf = rg.NurbsSurface.CreateFromPoints(pts, cloth.nu, cloth.nv, 3, 3)
        try:
            inner = m.Offset(-gap)
            thick = inner.Offset(thickness, True) if inner else None
        except Exception:
            thick = None

    state = "settled" if settled else ("running" if run else "paused")
    if held is not None:
        state = "grabbing"
    txt = ("{} | iteration {} | {:.0f} ms/frame | points {}x{} | folds {} | "
           "fold depth {:.1f} | thickness {:g} | max stretch {:.1%} | {}{}{}").format(
        state, cloth.total, ms, cloth.nu, cloth.nv, cloth.pleats, cloth.amp,
        thickness, cloth.max_strain(), st.get("ctext", "colliders: none"),
        st.get("ctext_build", "") if field is not None else "",
        " | grab: drag the cloth with the left mouse button" if grab_on else grab_err) + warn
    comp.Message = state
    st["running"] = run and not settled
    if st["running"]:
        _schedule_next(comp, 1)
    return m, thick, srf, pts, txt


if IN_RHINO:
    mesh, thick_mesh, surface, points, info = run_grasshopper()
