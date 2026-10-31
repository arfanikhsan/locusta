# r: numpy
"""
Drape Simulation (live) - Rhino 8 "Python 3 Script" component

A cloth hangs from a top rail. The fabric is wider than the rail
(fullness > 1), so it gathers into vertical folds, falls under gravity,
wraps over collider objects and can pool on a floor.

RAIL: a Curve gives a thin curtain. A planar strip Surface gives THICK
fabric: the strip's long direction is the rail (its centre line) and its
short width is the fabric thickness. Thick fabric can't bend tighter than
its thickness, so folds get bigger and rounder (cartoon fabric), folds stay
that far apart, and the `thick_mesh` output shows the fabric's thickness.

UNITS: everything scales with your model, so cm, m or mm all work. Gravity
and wind are relative to the drop length.

START: start = 1 or -1 drops the fabric from a flat, horizontal sheet at
rail height (to one side or the other), so you see it fall and swing down
under gravity. start = 0 starts it already hanging.

LIVE MODE: like Kangaroo, the cloth is stepped a few iterations per frame
and redrawn, so you watch it settle and can drag sliders while it runs.
  - Shape inputs (rail, drop, fullness, pleats, pleat_depth, variation,
    res_u, res_v, seed, start) rebuild the cloth from scratch.
  - Physics inputs (gravity, stiffness, bend, wind, colliders, thickness,
    floor_z, anchors) apply to the running cloth without a restart.

INPUTS  (right-click each input -> set Access and Type hint)
    run         bool          item   Toggle: True = simulate live            (default True)
    reset       bool          item   Button: restart the fabric
    grab        bool          item   Toggle: drag the cloth with the left mouse button
    anchors     Point3d       list   Optional points: the nearest fabric point is held at each
    rail        GeometryBase  item   Curve (thin) or planar strip Surface (thick)
                                     (default: 600 line at z=400)
    start       int           item   1 / -1 = fall from horizontal, 0 = hanging (default 1)
    drop        float         item   Hanging length                          (default 400)
    fullness    float         item   Fabric width / rail length, >1 = folds  (default 1.6)
    pleats      int           item   Number of folds along the rail          (default 8)
    pleat_depth float         item   Fold depth, 0 = auto from fullness      (default 0)
    variation   float         item   0-1, randomness of fold depth           (default 0.3)
    res_u       int           item   Points across                           (default 61)
    res_v       int           item   Points down                             (default 40)
    steps       int           item   Iterations per frame                    (default 5)
    iterations  int           item   Max iterations after a change           (default 1000)
    gravity     float         item   Pull per step, relative to drop         (default 0.1)
    stiffness   float         item   0-1, stretch resistance                 (default 1.0)
    bend        float         item   0-1, bending resistance                 (default 0.3)
    wind        Vector3d      item   Optional push, relative like gravity, e.g. (0, 0.03, 0)
    colliders   GeometryBase  list   Optional Mesh / Brep / SubD / Surface obstacles
    thickness   float         item   Gap kept between fabric and colliders/floor
                                     in model units                          (default drop / 80)
    floor_z     float         item   Optional floor height for pooling fabric
    seed        int           item   Random seed                             (default 1)

OUTPUTS
    mesh        Quad mesh of the fabric (middle surface; coarser while running)
    thick_mesh  Closed mesh with the fabric thickness (surface rail, when settled)
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
STATE_VERSION = 5

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
        self.R = np.eye(3)               # rigid move of the collider since
        self.t = np.zeros(3)             # the grid was built: x = R @ local + t

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
        Xl = self.R.T @ (X - self.t[:, None])            # into the grid's frame
        # cheap box test first: only points near the grid are sampled
        lo = self.origin
        hi = self.origin + (self.shape - 1) * self.cell
        for a in range(3):
            cand &= (Xl[a] > lo[a]) & (Xl[a] < hi[a])
        idx = np.nonzero(cand)[0]
        if idx.size == 0:
            return
        phi, g = self.sample(Xl[:, idx].T)
        hit = phi < gap
        if not hit.any():
            return
        idx, phi, n = idx[hit], phi[hit], g[hit] @ self.R.T  # normals to world
        Xi = X[:, idx].T
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
                 pleat_depth, variation, res_v, seed, start=0):
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
        off_d = _fold_offsets(amp, t, pleats, mult)
        dense = top_d + side_d * off_d[:, None]
        # ... then resampled so neighbouring points are equally far apart
        # along the fabric: the cloth starts at exactly its rest size
        seg = np.linalg.norm(np.diff(dense, axis=0), axis=1)
        cum = np.concatenate([[0.0], np.cumsum(seg)])
        at = np.linspace(0.0, cum[-1], nu)
        def resample(arr):
            return np.stack([np.interp(at, cum, arr[:, a]) for a in range(arr.shape[1])], 1)
        base = resample(dense)

        du = cum[-1] / (nu - 1)
        dv = drop / (nv - 1)
        P = np.repeat(base.T[:, :, None], nv, axis=2)
        if not start:
            P[2] -= np.arange(nv)[None, :] * dv           # hanging straight down
        else:
            # flat sheet sticking out sideways at rail height; the folds become
            # ripples up and down. Gravity swings it down from top to bottom.
            rail_r = resample(top_d)
            side_r = resample(side_d)
            side_r[:, 2] = 0.0
            side_r /= np.maximum(np.linalg.norm(side_r, axis=1, keepdims=True), 1e-12)
            off_r = resample(off_d[:, None])[:, 0]
            # start > 0: behind the rail (away from the side the fabric's
            # sideways normal points to), so it swings down without passing
            # through frames standing in front. start < 0: the other side.
            sign = -1.0 if start > 0 else 1.0
            j = np.arange(1, nv)[None, :]
            for a in range(3):
                P[a, :, 1:] = rail_r[:, a:a + 1] + sign * side_r[:, a:a + 1] * j * dv
            P[2, :, 1:] += off_r[:, None]

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
        self.age = 0                      # iterations since it was built

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

    def _bend(self, span, k, fade=None):
        """Bending: springs between points `span` apart, across and down.
        They only push apart (never pull), so they resist folds tighter than
        about span x spacing but leave big, round folds alone. The long span
        is 2 x thickness, so thicker fabric makes fewer, rounder folds.
        fade: optional (1, 1, nv) strength per row (weaker near the rail)."""
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
            if fade is not None:
                s = s * fade[ia]
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

    def step(self, n, accel, stiffness, bend, fabric=0.0, gap=0.0, floor_z=None,
             collide=None, damping=0.98, passes=4, max_speed=None, grabs=()):
        """accel is relative (x scale).
        fabric: fabric thickness (0 = thin cloth). gap: extra distance kept
        from colliders and the floor (on top of fabric / 2).
        grabs: list of ((i, j), target_xyz) points held in place."""
        P, W = self.P, self.W
        accel = (np.asarray(accel, float) * self.scale)[:, None, None]
        k_of = {"stretch": stiffness, "shear": stiffness * 0.5}
        span = max(2, min(int(round(2.0 * fabric / self.du)), (self.nu - 1) // 3))
        # thick-fabric bending fades in over 4 x thickness below the rail: up
        # there the fabric has to follow the pleats it hangs from, further
        # down it relaxes into big round folds (like a real thick curtain)
        fade = None
        if fabric > 0:
            depth = np.arange(self.nv) * self.dv
            fade = np.clip(depth / (4.0 * fabric), 0.0, 1.0)[None, None, :]
        gap = gap + 0.5 * fabric
        held = []
        for (gi, gj), target in grabs:
            held.append((gi, gj, W[0, gi, gj]))
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
                            self._bend(span, bend, fade)    # big, thick folds
                if fabric > 0:
                    self._self_collide(fabric, free)
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
                self.age += 1
                self.motion = float(np.abs(P - before).max())
        finally:
            for gi, gj, w in held:
                W[0, gi, gj] = w

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
            for m in rg.Mesh.CreateFromBrep(g, rg.MeshingParameters.Default) or []:
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


def _mesh_vertices(m):
    arr = m.Vertices.ToFloatArray()
    return np.fromiter(arr, dtype=float, count=len(arr)).reshape(-1, 3)


def rigid_fit(A, B):
    """Best rotation R and move t with B ~ A @ R.T + t, and the worst error."""
    ca, cb = A.mean(0), B.mean(0)
    U, _, Vt = np.linalg.svd((A - ca).T @ (B - cb))
    D = np.diag([1.0, 1.0, np.sign(np.linalg.det(Vt.T @ U.T))])
    R = Vt.T @ D @ U.T
    t = cb - R @ ca
    return R, t, float(np.abs(A @ R.T + t - B).max())


def _signed_push(cmesh, X, V, k, gap):
    """Pushes point k of X (3, N) out of the mesh to `gap`. Returns distance."""
    x, y, z = X[0, k], X[1, k], X[2, k]
    mp = cmesh.ClosestMeshPoint(rg.Point3d(x, y, z), 0.0)
    if mp is None:
        return float("inf")
    cp = mp.Point
    n = cmesh.NormalAt(mp)
    n.Unitize()
    dist = (x - cp.X) * n.X + (y - cp.Y) * n.Y + (z - cp.Z) * n.Z
    if dist < gap:
        X[0, k], X[1, k], X[2, k] = cp.X + n.X * gap, cp.Y + n.Y * gap, cp.Z + n.Z * gap
        V[:, k] = X[:, k]                              # friction
    return dist


def _make_direct(cmesh, reach, refresh=8):
    """Collider that asks Rhino directly (used right after the collider
    changes, while the distance grid is being built). Every `refresh` steps
    it finds the points that could reach the surface (they move at most
    reach / refresh per step); only those are tested every step."""
    bb = cmesh.GetBoundingBox(True)
    lo = np.array([bb.Min.X, bb.Min.Y, bb.Min.Z])
    hi = np.array([bb.Max.X, bb.Max.Y, bb.Max.Z])
    state = {"count": 0, "near": None}

    def collide(P, prev, free, gap):
        X, V, fr = P.reshape(3, -1), prev.reshape(3, -1), free.reshape(-1)
        pad = gap + reach
        if state["count"] % refresh == 0 or state["near"] is None:
            cand = fr.copy()
            for a in range(3):
                cand &= (X[a] > lo[a] - pad) & (X[a] < hi[a] + pad)
            near = []
            for k in np.nonzero(cand)[0]:
                pt = rg.Point3d(X[0, k], X[1, k], X[2, k])
                if cmesh.ClosestMeshPoint(pt, pad) is not None:   # fast reject
                    near.append(k)
            state["near"] = np.array(near, int)
        for k in state["near"]:
            if fr[k]:
                _signed_push(cmesh, X, V, k, gap)
        state["count"] += 1
    return collide


def _field_builder(cmesh, pad, max_nodes=30000, chunk=400):
    """Builds the distance grid a few hundred cells at a time (a generator),
    so a changed collider never freezes Grasshopper. Yields progress 0-1 and
    returns the FieldCollider."""
    bb = cmesh.GetBoundingBox(True)
    lo = np.array([bb.Min.X, bb.Min.Y, bb.Min.Z]) - pad
    hi = np.array([bb.Max.X, bb.Max.Y, bb.Max.Z]) + pad
    size = hi - lo
    cell = max((np.prod(size) / max_nodes) ** (1.0 / 3.0), size.max() / 200.0)
    nx, ny, nz = (int(v) for v in np.maximum(2, np.ceil(size / cell).astype(int) + 1))
    band = pad * 2.0
    phi = np.full((nx, ny, nz), band)
    total = nx * ny * nz
    for start in range(0, total, chunk):
        for flat in range(start, min(total, start + chunk)):
            i, rem = divmod(flat, ny * nz)
            j, k = divmod(rem, nz)
            x, y, z = lo[0] + i * cell, lo[1] + j * cell, lo[2] + k * cell
            mp = cmesh.ClosestMeshPoint(rg.Point3d(x, y, z), band)
            if mp is None:
                continue
            cp = mp.Point
            n = cmesh.NormalAt(mp)
            n.Unitize()
            d = math.sqrt((x - cp.X) ** 2 + (y - cp.Y) ** 2 + (z - cp.Z) ** 2)
            side = (x - cp.X) * n.X + (y - cp.Y) * n.Y + (z - cp.Z) * n.Z
            phi[i, j, k] = d if side >= 0 else -d
        yield float(start) / total
    field = FieldCollider(phi, lo, cell)
    field.pad = pad
    return field


def _rail_input(g):
    """Curve -> (curve, 0). Planar strip surface -> (centre line, width)."""
    g = _to_geometry(g)
    if g is None:
        return None, 0.0, ""
    if isinstance(g, rg.Curve):
        return g, 0.0, "thin (curve rail)"
    if isinstance(g, rg.Extrusion):
        g = g.ToBrep()
    if isinstance(g, rg.Brep):
        g = g.Faces[0] if g.Faces.Count else None
    if isinstance(g, rg.Surface):
        du, dv = g.Domain(0), g.Domain(1)
        cu = g.IsoCurve(0, dv.Mid)                    # runs along u
        cv = g.IsoCurve(1, du.Mid)                    # runs along v
        if cu is None or cv is None:
            return None, 0.0, ""
        if cu.GetLength() >= cv.GetLength():
            rail, width = cu, cv.GetLength()
        else:
            rail, width = cv, cu.GetLength()
        return rail, width, "thick (surface rail, fabric {:g})".format(width)
    return None, 0.0, "rail type {} not supported".format(type(g).__name__)


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
# What it sees is counted in GRAB_STATS and shown in `info`, to debug it.

def _grab_stats():
    return sc.sticky.setdefault("drape_sim_grab_stats",
                                {"events": 0, "downs": 0, "hits": 0, "error": ""})


def _mouse_callback():
    key = "drape_sim_mouse_v{}".format(STATE_VERSION)
    for k in list(sc.sticky.keys()):
        if str(k).startswith("drape_sim_mouse") and k != key:
            old = sc.sticky.pop(k, None)          # switch off hooks made by
            try:                                  # older versions of the script
                old.Enabled = False
            except Exception:
                pass
    cb = sc.sticky.get(key)
    if cb is not None:
        return cb

    def forward(kind, e):
        stats = sc.sticky.get("drape_sim_grab_stats")
        if stats is not None:
            stats["events"] += 1
        h = sc.sticky.get("drape_sim_handlers", {}).get(kind)
        if h is not None:
            try:
                h(e)
            except Exception as ex:
                if stats is not None:
                    stats["error"] = "{}: {}".format(kind, ex)
                Rhino.RhinoApp.WriteLine("drape grab: {}".format(ex))

    class DrapeMouse(Rhino.UI.MouseCallback):
        # a .NET type can only be defined once per namespace: new namespace
        # per script version
        __namespace__ = "LocustaDrape_v" + str(STATE_VERSION)

        def OnMouseDown(self, e):
            forward("down", e)

        def OnMouseMove(self, e):
            forward("move", e)

        def OnMouseUp(self, e):
            forward("up", e)

    cb = DrapeMouse()
    sc.sticky[key] = cb
    return cb


def _is_left(e):
    for attr in ("MouseButton", "Button"):
        try:
            b = getattr(e, attr)
        except Exception:
            continue
        try:
            if b == Rhino.UI.MouseButton.Left:
                return True
        except Exception:
            pass
        txt = str(b)
        # names, or the numbers pythonnet may show for the enums
        if txt.endswith("Left") or txt in ("1", "1048576"):
            return True
    return False


def _install_grab(comp, st):
    stats = _grab_stats()

    def wake():
        cloth = st.get("cloth")
        if cloth is not None:
            cloth.total, cloth.motion = 0, float("inf")
        if not st.get("running"):          # already looping otherwise
            comp.ExpireSolution(True)

    def down(e):
        stats["downs"] += 1
        cloth = st.get("cloth")
        if cloth is None or not _is_left(e):
            return
        vp = e.View.ActiveViewport
        xf = vp.GetTransform(rg.CoordinateSystem.World, rg.CoordinateSystem.Screen)
        try:
            M = np.array([[float(getattr(xf, "M{}{}".format(r, c))) for c in range(4)]
                          for r in range(4)])
        except Exception:
            M = np.array([[float(xf[r, c]) for c in range(4)] for r in range(4)])
        x0, y0 = float(e.ViewportPoint.X), float(e.ViewportPoint.Y)
        hit = cloth.nearest_on_screen(M, x0, y0, 20)
        if hit is None:
            return                         # missed the cloth: normal Rhino click
        stats["hits"] += 1
        x, y, z = (float(v) for v in cloth.points[hit])
        st["grab"] = {"idx": hit, "target": np.array([x, y, z]),
                      "plane": rg.Plane(rg.Point3d(x, y, z), vp.CameraDirection)}
        e.Cancel = True                    # don't start a selection window
        wake()

    def move(e):
        g = st.get("grab")
        if g is None:
            return
        res = e.View.ActiveViewport.GetFrustumLine(float(e.ViewportPoint.X),
                                                   float(e.ViewportPoint.Y))
        ok, line = res if isinstance(res, tuple) else (True, res)
        if ok:
            res = rg.Intersection.LinePlane(line, g["plane"])
            hit, t = res if isinstance(res, tuple) else (True, res)
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


def _anchor_grabs(st, cloth, anchors):
    """Each anchor point holds the fabric point nearest to it. The pairing
    is made once per set of anchors, so moving an anchor drags its point."""
    pts = []
    for a in anchors or []:
        a = _to_geometry(a) if not isinstance(a, rg.Point3d) else a
        if isinstance(a, rg.Point):
            a = a.Location
        if isinstance(a, rg.Point3d):
            pts.append((a.X, a.Y, a.Z))
    if not pts:
        st.pop("anchor_idx", None)
        return []
    if st.get("anchor_count") != len(pts) or st.get("anchor_idx") is None:
        flat = cloth.points.reshape(-1, 3)
        idx = []
        for p in pts:
            k = int(np.argmin(((flat - np.array(p)) ** 2).sum(1)))
            idx.append(divmod(k, cloth.nv))
        st["anchor_idx"], st["anchor_count"] = idx, len(pts)
    return [(ij, np.array(p)) for ij, p in zip(st["anchor_idx"], pts)]


def run_grasshopper():
    comp = ghenv.Component
    key = "drape_sim_state_v{}".format(STATE_VERSION)
    for k in list(sc.sticky.keys()):
        if (k == "drape_sim" or str(k).startswith("drape_sim_state")) and k != key:
            sc.sticky.pop(k, None)          # drop state saved by older versions
    store = sc.sticky.setdefault(key, {})
    sid = str(comp.InstanceGuid)
    st = store.get(sid) or {}

    # --- shape ---------------------------------------------------------------
    drop = float(_val("drop", 400.0))
    crv, fabric, mode = _rail_input(_val("rail", None))
    if crv is None:
        if not mode:
            mode = "thin (default rail)"
        crv = rg.LineCurve(rg.Point3d(0, 0, drop), rg.Point3d(1.5 * drop, 0, drop))
    nu = max(5, int(_val("res_u", 61)))
    nv = max(3, int(_val("res_v", 40)))
    top, side = _rail_frames(crv, nu)
    start = int(_val("start", 1))
    shape = (drop, float(_val("fullness", 1.6)),
             int(_val("pleats", 8)), float(_val("pleat_depth", 0.0)),
             float(_val("variation", 0.3)), nv, int(_val("seed", 1)), start)
    shape_key = (tuple(round(v, 6) for p in top for v in p), shape)

    cloth = st.get("cloth")
    if cloth is None or bool(_val("reset", False)) or st.get("shape_key") != shape_key:
        cloth = Cloth(top, side, crv.GetLength(), *shape)
        keep = {k: st[k] for k in ("field", "builder", "builder_v0", "direct",
                                   "cmesh", "ctext", "ckey") if k in st}
        st.clear()                          # keep the same dict: grab handlers use it
        st.update(keep)                     # colliders don't depend on the cloth
        st.update({"cloth": cloth, "shape_key": shape_key})
        store[sid] = st

    # --- physics -------------------------------------------------------------
    w = _val("wind", None)
    g = float(_val("gravity", 0.1))
    accel = (w.X if w else 0.0, w.Y if w else 0.0, (w.Z if w else 0.0) - g)
    stiffness = float(_val("stiffness", 1.0))
    bend = float(_val("bend", 0.3))
    floor_z = _val("floor_z", None)
    gap = max(0.0, float(_val("thickness", drop / 80.0)))
    total_gap = gap + 0.5 * fabric
    reach_step = 0.5 * cloth.dv             # max move per step near colliders

    # --- colliders -----------------------------------------------------------
    geos = _val("colliders", None)
    if geos is not None and not isinstance(geos, (list, tuple)) and \
            (isinstance(geos, rg.GeometryBase) or not hasattr(geos, "__iter__")):
        geos = [geos]                       # item access instead of list
    ckey = _geo_key(geos)
    cnote = ""
    if st.get("ckey") != ckey:
        st["ckey"] = ckey
        cmesh, st["ctext"] = _collider_mesh(geos)
        st["cmesh"] = cmesh
        if cmesh is None:
            st["field"] = st["builder"] = st["direct"] = None
        else:
            V1 = _mesh_vertices(cmesh)
            field = st.get("field")
            moved = False
            if field is not None and getattr(field, "V0", None) is not None \
                    and len(field.V0) == len(V1):
                R, t, err = rigid_fit(field.V0, V1)
                diag = float(np.linalg.norm(V1.max(0) - V1.min(0)))
                if err < 1e-4 * max(diag, 1e-9):
                    field.R, field.t = R, t     # just moved/rotated: reuse grid
                    moved = True
                    cnote = " | collider moved: grid reused"
            if not moved:
                # changed shape: use Rhino directly now, build the grid meanwhile
                st["field"] = None
                st["direct"] = _make_direct(cmesh, reach_step * 8)
                diag = cmesh.GetBoundingBox(True).Diagonal.Length
                st["builder"] = _field_builder(cmesh, max(0.1 * diag, 2.0 * (total_gap + 4 * cloth.dv)))
                st["builder_v0"] = V1
    field = st.get("field")
    if st.get("cmesh") is not None and field is not None and st.get("builder") is None \
            and field.pad < total_gap + 4 * cloth.dv:
        # gap grew past the grid's margin: rebuild it in the background
        cmesh = st["cmesh"]
        st["direct"] = _make_direct(cmesh, reach_step * 8)
        diag = cmesh.GetBoundingBox(True).Diagonal.Length
        V1 = _mesh_vertices(cmesh)
        st["builder"] = _field_builder(cmesh, max(0.1 * diag, 2.0 * (total_gap + 4 * cloth.dv)))
        st["builder_v0"] = V1
        st["field"] = field = None
    if st.get("builder") is not None:
        t_end = time.time() + 0.015          # a little grid building per frame
        try:
            while time.time() < t_end:
                st["build_progress"] = next(st["builder"])
        except StopIteration as done:
            field = done.value
            field.V0 = st.pop("builder_v0", None)
            st["field"], st["builder"], st["direct"] = field, None, None
    collide = field if field is not None else st.get("direct")
    if st.get("builder") is not None:
        cnote = " | building collider grid {:.0%} (direct mode meanwhile)".format(
            st.get("build_progress", 0.0))
    # speed limit: keeps the fall from folding the fabric through itself and
    # points from jumping through thin colliders between checks
    max_speed = reach_step

    phys_key = (accel, stiffness, bend, floor_z, gap, ckey)
    if st.get("phys_key") != phys_key:
        st["phys_key"] = phys_key
        cloth.total = 0
        cloth.motion = float("inf")

    warn = ""
    if field is not None:
        if "inside" not in st and cloth.age <= 10:      # only a fresh cloth
            phi0 = field.sample((field.R.T @ (cloth.P.reshape(3, -1) - field.t[:, None])).T)[0]
            st["inside"] = int((phi0 < 0).sum())
        if st.get("inside"):
            warn = (" | WARNING: fabric started inside the colliders at {} points, "
                    "it may get pierced. Move the rail away from the frame").format(st["inside"])

    # --- grabbing: mouse and anchor points ------------------------------------
    grab_on = bool(_val("grab", False))
    gnote = ""
    try:
        cb = _mouse_callback()
        if grab_on:
            _install_grab(comp, st)
        cb.Enabled = grab_on
        if grab_on:
            s_ = _grab_stats()
            gnote = (" | grab ON: left-drag the fabric (mouse events {}, clicks {}, "
                     "hits {}{})").format(s_["events"], s_["downs"], s_["hits"],
                                          ", error " + s_["error"] if s_["error"] else "")
    except Exception as ex:
        grab_on = False
        gnote = " | grab unavailable: {} (use the anchors input instead)".format(ex)
    if not grab_on:
        st.pop("grab", None)
    grabs = _anchor_grabs(st, cloth, _val("anchors", None))
    anchor_key = tuple(tuple(np.round(p, 6)) for _, p in grabs)
    if st.get("anchor_key") != anchor_key:
        st["anchor_key"] = anchor_key
        cloth.total, cloth.motion = 0, float("inf")
    held = st.get("grab")
    if held:
        grabs = grabs + [(held["idx"], held["target"])]

    # --- simulate --------------------------------------------------------------
    run = bool(_val("run", True))
    max_it = max(1, int(_val("iterations", 1000)))

    def is_settled():
        return held is None and (cloth.total >= max_it or cloth.motion < 1e-4 * cloth.dv)

    settled = is_settled()
    t0 = time.time()
    if run and not settled:
        cloth.step(max(1, int(_val("steps", 5))), accel, stiffness, bend,
                   fabric=fabric, gap=gap, floor_z=floor_z, collide=collide,
                   max_speed=max_speed, grabs=grabs)
        settled = is_settled()
    ms = (time.time() - t0) * 1000.0

    # --- outputs ---------------------------------------------------------------
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
        if fabric > 0:
            try:
                inner = m.Offset(-0.5 * fabric)
                thick = inner.Offset(fabric, True) if inner else None
            except Exception:
                thick = None

    state = "settled" if settled else ("running" if run else "paused")
    if held is not None:
        state = "grabbing"
    txt = ("{} | iteration {} | {:.0f} ms/frame | {} | points {}x{} | folds {} | "
           "fold depth {:.1f} | gap {:g} | max stretch {:.1%} | {}{}{}{}").format(
        state, cloth.total, ms, mode, cloth.nu, cloth.nv, cloth.pleats, cloth.amp,
        gap, cloth.max_strain(), st.get("ctext", "colliders: none"), cnote, gnote, warn)
    comp.Message = state
    st["running"] = (run and not settled) or st.get("builder") is not None
    if st["running"]:
        _schedule_next(comp, 1)
    return m, thick, srf, pts, txt


if IN_RHINO:
    mesh, thick_mesh, surface, points, info = run_grasshopper()
