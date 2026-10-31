# r: numpy
"""
Drape Simulation (live) - Rhino 8 "Python 3 Script" component

A cloth hangs from a top rail curve. The fabric is wider than the rail
(fullness > 1), so it gathers into vertical folds, falls under gravity,
wraps over collider objects and can pool on a floor.

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
    rail        Curve         item   Top edge the fabric hangs from (default: 6 m line at z=4)
    drop        float         item   Hanging length                          (default 4.0)
    fullness    float         item   Fabric width / rail length, >1 = folds  (default 1.6)
    pleats      int           item   Number of folds along the rail          (default 10)
    pleat_depth float         item   Fold depth, 0 = auto from fullness      (default 0)
    variation   float         item   0-1, randomness of fold depth           (default 0.3)
    res_u       int           item   Points across                           (default 61)
    res_v       int           item   Points down                             (default 40)
    steps       int           item   Iterations per frame                    (default 5)
    iterations  int           item   Max iterations after a change           (default 400)
    gravity     float         item   Pull per step                           (default 0.01)
    stiffness   float         item   0-1, stretch resistance                 (default 1.0)
    bend        float         item   0-1, bending resistance                 (default 0.1)
    wind        Vector3d      item   Optional constant push, e.g. (0, 0.003, 0)
    colliders   GeometryBase  list   Optional Mesh / Brep / SubD / Surface obstacles
    thickness   float         item   Gap kept between fabric and colliders   (default 0.05)
                                     (also the max move per step near colliders)
    floor_z     float         item   Optional floor height for pooling fabric
    seed        int           item   Random seed                             (default 1)

OUTPUTS
    mesh        Quad mesh of the fabric, updated every frame
    surface     NURBS surface, built only once the cloth has settled or is paused
    points      Simulated points
    info        Status text
"""
import math
import time
import numpy as np

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
        sets.append(((slice(r, nu - 1, 2), A), (slice(r + 1, nu, 2), A), du, "stretch"))
        sets.append(((A, slice(r, nv - 1, 2)), (A, slice(r + 1, nv, 2)), dv, "stretch"))
    for r in (0, 1):                                  # shear diagonals
        sets.append(((slice(r, nu - 1, 2), slice(0, nv - 1)),
                     (slice(r + 1, nu, 2), slice(1, nv)), diag, "shear"))
        sets.append(((slice(r, nu - 1, 2), slice(1, nv)),
                     (slice(r + 1, nu, 2), slice(0, nv - 1)), diag, "shear"))
    for r in range(4):                                # bending (skip one point)
        sets.append(((slice(r, nu - 2, 4), A), (slice(r + 2, nu, 4), A), 2 * du, "bend"))
        sets.append(((A, slice(r, nv - 2, 4)), (A, slice(r + 2, nv, 4)), 2 * dv, "bend"))
    return sets


class Cloth(object):
    def __init__(self, top, side, rail_len, drop, fullness, pleats,
                 pleat_depth, variation, res_v, seed):
        """top, side: (nu, 3) arrays of rail points and sideways unit normals."""
        top = np.asarray(top, float)
        side = np.asarray(side, float)
        rnd = np.random.default_rng(seed)
        nu, nv = len(top), int(res_v)
        pleats = max(0, min(int(pleats), max(1, (nu - 1) // 4)))  # >= 4 points per fold

        # smooth random variation of fold depth along the rail
        t = np.linspace(0.0, 1.0, nu)
        ph1, ph2 = rnd.uniform(0, 2 * np.pi, 2)
        mult = 1.0 + variation * 0.5 * (np.sin(2 * np.pi * 1.3 * t + ph1) +
                                        np.sin(2 * np.pi * 2.7 * t + ph2))
        if pleats == 0:
            amp = 0.0
        elif pleat_depth is None or pleat_depth <= 0:
            amp = auto_amplitude(rail_len, fullness, nu, pleats, mult)
        else:
            amp = float(pleat_depth)

        du = rail_len * fullness / (nu - 1)
        dv = drop / (nv - 1)
        base = top + side * _fold_offsets(amp, t, pleats, mult)[:, None]
        jitter = rnd.uniform(-0.02, 0.02, (nu, nv)) * dv
        jitter[:, 0] = 0.0
        P = np.repeat(base[:, None, :], nv, axis=1)
        P += side[:, None, :] * jitter[:, :, None]
        P[:, :, 2] -= np.arange(nv) * dv

        self.P = P                        # (nu, nv, 3) positions
        self.prev = P.copy()              # previous positions (Verlet)
        self.W = np.ones((nu, nv, 1))     # inverse mass, 0 = pinned
        self.W[:, 0] = 0.0                # top row hangs from the rail
        self.tether = (np.arange(nv) * dv)[None, :, None] * 1.0001
        self.sets = _constraint_sets(nu, nv, du, dv)
        self.nu, self.nv, self.du, self.dv = nu, nv, du, dv
        self.amp, self.pleats = amp, pleats
        self.total = 0                    # iterations since the last change
        self.motion = float("inf")        # largest move in the last step

    def _project(self, sa, sb, rest, k):
        a, b = self.P[sa], self.P[sb]     # views: edits write into self.P
        wa, wb = self.W[sa], self.W[sb]
        d = b - a
        L = np.sqrt((d * d).sum(-1, keepdims=True))
        ws = wa + wb
        s = k * (L - rest) / (np.maximum(L, 1e-12) * np.where(ws > 0, ws, 1.0))
        corr = d * s
        a += wa * corr
        b -= wb * corr

    def step(self, n, accel, stiffness, bend, floor_z=None, collide=None,
             damping=0.98, passes=4, max_speed=None, grab=None):
        """grab: ((i, j), target_xyz) holds one point at the cursor."""
        P, W = self.P, self.W
        accel = np.asarray(accel, float)
        k_of = {"stretch": stiffness, "shear": stiffness * 0.5, "bend": bend}
        held_w = None
        if grab is not None:
            (gi, gj), target = grab
            held_w = W[gi, gj].copy()
            W[gi, gj] = 0.0                       # held points act pinned
            self.prev[gi, gj] = P[gi, gj]         # so it flings when released
            P[gi, gj] = target
        free = W > 0
        try:
            for _ in range(n):
                before = P.copy()
                vel = (P - self.prev) * damping + accel
                if max_speed:
                    # never move further than max_speed per step, so points
                    # can't jump through a thin collider between checks
                    sp = np.sqrt((vel * vel).sum(-1, keepdims=True))
                    vel *= np.minimum(1.0, max_speed / np.maximum(sp, 1e-12))
                self.prev[:] = P
                P += vel * W
                for _ in range(passes):
                    for sa, sb, rest, kind in self.sets:
                        k = k_of[kind]
                        if k > 0:
                            self._project(sa, sb, rest, k)
                # tethers: no point gets further from its top point than the
                # fabric length between them (stops rubbery stretching)
                d = P - P[:, :1]
                L = np.sqrt((d * d).sum(-1, keepdims=True))
                over = (L > self.tether) & free
                if over.any():
                    P[:] = np.where(over, P[:, :1] + d * (self.tether / np.maximum(L, 1e-12)), P)
                if floor_z is not None:
                    low = (P[:, :, 2:3] < floor_z) & free
                    if low.any():
                        P[:, :, 2:3] = np.where(low, floor_z, P[:, :, 2:3])
                        self.prev[:] = np.where(low, P, self.prev)   # friction
                if collide is not None:
                    collide(P, self.prev, free)
                self.total += 1
                self.motion = float(np.abs(P - before).max())
        finally:
            if held_w is not None:
                W[gi, gj] = held_w

    def nearest_on_screen(self, M, x, y, radius):
        """Index of the cloth point closest to screen pixel (x, y), using the
        4x4 world->screen matrix M. None if nothing within radius pixels."""
        pts = self.P.reshape(-1, 3)
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

    def max_strain(self):
        """Largest stretch (tension) of the stretch springs, as a ratio."""
        worst = 0.0
        for sa, sb, rest, kind in self.sets:
            if kind != "stretch":
                continue
            d = self.P[sb] - self.P[sa]
            L = np.sqrt((d * d).sum(-1))
            live = (self.W[sa] + self.W[sb])[..., 0] > 0
            if live.any():
                worst = max(worst, float(((L[live] - rest) / rest).max()))
        return worst

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
            mp = rg.MeshingParameters.QualityRenderMesh
            part = rg.Mesh()
            for m in rg.Mesh.CreateFromBrep(g, mp) or []:
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


def _make_collider(cmesh, thickness, max_speed, refresh=8):
    """Pushes points out of the collider mesh, keeping a `thickness` gap.

    Broad phase: every `refresh` steps NumPy picks the points inside the
    collider's bounding box and Rhino measures how far they are from the
    surface. Only points close enough to touch it before the next refresh
    (they can move at most max_speed per step) are tested every step."""
    bb = cmesh.GetBoundingBox(True)
    pad = thickness + max_speed * refresh
    lo = np.array([bb.Min.X, bb.Min.Y, bb.Min.Z]) - pad
    hi = np.array([bb.Max.X, bb.Max.Y, bb.Max.Z]) + pad
    state = {"count": 0, "near": None}

    def signed_push(P, prev, i, j):
        x, y, z = P[i, j]
        # no distance limit: points that ended up deep inside still find
        # their way out (a limited search returns nothing for them)
        mp = cmesh.ClosestMeshPoint(rg.Point3d(x, y, z), 0.0)
        if mp is None:
            return float("inf")
        cp = mp.Point
        n = cmesh.NormalAt(mp)
        n.Unitize()
        dist = (x - cp.X) * n.X + (y - cp.Y) * n.Y + (z - cp.Z) * n.Z
        if dist < thickness:
            P[i, j] = (cp.X + n.X * thickness, cp.Y + n.Y * thickness,
                       cp.Z + n.Z * thickness)
            prev[i, j] = P[i, j]                      # friction
        return dist

    def collide(P, prev, free):
        if state["count"] % refresh == 0 or state["near"] is None:
            near = np.zeros(P.shape[:2], bool)
            inside_box = np.all((P >= lo) & (P <= hi), axis=-1) & free[..., 0]
            for i, j in zip(*np.nonzero(inside_box)):
                if signed_push(P, prev, i, j) < pad:
                    near[i, j] = True
            state["near"] = near
        else:
            for i, j in zip(*np.nonzero(state["near"] & free[..., 0])):
                signed_push(P, prev, i, j)
        state["count"] += 1
    return collide


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


def _quad_mesh(cloth):
    nu, nv = cloth.nu, cloth.nv
    m = rg.Mesh()
    for x, y, z in cloth.P.reshape(-1, 3).tolist():
        m.Vertices.Add(x, y, z)
    for i in range(nu - 1):
        for j in range(nv - 1):
            a, b = i * nv + j, (i + 1) * nv + j
            m.Faces.AddFace(a, b, b + 1, a + 1)
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
        x, y, z = cloth.P[hit]
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
    store = sc.sticky.setdefault("drape_sim", {})
    sid = str(comp.InstanceGuid)
    st = store.get(sid) or {}

    crv = _val("rail", None)
    if crv is None:
        crv = rg.LineCurve(rg.Point3d(0, 0, 4), rg.Point3d(6, 0, 4))
    nu = max(5, int(_val("res_u", 61)))
    nv = max(3, int(_val("res_v", 40)))
    top, side = _rail_frames(crv, nu)
    shape = (float(_val("drop", 4.0)), float(_val("fullness", 1.6)),
             int(_val("pleats", 10)), float(_val("pleat_depth", 0.0)),
             float(_val("variation", 0.3)), nv, int(_val("seed", 1)))
    shape_key = (tuple(round(v, 6) for p in top for v in p), shape)

    cloth = st.get("cloth")
    if cloth is None or bool(_val("reset", False)) or st.get("shape_key") != shape_key:
        cloth = Cloth(top, side, crv.GetLength(), *shape)
        st.clear()                          # keep the same dict: grab handlers use it
        st.update({"cloth": cloth, "shape_key": shape_key})
        store[sid] = st

    # physics: applied live, any change wakes the cloth up again
    w = _val("wind", None)
    g = float(_val("gravity", 0.01))
    accel = (w.X if w else 0.0, w.Y if w else 0.0, (w.Z if w else 0.0) - g)
    stiffness = float(_val("stiffness", 1.0))
    bend = float(_val("bend", 0.1))
    floor_z = _val("floor_z", None)
    thickness = max(1e-4, float(_val("thickness", 0.05)))
    geos = _val("colliders", None)
    if geos is not None and not isinstance(geos, (list, tuple)) and \
            (isinstance(geos, rg.GeometryBase) or not hasattr(geos, "__iter__")):
        geos = [geos]                       # item access instead of list
    max_speed = None
    ckey = (_geo_key(geos), thickness)
    if st.get("ckey") != ckey:
        cmesh, st["ctext"] = _collider_mesh(geos)
        st["cmesh"] = cmesh
        st["ckey"] = ckey
        st["collide"] = None
        if cmesh is not None:
            st["collide"] = _make_collider(cmesh, thickness, max(thickness, 0.5 * cloth.dv))
    if st.get("cmesh") is not None:
        # with colliders, points move at most max_speed per step (no tunnelling)
        max_speed = max(thickness, 0.5 * cloth.dv)
    phys_key = (accel, stiffness, bend, floor_z, ckey)
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
        return held is None and (cloth.total >= max_it or cloth.motion < 1e-5 * cloth.dv)
    settled = is_settled()
    t0 = time.time()
    if run and not settled:
        cloth.step(max(1, int(_val("steps", 5))), accel, stiffness, bend,
                   floor_z=floor_z, collide=st.get("collide"),
                   max_speed=max_speed, grab=grab)
        settled = is_settled()
    ms = (time.time() - t0) * 1000.0

    m = _quad_mesh(cloth)
    m.Normals.ComputeNormals()
    pts = [rg.Point3d(x, y, z) for x, y, z in cloth.P.reshape(-1, 3).tolist()]
    srf = None
    if settled or not run:          # surface fitting is slower: only when still
        srf = rg.NurbsSurface.CreateFromPoints(pts, cloth.nu, cloth.nv, 3, 3)

    state = "settled" if settled else ("running" if run else "paused")
    if held is not None:
        state = "grabbing"
    txt = ("{} | iteration {} | {:.0f} ms/frame | points {}x{} | folds {} | "
           "fold depth {:.3f} | max stretch {:.1%} | {}{}").format(
        state, cloth.total, ms, cloth.nu, cloth.nv, cloth.pleats, cloth.amp,
        cloth.max_strain(), st.get("ctext", "colliders: none"),
        " | grab: drag the cloth with the left mouse button" if grab_on else grab_err)
    comp.Message = state
    st["running"] = run and not settled
    if st["running"]:
        _schedule_next(comp, 1)
    return m, srf, pts, txt


if IN_RHINO:
    mesh, surface, points, info = run_grasshopper()
