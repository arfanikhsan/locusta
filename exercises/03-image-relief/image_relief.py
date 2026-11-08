# r: numpy, pillow
"""
Image Relief - Rhino 8 "Python 3 Script" component

Turns an image into real mesh geometry on a flat or curved surface: every
pixel's brightness moves the surface along its normal (white = out, black =
stays). No normal or displacement map at render time, the relief IS the
mesh, so you can bake it, 3D print it, mill it or use it as a facade panel.

Example: a pillar with a water texture. Draw a cylinder (or connect nothing
and get a default pillar), plug a picture of water into `image`, done.

BASE: Surfaces, Breps and Extrusions. Each untrimmed face is textured on its
own (a capped cylinder: the side gets the texture, the trimmed caps are
skipped). Image left->right follows the surface u direction, bottom->top
follows v. Use `swap` if the image runs the wrong way.

CLOSED SURFACES (pillar, tube): the image repeats a whole number of times
around, and the mesh is welded at the seam. If the image's left and right
edges don't match, turn on `mirror` for a seamless pattern.

INPUTS  (right-click each input -> set Access and Type hint)
    image       str           item   Path to a .png / .jpg / .bmp / .tif (File Path component)
    base        GeometryBase  list   Surface / Brep / Extrusion   (default: pillar r=25, h=300)
    depth       float         item   Relief depth for pure white in model units. Negative
                                     carves into the surface   (default: 3% of the surface size)
    tile_size   float         item   Width of one copy of the image in model units.
                                     0 = the image fits once across u, keeps its proportions (default 0)
    detail      int           item   Mesh vertices along the longest side   (default 300)
    blur        float         item   Extra smoothing in image pixels, removes photo noise (default 1)
    gamma       float         item   >1 sharper peaks, <1 flatter, rounder tops   (default 1)
    invert      bool          item   Dark = out instead of light = out
    mirror      bool          item   Mirror every other copy -> seamless tiling
    swap        bool          item   Run the image along v instead of u
    offset_u    float         item   Slide the image, in copies (0.5 = half an image)
    offset_v    float         item   Slide the image, in copies
    fade        float         item   Relief fades to 0 over this distance from open edges,
                                     so it meets the rest of the model cleanly   (default 0)
    flip        bool          item   Push to the other side of the surface
    thickness   float         item   >0: closed solid shell this thick behind the surface
    cap         bool          item   Close the ends of a pillar / tube with flat caps

OUTPUTS
    mesh        The relief mesh (a closed solid with thickness or cap)
    heights     Relief height per vertex, 0-1 (for Gradient colouring)
    info        Status text
"""
import math
import time
import numpy as np


# ----------------------------------------------------------------------------
# core (NumPy only, no Rhino needed)
# ----------------------------------------------------------------------------

LUMA = np.array([0.299, 0.587, 0.114])


def to_height(rgb):
    """(h, w, 3 or 4) uint8/float image -> (h, w) floats 0-1. Transparent
    pixels count as black."""
    a = np.asarray(rgb, dtype=float)
    if a.max() > 1.0:
        a = a / 255.0
    if a.ndim == 2:
        return a
    h = a[..., :3] @ LUMA
    if a.shape[2] == 4:
        h = h * a[..., 3]
    return h


def _pad_mode(mode):
    return {"repeat": "wrap", "mirror": "symmetric"}.get(mode, "edge")


def _box(a, r, axis, mode):
    """Box blur of radius r (int) along one axis, with the edge rule `mode`."""
    if r < 1:
        return a
    pad = [(0, 0)] * a.ndim
    pad[axis] = (r + 1, r)
    p = np.pad(a, pad, mode=_pad_mode(mode))
    c = np.cumsum(p, axis=axis)
    n = a.shape[axis]
    hi = np.take(c, np.arange(2 * r + 1, 2 * r + 1 + n), axis=axis)
    lo = np.take(c, np.arange(0, n), axis=axis)
    return (hi - lo) / (2 * r + 1)


def smooth_image(a, sigma_x, sigma_y=None, mode="repeat"):
    """Gaussian-like blur (3 box passes), sigmas in pixels."""
    if sigma_y is None:
        sigma_y = sigma_x
    for axis, s in ((1, sigma_x), (0, sigma_y)):
        if s <= 0.3:
            continue
        # three boxes of radius r have the variance of a gaussian with sigma s
        r = max(1, int(round((math.sqrt(4.0 * s * s + 1.0) - 1.0) / 2.0)))
        r = min(r, a.shape[axis] // 2)
        for _ in range(3):
            a = _box(a, r, axis, mode)
    return a


def normalize(h, clip=0.5):
    """Stretch to 0-1. The darkest and brightest `clip` percent are cut off,
    so a single glint (sun on water) doesn't flatten the rest of the relief."""
    lo, hi = np.percentile(h, [clip, 100.0 - clip]) if clip > 0 else (h.min(), h.max())
    lo, hi = float(lo), float(hi)
    if hi - lo < 1e-9:
        return np.zeros_like(h)
    return np.clip((h - lo) / (hi - lo), 0.0, 1.0)


def downsample(H, w, h):
    """Area-average H down to at most (h, w) pixels. Cheap anti-aliasing:
    the mesh never needs more than about 2 pixels per vertex."""
    hh, ww = H.shape
    w, h = int(min(max(w, 2), ww)), int(min(max(h, 2), hh))
    if w < ww:
        e = np.linspace(0, ww, w + 1).astype(int)
        H = np.add.reduceat(H, e[:-1], axis=1) / np.diff(e)
    if h < hh:
        e = np.linspace(0, hh, h + 1).astype(int)
        H = np.add.reduceat(H, e[:-1], axis=0) / np.diff(e)[:, None]
    return H


def bspline_basis(knots, p, t):
    """B-spline basis functions of degree p and their first derivatives at
    the parameters t, as dense (len(t), n_cv) arrays. `knots` is the full
    knot vector (n_cv + p + 1 values)."""
    k = np.asarray(knots, float)
    t = np.asarray(t, float)
    n_cv = len(k) - p - 1
    span = np.clip(np.searchsorted(k, t, side="right") - 1, p, n_cv - 1)
    N = np.zeros((len(t), len(k) - 1))
    N[np.arange(len(t)), span] = 1.0
    T = t[:, None]
    dN = np.zeros((len(t), n_cv))
    for d in range(1, p + 1):
        d1 = k[d:-1] - k[:-d - 1]
        d2 = k[d + 1:] - k[1:-d]
        s1 = np.where(d1 > 0, d1, 1.0)
        s2 = np.where(d2 > 0, d2, 1.0)
        if d == p:
            dN = d * (np.where(d1 > 0, N[:, :-1] / s1, 0.0)
                      - np.where(d2 > 0, N[:, 1:] / s2, 0.0))
        a = np.where(d1 > 0, (T - k[:-d - 1]) / s1, 0.0)
        b = np.where(d2 > 0, (k[d + 1:] - T) / s2, 0.0)
        N = a * N[:, :-1] + b * N[:, 1:]
    return N, dN


def nurbs_grid(C, W, ku, kv, pu, pv, us, vs):
    """Points and unnormalised normals (Su x Sv) of a NURBS surface on the
    grid us x vs, all at once. C (cu, cv, 3) control points, W (cu, cv)
    weights. Returns (nv, nu, 3) arrays."""
    Bu, dBu = bspline_basis(ku, pu, us)
    Bv, dBv = bspline_basis(kv, pv, vs)
    Cw = C * W[..., None]

    def ev(bv, bu, X):                      # (nv, nu) for each coordinate
        return np.einsum("jb,abc,ia->jic", bv, X, bu, optimize=True)
    A, w = ev(Bv, Bu, Cw), ev(Bv, Bu, W[..., None])
    Au, wu = ev(Bv, dBu, Cw), ev(Bv, dBu, W[..., None])
    Av, wv = ev(dBv, Bu, Cw), ev(dBv, Bu, W[..., None])
    P = A / w
    Su = (Au - wu * P) / w
    Sv = (Av - wv * P) / w
    return P, np.cross(Su, Sv)


def _wrap(t, mode):
    """Texture coordinate (any real number, 1 = one image) -> 0..1."""
    if mode == "repeat":
        return t - np.floor(t)
    if mode == "mirror":
        t = np.mod(t, 2.0)
        return np.where(t > 1.0, 2.0 - t, t)
    return np.clip(t, 0.0, 1.0)


def sample(H, tu, tv, mode="repeat"):
    """Bilinear lookup. tu runs left->right, tv bottom->top, 1 = one image."""
    h, w = H.shape
    x = _wrap(np.asarray(tu, float), mode) * w - 0.5
    y = (1.0 - _wrap(np.asarray(tv, float), mode)) * h - 0.5
    x0, y0 = np.floor(x), np.floor(y)
    fx, fy = x - x0, y - y0
    x0, y0 = x0.astype(int), y0.astype(int)

    def idx(i, n):
        if mode == "repeat":
            return np.mod(i, n)
        if mode == "mirror":            # pixel centres mirror around the edge
            i = np.mod(i, 2 * n)
            return np.where(i >= n, 2 * n - 1 - i, i)
        return np.clip(i, 0, n - 1)

    xa, xb = idx(x0, w), idx(x0 + 1, w)
    ya, yb = idx(y0, h), idx(y0 + 1, h)
    top = H[ya, xa] * (1 - fx) + H[ya, xb] * fx
    bot = H[yb, xa] * (1 - fx) + H[yb, xb] * fx
    return top * (1 - fy) + bot * fy


def _arc(P, axis, closed):
    """Cumulative distance along `axis` of the grid P (nv, nu, 3), and the
    total length of each line. With `closed`, the wrap-around segment counts."""
    d = np.diff(P, axis=axis)
    if closed:
        first = np.take(P, [0], axis=axis)
        last = np.take(P, [-1], axis=axis)
        d = np.concatenate([d, first - last], axis=axis)
    seg = np.linalg.norm(d, axis=2)
    cum = np.cumsum(seg, axis=axis)
    total = np.take(cum, [-1], axis=axis)
    shape = list(cum.shape)
    shape[axis] = 1
    cum = np.concatenate([np.zeros(shape), cum], axis=axis)
    if closed:                                     # drop the wrap point again
        cum = np.take(cum, np.arange(P.shape[axis]), axis=axis)
    return cum, total


def surface_lengths(P, closed_u, closed_v):
    """Average length of the grid along u and along v."""
    _, tu = _arc(P, 1, closed_u)
    _, tv = _arc(P, 0, closed_v)
    return float(tu.mean()), float(tv.mean())


def plan_repeats(len_u, len_v, img_w, img_h, tile_size, closed_u, closed_v, mirror):
    """How many copies of the image go across u and along v. Closed
    directions get a whole number (an even one when mirrored) so the seam
    matches."""
    aspect = img_h / float(img_w)
    tile = tile_size if tile_size > 0 else len_u
    ru = len_u / tile
    rv = len_v / (tile * aspect)

    def whole(r, closed):
        if not closed:
            return r
        step = 2 if mirror else 1
        return max(step, step * round(r / step))
    return whole(ru, closed_u), whole(rv, closed_v)


def tex_coords(P, closed_u, closed_v, rep_u, rep_v, off_u=0.0, off_v=0.0):
    """Texture coordinates by arc length: the image isn't stretched where the
    surface parameterisation is uneven."""
    cu, tu = _arc(P, 1, closed_u)
    cv, tv = _arc(P, 0, closed_v)
    TU = cu / np.maximum(tu, 1e-12) * rep_u + off_u
    TV = cv / np.maximum(tv, 1e-12) * rep_v + off_v
    return TU, TV, cu, tu, cv, tv


def edge_fade(cu, tu, cv, tv, closed_u, closed_v, fade):
    """1 inside, smoothly down to 0 at open edges within `fade` distance."""
    f = np.ones(cu.shape)
    if fade <= 0:
        return f

    def ramp(d):
        t = np.clip(d / fade, 0.0, 1.0)
        return t * t * (3.0 - 2.0 * t)
    if not closed_u:
        f *= ramp(np.minimum(cu, tu - cu))
    if not closed_v:
        f *= ramp(np.minimum(cv, tv - cv))
    return f


def grid_faces(nv, nu, closed_u, closed_v, start=0):
    """Quads (a, b, c, d) for a (nv, nu) vertex grid, normal along Su x Sv."""
    iu = np.arange(nu if closed_u else nu - 1)
    jv = np.arange(nv if closed_v else nv - 1)
    J, I = np.meshgrid(jv, iu, indexing="ij")
    I1, J1 = (I + 1) % nu, (J + 1) % nv
    q = np.stack([J * nu + I, J * nu + I1, J1 * nu + I1, J1 * nu + I], -1)
    return q.reshape(-1, 4) + start


def _boundary_rows(nv, nu, closed_u, closed_v):
    """The open edges of the grid as directed vertex-index chains, running
    the same way as the front faces' edges on that border."""
    chains = []
    idx = np.arange(nv * nu).reshape(nv, nu)
    if not closed_v:
        bottom = idx[0]                   # faces use bottom edges a->b (along +u)
        top = idx[-1][::-1]               # and top edges c->d (along -u)
        if closed_u:
            bottom = np.append(bottom, bottom[0])
            top = np.append(top, top[0])
        chains += [bottom, top]
    if not closed_u:
        right = idx[:, -1]                # b->c, along +v
        left = idx[:, 0][::-1]            # d->a, along -v
        if closed_v:
            right = np.append(right, right[0])
            left = np.append(left, left[0])
        chains += [right, left]
    return chains


def build_mesh(front, closed_u, closed_v, back=None, cap=False):
    """Vertex array (n, 3) and face array (m, 4) (triangles repeat the last
    index). `back` (same shape as front) closes it into a shell; `cap` closes
    the open ends of a tube with a fan around each end's centre."""
    nv, nu = front.shape[:2]
    n = nv * nu
    V = [front.reshape(-1, 3)]
    F = [grid_faces(nv, nu, closed_u, closed_v)]
    chains = _boundary_rows(nv, nu, closed_u, closed_v)
    if back is not None:
        V.append(back.reshape(-1, 3))
        F.append(grid_faces(nv, nu, closed_u, closed_v, n)[:, ::-1])
        for c in chains:                  # side walls: front edge a->b, back b'->a'
            a, b = c[:-1], c[1:]
            F.append(np.stack([b, a, a + n, b + n], -1))
    elif cap and closed_u != closed_v:
        start = n
        for c in chains:
            ring = c[:-1]
            V.append(front.reshape(-1, 3)[ring].mean(0)[None])
            a, b = c[:-1], c[1:]
            ctr = np.full_like(a, start)
            F.append(np.stack([b, a, ctr, ctr], -1))
            start += 1
    return np.concatenate(V), np.concatenate(F)


def relief(P, N, H, depth, rep_u, rep_v, closed_u, closed_v, mode="repeat",
           off_u=0.0, off_v=0.0, fade=0.0, gamma=1.0, aa=False):
    """Moves the base grid P (nv, nu, 3) along its normals N by the image H.
    Returns (displaced points, heights 0-1)."""
    TU, TV, cu, tu, cv, tv = tex_coords(P, closed_u, closed_v, rep_u, rep_v, off_u, off_v)
    if aa:                                 # don't let pixels finer than the mesh alias
        nv, nu = P.shape[:2]
        sx = H.shape[1] * rep_u / max(nu - 1, 1)
        sy = H.shape[0] * rep_v / max(nv - 1, 1)
        H = smooth_image(H, 0.5 * sx if sx > 1 else 0, 0.5 * sy if sy > 1 else 0, mode)
    h = sample(H, TU, TV, mode)
    if gamma != 1.0:
        h = np.clip(h, 0.0, 1.0) ** gamma
    h = h * edge_fade(cu, tu, cv, tv, closed_u, closed_v, fade)
    return P + N * (depth * h)[..., None], h


def mesh_is_closed(V, F):
    """True when every edge is shared by exactly two faces running opposite
    ways (a watertight, consistently oriented mesh)."""
    edges = {}
    for f in F.tolist():
        f = [f[0], f[1], f[2]] if f[2] == f[3] else f
        for a, b in zip(f, f[1:] + f[:1]):
            edges[(a, b)] = edges.get((a, b), 0) + 1
    return all(c == 1 and edges.get((b, a)) == 1 for (a, b), c in edges.items())


# ----------------------------------------------------------------------------
# grasshopper
# ----------------------------------------------------------------------------
try:
    import Rhino
    import Rhino.Geometry as rg
    import scriptcontext as sc
    IN_RHINO = True
except ImportError:
    IN_RHINO = False

DIAG = {}                  # which fast paths worked, and why not (for `info`)
MAX_PX = 2048              # bigger images are scaled down (enough for any mesh)
MAX_VERTS = 1500000


def _val(name, default):
    """A component input, or `default` when it's missing or not connected.
    (Functions and modules are skipped, in case a name ever clashes.)"""
    v = globals().get(name)
    if v is None or callable(v) or type(v).__name__ == "module":
        return default
    return v


def load_image(path):
    """Image file -> (h, w) heights 0-1. Pillow first, System.Drawing next."""
    try:
        from PIL import Image
        im = Image.open(path)
        im.load()
        if max(im.size) > MAX_PX:
            im.thumbnail((MAX_PX, MAX_PX))
        im = im.convert("RGBA")
        return to_height(np.asarray(im)), "Pillow"
    except ImportError:
        pass
    import base64
    import System
    import System.Drawing as sd
    from System.Runtime.InteropServices import Marshal
    bmp = sd.Bitmap(path)
    w, h = bmp.Width, bmp.Height
    if max(w, h) > MAX_PX:
        s = MAX_PX / float(max(w, h))
        bmp = sd.Bitmap(bmp, int(w * s), int(h * s))
        w, h = bmp.Width, bmp.Height
    rect = sd.Rectangle(0, 0, w, h)
    data = bmp.LockBits(rect, sd.Imaging.ImageLockMode.ReadOnly,
                        sd.Imaging.PixelFormat.Format32bppArgb)
    try:
        n = abs(data.Stride) * h
        buf = System.Array.CreateInstance(System.Byte, n)
        Marshal.Copy(data.Scan0, buf, 0, n)
        stride = abs(data.Stride)
    finally:
        bmp.UnlockBits(data)
    # base64 is the fastest way to move a big .NET byte array into Python
    raw = np.frombuffer(base64.b64decode(str(System.Convert.ToBase64String(buf))), np.uint8)
    px = raw.reshape(h, stride)[:, :w * 4].reshape(h, w, 4)
    return to_height(px[..., [2, 1, 0, 3]]), "System.Drawing"


def _image(st, path):
    import os
    path = str(path).strip().strip('"')
    if not os.path.isfile(path):
        raise IOError("image not found: " + path)
    key = (path, os.path.getmtime(path))
    if st.get("img_key") != key:
        st["img"], st["img_by"] = load_image(path)
        st["img_key"] = key
        st.pop("proc_key", None)
    return st["img"]


def _to_geometry(g):
    if g is None:
        return None
    if hasattr(g, "Value") and not isinstance(g, rg.GeometryBase):
        g = g.Value
    if not isinstance(g, rg.GeometryBase):
        try:
            obj = Rhino.RhinoDoc.ActiveDoc.Objects.FindId(g)
            g = obj.Geometry if obj else None
        except Exception:
            g = None
    return g


def _surfaces(geos):
    """Untrimmed surfaces to texture, with a flag for reversed Brep faces."""
    out, skipped = [], 0
    for g in geos or []:
        g = _to_geometry(g)
        if isinstance(g, rg.Extrusion):
            g = g.ToBrep()
        if isinstance(g, rg.BrepFace):
            g = g.DuplicateFace(False)
        if isinstance(g, rg.Brep):
            g = g.DuplicateBrep()
            # Boundary Surfaces / PlanarSrf make a rectangle as a trimmed,
            # oversized plane: shrinking fits it to its edges -> untrimmed
            g.Faces.ShrinkFaces()
            for f in g.Faces:
                if f.IsSurface:
                    out.append((f.UnderlyingSurface(), f.OrientationIsReversed))
                else:
                    skipped += 1
        elif isinstance(g, rg.Surface):
            out.append((g, False))
        elif g is not None:
            skipped += 1
    return out, skipped


def _default_base():
    c = rg.Circle(rg.Plane.WorldXY, 25.0)
    return rg.Cylinder(c, 300.0).ToNurbsSurface()


def _nurbs_arrays(ns):
    """Control points, weights and full knot vectors of a NurbsSurface."""
    cu, cv = ns.Points.CountU, ns.Points.CountV
    C = np.zeros((cu, cv, 3))
    W = np.ones((cu, cv))
    for i in range(cu):
        for j in range(cv):
            cp = ns.Points.GetControlPoint(i, j)
            p = cp.Location
            C[i, j] = (p.X, p.Y, p.Z)
            W[i, j] = cp.Weight
    ku = [ns.KnotsU[i] for i in range(ns.KnotsU.Count)]
    kv = [ns.KnotsV[i] for i in range(ns.KnotsV.Count)]
    # Rhino leaves out the first and last knot (they never matter)
    ku = np.array([ku[0]] + ku + [ku[-1]])
    kv = np.array([kv[0]] + kv + [kv[-1]])
    return C, W, ku, kv, ns.OrderU - 1, ns.OrderV - 1


def _evaluate_fast(srf, us, vs):
    """The whole grid in one go with NumPy, checked against Rhino at a few
    points. Returns None if anything doesn't match (then use Rhino)."""
    ns = srf if isinstance(srf, rg.NurbsSurface) else srf.ToNurbsSurface()
    if ns is None:
        return None
    C, W, ku, kv, pu, pv = _nurbs_arrays(ns)
    P, N = nurbs_grid(C, W, ku, kv, pu, pv, us, vs)
    ln = np.linalg.norm(N, axis=2)
    size = max(float(np.ptp(P.reshape(-1, 3), 0).max()), 1e-9)
    nv, nu = ln.shape
    for j, i in ((0, 0), (nv // 3, nu // 2), (nv - 1, nu - 1), (nv // 2, nu // 5)):
        p = srf.PointAt(float(us[i]), float(vs[j]))
        n = srf.NormalAt(float(us[i]), float(vs[j]))
        if np.linalg.norm(P[j, i] - (p.X, p.Y, p.Z)) > 1e-6 * size:
            return None
        if ln[j, i] > 1e-12 * size * size and np.dot(N[j, i] / ln[j, i], (n.X, n.Y, n.Z)) < 0.999:
            return None
    bad = ln <= 1e-9 * ln.max()                    # poles: ask Rhino there
    for j, i in zip(*np.nonzero(bad)):
        n = srf.NormalAt(float(us[i]), float(vs[j]))
        N[j, i] = (n.X, n.Y, n.Z)
        ln[j, i] = 1.0
    return P, N / ln[..., None]


def _evaluate(srf, nu, nv, closed_u, closed_v):
    """Points and unit normals on an (nv, nu) parameter grid."""
    du, dv = srf.Domain(0), srf.Domain(1)
    us = np.linspace(du.T0, du.T1, nu, endpoint=not closed_u)
    vs = np.linspace(dv.T0, dv.T1, nv, endpoint=not closed_v)
    try:
        fast = _evaluate_fast(srf, us, vs)
        if fast is None:
            DIAG["surface"] = "slow (NumPy didn't match Rhino)"
    except Exception as e:
        fast = None
        DIAG["surface"] = "slow (fast path error: {}: {})".format(type(e).__name__, str(e)[:120])
    if fast is not None:
        DIAG["surface"] = "fast"
        return fast
    P = np.zeros((nv, nu, 3))
    N = np.zeros((nv, nu, 3))
    for j, v in enumerate(vs.tolist()):
        for i, u in enumerate(us.tolist()):
            p = srf.PointAt(u, v)
            n = srf.NormalAt(u, v)
            P[j, i] = (p.X, p.Y, p.Z)
            N[j, i] = (n.X, n.Y, n.Z)
    ln = np.linalg.norm(N, axis=2, keepdims=True)
    return P, N / np.maximum(ln, 1e-12)


def _srf_key(srf):
    b = srf.GetBoundingBox(False)
    du, dv = srf.Domain(0), srf.Domain(1)
    pts = []
    for a, c in ((0.31, 0.17), (0.5, 0.5), (0.83, 0.71)):
        p = srf.PointAt(du.ParameterAt(a), dv.ParameterAt(c))
        pts += [p.X, p.Y, p.Z]
    return tuple(round(x, 6) for x in
                 [b.Min.X, b.Min.Y, b.Min.Z, b.Max.X, b.Max.Y, b.Max.Z,
                  du.T0, du.T1, dv.T0, dv.T1] + pts)


def _net_array(net_type, n, data):
    """A .NET array of n structs filled with the raw bytes of a NumPy array
    in one copy (base64 is the fastest way across from Python)."""
    import base64
    import System
    from System.Runtime.InteropServices import Marshal, GCHandle, GCHandleType
    import clr
    arr = System.Array.CreateInstance(clr.GetClrType(net_type), n)
    src = System.Convert.FromBase64String(base64.b64encode(data.tobytes()).decode("ascii"))
    h = GCHandle.Alloc(arr, GCHandleType.Pinned)
    try:
        Marshal.Copy(src, 0, h.AddrOfPinnedObject(), src.Length)
    finally:
        h.Free()
    return arr


def _check_mesh(m, V, F, how):
    """True if the Rhino mesh holds exactly V and F."""
    if m.Vertices.Count != len(V) or m.Faces.Count != len(F):
        DIAG["mesh_" + how] = "gave {} / {} vertices".format(m.Vertices.Count, len(V))
        return False
    for k in (0, len(V) // 2, len(V) - 1):
        p = m.Vertices[k]
        if abs(p.X - V[k, 0]) + abs(p.Y - V[k, 1]) + abs(p.Z - V[k, 2]) > 1e-3 * (1 + abs(V[k]).max()):
            DIAG["mesh_" + how] = "vertex mismatch"
            return False
    f = m.Faces[len(F) - 1]
    if (f.A, f.B, f.C, f.D) != tuple(int(x) for x in F[-1]):
        DIAG["mesh_" + how] = "face mismatch"
        return False
    return True


def _mesh_by_copy(V, F):
    """Fastest: copy the raw bytes into .NET Point3f / MeshFace arrays."""
    m = rg.Mesh()
    m.Vertices.AddVertices(_net_array(rg.Point3f, len(V), np.ascontiguousarray(V, np.float32)))
    m.Faces.AddFaces(_net_array(rg.MeshFace, len(F), np.ascontiguousarray(F, np.int32)))
    return m


def _mesh_by_list(V, F):
    """Fast enough: make the structs in Python, hand them over in one call
    each (the slow part of the old way was one Rhino call per vertex)."""
    import System
    P3f, MF = rg.Point3f, rg.MeshFace
    pts = System.Array[P3f]([P3f(x, y, z) for x, y, z in V.tolist()])
    fcs = System.Array[MF]([MF(a, b, c, d) for a, b, c, d in F.tolist()])
    m = rg.Mesh()
    m.Vertices.AddVertices(pts)
    m.Faces.AddFaces(fcs)
    return m


def _fast_mesh(V, F):
    """Tries the fast ways in turn; each result is checked against V and F.
    A triangle repeats its last index, which is also Rhino's convention."""
    for how, build in (("copy", _mesh_by_copy), ("list", _mesh_by_list)):
        try:
            m = build(V, F)
            if _check_mesh(m, V, F, how):
                m.Normals.ComputeNormals()
                m.FaceNormals.ComputeFaceNormals()
                DIAG["mesh"] = how
                return m
        except Exception as e:
            DIAG["mesh_" + how] = "{}: {}".format(type(e).__name__, str(e)[:100])
    DIAG["mesh"] = "slow (one vertex at a time)"
    return None


def _to_rhino_mesh(V, F, cache, key):
    """Fast path first. Fallback: faces are built once per topology and
    reused, and only the vertices are moved (slow for big meshes)."""
    m = _fast_mesh(V, F)
    if m is not None:
        return m
    m = cache.get(key)
    if m is None:
        m = rg.Mesh()
        for _ in range(len(V)):
            m.Vertices.Add(0.0, 0.0, 0.0)
        for a, b, c, d in F.tolist():
            if c == d:
                m.Faces.AddFace(a, b, c)
            else:
                m.Faces.AddFace(a, b, c, d)
        cache[key] = m
    m = m.DuplicateMesh()
    for k, (x, y, z) in enumerate(V.tolist()):
        m.Vertices.SetVertex(k, x, y, z)
    m.Normals.ComputeNormals()
    m.FaceNormals.ComputeFaceNormals()
    return m


def run_grasshopper():
    comp = ghenv.Component
    store = sc.sticky.setdefault("image_relief", {})
    st = store.setdefault(str(comp.InstanceGuid), {})
    t0 = time.time()
    DIAG.clear()

    path = _val("image", None)
    if not path:
        comp.Message = "no image"
        return None, None, "Connect an image file path to `image` (File Path component)."
    H0 = _image(st, path)

    blur_px = float(_val("blur", 1.0))
    invert = bool(_val("invert", False))
    mirror = bool(_val("mirror", False))
    mode = "mirror" if mirror else "repeat"
    swap = bool(_val("swap", False))
    H = np.rot90(H0, 1) if swap else H0    # image runs up the surface instead
    procs = st.setdefault("procs", {})

    def processed(nu, nv, ru, rv):
        """The image at the size this mesh needs (about 2 px per vertex),
        blurred and stretched. Cached: only blur / invert / mirror / swap
        or a big change in size recompute it."""
        tw = int(min(H.shape[1], 16 * math.ceil(2.0 * (nu - 1) / max(ru, 1e-9) / 16)))
        th = int(min(H.shape[0], 16 * math.ceil(2.0 * (nv - 1) / max(rv, 1e-9) / 16)))
        key = (st["img_key"], swap, tw, th, blur_px, invert, mode)
        if key not in procs:
            Hs = downsample(H, tw, th)
            s = Hs.shape[1] / float(H.shape[1])
            Hs = normalize(smooth_image(Hs, blur_px * s, blur_px * s, mode))
            if len(procs) > 8:
                procs.clear()
            procs[key] = 1.0 - Hs if invert else Hs
        return procs[key]
    t_img = time.time()

    geos = _val("base", None)
    surfaces, skipped = _surfaces(geos if isinstance(geos, (list, tuple)) else [geos])
    if not surfaces:
        surfaces = [(_default_base(), False)]

    detail = max(4, int(_val("detail", 300)))
    tile_size = float(_val("tile_size", 0.0))
    gamma = float(_val("gamma", 1.0))
    fade = float(_val("fade", 0.0))
    thickness = float(_val("thickness", 0.0))
    cap = bool(_val("cap", False))
    flip = bool(_val("flip", False))
    off_u = float(_val("offset_u", 0.0))
    off_v = float(_val("offset_v", 0.0))

    grids = st.setdefault("grids", {})
    meshes = st.setdefault("meshes", {})
    out_m, out_h, notes = [], [], []
    used_grids, used_meshes = set(), set()
    want_heights = any(o.NickName == "heights" and o.Recipients.Count > 0
                       for o in comp.Params.Output)
    t_srf = t_rel = t_mesh = 0.0
    for srf, reversed_face in surfaces:
        t1 = time.time()
        cu, cv = srf.IsClosed(0), srf.IsClosed(1)
        skey = _srf_key(srf)
        # --- grid size from the real surface size ---------------------------
        coarse = grids.get((skey, "coarse"))
        if coarse is None:
            coarse = surface_lengths(_evaluate(srf, 24, 24, cu, cv)[0], cu, cv)
            grids[(skey, "coarse")] = coarse
        used_grids.add((skey, "coarse"))
        len_u, len_v = coarse
        spacing = max(len_u, len_v) / detail
        nu = max(4, int(math.ceil(len_u / spacing)) + 1)
        nv = max(4, int(math.ceil(len_v / spacing)) + 1)
        if nu * nv > MAX_VERTS:
            s = math.sqrt(MAX_VERTS / float(nu * nv))
            nu, nv = max(4, int(nu * s)), max(4, int(nv * s))
        gkey = (skey, nu, nv)
        if gkey not in grids:              # the slow part (Rhino calls): cached
            grids[gkey] = _evaluate(srf, nu, nv, cu, cv)
        used_grids.add(gkey)
        P, N = grids[gkey]
        if reversed_face != flip:
            N = -N

        # --- relief ----------------------------------------------------------
        depth = _val("depth", None)
        depth = 0.03 * min(len_u, len_v) if depth is None else float(depth)
        ru, rv = plan_repeats(len_u, len_v, H.shape[1], H.shape[0], tile_size, cu, cv, mirror)
        t_srf += time.time() - t1
        t1 = time.time()
        Pd, h = relief(P, N, processed(nu, nv, ru, rv), depth, ru, rv, cu, cv,
                       mode, off_u, off_v, fade, gamma)
        back = P - N * thickness if thickness > 0 else None
        V, F = build_mesh(Pd, cu, cv, back, cap)
        mkey = (nu, nv, cu, cv, back is not None, cap)
        if (reversed_face != flip):        # keep faces facing the relief side
            F = F[:, [0, 3, 2, 1]]
            mkey += ("flipped",)
        used_meshes.add(mkey)
        t_rel += time.time() - t1
        t1 = time.time()
        m = _to_rhino_mesh(V, F, meshes, mkey)
        if back is not None or (cap and cu != cv):
            if m.SolidOrientation() == -1:
                m.Flip(True, True, True)
        t_mesh += time.time() - t1
        out_m.append(m)
        if want_heights:                    # thousands of numbers: only if used
            out_h.extend(h.reshape(-1).tolist())
        notes.append("{}x{} vertices, {:.3g} x {:.3g} copies, depth {:g}".format(
            nu, nv, ru, rv, depth))

    for k in list(grids):                  # forget grids of deleted surfaces
        if k not in used_grids:
            grids.pop(k)
    for k in list(meshes):
        if k not in used_meshes:
            meshes.pop(k)

    ms = (time.time() - t0) * 1000.0
    txt = ("image {}x{} ({}) | {} | {:.0f} ms (image {:.0f}, surface {:.0f}, "
           "relief {:.0f}, mesh {:.0f})").format(
        H0.shape[1], H0.shape[0], st.get("img_by"), " ; ".join(notes), ms,
        (t_img - t0) * 1000, t_srf * 1000, t_rel * 1000, t_mesh * 1000)
    if not want_heights:
        txt += " | heights: connect it to get values"
    if skipped:
        txt += " | skipped {} trimmed / unsupported face(s)".format(skipped)
    bad = [m for m in out_m if not m.IsValid]
    if bad:
        ok, log = bad[0].IsValidWithLog()
        txt += " | INVALID MESH: " + str(log).strip().replace("\n", " ")[:300]
    names = [o.NickName for o in comp.Params.Output]
    for need in ("mesh", "heights", "info"):
        if need not in names:
            txt += " | WARNING: no output named '{}' (outputs are: {})".format(need, ", ".join(names))
    txt += " | paths: surface {}, mesh {}".format(DIAG.get("surface", "cached"), DIAG.get("mesh", "-"))
    errs = ["{} failed ({})".format(k[5:], v) for k, v in sorted(DIAG.items()) if k.startswith("mesh_")]
    if errs:
        txt += " | " + "; ".join(errs)
    comp.Message = "{} vertices, {:.0f} ms".format(sum(m.Vertices.Count for m in out_m), ms)
    return (out_m[0] if len(out_m) == 1 else out_m), (out_h or None), txt


if IN_RHINO:
    mesh, heights, info = run_grasshopper()
