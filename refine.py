"""Refinement of a line from the mesh itself, at a resolution far below the walk's cubes (session 2026-10-21).

Stage 3 of the pipeline: the walk (trace4) gives a coarse line on each crease; this module moves it onto the crease
exactly as the surface shows it, then draws it as a clean curve.

Measurement (refine_line): the line is resampled every STEP mm. At each point a slice is walked across the line on the
surface in FINE steps, the facing read as the render blends it (corner normals interpolated across the triangle). The
surface on each side is read as a plane: a weighted average of position and facing where every weight is a smooth
function of distance (rises in from the point so the strip of connecting triangles counts little, fades out so the
rounding further away counts little) and the planes are averaged gradually along the line. The point moves to where
the two planes cross. Repeated from the new point until it settles. No cut-offs, no picks, no boxes anywhere:
neighbouring points 0.02 mm apart must land within a hair of each other, since the surface is continuous.

Curve (fit_curve): the measured path is ironed; where it still turns more than CORNER it breaks (two lines meeting);
where it bends more than BEND it gets a joint (one continuous curve through it). One cubic spline per stretch between
corners, with a knot at each joint; where the curve misses the measurements by more than TOL a joint is added there
and the stretch refitted. The result is a few control points per line and a smooth sweep between them.
"""
import numpy as np, trimesh
from scipy.spatial import cKDTree
from scipy.interpolate import make_lsq_spline

DEFAULTS = dict(step=0.02, fine=0.0087, range=0.6, sin=0.06, sout=0.2, salong=0.05, passes=6,
                iron=0.5, chord=0.8, bend=15.0, corner=45.0, tol=0.05)


class Band:
    """the part of the mesh near a set of lines, with shared corners so facings blend across triangles"""
    def __init__(self, P, F, near_pts, radius):
        cen = P[F].mean(1); keep = cKDTree(near_pts).query(cen)[0] <= radius; Fp = F[keep]
        used = np.unique(Fp); rm = -np.ones(len(P), int); rm[used] = np.arange(len(used))
        self.m = trimesh.Trimesh(P[used], rm[Fp], process=False); self.m.merge_vertices()
        self.FN = self.m.face_normals; self.VN = self.m.vertex_normals; self.n_faces = len(self.m.faces)

    def drop(self, Q):
        """nearest surface point and its triangle"""
        C, _, t = trimesh.proximity.closest_point(self.m, Q); return C, t

    def facing(self, Q, t):
        """the facing the render shows at surface point Q on triangle t: the corner facings blended by position"""
        tri = self.m.triangles[t]; A, B, C = tri[:, 0], tri[:, 1], tri[:, 2]; v0 = B-A; v1 = C-A; v2 = Q-A
        d00 = (v0*v0).sum(1); d01 = (v0*v1).sum(1); d11 = (v1*v1).sum(1); d20 = (v2*v0).sum(1); d21 = (v2*v1).sum(1)
        den = d00*d11-d01*d01; den = np.where(np.abs(den) > 1e-18, den, 1); v = (d11*d20-d01*d21)/den; w = (d00*d21-d01*d20)/den; u = 1-v-w
        fv = self.m.faces[t]; n = self.VN[fv[:, 0]]*u[:, None]+self.VN[fv[:, 1]]*v[:, None]+self.VN[fv[:, 2]]*w[:, None]
        return n/(np.linalg.norm(n, axis=1, keepdims=True)+1e-12)

    def push(self, Q, n):
        """back onto the surface along direction n (continuous), nearest point only where the ray misses"""
        out = Q.copy(); hit = np.zeros(len(Q), bool)
        for sgn in (1, -1):
            if hit.all(): break
            loc, idr, _ = self.m.ray.intersects_location(Q[~hit]-sgn*0.3*n[~hit], sgn*n[~hit], multiple_hits=False)
            if len(idr): ii = np.flatnonzero(~hit)[idr]; out[ii] = loc.reshape(-1, 3); hit[ii] = True
        if (~hit).any(): out[~hit] = self.drop(Q[~hit])[0]
        return out


def resample(L, step):
    d = np.r_[0, np.cumsum(np.linalg.norm(np.diff(L, axis=0), axis=1))]
    if d[-1] < 2*step: return L.copy()
    t = np.arange(0, d[-1], step); return np.stack([np.interp(t, d, L[:, i]) for i in range(3)], 1)


def refine_line(L, band, step=0.02, fine=0.0087, sin=0.06, sout=0.2, salong=0.05, passes=6, **kw):
    """L: coarse line (mm). Returns the measured line (same spacing as STEP) and the mean facing at each point."""
    rng = kw.get('range', 0.6); L = resample(L, step); K = len(L)
    if K < 4: C, t = band.drop(L); return C, band.facing(C, t)
    T = np.gradient(L, axis=0); T /= np.linalg.norm(T, axis=1, keepdims=True)+1e-12
    dl = np.r_[0, np.cumsum(np.linalg.norm(np.diff(L, axis=0), axis=1))]
    WA = np.exp(-((dl[:, None]-dl[None, :])/salong)**2); WA /= WA.sum(1, keepdims=True)
    nst = int(round(rng/fine)); S = np.arange(-nst, nst+1)*fine
    say = (1-np.exp(-(S/sin)**2))*np.exp(-(S/sout)**2); wA = np.where(S < 0, say, 0.0); wB = np.where(S > 0, say, 0.0)

    def march(C0, N0, sign):
        q = C0.copy(); u = np.cross(N0, T)*sign; u /= np.linalg.norm(u, axis=1, keepdims=True)+1e-12; pos = [q.copy()]; fac = [N0.copy()]
        for _ in range(nst):
            q2, tq = band.drop(q+fine*u); n2 = band.facing(q2, tq); u = u-(u*n2).sum(1, keepdims=True)*n2; u /= np.linalg.norm(u, axis=1, keepdims=True)+1e-12
            q = q2; pos.append(q.copy()); fac.append(n2.copy())
        return np.stack(pos, 1), np.stack(fac, 1)

    start = L.copy(); nA = nB = None
    for it in np.arange(passes):
        C0, T0 = band.drop(start); N0 = WA @ band.facing(C0, T0)
        if nA is not None: N0 = nA+nB
        N0 /= np.linalg.norm(N0, axis=1, keepdims=True)+1e-12
        pp, fp = march(C0, N0, 1); pm, fm = march(C0, N0, -1)
        pos = np.concatenate([pm[:, :0:-1], pp], 1); fac = np.concatenate([fm[:, :0:-1], fp], 1)
        nA = (fac*wA[None, :, None]).sum(1)/wA.sum(); pA = (pos*wA[None, :, None]).sum(1)/wA.sum()
        nB = (fac*wB[None, :, None]).sum(1)/wB.sum(); pB = (pos*wB[None, :, None]).sum(1)/wB.sum()
        nA, pA, nB, pB = WA @ nA, WA @ pA, WA @ nB, WA @ pB
        nA /= np.linalg.norm(nA, axis=1, keepdims=True)+1e-12; nB /= np.linalg.norm(nB, axis=1, keepdims=True)+1e-12
        R = C0.copy()
        for k in np.arange(K):
            A = np.stack([nA[k], nB[k], T[k]]); b = np.array([nA[k] @ pA[k], nB[k] @ pB[k], T[k] @ C0[k]])
            if abs(np.linalg.det(A)) > 1e-6: R[k] = np.linalg.solve(A, b)
        nn = nA+nB; nn /= np.linalg.norm(nn, axis=1, keepdims=True)+1e-12
        new = band.push(R, nn); shift = np.linalg.norm(new-start, axis=1); start = new
        if np.percentile(shift, 90) < 0.002: break
    return R, nn


def fit_curve(R, iron=0.5, chord=0.8, bend=15.0, corner=45.0, tol=0.05, **_):
    """R: measured line. Returns (curve points every 0.02 mm, corner points, joint points)."""
    d = np.r_[0, np.cumsum(np.linalg.norm(np.diff(R, axis=0), axis=1))]
    if len(R) < 8 or d[-1] < 2*chord: return R.copy(), R[:0], R[:0]
    W = np.exp(-((d[:, None]-d[None, :])/iron)**2); I = (W @ R)/W.sum(1, keepdims=True)
    at = lambda Q, s: np.stack([np.interp(s, d, Q[:, i]) for i in range(3)], 1)
    a = I-at(I, np.maximum(d-chord, 0)); b = at(I, np.minimum(d+chord, d[-1]))-I
    a /= np.linalg.norm(a, axis=1, keepdims=True)+1e-12; b /= np.linalg.norm(b, axis=1, keepdims=True)+1e-12
    turn = np.degrees(np.arccos(np.clip((a*b).sum(1), -1, 1))); turn[(d < chord) | (d > d[-1]-chord)] = 0
    cand = np.flatnonzero((turn >= bend) & (turn >= np.r_[0, turn[:-1]]) & (turn >= np.r_[turn[1:], 0])); picked = []
    for i in cand[np.argsort(-turn[cand])]:
        if all(abs(d[i]-d[j]) > chord for j in picked): picked.append(i)
    corners = sorted(i for i in picked if turn[i] >= corner); joints = sorted(i for i in picked if turn[i] < corner)

    def fit(s, e, knots):
        u = d[s:e+1]-d[s]; P = R[s:e+1]; w = np.ones(len(u)); w[0] = w[-1] = 100; k = 3
        inner = sorted(d[j]-d[s] for j in knots if s < j < e and d[j]-d[s] > 0.1 and d[e]-d[j] > 0.1)
        t = np.r_[[0]*(k+1), inner, [u[-1]]*(k+1)]
        if len(u) <= len(t)-k-1: return P, np.zeros(len(u))
        spl = make_lsq_spline(u, P, t, k=k, w=w); return spl(np.arange(0, u[-1], 0.02)), np.linalg.norm(spl(u)-P, axis=1)
    bounds = [0]+corners+[len(R)-1]; added = 0
    while added < 40:
        worst = None
        for s, e in zip(bounds[:-1], bounds[1:]):
            _, miss = fit(s, e, joints); ok = np.ones(len(miss), bool); ok[:6] = ok[-6:] = False
            for q in joints: ok &= np.abs(d[s:e+1]-d[q]) > 0.15
            if not ok.any(): continue
            j = int(np.argmax(np.where(ok, miss, -1)))
            if miss[j] > tol and (worst is None or miss[j] > worst[0]): worst = (miss[j], s+j)
        if worst is None: break
        joints.append(worst[1]); added += 1
    joints = sorted(joints); pieces = [fit(s, e, joints)[0] for s, e in zip(bounds[:-1], bounds[1:])]
    return np.vstack(pieces+[R[-1:]]), R[corners], R[joints]


def refine_all(curves, P, F, log=None, only=None, **kw):
    """curves: the walk's lines (mm). only: index of one line, or None for all. Returns dict(lines, measured, corners, joints)."""
    import time
    kw = {**DEFAULTS, **kw}; idx = [only] if only is not None else list(np.arange(len(curves)))
    t0 = time.time(); band = Band(P, F, np.vstack([curves[i] for i in idx]), kw['range']+0.6)
    if log: log(f'refine: {len(idx)} lines, band of {band.n_faces} triangles, {time.time()-t0:.1f} s')
    lines, measured, corners, joints = [], [], [], []
    for n, i in enumerate(idx):
        R, nn = refine_line(np.asarray(curves[i], float), band, **kw); C, cp, jp = fit_curve(R, **kw)
        lines.append(C); measured.append(R); corners.append(cp); joints.append(jp)
        if log and (n+1) % 50 == 0: log(f'  refined {n+1} of {len(idx)}, {time.time()-t0:.0f} s')
    if log: log(f'refine done: {len(lines)} lines, {sum(len(c) for c in corners)} corners, {sum(len(j) for j in joints)} joints, {time.time()-t0:.1f} s')
    return dict(lines=lines, measured=measured, corners=corners, joints=joints, index=idx)


def nearest_line(curves, p):
    """index of the line with a point nearest p"""
    p = np.asarray(p, float); return int(np.argmin([np.linalg.norm(np.asarray(c)-p, axis=1).min() for c in curves]))
