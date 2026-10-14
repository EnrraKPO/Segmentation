"""Shape features on a single grid fixed to the model: a shell of voxels over the surface.

Each surface cell holds the surface direction at the surface point nearest its centre (interpolated, as a texel does). A line is where the
direction turns fastest across the surface, measured over the whole neighbourhood at once (shape operator fit),
judged per cell against `sensitivity`
degrees per cell. Nothing else: no averaging, no viewpoint.

run(P, F, cells=300, sensitivity=10) -> dict(idx (k,3) cell indices, N (k,3) directions, turn (k,), line (k,) bool,
                                            face (k,) a mesh face in the cell, size (mm))
"""
import numpy as np
from scipy.spatial import cKDTree

OFFS = np.array([[i, j, k] for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1) if (i, j, k) != (0, 0, 0)])


def _cell_samples(P, F, VN, size, spp, flat=False):
    """flat=True: each sample carries the facing of its own triangle (the print's surface), no corner blending"""
    """sample points on every face densely enough that every cell the surface passes through gets hit"""
    A, B, C = P[F[:, 0]], P[F[:, 1]], P[F[:, 2]]
    longest = np.maximum.reduce([np.linalg.norm(B-A, axis=1), np.linalg.norm(C-B, axis=1), np.linalg.norm(A-C, axis=1)])
    n_sub = np.maximum(1, np.ceil(spp*longest/size)).astype(int)         # subdivisions per face edge
    pts, nrm, fid = [], [], []
    for k in np.unique(n_sub):
        sel = np.nonzero(n_sub == k)[0]
        u = np.arange(k)+0.5
        bu, bv = np.meshgrid(u, u); m = (bu+bv) < k; bu, bv = bu[m]/k, bv[m]/k                # points in the triangle
        S = A[sel][:, None]+(B[sel]-A[sel])[:, None]*bu[None, :, None]+(C[sel]-A[sel])[:, None]*bv[None, :, None]
        if flat:
            fn = np.cross(B[sel]-A[sel], C[sel]-A[sel]); fn /= np.linalg.norm(fn, axis=1, keepdims=True)+1e-12
            Ns = np.repeat(fn[:, None, :], len(bu), 1)
        else:
            Ns = VN[F[sel, 0]][:, None]*(1-bu-bv)[None, :, None]+VN[F[sel, 1]][:, None]*bu[None, :, None]+VN[F[sel, 2]][:, None]*bv[None, :, None]
        pts.append(S.reshape(-1, 3)); nrm.append(Ns.reshape(-1, 3)); fid.append(np.repeat(sel, len(bu)))
    return np.vstack(pts), np.vstack(nrm), np.concatenate(fid)


KB = 1 << 16   # index packing base (cells per axis must stay below this)


def run(P, F, cells=300, sensitivity=10.0, log=None, height=None, spp=2.0, flat=False):
    """spp: sample points per cube edge on each triangle; the cube reads its arrow from the sample nearest its centre,
    so higher spp reads closer to the centre (error ~ 1/(2 spp) cube)"""
    """height: reference model height (mm) for the cell size and the turn scale; default the mesh's own extent, so a
    band cut out of a model can be gridded at the model's density by passing the full model height"""
    lo = P.min(0); ext = P.max(0)-lo; height = ext.max() if height is None else height; size = height/cells
    import trimesh
    VN = trimesh.Trimesh(P, F, process=False).vertex_normals
    pts, nrm, fid = _cell_samples(P, F, VN, size, spp, flat)
    ijk = np.floor((pts-lo)/size).astype(np.int64)
    key = (ijk[:, 0]*KB+ijk[:, 1])*KB+ijk[:, 2]
    uk, inv = np.unique(key, return_inverse=True); k = len(uk)
    idx = np.stack([uk // (KB*KB), (uk // KB) % KB, uk % KB], 1)
    # direction per cell = the interpolated surface direction at the surface point nearest the cell centre (a texel)
    centre = lo+(idx[inv]+0.5)*size; dc = np.linalg.norm(pts-centre, axis=1)
    order = np.lexsort((dc, inv)); first = np.r_[True, inv[order][1:] != inv[order][:-1]]; pick = order[first]
    N = nrm[pick]; N /= np.linalg.norm(N, axis=1, keepdims=True)+1e-12
    face = fid[pick]
    if log: log(f'grid {cells} cells/height, cell {size:.3f} mm: {k} surface cells from {len(pts)} samples')
    # neighbours: the 26 surrounding cells that are also surface cells
    lut = dict(zip(uk.tolist(), range(k)))
    nb = np.full((k, 26), -1, np.int64)
    for q, o in enumerate(OFFS):
        kk = ((idx[:, 0]+o[0])*KB+(idx[:, 1]+o[1]))*KB+(idx[:, 2]+o[2])
        nb[:, q] = [lut.get(int(x), -1) for x in kk]
    has = nb >= 0
    # turn of the direction across the whole neighbourhood at once: least-squares fit of dN = S . dx over all present
    # neighbours (dx in the cell's tangent plane, dN the change of direction). S is the shape operator; its largest
    # eigenvalue is the strongest rate of turn, its eigenvector the direction across which the surface turns.
    ref = cells/1000.0                                                  # degrees per reference cell (1/1000 height)
    a = np.where((np.abs(N[:, 0]) < 0.9)[:, None], np.array([1., 0, 0]), np.array([0, 1., 0]))
    t1 = np.cross(N, a); t1 /= np.linalg.norm(t1, axis=1, keepdims=True)+1e-12; t2 = np.cross(N, t1)
    A11 = np.zeros(k); A12 = np.zeros(k); A22 = np.zeros(k); B11 = np.zeros(k); B12 = np.zeros(k); B21 = np.zeros(k); B22 = np.zeros(k)
    for q, o in enumerate(OFFS):
        h = has[:, q]; j = nb[h, q]
        dx = (o.astype(float))[None, :]*np.ones((h.sum(), 1))               # offset in cells
        dx -= (dx*N[h]).sum(1, keepdims=True)*N[h]                         # in the tangent plane
        dn = N[j]-N[h]
        ex, ey = (dx*t1[h]).sum(1), (dx*t2[h]).sum(1); nx, ny = (dn*t1[h]).sum(1), (dn*t2[h]).sum(1)
        A11[h] += ex*ex; A12[h] += ex*ey; A22[h] += ey*ey
        B11[h] += nx*ex; B12[h] += nx*ey; B21[h] += ny*ex; B22[h] += ny*ey
    det = A11*A22-A12*A12+1e-12; I11, I12, I22 = A22/det, -A12/det, A11/det
    S11 = B11*I11+B12*I12; S12 = B11*I12+B12*I22; S21 = B21*I11+B22*I12; S22 = B21*I12+B22*I22
    s12 = 0.5*(S12+S21); tr = S11+S22; dd = np.sqrt(np.maximum(((S11-S22)/2)**2+s12**2, 0))
    l1, l2 = tr/2+dd, tr/2-dd
    big = np.abs(l1) >= np.abs(l2); k1 = np.where(big, l1, l2)
    turn = np.degrees(np.abs(k1))*ref                                      # rate of turn, degrees per reference cell
    vx = np.where(np.abs(s12) > 1e-12, s12, np.where(np.abs(S11-k1) < np.abs(S22-k1), 0.0, 1.0))
    vy = np.where(np.abs(s12) > 1e-12, k1-S11, np.where(np.abs(S11-k1) < np.abs(S22-k1), 1.0, 0.0))
    nv = np.sqrt(vx*vx+vy*vy)+1e-18; d1 = (vx/nv)[:, None]*t1+(vy/nv)[:, None]*t2
    nnb = has.sum(1); turn[nnb < 6] = 0                                    # a cell needs a neighbourhood to be judged
    # crest: the turn must peak across d1 (nearest cells along +-d1)
    def along(sign):
        step = np.rint(sign*d1).astype(int)
        kk = ((idx[:, 0]+step[:, 0])*KB+(idx[:, 1]+step[:, 1]))*KB+(idx[:, 2]+step[:, 2])
        j = np.array([lut.get(int(x), -1) for x in kk]); return np.where(j >= 0, turn[np.maximum(j, 0)], 0.0)
    tp, tm = along(1.0), along(-1.0)
    line = (turn >= sensitivity) & (turn >= tp) & (turn > tm)
    if log: log(f'turn per cell: p50 {np.percentile(turn, 50):.1f} p90 {np.percentile(turn, 90):.1f} p99 {np.percentile(turn, 99):.1f} deg; line cells {line.sum()}')
    return dict(idx=idx, N=N, turn=turn, line=line, face=face, size=size, lo=lo, d1=d1)


def _closest_on_triangles(Q, A, B, C):
    """closest point on triangle ABC to Q, per row (Ericson, Real-Time Collision Detection 5.1.5), vectorised"""
    ab = B-A; ac = C-A; ap = Q-A
    d1 = (ab*ap).sum(1); d2 = (ac*ap).sum(1)
    bp = Q-B; d3 = (ab*bp).sum(1); d4 = (ac*bp).sum(1)
    cp = Q-C; d5 = (ab*cp).sum(1); d6 = (ac*cp).sum(1)
    vc = d1*d4-d3*d2; vb = d5*d2-d1*d6; va = d3*d6-d5*d4
    out = np.empty_like(Q); done = np.zeros(len(Q), bool)
    def put(m, pts): m = m & ~done; out[m] = pts[m] if pts.ndim == 2 else pts; done[m] = True
    put((d1 <= 0) & (d2 <= 0), A)
    put((d3 >= 0) & (d4 <= d3), B)
    put((d6 >= 0) & (d5 <= d6), C)
    v = d1/np.where(d1-d3 != 0, d1-d3, 1); put((vc <= 0) & (d1 >= 0) & (d3 <= 0), A+v[:, None]*ab)
    w = d2/np.where(d2-d6 != 0, d2-d6, 1); put((vb <= 0) & (d2 >= 0) & (d6 <= 0), A+w[:, None]*ac)
    den = (d4-d3)+(d5-d6); w = (d4-d3)/np.where(den != 0, den, 1); put((va <= 0) & (d4-d3 >= 0) & (d5-d6 >= 0), B+w[:, None]*(C-B))
    s = va+vb+vc; s = np.where(s != 0, s, 1); v = vb/s; w = vc/s; put(~done, A+v[:, None]*ab+w[:, None]*ac)
    return out


def _tri_box_overlap(Q, A, B, C, half):
    """exact triangle / cube overlap test (separating axes, Akenine-Moller), per row; cube centre Q, half-size half"""
    v0, v1, v2 = A-Q, B-Q, C-Q
    e0, e1, e2 = v1-v0, v2-v1, v0-v2
    ok = np.ones(len(Q), bool)
    # the 3 cube axes
    for a in range(3):
        mn = np.minimum(np.minimum(v0[:, a], v1[:, a]), v2[:, a]); mx = np.maximum(np.maximum(v0[:, a], v1[:, a]), v2[:, a])
        ok &= ~((mn > half) | (mx < -half))
    # the triangle plane
    n = np.cross(e0, e1); r = half*np.abs(n).sum(1); d = (n*v0).sum(1); ok &= np.abs(d) <= r
    # the 9 edge x axis directions
    for e in (e0, e1, e2):
        for a in range(3):
            ax = np.zeros_like(e); ax[:, a] = 1.0; L = np.cross(e, ax)
            p0 = (L*v0).sum(1); p1 = (L*v1).sum(1); p2 = (L*v2).sum(1)
            r = half*np.abs(L).sum(1)
            ok &= ~((np.minimum(np.minimum(p0, p1), p2) > r) | (np.maximum(np.maximum(p0, p1), p2) < -r))
    return ok


def run_direct(P, F, cells=300, log=None, height=None, flat=True, batch=4_000_000):
    """Surface cubes straight from the triangles, no sample points. For each triangle, every cube in the box its
    corners span is tested exactly for overlap with the triangle (separating axes). A cube overlapped by several
    triangles keeps the one whose nearest point is closest to the cube centre
    and reads its facing (flat: the triangle's own; else the corner-blended facing at the nearest point).
    Returns the same fields as run(); turn/line are zeros (creases.fold_change computes the fold change)."""
    import trimesh
    lo = P.min(0); ext = P.max(0)-lo; height = ext.max() if height is None else height; size = height/cells
    A, B, C = P[F[:, 0]], P[F[:, 1]], P[F[:, 2]]
    T = np.stack([A, B, C], 1); g = (T-lo)/size
    imin = np.floor(g.min(1)).astype(np.int64); imax = np.floor(g.max(1)).astype(np.int64); span = imax-imin+1
    VN = trimesh.Trimesh(P, F, process=False).vertex_normals if not flat else None
    fn = np.cross(B-A, C-A); fn /= np.linalg.norm(fn, axis=1, keepdims=True)+1e-12
    half = size*0.5
    keys, dists, fids, nrms = [], [], [], []
    skey = (span[:, 0]*KB+span[:, 1])*KB+span[:, 2]; order = np.argsort(skey, kind='stable'); skey_s = skey[order]
    starts = np.r_[0, np.nonzero(skey_s[1:] != skey_s[:-1])[0]+1, len(order)]
    ntest = 0
    for a, b in zip(starts[:-1], starts[1:]):
        grp = order[a:b]; sp = span[grp[0]]; ncand = int(sp.prod())
        offs = np.stack(np.meshgrid(np.arange(sp[0]), np.arange(sp[1]), np.arange(sp[2]), indexing='ij'), -1).reshape(-1, 3)
        step = max(1, batch//ncand)
        for c0 in range(0, len(grp), step):
            tri = grp[c0:c0+step]; n = len(tri)
            cid = imin[tri][:, None, :]+offs[None]                          # (n, ncand, 3) cube indices
            centre = lo+(cid+0.5)*size
            Q = centre.reshape(-1, 3); rep = np.repeat(tri, ncand)
            inside = _tri_box_overlap(Q, A[rep], B[rep], C[rep], half); ntest += len(Q)
            if not inside.any(): continue
            Q = Q[inside]; rep = rep[inside]; cid = cid.reshape(-1, 3)[inside]
            X = _closest_on_triangles(Q, A[rep], B[rep], C[rep]); d = np.linalg.norm(X-Q, axis=1)
            keys.append((cid[:, 0]*KB+cid[:, 1])*KB+cid[:, 2]); dists.append(d.astype(np.float32)); fids.append(rep.astype(np.int32))
            if not flat:
                Xi = X; Ai, Bi, Ci = A[rep], B[rep], C[rep]
                v0 = Bi-Ai; v1 = Ci-Ai; v2 = Xi-Ai
                d00 = (v0*v0).sum(1); d01 = (v0*v1).sum(1); d11 = (v1*v1).sum(1); d20 = (v2*v0).sum(1); d21 = (v2*v1).sum(1)
                den = d00*d11-d01*d01; den = np.where(den != 0, den, 1); bv = (d11*d20-d01*d21)/den; bw = (d00*d21-d01*d20)/den; bu = 1-bv-bw
                Nn = VN[F[rep, 0]]*bu[:, None]+VN[F[rep, 1]]*bv[:, None]+VN[F[rep, 2]]*bw[:, None]; nrms.append(Nn.astype(np.float32))
    key = np.concatenate(keys); dist = np.concatenate(dists); fid = np.concatenate(fids)
    o = np.lexsort((dist, key)); key, dist, fid = key[o], dist[o], fid[o]
    first = np.r_[True, key[1:] != key[:-1]]; uk = key[first]; face = fid[first]
    idx = np.stack([uk // (KB*KB), (uk // KB) % KB, uk % KB], 1)
    N = fn[face].copy() if flat else np.concatenate(nrms)[o][first].astype(np.float64)
    N /= np.linalg.norm(N, axis=1, keepdims=True)+1e-12
    k = len(idx)
    if log: log(f'grid {cells} cells/height, cell {size:.3f} mm: {k} surface cells, {ntest/1e6:.0f}M cube-triangle tests (direct)')
    z = np.zeros(k)
    return dict(idx=idx, N=N, turn=z, line=np.zeros(k, bool), face=face, size=size, lo=lo, d1=N.copy())
