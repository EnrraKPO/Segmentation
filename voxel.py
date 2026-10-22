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
    # the triangle plane first: it is cheap and alone rejects most of the pairs (a cube the plane misses cannot hold the
    # triangle); the other twelve checks run only on its survivors. Same hits as testing everything on everyone.
    n = np.cross(e0, e1); r = half*np.abs(n).sum(1); d = (n*v0).sum(1); ok = np.abs(d) <= r
    s = np.nonzero(ok)[0]
    if len(s) == 0: return ok
    v0, v1, v2, e0, e1, e2 = v0[s], v1[s], v2[s], e0[s], e1[s], e2[s]; oks = np.ones(len(s), bool)
    # the 3 cube axes
    for a in range(3):
        mn = np.minimum(np.minimum(v0[:, a], v1[:, a]), v2[:, a]); mx = np.maximum(np.maximum(v0[:, a], v1[:, a]), v2[:, a])
        oks &= ~((mn > half) | (mx < -half))
    # the 9 edge x axis directions
    for e in (e0, e1, e2):
        for a in range(3):
            ax = np.zeros_like(e); ax[:, a] = 1.0; L = np.cross(e, ax)
            p0 = (L*v0).sum(1); p1 = (L*v1).sum(1); p2 = (L*v2).sum(1)
            r = half*np.abs(L).sum(1)
            oks &= ~((np.minimum(np.minimum(p0, p1), p2) > r) | (np.maximum(np.maximum(p0, p1), p2) < -r))
    ok[s] = oks
    return ok


def _tri_in_cube_area(A, B, C, half):
    """area of triangle ABC (given relative to the cube centre) inside the cube [-half, half]^3, per row: the triangle is
    clipped against the six cube planes (Sutherland-Hodgman, a convex polygon of at most 9 corners), then its area taken"""
    n = len(A); V = np.zeros((n, 10, 3)); V[:, 0] = A; V[:, 1] = B; V[:, 2] = C; cnt = np.full(n, 3); slot = np.arange(10)
    for axis in range(3):
        for sign in (1.0, -1.0):                                                 # keep sign*v[axis] <= half
            outside = (sign*V[:, :, axis] > half) & (slot[None, :] < cnt[:, None])
            sel = np.nonzero(outside.any(1))[0]                                  # only the polygons this plane cuts; the others are left as they are
            if len(sel) == 0: continue
            Vs = V[sel]; cs = cnt[sel]; rows = np.arange(len(sel)); out = np.zeros_like(Vs); oc = np.zeros(len(sel), int)
            def emit(m, pt):
                out[rows[m], oc[m]] = pt[m]; oc[m] += 1
            for k in range(9):
                valid = k < cs
                if not valid.any(): break
                i = np.minimum(k, cs-1); j = (k+1) % np.maximum(cs, 1)
                Pi = Vs[rows, i]; Pj = Vs[rows, j]; di = sign*Pi[:, axis]-half; dj = sign*Pj[:, axis]-half
                ini = di <= 0; inj = dj <= 0; t = di/np.where(di-dj != 0, di-dj, 1); cross = Pi+t[:, None]*(Pj-Pi)
                emit(valid & ini & ~inj, cross); emit(valid & ~ini & inj, cross); emit(valid & inj, Pj)
            V[sel] = out; cnt[sel] = oc
    area = np.zeros(n)
    for k in range(1, 9):                                                        # fan from corner 0
        m = k+1 < cnt
        if not m.any(): break
        area[m] += 0.5*np.linalg.norm(np.cross(V[m, k]-V[m, 0], V[m, k+1]-V[m, 0]), axis=1)
    return area


_W = {}   # per-process data for the direct builder (set by _direct_init)


def _direct_init(spec, lo, size):
    """attach the shared arrays (name, shape, dtype) published by run_direct; no copies per process"""
    from multiprocessing import shared_memory
    _W.clear(); _W['_shm'] = []
    for k, (name, shape, dtype) in spec.items():
        if name is None: _W[k] = None; continue
        shm = shared_memory.SharedMemory(name=name); _W['_shm'].append(shm)
        _W[k] = np.ndarray(shape, dtype=dtype, buffer=shm.buf)
    _W['lo'] = lo; _W['size'] = size


def _direct_task(args):
    """one batch of triangles that all span the same box of cubes: returns the cubes they overlap"""
    tri, sp, imin_tri, flat, blend = args
    A, B, C, VN, F, lo, size = (_W[k] for k in ('A', 'B', 'C', 'VN', 'F', 'lo', 'size'))
    ncand = int(np.prod(sp)); half = size*0.5
    offs = np.stack(np.meshgrid(np.arange(sp[0]), np.arange(sp[1]), np.arange(sp[2]), indexing='ij'), -1).reshape(-1, 3)
    cid = imin_tri[:, None, :]+offs[None]                                     # (n, ncand, 3) cube indices
    Q = (lo+(cid+0.5)*size).reshape(-1, 3); rep = np.repeat(tri, ncand)
    inside = _tri_box_overlap(Q, A[rep], B[rep], C[rep], half); ntest = len(Q)
    if not inside.any(): return None, None, None, None, ntest, None
    Q = Q[inside]; rep = rep[inside]; cid = cid.reshape(-1, 3)[inside]
    X = _closest_on_triangles(Q, A[rep], B[rep], C[rep]); d = np.linalg.norm(X-Q, axis=1)
    key = (cid[:, 0]*KB+cid[:, 1])*KB+cid[:, 2]; nrm = None
    # side score for ties (cube centre equally near two triangles = nearest point on their shared edge): the triangle
    # whose facing points most directly at the cube centre is the side the cube is over; lower score = better
    fnr = _W['FN'][rep]; side = (-np.abs((fnr*(Q-X)).sum(1))/np.maximum(d, 1e-12)).astype(np.float32)
    o = np.lexsort((side, d, key))
    if not blend: first = np.r_[True, key[o][1:] != key[o][:-1]]; o = o[first]   # one hit per cube: the nearest (blend: every hit is kept, the merge needs all the pairs)
    key, d, rep, X, side = key[o], d[o], rep[o], X[o], side[o]
    if not flat:
        Ai, Bi, Ci = A[rep], B[rep], C[rep]; v0 = Bi-Ai; v1 = Ci-Ai; v2 = X-Ai
        d00 = (v0*v0).sum(1); d01 = (v0*v1).sum(1); d11 = (v1*v1).sum(1); d20 = (v2*v0).sum(1); d21 = (v2*v1).sum(1)
        den = d00*d11-d01*d01; den = np.where(den != 0, den, 1); bv = (d11*d20-d01*d21)/den; bw = (d00*d21-d01*d20)/den; bu = 1-bv-bw
        nrm = (VN[F[rep, 0]]*bu[:, None]+VN[F[rep, 1]]*bv[:, None]+VN[F[rep, 2]]*bw[:, None]).astype(np.float32)
    return key, d.astype(np.float32), rep.astype(np.int32), nrm, ntest, side


def _blend_sums(pk, pr, cu):
    """(cube, triangle) pairs sorted by cube, cu = each pair's cube number from 0: each triangle clipped to its cube and its
    facing weighed by the area inside; returns the summed weighted facing per cube number"""
    A, B, C, FN, lo, size = (_W[k] for k in ('A', 'B', 'C', 'FN', 'lo', 'size')); n = int(cu[-1])+1; tot = np.zeros((n, 3))
    for s in range(0, len(pk), _BLEND_SUB):                                    # small pieces: the clipping's working arrays stay in cache
        k, r, c = pk[s:s+_BLEND_SUB], pr[s:s+_BLEND_SUB], cu[s:s+_BLEND_SUB]
        cid = np.stack([k // (KB*KB), (k // KB) % KB, k % KB], 1); Q = lo+(cid+0.5)*size
        wa = _tri_in_cube_area(A[r]-Q, B[r]-Q, C[r]-Q, size*0.5); fa = FN[r]*wa[:, None]
        for c_ in range(3): tot[:, c_] += np.bincount(c, fa[:, c_], minlength=n)
    return tot


_BLEND_SUB = 16384   # pairs per clipping piece (16k: the working arrays stay in cache; measured 5 us per pair against 8 unchunked)


def _blend_task(args):
    """one chunk of the pairs, read from shared memory (nothing is copied to the core); the chunk holds whole cubes, so
    its sums are written straight into their rows of the shared result"""
    from multiprocessing import shared_memory
    spec, s, e = args; shm = {k: shared_memory.SharedMemory(name=v[0]) for k, v in spec.items()}
    arr = {k: np.ndarray(v[1], dtype=v[2], buffer=shm[k].buf) for k, v in spec.items()}
    import time; t0 = time.time()
    cu = arr['CU'][s:e]; c0 = int(cu[0]); arr['OUT'][c0:int(cu[-1])+1] = _blend_sums(arr['PK'][s:e], arr['PR'][s:e], cu-c0)
    for m in shm.values(): m.close()
    return e-s, time.time()-t0


def run_direct(P, F, cells=300, log=None, height=None, flat=True, batch=2_000_000, workers=None, box=None, blend=False):
    """Surface cubes straight from the triangles, no sample points. box=(lo, hi) in mm: cubes are laid only inside it. For each triangle, every cube in the box its
    corners span is tested exactly for overlap with the triangle (separating axes). A cube overlapped by several
    triangles keeps the one whose nearest point is closest to the cube centre and reads its facing there (flat: the
    triangle's own; else the corner-blended facing at that point). blend=True: the cube's facing is instead the facings
    of ALL the triangles inside it, each weighted by the area it covers inside the cube (the triangle clipped to the
    cube), so the arrow changes gradually as a crease crosses from one cube to the next; `face` stays the nearest triangle. The triangles are split into batches and spread
    over `workers` processes (default: all cores). Returns the same fields as run(); turn/line are zeros
    (creases.fold_change computes the fold change)."""
    import trimesh, time, os
    from concurrent.futures import ProcessPoolExecutor
    lo = P.min(0); ext = P.max(0)-lo; height = ext.max() if height is None else height; size = height/cells
    A, B, C = P[F[:, 0]], P[F[:, 1]], P[F[:, 2]]; tri_id = np.arange(len(F))
    VN = trimesh.Trimesh(P, F, process=False).vertex_normals if not flat else None
    T = np.stack([A, B, C], 1); g = (T-lo)/size
    imin = np.floor(g.min(1)).astype(np.int64); imax = np.floor(g.max(1)).astype(np.int64)
    if box is not None:                                                        # cubes are laid only inside this box (mm); a triangle's candidate cubes are clipped to it
        bmin = np.floor((np.asarray(box[0])-lo)/size).astype(np.int64); bmax = np.floor((np.asarray(box[1])-lo)/size).astype(np.int64)
        imin = np.maximum(imin, bmin); imax = np.minimum(imax, bmax); keep = (imax >= imin).all(1)
        A, B, C, imin, imax = A[keep], B[keep], C[keep], imin[keep], imax[keep]; F = F[keep]; tri_id = tri_id[keep]
    span = imax-imin+1
    fn = np.cross(B-A, C-A); fn /= np.linalg.norm(fn, axis=1, keepdims=True)+1e-12
    skey = (span[:, 0]*KB+span[:, 1])*KB+span[:, 2]; order = np.argsort(skey, kind='stable'); skey_s = skey[order]
    starts = np.r_[0, np.nonzero(skey_s[1:] != skey_s[:-1])[0]+1, len(order)]
    tasks = []
    for a, b in zip(starts[:-1], starts[1:]):
        grp = order[a:b]; sp = span[grp[0]]; step = max(1, batch//int(sp.prod()))
        for c0 in range(0, len(grp), step):
            tri = grp[c0:c0+step]; tasks.append((tri, sp, imin[tri], flat, blend))
    tasks.sort(key=lambda t: -len(t[0])*int(np.prod(t[1])))                     # biggest batches first: balanced finish
    workers = workers or os.cpu_count() or 1
    keys, dists, fids, nrms, sides = [], [], [], [], []; ntest = 0; done = 0; t0 = time.time(); nxt = 0.1
    def take(r):
        nonlocal ntest, done, nxt
        key, d, rep, nrm, nt, side = r; ntest += nt; done += 1
        if key is not None: keys.append(key); dists.append(d); fids.append(rep); nrms.append(nrm); sides.append(side)
        if log and done >= nxt*len(tasks): log(f'  cubes: {100*done/len(tasks):.0f} % of triangle batches, {time.time()-t0:.0f} s'); nxt += 0.1
    def merge():
        """one hit per cube: the nearest triangle (ties: the one facing the cube centre most directly)"""
        t1 = time.time(); key = np.concatenate(keys); dist = np.concatenate(dists); fid = np.concatenate(fids); side = np.concatenate(sides)
        o = np.argsort(fid, kind='stable'); o = o[np.argsort(side[o], kind='stable')]; o = o[np.argsort(key[o], kind='stable')]   # by cube, then best side, then triangle number
        key, dist, fid = key[o], dist[o], fid[o]
        first = np.r_[True, key[1:] != key[:-1]]; grp = np.cumsum(first)-1                         # hits of one cube are now adjacent
        best = np.full(first.sum(), np.inf, np.float32); np.minimum.at(best, grp, dist)             # nearest hit per cube
        pick = np.nonzero(dist == best[grp])[0]; pick = pick[np.r_[True, grp[pick][1:] != grp[pick][:-1]]]   # nearest; ties (centre nearest to a shared edge) -> the triangle facing the cube centre most directly
        if log: log(f'  cubes: merging {len(grp)/1e6:.0f}M hits, {time.time()-t1:.0f} s')
        multi = None
        if blend:                                                                                  # the hits of a cube are adjacent: a run of two or more
            m = np.bincount(grp)[grp] >= 2; g = grp[m]; fm = np.r_[True, g[1:] != g[:-1]]            # is a cube holding two or more triangles, to be blended
            multi = (key[m], fid[m], np.cumsum(fm)-1, g[fm])                                       # pairs, their cube number from 0, that cube's row in the result
        return key[pick], fid[pick], o[pick], multi

    def blend_chunks(cu):
        """cut the pairs into chunks of whole cubes, a few per core"""
        n = max(1, workers*4); b = np.searchsorted(cu, cu[np.linspace(0, len(cu)-1, n+1).astype(int)]); b[-1] = len(cu); b = np.unique(b)
        return list(zip(b[:-1], b[1:]))

    if workers > 1 and len(tasks) > 1:
        from multiprocessing import shared_memory
        shms = []; spec = {}
        for k, arr in dict(A=A, B=B, C=C, VN=VN, F=F, FN=fn).items():
            if arr is None: spec[k] = (None, None, None); continue
            arr = np.ascontiguousarray(arr); shm = shared_memory.SharedMemory(create=True, size=arr.nbytes); shms.append(shm)
            np.ndarray(arr.shape, dtype=arr.dtype, buffer=shm.buf)[:] = arr; spec[k] = (shm.name, arr.shape, arr.dtype)
        try:
            with ProcessPoolExecutor(workers, initializer=_direct_init, initargs=(spec, lo, size)) as ex:
                for r in ex.map(_direct_task, tasks, chunksize=max(1, len(tasks)//(workers*8))): take(r)
                uk, fid, o, multi = merge()
                if blend:                                                        # the pairs go to the cores through shared memory, the sums come back the same way
                    t1 = time.time(); pk, pr, cu, row = multi; bspec = {}; tot = np.zeros((len(row), 3))
                    for k, arr in dict(PK=pk, PR=pr, CU=cu, OUT=tot).items():
                        arr = np.ascontiguousarray(arr); shm = shared_memory.SharedMemory(create=True, size=arr.nbytes); shms.append(shm)
                        np.ndarray(arr.shape, dtype=arr.dtype, buffer=shm.buf)[:] = arr; bspec[k] = (shm.name, arr.shape, arr.dtype)
                    rr = list(ex.map(_blend_task, [(bspec, s, e) for s, e in blend_chunks(cu)])); npairs = sum(r[0] for r in rr); busy = sum(r[1] for r in rr)
                    if log: log(f'  blend: {len(rr)} chunks, worker time {busy:.0f} s in all, longest {max(r[1] for r in rr):.1f} s')
                    tot = np.ndarray(tot.shape, dtype=tot.dtype, buffer=shms[-1].buf).copy()
                    if log: log(f'  blend: clipped on {workers} cores in {time.time()-t1:.0f} s')
        finally:
            for shm in shms: shm.close(); shm.unlink()
    else:
        _W.update(A=A, B=B, C=C, VN=VN, F=F, FN=fn, lo=lo, size=size)
        for t in tasks: take(_direct_task(t))
        uk, fid, o, multi = merge()
        if blend: pk, pr, cu, row = multi; tot = _blend_sums(pk, pr, cu); npairs = len(pk)
    face = tri_id[fid]
    idx = np.stack([uk // (KB*KB), (uk // KB) % KB, uk % KB], 1)
    N = fn[fid].copy() if flat else np.concatenate(nrms)[o].astype(np.float64)
    if blend:                                                                   # only the cubes holding two or more triangles need the blend: there each
        ok = np.linalg.norm(tot, axis=1) > 1e-12; N[row[ok]] = tot[ok]        # triangle was clipped to the cube and weighed by the area it covers inside it;
        if log: log(f'  blend: {len(row)} of {len(uk)} cubes hold two or more triangles ({npairs} pairs clipped on {workers} cores)')   # a cube with no area inside it (triangles only touching it) keeps the nearest facing
    N /= np.linalg.norm(N, axis=1, keepdims=True)+1e-12
    k = len(idx)
    if log: log(f'grid {cells} cells/height, cell {size:.3f} mm: {k} surface cells, {ntest/1e6:.0f}M cube-triangle tests, {len(tasks)} batches on {workers} cores, {time.time()-t0:.0f} s')
    z = np.zeros(k)
    return dict(idx=idx, N=N, turn=z, line=np.zeros(k, bool), face=face, size=size, lo=lo, d1=N.copy())
