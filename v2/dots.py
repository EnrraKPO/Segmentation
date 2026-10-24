"""Stage 1, dots: where the surface folds.

Cubes of one size are laid over the model; a cube is a surface cube when a triangle passes through it (exact
triangle/cube overlap test, no sample points). Each cube's arrow is the facing of the flat triangles inside it, the
surface the print has (no corner blending): with blend on, the facings of all the triangles in the cube weighted by the
area each covers inside it, so the arrow changes gradually as a crease crosses from one cube to the next; with blend
off, the facing of the triangle nearest the cube centre alone (cheaper: no clipping). Per cube the fold angle is read
from its 26 neighbours at the reach (reach 1: the cubes next door, reach r: the cubes r steps away in those directions):
the difference between the cube's arrow and each neighbour's, averaged with each difference weighted by itself, in
degrees. A neighbour that matches the cube's arrow, faces away (the far side of a thin wall) or does not exist weighs
nothing, so one neighbour across the crease is enough to read the whole fold. The reach is the handle for soft bends: a
bend spread over several cubes reads its full angle when the arrows compared sit that far apart. A cube is a dot when
its fold angle is at least the cut; the cut is applied afterwards, as a filter on the field, so it can be moved without
recomputing anything. (Until 2026-10-24 the fold angle was read between the cube one behind and one ahead on a straight
line through the cube, which at a crease leaves the surface on one side and read zero on 83 % of sharp edges.)

field(P, F, cube_mm, blend, log, reach) -> dict(idx, N, turn, lo, size): every surface cube's grid index, arrow and fold angle,
the grid origin and cube size (mm). dots(field, cut) -> mask of the dots. Copied 2026-10-24 from voxel.py and creases.py
of the first app, trimmed to the path that ran there (flat facing always; the sample-point grid and the thick-stroke
method are gone)."""
import os, time, numpy as np

KB = 1 << 16   # index packing base (cubes per axis must stay below this)
OFFS = np.array([[i, j, k] for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1) if (i, j, k) != (0, 0, 0)])


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
        shm = shared_memory.SharedMemory(name=name); _W['_shm'].append(shm)
        _W[k] = np.ndarray(shape, dtype=dtype, buffer=shm.buf)
    _W['lo'] = lo; _W['size'] = size


def _direct_task(args):
    """one batch of triangles that all span the same box of cubes: returns the cubes they overlap"""
    tri, sp, imin_tri, blend = args
    A, B, C, lo, size = (_W[k] for k in ('A', 'B', 'C', 'lo', 'size'))
    ncand = int(np.prod(sp)); half = size*0.5
    offs = np.stack(np.meshgrid(np.arange(sp[0]), np.arange(sp[1]), np.arange(sp[2]), indexing='ij'), -1).reshape(-1, 3)
    cid = imin_tri[:, None, :]+offs[None]                                     # (n, ncand, 3) cube indices
    Q = (lo+(cid+0.5)*size).reshape(-1, 3); rep = np.repeat(tri, ncand)
    inside = _tri_box_overlap(Q, A[rep], B[rep], C[rep], half); ntest = len(Q)
    if not inside.any(): return None, None, None, ntest, None
    Q = Q[inside]; rep = rep[inside]; cid = cid.reshape(-1, 3)[inside]
    X = _closest_on_triangles(Q, A[rep], B[rep], C[rep]); d = np.linalg.norm(X-Q, axis=1)
    key = (cid[:, 0]*KB+cid[:, 1])*KB+cid[:, 2]
    # side score for ties (cube centre equally near two triangles = nearest point on their shared edge): the triangle
    # whose facing points most directly at the cube centre is the side the cube is over; lower score = better
    fnr = _W['FN'][rep]; side = (-np.abs((fnr*(Q-X)).sum(1))/np.maximum(d, 1e-12)).astype(np.float32)
    o = np.lexsort((side, d, key))
    if not blend: first = np.r_[True, key[o][1:] != key[o][:-1]]; o = o[first]   # one hit per cube: the nearest (blend: every hit is kept, the merge needs all the pairs)
    key, d, rep, X, side = key[o], d[o], rep[o], X[o], side[o]
    return key, d.astype(np.float32), rep.astype(np.int32), ntest, side


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


def cubes(P, F, cells, log=None, height=None, batch=2_000_000, workers=None, box=None, blend=True):
    """Surface cubes straight from the triangles. For each triangle, every cube in the box its corners span is tested
    exactly for overlap with the triangle (separating axes). A cube's `face` is the triangle whose nearest point is closest
    to the cube centre. Its arrow N: with blend, the facings of ALL the triangles inside the cube, each weighted by the
    area it covers inside (the triangle clipped to the cube); without, the facing of `face` alone. box=(lo, hi) in mm:
    cubes are laid only inside it. The triangles are split into batches spread over `workers` processes (default: all
    cores). Returns dict(idx, N, face, size, lo)."""
    from concurrent.futures import ProcessPoolExecutor
    lo = P.min(0); ext = P.max(0)-lo; height = ext.max() if height is None else height; size = height/cells
    A, B, C = P[F[:, 0]], P[F[:, 1]], P[F[:, 2]]; tri_id = np.arange(len(F))
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
            tri = grp[c0:c0+step]; tasks.append((tri, sp, imin[tri], blend))
    tasks.sort(key=lambda t: -len(t[0])*int(np.prod(t[1])))                     # biggest batches first: balanced finish
    workers = workers or os.cpu_count() or 1
    keys, dists, fids, sides = [], [], [], []; ntest = 0; done = 0; t0 = time.time(); nxt = 0.1
    def take(r):
        nonlocal ntest, done, nxt
        key, d, rep, nt, side = r; ntest += nt; done += 1
        if key is not None: keys.append(key); dists.append(d); fids.append(rep); sides.append(side)
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
        for k, arr in dict(A=A, B=B, C=C, FN=fn).items():
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
        _W.update(A=A, B=B, C=C, FN=fn, lo=lo, size=size)
        for t in tasks: take(_direct_task(t))
        uk, fid, o, multi = merge()
        if blend: pk, pr, cu, row = multi; tot = _blend_sums(pk, pr, cu); npairs = len(pk)
    face = tri_id[fid]
    idx = np.stack([uk // (KB*KB), (uk // KB) % KB, uk % KB], 1)
    N = fn[fid].copy()
    if blend:                                                                   # only the cubes holding two or more triangles need the blend: there each
        ok = np.linalg.norm(tot, axis=1) > 1e-12; N[row[ok]] = tot[ok]        # triangle was clipped to the cube and weighed by the area it covers inside it;
        if log: log(f'  blend: {len(row)} of {len(uk)} cubes hold two or more triangles ({npairs} pairs clipped on {workers} cores)')   # a cube with no area inside it (triangles only touching it) keeps the nearest facing
    N /= np.linalg.norm(N, axis=1, keepdims=True)+1e-12
    k = len(idx)
    if log: log(f'grid {cells} cubes/height, cube {size:.3f} mm: {k} surface cubes, {ntest/1e6:.0f}M cube-triangle tests, {len(tasks)} batches on {workers} cores, {time.time()-t0:.0f} s')
    return dict(idx=idx, N=N, face=face, size=size, lo=lo)


def fold_angle(idx, N, reach=1):
    """per cube: the fold angle, in degrees: the differences between the cube's arrow and the arrows of its 26 neighbours
    at the reach, averaged with each difference as its own weight (sum of a*a over sum of a). A missing or facing-away
    neighbour contributes nothing; a cube with no difference anywhere reads 0."""
    from concurrent.futures import ThreadPoolExecutor
    k = len(idx); key = (idx[:, 0]*KB+idx[:, 1])*KB+idx[:, 2]; srt = np.argsort(key, kind='stable'); skey = key[srt]

    def look(o):
        """the arrow difference to the neighbour reach steps away in one of the 26 directions (0 where there is none)"""
        o = o*reach; kk = ((idx[:, 0]+o[0])*KB+(idx[:, 1]+o[1]))*KB+(idx[:, 2]+o[2])
        pos = np.minimum(np.searchsorted(skey, kk), k-1); j = np.where(skey[pos] == kk, srt[pos], -1)
        h = j >= 0; ang = np.zeros(k)
        dot = (N[h]*N[j[h]]).sum(1); a = np.degrees(np.arccos(np.clip(dot, -1, 1))); a[dot < 0] = 0   # facing away = far side of a thin wall
        ang[h] = a; return ang
    num = np.zeros(k); den = np.zeros(k)
    with ThreadPoolExecutor(min(26, os.cpu_count() or 1)) as ex:
        for a in ex.map(look, OFFS): num += a*a; den += a
    return np.where(den > 0, num/np.maximum(den, 1e-12), 0.0)


def field(P, F, cube_mm, blend=True, log=None, workers=None, reach=1):
    """The stage: cubes of cube_mm over the mesh, each with its arrow and fold angle read at the reach (cubes). The cut is not
    part of it: a dot is a cube whose fold angle is at least the cut, chosen afterwards (dots(field, cut)). Returns
    dict(idx, N, turn, lo, size)."""
    height = float((P.max(0)-P.min(0)).max()); cells = max(1, int(round(height/cube_mm)))
    if cells >= KB: raise ValueError(f'{cells} cubes along the height of {height:.2f} mm: the grid index holds fewer than {KB} per axis, so the cube must be at least {height/(KB-1):.5f} mm')
    r = cubes(P, F, cells, log, height=height, blend=blend, workers=workers)
    t0 = time.time(); span = fold_angle(r['idx'], r['N'], reach)
    if log: log(f'fold angle per cube at reach {reach}: p50 {np.percentile(span, 50):.1f} p90 {np.percentile(span, 90):.1f} p99 {np.percentile(span, 99):.1f} deg, {time.time()-t0:.0f} s')
    return dict(idx=r['idx'].astype(np.int32), N=r['N'].astype(np.float32), turn=span.astype(np.float32), lo=r['lo'], size=float(r['size']))


def dots(fld, cut):
    """the dots of a field at a cut: boolean mask over its cubes"""
    return fld['turn'] >= cut
