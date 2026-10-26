"""Stage 1, dots: where the surface folds.

Cubes of one size are laid over the model; a cube is a surface cube when a triangle passes through it (exact
triangle/cube overlap test, no sample points). Each cube's arrow is the facing of the flat triangles inside it, the
surface the print has (no corner blending): with blend on, the facings of all the triangles in the cube weighted by the
area each covers inside it, so the arrow changes gradually as a crease crosses from one cube to the next; with blend
off, the facing of the triangle nearest the cube centre alone (cheaper: no clipping). That is the field, computed once.

The strength of a cube is read afterwards from the cubes around it, each compared with the cube itself (designed
2026-10-24 night with the user, after the pairwise variety x closeness reading demoted a real crease below a button):
  for every surface cube within the reach (all of them, not 26 rays), the angle between its arrow and the cube's
  arrow, a. A cube facing away (the far side of a thin wall) is left out.
  distance weight: 1 next door, fading in a straight line to nothing at one step past the reach, so near cubes carry
             the reading: a bend concentrated in one step reads its full angle, the same bend spread over the reach
             reads less, and a bump that bends evenly in every direction reads its small near-neighbour difference.
  difference weight: a^steep, steep a handle (0: every cube counts the same, a plain weighted mean; 1: a cube counts
             as much as it differs, so the cubes across a crease decide and the many that agree weigh nothing).
  strength = sum(w_dist * a^steep * a) / sum(w_dist * a^steep), in degrees.
The reach is the handle for soft bends: a bend spread over several cubes reads its full angle when the arrows compared
sit that far apart.
Concentration (2026-10-24, agreed with the user, built last): a curved surface is never a crease. By definition a crease
holds its turn in one place and a curved surface spreads it, so the reading is the cube's own turn against the average
turn of the surface cubes around it. own turn = the through turn: the angle between the arrows of the cube's two
opposite-side neighbours (first side: the next-door cube whose arrow differs most from the cube's; second side: the
existing next-door cube most nearly opposite the first in position), the cube's own arrow unused. On the box it reads
90 at the crest, 45 on the cube beside it, 0 beyond; where a slope cube touches the far face at a corner it reads 90
too, so the full reading sits on one to three cubes. (The first build used the largest next-door difference, which by
the same box reads 90 on the cube beside the crest and 45 on the crest itself: wrong by derivation.) region average = the fade-weighted mean of own turn over the surface
cubes within the reach, the cube itself included. concentration = own / average, a ratio: 1 where every cube turns
alike (a button, fur, a plane), several at a crest whose region is mostly flank. It is silent on a crease whose turn is
wider than the reach (the region is all turn) and dumbed down when the reach spans several features; on flat surface
both numbers are noise, so the strength cut stays underneath. Two earlier forms failed in the doing, not the idea:
contrast (an averaged turn against a ring) and a concentration read along a straight grid row that leaves the surface
after one cube. A dot needs concentration >= the concentration cut (panel slider, 0 = off). Contrast stays in the app
as its own slider (strength against the ring beyond the reach), computed only once its cut is raised.
Spread (2026-10-24, last): the centre-compared strength cannot peak on the crest (the crest's arrow is the middle of the
turn, its flanks see the whole fold) and its band is as wide as the reach. Spread compares the region's cubes with each
other: every pair of surface cubes within the reach (the cube itself among them), pair weight = the product of the two
distance fades, difference = 1 - cos(angle), weighted by difference^steep (steep 0, 1 or 2), the weighted mean turned
back into degrees (arccos(1 - mean)). A pair that straddles the crest is near on both sides at the crest and far on one
side at a flank, so the peak sits on the crest; a 90-degree edge reads from the pair across it even where no straight
line through the crest stays on the surface. Computed without touching any pair: the sums of pair cosines to the power
q come from the sums of w*n, w*n*n and w*n*n*n over the members (the "shortcut": small differences count as their
square, about difference weight 1 built in, which the user accepted). A dot needs spread >= the spread cut as well.
Line (2026-10-24, night, the user's idea): a size-free form of direction. At a crease every arrow is tilted about one
axis, the crease line, so all the arrows are at right angles to that line: the arrow cloud lies along one great circle.
On a button the tilt axes point every way and no line is at right angles to all the arrows. line = 1 - (short width /
long width) of the cloud, 0 a blob, 1 a perfect line, the size of the fold divided out, so a subtle crease scores as a
deep one ("lifting" rather than filtering). Known by derivation: any cylinder-shaped surface (a log, a fur stroke, a
rounded limb) scores as a line over its whole surface, and where there is no bend the shape is noise, so the strength
cut stays underneath as a floor. From the same covariance as direction, free.
Lift (same night): the user's preferred use of line is a boost, not a gate: strength = raw strength x (1 + lift x line),
lift a panel handle (0 = off). A perfect line at lift 1 doubles its strength, a blob keeps its own, so a subtle crease
clears the same cut as a deep one. The cut, the crest share on strength, contrast and the saved files see the lifted
strength; the raw one is kept as turn_raw.
Crest test (2026-10-24, last): every reading lights the whole slope, since any cube within reach sees the crease, so a
cut cannot thin a band into a line. crest(field, value, reach): for each cube, the across direction is the next-door
neighbour whose arrow differs most from the cube's; the other side is the existing next-door cube lying most nearly
opposite to it. The cube is a crest where no cube within the reach along either side reads higher (ties kept). With
the test on, a dot must be a crest of the reading the cuts are judged by (spread when its cut is on, else strength).
Crest share (2026-10-24, last): the strict test kept one cube and broke lines where two or three cubes share the top.
crest_share(): each cube's reading as a share of the highest reading within reach on its across line (itself included),
and the panel handle "Crest (%)" keeps cubes at or above that share; 100 is the strict test, 0 is off.
Speed: readings() builds strength, spread and the across direction from one pass over the neighbour lookups (they are
the cost), the crest walk reuses the direction, and contrast is only computed when its cut is on.
Direction (2026-10-24, last): what a crease is, as opposed to a bump: the facing turns across one direction and not
along it. The arrows within the reach (fade weights, the cube itself at 1) form a cloud on the unit sphere: a point on
a plane, a short line at a crease, a round blob on a bump. Its widths are the square roots of the eigenvalues of the
weighted covariance of the arrows, as angles (asin); direction = widest width minus the width at right angles, in
degrees. Crease high, bump near zero whatever its size, plane zero. Measured on the user's dump at reach 2: pocket edge
17, button 0.9, fur strokes 3, flat 0.2. The covariance comes from the same sums as spread, so it is free. A dot needs
direction >= the direction cut as well (panel slider, 0 = off).
Crest test (2026-10-24, last): every reading lights the whole slope, since any cube within reach sees the crease, so a
cut cannot thin a band into a line. crest(field, value, reach): for each cube, the across direction is the next-door
neighbour whose arrow differs most from the cube's; the other side is the existing next-door cube lying most nearly
opposite to it. The cube is a crest where no cube within the reach along either side reads higher (ties kept). With
the test on, a dot must be a crest of the reading the cuts are judged by (spread when its cut is on, else strength).
Crest share (2026-10-24, last): the strict test kept one cube and broke lines where two or three cubes share the top.
crest_share(): each cube's reading as a share of the highest reading within reach on its across line (itself included),
and the panel handle "Crest (%)" keeps cubes at or above that share; 100 is the strict test, 0 is off.
Speed: readings() builds strength, spread and the across direction from one pass over the neighbour lookups (they are
the cost), the crest walk reuses the direction, and contrast is only computed when its cut is on.
Direction (2026-10-24, last): what a crease is, as opposed to a bump: the facing turns across one direction and not
along it. The arrows within the reach (fade weights, the cube itself at 1) form a cloud on the unit sphere: a point on
a plane, a short line at a crease, a round blob on a bump. Its widths are the square roots of the eigenvalues of the
weighted covariance of the arrows, as angles (asin); direction = widest width minus the width at right angles, in
degrees. Crease high, bump near zero whatever its size, plane zero. Measured on the user's dump at reach 2: pocket edge
17, button 0.9, fur strokes 3, flat 0.2. The covariance comes from the same sums as spread, so it is free. A dot needs
direction >= the direction cut as well (panel slider, 0 = off).
Concentration (2026-10-24, last): the one property a curved surface never has. Every reading above measures how much
the surface turns, so a button sits on the same scale as the creases and a cut that drops it drops creases too.
Concentration measures how the turn is distributed along the across line: the turn at the cube's own step over the
sum of the turns at every step out to one step past the reach on both sides. A sphere spreads its turn evenly (1/(r+1)), a
crease holds nearly all of it in one step (near 1) whatever its size, so a 5-degree crease and a 90-degree edge read
alike and a button never reads high. A crest on a cube wall splits it between two cubes. Panel slider
"Concentration (%)", 0 = off; computed with the other readings (a few lookups per cube). Earlier readings, each measured on the chimp and dropped (v2/README.md): a fold angle
on a straight line through the cube (read zero on 83 % of sharp edges), the self-weighted mean of the differences to
the 26 next-door cubes (band three wide, crest not the peak), pairwise variety x closeness over 26 rays at the reach
(closeness demoted a crest on a cube wall below a sphere; the rays found 6 cubes of 26 at reach 3).

field(P, F, cube_mm, blend, log) -> dict(idx, N, lo, size): every surface cube's grid index and arrow, the grid origin
and cube size (mm). strength(field, reach, steep) -> the strength per cube; contrast(field, turn, reach) -> the contrast per
cube; readings(field, reach, steep) -> strength, spread, direction, concentration, own turn, across direction per cube.
dots(field, cut, scut, crest_pct, dcut, conc) -> mask of the dots. Cube laying copied 2026-10-24 from voxel.py
and creases.py of the first app, trimmed to the path that ran there."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')   # the strength runs its own threads; BLAS must not add more per call
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


def strength(fld, reach=1, steep=1.0, log=None):
    """per cube: the distance- and difference-weighted mean of the angle between its arrow and the arrows of the surface
    cubes within the reach (see the module header). Returns float32 per cube. One neighbour offset at a time, the offsets
    shared out over all cores; facing-away cubes count nothing."""
    from concurrent.futures import ThreadPoolExecutor
    idx = fld['idx'].astype(np.int64); N = fld['N'].astype(np.float32); k = len(idx); t0 = time.time(); r = int(reach); p = float(steep)
    key = (idx[:, 0]*KB+idx[:, 1])*KB+idx[:, 2]; srt = np.argsort(key, kind='stable'); skey = key[srt]
    offs = np.array([[i, j, l] for i in range(-r, r+1) for j in range(-r, r+1) for l in range(-r, r+1) if (i, j, l) != (0, 0, 0)])
    dist = np.linalg.norm(offs, axis=1); w = (r+1-dist)/r; keep = w > 0; offs, w = offs[keep], w[keep]           # 1 next door, 0 one step past the reach

    def part(q):
        num = np.zeros(k, np.float64); den = np.zeros(k, np.float64)
        for o, wd in zip(offs[q::nthreads], w[q::nthreads]):
            kk = ((idx[:, 0]+o[0])*KB+(idx[:, 1]+o[1]))*KB+(idx[:, 2]+o[2])
            pos = np.minimum(np.searchsorted(skey, kk), k-1); j = np.where(skey[pos] == kk, srt[pos], -1); h = j >= 0
            d = (N[h]*N[j[h]]).sum(1); use = d >= 0; hh = np.where(h)[0][use]
            a = np.degrees(np.arccos(np.clip(d[use], -1, 1))).astype(np.float64); ap = np.power(a, p) if p else np.ones_like(a)
            num[hh] += wd*ap*a; den[hh] += wd*ap
        return num, den
    nthreads = min(len(offs), os.cpu_count() or 1); num = np.zeros(k); den = np.zeros(k)
    with ThreadPoolExecutor(nthreads) as ex:
        for n_, d_ in ex.map(part, range(nthreads)): num += n_; den += d_
    out = np.where(den > 0, num/np.maximum(den, 1e-12), 0.0).astype(np.float32)
    if log: log(f'strength at reach {r}, steepness {p:g} ({len(offs)} offsets): p50 {np.percentile(out, 50):.1f} p90 {np.percentile(out, 90):.1f} p99 {np.percentile(out, 99):.1f}, {time.time()-t0:.1f} s')
    return out


def field(P, F, cube_mm, blend=True, log=None, workers=None):
    """The field: cubes of cube_mm over the mesh, each with its arrow. Reach and cut are not part of it (strength(), dots()).
    Returns dict(idx, N, lo, size)."""
    height = float((P.max(0)-P.min(0)).max()); cells = max(1, int(round(height/cube_mm)))
    if cells >= KB: raise ValueError(f'{cells} cubes along the height of {height:.2f} mm: the grid index holds fewer than {KB} per axis, so the cube must be at least {height/(KB-1):.5f} mm')
    r = cubes(P, F, cells, log, height=height, blend=blend, workers=workers)
    return dict(idx=r['idx'].astype(np.int32), N=r['N'].astype(np.float32), lo=r['lo'], size=float(r['size']))


def spread(fld, reach=1, steep=1.0, log=None):
    """per cube: the pairwise reading (module header). Members: the cube (weight 1) and the surface cubes within the reach
    (weight fading from 1 next door to 0 one step past the reach; facing-away cubes left out). Over all pairs i != j:
    mean of d_ij weighted by w_i*w_j*d_ij^p, d = 1 - cos(angle), p = round(steep) capped at 2; result arccos(1 - mean) in
    degrees. The pair sums come from moment tensors of the members, so the cost is one pass over the offsets."""
    from concurrent.futures import ThreadPoolExecutor
    from math import comb
    idx = fld['idx'].astype(np.int64); N = fld['N'].astype(np.float64); k = len(idx); t0 = time.time(); r = int(reach); p = int(min(2, max(0, round(steep))))
    key = (idx[:, 0]*KB+idx[:, 1])*KB+idx[:, 2]; srt = np.argsort(key, kind='stable'); skey = key[srt]
    offs = np.array([[i, j, l] for i in range(-r, r+1) for j in range(-r, r+1) for l in range(-r, r+1) if (i, j, l) != (0, 0, 0)])
    dist = np.linalg.norm(offs, axis=1); w = (r+1-dist)/r; keep = w > 0; offs, w = offs[keep], w[keep]
    # symmetric tensor slots: order-2 (6 unique with multiplicity) and order-3 (10 unique), so ||T||^2 = sum mult * T_slot^2
    s2 = [(0, 0, 1), (1, 1, 1), (2, 2, 1), (0, 1, 2), (0, 2, 2), (1, 2, 2)]
    s3 = [(0, 0, 0, 1), (1, 1, 1, 1), (2, 2, 2, 1), (0, 0, 1, 3), (0, 0, 2, 3), (1, 1, 0, 3), (1, 1, 2, 3), (2, 2, 0, 3), (2, 2, 1, 3), (0, 1, 2, 6)]
    need3 = p >= 2

    def moments(n, wv):
        """add weighted members: n (m, 3) arrows, wv (m,) weights -> S0, S0sq, V (3), T2 (6), T3 (10)"""
        out = [wv, wv*wv, wv[:, None]*n, np.stack([wv*n[:, a]*n[:, b] for a, b, _ in s2], 1)]
        out.append(np.stack([wv*n[:, a]*n[:, b]*n[:, c] for a, b, c, _ in s3], 1) if need3 else None); return out
    def part(q):
        S0 = np.ones(k); S0sq = np.ones(k); V = N.copy(); T2 = np.stack([N[:, a]*N[:, b] for a, b, _ in s2], 1); T3 = np.stack([N[:, a]*N[:, b]*N[:, c] for a, b, c, _ in s3], 1) if need3 else None
        if q: S0[:] = 0; S0sq[:] = 0; V[:] = 0; T2[:] = 0; T3 is not None and T3.fill(0)      # the centre member is added once, by thread 0
        for o, wd in zip(offs[q::nthreads], w[q::nthreads]):
            kk = ((idx[:, 0]+o[0])*KB+(idx[:, 1]+o[1]))*KB+(idx[:, 2]+o[2])
            pos = np.minimum(np.searchsorted(skey, kk), k-1); j = np.where(skey[pos] == kk, srt[pos], -1); h = j >= 0
            use = (N[h]*N[j[h]]).sum(1) >= 0; hh = np.where(h)[0][use]; n = N[j[hh]]
            S0[hh] += wd; S0sq[hh] += wd*wd; V[hh] += wd*n; T2[hh] += wd*np.stack([n[:, a]*n[:, b] for a, b, _ in s2], 1)
            if need3: T3[hh] += wd*np.stack([n[:, a]*n[:, b]*n[:, c] for a, b, c, _ in s3], 1)
        return S0, S0sq, V, T2, T3
    nthreads = min(len(offs), os.cpu_count() or 1)
    S0 = np.zeros(k); S0sq = np.zeros(k); V = np.zeros((k, 3)); T2 = np.zeros((k, 6)); T3 = np.zeros((k, 10)) if need3 else None
    with ThreadPoolExecutor(nthreads) as ex:
        for a_, b_, c_, d_, e_ in ex.map(part, range(nthreads)):
            S0 += a_; S0sq += b_; V += c_; T2 += d_
            if need3: T3 += e_
    # M_q = sum over ordered pairs i != j of w_i w_j cos^q = (full double sum) - (diagonal: sum w_i^2, cos = 1)
    M = [S0*S0-S0sq, (V*V).sum(1)-S0sq, (T2*T2*np.array([m for *_, m in s2])).sum(1)-S0sq]
    if need3: M.append((T3*T3*np.array([m for *_, m in s3])).sum(1)-S0sq)
    def delta_pow(q):   # sum w_i w_j (1 - cos)^q = sum_m C(q, m) (-1)^m M_m
        return sum(comb(q, m)*(-1)**m*M[m] for m in range(q+1))
    num = delta_pow(p+1); den = delta_pow(p); mean = np.where(den > 1e-12, num/np.maximum(den, 1e-12), 0.0)
    out = np.degrees(np.arccos(np.clip(1-mean, -1, 1))).astype(np.float32)
    if log: log(f'spread at reach {r}, steepness {p} ({len(offs)} offsets): p50 {np.percentile(out, 50):.1f} p90 {np.percentile(out, 90):.1f} p99 {np.percentile(out, 99):.1f}, {time.time()-t0:.1f} s')
    return out


def _eig3(C):
    """eigenvalues of symmetric 3x3 matrices (n, 3, 3), descending, by the trigonometric closed form (vectorised; LAPACK
    per matrix is 50x slower on millions of them)"""
    q = np.trace(C, axis1=1, axis2=2)/3; B = C-q[:, None, None]*np.eye(3)
    p = np.sqrt(np.maximum((B*B).sum((1, 2))/6, 0)); safe = np.where(p > 1e-15, p, 1.0)
    Bn = B/safe[:, None, None]; r = np.linalg.det(Bn)/2; phi = np.arccos(np.clip(r, -1, 1))/3
    l1 = q+2*p*np.cos(phi); l3 = q+2*p*np.cos(phi+2*np.pi/3); l2 = 3*q-l1-l3
    out = np.stack([l1, l2, l3], 1); out[p <= 1e-15] = q[p <= 1e-15, None]; return out


def readings(fld, reach=1, steep=1.0, log=None):
    """one pass over the neighbour lookups within the reach (plus a pre-pass next door) -> dict(turn, spread, direction,
    concentration, own, across, far). turn and spread equal strength() and spread(); direction is the widest width of the
    arrow cloud minus the width at right angles (degrees); own is the cube's through turn (the angle between its two
    opposite-side neighbours' arrows, degrees) and concentration its ratio to the fade-weighted mean own turn of the
    region (module header);
    across is the index (0..25) of the next-door direction whose arrow differs most from the cube's (-1 when none
    differs), far the index of the existing next-door direction most nearly opposite (-1 if none)."""
    from concurrent.futures import ThreadPoolExecutor
    from math import comb
    idx = fld['idx'].astype(np.int64); N = fld['N'].astype(np.float64); k = len(idx); t0 = time.time(); r = int(reach); p = float(steep); pi = int(min(2, max(0, round(steep))))
    key = (idx[:, 0]*KB+idx[:, 1])*KB+idx[:, 2]; srt = np.argsort(key, kind='stable'); skey = key[srt]
    offs = np.array([[i, j, l] for i in range(-r, r+1) for j in range(-r, r+1) for l in range(-r, r+1) if (i, j, l) != (0, 0, 0)])
    dist = np.linalg.norm(offs, axis=1); w = (r+1-dist)/r; keep = w > 0; offs, w = offs[keep], w[keep]
    s2 = [(0, 0, 1), (1, 1, 1), (2, 2, 1), (0, 1, 2), (0, 2, 2), (1, 2, 2)]
    s3 = [(0, 0, 0, 1), (1, 1, 1, 1), (2, 2, 2, 1), (0, 0, 1, 3), (0, 0, 2, 3), (1, 1, 0, 3), (1, 1, 2, 3), (2, 2, 0, 3), (2, 2, 1, 3), (0, 1, 2, 6)]
    need3 = pi >= 2
    # pre-pass over the 26 next-door offsets: the two opposite-side neighbours and the through turn (the box derivation)
    nb = np.full((k, 26), -1, np.int64); best = np.zeros(k); barg = np.full(k, -1, np.int32)
    for q, o in enumerate(OFFS):
        kk = ((idx[:, 0]+o[0])*KB+(idx[:, 1]+o[1]))*KB+(idx[:, 2]+o[2]); pos = np.minimum(np.searchsorted(skey, kk), k-1); j = np.where(skey[pos] == kk, srt[pos], -1); h = j >= 0
        d = np.ones(k); d[h] = (N[h]*N[j[h]]).sum(1); ok = h & (d >= 0); nb[ok, q] = j[ok]
        dd = np.where(ok, 1-d, -1.0); upd = dd > best; best[upd] = dd[upd]; barg[upd] = q
    exist = nb >= 0; unit = OFFS/np.linalg.norm(OFFS, axis=1)[:, None]; cosd = unit @ unit.T
    cand = cosd[np.maximum(barg, 0)]; cand[~exist] = 2.0; far = np.where((barg >= 0) & (cand.min(1) < 0), cand.argmin(1), -1).astype(np.int32)
    across = np.where(far >= 0, barg, -1).astype(np.int32); okt = across >= 0
    ja = nb[np.arange(k), np.maximum(across, 0)]; jf = nb[np.arange(k), np.maximum(far, 0)]
    own = np.zeros(k); own[okt] = np.degrees(np.arccos(np.clip((N[ja[okt]]*N[jf[okt]]).sum(1), -1, 1)))   # the through turn: the two sides against each other, the cube's own arrow unused
    del nb

    def part(q):
        num = np.zeros(k); den = np.zeros(k); S0 = np.zeros(k); S0sq = np.zeros(k); V = np.zeros((k, 3)); T2 = np.zeros((k, 6)); T3 = np.zeros((k, 10)) if need3 else None; aown = np.zeros(k)
        for o, wd in zip(offs[q::nthreads], w[q::nthreads]):
            kk = ((idx[:, 0]+o[0])*KB+(idx[:, 1]+o[1]))*KB+(idx[:, 2]+o[2])
            pos = np.minimum(np.searchsorted(skey, kk), k-1); j = np.where(skey[pos] == kk, srt[pos], -1); h = j >= 0
            d = (N[h]*N[j[h]]).sum(1); use = d >= 0; hh = np.where(h)[0][use]; d = d[use]; n = N[j[hh]]
            a = np.degrees(np.arccos(np.clip(d, -1, 1))); ap = np.power(a, p) if p else np.ones_like(a); num[hh] += wd*ap*a; den[hh] += wd*ap
            S0[hh] += wd; S0sq[hh] += wd*wd; V[hh] += wd*n; T2[hh] += wd*np.stack([n[:, x]*n[:, y] for x, y, _ in s2], 1); aown[hh] += wd*own[j[hh]]
            if need3: T3[hh] += wd*np.stack([n[:, x]*n[:, y]*n[:, z] for x, y, z, _ in s3], 1)
        return num, den, S0, S0sq, V, T2, T3, aown
    nthreads = min(len(offs), os.cpu_count() or 1)
    num = np.zeros(k); den = np.zeros(k); S0 = np.ones(k); S0sq = np.ones(k); V = N.copy(); T2 = np.stack([N[:, x]*N[:, y] for x, y, _ in s2], 1); T3 = np.stack([N[:, x]*N[:, y]*N[:, z] for x, y, z, _ in s3], 1) if need3 else None
    aown = own.copy()   # the cube itself, weight 1
    with ThreadPoolExecutor(nthreads) as ex:
        for n_, d_, a_, b_, v_, t_, u_, ao_ in ex.map(part, range(nthreads)):
            num += n_; den += d_; S0 += a_; S0sq += b_; V += v_; T2 += t_; aown += ao_
            if need3: T3 += u_
    turn = np.where(den > 0, num/np.maximum(den, 1e-12), 0.0).astype(np.float32)
    M = [S0*S0-S0sq, (V*V).sum(1)-S0sq, (T2*T2*np.array([m for *_, m in s2])).sum(1)-S0sq]
    if need3: M.append((T3*T3*np.array([m for *_, m in s3])).sum(1)-S0sq)
    dp = lambda q: sum(comb(q, m)*(-1)**m*M[m] for m in range(q+1))
    mean = np.where(dp(pi) > 1e-12, dp(pi+1)/np.maximum(dp(pi), 1e-12), 0.0); spread = np.degrees(np.arccos(np.clip(1-mean, -1, 1))).astype(np.float32)
    # direction: weighted covariance of the arrows from the same sums (T2 holds the 6 unique entries of sum w n n^T)
    mu = V/np.maximum(S0, 1e-12)[:, None]; C = np.empty((k, 3, 3))
    for (x, y, _), col in zip(s2, range(6)): C[:, x, y] = T2[:, col]/np.maximum(S0, 1e-12)-mu[:, x]*mu[:, y]; C[:, y, x] = C[:, x, y]
    lam = _eig3(C); width = np.degrees(np.arcsin(np.sqrt(np.clip(lam[:, :2], 0, 1)))); direction = (width[:, 0]-width[:, 1]).astype(np.float32)
    line = np.where(width[:, 0] > 1e-9, 1-width[:, 1]/np.maximum(width[:, 0], 1e-9), 0.0).astype(np.float32)   # size-free: short width over long width
    avg = aown/np.maximum(S0, 1e-12); concentration = np.where(avg > 1e-9, own/np.maximum(avg, 1e-9), 0.0).astype(np.float32)
    if log: log(f'readings at reach {r}, weight {p:g} ({len(offs)} offsets): strength p50 {np.percentile(turn, 50):.1f} p90 {np.percentile(turn, 90):.1f}; spread p50 {np.percentile(spread, 50):.1f} p90 {np.percentile(spread, 90):.1f}; direction p50 {np.percentile(direction, 50):.1f} p90 {np.percentile(direction, 90):.1f}; concentration p50 {np.percentile(concentration, 50):.2f} p90 {np.percentile(concentration, 90):.2f}; line p50 {np.percentile(line, 50):.2f} p90 {np.percentile(line, 90):.2f}; {time.time()-t0:.1f} s')
    return dict(turn=turn, spread=spread, direction=direction, line=line, concentration=concentration, own=own.astype(np.float32), across=across, far=far)


def contrast(fld, turn, reach=1, log=None):
    """per cube: its strength minus the mean strength of the surface cubes in the ring just beyond the reach (distance over
    the reach up to twice the reach, weight 1 at the inner edge fading to 0 at the outer; facing-away cubes left out). The
    ring starts where the crest's own band ends, since strength spreads exactly as far as the reach. Zero on uniform bend,
    positive on a crest."""
    from concurrent.futures import ThreadPoolExecutor
    idx = fld['idx'].astype(np.int64); N = fld['N'].astype(np.float32); turn = np.asarray(turn, np.float64); k = len(idx); t0 = time.time(); r = int(reach)
    key = (idx[:, 0]*KB+idx[:, 1])*KB+idx[:, 2]; srt = np.argsort(key, kind='stable'); skey = key[srt]
    R2 = 2*r+1; offs = np.array([[i, j, l] for i in range(-R2, R2+1) for j in range(-R2, R2+1) for l in range(-R2, R2+1) if (i, j, l) != (0, 0, 0)])
    dist = np.linalg.norm(offs, axis=1); w = (R2-dist)/r; keep = (dist > r) & (w > 0); offs, w = offs[keep], w[keep]   # the ring: past the reach, fading to 0 at 2r+1

    def part(q):
        num = np.zeros(k); den = np.zeros(k)
        for o, wd in zip(offs[q::nthreads], w[q::nthreads]):
            kk = ((idx[:, 0]+o[0])*KB+(idx[:, 1]+o[1]))*KB+(idx[:, 2]+o[2])
            pos = np.minimum(np.searchsorted(skey, kk), k-1); j = np.where(skey[pos] == kk, srt[pos], -1); h = j >= 0
            use = (N[h]*N[j[h]]).sum(1) >= 0; hh = np.where(h)[0][use]; num[hh] += wd*turn[j[hh]]; den[hh] += wd
        return num, den
    nthreads = min(len(offs), os.cpu_count() or 1); num = np.zeros(k); den = np.zeros(k)
    with ThreadPoolExecutor(nthreads) as ex:
        for n_, d_ in ex.map(part, range(nthreads)): num += n_; den += d_
    around = np.where(den > 0, num/np.maximum(den, 1e-12), 0.0); out = (turn-around).astype(np.float32)
    if log: log(f'contrast at reach {r}: p50 {np.percentile(out, 50):.1f} p90 {np.percentile(out, 90):.1f} p99 {np.percentile(out, 99):.1f}, {time.time()-t0:.1f} s')
    return out



def crest_share(fld, value, across, far, reach=1, log=None):
    """per cube: its value as a share (0..1) of the highest value within the reach on its across line, itself included;
    0 where the cube has no across direction. 1 is the strict crest."""
    idx = fld['idx'].astype(np.int64); value = np.asarray(value, np.float64); k = len(idx); t0 = time.time(); r = int(reach)
    key = (idx[:, 0]*KB+idx[:, 1])*KB+idx[:, 2]; srt = np.argsort(key, kind='stable'); skey = key[srt]
    def at(ijk):
        kk = (ijk[:, 0]*KB+ijk[:, 1])*KB+ijk[:, 2]; pos = np.minimum(np.searchsorted(skey, kk), k-1); return np.where(skey[pos] == kk, srt[pos], -1)
    top = value.copy(); ok0 = across >= 0
    for side in (across, far):
        alive = ok0.copy(); o = OFFS[np.maximum(side, 0)]
        for step in range(1, r+1):
            j = at(idx+step*o); alive &= j >= 0; sel = np.where(alive)[0]; top[sel] = np.maximum(top[sel], value[j[sel]])
    out = np.where(ok0 & (top > 0), value/np.maximum(top, 1e-12), 0.0).astype(np.float32)
    if log: log(f'crest share at reach {r}: {int((out >= 1).sum())} strict crests of {k}, {time.time()-t0:.1f} s')
    return out


def dots(fld, cut, scut=0.0, crest_pct=0.0, dcut=0.0, conc=0.0, ccut=0.0, lcut=0.0):
    """the dots of a field: strength at least the cut, spread at least the spread cut, direction at least the direction
    cut, concentration at least conc (a ratio), contrast at least the contrast cut, crest share at least crest_pct/100
    (of spread when its cut is on, else of strength); a cut of 0 turns its filter off"""
    m = fld['turn'] >= cut
    if ccut > 0 and fld.get('contrast') is not None: m &= fld['contrast'] >= ccut
    if scut > 0 and fld.get('spread') is not None: m &= fld['spread'] >= scut
    if dcut > 0 and fld.get('direction') is not None: m &= fld['direction'] >= dcut
    if conc > 0 and fld.get('concentration') is not None: m &= fld['concentration'] >= conc
    if lcut > 0 and fld.get('line') is not None: m &= fld['line'] >= lcut
    if crest_pct > 0:
        c = fld.get('crest_spread' if scut > 0 else 'crest_turn')
        if c is not None: m &= c >= crest_pct/100-1e-6
    return m
