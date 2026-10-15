"""The crease pipeline as functions (the CLI scripts shift.py / trace3.py do the same step by step).

Stage 1  find_dots(P, F, cells, cut, spp, flat, log) -> res
   Cubes over the surface (voxel.run_direct: every cube a triangle passes through, exact test); each cube's arrow is
   the facing of its own triangle (flat=True, the print's surface), the triangle nearest the cube centre. Fold change per cube: along the
   neighbour direction with the largest arrow angle, |shift ahead - shift behind| in degrees. A cube is a dot when
   that change is >= cut and not smaller than the cubes behind/ahead.
Stage 2  strokes(res, W, L, K, align, agree, minlen, log) -> list of curves (mm)
   Thick stroke of half-width W mm through the dots: at each step the dots inside a window (L mm long, W mm wide) give
   the position (weighted centre) and the heading (the axis they are arranged along); stops when less than K mm of
   dot row lies ahead, on a sharp turn, or on its own tail. Dots are every cube whose turn (facing one cube behind vs
   one cube ahead, across the feature) is >= cut, weighted by turn above the cut, so the stroke sits on the centre of
   the band; a stroke does not stop at dots another stroke owns, but no stroke starts where its window is mostly owned.
Stage 3  join(curves, G, TOL, ANG): stroke ends meeting ahead (corners), facing each other (gaps) or hitting another
   stroke (tees) are joined. All lengths in mm, so the cube size only sets the resolution."""
import os, sys, time, numpy as np
from scipy.spatial import cKDTree
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import voxel

OFFS = np.array([[i, j, l] for i in (-1, 0, 1) for j in (-1, 0, 1) for l in (-1, 0, 1) if (i, j, l) != (0, 0, 0)])
KB = 1 << 16


def fold_change(idx, N, cut, log=None):
    """per cube: |shift ahead - shift behind| (deg) along the direction of the largest neighbour angle, and the crest flag"""
    from concurrent.futures import ThreadPoolExecutor
    k = len(idx); key = (idx[:, 0]*KB+idx[:, 1])*KB+idx[:, 2]; srt = np.argsort(key, kind='stable'); skey = key[srt]

    def look(o):
        """neighbour index (or -1) and the arrow angle to it, for one of the 26 directions"""
        kk = ((idx[:, 0]+o[0])*KB+(idx[:, 1]+o[1]))*KB+(idx[:, 2]+o[2])
        pos = np.minimum(np.searchsorted(skey, kk), k-1); j = np.where(skey[pos] == kk, srt[pos], -1).astype(np.int32)
        h = j >= 0; ang = np.zeros(k)
        dot = (N[h]*N[j[h]]).sum(1); a = np.degrees(np.arccos(np.clip(dot, -1, 1))); a[dot < 0] = 0   # facing away = far side of a thin wall
        ang[h] = a; return j, ang
    with ThreadPoolExecutor(min(26, os.cpu_count() or 1)) as ex: res = list(ex.map(look, OFFS))
    pmax = np.zeros(k); pdir = np.zeros(k, np.int32); nb = np.empty((k, 26), np.int32)
    for q, (j, ang) in enumerate(res):
        nb[:, q] = j; upd = ang > pmax; pmax[upd] = ang[upd]; pdir[upd] = q
    del res
    opp = np.array([np.nonzero((OFFS == -o).all(1))[0][0] for o in OFFS])
    fwd = nb[np.arange(k), pdir]; bwd = nb[np.arange(k), opp[pdir]]; ok = (fwd >= 0) & (bwd >= 0)
    near = N-N[np.maximum(bwd, 0)]; far = N[np.maximum(fwd, 0)]-N
    change = np.degrees(np.linalg.norm(far-near, axis=1)); change[~ok] = 0
    span = np.zeros(k); span[ok] = np.degrees(np.arccos(np.clip((N[fwd[ok]]*N[bwd[ok]]).sum(1), -1, 1)))   # angle between the facings one cube behind and one cube ahead
    cf = np.where(fwd >= 0, change[np.maximum(fwd, 0)], 0); cb = np.where(bwd >= 0, change[np.maximum(bwd, 0)], 0)
    crest = (change >= cut) & (change >= cf) & (change > cb)
    return change, crest, span


def find_dots(P, F, cells=1200, cut=17.0, spp=0.0, flat=True, log=None, height=None):
    """spp=0: cubes straight from the triangles (voxel.run_direct, exact overlap test, facing read at the cube centre);
    spp>0: the old sample-point sprinkle with spp samples per cube edge (slow, memory-heavy; kept for comparison)"""
    r = voxel.run_direct(P, F, cells, log, height=height, flat=flat) if spp <= 0 else voxel.run(P, F, cells, 0.0, log, height=height, spp=spp, flat=flat)
    t0 = time.time(); change, crest, span = fold_change(r['idx'], r['N'], cut)
    if log: log(f'fold change {time.time()-t0:.0f} s')
    band = span >= cut
    if log: log(f'turn per cube: p50 {np.percentile(span, 50):.1f} p90 {np.percentile(span, 90):.1f} p99 {np.percentile(span, 99):.1f} deg; dots (turn >= {cut} deg): {int(band.sum())}, of which crest {int(crest.sum())}')
    return dict(idx=r['idx'], N=r['N'], lo=r['lo'], size=r['size'], turn=span, line=band, crest=crest, change=change, d1=r['d1'], cut=cut)


def _trace(args):
    """strokes through one cluster of dots; X in cubes, W L step in cubes, K in dots. Returns curves in cube units."""
    X, wgt, W, L, K, align, agree, minlen, step, stop_foreign = args
    tree = cKDTree(X); owner = np.full(len(X), -1); R = np.hypot(L, W); cos_agree = np.cos(np.radians(agree))

    def window(p, h):
        nb = np.array(tree.query_ball_point(p, R), int)
        if len(nb) == 0: return nb, None
        off = X[nb]-p; al = off @ h; sd = np.linalg.norm(off-al[:, None]*h, axis=1)
        keep = (np.abs(al) <= L) & (sd <= W); return nb[keep], al[keep]

    def fit(p, h, tid):
        nb, al = window(p, h)
        if len(nb) < K: return None
        w = wgt[nb]; c = (X[nb]*w[:, None]).sum(0)/w.sum()
        d = X[nb]-c; ev, vec = np.linalg.eigh((d*w[:, None]).T @ d); ax = vec[:, -1]
        if ev[-1] < align*align*max(ev[-2], 1e-9): return None
        if ax @ h < 0: ax = -ax
        q = p+(c-p)-((c-p) @ h)*h
        ow = owner[nb]; return q, ax, int((al > 0).sum()), int(((ow != -1) & (ow != tid)).sum()), len(nb)

    def march(p, h, tid, claimed):
        out = []; back = int(2*L)+5
        while True:
            q = p+step*h; r = fit(q, h, tid)
            if r is None: break
            q, hq, ahead, foreign, _ = r
            if ahead < K or hq @ h < cos_agree or (stop_foreign and foreign > K): break   # foreign dots no longer stop a stroke (2026-10-15); new strokes still only seed on unowned dots
            if len(out) > back and np.linalg.norm(np.array(out[:-back])-q, axis=1).min() < W: break
            kk = tree.query_ball_point(q, W); owner[kk] = tid; claimed.extend(kk)
            out.append(q); p, h = q, hq
            if len(out) > 20000: break
        return out

    curves = []
    for s in np.argsort(-wgt):
        if owner[s] != -1: continue
        nb = np.array(tree.query_ball_point(X[s], L), int)
        if len(nb) < K: continue
        d = X[nb]-X[nb].mean(0); ev, vec = np.linalg.eigh(d.T @ d); h0 = vec[:, -1]
        r = fit(X[s], h0, len(curves))
        if r is None: continue
        if r[3] > 0.5*r[4]: continue      # seed window mostly owned by another stroke: that line is drawn already
        p0, h0 = r[0], r[1]; tid = len(curves)
        kk = tree.query_ball_point(p0, W); owner[kk] = tid; claimed = list(kk)
        ch = march(p0, -h0, tid, claimed)[::-1]+[p0]+march(p0, h0, tid, claimed)
        if len(ch) < minlen: owner[claimed] = -1; curves.append(None); continue
        curves.append(np.array(ch))
    return [c for c in curves if c is not None], int((owner == -1).sum())


def clusters(X, r):
    """label dots by connectivity: dots within r cubes of each other are in one cluster"""
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    pairs = cKDTree(X).query_pairs(r, output_type='ndarray')
    A = coo_matrix((np.ones(len(pairs), bool), (pairs[:, 0], pairs[:, 1])), shape=(len(X), len(X)))
    return connected_components(A, directed=False)[1]


def strokes(res, W=0.4, L=1.0, K=0.26, align=3.0, agree=30.0, minlen=10, log=None, step=0.087, stop_foreign=False, workers=None):
    """W, L, K, step in mm (independent of the cube size): half-width and half-length of the window, the length of dot
    row needed ahead to continue (K/cube = dot count), and the node spacing. minlen in nodes.
    The dots are first grouped into clusters (connected within 3 cubes); each cluster is traced on its own, in
    parallel over the cores, biggest clusters first. Lines never cross between clusters (the join step bridges gaps)."""
    import time
    from concurrent.futures import ProcessPoolExecutor
    line = res['line']; lo, size = res['lo'], float(res['size']); X = (res['idx'][line]+0.5).astype(float)
    wgt = np.maximum(res['turn'][line]-float(res.get('cut', 0)), 0)+1e-3     # weight = turn above the cut, so the stroke centres on the band
    W, L, step = W/size, L/size, step/size; K = max(1, int(round(K/size)))
    t0 = time.time(); lab = clusters(X, 3.0); order = np.argsort(lab, kind='stable'); bounds = np.flatnonzero(np.diff(lab[order]))+1
    groups = [g for g in np.split(order, bounds) if len(g) >= K]
    CAP, B = 20000, 250      # a cluster above CAP dots is cut into blocks of B cubes (10 mm at 0.087) so no single core holds the whole plate outline; join bridges the cuts
    split = []
    for g in groups:
        if len(g) <= CAP: split.append(g); continue
        blk = np.floor(X[g]/B).astype(int); key = (blk[:, 0]*4096+blk[:, 1])*4096+blk[:, 2]; o = np.argsort(key, kind='stable'); b = np.flatnonzero(np.diff(key[o]))+1
        split.extend(g[o][p] for p in np.split(np.arange(len(g)), b) if len(p) >= K)
    groups = sorted(split, key=len, reverse=True)
    if log: log(f'{len(X)} dots in {len(groups)} clusters (largest {len(groups[0]) if groups else 0}), {time.time()-t0:.0f} s')
    tasks = [(X[g], wgt[g], W, L, K, align, agree, minlen, step, stop_foreign) for g in groups]
    workers = workers or os.cpu_count() or 1
    if workers > 1 and len(tasks) > 1:
        with ProcessPoolExecutor(workers) as ex: results = list(ex.map(_trace, tasks, chunksize=1))
    else: results = [_trace(t) for t in tasks]
    curves = [lo+c*size for r in results for c in r[0]]; unclaimed = sum(r[1] for r in results)
    if log: log(f'strokes {len(curves)}, nodes {sum(len(c) for c in curves)}, dots unclaimed {unclaimed} of {len(X)}')
    return curves


def join(curves, G=1.0, TOL=0.15, ANG=25.0, log=None):
    """Join stroke ends (mm). corner: two ends whose headings meet ahead of both within G, rays passing within TOL, are
    extended to the meeting point and merged. gap: two ends facing each other (within ANG deg) closer than G are merged.
    tee: an end whose extrapolation hits another stroke's body within TOL, ahead within G, is extended to the hit point."""
    S = [np.asarray(c, float) for c in curves]; NH = 6
    def head(c, end):
        q = c[-NH:] if end else c[:NH][::-1]; h = q[-1]-q[0]; return q[-1], h/(np.linalg.norm(h)+1e-12)
    def ray_meet(p1, h1, p2, h2):
        w = p1-p2; a = h1 @ h1; b = h1 @ h2; c = h2 @ h2; d = h1 @ w; e = h2 @ w; den = a*c-b*b
        if den < 1e-9: return None
        t = (b*e-c*d)/den; u = (a*e-b*d)/den; return t, u, np.linalg.norm((p1+t*h1)-(p2+u*h2))
    cos_ang = np.cos(np.radians(ANG)); joins = dict(corner=0, gap=0, tee=0); alive = [True]*len(S); changed = True
    while changed:
        changed = False
        ends = [(i, e) for i in range(len(S)) if alive[i] for e in (0, 1)]
        if not ends: break
        E = np.array([head(S[i], e)[0] for i, e in ends]); tree = cKDTree(E)
        used = set()
        for a, (i, ei) in enumerate(ends):
            if not alive[i] or a in used: continue
            pa, ha = head(S[i], ei); best = None
            for b in tree.query_ball_point(pa, G):
                j, ej = ends[b]
                if b == a or b in used or j == i or not alive[j]: continue
                pb, hb = head(S[j], ej)
                if ha @ hb < -cos_ang:
                    d = pb-pa; Ld = np.linalg.norm(d); d = d/(Ld+1e-12)
                    if ha @ d > cos_ang and -hb @ d > cos_ang: cand = ('gap', Ld, j, ej, None)
                    else: continue
                else:
                    r = ray_meet(pa, ha, pb, hb)
                    if r is None: continue
                    t, u, dist = r
                    if 0 < t < G and 0 < u < G and dist < TOL: cand = ('corner', t+u, j, ej, (pa+t*ha+pb+u*hb)/2)
                    else: continue
                if best is None or cand[1] < best[1]: best = cand
            if best is None: continue
            kind, _, j, ej, X = best
            A = S[i] if ei else S[i][::-1]; B = S[j][::-1] if ej else S[j]
            S[i] = np.vstack([A]+([X] if X is not None else [])+[B]); alive[j] = False; joins[kind] += 1; changed = True
            used.update(k for k, (jj, ee) in enumerate(ends) if jj in (i, j))   # both strokes' ends have moved: next pass
    S = [c for c, al in zip(S, alive) if al]
    if S:
        allp = np.vstack(S); owner = np.concatenate([np.full(len(c), k) for k, c in enumerate(S)]); tree = cKDTree(allp)
        for i in range(len(S)):
            for e in (0, 1):
                p, h = head(S[i], e); hit = None
                for t in np.arange(0.05, G, 0.05):
                    q = p+t*h; dd, k = tree.query(q)
                    if owner[k] != i and dd < TOL: hit = allp[k]; break
                if hit is not None: S[i] = np.vstack([S[i], hit[None]]) if e else np.vstack([hit[None], S[i]]); joins['tee'] += 1
    if log: log(f'join: {len(curves)} strokes -> {len(S)}; corners {joins["corner"]}, gaps {joins["gap"]}, tees {joins["tee"]}')
    return S


if __name__ == '__main__':
    """usage: python creases.py mesh out.npz [cells cut W]  (whole model; env SPP, FLAT=0 to blend)"""
    import trimesh, time
    path, outp = sys.argv[1], sys.argv[2]; a = sys.argv[3:]
    cells = int(a[0]) if a else 1200; cut = float(a[1]) if len(a) > 1 else 17.0; W = float(a[2]) if len(a) > 2 else 0.35
    m = trimesh.load(path, process=True); P, F = np.asarray(m.vertices), np.asarray(m.faces)
    t0 = time.time(); res = find_dots(P, F, cells, cut, float(os.environ.get('SPP', 0)), os.environ.get('FLAT', '1') == '1', print); print(f'dots {time.time()-t0:.0f} s')
    t0 = time.time(); cv = join(strokes(res, W=W, log=print), log=print); print(f'strokes {time.time()-t0:.0f} s')
    np.savez(outp, n=np.array([len(c) for c in cv]), pts=np.vstack(cv)); print('saved', outp)
