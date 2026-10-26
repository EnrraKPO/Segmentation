"""one-off: regions.py gains (a) bisection of stretches that found no strip, (b) a big-triangle pass in project(), (c) absorb()
for regions below an area handle, (d) run(min_area=...)"""
import os
root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); path = os.path.join(root, 'regions.py')
src = open(path, encoding='utf-8').read()

# (b) project: big triangles get their own candidate pass before the exact fallback
old = src[src.index('def project('):src.index('# ------------------------------------------------------------------ 1. strokes across the surface')]
new = '''def project(P, F, pts, k=16, log=None):
    """closest surface point and its triangle for every point. Candidates: the k triangles with the nearest centres, plus
    the nearest of the large triangles (radius above the 95th percentile) whose reach allows them; a point is trusted
    when no other triangle can come closer (the next centre is farther than the best distance plus the radius bound);
    the rest go to trimesh's exact query."""
    from scipy.spatial import cKDTree
    T = P[F]; cen = T.mean(1); R = np.linalg.norm(T-cen[:, None, :], axis=2).max(1)
    n = len(pts); best = np.full(n, np.inf); bestq = np.zeros((n, 3)); bestf = np.zeros(n, np.int64)
    def take(f):
        q = _closest_on_triangles(pts, T[f, 0], T[f, 1], T[f, 2]); d = np.linalg.norm(q-pts, axis=1)
        m = d < best; best[m] = d[m]; bestq[m] = q[m]; bestf[m] = f[m]
    Rbig = float(np.percentile(R, 95)); big = np.where(R > Rbig)[0]
    tree = cKDTree(cen); dk, ik = tree.query(pts, k=k)
    for j in range(k): take(ik[:, j])
    sure = dk[:, -1]-Rbig > best                                                      # no ordinary triangle beyond the k can come closer
    if len(big):
        kb = min(8, len(big)); dkb, ikb = cKDTree(cen[big]).query(pts, k=kb); dkb = dkb.reshape(n, -1); ikb = ikb.reshape(n, -1)
        for j in range(kb): take(big[ikb[:, j]])
        sure &= dkb[:, -1]-float(R[big].max()) > best
    unsure = np.where(~sure)[0]
    if len(unsure):
        if log: log(f'  {len(unsure)} of {n} nodes need the exact query')
        import trimesh
        mesh = trimesh.Trimesh(P, F, process=False); q, d, f = trimesh.proximity.closest_point(mesh, pts[unsure])
        bestq[unsure] = q; bestf[unsure] = f; best[unsure] = d
    return bestq, bestf, best


'''
src = src.replace(old, new)

# (a) embed: bisection rounds
old = src[src.index('def embed('):src.index('# ------------------------------------------------------------------ 2. partition')]
new = '''def embed(P, F, lines, log=None, depth=8, rounds=8):
    """lines (list of (k,3) arrays) onto the surface: see _walk. A stretch that finds no strip (two nodes projected far
    apart on the surface, or across a sharp fold) is bisected: the midpoint of its two surface points is dropped onto the
    surface as a new node, and the walk runs again, up to `rounds` times. Returns dict(strokes, edge points, nodes, ...)."""
    t0 = time.time(); P = np.ascontiguousarray(P, np.float64); F = np.ascontiguousarray(F, np.int64)
    lines = [np.asarray(L, float).reshape(-1, 3) for L in lines if len(L) > 1]
    nodes = np.vstack(lines); start = np.r_[0, np.cumsum([len(L) for L in lines])].astype(np.int64)
    E, fedge, nbr = edges_of(F); v0, e1, e2, C2 = frames(P, F)
    Q, nf, dist = project(P, F, nodes, log=log)
    if log: log(f'  {len(nodes)} nodes onto the surface in {time.time()-t0:.1f} s (distance median {np.median(dist):.3f}, max {dist.max():.3f} mm)')
    def frame2(Q, nf): d = Q-v0[nf]; return np.stack([(d*e1[nf]).sum(1), (d*e2[nf]).sum(1)], 1)
    node2 = frame2(Q, nf); added = 0
    for rnd in range(rounds+1):
        cap = 4*len(Q)+1000
        while True:
            sf, sp, sq, s0, s1, sl, epe, ept, ep3, fid, frs, over = _walk(P, F, Q, nf, node2, start, C2, nbr, fedge, E, depth, cap)
            if not over: break
            cap *= 2
        if len(fid) == 0 or rnd == rounds: break
        fid = np.unique(fid); mid = 0.5*(Q[fid]+Q[fid+1]); mq, mf, _ = project(P, F, mid)
        Q = np.insert(Q, fid+1, mq, axis=0); nf = np.insert(nf, fid+1, mf); node2 = np.insert(node2, fid+1, frame2(mq, mf), axis=0)
        start = start+np.searchsorted(fid, start, side='left'); added += len(fid)                  # a line's offset grows by the insertions before it
    if log: log(f'  {len(sf)} strokes in {len(np.unique(sf))} crossed triangles, {len(epe)} edge crossings, {added} nodes added by bisection, {len(fid)} stretches not drawn through, {time.time()-t0:.1f} s')
    return dict(P=P, F=F, E=E, fedge=fedge, nbr=nbr, v0=v0, e1=e1, e2=e2, C2=C2, nodes=Q, nface=nf, node2=node2, start=start,
                sf=sf, sp=sp, sq=sq, s0=s0, s1=s1, sl=sl, epe=epe, ept=ept, ep3=ep3, fails=len(fid), fid=fid, added=added)


'''
src = src.replace(old, new)

# (c) absorb
old = "PALETTE = ["
new = '''def absorb(con, min_area, log=None):
    """regions below min_area (mm2) are merged into their largest neighbour, smallest first; regions without a neighbour
    stay. Rewrites con (region, nreg, pairs, areas) in place; returns the number absorbed."""
    t0 = time.time(); nreg = con['nreg']; areas = con['areas'].copy(); root = np.arange(nreg)
    def find(r):
        while root[r] != r: root[r] = root[root[r]]; r = root[r]
        return r
    nb = [set() for _ in range(nreg)]
    for a, b in con['pairs']: nb[a].add(b); nb[b].add(a)
    done = 0
    for r in sorted(range(nreg), key=lambda r: areas[r]):
        if areas[r] >= min_area: break
        if find(r) != r: continue
        cand = {find(q) for q in nb[r]}; cand.discard(r)
        if not cand: continue
        t = max(cand, key=lambda q: areas[q]); root[r] = t; areas[t] += areas[r]; nb[t] |= nb[r]; done += 1
    roots = np.array([find(r) for r in range(nreg)]); keep = np.unique(roots); remap = np.full(nreg, -1, np.int64); remap[keep] = np.arange(len(keep))
    reg = con['region']; m = reg >= 0; reg[m] = remap[roots[reg[m]]]
    con['nreg'] = len(keep); con['areas'] = np.bincount(remap[roots], weights=areas*(roots == np.arange(nreg)), minlength=len(keep))
    con['pairs'] = sorted({(min(remap[roots[a]], remap[roots[b]]), max(remap[roots[a]], remap[roots[b]])) for a, b in con['pairs'] if roots[a] != roots[b]})
    if log: log(f'  {done} regions below {min_area} mm2 absorbed into a neighbour: {con["nreg"]} regions, {time.time()-t0:.1f} s')
    return done


PALETTE = ['''
src = src.replace(old, new, 1)

# (d) run
old = '''def run(P, F, lines, log=None, k=9, workers=None):
    t0 = time.time()
    if log: log(f'embedding {len(lines)} lines on {len(F)} triangles')
    emb = embed(P, F, lines, log=log)
    part = partition(emb, log=log, workers=workers)
    con = connect(emb, part, log=log)
'''
new = '''def run(P, F, lines, log=None, k=9, workers=None, min_area=0.0, depth=8):
    t0 = time.time()
    if log: log(f'embedding {len(lines)} lines on {len(F)} triangles')
    emb = embed(P, F, lines, log=log, depth=depth)
    part = partition(emb, log=log, workers=workers)
    con = connect(emb, part, log=log)
    if min_area > 0: absorb(con, min_area, log=log)
'''
assert old in src; src = src.replace(old, new)
open(path, 'w', encoding='utf-8').write(src); print('patched')
