"""Stage 2b: merge the redundant walk lines (rule agreed 2026-10-23, late).

Two lines are one line wherever averaging them would barely move either: where a line's nodes lie within MARGIN of the
other line (the walk node's own scatter about its crease, measured 0.07 mm median), the stretch is replaced by one shared
stretch, the average of the two, and the rest of each line continues from it. The length of the stretch is not a
condition: a long stretch is one crease walked twice, now one; a stretch of a node or two is a crossing or a touch, and
the shared node is the junction record. The longer line keeps the stretch and stays whole; the shorter line's parts
outside the stretch become lines of their own, each ending on the shared node at the stretch's end. Repeated in passes
until no pair qualifies (a merge can bring a third line within the margin).

run(curves, margin=0.07, labels=None) -> dict(lines, labels, shared, report)   lengths in mm. labels: one per line (e.g. the
walk line a refined piece came from), carried through: a kept line keeps its label, a tail keeps its line's.
The same stage serves the refined pieces (3b): several fine pieces of one tube retrace one crease."""
import numpy as np
from scipy.spatial import cKDTree

DEFAULTS = dict(margin=0.07)


def _seg_dist(p, L, k):
    """distance from point p to line L near its node k: the nearest of the two segments at node k"""
    best = np.linalg.norm(L[k]-p); hit = L[k]
    for a, b in ((k-1, k), (k, k+1)):
        if a < 0 or b >= len(L): continue
        d = L[b]-L[a]; n2 = d @ d
        if n2 < 1e-18: continue
        t = np.clip(((p-L[a]) @ d)/n2, 0, 1); q = L[a]+t*d; dd = np.linalg.norm(q-p)
        if dd < best: best, hit = dd, q
    return best, hit


def _runs(mask):
    """maximal runs of True as (start, end) inclusive"""
    out = []; k = 0; n = len(mask)
    while k < n:
        if mask[k]:
            s = k
            while k+1 < n and mask[k+1]: k += 1
            out.append((s, k))
        k += 1
    return out


def one_pass(lines, margin, labels=None):
    """one pass over all pairs: returns the new line list, the shared nodes made, and the number of merges"""
    lines = [np.asarray(l, float) for l in lines]; n = len(lines); labels = list(range(n)) if labels is None else list(labels)
    allp = np.vstack(lines); owner = np.repeat(np.arange(n), [len(l) for l in lines]); local = np.concatenate([np.arange(len(l)) for l in lines]); tree = cKDTree(allp)
    order = sorted(range(n), key=lambda i: -len(lines[i]))                      # the longer line keeps the stretch
    rank = np.empty(n, int); rank[order] = np.arange(n)
    touched = np.zeros(n, bool); out = {i: [lines[i]] for i in range(n)}; shared = []; merges = 0
    for j in order:                                                              # j = the line whose stretch is given up
        if touched[j] or len(lines[j]) < 2: continue
        Lj = lines[j]; nb = tree.query_ball_point(Lj, margin*2)                  # candidates: the other lines near each node of j
        cand = {}
        for k, lst in enumerate(nb):
            for m in lst:
                i = int(owner[m])
                if i == j or touched[i] or rank[i] > rank[j]: continue            # only onto a longer (or equal, earlier) line
                cand.setdefault(i, set()).add(k)
        if not cand: continue
        i = max(cand, key=lambda i: len(cand[i]))                                # the line j runs along the most
        Li = lines[i]; ti = cKDTree(Li); near = np.zeros(len(Lj), bool); proj = np.zeros((len(Lj), 3)); idx = np.zeros(len(Lj), int)
        for k in range(len(Lj)):
            d0, m = ti.query(Lj[k]); dd, q = _seg_dist(Lj[k], Li, int(m)); near[k] = dd <= margin; proj[k] = q; idx[k] = m
        done = np.array([np.linalg.norm(proj[k]-Lj[k]) < 0.1*margin for k in range(len(Lj))])   # nodes already on line i within a tenth of the margin (a tail's end on its shared node, moved a hair by a later merge): merged before
        runs = [r for r in _runs(near) if not done[r[0]:r[1]+1].all()]
        if not runs: continue
        # apply the longest run (the others wait for the next pass, since i changes here)
        c, d = max(runs, key=lambda r: r[1]-r[0]); a, b = int(idx[c:d+1].min()), int(idx[c:d+1].max())
        Li2 = Li.copy(); tj = cKDTree(Lj[c:d+1])
        for k in range(a, b+1):                                                  # the stretch on i: each node averaged with its foot on j
            _, m = tj.query(Li[k]); dd, q = _seg_dist(Li[k], Lj[c:d+1], int(m)); Li2[k] = 0.5*(Li[k]+q)
        pieces = []
        if c >= 1: head = np.vstack([Lj[:c], Li2[a if idx[c] <= idx[d] else b][None]]); pieces.append(head)         # j's tail before the stretch ends on the shared node
        if d <= len(Lj)-2: tail = np.vstack([Li2[b if idx[c] <= idx[d] else a][None], Lj[d+1:]]); pieces.append(tail)
        out[i] = [Li2]; out[j] = [p for p in pieces if len(p) >= 2]; touched[i] = touched[j] = True; merges += 1
        shared.append((i, j, Li2[a].copy(), Li2[b].copy(), b-a+1, d-c+1))
    new = [p for i in range(n) for p in out[i]]; newlab = [labels[i] for i in range(n) for _ in out[i]]
    return new, newlab, shared, merges


def run(curves, margin=0.07, log=None, max_passes=50, labels=None):
    lines = [np.asarray(c, float) for c in curves]; labels = list(range(len(lines))) if labels is None else list(labels); n0 = len(lines); nodes0 = sum(len(l) for l in lines); shared_all = []; passes = 0; total = 0
    while passes < max_passes:
        lines, labels, shared, merges = one_pass(lines, margin, labels); passes += 1; total += merges; shared_all += shared
        if log: log(f'  merge pass {passes}: {merges} stretches merged, {len(lines)} lines')
        if merges == 0: break
    long = sum(1 for s in shared_all if s[4] >= 3); short = len(shared_all)-long
    rep = dict(merges=total, long_stretches=long, touches=short, lines_before=n0, lines_after=len(lines), nodes_before=nodes0, nodes_after=sum(len(l) for l in lines), passes=passes)
    if log: log(f'merge: {n0} lines / {nodes0} nodes -> {len(lines)} lines / {rep["nodes_after"]} nodes; {total} merges ({long} stretches of 3+ nodes, {short} touches) in {passes} passes')
    return dict(lines=lines, labels=labels, shared=shared_all, report=rep)
