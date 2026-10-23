"""how shy are the refined lines: on the dump, per walk line with refined pieces, number of pieces; for each walk END that
sits on another walk line (within one cube), the distance from the nearest refined piece end to the refined lines of that other walk line"""
import sys, numpy as np; from scipy.spatial import cKDTree
d = sys.argv[1]; cube = 0.087
w = np.load(d+'/lines.npz'); W = [w['pts'][a:a+k] for a, k in zip(np.cumsum(np.r_[0, w['n'][:-1]]), w['n'])]
r = np.load(d+'/refined.npz'); R = [r['pts'][a:a+k] for a, k in zip(np.cumsum(np.r_[0, r['n'][:-1]]), r['n'])]; ri = r['index']
print(f'walk lines {len(W)}, refined pieces {len(R)} from {len(set(ri.tolist()))} walk lines')
pieces = np.bincount(ri); pc = pieces[pieces > 0]; print('pieces per refined walk line: 1:', (pc == 1).sum(), '2:', (pc == 2).sum(), '3+:', (pc >= 3).sum(), 'max', pc.max())
# gaps between the pieces of one walk line: distance between consecutive piece ends along the walk line
gaps = []
for i in np.flatnonzero(pieces >= 2):
    P = [R[k] for k in np.flatnonzero(ri == i)]; L = W[i]; t = cKDTree(L)
    pos = sorted([(t.query(p[0])[1], t.query(p[-1])[1], k) for k, p in enumerate(P)])   # pieces ordered along the walk line
    for (a0, a1, k), (b0, b1, m) in zip(pos[:-1], pos[1:]):
        e = P[k][-1] if a1 >= a0 else P[k][0]; s = P[m][0] if b1 >= b0 else P[m][-1]; gaps.append(np.linalg.norm(s-e))
gaps = np.array(gaps); print(f'gaps between pieces of one line: n {len(gaps)}, p50 {np.percentile(gaps,50):.3f} p90 {np.percentile(gaps,90):.3f} max {gaps.max():.2f} mm')
# ends: walk end on another walk line -> refined end distance to that line's refined pieces (or its walk line if unrefined)
wt = cKDTree(np.vstack(W)); wo = np.repeat(np.arange(len(W)), w['n'])
rt = cKDTree(np.vstack(R)); ro = ri[np.repeat(np.arange(len(R)), r['n'])]
walk_on, ref_on, nref_free = [], [], 0
for i in np.flatnonzero(pieces > 0):
    L = W[i]
    for e in (L[0], L[-1]):
        dd, j = wt.query(e, k=50); other = wo[j][(wo[j] != i) & (dd <= 1.0)]
        if len(other) == 0 or dd[(wo[j] != i)][0] > cube: continue   # this walk end does not sit on another line
        o = other[0]; walk_on.append(dd[wo[j] == o][0])
        # nearest refined piece end of line i to this walk end
        ends = np.array([q for k in np.flatnonzero(ri == i) for q in (R[k][0], R[k][-1])]); re = ends[np.argmin(np.linalg.norm(ends-e, axis=1))]
        tgt = np.vstack([R[k] for k in np.flatnonzero(ri == o)]) if pieces[o] > 0 else W[o]
        ref_on.append(np.linalg.norm(tgt-re, axis=1).min())
walk_on, ref_on = np.array(walk_on), np.array(ref_on)
print(f'walk ends sitting on another line (within one cube): {len(walk_on)}')
print(f'  their refined ends to that line: p50 {np.percentile(ref_on,50):.3f} p90 {np.percentile(ref_on,90):.3f} mm; within one cube {(ref_on<=cube).mean()*100:.0f} %, beyond 0.3 mm {(ref_on>0.3).mean()*100:.0f} %')
# how much of the walk line's length the pieces cover
cov = []
for i in np.flatnonzero(pieces > 0):
    lw = np.linalg.norm(np.diff(W[i], axis=0), axis=1).sum(); lr = sum(np.linalg.norm(np.diff(R[k], axis=0), axis=1).sum() for k in np.flatnonzero(ri == i)); cov.append(lr/max(lw, 1e-9))
cov = np.array(cov); print(f'refined length / walk length: p10 {np.percentile(cov,10):.2f} p50 {np.percentile(cov,50):.2f} p90 {np.percentile(cov,90):.2f}')
