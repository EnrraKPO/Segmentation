"""Thick-line walk: a smooth stroke of half-width W that covers the dots (line cubes). At each step the stroke looks at
the dots inside a window L long (ahead and behind) and W wide, takes their weighted centre as its position and the
direction they are arranged along (principal axis) as its heading. No per-cube direction is used. It stops when the
window ahead holds fewer than K dots, when the heading turns more than `agree` degrees in one step, or when it meets
dots already claimed by another stroke. Seeds: strongest unclaimed dot that has at least K dots in its window.
usage: python trace3.py grid.npz out.npz      env: W=3 L=12 K=3 ALIGN=3 STEP=1 AGREE=30 MINLEN=10 (cubes / degrees / nodes; ALIGN = how much longer than wide the dot cloud in the window must be)"""
import sys, os, numpy as np
from scipy.spatial import cKDTree
g = np.load(sys.argv[1]); outp = sys.argv[2]
W = float(os.environ.get('W', 3)); L = float(os.environ.get('L', 12)); K = int(os.environ.get('K', 3)); ALIGN = float(os.environ.get('ALIGN', 3)); STEP = float(os.environ.get('STEP', 1))
AGREE = np.cos(np.radians(float(os.environ.get('AGREE', 30)))); MINLEN = int(os.environ.get('MINLEN', 10))
line = g['line']; X = (g['idx'][line]+0.5).astype(float); wgt = g['turn'][line]; N = g['N'][line]; lo, size = g['lo'], g['size']
tree = cKDTree(X); owner = np.full(len(X), -1); R = np.hypot(L, W)
print(f'{len(X)} dots')


def window(p, h):
    """dots inside the stroke window at p heading h: indices, along- and across-offsets"""
    nb = np.array(tree.query_ball_point(p, R), int)
    if len(nb) == 0: return nb, None, None
    off = X[nb]-p; al = off @ h; side = off-al[:, None]*h; sd = np.linalg.norm(side, axis=1)
    keep = (np.abs(al) <= L) & (sd <= W); return nb[keep], al[keep], side[keep]


def fit(p, h, tid):
    """position and heading of the stroke at p: weighted centre and principal axis of the window's dots"""
    nb, al, side = window(p, h)
    if len(nb) < K: return None
    w = wgt[nb]; c = (X[nb]*w[:, None]).sum(0)/w.sum()
    d = X[nb]-c; cov = (d*w[:, None]).T @ d; ev, vec = np.linalg.eigh(cov); ax = vec[:, -1]
    if ev[-1] < ALIGN*ALIGN*max(ev[-2], 1e-9): return None                       # dots not arranged along a line
    if ax @ h < 0: ax = -ax
    q = p+(c-p)-((c-p) @ h)*h                                                # slide sideways only
    ahead = (al > 0).sum(); foreign = [k for k in nb if owner[k] not in (-1, tid)]
    return q, ax, ahead, foreign


def march(p, h, tid):
    outp = []
    while True:
        q = p+STEP*h; r = fit(q, h, tid)
        if r is None: break
        q, hq, ahead, foreign = r
        if ahead < K: break                                                   # stroke ahead is empty: line ended
        if hq @ h < AGREE: break                                              # sharp corner
        if len(foreign) > K: break                                            # ran into another stroke
        back = int(2*L)+5
        if len(outp) > back and np.linalg.norm(np.array(outp[:-back])-q, axis=1).min() < W: break   # met its own tail (closed loop)
        for k in tree.query_ball_point(q, W): owner[k] = tid
        outp.append(q); p, h = q, hq
        if len(outp) > 20000: break
    return outp


strokes = []; order = np.argsort(-wgt)
for s in order:
    if owner[s] != -1: continue
    p0 = X[s]
    # initial heading: principal axis of all dots around the seed
    nb = np.array(tree.query_ball_point(p0, L), int)
    if len(nb) < K: continue
    d = X[nb]-X[nb].mean(0); ev, vec = np.linalg.eigh(d.T @ d); h0 = vec[:, -1]
    r = fit(p0, h0, len(strokes))
    if r is None: continue
    p0, h0 = r[0], r[1]; tid = len(strokes)
    for k in tree.query_ball_point(p0, W): owner[k] = tid
    fw = march(p0, h0, tid); bw = march(p0, -h0, tid); ch = bw[::-1]+[p0]+fw
    if len(ch) < MINLEN: owner[owner == tid] = -1; strokes.append(None); continue
    strokes.append(np.array(ch))
strokes = [c for c in strokes if c is not None]; Ls = [len(c) for c in strokes]
print(f'strokes {len(strokes)}, nodes {sum(Ls)}, longest {max(Ls) if Ls else 0}, dots unclaimed {(owner == -1).sum()} of {len(X)}')
np.savez(outp, n=np.array(Ls), pts=np.vstack([lo+c*size for c in strokes]))
print('saved', outp)
