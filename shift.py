"""Direction-aware test. For each cube, take the neighbour direction of its largest pair angle (across the feature).
Shift on the near side: arrow(cube) - arrow(cube behind). Shift on the far side: arrow(cube ahead) - arrow(cube).
Store how much the shift CHANGES from one side to the other, |far - near|, in degrees. A steady tilt (toe flank)
gives ~0; a wrinkle (tilt one way then back) gives twice the tilt; a step edge gives the step.
Crest = above thr and not smaller than the cubes behind/ahead. Same npz fields as the grid, so trace2.py runs on it.
usage: python shift.py grid.npz out.npz [thr_deg]"""
import sys, numpy as np
g = np.load(sys.argv[1]); outp = sys.argv[2]; thr = float(sys.argv[3]) if len(sys.argv) > 3 else 6.0
idx, N = g['idx'], g['N']; k = len(idx); KB = 1 << 16
key = (idx[:, 0]*KB+idx[:, 1])*KB+idx[:, 2]; lut = dict(zip(key.tolist(), range(k)))
OFFS = np.array([[i, j, l] for i in (-1, 0, 1) for j in (-1, 0, 1) for l in (-1, 0, 1) if (i, j, l) != (0, 0, 0)])
pmax = np.zeros(k); pdir = np.zeros(k, int); nb = np.full((k, 26), -1)
for q, o in enumerate(OFFS):
    kk = ((idx[:, 0]+o[0])*KB+(idx[:, 1]+o[1]))*KB+(idx[:, 2]+o[2])
    j = np.array([lut.get(int(x), -1) for x in kk]); nb[:, q] = j; h = j >= 0
    dot = (N[h]*N[j[h]]).sum(1); ang = np.degrees(np.arccos(np.clip(dot, -1, 1))); ang[dot < 0] = 0
    upd = np.zeros(k, bool); upd[h] = ang > pmax[h]; pmax[upd] = ang[upd[h]]; pdir[upd] = q
opp = np.array([np.nonzero((OFFS == -o).all(1))[0][0] for o in OFFS])
fwd = nb[np.arange(k), pdir]; bwd = nb[np.arange(k), opp[pdir]]; ok = (fwd >= 0) & (bwd >= 0)
near = N-N[np.maximum(bwd, 0)]; far = N[np.maximum(fwd, 0)]-N
change = np.degrees(np.linalg.norm(far-near, axis=1)); change[~ok] = 0          # |far - near| in degrees (small-angle)
cf = np.where(fwd >= 0, change[np.maximum(fwd, 0)], 0); cb = np.where(bwd >= 0, change[np.maximum(bwd, 0)], 0)
crest = (change >= thr) & (change >= cf) & (change > cb)
print(f'{k} cubes; shift change p50 {np.percentile(change, 50):.1f} p90 {np.percentile(change, 90):.1f} p99 {np.percentile(change, 99):.1f} deg; above {thr}: {(change >= thr).sum()}; crest {crest.sum()}')
np.savez(outp, idx=idx, N=N, lo=g['lo'], size=g['size'], turn=change, line=crest, d1=g['d1'])
