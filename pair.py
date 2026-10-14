"""Pair angles: for every surface cube, the angle between its arrow and each of its 26 neighbours' arrows, NO blending.
Stores pmax (largest neighbour angle, deg), pdir (offset of that neighbour) and a crest flag: pmax above thr and not
smaller than pmax of the cubes on either side along pdir.
usage: python pair.py grid.npz out.npz [thr_deg]"""
import sys, numpy as np
g = np.load(sys.argv[1]); outp = sys.argv[2]; thr = float(sys.argv[3]) if len(sys.argv) > 3 else 15.0
idx, N = g['idx'], g['N']; k = len(idx); KB = 1 << 16
key = (idx[:, 0]*KB+idx[:, 1])*KB+idx[:, 2]; lut = dict(zip(key.tolist(), range(k)))
OFFS = np.array([[i, j, l] for i in (-1, 0, 1) for j in (-1, 0, 1) for l in (-1, 0, 1) if (i, j, l) != (0, 0, 0)])
pmax = np.zeros(k); pdir = np.zeros(k, int); nb = np.full((k, 26), -1)
for q, o in enumerate(OFFS):
    kk = ((idx[:, 0]+o[0])*KB+(idx[:, 1]+o[1]))*KB+(idx[:, 2]+o[2])
    j = np.array([lut.get(int(x), -1) for x in kk]); nb[:, q] = j; h = j >= 0
    dot = (N[h]*N[j[h]]).sum(1); ang = np.degrees(np.arccos(np.clip(dot, -1, 1))); ang[dot < 0] = 0      # a neighbour facing away is the far side of a thin wall, not a crease
    upd = np.zeros(k, bool); upd[h] = ang > pmax[h]; pmax[upd] = ang[upd[h]]; pdir[upd] = q
# crest along the direction of the largest pair: the cube must be at least as sharp as the cube on either side
opp = np.array([np.nonzero((OFFS == -o).all(1))[0][0] for o in OFFS])
fwd = nb[np.arange(k), pdir]; bwd = nb[np.arange(k), opp[pdir]]
pf = np.where(fwd >= 0, pmax[np.maximum(fwd, 0)], 0); pb = np.where(bwd >= 0, pmax[np.maximum(bwd, 0)], 0)
crest = (pmax >= thr) & (pmax >= pf) & (pmax > pb)
print(f'{k} cubes; pair angle p50 {np.percentile(pmax, 50):.1f} p90 {np.percentile(pmax, 90):.1f} p99 {np.percentile(pmax, 99):.1f} deg; above {thr}: {(pmax >= thr).sum()}; crest cubes {crest.sum()}')
np.savez(outp, idx=idx, N=N, lo=g['lo'], size=g['size'], turn=pmax, line=crest, d1=g['d1'])
