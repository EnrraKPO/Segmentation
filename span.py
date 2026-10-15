"""Turn over a span: per cube, the angle between the facings s cubes behind and s cubes ahead along the across
direction (the neighbour direction with the largest facing angle). Saves a grid npz with turn=that angle and
line=crest (>= cut and a local max along the across direction, ties kept).
usage: python span.py grid.npz out.npz s [cut]"""
import sys, numpy as np
g = np.load(sys.argv[1]); outp = sys.argv[2]; s = int(sys.argv[3]); cut = float(sys.argv[4]) if len(sys.argv) > 4 else 17.0
idx, N = g['idx'], g['N']; k = len(idx); KB = 1 << 16
key = (idx[:, 0]*KB+idx[:, 1])*KB+idx[:, 2]; srt = np.argsort(key, kind='stable'); skey = key[srt]
def look(o):
    kk = ((idx[:, 0]+o[0])*KB+(idx[:, 1]+o[1]))*KB+(idx[:, 2]+o[2])
    pos = np.minimum(np.searchsorted(skey, kk), k-1); return np.where(skey[pos] == kk, srt[pos], -1)
OFFS = np.array([[i, j, l] for i in (-1, 0, 1) for j in (-1, 0, 1) for l in (-1, 0, 1) if (i, j, l) != (0, 0, 0)])
def ang(a, b): return np.degrees(np.arccos(np.clip((a*b).sum(1), -1, 1)))
pmax = np.zeros(k); pdir = np.zeros(k, int)
for q, o in enumerate(OFFS):
    j = look(o); h = j >= 0; a = np.zeros(k); a[h] = ang(N[h], N[j[h]]); a[h & ((N*N[np.maximum(j, 0)]).sum(1) < 0)] = 0
    upd = a > pmax; pmax[upd] = a[upd]; pdir[upd] = q
# s cubes each way along the across direction, on the surface: the across direction is the best-neighbour offset
# projected onto the cube's tangent plane; the target point p +- s*dir is matched to the nearest cube centre (within 0.9 cube)
from scipy.spatial import cKDTree
Xc = idx.astype(float); o = OFFS[pdir].astype(float); d = o-(o*N).sum(1, keepdims=True)*N; nd = np.linalg.norm(d, axis=1); d /= np.maximum(nd, 1e-9)[:, None]
tree = cKDTree(Xc)
def at(sign):
    dist, j = tree.query(Xc+sign*s*d, distance_upper_bound=0.9); okk = np.isfinite(dist) & (nd > 1e-6); return np.where(okk, j, 0), okk
fw, af = at(1); bw, ab = at(-1); ok = af & ab
turn = np.zeros(k); turn[ok] = ang(N[fw[ok]], N[bw[ok]])
# crest along the across direction: compare with the immediate neighbours ahead/behind
o = OFFS[pdir]
def nb1(sign):
    kk = ((idx[:, 0]+sign*o[:, 0])*KB+(idx[:, 1]+sign*o[:, 1]))*KB+(idx[:, 2]+sign*o[:, 2])
    pos = np.minimum(np.searchsorted(skey, kk), k-1); return np.where(skey[pos] == kk, srt[pos], -1)
nf, nb = nb1(1), nb1(-1)
tf = np.where(nf >= 0, turn[np.maximum(nf, 0)], 0); tb = np.where(nb >= 0, turn[np.maximum(nb, 0)], 0)
crest = (turn >= cut) & (turn >= tf) & (turn >= tb)
print(f'{k} cubes, span {s}: turn p50 {np.percentile(turn, 50):.1f} p90 {np.percentile(turn, 90):.1f} p99 {np.percentile(turn, 99):.1f}; above {cut}: {(turn >= cut).sum()}; crest {crest.sum()}')
np.savez(outp, idx=idx, N=N, lo=g['lo'], size=g['size'], turn=turn, line=crest, d1=g['d1'])
