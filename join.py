"""Join stroke ends. Two rules, in mm:
  corner: ends A and B whose extrapolated headings meet ahead of both within G, rays passing within TOL -> both extended
          to the meeting point and merged into one stroke.
  gap:    ends A and B facing each other (A's heading points at B and B's at A within ANG deg), distance < G -> merged.
  tee:    end A whose extrapolation hits another stroke's body within TOL, ahead within G -> A extended to the hit point.
usage: python join.py in.npz out.npz [G TOL ANG]   (defaults 1.0 0.15 25)"""
import sys, numpy as np
from scipy.spatial import cKDTree
inp, outp = sys.argv[1:3]; G = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0; TOL = float(sys.argv[4]) if len(sys.argv) > 4 else 0.15; ANG = float(sys.argv[5]) if len(sys.argv) > 5 else 25.0
z = np.load(inp); pts, n = z['pts'], z['n']; S = [pts[a:a+k].copy() for a, k in zip(np.cumsum(np.r_[0, n[:-1]]), n)]
NH = 6
def head(c, end):   # unit heading pointing OUT of the stroke at that end, and the end point
    q = c[-NH:] if end else c[:NH][::-1]
    h = q[-1]-q[0]; return q[-1], h/(np.linalg.norm(h)+1e-12)
def ray_meet(p1, h1, p2, h2):
    """closest points of rays p1+t h1, p2+u h2: returns t, u, distance between closest points"""
    w = p1-p2; a = h1 @ h1; b = h1 @ h2; c = h2 @ h2; d = h1 @ w; e = h2 @ w; den = a*c-b*b
    if den < 1e-9: return None
    t = (b*e-c*d)/den; u = (a*e-b*d)/den; return t, u, np.linalg.norm((p1+t*h1)-(p2+u*h2))
cos_ang = np.cos(np.radians(ANG)); joins = dict(corner=0, gap=0, tee=0)
alive = [True]*len(S)
changed = True
while changed:
    changed = False
    ends = [(i, e) for i in range(len(S)) if alive[i] for e in (0, 1)]
    E = np.array([head(S[i], e)[0] for i, e in ends]); tree = cKDTree(E)
    used = set()
    for a, (i, ei) in enumerate(ends):
        if a in used or not alive[i]: continue
        pa, ha = head(S[i], ei); best = None
        for b in tree.query_ball_point(pa, G):
            j, ej = ends[b]
            if b == a or b in used or j == i or not alive[j]: continue
            pb, hb = head(S[j], ej)
            if ha @ hb < -cos_ang:   # facing each other: gap
                d = pb-pa; L = np.linalg.norm(d); d = d/(L+1e-12)
                if ha @ d > cos_ang and -hb @ d > cos_ang: cand = ('gap', L, j, ej, None)
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
        A = S[i] if ei else S[i][::-1]; B = S[j][::-1] if ej else S[j]      # A ends at its joining end, B starts at its joining end
        mid = [X] if X is not None else []
        S[i] = np.vstack([A]+mid+[B]); alive[j] = False; joins[kind] += 1; changed = True
        used.add(a); used.add([k for k, (jj, ee) in enumerate(ends) if jj == j and ee == ej][0]); break   # rebuild ends after each merge
S = [c for c, al in zip(S, alive) if al]
# tee: extend ends to the body of another stroke
allp = np.vstack(S); owner = np.concatenate([np.full(len(c), k) for k, c in enumerate(S)]); tree = cKDTree(allp)
for i in range(len(S)):
    for e in (0, 1):
        p, h = head(S[i], e); hit = None
        for t in np.arange(0.05, G, 0.05):
            q = p+t*h; dd, k = tree.query(q)
            if owner[k] != i and dd < TOL: hit = allp[k]; break
        if hit is not None:
            S[i] = np.vstack([S[i], hit[None]]) if e else np.vstack([hit[None], S[i]]); joins['tee'] += 1
print(f'{len(n)} strokes -> {len(S)}; joins {joins}')
np.savez(outp, n=np.array([len(c) for c in S]), pts=np.vstack(S))
