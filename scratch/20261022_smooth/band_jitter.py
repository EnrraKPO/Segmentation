"""is the wiggle in the dots or in our band mean? Ring box, walk with the trend's vote off (pull 0, look 0), so every node is
the bare band mean. Fit a circle to the longest line's nodes; across-residual of the nodes = band-mean jitter; across-residual
of the dots within the band width of the circle = dot scatter. Repeated with LEVP (steep dot-say exponent) off."""
import sys, os, numpy as np
sys.path.insert(0, r'D:\Stuff\PrintSuit\SegmentationTool'); os.chdir(r'D:\Stuff\PrintSuit\SegmentationTool'); sys.path.insert(0, 'scratch/20261022_smooth'); import smooth_test as T, trace4
z = np.load(T.DOTS, allow_pickle=True); lo = z['lo']; size = float(z['size']); cut = float(z['cut']); idx, N, turn = z['idx'], z['N'], z['turn']; C = lo+(idx+0.5)*size
focus = np.array(T.SPOTS['ring'][0]); sel = (np.abs(C-focus) < T.BOX).all(1); X = (idx[sel]+0.5).astype(float); wgt = np.maximum(turn[sel]-cut, 0)+1e-3
def circle_fit(P):
    c0 = P.mean(0); u, s, vt = np.linalg.svd(P-c0); n = vt[2]; Q = (P-c0) @ vt[:2].T          # plane, then 2-D circle (algebraic fit)
    A = np.c_[2*Q, np.ones(len(Q))]; b = (Q**2).sum(1); x = np.linalg.lstsq(A, b, rcond=None)[0]; cc = x[:2]; r = np.sqrt(x[2]+cc @ cc)
    centre = c0+vt[:2].T @ cc; return centre, n, r
for name, extra in [('today', {}), ('dots only', dict(pull=0, look=0)), ('dots only, LEVP off', dict(pull=0, look=0, levp=0)), ('dots only, LEVP off, LEAD off', dict(pull=0, look=0, levp=0, lead=0))]:
    lines, _ = trace4.walk(X, wgt, N[sel], size, log=lambda *_: None, **{**T.walk_params(), 'curve': 0, 'trend': 6, **extra}); L = max(lines, key=len)
    centre, n, r = circle_fit(L); rad = lambda P: np.linalg.norm((P-centre)-((P-centre) @ n)[:, None]*n, axis=1)-r; up = lambda P: (P-centre) @ n
    res = np.hypot(rad(L), up(L)); near = np.hypot(rad(X), up(X)) < 3; dres = np.hypot(rad(X[near]), up(X[near]))
    v = np.diff(L, axis=0); v /= np.linalg.norm(v, axis=1, keepdims=True); tr = np.degrees(np.arccos(np.clip((v[1:]*v[:-1]).sum(1), -1, 1)))
    print(f'{name:28} longest line {len(L):4} nodes, radius {r*size:.2f} mm; node off the circle (cubes): p50 {np.median(res):.2f} p90 {np.percentile(res, 90):.2f}; dots within 3 cubes: {near.sum()}, off the circle p50 {np.median(dres):.2f} p90 {np.percentile(dres, 90):.2f}; turn p50 {np.median(tr):.1f} p90 {np.percentile(tr, 90):.1f}')
