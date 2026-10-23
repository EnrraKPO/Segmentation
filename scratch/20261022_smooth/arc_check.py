"""how well does each trend predict the next node, on today's lines? For every node k of every line (ring box), fit the trend
over nodes 0..k and compare its direction with the true step k->k+1. Chord over 6 vs arc with reach R."""
import numpy as np
z = np.load('out/trend_ring.npz'); n = z['off_n']; pts = z['off_pts']; lines = np.split(pts, np.cumsum(n)[:-1])
def arc_dir(P, reach):
    k = len(P)-1; seg = np.linalg.norm(np.diff(P, axis=0), axis=1); t = np.concatenate([[0.0], np.cumsum(seg)]); t = t-t[-1]
    w = np.exp(-(k-np.arange(k+1))/reach); A = np.stack([np.ones(k+1), t, t*t], 1)*w[:, None]
    a0, b, c = np.linalg.lstsq(A, P*w[:, None], rcond=None)[0]; return b/np.linalg.norm(b)
for name, f in [('chord6', lambda P: (P[-1]-P[-min(6, len(P)-1)-1])), ('arc6', lambda P: arc_dir(P, 6)), ('arc20', lambda P: arc_dir(P, 20)), ('arc100', lambda P: arc_dir(P, 100))]:
    errs = []
    for L in lines:
        for k in range(2, len(L)-1):
            d = f(L[:k+1]); d = d/np.linalg.norm(d); s = L[k+1]-L[k]; s /= np.linalg.norm(s)+1e-12
            errs.append(np.degrees(np.arccos(np.clip(d @ s, -1, 1))))
    e = np.array(errs); print(f'{name:7} angle to the true next step: p50 {np.median(e):5.1f}  p90 {np.percentile(e, 90):5.1f}  max {e.max():5.1f}')
