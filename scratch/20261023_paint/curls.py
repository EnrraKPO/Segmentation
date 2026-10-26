"""curls at line ends (turn over the last millimetre) and wobble along the body (lateral deviation from a 1 mm running mean)"""
import sys, os, pickle; root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); sys.path.insert(0, root)
import numpy as np, pyvista as pv
res = pickle.load(open(os.path.join(root, 'out', 'paint_dump.pkl'), 'rb')); emb = res['emb']; Q = emb['nodes']; start = emb['start']; N = np.cross(emb['e1'], emb['e2']); P = emb['P']
def turn_last(L, mm=1.0):
    """total heading change over the last mm of a chain (degrees)"""
    s = np.r_[0, np.cumsum(np.linalg.norm(np.diff(L, axis=0), axis=1))]
    k = np.searchsorted(s, s[-1]-mm); seg = L[k:]
    if len(seg) < 3: return 0.0
    d = np.diff(seg, axis=0); d /= np.linalg.norm(d, axis=1)[:, None]+1e-12
    return float(np.degrees(np.arccos(np.clip((d[1:]*d[:-1]).sum(1), -1, 1))).sum())
def wobble(L, mm=1.0):
    """lateral deviation of the nodes from a running mean over mm (median, mm)"""
    s = np.r_[0, np.cumsum(np.linalg.norm(np.diff(L, axis=0), axis=1))]
    if s[-1] < 2*mm: return np.nan
    dev = []
    for i in range(len(L)):
        m = (s >= s[i]-mm/2) & (s <= s[i]+mm/2); dev.append(np.linalg.norm(L[i]-L[m].mean(0)))
    return float(np.median(dev))
turns = []; wob = []; ends = []
for li in range(len(start)-1):
    L = Q[start[li]:start[li+1]]
    if len(L) < 12: continue
    for e, LL in ((0, L[::-1]), (1, L)):
        t = turn_last(LL); turns.append(t); ends.append((li, e, t))
    wob.append(wobble(L))
turns = np.array(turns); wob = np.array([w for w in wob if not np.isnan(w)])
print(f'{len(turns)} ends of lines with 12+ nodes: turn over the last 1 mm  >90 deg: {np.mean(turns>90):.0%}, >180 deg: {np.mean(turns>180):.0%}, median {np.median(turns):.0f} deg')
print(f'wobble (median lateral deviation from a 1 mm running mean): median over lines {np.median(wob):.3f} mm, p90 {np.percentile(wob, 90):.3f} mm')
# gallery of the strongest curls
ends.sort(key=lambda x: -x[2]); pick = ends[:8]
pl = pv.Plotter(off_screen=True, shape=(2, 4), window_size=(3200, 1600), border=True)
mesh = pv.PolyData(P, np.hstack([np.full((len(emb['F']), 1), 3), emb['F']]).ravel())
for k, (li, e, t) in enumerate(pick):
    L = Q[start[li]:start[li+1]]; p = L[0] if e == 0 else L[-1]; R = 1.5; pl.subplot(k//4, k % 4)
    box = [p[0]-R, p[0]+R, p[1]-R, p[1]+R, p[2]-R, p[2]+R]; pl.add_mesh(mesh.clip_box(box, invert=False), color='lightgrey', smooth_shading=False)
    LL = L+0.01*N[emb['nface'][start[li]:start[li+1]]]; pl.add_mesh(pv.PolyData(LL, lines=np.r_[len(LL), np.arange(len(LL))]), color='blue', line_width=4)
    pl.add_mesh(pv.PolyData(LL), color='blue', point_size=6, render_points_as_spheres=True)
    pl.add_mesh(pv.PolyData(p.reshape(1, 3)), color='red', point_size=14, render_points_as_spheres=True)
    n = N[emb['nface'][start[li] if e == 0 else start[li+1]-1]]; pl.camera_position = [p+3*R*n, p, [0, 0, 1] if abs(n[2]) < 0.9 else [1, 0, 0]]
    pl.add_text(f'line {li} end {e}: {t:.0f} deg over the last 1 mm', font_size=11)
pl.screenshot(os.path.join(root, 'out', 'curls_gallery.png')); print('wrote out/curls_gallery.png')
