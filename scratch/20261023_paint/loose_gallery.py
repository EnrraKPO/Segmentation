"""close-ups of loose ends of the walk set on the orc: the end (red), the lines around it (black), its own line (blue)"""
import sys, os, pickle; root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); sys.path.insert(0, root)
import numpy as np, pyvista as pv, regions
from scipy.spatial import cKDTree
res = pickle.load(open(os.path.join(root, 'out', 'paint_dump.pkl'), 'rb')); emb = res['emb']; P, F = emb['P'], emb['F']
Q = emb['nodes']; start = emb['start']; line = np.repeat(np.arange(len(start)-1), np.diff(start)); t = cKDTree(Q)
loose, _ = regions.loose_ends(emb); rng = np.random.default_rng(3); pick = rng.choice(len(loose), 12, replace=False)
mesh = pv.PolyData(P, np.hstack([np.full((len(F), 1), 3), F]).ravel())
N = np.cross(emb['e1'], emb['e2'])
pl = pv.Plotter(off_screen=True, shape=(3, 4), window_size=(2400, 1800), border=True)
for k, s in enumerate(pick):
    p = loose[s]; i0 = t.query(p)[1]; li = line[i0]
    pl.subplot(k//4, k % 4); R = 2.0
    box = mesh.clip_box([p[0]-R, p[0]+R, p[1]-R, p[1]+R, p[2]-R, p[2]+R], invert=False)
    pl.add_mesh(box, color='lightgrey', smooth_shading=False)
    near = t.query_ball_point(p, R); lines_here = set(line[near])
    for lj in lines_here:
        a, b = start[lj], start[lj+1]; pts = Q[a:b]+0.01*N[emb['nface'][a:b]]
        pl.add_mesh(pv.PolyData(pts, lines=np.r_[b-a, np.arange(b-a)]), color='blue' if lj == li else 'black', line_width=3 if lj == li else 2)
    pl.add_mesh(pv.PolyData(p.reshape(1, 3)), color='red', point_size=14, render_points_as_spheres=True)
    others = [j for j in near if line[j] != li]
    d = np.linalg.norm(Q[others]-p, axis=1).min() if others else np.inf
    pl.add_text(f'end {k}: nearest other line {d:.2f} mm', font_size=10)
    n = N[emb['nface'][i0]]; pl.camera_position = [p+3.5*R*n, p, [0, 0, 1] if abs(n[2]) < 0.9 else [1, 0, 0]]
out = os.path.join(root, 'out', 'loose_ends_gallery.png'); pl.screenshot(out); print('wrote', out)
