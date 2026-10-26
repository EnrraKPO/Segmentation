import sys, os, pickle; root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); sys.path.insert(0, root)
import numpy as np, pyvista as pv, regions
res = pickle.load(open(os.path.join(root, 'out', 'paint_dump.pkl'), 'rb')); emb = res['emb']; con = res['con']; P = emb['P']
Q = emb['nodes']; start = emb['start']; N = np.cross(emb['e1'], emb['e2'])
DV, DF, pid = res['DV'], res['DF'], res['pid']; mesh = pv.PolyData(DV, np.hstack([np.full((len(DF), 1), 3), DF]).ravel()); mesh.cell_data['c'] = res['colour'][con['region'][pid]].astype(float)
loose, gaps = regions.loose_ends(emb)
pts_all = []; cells = []; k = 0
for a, b in zip(start[:-1], start[1:]): pts_all.append(Q[a:b]+0.01*N[emb['nface'][a:b]]); cells.append(np.r_[b-a, k+np.arange(b-a)]); k += b-a
lines_pd = pv.PolyData(np.vstack(pts_all), lines=np.concatenate(cells))
top = P[P[:, 2] > np.percentile(P[:, 2], 97)]; c = top.mean(0); print('top centre', np.round(c, 1), 'z range', P[:, 2].min().round(1), P[:, 2].max().round(1))
R = 30
box = [c[0]-R, c[0]+R, c[1]-R, c[1]+R, c[2]-R, c[2]+R]
pl = pv.Plotter(off_screen=True, shape=(2, 2), window_size=(2400, 2400), border=True)
views = [('above', [0, 0, 1], [0, 1, 0]), ('front', [0, -1, 0.3], [0, 0, 1]), ('back', [0, 1, 0.3], [0, 0, 1]), ('side', [1, 0, 0.3], [0, 0, 1])]
for k, (name, d, up) in enumerate(views):
    pl.subplot(k//2, k % 2); d = np.array(d, float); d /= np.linalg.norm(d)
    pl.add_mesh(mesh.clip_box(box, invert=False), scalars='c', cmap=list(regions.PALETTE), clim=(-0.5, 8.5), show_scalar_bar=False)
    pl.add_mesh(lines_pd.clip_box(box, invert=False), color='black', line_width=2)
    near = loose[np.all(np.abs(loose-c) < R, axis=1)]
    if len(near): pl.add_mesh(pv.PolyData(near), color='red', point_size=10, render_points_as_spheres=True)
    pl.camera_position = [c+3.2*R*d, c, up]; pl.add_text(f'hair, {name}', font_size=12)
out = os.path.join(root, 'out', 'hair.png'); pl.screenshot(out); print('wrote', out)
