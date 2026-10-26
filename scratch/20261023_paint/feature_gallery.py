"""close-ups of the painted orc around the strongest (longest) walk lines: regions coloured, lines black, loose ends red"""
import sys, os, pickle; root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); sys.path.insert(0, root)
import numpy as np, pyvista as pv, regions
res = pickle.load(open(os.path.join(root, 'out', 'paint_dump.pkl'), 'rb')); emb = res['emb']; con = res['con']
Q = emb['nodes']; start = emb['start']; n = np.diff(start); N = np.cross(emb['e1'], emb['e2'])
DV, DF, pid = res['DV'], res['DF'], res['pid']; mesh = pv.PolyData(DV, np.hstack([np.full((len(DF), 1), 3), DF]).ravel()); mesh.cell_data['c'] = res['colour'][con['region'][pid]].astype(float)
loose, gaps = regions.loose_ends(emb)
order = np.argsort(n)[::-1]; R = float(sys.argv[1]) if len(sys.argv) > 1 else 8.0
pts_all = []; cells = []; k = 0
for a, b in zip(start[:-1], start[1:]): pts_all.append(Q[a:b]+0.01*N[emb['nface'][a:b]]); cells.append(np.r_[b-a, k+np.arange(b-a)]); k += b-a
lines_pd = pv.PolyData(np.vstack(pts_all), lines=np.concatenate(cells))
pl = pv.Plotter(off_screen=True, shape=(2, 3), window_size=(2700, 1800), border=True)
for k, li in enumerate(order[:6]):
    a, b = start[li], start[li+1]; p = Q[(a+b)//2]; pl.subplot(k//3, k % 3)
    box = [p[0]-R, p[0]+R, p[1]-R, p[1]+R, p[2]-R, p[2]+R]
    pl.add_mesh(mesh.clip_box(box, invert=False), scalars='c', cmap=list(regions.PALETTE), clim=(-0.5, 8.5), show_scalar_bar=False)
    pl.add_mesh(lines_pd.clip_box(box, invert=False), color='black', line_width=2)
    near = loose[np.all(np.abs(loose-p) < R, axis=1)]
    if len(near): pl.add_mesh(pv.PolyData(near), color='red', point_size=12, render_points_as_spheres=True)
    nn = N[emb['nface'][(a+b)//2]]; pl.camera_position = [p+3.0*R*nn, p, [0, 0, 1] if abs(nn[2]) < 0.9 else [1, 0, 0]]
    pl.add_text(f'line {li}, {n[li]} nodes; {len(near)} loose ends in view', font_size=11)
out = os.path.join(root, 'out', f'feature_gallery_{int(R)}.png'); pl.screenshot(out); print('wrote', out)
