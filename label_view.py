"""Close-up with every line numbered at both ends (and dots coloured by strength), lines snapped to the surface.
usage: FOCUS=x,y,z VIEW=xz ZOOM=30 AZ=0 python label_view.py mesh grid.npz chains.npz out.png"""
import sys, os, numpy as np, pyvista as pv, trimesh
M, G, C, outp = sys.argv[1:5]
focus = np.array([float(v) for v in os.environ.get('FOCUS', '0,0,0').split(',')]); zoom = float(os.environ.get('ZOOM', 30)); VIEW = os.environ.get('VIEW', 'xz')
m = trimesh.load(M, process=True); P, F = np.asarray(m.vertices), np.asarray(m.faces)
mesh = pv.PolyData(P, np.hstack([np.full((len(F), 1), 3), F]).ravel())
g = np.load(G); X = g['lo']+(g['idx']+0.5)*g['size']; line = g['line']; turn = g['turn']; X = X[line]; turn = turn[line]
z = np.load(C); pts, n = z['pts'], z['n']; chains = np.split(np.arange(len(pts)), np.cumsum(n)[:-1])
R = 120/zoom; near = np.linalg.norm(X-focus, axis=1) < R
pl = pv.Plotter(off_screen=True, window_size=(1500, 1500)); pl.add_mesh(mesh, color='lightgrey', smooth_shading=True)
pl.camera_position = VIEW; pl.camera.focal_point = focus.tolist(); pl.camera.position = (focus+(np.array(pl.camera.position)-np.array(pl.camera.focal_point))).tolist(); pl.camera.zoom(zoom); pl.camera.Azimuth(float(os.environ.get('AZ', 0)))
d = pv.PolyData(X[near]); d['turn'] = turn[near]
pl.add_mesh(d, scalars='turn', cmap='plasma', point_size=6, render_points_as_spheres=True, reset_camera=False, show_scalar_bar=False)
rng = np.random.default_rng(1); lab_pts = []; labs = []
for i, ch in enumerate(chains):
    col = rng.uniform(0, 0.8, 3)
    if len(ch) == 0 or np.linalg.norm(pts[ch]-focus, axis=1).min() > R: continue
    Q = pts[ch]; Cq, dq, tq = trimesh.proximity.closest_point(m, Q); Qs = Cq+0.01*m.face_normals[tq]
    pl.add_mesh(pv.PolyData(Qs, lines=np.r_[len(ch), np.arange(len(ch))]), color=col, line_width=4, reset_camera=False)
    for k in (0, len(ch)-1):
        if np.linalg.norm(Q[k]-focus) < R: lab_pts.append(Qs[k]+0.05*m.face_normals[tq[k]]); labs.append(f'{i}{"s" if k == 0 else "e"}')
    mid = len(ch)//2
    if np.linalg.norm(Q[mid]-focus) < R: lab_pts.append(Qs[mid]+0.05*m.face_normals[tq[mid]]); labs.append(f'{i}')
if lab_pts: pl.add_point_labels(np.array(lab_pts), labs, font_size=22, text_color='black', point_color='red', point_size=8, shape=None, always_visible=True)
pl.screenshot(outp); print('saved', outp)
