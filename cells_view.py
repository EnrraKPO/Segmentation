"""Diagnostic: at FOCUS, draw the grid's line cells as dots (coloured by turn) and every trail (any length) at LIFT 0.
usage: SNAP=1 snaps nodes onto the surface. FOCUS=x,y,z VIEW=xz ZOOM=25 python cells_view.py mesh grid.npz chains.npz out.png"""
import sys, os, numpy as np, pyvista as pv, trimesh
M, G, C, outp = sys.argv[1:5]
focus = np.array([float(v) for v in os.environ.get('FOCUS', '0,0,0').split(',')]); zoom = float(os.environ.get('ZOOM', 25)); VIEW = os.environ.get('VIEW', 'xz')
m = trimesh.load(M, process=True); P, F = np.asarray(m.vertices), np.asarray(m.faces)
mesh = pv.PolyData(P, np.hstack([np.full((len(F), 1), 3), F]).ravel())
g = np.load(G); X = g['lo']+(g['idx']+0.5)*g['size']; line = g['line']; turn = g['turn']
R = 120/zoom; near = (np.linalg.norm(X-focus, axis=1) < R) & line
z = np.load(C); pts, n = z['pts'], z['n']; chains = np.split(np.arange(len(pts)), np.cumsum(n)[:-1])
pl = pv.Plotter(off_screen=True, window_size=(1500, 1500)); pl.add_mesh(mesh, color='lightgrey', smooth_shading=True)
pl.camera_position = VIEW; pl.camera.focal_point = focus.tolist(); pl.camera.position = (focus+(np.array(pl.camera.position)-np.array(pl.camera.focal_point))).tolist(); pl.camera.zoom(zoom); pl.camera.Azimuth(float(os.environ.get('AZ', 0)))
d = pv.PolyData(X[near]); d['turn'] = turn[near]
pl.add_mesh(d, scalars='turn', cmap='plasma', point_size=6, render_points_as_spheres=True, reset_camera=False, show_scalar_bar=False)
rng = np.random.default_rng(1); cnt = 0; SNAP = os.environ.get('SNAP', '0') == '1'
for ch in chains:
    if np.linalg.norm(pts[ch]-focus, axis=1).min() > R: continue
    Q = pts[ch]
    if SNAP: Cq, dq, tq = trimesh.proximity.closest_point(m, Q); Q = Cq+0.01*m.face_normals[tq]
    pl.add_mesh(pv.PolyData(Q, lines=np.r_[len(ch), np.arange(len(ch))]), color=rng.uniform(0, 0.8, 3), line_width=4, reset_camera=False); cnt += 1
print('line cells shown', near.sum(), 'trails shown', cnt)
pl.screenshot(outp); print('saved', outp)
