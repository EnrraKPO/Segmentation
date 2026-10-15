"""Close-up at FOCUS with dots (plasma by turn) and strokes drawn thick, saves PNG and prints wobble stats of strokes in view.
usage: FOCUS=x,y,z VIEW=xz ZOOM=40 python closeup2.py mesh grid.npz chains.npz out.png"""
import sys, os, numpy as np, pyvista as pv, trimesh
M, G, C, outp = sys.argv[1:5]
focus = np.array([float(v) for v in os.environ.get('FOCUS', '0,0,0').split(',')]); zoom = float(os.environ.get('ZOOM', 40)); VIEW = os.environ.get('VIEW', 'xz')
m = trimesh.load(M, process=True); P, F = np.asarray(m.vertices), np.asarray(m.faces)
mesh = pv.PolyData(P, np.hstack([np.full((len(F), 1), 3), F]).ravel())
g = np.load(G); X = g['lo']+(g['idx']+0.5)*g['size']; line = g['line']; turn = g['turn']; X = X[line]; turn = turn[line]; size = float(g['size'])
z = np.load(C); pts, n = z['pts'], z['n']; chains = np.split(np.arange(len(pts)), np.cumsum(n)[:-1])
R = 120/zoom; near = np.linalg.norm(X-focus, axis=1) < R
pl = pv.Plotter(off_screen=True, window_size=(1500, 1500)); pl.add_mesh(mesh, color='lightgrey', smooth_shading=True)
pl.camera_position = VIEW; pl.camera.focal_point = focus.tolist(); pl.camera.position = (focus+(np.array(pl.camera.position)-np.array(pl.camera.focal_point))).tolist(); pl.camera.zoom(zoom)
d = pv.PolyData(X[near]); d['turn'] = turn[near]
pl.add_mesh(d, scalars='turn', cmap='plasma', point_size=7, render_points_as_spheres=True, reset_camera=False, show_scalar_bar=False)
rng = np.random.default_rng(1)
for ch in chains:
    if len(ch) == 0: continue
    if np.linalg.norm(pts[ch]-focus, axis=1).min() > R: continue
    Q = pts[ch]; Cq, dq, tq = trimesh.proximity.closest_point(m, Q); Qs = Cq+0.01*m.face_normals[tq]
    seg = np.diff(Q, axis=0); seg /= np.linalg.norm(seg, axis=1, keepdims=True)+1e-12
    ang = np.degrees(np.arccos(np.clip((seg[:-1]*seg[1:]).sum(1), -1, 1))) if len(seg) > 1 else np.zeros(1)
    print(f'stroke {len(ch):4d} nodes, heading change per node: median {np.median(ang):.1f} p90 {np.percentile(ang, 90):.1f} deg; off surface median {np.median(dq)/size:.2f} cubes')
    pl.add_mesh(pv.PolyData(Qs, lines=np.r_[len(ch), np.arange(len(ch))]), color=rng.uniform(0, 0.8, 3), line_width=5, reset_camera=False)
print('dots in view', int(near.sum()), 'cube', size); pl.screenshot(outp); print('saved', outp)
