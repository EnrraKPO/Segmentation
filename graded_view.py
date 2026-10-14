"""All surface cubes near FOCUS as dots coloured by a field (turn), like a shade. usage: FOCUS VIEW ZOOM AZ MIN=deg python graded_view.py mesh grid.npz out.png"""
import sys, os, numpy as np, pyvista as pv, trimesh
M, G, outp = sys.argv[1:4]; focus = np.array([float(v) for v in os.environ['FOCUS'].split(',')]); zoom = float(os.environ.get('ZOOM', 10)); VIEW = os.environ.get('VIEW', 'xz'); MIN = float(os.environ.get('MIN', 0))
m = trimesh.load(M, process=True); P, F = np.asarray(m.vertices), np.asarray(m.faces); mesh = pv.PolyData(P, np.hstack([np.full((len(F), 1), 3), F]).ravel())
g = np.load(G); X = g['lo']+(g['idx']+0.5)*g['size']; t = g['turn']; near = (np.linalg.norm(X-focus, axis=1) < 120/zoom) & (t >= MIN)
pl = pv.Plotter(off_screen=True, window_size=(1500, 1500)); pl.add_mesh(mesh, color='lightgrey', smooth_shading=True)
pl.camera_position = VIEW; pl.camera.focal_point = focus.tolist(); pl.camera.position = (focus+(np.array(pl.camera.position)-np.array(pl.camera.focal_point))).tolist(); pl.camera.zoom(zoom); pl.camera.Azimuth(float(os.environ.get('AZ', 0)))
d = pv.PolyData(X[near]); d['v'] = t[near]; pl.add_mesh(d, scalars='v', cmap='plasma', clim=[MIN, max(MIN+1, np.percentile(t[near], 98))], point_size=7, render_points_as_spheres=True, reset_camera=False, scalar_bar_args=dict(title='deg'))
pl.screenshot(outp); print('saved', outp, near.sum(), 'cubes')
