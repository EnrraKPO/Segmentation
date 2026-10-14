"""Close-up with the mesh triangle edges drawn, plus the line cubes (dots by value) and trails. Shows where the dots
sit relative to the triangulation. usage: FOCUS VIEW ZOOM python wire_view.py mesh grid.npz chains.npz out.png"""
import sys, os, numpy as np, pyvista as pv, trimesh
M, G, C, outp = sys.argv[1:5]; focus = np.array([float(v) for v in os.environ['FOCUS'].split(',')]); zoom = float(os.environ.get('ZOOM', 30)); VIEW = os.environ.get('VIEW', 'xz')
m = trimesh.load(M, process=True); P, F = np.asarray(m.vertices), np.asarray(m.faces)
near_f = np.linalg.norm(P[F].mean(1)-focus, axis=1) < 150/zoom; Fn = F[near_f]
mesh = pv.PolyData(P, np.hstack([np.full((len(F), 1), 3), F]).ravel()); sub = pv.PolyData(P, np.hstack([np.full((len(Fn), 1), 3), Fn]).ravel())
g = np.load(G); X = g['lo']+(g['idx']+0.5)*g['size']; R = 120/zoom; near = (np.linalg.norm(X-focus, axis=1) < R) & g['line']
z = np.load(C); pts, n = z['pts'], z['n']; chains = np.split(np.arange(len(pts)), np.cumsum(n)[:-1])
pl = pv.Plotter(off_screen=True, window_size=(1500, 1500)); pl.add_mesh(mesh, color='lightgrey', smooth_shading=True)
pl.camera_position = VIEW; pl.camera.focal_point = focus.tolist(); pl.camera.position = (focus+(np.array(pl.camera.position)-np.array(pl.camera.focal_point))).tolist(); pl.camera.zoom(zoom)
pl.add_mesh(sub, style='wireframe', color='black', line_width=1, opacity=0.5, reset_camera=False)
Cq, dq, tq = trimesh.proximity.closest_point(m, X[near]); d = pv.PolyData(Cq+0.01*m.face_normals[tq]); d['v'] = g['turn'][near]
pl.add_mesh(d, scalars='v', cmap='plasma', point_size=9, render_points_as_spheres=True, reset_camera=False, show_scalar_bar=False)
for ch in chains:
    if np.linalg.norm(pts[ch]-focus, axis=1).min() > R: continue
    Cq, dq, tq = trimesh.proximity.closest_point(m, pts[ch]); pl.add_mesh(pv.PolyData(Cq+0.01*m.face_normals[tq], lines=np.r_[len(ch), np.arange(len(ch))]), color='lime', line_width=4, reset_camera=False)
pl.screenshot(outp); print('saved', outp, 'triangles drawn', len(Fn), 'dots', near.sum())
