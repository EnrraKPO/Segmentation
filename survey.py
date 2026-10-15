"""Survey: render dots + strokes at several spots, loading the grid once.
usage: python survey.py mesh grid.npz chains.npz out_prefix  (spots and views listed below)"""
import sys, os, numpy as np, pyvista as pv, trimesh
M, G, C, out = sys.argv[1:5]
SPOTS = [('mouth', 'xz', 168.87, 121.33, 177.72), ('pauldron', 'xz', 203.49, 118.88, 167.79), ('buckle', 'xz', 165.0, 138.41, 112.24),
         ('bracer', 'xz', 249.44, 104.57, 126.8), ('axe', 'xz', 85.12, 140.53, 140.23), ('knee', 'xz', 129.77, 163.08, 56.62),
         ('toes', 'xz', 225.43, 108.11, 6.09), ('skirt', 'xz', 184.84, 136.63, 84.44), ('belt', 'xz', 228.11, 115.31, 119.58),
         ('backbelt', 'yz', 213.23, 169.5, 108.08), ('backplate', 'yz', 232.65, 145.32, 177.7)]
ZOOMS = [float(z) for z in os.environ.get('ZOOMS', '6,15').split(',')]
m = trimesh.load(M, process=True); P, F = np.asarray(m.vertices), np.asarray(m.faces)
mesh = pv.PolyData(P, np.hstack([np.full((len(F), 1), 3), F]).ravel())
g = np.load(G); X = g['lo']+(g['idx']+0.5)*g['size']; line = g['line']; turn = g['turn']; X = X[line]; turn = turn[line]
z = np.load(C); pts, n = z['pts'], z['n']; chains = np.split(np.arange(len(pts)), np.cumsum(n)[:-1])
Cq, dq, tq = trimesh.proximity.closest_point(m, pts); spts = Cq+0.01*m.face_normals[tq]
for name, VIEW, *f in SPOTS:
    focus = np.array(f)
    for zoom in ZOOMS:
        R = 120/zoom; near = np.linalg.norm(X-focus, axis=1) < R
        pl = pv.Plotter(off_screen=True, window_size=(1500, 1500)); pl.add_mesh(mesh, color='lightgrey', smooth_shading=True)
        pl.camera_position = VIEW; pl.camera.focal_point = focus.tolist(); pl.camera.position = (focus+(np.array(pl.camera.position)-np.array(pl.camera.focal_point))).tolist(); pl.camera.zoom(zoom)
        d = pv.PolyData(X[near]); d['turn'] = turn[near]
        pl.add_mesh(d, scalars='turn', cmap='plasma', point_size=4, render_points_as_spheres=True, reset_camera=False, show_scalar_bar=False)
        rng = np.random.default_rng(1); cnt = 0
        for ch in chains:
            if np.linalg.norm(pts[ch]-focus, axis=1).min() > R: continue
            pl.add_mesh(pv.PolyData(spts[ch], lines=np.r_[len(ch), np.arange(len(ch))]), color=rng.uniform(0, 0.8, 3), line_width=3, reset_camera=False); cnt += 1
        o = f'{out}_{name}_z{zoom:g}.png'; pl.screenshot(o); print(name, zoom, 'dots', int(near.sum()), 'strokes', cnt, o, flush=True)
