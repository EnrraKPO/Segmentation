"""Clean render of chains over the model: chains lifted off the surface along the normal so chords on convex surface
are not hidden, no dots. usage: python view_chains.py mesh chains.npz out_prefix [lift_mm zooms...]"""
import sys, numpy as np, trimesh, pyvista as pv
from scipy.spatial import cKDTree
path, npz, out = sys.argv[1:4]; a = sys.argv[4:]; lift = float(a[0]) if a else 0.06; zooms = [int(x) for x in a[1:]] or [6, 12, 25]
z = np.load(npz); pts, n = z['pts'], z['n']; chains = np.split(np.arange(len(pts)), np.cumsum(n)[:-1])
m = trimesh.load(path, process=True); P, F = np.asarray(m.vertices), np.asarray(m.faces)
VN = m.vertex_normals; pts = pts+lift*VN[cKDTree(P).query(pts)[1]]
mesh = pv.PolyData(P, np.hstack([np.full((len(F), 1), 3), F]).ravel()); rng = np.random.default_rng(0); cols = rng.uniform(0.1, 0.95, (len(chains), 3))
sc = np.vstack([np.repeat(cols[ci][None], len(ch), 0) for ci, ch in enumerate(chains)]); cl = np.concatenate([np.r_[len(ch), ch] for ch in chains])
focus = np.array([12, 135, 12.])
for zoom in zooms:
    pl = pv.Plotter(off_screen=True, window_size=(1500, 1500)); pl.add_mesh(mesh, color='lightgrey', smooth_shading=True)
    poly = pv.PolyData(pts, lines=cl); poly['c'] = sc; pl.add_mesh(poly, scalars='c', rgb=True, line_width=5, render_lines_as_tubes=True)
    pl.add_text(f'{npz} lifted {lift} mm, zoom {zoom}', font_size=10)
    pl.camera_position = 'zy'; pl.camera.focal_point = focus.tolist(); pl.camera.position = (focus+(np.array(pl.camera.position)-np.array(pl.camera.focal_point))).tolist(); pl.camera.zoom(zoom)
    pl.screenshot(out+f'_zoom{zoom}.png'); print('saved', out+f'_zoom{zoom}.png')
