"""bands in colour: each band's cells as dots and its curve as a line, same colour. usage: bands_view.py mesh npz out.png view zoom [fx fy fz]"""
import sys, numpy as np, trimesh, pyvista as pv
path, npz, out, view, zoom = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], float(sys.argv[5])
focus = np.array([float(v) for v in sys.argv[6:9]]) if len(sys.argv) >= 9 else None
m = trimesh.load(path, process=True); P, F = np.asarray(m.vertices), np.asarray(m.faces)
z = np.load(npz); n, pts, bn, bpts = z['n'], z['pts'], z['bn'], z['bpts']
rng = np.random.default_rng(0); cols = rng.uniform(0.1, 0.9, (len(n), 3))
mesh = pv.PolyData(P, np.hstack([np.full((len(F), 1), 3), F]).ravel())
pl = pv.Plotter(off_screen=True, shape=(1, 2), window_size=(3000, 1500))
boff = np.cumsum(np.r_[0, bn[:-1]]); coff = np.cumsum(np.r_[0, n[:-1]])
bc = np.repeat(cols, bn, axis=0); cc = np.repeat(cols, n, axis=0)
for i in (0, 1):
    pl.subplot(0, i); pl.add_mesh(mesh, color='lightgrey', smooth_shading=True)
    if i == 0:
        pd = pv.PolyData(bpts); pd['c'] = bc; pl.add_mesh(pd, scalars='c', rgb=True, point_size=6, render_points_as_spheres=True)
        pl.add_text('band cells, one colour per band', font_size=10)
    else:
        cl = np.concatenate([np.r_[k, np.arange(o, o+k)] for k, o in zip(n, coff)]); pd = pv.PolyData(pts, lines=cl); pd['c'] = cc
        pl.add_mesh(pd, scalars='c', rgb=True, line_width=4); pl.add_text('curves, same colours', font_size=10)
    pl.camera_position = view
    if focus is not None: pl.camera.focal_point = focus.tolist(); pl.camera.position = (focus+(np.array(pl.camera.position)-np.array(pl.camera.focal_point))).tolist()
    pl.camera.zoom(zoom)
pl.screenshot(out); print('saved', out)
