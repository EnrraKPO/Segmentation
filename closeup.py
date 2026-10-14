"""Close-up of saved curves on the model: usage python closeup.py mesh.stl curves.npz out.png view zoom [focus_x focus_y focus_z]
view: xy yx xz zx yz zy (pyvista camera positions). zoom: e.g. 3. focus: world point to centre on (default model centre)."""
import sys, os, numpy as np, trimesh, pyvista as pv
path, npz, out, view, zoom = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], float(sys.argv[5])
focus = np.array([float(v) for v in sys.argv[6:9]]) if len(sys.argv) >= 9 else None
m = trimesh.load(path, process=True); P, F = np.asarray(m.vertices), np.asarray(m.faces)
z = np.load(npz); n, pts = z['n'], z['pts']; off = np.cumsum(np.r_[0, n[:-1]])
cl = np.concatenate([np.r_[k, np.arange(o, o+k)] for k, o in zip(n, off)]); poly = pv.PolyData(pts, lines=cl)
mesh = pv.PolyData(P, np.hstack([np.full((len(F), 1), 3), F]).ravel())
pl = pv.Plotter(off_screen=True, shape=(1, 2), window_size=(3000, 1500))
for i in (0, 1):
    pl.subplot(0, i); pl.add_mesh(mesh, color='lightgrey', smooth_shading=True)
    if i == 1: pl.add_mesh(poly, color='red', line_width=3)
    pl.camera_position = view
    if focus is not None: pl.camera.focal_point = focus.tolist(); pl.camera.position = (focus+(np.array(pl.camera.position)-np.array(pl.camera.focal_point))).tolist()
    pl.camera.zoom(zoom)
pl.screenshot(out); print('saved', out)
