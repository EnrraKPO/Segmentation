"""Verification without the window: run the method on a model and save full views, shaded next to the curves.
usage: python check.py mesh.stl out.png [cells sensitivity reach agree beside merge_reach pull gap cycles min_len]"""
import sys, os, time, numpy as np, trimesh
import pyvista as pv
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import voxel, lines
path, out = sys.argv[1], sys.argv[2]
cells, sens = (int(sys.argv[3]) if len(sys.argv) > 3 else 600), (float(sys.argv[4]) if len(sys.argv) > 4 else 6.0)
kw = dict(zip(['reach', 'agree', 'beside', 'merge_reach', 'pull', 'gap', 'cycles', 'min_len'], [float(v) for v in sys.argv[5:13]]))
t0 = time.time()
m = trimesh.load(path, process=True); P, F = np.asarray(m.vertices, np.float64), np.asarray(m.faces, np.int64)
res = voxel.run(P, F, cells, sens, print); curves = lines.run(P, F, res, log=print, **kw); print(f'{time.time()-t0:.0f} s')
np.savez(os.path.splitext(out)[0]+'_curves.npz', n=np.array([len(c) for c in curves]), pts=np.vstack(curves), bn=np.array([len(b) for b in lines.run.groups]), bpts=np.vstack(lines.run.groups))
mesh = pv.PolyData(P, np.hstack([np.full((len(F), 1), 3), F]).ravel())
pts = np.vstack(curves); n = [len(c) for c in curves]; off = np.cumsum([0]+n[:-1])
cl = np.concatenate([np.r_[k, np.arange(o, o+k)] for k, o in zip(n, off)]); poly = pv.PolyData(pts, lines=cl)
ext = P.max(0)-P.min(0); up = 'y' if ext[1] >= ext[2] else 'z'
views = ['xy', 'yx', 'xz', 'zx'] if up == 'z' else ['xy', 'yx', 'zy', 'yz']
pl = pv.Plotter(off_screen=True, shape=(2, 4), window_size=(4000, 2000))
for i, v in enumerate(views):
    for row in (0, 1):
        pl.subplot(row, i); pl.add_mesh(mesh, color='lightgrey', smooth_shading=True)
        if row == 1: pl.add_mesh(poly, color='red', line_width=4)
        pl.add_text(('shaded ' if row == 0 else 'curves ')+v, font_size=10)
        pl.camera_position = v; pl.camera.zoom(1.3)
pl.screenshot(out); print('saved', out)
