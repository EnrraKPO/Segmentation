import sys, os, numpy as np, trimesh, pyvista as pv
os.chdir(r'D:\Stuff\PrintSuit\SegmentationTool')
npz = sys.argv[1]; outp = sys.argv[2]; show_dots = len(sys.argv) > 3 and sys.argv[3] == 'dots'; zoom = float(sys.argv[4]) if len(sys.argv) > 4 else 0.6; at = float(sys.argv[5]) if len(sys.argv) > 5 else -1
g = np.load(npz); L = g['coarse']; n = g['n']; pts = g['pts']; off = np.r_[0, np.cumsum(n)]
m = trimesh.load('out/orcfull_patch.stl', process=False); P = np.asarray(m.vertices); F = np.asarray(m.faces)
cen = P[F].mean(1); near = (np.abs(cen-L.mean(0)) < 5).all(1); Fp = F[near]
mesh = pv.PolyData(P, np.c_[np.full(len(Fp), 3), Fp].ravel())
pl = pv.Plotter(off_screen=True, window_size=(1600, 1200)); pl.add_mesh(mesh, color='lightgrey', smooth_shading=True)
def lift(Q): C, _, t = trimesh.proximity.closest_point(m, Q); return C+0.01*m.face_normals[t]
Lq = lift(L); pl.add_mesh(pv.PolyData(Lq, lines=np.r_[len(Lq), np.arange(len(Lq))]), color='blue', line_width=4, reset_camera=False)
for i in range(len(n)):
    ch = lift(pts[off[i]:off[i+1]]); pl.add_mesh(pv.PolyData(ch, lines=np.r_[len(ch), np.arange(len(ch))]), color='red', line_width=4, reset_camera=False)
    pl.add_mesh(pv.PolyData(ch[[0, -1]]), color='black', point_size=10, render_points_as_spheres=True, reset_camera=False)
if show_dots:
    X = g['lo']+(g['idx'][g['line']]+0.5)*float(g['size']); d = pv.PolyData(lift(X)); d['turn'] = g['turn'][g['line']]
    pl.add_mesh(d, scalars='turn', cmap='plasma', point_size=3, render_points_as_spheres=True, reset_camera=False, show_scalar_bar=False)
# view: from the mean facing of the surface near the line, looking at the line's centre
nn = m.face_normals[near]; v = nn.mean(0); v /= np.linalg.norm(v); d = np.r_[0, np.cumsum(np.linalg.norm(np.diff(L, axis=0), axis=1))]; focus = L.mean(0) if at < 0 else np.array([np.interp(at*d[-1], d, L[:, i]) for i in range(3)]); ext = np.linalg.norm(L.max(0)-L.min(0))
pl.camera.focal_point = focus.tolist(); pl.camera.position = (focus+v*3*ext).tolist(); pl.camera.up = (0, 0, 1); pl.camera.parallel_projection = True; pl.camera.parallel_scale = ext*zoom
pl.screenshot(outp); print('saved', outp)
