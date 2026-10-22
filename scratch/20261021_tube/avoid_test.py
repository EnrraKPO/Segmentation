import sys, os, numpy as np, trimesh, pyvista as pv
from scipy.spatial import cKDTree
sys.path.insert(0, r'D:\Stuff\PrintSuit\SegmentationTool'); os.chdir(r'D:\Stuff\PrintSuit\SegmentationTool'); import trace4, refine
reach = float(sys.argv[1]) if len(sys.argv) > 1 else 0.174
g = np.load('out/tube_10_0.3.npz'); size = float(g['size']); L = g['coarse']; band = g['line']; X = g['lo']+(g['idx'][band]+0.5)*size; turn = g['turn'][band]; N = g['N'][band]
dmp = np.load('out/debug/20261020_144432/lines.npz'); n = dmp['n']; pts = dmp['pts']; off = np.r_[0, np.cumsum(n)]
curves = [pts[off[i]:off[i+1]] for i in range(len(n))]; me = refine.nearest_line(curves, L[len(L)//2])
others = np.vstack([refine.resample(c, 0.02) for i, c in enumerate(curves) if i != me and cKDTree(L).query(c)[0].min() < 1.0])
dist = cKDTree(others).query(X)[0]; keep = dist > reach
print(f'line {me}; other lines within 3 mm: resampled to {len(others)} points; fine dots {len(X)}, within {reach} mm of another line: {(~keep).sum()}, kept {keep.sum()}')
Xk, tk, Nk = X[keep], turn[keep], N[keep]; cut = 17.0
lines, _ = trace4.walk((Xk-g['lo'])/size, np.maximum(tk-cut, 0)+1e-3, Nk, size, log=lambda *a: None, step=0.087, slice=0.7, reach=0.087, pull=1.0, look=0.45, levp=2.0, minlen=8)
Rs = [g['lo']+l*size for l in lines]; R = max(Rs, key=len); dd = cKDTree(Xk).query(R)[0]
v = np.diff(R, axis=0); v /= np.linalg.norm(v, axis=1, keepdims=True)+1e-12; tr = np.degrees(np.arccos(np.clip((v[1:]*v[:-1]).sum(1), -1, 1)))
d = np.r_[0, np.cumsum(np.linalg.norm(np.diff(R, axis=0), axis=1))]
print(f'walk: pieces {[len(l) for l in lines]}, length {d[-1]:.1f} mm; turn p50 {np.percentile(tr, 50):.1f} p90 {np.percentile(tr, 90):.1f} max {tr.max():.1f}; nodes off the band {(dd > 0.015).sum()}, farthest {dd.max():.3f} mm')
np.savez('out/tube_avoid.npz', keep=keep, line=R, others=others)
m = trimesh.load('out/orcfull_patch.stl', process=False); P = np.asarray(m.vertices); F = np.asarray(m.faces)
lift = lambda Q: (lambda r: r[0]+0.01*m.face_normals[r[2]])(trimesh.proximity.closest_point(m, Q)); Xl = lift(X); Rl = lift(R); Ll = lift(L)
for name, at in [('A', 0.6), ('B', d[-1]-0.6)]:
    focus = np.array([np.interp(at, d, R[:, i]) for i in range(3)]); cen = P[F].mean(1); near = (np.abs(cen-focus) < 2).all(1); Fp = F[near]
    mesh = pv.PolyData(P, np.c_[np.full(len(Fp), 3), Fp].ravel()); pl = pv.Plotter(off_screen=True, window_size=(1600, 1200)); pl.add_mesh(mesh, color='white', smooth_shading=True); pl.add_mesh(mesh, style='wireframe', color='black', line_width=1, opacity=0.25)
    pl.add_mesh(pv.PolyData(Xl[~keep]), color='black', point_size=4, render_points_as_spheres=True, reset_camera=False)
    pts_ = pv.PolyData(Xl[keep]); pts_['turn'] = turn[keep]; pl.add_mesh(pts_, scalars='turn', cmap='viridis', clim=(17, 45), point_size=4, render_points_as_spheres=True, show_scalar_bar=False, reset_camera=False)
    for i, c in enumerate(curves):
        if i != me and np.linalg.norm(c-focus, axis=1).min() < 2: cl = lift(c); pl.add_mesh(pv.PolyData(cl, lines=np.r_[len(cl), np.arange(len(cl))]), color='blue', line_width=3, reset_camera=False)
    pl.add_mesh(pv.PolyData(Rl, lines=np.r_[len(Rl), np.arange(len(Rl))]), color='red', line_width=5, reset_camera=False); pl.add_mesh(pv.PolyData(Rl), color='white', point_size=7, render_points_as_spheres=True, reset_camera=False)
    vv = m.face_normals[near].mean(0); vv /= np.linalg.norm(vv); pl.camera.focal_point = focus.tolist(); pl.camera.position = (focus+vv*30).tolist(); pl.camera.up = (0, 0, 1); pl.camera.parallel_projection = True; pl.camera.parallel_scale = 0.9
    pl.add_text(f'end {name}: black = dots within {reach} mm of another line (blue), removed', font_size=12); pl.screenshot(f'out/avoid_{name}.png')
print('saved out/avoid_A.png out/avoid_B.png')
