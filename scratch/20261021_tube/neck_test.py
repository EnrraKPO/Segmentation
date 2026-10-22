import sys, os, numpy as np, trimesh, pyvista as pv
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
sys.path.insert(0, r'D:\Stuff\PrintSuit\SegmentationTool'); os.chdir(r'D:\Stuff\PrintSuit\SegmentationTool'); import trace4
g = np.load('out/tube_10_0.3.npz'); size = float(g['size']); L = g['coarse']; band = g['line']; Xall = g['lo']+(g['idx']+0.5)*size
X = Xall[band]; turn = g['turn'][band]; N = g['N'][band]; n0 = len(X)
rho = float(sys.argv[1]) if len(sys.argv) > 1 else 0.05; fill = float(sys.argv[2]) if len(sys.argv) > 2 else 0.6
tree = cKDTree(X); tall = cKDTree(Xall)
cnt = np.array([len(b) for b in tree.query_ball_point(X, rho)]); full = np.array([len(b) for b in tall.query_ball_point(X, rho)])   # band dots vs all surface cubes in the disc
body = cnt >= fill*full
print(f'band {n0} dots; disc radius {rho} mm; body (band fills >= {fill:.0%} of the disc): {body.sum()}')
Xb = X[body]; pr = cKDTree(Xb).query_pairs(1.8*size, output_type='ndarray')
lab = connected_components(coo_matrix((np.ones(len(pr), bool), (pr[:, 0], pr[:, 1])), shape=(len(Xb), len(Xb))), directed=False)[1]
mid = L[len(L)//2]; seed = lab[np.argmin(np.linalg.norm(Xb-mid, axis=1))]; core = lab == seed
sizes = np.bincount(lab); print(f'body pieces: {len(sizes)}, sizes of the largest 6: {sorted(sizes, reverse=True)[:6]}; the piece under the line middle: {core.sum()} dots')
keep = cKDTree(Xb[core]).query(X)[0] <= rho
print(f'kept after taking the edge back: {keep.sum()} of {n0} ({n0-keep.sum()} severed)')
Xk, tk, Nk = X[keep], turn[keep], N[keep]; cut = 17.0
lines, _ = trace4.walk((Xk-g['lo'])/size, np.maximum(tk-cut, 0)+1e-3, Nk, size, log=lambda *a: None, step=0.087, slice=0.7, reach=0.087, pull=1.0, look=0.45, levp=2.0, minlen=8)
Rs = [g['lo']+l*size for l in lines]; R = max(Rs, key=len); dd = cKDTree(Xk).query(R)[0]
v = np.diff(R, axis=0); v /= np.linalg.norm(v, axis=1, keepdims=True)+1e-12; tr = np.degrees(np.arccos(np.clip((v[1:]*v[:-1]).sum(1), -1, 1)))
print(f'walk: pieces {[len(l) for l in lines]}; turn p50 {np.percentile(tr, 50):.1f} p90 {np.percentile(tr, 90):.1f} max {tr.max():.1f}; nodes off the band {(dd > 0.015).sum()}, farthest {dd.max():.3f} mm')
np.savez('out/tube_neck.npz', keep=keep, line=R)
m = trimesh.load('out/orcfull_patch.stl', process=False); P = np.asarray(m.vertices); F = np.asarray(m.faces)
lift = lambda Q: (lambda r: r[0]+0.01*m.face_normals[r[2]])(trimesh.proximity.closest_point(m, Q)); Xl = lift(X); Rl = lift(R); Ll = lift(L)
d = np.r_[0, np.cumsum(np.linalg.norm(np.diff(R, axis=0), axis=1))]
for name, at in [('A', 0.8), ('B', d[-1]-0.8)]:
    focus = np.array([np.interp(at, d, R[:, i]) for i in range(3)]); cen = P[F].mean(1); near = (np.abs(cen-focus) < 2).all(1); Fp = F[near]
    mesh = pv.PolyData(P, np.c_[np.full(len(Fp), 3), Fp].ravel()); pl = pv.Plotter(off_screen=True, window_size=(1600, 1200)); pl.add_mesh(mesh, color='lightgrey', smooth_shading=True); pl.add_mesh(mesh, style='wireframe', color='black', line_width=1, opacity=0.3)
    pl.add_mesh(pv.PolyData(Xl[~keep]), color='black', point_size=4, render_points_as_spheres=True, reset_camera=False)       # severed dots black
    pts = pv.PolyData(Xl[keep]); pts['turn'] = turn[keep]; pl.add_mesh(pts, scalars='turn', cmap='plasma', point_size=3, render_points_as_spheres=True, show_scalar_bar=False, reset_camera=False)
    pl.add_mesh(pv.PolyData(Ll, lines=np.r_[len(Ll), np.arange(len(Ll))]), color='blue', line_width=3, reset_camera=False)
    pl.add_mesh(pv.PolyData(Rl, lines=np.r_[len(Rl), np.arange(len(Rl))]), color='red', line_width=5, reset_camera=False); pl.add_mesh(pv.PolyData(Rl), color='white', point_size=7, render_points_as_spheres=True, reset_camera=False)
    vv = m.face_normals[near].mean(0); vv /= np.linalg.norm(vv); pl.camera.focal_point = focus.tolist(); pl.camera.position = (focus+vv*30).tolist(); pl.camera.up = (0, 0, 1); pl.camera.parallel_projection = True; pl.camera.parallel_scale = 0.9
    pl.screenshot(f'out/neck_{name}.png')
print('saved out/neck_A.png out/neck_B.png (black = severed dots)')
