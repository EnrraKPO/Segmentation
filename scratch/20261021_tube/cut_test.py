import os, numpy as np, trimesh, pyvista as pv
from scipy.spatial import cKDTree
os.chdir(r'D:\Stuff\PrintSuit\SegmentationTool'); g = np.load('out/tube_10_0.3.npz'); size = float(g['size']); R = g['pts']; L = g['coarse']
K = len(R); d = np.r_[0, np.cumsum(np.linalg.norm(np.diff(R, axis=0), axis=1))]
H = 3; T = 6                                                   # heading over 3 nodes (0.26 mm); trend = mean heading of the 6 nodes before (0.5 mm)
head = np.zeros((K, 3)); head[H:] = R[H:]-R[:-H]; head[:H] = head[H]; head /= np.linalg.norm(head, axis=1, keepdims=True)+1e-12
def dev_out(i, sgn):                                           # heading at i against the trend of the T nodes behind it (toward the middle)
    back = head[i-sgn*T:i:sgn] if sgn > 0 else head[i+1:i+1+T]; t = back.mean(0); t /= np.linalg.norm(t)+1e-12
    return np.degrees(np.arccos(np.clip(head[i] @ t, -1, 1)))
mid = np.arange(int(0.2*K), int(0.8*K)); dm = np.array([dev_out(i, 1) for i in mid]); allow = dm.max()
print(f'{K} nodes, {d[-1]:.1f} mm. Middle 60 %: heading departs from its trend by p50 {np.median(dm):.1f}, max {allow:.1f} deg -> allowed deviation {allow:.1f} deg')
cuts = {}
for name, rng, sgn in (('end A (node 0)', range(mid[0], T, -1), -1), ('end B (last node)', range(mid[-1], K-H), 1)):
    dv = [(i, dev_out(i, sgn)) for i in rng]; bad = [i for i, v in dv if v > allow]
    cuts[name] = bad[0] if bad else None
    print(f'{name}: first node beyond the allowance: {bad[0] if bad else "none"} at {d[bad[0]] if bad else -1:.2f} mm (deviation {dict(dv)[bad[0]] if bad else 0:.0f} deg); nodes cut: {(bad[0] if bad else 0) if sgn < 0 else (K-bad[0] if bad else 0)}')
ca = cuts['end A (node 0)']; cb = cuts['end B (last node)']; keep = R[(ca if ca is not None else 0):(cb if cb is not None else K)]
np.savez('out/tube_cut.npz', keep=keep, cutA=ca if ca is not None else -1, cutB=cb if cb is not None else -1); print('kept', len(keep), 'nodes')
m = trimesh.load('out/orcfull_patch.stl', process=False); P = np.asarray(m.vertices); F = np.asarray(m.faces)
lift = lambda Q: (lambda r: r[0]+0.01*m.face_normals[r[2]])(trimesh.proximity.closest_point(m, Q))
X = g['lo']+(g['idx'][g['line']]+0.5)*size; Xl = lift(X); Rl = lift(R); Kl = lift(keep); Ll = lift(L)
for name, at in [('A', 0.8), ('B', d[-1]-0.8)]:
    focus = np.array([np.interp(at, d, R[:, i]) for i in range(3)]); cen = P[F].mean(1); near = (np.abs(cen-focus) < 2).all(1); Fp = F[near]
    mesh = pv.PolyData(P, np.c_[np.full(len(Fp), 3), Fp].ravel()); pl = pv.Plotter(off_screen=True, window_size=(1600, 1200)); pl.add_mesh(mesh, color='lightgrey', smooth_shading=True); pl.add_mesh(mesh, style='wireframe', color='black', line_width=1, opacity=0.3)
    pts = pv.PolyData(Xl); pts['turn'] = g['turn'][g['line']]; pl.add_mesh(pts, scalars='turn', cmap='plasma', point_size=3, render_points_as_spheres=True, show_scalar_bar=False, reset_camera=False)
    pl.add_mesh(pv.PolyData(Rl, lines=np.r_[len(Rl), np.arange(len(Rl))]), color='orange', line_width=6, reset_camera=False)   # cut-away part shows as orange under the red
    pl.add_mesh(pv.PolyData(Kl, lines=np.r_[len(Kl), np.arange(len(Kl))]), color='red', line_width=6, reset_camera=False); pl.add_mesh(pv.PolyData(Kl[[0, -1]]), color='black', point_size=14, render_points_as_spheres=True, reset_camera=False)
    vv = m.face_normals[near].mean(0); vv /= np.linalg.norm(vv); pl.camera.focal_point = focus.tolist(); pl.camera.position = (focus+vv*30).tolist(); pl.camera.up = (0, 0, 1); pl.camera.parallel_projection = True; pl.camera.parallel_scale = 0.9
    pl.screenshot(f'out/cut_{name}.png')
print('saved out/cut_A.png out/cut_B.png (orange = cut away, red = kept, black = new ends)')
