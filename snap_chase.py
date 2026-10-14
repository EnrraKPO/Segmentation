"""Snap refined chains onto the crease using a FINE grid laid only along the trails.

The coarse grid finds the trails; its dots carry up to half a cell of grid jitter, so the chains zigzag. Here the mesh
faces within `band` mm of any node are gridded again at `fine` cells per model height (same turn measure as
voxel.py), and each node moves to the fine crest cell (turn peaks across the line) with the strongest turn within
`reach` mm of it, looking only sideways to the chain (not along it). Nodes with no crest cell in reach stay put.
usage: python snap_chase.py mesh refined.npz chase_nodes.npz out_prefix [fine sens band reach]"""
import sys, os, time, numpy as np, trimesh, pyvista as pv
from scipy.spatial import cKDTree
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import voxel

path, ref_npz, nodes_npz, out = sys.argv[1:5]; a = sys.argv[5:]
fine, sens, band, reach = (int(a[0]) if a else 2400), (float(a[1]) if len(a) > 1 else 6.0), (float(a[2]) if len(a) > 2 else 0.6), (float(a[3]) if len(a) > 3 else 0.12)

z = np.load(ref_npz); pts, n = z['pts'], z['n']; chains = np.split(np.arange(len(pts)), np.cumsum(n)[:-1])
zn = np.load(nodes_npz); X, lo, size, focus = zn['X'], zn['lo'], zn['size'], zn['focus']
m = trimesh.load(path, process=True); P, F = np.asarray(m.vertices), np.asarray(m.faces); height = (P.max(0)-P.min(0)).max()
print(f'{len(chains)} chains, {len(pts)} nodes; coarse cell {size:.4f} mm, fine cell {height/fine:.4f} mm')

# band of faces along the trails
t0 = time.time(); cen = P[F].mean(1); near = cKDTree(pts).query(cen, distance_upper_bound=band)[0] < band
Fb = F[near]; vid, Fb = np.unique(Fb, return_inverse=True); Pb = P[vid]; Fb = Fb.reshape(-1, 3)
print(f'band: {near.sum()} of {len(F)} faces ({100*near.mean():.1f}%)')
cache = out+f'_fine{fine}.npz'
if os.path.exists(cache):
    zc = np.load(cache); G, GT, Gline = zc['G'], zc['GT'], zc['line']; print('fine grid from cache')
else:
    res = voxel.run(Pb, Fb, fine, sens, print, height=height); print(f'fine grid {time.time()-t0:.0f} s')
    G = res['lo']+(res['idx']+0.5)*res['size']; GT = res['turn']; Gline = res['line']; np.savez(cache, G=G, GT=GT, line=Gline)
C = G[Gline]; gtree = cKDTree(G); print(f'{len(G)} fine cells, {Gline.sum()} crest')

# tangent per node from its chain neighbours
T = np.zeros_like(pts)
for ch in chains:
    if len(ch) < 2: continue
    Q = pts[ch]; T[ch[1:-1]] = Q[2:]-Q[:-2]; T[ch[0]] = Q[1]-Q[0]; T[ch[-1]] = Q[-1]-Q[-2]
T /= np.linalg.norm(T, axis=1, keepdims=True)+1e-12

# each node moves sideways to the turn-weighted centre of the fine cells in a strip across the chain
new = pts.copy(); moved = 0; shift = []
for i in range(len(pts)):
    cand = gtree.query_ball_point(pts[i], reach)
    if not cand: continue
    cand = np.array(cand); off = G[cand]-pts[i]; along = off @ T[i]
    ok = np.abs(along) <= 0.4*reach                                      # strip across the chain, not along it
    w = GT[cand]**2*ok
    if w.sum() <= 0: continue
    d = (off*w[:, None]).sum(0)/w.sum(); d -= (d @ T[i])*T[i]           # sideways move only
    new[i] = pts[i]+d; moved += 1; shift.append(np.linalg.norm(d))
shift = np.array(shift)
print(f'snapped {moved} of {len(pts)} nodes; shift mm p50 {np.median(shift):.4f} p90 {np.percentile(shift, 90):.4f} max {shift.max():.4f}')

np.savez(out+'.npz', n=n, pts=new)
mesh = pv.PolyData(P, np.hstack([np.full((len(F), 1), 3), F]).ravel()); rng = np.random.default_rng(0); cols = rng.uniform(0.1, 0.95, (len(chains), 3))
L = list(n); sc = np.vstack([np.repeat(cols[ci][None], len(ch), 0) for ci, ch in enumerate(chains)])
cl = np.concatenate([np.r_[len(ch), ch] for ch in chains])
for zoom in (6, 12, 25):
    pl = pv.Plotter(off_screen=True, window_size=(1500, 1500)); pl.add_mesh(mesh, color='lightgrey', smooth_shading=True)
    pl.add_mesh(pv.PolyData(C), color='yellow', point_size=3, render_points_as_spheres=True)
    pl.add_mesh(pv.PolyData(lo+X*size), color='red', point_size=4, render_points_as_spheres=True)
    poly = pv.PolyData(new, lines=cl); poly['c'] = sc; pl.add_mesh(poly, scalars='c', rgb=True, line_width=4)
    pl.add_mesh(pv.PolyData(new), scalars=sc, rgb=True, point_size=9, render_points_as_spheres=True)
    pl.add_text(f'snapped (fine {fine}), zoom {zoom}: coarse dots red, fine crest yellow, chains coloured', font_size=10)
    pl.camera_position = 'zy'; pl.camera.focal_point = focus.tolist(); pl.camera.position = (focus+(np.array(pl.camera.position)-np.array(pl.camera.focal_point))).tolist(); pl.camera.zoom(zoom)
    pl.screenshot(out+f'_zoom{zoom}.png'); print('saved', out+f'_zoom{zoom}.png')
