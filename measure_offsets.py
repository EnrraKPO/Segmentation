"""How far are the trails from the dark line the render shows? Renders the shaded view (same camera as view_chains),
maps the darkest pixel per row back to surface y,z through coordinate images, and reports each trail's sideways
offset there. usage: python measure_offsets.py mesh chains.npz [lean.npz]"""
import sys, numpy as np, pyvista as pv, trimesh
from scipy.spatial import cKDTree
M, npz = sys.argv[1], sys.argv[2]; leanf = sys.argv[3] if len(sys.argv) > 3 else None
m = trimesh.load(M, process=True); P, F = np.asarray(m.vertices), np.asarray(m.faces)
mesh = pv.PolyData(P, np.hstack([np.full((len(F), 1), 3), F]).ravel()); focus = np.array([12, 135, 12.]); zoom = 25; W = 1500
def cam(pl):
    pl.camera_position = 'zy'; pl.camera.focal_point = focus.tolist(); pl.camera.position = (focus+(np.array(pl.camera.position)-np.array(pl.camera.focal_point))).tolist(); pl.camera.zoom(zoom)
pl = pv.Plotter(off_screen=True, window_size=(W, W)); pl.add_mesh(mesh, color='lightgrey', smooth_shading=True); cam(pl); shade = pl.screenshot(return_img=True)[:, :, 0].astype(float)
coord = {}
for ax, name in ((1, 'y'), (2, 'z')):
    pl2 = pv.Plotter(off_screen=True, window_size=(W, W)); pl2.set_background('white'); mesh[name] = P[:, ax]; lo_, hi_ = focus[ax]-2.5, focus[ax]+2.5
    pl2.add_mesh(mesh, scalars=name, cmap='gray', clim=[lo_, hi_], lighting=False, show_scalar_bar=False); cam(pl2)
    coord[name] = lo_+(hi_-lo_)*pl2.screenshot(return_img=True)[:, :, 0].astype(float)/255
zc = np.load(npz); pts, n = zc['pts'], zc['n']; chains = np.split(np.arange(len(pts)), np.cumsum(n)[:-1])
near = [c for c in chains if np.linalg.norm(pts[c]-focus, axis=1).min() < 4 and len(c) > 40]
if leanf:
    r = np.load('out/trace_grid600.npz'); Xmm = r['lo']+(r['idx']+0.5)*r['size']; lean = np.load(leanf)['lean']; t = cKDTree(Xmm)
offs = {i: [] for i in range(len(near))}
print('row | darkest y z' + (' lean' if leanf else '') + ' | per trail: offset (um)' + (', lean at trail' if leanf else ''))
for row in range(200, 1400, 150):
    line = shade[row, 300:1300]; c = 300+np.argmin(line); y, z = coord['y'][row, c], coord['z'][row, c]; s = f'{row:4d} | {y:7.3f} {z:6.3f}'
    if leanf:
        cand = np.nonzero((np.abs(Xmm[:, 1]-y) < 0.08) & (np.abs(Xmm[:, 2]-z) < 0.08))[0]; s += f' {lean[cand].max() if len(cand) else float("nan"):5.1f}'
    s += ' |'
    for i, ch in enumerate(near):
        Q = pts[ch]; k = np.argmin(np.abs(Q[:, 1]-y)); d = (Q[k, 2]-z)*1000; offs[i].append(d); s += f'  {d:+6.0f}' + (f' ({lean[t.query(Q[k])[1]]:4.1f})' if leanf else '')
    print(s)
for i, ch in enumerate(near): print(f'trail {i} ({len(ch)} nodes): median |offset| {np.median(np.abs(offs[i])):.0f} um')
