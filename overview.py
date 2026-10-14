"""Whole-model render of trails from several sides, to pick seams to judge on a new model.
usage: python overview.py mesh chains.npz out_prefix [min_nodes]
Writes out_prefix_<view>.png for front/back/left/right, plus prints the model's bounds and the longest trails' centres."""
import sys, numpy as np, trimesh, pyvista as pv
path, npz, out = sys.argv[1:4]; min_nodes = int(sys.argv[4]) if len(sys.argv) > 4 else 20
z = np.load(npz); pts, n = z['pts'], z['n']; chains = [c for c in np.split(np.arange(len(pts)), np.cumsum(n)[:-1]) if len(c) >= min_nodes]
m = trimesh.load(path, process=True); P, F = np.asarray(m.vertices), np.asarray(m.faces)
print('bounds lo', P.min(0).round(1), 'hi', P.max(0).round(1), 'faces', len(F), 'trails shown', len(chains))
for ch in sorted(chains, key=len, reverse=True)[:12]:
    print(f'  trail {len(ch):4d} nodes, centre {pts[ch].mean(0).round(1)}, from {pts[ch[0]].round(1)} to {pts[ch[-1]].round(1)}')
mesh = pv.PolyData(P, np.hstack([np.full((len(F), 1), 3), F]).ravel()); rng = np.random.default_rng(0); cols = rng.uniform(0.1, 0.95, (len(chains), 3))
sc = np.zeros((len(pts), 3)); cl = np.concatenate([np.r_[len(ch), ch] for ch in chains])
for ci, ch in enumerate(chains): sc[ch] = cols[ci]
poly = pv.PolyData(pts, lines=cl); poly['c'] = sc
for name, cp, az in [('front', 'xz', 0), ('back', 'xz', 180), ('left', 'yz', 0), ('right', 'yz', 180)]:
    pl = pv.Plotter(off_screen=True, window_size=(1600, 1600)); pl.add_mesh(mesh, color='lightgrey', smooth_shading=True)
    pl.add_mesh(poly, scalars='c', rgb=True, line_width=3, render_lines_as_tubes=True)
    pl.camera_position = cp
    if az: pl.camera.Azimuth(az)
    pl.add_text(f'{npz} {name}', font_size=10); pl.screenshot(out+f'_{name}.png'); print('saved', out+f'_{name}.png')
