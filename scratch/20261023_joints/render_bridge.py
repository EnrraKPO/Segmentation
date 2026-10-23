"""render out/bridge_<dump>.npz: measured nodes dark red, bridges magenta, loose ends black, shared nodes green, walk lines blue.
usage: python render_bridge.py dumpdir [x,y,z half] ...  (each extra arg = one close-up; no arg = the dump's own view)"""
import sys, os, json, numpy as np, pyvista as pv, trimesh
os.chdir(r'D:\Stuff\PrintSuit\SegmentationTool'); d = sys.argv[1]; tag = os.path.basename(d)+os.environ.get('TAG', '')
z = np.load(f'out/bridge_{tag}.npz'); n, pts, meas = z['n'], z['pts'], z['measured']; off = np.cumsum(np.r_[0, n[:-1]])
w = np.load(d+'/lines.npz'); prm = json.load(open(d+'/settings.json'))
m = trimesh.load('out/orcfull_patch.stl', process=False); P, F = np.asarray(m.vertices), np.asarray(m.faces); mesh = pv.PolyData(P, np.hstack([np.full((len(F), 1), 3), F]).ravel())
def snap(p): C, _, f = trimesh.proximity.closest_point(m, p); return C+0.01*m.face_normals[f]
S = snap(pts); Wp = snap(w['pts'])
def polylines(sel_runs):
    cells = np.concatenate([np.r_[len(rn), rn] for rn in sel_runs]) if sel_runs else None; return cells
mruns, bruns = [], []
for a, k in zip(off, n):
    f = meas[a:a+k]; idx = np.arange(a, a+k); cut = np.flatnonzero(np.diff(f.astype(int)) != 0)+1; runs = np.split(idx, cut)
    for rn in runs:
        if f[rn[0]-a]: mruns.append(rn)
        else: bruns.append(np.r_[max(rn[0]-1, a), rn, min(rn[-1]+1, a+k-1)])
views = [('view', None)]+[(f'zoom{i}', (np.array([float(v) for v in s.split(',')[:3]]), float(s.split(',')[3]))) for i, s in enumerate(sys.argv[2:])]
for name, zoom in views:
    pl = pv.Plotter(off_screen=True, window_size=(1600, 1300)); pl.set_background('white'); pl.add_mesh(mesh, color='lightgrey', smooth_shading=True)
    pl.add_mesh(pv.PolyData(Wp, lines=np.concatenate([np.r_[k, np.arange(a, a+k)] for a, k in zip(np.cumsum(np.r_[0, w['n'][:-1]]), w['n'])])), color='blue', line_width=2)
    pl.add_mesh(pv.PolyData(S, lines=polylines(mruns)), color='darkred', line_width=4)
    if bruns: pl.add_mesh(pv.PolyData(S, lines=polylines(bruns)), color='magenta', line_width=5)
    if len(z['loose']): pl.add_mesh(pv.PolyData(snap(z['loose'])), color='black', point_size=9, render_points_as_spheres=True)
    if len(z['shared']): pl.add_mesh(pv.PolyData(snap(z['shared'])), color='lime', point_size=9, render_points_as_spheres=True)
    cam = prm['camera']; pl.camera_position = cam
    if zoom is not None:
        c, half = zoom; pos = np.array(cam[0]); foc = np.array(cam[1]); dirn = (pos-foc)/np.linalg.norm(pos-foc); pl.camera_position = [c+dirn*20, c, cam[2]]; pl.camera.parallel_projection = True; pl.camera.parallel_scale = half
    pl.add_text('walk blue, measured dark red, bridges magenta, loose ends black, shared nodes green', font_size=10); pl.screenshot(f'out/bridge_{tag}_{name}.png'); pl.close(); print('out/bridge_'+tag+'_'+name+'.png')
