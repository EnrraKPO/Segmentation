"""the rivet ring bridged with the parabola trajectory and with the helix: out/rivet_<model>.png"""
import sys, os, json, numpy as np, pyvista as pv, trimesh; from scipy.spatial import cKDTree
sys.path.insert(0, r'D:\Stuff\PrintSuit\SegmentationTool'); os.chdir(r'D:\Stuff\PrintSuit\SegmentationTool'); import bridge
def main():
    d = 'out/debug/20261023_114229'; prm = json.load(open(d+'/settings.json'))
    w = np.load(d+'/lines.npz'); W = [w['pts'][a:a+k] for a, k in zip(np.cumsum(np.r_[0, w['n'][:-1]]), w['n'])]
    r = np.load(d+'/refined.npz'); R = [r['pts'][a:a+k] for a, k in zip(np.cumsum(np.r_[0, r['n'][:-1]]), r['n'])]; ri = r['index']
    c47 = np.array([173.08, 143.54, 96.98])
    cam = prm['camera']; pos = np.array(cam[0]); foc = np.array(cam[1]); up = np.array(cam[2]); dirn = (pos-foc)/np.linalg.norm(pos-foc)
    right = np.cross(up, dirn); right /= np.linalg.norm(right); upv = np.cross(dirn, right)      # the rivet's place read off out/bridge_20261023_114229_zoom0.png
    c = c47+0.32*right+3.0*upv; half = 1.3; box = (np.abs(np.vstack(W)-c) < 4.0).all(1); print('rivet centre', np.round(c, 2))
    sel = sorted(set(np.repeat(np.arange(len(W)), w['n'])[box])); ks = [k for k in range(len(R)) if ri[k] in sel]
    pieces = [R[k] for k in ks]; index = [sel.index(int(ri[k])) for k in ks]; curves = [W[j] for j in sel]
    m = trimesh.load('out/orcfull_patch.stl', process=False); P, F = np.asarray(m.vertices), np.asarray(m.faces); mesh = pv.PolyData(P, np.hstack([np.full((len(F), 1), 3), F]).ravel())
    def snap(p): C, _, f = trimesh.proximity.closest_point(m, p); return C+0.01*m.face_normals[f]
    for model in ('para', 'helix'):
        res = bridge.run(pieces, index, curves, log=None, workers=1, model=model)
        pl = pv.Plotter(off_screen=True, window_size=(1400, 1200)); pl.set_background('white'); pl.add_mesh(mesh, color='lightgrey', smooth_shading=True)
        for L, mm in zip(res['lines'], res['measured']):
            S = snap(L); cut = np.flatnonzero(np.diff(mm.astype(int)) != 0)+1
            for rn in np.split(np.arange(len(L)), cut):
                if mm[rn[0]]: pl.add_mesh(pv.PolyData(S[rn], lines=np.r_[len(rn), np.arange(len(rn))]), color='darkred', line_width=5)
                else: rr = np.arange(max(rn[0]-1, 0), min(rn[-1]+2, len(L))); pl.add_mesh(pv.PolyData(S[rr], lines=np.r_[len(rr), np.arange(len(rr))]), color='magenta', line_width=6)
        for L in curves: S = snap(L); pl.add_mesh(pv.PolyData(S, lines=np.r_[len(S), np.arange(len(S))]), color='blue', line_width=2)
        if res['loose']: pl.add_mesh(pv.PolyData(snap(np.array([p for *_, p in res['loose']]))), color='black', point_size=10, render_points_as_spheres=True)
        if res['shared']: pl.add_mesh(pv.PolyData(snap(np.array([p for *_, p in res['shared']]))), color='lime', point_size=10, render_points_as_spheres=True)
        cam = prm['camera']; pos = np.array(cam[0]); foc = np.array(cam[1]); dirn = (pos-foc)/np.linalg.norm(pos-foc); pl.camera_position = [c+dirn*20, c, cam[2]]; pl.camera.parallel_projection = True; pl.camera.parallel_scale = half
        pl.add_text(f'rivet, trajectory = {"parabola over 0.5 mm" if model == "para" else "helix over the grown segment"}: walk blue, measured dark red, bridges magenta', font_size=10); pl.screenshot(f'out/rivet_{model}.png'); pl.close(); print(f'out/rivet_{model}.png', res['report'])
if __name__ == '__main__': main()
