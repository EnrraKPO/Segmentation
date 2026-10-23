"""stage 2b on the dump's walk lines: counts, timing, out/merge_<dump>.npz, and renders before/after at a spot"""
import sys, os, time, json, numpy as np, pyvista as pv, trimesh
sys.path.insert(0, r'D:\Stuff\PrintSuit\SegmentationTool'); os.chdir(r'D:\Stuff\PrintSuit\SegmentationTool'); import merge
def main():
    d = sys.argv[1]; margin = float(sys.argv[2]) if len(sys.argv) > 2 else 0.07; prm = json.load(open(d+'/settings.json'))
    w = np.load(d+'/lines.npz'); W = [w['pts'][a:a+k] for a, k in zip(np.cumsum(np.r_[0, w['n'][:-1]]), w['n'])]
    t0 = time.time(); res = merge.run(W, margin, log=print); print(f'{time.time()-t0:.1f} s', res['report'])
    M = res['lines']; np.savez(f'out/merge_{os.path.basename(d)}.npz', n=np.array([len(l) for l in M]), pts=np.vstack(M))
    # how close do lines still run: for each node, distance to the nearest node of another line, before and after
    from scipy.spatial import cKDTree
    def closeness(S):
        pts = np.vstack(S); o = np.repeat(np.arange(len(S)), [len(c) for c in S]); t = cKDTree(pts); dd, j = t.query(pts, k=40); out = np.full(len(pts), np.inf)
        for q in range(40):
            m = (o[j[:, q]] != o) & np.isinf(out); out[m] = dd[m, q]
        return f'nodes with another line within {margin} mm: {(out <= margin).mean()*100:.1f} %, within 0.2 mm: {(out <= 0.2).mean()*100:.1f} %'
    print('before:', closeness(W)); print('after: ', closeness(M))
    m = trimesh.load('out/orcfull_patch.stl', process=False); P, F = np.asarray(m.vertices), np.asarray(m.faces); mesh = pv.PolyData(P, np.hstack([np.full((len(F), 1), 3), F]).ravel())
    def snap(p): C, _, f = trimesh.proximity.closest_point(m, p); return C+0.01*m.face_normals[f]
    cam = prm['camera']; pos = np.array(cam[0]); foc = np.array(cam[1]); dirn = (pos-foc)/np.linalg.norm(pos-foc)
    for spot in sys.argv[3:]:
        c = np.array([float(v) for v in spot.split(',')[:3]]); half = float(spot.split(',')[3])
        for name, S, col in (('before', W, 'blue'), ('after', M, 'darkgreen')):
            pl = pv.Plotter(off_screen=True, window_size=(1400, 1200)); pl.set_background('white'); pl.add_mesh(mesh, color='lightgrey', smooth_shading=True)
            for L in S:
                if (np.abs(L-c) < half*1.5).all(1).any(): p = snap(L); pl.add_mesh(pv.PolyData(p, lines=np.r_[len(p), np.arange(len(p))]), color=col, line_width=3); pl.add_mesh(pv.PolyData(p[[0, -1]]), color='black', point_size=8, render_points_as_spheres=True)
            if name == 'after' and res['shared']:
                sh = np.array([s[2] for s in res['shared']]+[s[3] for s in res['shared']]); sh = sh[(np.abs(sh-c) < half*1.5).all(1)]
                if len(sh): pl.add_mesh(pv.PolyData(snap(sh)), color='lime', point_size=10, render_points_as_spheres=True)
            pl.camera_position = [c+dirn*20, c, cam[2]]; pl.camera.parallel_projection = True; pl.camera.parallel_scale = half
            pl.add_text(f'{name} merge (margin {margin} mm): lines {col}, ends black' + (', shared nodes green' if name == 'after' else ''), font_size=10); pl.screenshot(f'out/merge_{name}_{spot.replace(",", "_")}.png'); pl.close()
        print('rendered', spot)
if __name__ == '__main__': main()
