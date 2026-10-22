import numpy as np, sys, os, time, trimesh
sys.path.insert(0, r'D:\Stuff\PrintSuit\SegmentationTool'); os.chdir(r'D:\Stuff\PrintSuit\SegmentationTool'); import voxel, refine
from scipy.spatial import cKDTree
def main():
    m = trimesh.load('out/orcfull_patch.stl', process=False); P = np.asarray(m.vertices, float); F = np.asarray(m.faces); h = float((P.max(0)-P.min(0)).max())
    d = np.load('out/debug/20261020_144432/lines.npz'); n = d['n']; pts = d['pts']; off = np.r_[0, np.cumsum(n)]; L = pts[off[318]:off[319]]
    size = 0.087/3; box = (L.min(0)-0.3-size, L.max(0)+0.3+size)
    t0 = time.time(); r0 = voxel.run_direct(P, F, h/size, None, height=h, flat=True, box=box, workers=1); t_off = time.time()-t0
    t0 = time.time(); r1 = voxel.run_direct(P, F, h/size, None, height=h, flat=True, box=box, workers=1, blend=True); t_on = time.time()-t0
    r8 = voxel.run_direct(P, F, h/size, None, height=h, flat=True, box=box, workers=8, blend=True)
    print(f'box cubes 1 core: blend off {t_off:.1f} s, on {t_on:.1f} s')
    print('same cubes:', (r0['idx'] == r1['idx']).all(), (r1['idx'] == r8['idx']).all(), 'same faces:', (r0['face'] == r1['face']).all(), '1 vs 8 cores max N diff:', np.abs(r1['N']-r8['N']).max())
    t = np.load('out/blend_test.npz'); C1 = r1['lo']+(r1['idx']+0.5)*r1['size']; j = cKDTree(C1).query(t['C'])[1]
    ang = np.degrees(np.arccos(np.clip((r1['N'][j]*t['N2']).sum(1), -1, 1)))
    print(f'exact clipped area vs the sample-point blend of the test: angle p50 {np.median(ang):.2f} p90 {np.percentile(ang,90):.2f} max {ang.max():.2f} deg')
    ang0 = np.degrees(np.arccos(np.clip((r0['N'][j]*t['N1']).sum(1), -1, 1))); print('blend off vs the test nearest arrows, max deg:', ang0.max())
    # the real stage on line 318 against the saved result
    curves = [pts[off[i]:off[i+1]] for i in range(len(n))]
    t0 = time.time(); res = refine.refine_all(curves, P, F, 0.087, h, log=print, only=318); print(f'refine_all line 318: {time.time()-t0:.1f} s')
    old = np.load('out/app_refine_test.npz'); on = old['n']; op = old['pts']
    new = np.vstack(res['lines']); dd = cKDTree(op).query(new)[0]
    print(f'refined line 318 with the blend: {len(res["lines"])} line(s), nodes {[len(l) for l in res["lines"]]}; vs this morning (nodes {list(on)}): median {np.median(dd):.3f} mm, max {dd.max():.3f} mm')
if __name__ == '__main__': main()
