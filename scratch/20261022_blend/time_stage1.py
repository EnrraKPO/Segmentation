"""time the cube stage alone on the whole patch, blend off and on, and check the dots are the same as stage1_blend.py's"""
import sys, os, time, numpy as np, trimesh
sys.path.insert(0, r'D:\Stuff\PrintSuit\SegmentationTool'); os.chdir(r'D:\Stuff\PrintSuit\SegmentationTool'); import voxel

def main():
    m = trimesh.load('out/orcfull_patch.stl', process=False); P = np.asarray(m.vertices, float); F = np.asarray(m.faces); h = float((P.max(0)-P.min(0)).max())
    cells = max(1, int(round(h/0.087))); ref = np.load('out/blend1_dots.npz')
    for name, blend in (('off', False), ('on', True)):
        t0 = time.time(); r = voxel.run_direct(P, F, cells, None, height=h, flat=True, blend=blend); t = time.time()-t0
        same = (r['idx'] == ref['idx']).all() and np.abs(r['N']-ref['N1' if blend else 'N0']).max() < 1e-9
        print(f'cube stage blend {name}: {t:.0f} s on {os.cpu_count()} cores; same cubes and arrows as the saved run: {same}')
if __name__ == '__main__': main()
