import sys, os, time, numpy as np, trimesh
sys.path.insert(0, r'D:\Stuff\PrintSuit\SegmentationTool'); os.chdir(r'D:\Stuff\PrintSuit\SegmentationTool'); import voxel
def main():
    m = trimesh.load('out/orcfull_patch.stl', process=False); P = np.asarray(m.vertices, float); F = np.asarray(m.faces); h = float((P.max(0)-P.min(0)).max())
    t0 = time.time(); log = lambda s: print(f'{time.time()-t0:6.1f} s  {s}', flush=True)
    voxel.run_direct(P, F, int(round(h/0.087)), log, height=h, flat=True, blend=True)
if __name__ == '__main__': main()
