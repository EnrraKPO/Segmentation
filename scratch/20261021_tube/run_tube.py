import sys, os, time, numpy as np
sys.path.insert(0, r'D:\Stuff\PrintSuit\SegmentationTool'); os.chdir(r'D:\Stuff\PrintSuit\SegmentationTool')
import refine, trimesh

def main():
    finer = float(sys.argv[1]) if len(sys.argv) > 1 else 10; radius = float(sys.argv[2]) if len(sys.argv) > 2 else 0.3
    m = trimesh.load('out/orcfull_patch.stl', process=False); P = np.asarray(m.vertices, float); F = np.asarray(m.faces)
    h = float((P.max(0)-P.min(0)).max())
    d = np.load('out/debug/20261020_144432/lines.npz'); n = d['n']; pts = d['pts']; off = np.r_[0, np.cumsum(n)]
    curves = [pts[off[i]:off[i+1]] for i in range(len(n))]; i = refine.nearest_line(curves, [219.51, 167.43, 165.39]); print('line', i, len(curves[i]), 'nodes')
    walk = dict(step=0.087, slice=0.7, reach=0.087, pull=1.0, look=0.45, levp=2.0, minlen=8)
    t0 = time.time(); lines, t = refine.refine_line(curves[i], P, F, 0.087, h, finer=finer, radius=radius, cut=17.0, log=print, walk=walk)
    print(f'{len(lines)} fine lines, nodes {[len(l) for l in lines]}, {time.time()-t0:.1f} s')
    out = f'out/tube_{finer:g}_{radius:g}.npz'
    np.savez(out, coarse=curves[i], n=np.array([len(l) for l in lines]), pts=np.vstack(lines) if lines else np.zeros((0, 3)),
             idx=t['idx'], N=t['N'], turn=t['turn'], line=t['line'], crest=t['crest'], lo=t['lo'], size=t['size']); print('saved', out)
if __name__ == '__main__': main()
