"""Whole pipeline on a patch: cut a box (+-R mm) around CENTRE, find dots at the given cube size (mm), build strokes, save.
usage: python patch.py mesh out_prefix cube_mm x,y,z R [cut]   env: W L K ALIGN AGREE MINLEN (stroke params, mm)
Writes out_prefix_grid.npz (dots) and out_prefix.npz (strokes)."""
import sys, os, time, numpy as np, trimesh
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import creases


def main():
    path, out, cube = sys.argv[1], sys.argv[2], float(sys.argv[3]); c = np.array([float(v) for v in sys.argv[4].split(',')]); R = float(sys.argv[5]); cut = float(sys.argv[6]) if len(sys.argv) > 6 else 17.0
    m = trimesh.load(path, process=True); P, F = np.asarray(m.vertices), np.asarray(m.faces); height = float(m.extents.max())
    keep = (np.abs(P[F]-c).max(axis=(1, 2)) < R); Fp = F[keep]; used = np.unique(Fp); remap = -np.ones(len(P), int); remap[used] = np.arange(len(used)); Pp, Fp = P[used], remap[Fp]
    cells = int(round(height/cube)); print(f'patch {len(Fp)} faces, cube {height/cells:.4f} mm ({cells} cells)')
    t0 = time.time(); res = creases.find_dots(Pp, Fp, cells, cut, 0.0, True, print, height=height); print(f'dots {time.time()-t0:.0f} s')
    np.savez(out+'_grid.npz', **res)
    e = os.environ; kw = {k: float(e[K]) for k, K in [('W', 'W'), ('L', 'L'), ('K', 'K'), ('align', 'ALIGN'), ('agree', 'AGREE')] if K in e}
    if 'MINLEN' in e: kw['minlen'] = int(e['MINLEN'])
    t0 = time.time(); cv = creases.join(creases.strokes(res, log=print, **kw), log=print); print(f'strokes {time.time()-t0:.0f} s')
    np.savez(out+'.npz', n=np.array([len(c) for c in cv]), pts=np.vstack(cv) if cv else np.zeros((0, 3))); print('saved', out+'.npz')


if __name__ == '__main__': main()
