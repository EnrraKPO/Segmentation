"""Fine cubes over a patch only: cut the mesh to a box around CENTRE (+-R mm), grid it at `cells` per FULL model
height (so cell size = height/cells, comparable across runs), save the grid cache trace2.py expects.
usage: python patch_grid.py mesh out_prefix cells x,y,z R"""
import sys, os, time, numpy as np, trimesh
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import voxel
path, out, cells = sys.argv[1], sys.argv[2], int(sys.argv[3]); c = np.array([float(v) for v in sys.argv[4].split(',')]); R = float(sys.argv[5])
m = trimesh.load(path, process=True); P, F = np.asarray(m.vertices), np.asarray(m.faces); height = m.extents.max()
keep = (np.abs(P[F]-c).max(axis=(1, 2)) < R); Fp = F[keep]; used = np.unique(Fp); remap = -np.ones(len(P), int); remap[used] = np.arange(len(used))
Pp, Fp = P[used], remap[Fp]; trimesh.Trimesh(Pp, Fp, process=False).export(out+'_patch.stl')
print(f'patch: {len(Fp)} faces of {len(F)}, cell {height/cells:.3f} mm')
t0 = time.time(); spp = float(os.environ.get('SPP', 2)); flat = os.environ.get('FLAT', '0') == '1'; r = voxel.run(Pp, Fp, cells, 6.0, print, height=height, spp=spp, flat=flat); print('spp', spp, 'flat', flat); print(f'grid {time.time()-t0:.0f} s')
np.savez(out+f'_grid{cells}.npz', idx=r['idx'], N=r['N'], turn=r['turn'], line=r['line'], d1=r['d1'], lo=r['lo'], size=r['size'])
print('saved', out+f'_grid{cells}.npz')
