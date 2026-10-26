import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import numpy as np, trimesh, regions
exec(open(os.path.join(os.path.dirname(__file__), 'synthetic.py')).read().split('log = print')[0])
emb = regions.embed(P, F, [circle(1500, 3.0)], log=print)
for i in emb['fid']:
    fa, fb = emb['nface'][i], emb['nface'][i+1]; a2 = emb['node2'][i]; C = emb['C2'][fa]
    shared = set(F[fa]) & set(F[fb])
    print(f'node {i}: faces {fa}->{fb}, shared verts {len(shared)}, a2 {a2}, corners {C.tolist()}, nbr {emb["nbr"][fa]}')
    # barycentric of a in fa
    A = np.c_[C[1]-C[0], C[2]-C[0]]; l = np.linalg.solve(A, a2-C[0]); print('   bary', 1-l.sum(), l)
    strokes = np.where(emb['sl'] == 0)[0]
print('--- exit search replay')
for i in emb['fid']:
    f = emb['nface'][i]; B = emb['nodes'][i+1]; pu, pv = emb['node2'][i]; C = emb['C2'][f]
    d = B-emb['v0'][f]; du = d@emb['e1'][f]-pu; dv = d@emb['e2'][f]-pv; dn = np.hypot(du, dv)
    print(f'node {i} face {f} fb {emb["nface"][i+1]} p ({pu:.6f},{pv:.6f}) d ({du:.3e},{dv:.3e}) |d| {dn:.3e}  normal part {np.cross(emb["e1"][f], emb["e2"][f])@d:.3e}')
    for k in range(3):
        k1 = (k+1)%3; cx, cy = C[k]; ex, ey = C[k1]-C[k]; den = du*ey-dv*ex; rx, ry = cx-pu, cy-pv
        s = (rx*ey-ry*ex)/den; t = (rx*dv-ry*du)/den
        print(f'   edge {k}: den {den:.3e} s {s:.4e} t {t:.6f} ttol {1e-6/np.hypot(ex,ey):.1e} smin {-1e-6/dn:.1e}')
