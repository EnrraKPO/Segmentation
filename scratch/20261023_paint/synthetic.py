"""regions.py on a sphere: one closed circle -> 2 regions, two crossing circles -> 4, an open arc -> 1."""
import sys, os, time; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import numpy as np, trimesh, regions
m = trimesh.creation.icosphere(subdivisions=5, radius=10.0); P = np.asarray(m.vertices); F = np.asarray(m.faces)
print('sphere', len(F), 'faces, edge median', np.median(m.edges_unique_length))
def circle(n, z, tilt=0.0, closed=True, frac=1.0):
    th = np.linspace(0, 2*np.pi*frac, n, endpoint=False) if closed else np.linspace(0, 2*np.pi*frac, n)
    r = np.sqrt(100-z*z); pts = np.stack([r*np.cos(th), r*np.sin(th), np.full(n, z)], 1)
    c, s = np.cos(tilt), np.sin(tilt); R = np.array([[1, 0, 0], [0, c, -s], [0, s, c]]); pts = pts @ R.T
    pts *= 1.004                                            # a hair off the surface, as walk nodes are
    return np.vstack([pts, pts[:1]]) if closed else pts
log = print
for name, lines, want in [('one circle', [circle(1500, 3.0)], 2), ('two crossing circles', [circle(1500, 3.0), circle(1500, -2.0, tilt=0.9)], 4), ('open arc', [circle(800, 3.0, closed=False, frac=0.6)], 1), ('circle + arc into it', [circle(1500, 3.0), circle(600, -2.0, tilt=0.9, closed=False, frac=0.3)], 2)]:
    t0 = time.time(); r = regions.run(P, F, lines, log=log, workers=1); con = r['con']
    tot = con['areas'].sum(); A = m.area
    print(f'{name}: {con["nreg"]} regions (want {want}), fails {r["emb"]["fails"]}, areas {np.round(np.sort(con["areas"])[::-1][:6], 2)} sum {tot:.3f} vs mesh {A:.3f}, colours {r["colour"][:6]}, {time.time()-t0:.1f} s')
    assert con['nreg'] == want, name
    assert abs(tot-A) < 1e-6*A
print('all good')
