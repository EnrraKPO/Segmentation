"""Trace the lines the light sees: step along the surface where its normal LEANS away from its surroundings.

Field (per surface cell of the voxel grid): lean = angle between the cell's normal and the mean normal of all surface
cells within `radius` cells. Flat or gently curved surface: ~0. The flank of a groove or cord: large, peaking in the
middle of the flank, which is where a shader darkens, for any light. The lean direction (the normal's sideways
component) points across the line, so the line direction is at right angles to it on the surface: a geometric truth
per cell, no neighbour statistics.

Tracer: from the strongest unused cell, step one cell along the line direction; balance the new node sideways to the
lean-weighted centre of the cells in a strip (`along` cells long, `strip` cells wide, gaussian weights); read the
direction there from the cells' own line directions, with momentum. Stop when the lean in the strip falls under
`thresh` degrees (line ends), on a sharp turn (over `agree` degrees), on running into another trail or itself.
usage: python trace_lean.py mesh out_prefix [cells radius thresh agree min_len]"""
import sys, os, time, numpy as np, trimesh
from scipy.spatial import cKDTree
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import voxel

path, out = sys.argv[1], sys.argv[2]; a = sys.argv[3:]
cells, radius, thresh, agree, min_len = (int(a[0]) if a else 600), (float(a[1]) if len(a) > 1 else 4.0), (float(a[2]) if len(a) > 2 else 8.0), (float(a[3]) if len(a) > 3 else 25.0), (int(a[4]) if len(a) > 4 else 10)
STEP = 1.0; STRIP = float(os.environ.get('STRIP', 3.5)); ALONG = float(os.environ.get('ALONG', 3.0)); SIGMA = 1.0; MOM = 0.7; CLAIM = float(os.environ.get('CLAIM', 2.0))
stops = dict(faded=0, corner=0, other=0, self=0, none=0)

# grid (cached)
cache = out+f'_grid{cells}.npz'
if os.path.exists(cache):
    r = dict(np.load(cache)); print('grid from cache')
else:
    m = trimesh.load(path, process=True); P, F = np.asarray(m.vertices), np.asarray(m.faces)
    t0 = time.time(); res = voxel.run(P, F, cells, 6.0, print); print(f'grid {time.time()-t0:.0f} s')
    r = dict(idx=res['idx'], N=res['N'], lo=res['lo'], size=res['size']); np.savez(cache, **r)
X = (r['idx']+0.5).astype(float); N = r['N']; lo, size = r['lo'], r['size']; n = len(X); tree = cKDTree(X)

# lean field (cached): angle between a cell's normal and the SMOOTH BODY normal there. The smooth body normal is the
# normal field blurred over `radius` mm (gaussian, on a coarse 3D grid), so a seam, cord or fold of any width leans
# against the body it sits on, and its crest is the fold itself. With SMOOTH unset the old local mean (radius in cells) is used.
smooth = os.environ.get('SMOOTH')          # mm: gaussian sigma of the body blur
lcache = out+(f'_lean{cells}_s{float(smooth):g}.npz' if smooth else f'_lean{cells}_r{radius:g}.npz')
if os.path.exists(lcache):
    z = np.load(lcache); lean, T = z['lean'], z['T']; print('lean from cache')
else:
    t0 = time.time(); lean = np.zeros(n); T = np.zeros((n, 3))
    if smooth:
        from scipy.ndimage import gaussian_filter, map_coordinates
        sig = float(smooth); h = sig/3.0                                         # coarse voxel = sigma/3 mm
        Xmm = lo+X*size; g = ((Xmm-lo)/h).astype(int); dims = g.max(0)+3
        cnt = np.zeros(dims); vec = np.zeros((3,)+tuple(dims))
        np.add.at(cnt, tuple(g.T), 1.0)
        for a in range(3): np.add.at(vec[a], tuple(g.T), N[:, a])
        cnt = gaussian_filter(cnt, 3.0); vec = np.stack([gaussian_filter(v, 3.0) for v in vec])
        crd = ((Xmm-lo)/h).T
        M = np.stack([map_coordinates(v, crd, order=1) for v in vec], 1)/(map_coordinates(cnt, crd, order=1)[:, None]+1e-9)
        M /= np.linalg.norm(M, axis=1, keepdims=True)+1e-12
        c = np.clip((N*M).sum(1), -1, 1); lean[:] = np.degrees(np.arccos(c)); tilt = N-c[:, None]*M
        t = np.cross(N, tilt); T[:] = t/(np.linalg.norm(t, axis=1, keepdims=True)+1e-12)
    else:
        K = int(10*radius**2)+40
        for s in range(0, n, 20000):
            e = min(n, s+20000); d, j = tree.query(X[s:e], k=K, distance_upper_bound=radius)
            ok = np.isfinite(d); jj = np.where(ok, j, 0)
            M = (N[jj]*ok[:, :, None]).sum(1); M /= np.linalg.norm(M, axis=1, keepdims=True)+1e-12
            Ns = N[s:e]; c = np.clip((Ns*M).sum(1), -1, 1); lean[s:e] = np.degrees(np.arccos(c))
            tilt = Ns-c[:, None]*M; t = np.cross(Ns, tilt); T[s:e] = t/(np.linalg.norm(t, axis=1, keepdims=True)+1e-12)
    np.savez(lcache, lean=lean, T=T); print(f'lean field {time.time()-t0:.0f} s')
print(f'{n} surface cells; lean deg p50 {np.median(lean):.1f} p90 {np.percentile(lean, 90):.1f} p99 {np.percentile(lean, 99):.1f}; cells over {thresh:g}: {(lean >= thresh).sum()}')
A = np.cross(T, N); A /= np.linalg.norm(A, axis=1, keepdims=True)+1e-12                    # across-the-line direction on the surface
def side_lean(sign):
    d, j = tree.query(X+sign*A, distance_upper_bound=1.0); return np.where(np.isfinite(d), lean[np.minimum(j, n-1)], 0.0)
crest = (lean >= thresh) & (lean >= side_lean(1.0)) & (lean > side_lean(-1.0))
print(f'crest cells {crest.sum()} of {(lean >= thresh).sum()} strong')
W = lean**2; li = np.nonzero(crest if os.environ.get('SEED', 'crest') == 'crest' else lean >= thresh)[0]; ltree = cKDTree(X[li]); owner = np.full(len(li), -1); cos_agree = np.cos(np.radians(agree))


def local(p, h):
    nb = np.array(tree.query_ball_point(p, STRIP))
    if len(nb) == 0: return None
    off = X[nb]-p; al = off @ h; side = off-al[:, None]*h
    w = W[nb]*np.exp(-0.5*(np.linalg.norm(side, axis=1)/SIGMA)**2)*np.exp(-0.5*(al/ALONG)**2)
    if w.sum() <= 0: return None
    d = (side*w[:, None]).sum(0)/w.sum()
    t = T[nb]*np.sign(T[nb] @ h)[:, None]; hd = (t*w[:, None]).sum(0); hd /= np.linalg.norm(hd)+1e-12
    return p+d, hd, lean[nb].max()


def march(p, h, tid):
    outp = []
    while True:
        loc = local(p+STEP*h, h)
        if loc is None: stops['none'] += 1; break
        q, hq, lmax = loc; hq = MOM*h+(1-MOM)*hq; hq /= np.linalg.norm(hq)+1e-12
        if lmax < thresh: stops['faded'] += 1; break
        if hq @ h < cos_agree: stops['corner'] += 1; break
        near = ltree.query_ball_point(q, CLAIM)
        if any(owner[k] not in (-1, tid) for k in near): stops['other'] += 1; break
        if len(outp) > 5 and np.linalg.norm(np.array(outp[:-5])-q, axis=1).min() < CLAIM: stops['self'] += 1; break
        for k in near: owner[k] = tid
        outp.append(q); p, h = q, hq
        if len(outp) > 20000: break
    return outp


chains = []
for c in li[np.argsort(-lean[li])]:
    k = np.searchsorted(li, c)
    if owner[k] != -1: continue
    tid = len(chains); loc = local(X[c], T[c])
    if loc is None: continue
    p0, h0, _ = loc
    for k2 in ltree.query_ball_point(p0, CLAIM): owner[k2] = tid
    ch = march(p0, -h0, tid)[::-1]+[p0]+march(p0, h0, tid)
    if len(ch) < min_len: owner[owner == tid] = -1; chains.append(None); continue
    chains.append(np.array(ch))
chains = [c for c in chains if c is not None]; L = [len(c) for c in chains]
print('stop reasons (trail ends):', stops)
print(f'trails {len(chains)}, nodes {sum(L)}, longest {max(L)}, >=20 nodes: {sum(l >= 20 for l in L)}, strong cells unclaimed {(owner == -1).sum()} of {len(li)}')
np.savez(out+'.npz', n=np.array(L), pts=np.vstack([lo+c*size for c in chains])); print('saved', out+'.npz')
