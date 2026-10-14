"""Trace lines by STEPPING along the surface's own line direction (not from dot to dot).

The grid (voxel.py) gives every surface cell its turn (how fast the surface bends) and d1 (the direction across
which it bends). The line direction at a cell is at right angles to d1 on the surface. A trail starts at the
strongest unused line cell and marches: step one cell along the current direction, then slide the new point
SIDEWAYS to the turn-weighted balance point of the surface cells in a strip across the line (both rims of a groove
balance to its bottom), then read the direction there from the cells' own line directions. Nodes are thus generated
by the direction; the dots only confirm the line goes on. A trail stops where the bending in the strip drops below
`sens` (line ends), where the direction turns more than `agree` degrees in one step (corner), or where it runs
into cells already claimed by another trail (junction / same line again). Traced both ways from each seed.
usage: python trace.py mesh out_prefix [cells sens agree min_len]"""
import sys, os, time, numpy as np, trimesh
from scipy.spatial import cKDTree
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import voxel

path, out = sys.argv[1], sys.argv[2]; a = sys.argv[3:]
cells, sens, agree, min_len = (int(a[0]) if a else 600), (float(a[1]) if len(a) > 1 else 6.0), (float(a[2]) if len(a) > 2 else 25.0), (int(a[3]) if len(a) > 3 else 10)
STEP = 1.0      # cells: advance per node
STRIP = float(os.environ.get('STRIP', 2.0))     # cells: half-width of the strip across the line used to balance the node and read the direction
ALONG = float(os.environ.get('ALONG', 0.6))     # cells: half-length of the strip along the line
SIGMA = float(os.environ.get('SIGMA', 1.0))     # cells: gaussian falloff of a cell's weight with sideways distance (no hard strip edge)
MOM = float(os.environ.get('MOM', 0.7))       # heading momentum: new heading = MOM*old + (1-MOM)*read (0 = read fresh each step)
CLAIM = float(os.environ.get('CLAIM', 1.2))     # cells: line cells this close to a node belong to its trail

cache = out+f'_grid{cells}.npz'
if os.path.exists(cache):
    r = dict(np.load(cache)); print('grid from cache')
else:
    m = trimesh.load(path, process=True); P, F = np.asarray(m.vertices), np.asarray(m.faces)
    t0 = time.time(); res = voxel.run(P, F, cells, sens, print); print(f'grid {time.time()-t0:.0f} s')
    r = dict(idx=res['idx'], N=res['N'], turn=res['turn'], line=res['line'], d1=res['d1'], lo=res['lo'], size=res['size']); np.savez(cache, **r)
X = r['idx']+0.5; N = r['N']; turn = r['turn']; line = r['line']; lo, size = r['lo'], r['size']
T = np.cross(N, r['d1']); T /= np.linalg.norm(T, axis=1, keepdims=True)+1e-12     # line direction per cell (unsigned)
W = turn**2; tree = cKDTree(X); li = np.nonzero(line)[0]; ltree = cKDTree(X[li])
owner = np.full(len(li), -1)                                                      # trail id per line cell
cos_agree = np.cos(np.radians(agree))
print(f'{len(X)} surface cells, {len(li)} line cells')


def local(p, h):
    """at point p heading h: balanced point across the line, direction there, strongest turn in the strip"""
    nb = np.array(tree.query_ball_point(p, max(STRIP, ALONG)))
    if len(nb) == 0: return None
    off = X[nb]-p; al = off @ h; keep = np.abs(al) <= ALONG
    nb, off, al = nb[keep], off[keep], al[keep]
    if len(nb) == 0: return None
    side = off-al[:, None]*h; w = W[nb]*np.exp(-0.5*(np.linalg.norm(side, axis=1)/SIGMA)**2)*np.exp(-0.5*(al/ALONG)**2)
    if w.sum() <= 0: return None
    d = (off*w[:, None]).sum(0)/w.sum(); d -= (d @ h)*h                           # sideways only
    t = T[nb]*np.sign(T[nb] @ h)[:, None]; hd = (t*w[:, None]).sum(0); hd /= np.linalg.norm(hd)+1e-12
    return p+d, hd, turn[nb].max()


def march(p, h, tid):
    """from node p heading h, generate nodes until the line ends; returns list of positions"""
    outp = []
    while True:
        q = p+STEP*h; loc = local(q, h)
        if loc is None: break
        q, hq, tmax = loc; hq = MOM*h+(1-MOM)*hq; hq /= np.linalg.norm(hq)+1e-12
        if tmax < sens: break                                                     # no bending here: line ended
        if hq @ h < cos_agree: break                                              # sharp corner
        near = ltree.query_ball_point(q, CLAIM)
        if any(owner[k] not in (-1, tid) for k in near): break                   # another trail owns this
        back = int(2*CLAIM)+5
        if len(outp) > back and np.linalg.norm(np.array(outp[:-back])-q, axis=1).min() < CLAIM: break   # ran into itself (closed loop); ignore the last 2*CLAIM nodes, which are always within CLAIM
        for k in near: owner[k] = tid
        outp.append(q); p, h = q, hq
        if len(outp) > 20000: break
    return outp


chains = []; order = li[np.argsort(-turn[li])]
for c in order:
    k = np.searchsorted(li, c)
    if owner[k] != -1: continue
    tid = len(chains); p0 = X[c]; h0 = T[c]
    loc = local(p0, h0)
    if loc is None: continue
    p0, h0, _ = loc
    for k2 in ltree.query_ball_point(p0, CLAIM): owner[k2] = tid
    fw = march(p0, h0, tid); bw = march(p0, -h0, tid)
    ch = bw[::-1]+[p0]+fw
    if len(ch) < min_len:
        owner[owner == tid] = -1; chains.append(None); continue                   # keep ids stable; drop trail
    chains.append(np.array(ch))
chains = [c for c in chains if c is not None]
L = [len(c) for c in chains]
print(f'trails {len(chains)}, nodes {sum(L)}, longest {max(L)}, >=20 nodes: {sum(l >= 20 for l in L)}, line cells unclaimed {(owner == -1).sum()} of {len(li)}')
np.savez(out+'.npz', n=np.array(L), pts=np.vstack([lo+c*size for c in chains]))
print('saved', out+'.npz')
