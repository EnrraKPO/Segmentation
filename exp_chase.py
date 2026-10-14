"""Test flight: dots -> nodes -> chains by CHASING.

Dots are information only (line cells with a tangent direction, never changed). Nodes are the mutable objects: a
node is spawned at the mean of a dot and the unused dots beside it that share its direction. A chain starts from a
seed dot and is chased in both directions: from the end node, among the unused dots that lie along its direction
(divergence under `max_div`), the one with the MOST SIMILAR direction is taken next, and only if the direction
shift is under `agree`; an abrupt shift breaks the chain. Seeds are tried strongest (most definite direction)
first. A chain with fewer than three nodes is not a trail and its dots are released.
usage: python exp_chase.py mesh out_prefix [cells sens reach agree max_div min_len]"""
import sys, os, numpy as np, trimesh, pyvista as pv
from scipy.spatial import cKDTree
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import voxel

path, out = sys.argv[1], sys.argv[2]; a = sys.argv[3:]
cells, sens, reach, agree, max_div, min_len = (int(a[0]) if a else 600), (float(a[1]) if len(a) > 1 else 6.0), (float(a[2]) if len(a) > 2 else 4.0), (float(a[3]) if len(a) > 3 else 25.0), (float(a[4]) if len(a) > 4 else 45.0), (int(a[5]) if len(a) > 5 else 3)
STEP = 2.0      # cells: how far the chase looks for the next dot (adjacent cells incl. diagonals are within 1.74)
BESIDE = 1.5    # cells: dots this close, sharing the direction and not ahead/behind, are the same node

m = trimesh.load(path, process=True); P, F = np.asarray(m.vertices), np.asarray(m.faces)
res = voxel.run(P, F, cells, sens, print)
li = np.nonzero(res['line'])[0]; X = (res['idx'][li]+0.5).astype(float); N = res['N'][li]; n = len(X)

# directions (unsigned tangents) and how definite they are, from all dots within reach
tree = cKDTree(X); pairs = tree.query_pairs(reach, output_type='ndarray'); i, j = pairs[:, 0], pairs[:, 1]
v = X[j]-X[i]; v /= np.linalg.norm(v, axis=1, keepdims=True)+1e-12
T = np.zeros((n, 3, 3))
for s_ in (i, j): np.add.at(T, s_, v[:, :, None]*v[:, None, :])
w, e = np.linalg.eigh(T); D = e[:, :, -1]; D -= (D*N).sum(1, keepdims=True)*N; D /= np.linalg.norm(D, axis=1, keepdims=True)+1e-12
aniso = (w[:, -1]-w[:, -2])/(w[:, -1]+1e-12); nn = np.bincount(i, minlength=n)+np.bincount(j, minlength=n)
print(f'{n} dots, definite direction {(aniso > 0.5).sum()}')

cos_agree, cos_div = np.cos(np.radians(agree)), np.cos(np.radians(max_div))
used = np.zeros(n, bool); node_pos, node_dir, node_dots = [], [], []


def spawn(k, d):
    """node at dot k: k plus the unused dots beside it (same direction, sideways not along); returns node index"""
    nb = np.array(tree.query_ball_point(X[k], BESIDE)); nb = nb[~used[nb]]
    off = X[nb]-X[k]; along = np.abs(off @ d); same = np.abs(D[nb] @ d) >= cos_agree
    mem = nb[(along < 0.5) & same]; mem = mem if k in mem else np.r_[mem, k]
    used[mem] = True
    dd = D[mem]*np.sign(D[mem] @ d)[:, None]; dm = dd.mean(0); dm /= np.linalg.norm(dm)+1e-12
    node_pos.append(X[mem].mean(0)); node_dir.append(dm); node_dots.append(mem); return len(node_pos)-1


def chase(node, sgn):
    """from node, step in direction sgn*node_dir; returns list of new node indices in order"""
    outn = []
    while True:
        p, d = node_pos[node], node_dir[node]
        nb = np.array(tree.query_ball_point(p, STEP)); nb = nb[~used[nb]]
        if len(nb) == 0: break
        off = X[nb]-p; dist = np.linalg.norm(off, axis=1); u = off/(dist[:, None]+1e-12)
        fwd = (u @ d)*sgn                                   # forward side of the node, along its direction
        c_dir = np.abs(D[nb] @ d)                           # direction similarity (unsigned)
        ok = (fwd >= cos_div) & (c_dir >= cos_agree)
        if not ok.any(): break
        cand = nb[ok]; best = cand[np.argmax(c_dir[ok]-1e-3*dist[ok])]   # most similar direction; nearest breaks ties
        nxt = spawn(best, node_dir[node]*sgn)
        # keep the chase sign consistent with the new node's direction
        sgn = 1.0 if node_dir[nxt] @ (node_dir[node]*sgn) >= 0 else -1.0
        outn.append(nxt); node = nxt
    return outn


chains = []; order = np.argsort(-aniso); tried = np.zeros(n, bool); freckles = 0
for k in order:
    if used[k] or tried[k] or nn[k] < 2: continue
    tried[k] = True; n0 = len(node_pos)
    s = spawn(k, D[k]); fw = chase(s, 1.0); bw = chase(s, -1.0)
    ch = bw[::-1]+[s]+fw
    if len(ch) < min_len:
        for t in ch: used[node_dots[t]] = False           # release: not a trail
        del node_pos[n0:], node_dir[n0:], node_dots[n0:]; freckles += 1; continue
    chains.append(ch)
L = [len(c) for c in chains]
print(f'chains {len(chains)}, nodes {sum(L)}, longest {max(L) if L else 0}, >=20 nodes: {sum(l >= 20 for l in L)}, freckle seeds {freckles}, dots unused {(~used).sum()} of {n}')

lo, size = res['lo'], res['size']; NP = np.array(node_pos) if node_pos else np.zeros((0, 3))
np.savez(out+'.npz', n=np.array(L), pts=np.vstack([lo+NP[c]*size for c in chains]) if chains else np.zeros((0, 3)))
np.savez(out+'_nodes.npz', NP=NP, ND=np.array(node_dir) if node_dir else np.zeros((0, 3)), n=np.array(L), chain=np.concatenate(chains) if chains else np.zeros(0, int), X=X, lo=lo, size=size, focus=np.array([12, 135, 12.]))

# pictures: waist close-up at two zooms. dots red, chains coloured, nodes as spheres on the chains
mesh = pv.PolyData(P, np.hstack([np.full((len(F), 1), 3), F]).ravel()); rng = np.random.default_rng(0)
cols = rng.uniform(0.1, 0.95, (max(len(chains), 1), 3))
for zoom in (6, 12):
    pl = pv.Plotter(off_screen=True, window_size=(1500, 1500)); pl.add_mesh(mesh, color='lightgrey', smooth_shading=True)
    pl.add_mesh(pv.PolyData(lo+X*size), color='red', point_size=4, render_points_as_spheres=True)
    if chains:
        segs = [lo+NP[c]*size for c in chains]; sc = [np.repeat(cols[ci][None], len(c), 0) for ci, c in enumerate(chains)]
        cl = np.concatenate([np.r_[len(c), np.arange(o, o+len(c))] for c, o in zip(chains, np.cumsum([0]+L[:-1]))])
        poly = pv.PolyData(np.vstack(segs), lines=cl); poly['c'] = np.vstack(sc)
        pl.add_mesh(poly, scalars='c', rgb=True, line_width=4)
        pl.add_mesh(pv.PolyData(np.vstack(segs)), scalars=np.vstack(sc), rgb=True, point_size=9, render_points_as_spheres=True)
    pl.add_text(f'zoom {zoom}: dots red, chains coloured, nodes as spheres', font_size=10)
    focus = np.array([12, 135, 12.]); pl.camera_position = 'zy'; pl.camera.focal_point = focus.tolist(); pl.camera.position = (focus+(np.array(pl.camera.position)-np.array(pl.camera.focal_point))).tolist(); pl.camera.zoom(zoom)
    pl.screenshot(out+f'_zoom{zoom}.png'); print('saved', out+f'_zoom{zoom}.png')
