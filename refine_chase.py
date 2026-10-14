"""Refine the chains from exp_chase.py by manipulating their NODES (the dots stay as they are).

Three operations, in cycles: TRIM branching (a chain running beside a longer chain loses that stretch), BRIDGE gaps
(two chain ends that point at each other within `gap` cells, with agreeing directions, are joined), SMOOTH (each
inner node moves halfway to the midpoint of its neighbours). Then chains under `min_len` nodes are dropped.
usage: python refine_chase.py mesh nodes.npz out_prefix [gap agree max_div cycles min_len]"""
import sys, os, numpy as np, trimesh, pyvista as pv
from scipy.spatial import cKDTree

path, npz, out = sys.argv[1], sys.argv[2], sys.argv[3]; a = sys.argv[4:]
gap, agree, max_div, cycles, min_len = (float(a[0]) if a else 6.0), (float(a[1]) if len(a) > 1 else 25.0), (float(a[2]) if len(a) > 2 else 45.0), (int(a[3]) if len(a) > 3 else 3), (int(a[4]) if len(a) > 4 else 3)
BESIDE = 1.0    # cells: a node this close to another chain's segment is on the same line
cos_agree, cos_div = np.cos(np.radians(agree)), np.cos(np.radians(max_div))

z = np.load(npz); NP = z['NP'].copy(); X, lo, size, focus = z['X'], z['lo'], z['size'], z['focus']
chains = list(np.split(z['chain'], np.cumsum(z['n'])[:-1])); chains = [list(c) for c in chains]


def end_dir(ch, at_start):
    """outward direction at a chain end, from the chord of the last three nodes"""
    seg = ch[:3] if at_start else ch[-3:][::-1]; v = NP[seg[0]]-NP[seg[-1]]; return v/(np.linalg.norm(v)+1e-12)


def trim(chains):
    """drop nodes that lie beside a longer chain; a chain may split into pieces"""
    order = sorted(range(len(chains)), key=lambda c: -len(chains[c])); out = []; kept = []   # kept: (chain nodes) of longer chains accepted so far
    seg_tree = None; segs = []
    for ci in order:
        ch = chains[ci]; keep = np.ones(len(ch), bool)
        if segs:
            A, B = np.array([s[0] for s in segs]), np.array([s[1] for s in segs]); mid = (A+B)/2; tree = cKDTree(mid)
            for t, nd in enumerate(ch):
                p = NP[nd]; cand = tree.query_ball_point(p, gap)
                for s in cand:
                    ab = B[s]-A[s]; L2 = ab @ ab
                    u = np.clip(((p-A[s]) @ ab)/(L2+1e-12), 0, 1); d = np.linalg.norm(p-(A[s]+u*ab))
                    if d <= BESIDE: keep[t] = False; break
        # split into kept runs
        run = []
        for t, nd in enumerate(ch):
            if keep[t]: run.append(nd)
            elif run: out.append(run); run = []
        if run: out.append(run)
        for t in range(len(ch)-1): segs.append((NP[ch[t]], NP[ch[t+1]]))
    return out


def bridge(chains):
    """join chain ends that point at each other within gap"""
    ends = [(ci, 0) for ci in range(len(chains))]+[(ci, 1) for ci in range(len(chains))]
    E = np.array([NP[chains[ci][0] if e == 0 else chains[ci][-1]] for ci, e in ends]); Dv = np.array([end_dir(chains[ci], e == 0) for ci, e in ends])
    tree = cKDTree(E); pairs = tree.query_pairs(gap, output_type='ndarray'); cand = []
    for a_, b_ in pairs:
        if ends[a_][0] == ends[b_][0]: continue
        v = E[b_]-E[a_]; L = np.linalg.norm(v); u = v/(L+1e-12)
        if (u @ Dv[a_]) >= cos_div and (-u @ Dv[b_]) >= cos_div and (-Dv[a_] @ Dv[b_]) >= cos_agree: cand.append((L, a_, b_))
    cand.sort(); used_end = set(); joins = []
    for L, a_, b_ in cand:
        if a_ in used_end or b_ in used_end: continue
        used_end.add(a_); used_end.add(b_); joins.append((ends[a_], ends[b_]))
    # apply joins via union of chains
    parent = list(range(len(chains))); cur = [list(c) for c in chains]

    def find(x):
        while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for (ca, ea), (cb, eb) in joins:
        ra, rb = find(ca), find(cb)
        if ra == rb: continue
        A, B = cur[ra], cur[rb]
        # orient: A's joined end must be its last node, B's joined end its first
        if A[0] == chains[ca][0] if ea == 0 else A[-1] != chains[ca][-1]: A = A[::-1]
        if B[-1] == chains[cb][-1] if eb == 1 else B[0] != chains[cb][0]: B = B[::-1]
        cur[ra] = A+B; parent[rb] = ra; cur[rb] = None
    return [c for c in cur if c], len(joins)


def smooth(chains):
    for ch in chains:
        if len(ch) < 3: continue
        P = NP[ch]; mid = (P[:-2]+P[2:])/2; NP[ch[1:-1]] = 0.5*(P[1:-1]+mid)


for c in range(1, cycles+1):
    n0 = len(chains); chains = trim(chains); n1 = len(chains); chains, nj = bridge(chains); smooth(chains)
    chains = [ch for ch in chains if len(ch) >= min_len]
    L = [len(ch) for ch in chains]
    print(f'cycle {c}: {n0} -> trim {n1} -> bridge {nj} joins -> {len(chains)} chains, longest {max(L)}, >=20: {sum(l >= 20 for l in L)}, nodes {sum(L)}')

L = [len(ch) for ch in chains]
np.savez(out+'.npz', n=np.array(L), pts=np.vstack([lo+NP[ch]*size for ch in chains]))
m = trimesh.load(path, process=True); P, F = np.asarray(m.vertices), np.asarray(m.faces)
mesh = pv.PolyData(P, np.hstack([np.full((len(F), 1), 3), F]).ravel()); rng = np.random.default_rng(0); cols = rng.uniform(0.1, 0.95, (len(chains), 3))
for zoom in (6, 12):
    pl = pv.Plotter(off_screen=True, window_size=(1500, 1500)); pl.add_mesh(mesh, color='lightgrey', smooth_shading=True)
    pl.add_mesh(pv.PolyData(lo+X*size), color='red', point_size=4, render_points_as_spheres=True)
    segs = [lo+NP[ch]*size for ch in chains]; sc = [np.repeat(cols[ci][None], len(ch), 0) for ci, ch in enumerate(chains)]
    cl = np.concatenate([np.r_[len(ch), np.arange(o, o+len(ch))] for ch, o in zip(chains, np.cumsum([0]+L[:-1]))])
    poly = pv.PolyData(np.vstack(segs), lines=cl); poly['c'] = np.vstack(sc); pl.add_mesh(poly, scalars='c', rgb=True, line_width=4)
    pl.add_mesh(pv.PolyData(np.vstack(segs)), scalars=np.vstack(sc), rgb=True, point_size=9, render_points_as_spheres=True)
    pl.add_text(f'refined, zoom {zoom}: dots red, chains coloured, nodes as spheres', font_size=10)
    pl.camera_position = 'zy'; pl.camera.focal_point = focus.tolist(); pl.camera.position = (focus+(np.array(pl.camera.position)-np.array(pl.camera.focal_point))).tolist(); pl.camera.zoom(zoom)
    pl.screenshot(out+f'_zoom{zoom}.png'); print('saved', out+f'_zoom{zoom}.png')
