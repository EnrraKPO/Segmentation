"""Test flight: chains by DIVERGENCE (the user's definition).

A dot's direction is a tangent, so it never lands exactly on the next dot; it misses by an angle, the divergence.
Two dots in a trail: the divergence from each to the line joining them is below `max_div` degrees (< 90), and a dot
in the middle of a trail has about the same divergence towards both neighbours. Directions are set FROM the
neighbours (the chord through them), not estimated and refined. A dot moves only when it cannot meet the threshold
to two trail dots from where it is: then it moves onto the chord of its would-be neighbours.

Cycle: merge beside dots (as before) -> link by divergence -> set directions from neighbours -> move off-trail dots.
usage: python exp_diverge.py mesh out_prefix [cells sens reach agree merge_reach max_div cycles min_len]"""
import sys, os, numpy as np, trimesh, pyvista as pv
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import voxel, lines

path, out = sys.argv[1], sys.argv[2]; a = sys.argv[3:]
cells, sens, reach, agree, mreach, max_div, cycles, min_len = (int(a[0]) if a else 600), (float(a[1]) if len(a) > 1 else 6.0), (float(a[2]) if len(a) > 2 else 4.0), (float(a[3]) if len(a) > 3 else 25.0), (float(a[4]) if len(a) > 4 else 2.0), (float(a[5]) if len(a) > 5 else 45.0), (int(a[6]) if len(a) > 6 else 6), (int(a[7]) if len(a) > 7 else 10)
LINK_REACH = 2.0

m = trimesh.load(path, process=True); P, F = np.asarray(m.vertices), np.asarray(m.faces)
res = voxel.run(P, F, cells, sens, print)
li = np.nonzero(res['line'])[0]; X = (res['idx'][li]+0.5).astype(float); N = res['N'][li]
D = lines.directions(X, N, reach); alive = np.ones(len(X), bool)
prm = lines.make_ops(mreach, agree, 0.6, 0.5, 4.0)


def umean(d1, d2):
    d2 = d2 if d1 @ d2 >= 0 else -d2; v = d1+d2; return v/(np.linalg.norm(v)+1e-12)


def link_by_divergence(X, D, alive, max_div):
    """ahead/behind per dot: among dots within LINK_REACH, the one on each side with the smallest divergence
    (angle between the dot's direction and the line to it), if under max_div. Sides are the two half-spaces of D."""
    cosd = np.cos(np.radians(max_div)); ia = np.nonzero(alive)[0]; tree = cKDTree(X[ia])
    ahead = np.full(len(X), -1); behind = np.full(len(X), -1); div_a = np.full(len(X), np.nan); div_b = np.full(len(X), np.nan)
    for k in ia:
        nb = ia[tree.query_ball_point(X[k], LINK_REACH)]; nb = nb[nb != k]
        if len(nb) == 0: continue
        v = X[nb]-X[k]; d = np.linalg.norm(v, axis=1); u = v/(d[:, None]+1e-12); c = u @ D[k]
        # divergence from k towards j, and from j back towards k (both must be under the threshold: they point at each other)
        c_back = np.abs((u*D[nb]).sum(1))
        for sgn, slot, dv in ((1, ahead, div_a), (-1, behind, div_b)):
            side = (sgn*c >= cosd) & (c_back >= cosd)
            if side.any():
                j = np.nonzero(side)[0]; best = j[np.argmax(sgn*c[j]+c_back[j]-0.1*d[j])]   # most aligned, nearest breaks ties
                slot[k] = nb[best]; dv[k] = np.degrees(np.arccos(np.clip(sgn*c[best], -1, 1)))
    return ahead, behind, div_a, div_b


def set_directions(X, D, alive, ahead, behind):
    for k in np.nonzero(alive & ((ahead >= 0) | (behind >= 0)))[0]:
        a_ = X[ahead[k]] if ahead[k] >= 0 else X[k]; b_ = X[behind[k]] if behind[k] >= 0 else X[k]
        ch = a_-b_; n = np.linalg.norm(ch)
        if n > 1e-9: D[k] = ch/n


def move_off_trail(X, D, alive, ahead, behind, div_a, div_b, max_div):
    """a dot with links on both sides whose two divergences are unequal or over the threshold moves onto the chord of
    its neighbours (halfway to it); a dot with no link at all stays"""
    k = np.nonzero(alive & (ahead >= 0) & (behind >= 0))[0]; n = 0
    for kk in k:
        a_, b_ = X[ahead[kk]], X[behind[kk]]; ch = a_-b_; L = np.linalg.norm(ch)
        if L < 1e-9: continue
        t = np.clip(((X[kk]-b_) @ ch)/(L*L), 0, 1); foot = b_+t*ch; off = np.linalg.norm(X[kk]-foot)
        if off > 0.25 or abs(div_a[kk]-div_b[kk]) > 0.5*max_div: X[kk] = X[kk]+0.5*(foot-X[kk]); n += 1
    return n


def chains(X, alive, ahead, behind, min_len):
    k = len(X); a = np.nonzero(alive & (ahead >= 0))[0]; b = np.nonzero(alive & (behind >= 0))[0]
    src = np.r_[a, b]; dst = np.r_[ahead[a], behind[b]]
    _, comp = connected_components(coo_matrix((np.ones(len(src)), (src, dst)), shape=(k, k)), directed=False)
    comp = np.where(alive, comp, -1); out = []
    for c in np.unique(comp[comp >= 0]):
        mem = np.nonzero(comp == c)[0]
        if len(mem) < min_len: continue
        ms = set(mem.tolist()); nbrs = {mm: set() for mm in mem}
        for mm in mem:
            for t in (ahead[mm], behind[mm]):
                if t in ms: nbrs[mm].add(t); nbrs[t].add(mm)
        ends = [mm for mm in mem if len(nbrs[mm]) <= 1]; start = ends[0] if ends else mem[0]
        path = [start]; seen = {start}; cur = start
        while True:
            nxt = [t for t in nbrs[cur] if t not in seen]
            if not nxt: break
            cur = min(nxt, key=lambda t: np.linalg.norm(X[t]-X[path[-1]])); path.append(cur); seen.add(cur)
        if len(path) >= min_len: out.append(path)
    return out


for c in range(1, cycles+1):
    nm = lines.merge(X, D, alive, prm)
    ahead, behind, div_a, div_b = link_by_divergence(X, D, alive, max_div)
    set_directions(X, D, alive, ahead, behind)
    ahead, behind, div_a, div_b = link_by_divergence(X, D, alive, max_div)
    nmv = move_off_trail(X, D, alive, ahead, behind, div_a, div_b, max_div)
    ahead, behind, div_a, div_b = link_by_divergence(X, D, alive, max_div)
    k = np.nonzero(alive)[0]; ends = ((ahead[k] < 0) | (behind[k] < 0)).sum(); linked = (ahead[k] >= 0).sum()
    print(f'cycle {c}: merged {nm}, moved {nmv} -> {len(k)} dots, {linked} linked ahead, {ends} ends, divergence p50 {np.nanmedian(div_a[k]):.1f} p90 {np.nanpercentile(div_a[k], 90):.1f} deg')
paths = chains(X, alive, ahead, behind, min_len); L = [len(p) for p in paths]
print(f'chains {len(paths)}, longest {max(L) if L else 0}, dots in chains {sum(L)}')
lo, size = res['lo'], res['size']
np.savez(out+'.npz', n=np.array(L), pts=np.vstack([lo+X[p]*size for p in paths]) if paths else np.zeros((0, 3)))
# close-up with directions and links
k = np.nonzero(alive)[0]; W = lo+X[k]*size; mesh = pv.PolyData(P, np.hstack([np.full((len(F), 1), 3), F]).ravel())
pl = pv.Plotter(off_screen=True, window_size=(1500, 1500)); pl.add_mesh(mesh, color='lightgrey', smooth_shading=True)
pl.add_mesh(pv.PolyData(W), color='red', point_size=5, render_points_as_spheres=True)
seg = np.stack([W-0.6*size*D[k], W+0.6*size*D[k]], 1).reshape(-1, 3); n = len(k)
pl.add_mesh(pv.PolyData(seg, lines=np.c_[np.full(n, 2), np.arange(0, 2*n, 2), np.arange(1, 2*n, 2)].ravel()), color='green', line_width=2)
lk = k[ahead[k] >= 0]
if len(lk):
    seg = np.stack([lo+X[lk]*size, lo+X[ahead[lk]]*size], 1).reshape(-1, 3); n = len(lk)
    pl.add_mesh(pv.PolyData(seg, lines=np.c_[np.full(n, 2), np.arange(0, 2*n, 2), np.arange(1, 2*n, 2)].ravel()), color='blue', line_width=3)
focus = np.array([12, 135, 12.]); pl.camera_position = 'zy'; pl.camera.focal_point = focus.tolist(); pl.camera.position = (focus+(np.array(pl.camera.position)-np.array(pl.camera.focal_point))).tolist(); pl.camera.zoom(12)
pl.screenshot(out+'_zoom.png'); print('saved')
