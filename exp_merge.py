"""Test flight of the merge / align / close-gaps cycle on the corset's line cells.

Dots carry a position (cells) and a direction (unit, unsigned). One cycle:
  merge : two dots within `reach` whose directions agree within `agree` deg and that lie BESIDE each other (the
          vector between them within `beside` deg of the perpendicular to their shared direction) become one dot at
          the average position with the average direction. Merging is done greedily on the closest pairs first.
  align : each dot with a neighbour ahead and one behind (within `reach` along its direction) moves a fraction
          `pull` towards the midpoint of the two; its direction is averaged with theirs.
  gaps  : a dot with no neighbour ahead (an end) whose continuation along its direction reaches, within `gap`
          cells, a dot pointing back at it, is linked to it (a new dot is placed halfway if the gap is > 1 cell).
Repeat `cycles` times. Output: the dots and their links, drawn per cycle on the side cord close-up.
usage: python exp_merge.py mesh out_prefix [cells sens reach agree beside pull gap cycles]"""
import sys, os, numpy as np, trimesh, pyvista as pv
from scipy.spatial import cKDTree
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import voxel, lines

path, out = sys.argv[1], sys.argv[2]
a = sys.argv[3:]; cells, sens = int(a[0]) if a else 1000, float(a[1]) if len(a) > 1 else 6.0
reach, agree, beside, pull, gap, cycles = (float(a[2]) if len(a) > 2 else 2.0), (float(a[3]) if len(a) > 3 else 25.0), (float(a[4]) if len(a) > 4 else 35.0), (float(a[5]) if len(a) > 5 else 0.5), (float(a[6]) if len(a) > 6 else 4.0), (int(a[7]) if len(a) > 7 else 6)

m = trimesh.load(path, process=True); P, F = np.asarray(m.vertices), np.asarray(m.faces)
res = voxel.run(path and P, F, cells, sens, print)
li = np.nonzero(res['line'])[0]; X = (res['idx'][li]+0.5).astype(float); N = res['N'][li]
D, an, nn, tree, i, j = lines.directions(X, N, 4.0)            # direction per dot from its neighbourhood (verified)
alive = np.ones(len(X), bool)
cos_agree, cos_beside = np.cos(np.radians(agree)), np.cos(np.radians(90-beside))   # beside: angle to the direction >= 90-beside


def unsign_mean(d1, d2):
    d2 = d2 if d1 @ d2 >= 0 else -d2; v = d1+d2; return v/(np.linalg.norm(v)+1e-12)


def merge(X, D, alive):
    idx_alive = np.nonzero(alive)[0]; tree = cKDTree(X[idx_alive]); pairs = tree.query_pairs(reach, output_type='ndarray')
    if len(pairs) == 0: return 0
    p, q = idx_alive[pairs[:, 0]], idx_alive[pairs[:, 1]]
    v = X[q]-X[p]; dist = np.linalg.norm(v, axis=1); u = v/(dist[:, None]+1e-12)
    par = np.abs((D[p]*D[q]).sum(1)) >= cos_agree                     # same direction
    dm = D[p]+np.where((D[p]*D[q]).sum(1)[:, None] < 0, -D[q], D[q]); dm /= np.linalg.norm(dm, axis=1, keepdims=True)+1e-12
    bes = np.abs((u*dm).sum(1)) <= cos_beside                           # beside: across the shared direction
    ok = par & bes; p, q, dist = p[ok], q[ok], dist[ok]
    order = np.argsort(dist); taken = np.zeros(len(X), bool); n = 0
    for a_, b_ in zip(p[order], q[order]):
        if taken[a_] or taken[b_]: continue
        X[a_] = 0.5*(X[a_]+X[b_]); D[a_] = unsign_mean(D[a_], D[b_]); alive[b_] = False; taken[a_] = taken[b_] = True; n += 1
    return n


def links(X, D, alive):
    """for each live dot: nearest live dot ahead and behind along its direction within reach (and beside-test passed as NOT beside)"""
    idx_alive = np.nonzero(alive)[0]; tree = cKDTree(X[idx_alive]); ahead = np.full(len(X), -1); behind = np.full(len(X), -1)
    for k in idx_alive:
        nb = idx_alive[tree.query_ball_point(X[k], reach)]; nb = nb[nb != k]
        if len(nb) == 0: continue
        v = X[nb]-X[k]; dist = np.linalg.norm(v, axis=1); u = v/(dist[:, None]+1e-12); along = u @ D[k]
        par = np.abs((D[nb]*D[k]).sum(1)) >= cos_agree; on_line = (np.abs(along) > cos_beside) & par
        fw = nb[on_line & (along > 0)]; bw = nb[on_line & (along < 0)]
        if len(fw): ahead[k] = fw[np.argmin(dist[on_line & (along > 0)])]
        if len(bw): behind[k] = bw[np.argmin(dist[on_line & (along < 0)])]
    return ahead, behind


def align(X, D, alive, ahead, behind):
    both = alive & (ahead >= 0) & (behind >= 0); k = np.nonzero(both)[0]
    mid = 0.5*(X[ahead[k]]+X[behind[k]]); Xn = X.copy(); Xn[k] = X[k]+pull*(mid-X[k])
    for kk in k: D[kk] = unsign_mean(unsign_mean(D[kk], D[ahead[kk]]), D[behind[kk]])
    return Xn, len(k)


def close_gaps(X, D, alive, ahead, behind):
    idx_alive = np.nonzero(alive)[0]; tree = cKDTree(X[idx_alive]); n = 0; newX, newD = [], []
    ends = [k for k in idx_alive if ahead[k] < 0 or behind[k] < 0]
    for k in ends:
        for sgn, slot in ((1, ahead), (-1, behind)):
            if slot[k] >= 0: continue
            probe = X[k]+sgn*D[k]*np.arange(1, int(gap)+1)[:, None]
            for pp in probe:
                nb = idx_alive[tree.query_ball_point(pp, 1.0)]; nb = nb[nb != k]
                if len(nb) == 0: continue
                par = np.abs((D[nb]*D[k]).sum(1)) >= cos_agree
                if not par.any(): continue
                t = nb[par][0]; dist = np.linalg.norm(X[t]-X[k])
                if dist > 1.5: newX.append(0.5*(X[t]+X[k])); newD.append(unsign_mean(D[k], D[t]))
                n += 1; break
    if newX:
        X = np.vstack([X, np.array(newX)]); D = np.vstack([D, np.array(newD)]); alive = np.r_[alive, np.ones(len(newX), bool)]
    return X, D, alive, n


focus = np.array([12, 135, 12.]); lo, size = res['lo'], res['size']
def draw(X, D, alive, ahead, tag, out_png):
    Wd = lo+X*size; mesh = pv.PolyData(P, np.hstack([np.full((len(F), 1), 3), F]).ravel())
    pl = pv.Plotter(off_screen=True, window_size=(1500, 1500)); pl.add_mesh(mesh, color='lightgrey', smooth_shading=True)
    k = np.nonzero(alive)[0]; pl.add_mesh(pv.PolyData(Wd[k]), color='red', point_size=6, render_points_as_spheres=True)
    lk = k[ahead[k] >= 0]
    if len(lk):
        seg = np.stack([Wd[lk], Wd[ahead[lk]]], 1).reshape(-1, 3); ln = pv.PolyData(seg, lines=np.c_[np.full(len(lk), 2), np.arange(0, 2*len(lk), 2), np.arange(1, 2*len(lk), 2)].ravel()); pl.add_mesh(ln, color='blue', line_width=3)
    pl.add_text(tag, font_size=12); pl.camera_position = 'zy'; pl.camera.focal_point = focus.tolist(); pl.camera.position = (focus+(np.array(pl.camera.position)-np.array(pl.camera.focal_point))).tolist(); pl.camera.zoom(6)
    pl.screenshot(out_png)

ahead, behind = links(X, D, alive); draw(X, D, alive, ahead, f'cycle 0: {alive.sum()} dots, {(ahead >= 0).sum()} links', out+'_c0.png')
for c in range(1, cycles+1):
    nm = merge(X, D, alive); ahead, behind = links(X, D, alive); X, na = align(X, D, alive, ahead, behind)
    ahead, behind = links(X, D, alive); X, D, alive, ng = close_gaps(X, D, alive, ahead, behind); ahead, behind = links(X, D, alive)
    ends = (alive & ((ahead < 0) | (behind < 0))).sum()
    print(f'cycle {c}: merged {nm}, aligned {na}, gaps closed {ng} -> {alive.sum()} dots, {(ahead >= 0).sum()} links, {ends} ends')
    if c in (1, 2, cycles): draw(X, D, alive, ahead, f'cycle {c}: {alive.sum()} dots, {ends} ends', out+f'_c{c}.png')
np.savez(out+'.npz', X=lo+X[alive]*size, D=D[alive])
print('saved')
