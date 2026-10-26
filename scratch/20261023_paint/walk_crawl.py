"""The crawl inside the walk, on the hair patch. The walk of the app (trace4.walk) with one change: where it finds no
anchor ahead and would end, it first spawns fine dots (refine's tube machinery, finer 3, no cut) in a cloud ahead of its
last node, adds them to the dots it walks on, and looks again. Fine dots are dots like any other. Renders
out/walk_crawl_hair.png: the walk without the spawn (left) and with it (right)."""
import sys, os, json, time, re; root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); sys.path.insert(0, root)
import numpy as np
from scipy.spatial import cKDTree
import linesio, creases, trace4, refine

# the walk's source with the spawn hook, loaded as a module of its own
src = open(os.path.join(root, 'trace4.py'), encoding='utf-8').read()
old = "                if len(cand) == 0: stats['no_anchor'] += 1; reasons.append('no_anchor'); ends_at.append((X[a] if 'a' in dir() else q).copy()); break"
new = """                if len(cand) == 0 and spawn is not None and nspawn < spawn_max:
                    nspawn += 1; newd = spawn(q, d, nodes)
                    if newd is not None and len(newd[0]):
                        Xn, wn, Nn = newd; X = np.vstack([X, Xn]); wgt = np.concatenate([wgt, wn]); N = np.vstack([N, Nn]); owner = np.concatenate([owner, np.full(len(Xn), -1)]); tree = cKDTree(X); n = len(X); stats['spawned'] += len(Xn); stats['spawns'] += 1
                        continue
                if len(cand) == 0: stats['no_anchor'] += 1; reasons.append('no_anchor'); ends_at.append((X[a] if 'a' in dir() else q).copy()); break"""
assert old in src; src = src.replace(old, new)
old2 = "    walk.probe = []\n"
new2 = "    walk.probe = []; spawn = kw.get('spawn'); spawn_max = int(kw.get('spawn_max', 500))\n"
assert old2 in src; src = src.replace(old2, new2, 1)
old3 = "        for sign in (1, -1):\n            nodes = [halves[0][0] if halves else p0]"
new3 = "        for sign in (1, -1):\n            nspawn = 0\n            nodes = [halves[0][0] if halves else p0]"
assert old3 in src; src = src.replace(old3, new3, 1)
src = src.replace("stats = dict(no_anchor=0,", "stats = dict(spawned=0, spawns=0, no_anchor=0,", 1)
mod = {}; exec(compile(src, 'trace4_crawl', 'exec'), mod); walk_crawl = mod['walk']

if __name__ == '__main__':
    import pyvista as pv
    c = np.array([180., 150., 200.]); R = 22.0; Rv = 18.0
    P, F = linesio.load_mesh(os.path.join(root, 'out', 'orcfull_patch.stl'))
    inside = np.all(np.abs(P-c) < R, axis=1); keep = inside[F].all(1); F = F[keep]; used, inv = np.unique(F, return_inverse=True); P = P[used]; F = inv.reshape(-1, 3)
    h = float((P.max(0)-P.min(0)).max())
    S = json.load(open(os.path.join(root, 'out', 'debug', '20261023_114229', 'settings.json'))); cube, cut = S['cube_mm'], S['cut']
    wp = {k: v for k, v in S.items() if isinstance(v, (int, float)) and not isinstance(v, bool) and k not in ('faces', 'cube_mm', 'cut', 'dots', 'lines', 'refined', 'joined')}
    wp = dict(wp, reach=wp['step'])                                                # the walk no longer reaches past a hole
    res = creases.find_dots(P, F, max(1, int(round(h/cube))), cut, 0.0, S['flat'], None, height=h, blend=S['blend'])
    line = res['line']; size = float(res['size']); lo = res['lo']
    X = (res['idx'][line]+0.5).astype(float); wgt = np.maximum(res['turn'][line]-cut, 0)+1e-3; N = res['N'][line]
    print(f'patch {len(F)} triangles, {len(X)} coarse dots')
    finer, area, unit = refine.DEFAULTS['finer'], refine.DEFAULTS['area'], refine.DEFAULTS['unit']; ahead = 1.0; radius = 0.5
    spawned = {'pts': np.zeros((0, 3)), 'tree': None, 'calls': 0, 'time': 0.0}
    def spawn(q, d, nodes):
        """fine dots in a cloud ahead of q (cube units) along d: positions in cube units, weights, facings. The fine bend reads
        about half the coarse span on one crease (measured 2026-10-21), so a fine dot's say is 2*bend-cut like a coarse dot's, floored at a hair."""
        t0 = time.time(); qm = lo+q*size; tail = lo+np.array(nodes[-12:])*size
        path = np.vstack([tail, qm+np.outer(np.linspace(0, ahead, 12)[1:], d)])
        try: t = refine.tube(path, P, F, cube, h, finer=finer, radius=radius, cut=1.0, flat=S['flat'], log=None, workers=1, area=area, unit=unit, rungs=False)
        except Exception: return None
        C = t['lo']+(t['idx']+0.5)*t['size']
        if spawned['tree'] is not None:
            dd = spawned['tree'].query(C)[0]; C_ = C[dd > 0.5*t['size']]; turn = t['turn'][dd > 0.5*t['size']]; Nn = t['N'][dd > 0.5*t['size']]
        else: C_, turn, Nn = C, t['turn'], t['N']
        spawned['pts'] = np.vstack([spawned['pts'], C_]); spawned['tree'] = cKDTree(spawned['pts']); spawned['calls'] += 1; spawned['time'] += time.time()-t0
        return (C_-lo)/size, np.maximum(2*turn-cut, 0)+1e-3, Nn
    t0 = time.time(); lines0, _ = trace4.walk(X, wgt, N, size, log=print, **wp); print(f'walk without the spawn: {len(lines0)} lines, {time.time()-t0:.0f} s')
    t0 = time.time(); lines1, owner1 = walk_crawl(X.copy(), wgt.copy(), N.copy(), size, log=print, spawn=spawn, **wp); print(f'walk with the spawn: {len(lines1)} lines, {time.time()-t0:.0f} s ({spawned["calls"]} spawns, {len(spawned["pts"])} fine dots, {spawned["time"]:.0f} s in the spawns)')
    L0 = [lo+l*size for l in lines0]; L1 = [lo+l*size for l in lines1]
    def dead_ends(lines):
        allp = np.vstack(lines); own = np.repeat(np.arange(len(lines)), [len(l) for l in lines]); tree = cKDTree(allp); out = []
        for i, L in enumerate(lines):
            if np.linalg.norm(L[0]-L[-1]) < 1e-9: continue
            for p in (L[0], L[-1]):
                if all(own[j] == i for j in tree.query_ball_point(p, 1e-6)) and np.all(np.abs(p-c) < Rv): out.append(p)
        return np.array(out).reshape(-1, 3)
    e0, e1 = dead_ends(L0), dead_ends(L1); n0 = sum(np.linalg.norm(np.diff(l, axis=0), axis=1).sum() for l in L0); n1 = sum(np.linalg.norm(np.diff(l, axis=0), axis=1).sum() for l in L1)
    print(f'dead ends in view: {len(e0)} -> {len(e1)}; total line length {n0:.0f} -> {n1:.0f} mm; longest line {max(len(l) for l in L0)} -> {max(len(l) for l in L1)} nodes')
    np.savez(os.path.join(root, 'out', 'walk_crawl_hair.npz'), n0=[len(l) for l in L0], p0=np.vstack(L0), n1=[len(l) for l in L1], p1=np.vstack(L1), fine=spawned['pts'])
    mesh = pv.PolyData(P, np.hstack([np.full((len(F), 1), 3), F]).ravel()); box = [c[0]-Rv, c[0]+Rv, c[1]-Rv, c[1]+Rv, c[2]-Rv, c[2]+Rv]
    def polyline(L): return pv.PolyData(L, lines=np.r_[len(L), np.arange(len(L))])
    pl = pv.Plotter(off_screen=True, shape=(1, 2), window_size=(3000, 1500), border=True)
    for k, (LL, ee) in enumerate(((L0, e0), (L1, e1))):
        pl.subplot(0, k); pl.add_mesh(mesh.clip_box(box, invert=False), color='lightgrey', smooth_shading=False)
        for L in LL: pl.add_mesh(polyline(L), color='blue', line_width=3)
        if len(ee): pl.add_mesh(pv.PolyData(ee), color='red', point_size=9, render_points_as_spheres=True)
        pl.camera_position = [c+np.array([-1, 0.3, 2.2])*Rv*1.6, c, [0, 1, 0]]
        pl.add_text([f'walk, reach = one step: {len(L0)} lines, {len(e0)} dead ends', f'the same walk spawning fine dots where it runs out: {len(L1)} lines, {len(e1)} dead ends'][k], font_size=12)
    pl.screenshot(os.path.join(root, 'out', 'walk_crawl_hair.png')); print('wrote out/walk_crawl_hair.png')
