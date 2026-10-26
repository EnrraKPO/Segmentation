"""The crawl inside the walk, second run: fine dots cached per block of space (computed once, reused by every spawn in
the block), and at every dead end a record of the dots ahead and which criterion threw each one out.
Renders out/walk_crawl2_hair.png; writes out/walk_crawl2_hair.npz (lines, fine dots, end diagnostics)."""
import sys, os, json, time; root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); sys.path.insert(0, root)
import numpy as np
from scipy.spatial import cKDTree
import linesio, creases, trace4, refine, voxel

src = open(os.path.join(root, 'trace4.py'), encoding='utf-8').read()
old = "                if len(cand) == 0: stats['no_anchor'] += 1; reasons.append('no_anchor'); ends_at.append((X[a] if 'a' in dir() else q).copy()); break"
new = """                if len(cand) == 0 and spawn is not None and nspawn < spawn_max:
                    nspawn += 1; newd = spawn(q, d, nodes)
                    if newd is not None and len(newd[0]):
                        Xn, wn, Nn = newd; X = np.vstack([X, Xn]); wgt = np.concatenate([wgt, wn]); N = np.vstack([N, Nn]); owner = np.concatenate([owner, np.full(len(Xn), -1)]); tree = cKDTree(X); n = len(X); stats['spawned'] += len(Xn); stats['spawns'] += 1
                        continue
                if len(cand) == 0:
                    cc = np.array(tree.query_ball_point(q, 3*reach+width), int); rec = dict(q=q.copy(), d=d.copy(), width=width, slack=slack, reach=reach, step=step, ball=len(cc), ahead=0, along=0, across=0, cone=0, wmax=0.0, wmed=0.0, spawns=nspawn, nodes=len(nodes))
                    if len(cc):
                        off = X[cc]-q; al = off @ d; ac = np.linalg.norm(off-al[:, None]*d, axis=1); a0 = (al > 0) & (al <= 3*reach); a1 = (al >= 0.5*step) & (al <= reach); a2 = a1 & (ac <= 0.5*width+slack); a3 = a2 & (al/np.maximum(np.linalg.norm(off, axis=1), 1e-9) >= cos_cone)
                        rec.update(ahead=int(a0.sum()), along=int(a1.sum()), across=int(a2.sum()), cone=int(a3.sum()), wmax=float(wgt[cc].max()), wmed=float(np.median(wgt[cc])), foreign=int(((owner[cc] != -1) & (owner[cc] != tid)).sum()))
                    walk.enddiag.append(rec)
                    stats['no_anchor'] += 1; reasons.append('no_anchor'); ends_at.append((X[a] if 'a' in dir() else q).copy()); break"""
assert old in src; src = src.replace(old, new)
src = src.replace("    walk.probe = []\n", "    walk.probe = []; walk.enddiag = []; spawn = kw.get('spawn'); spawn_max = int(kw.get('spawn_max', 500))\n", 1)
src = src.replace("        for sign in (1, -1):\n            nodes = [halves[0][0] if halves else p0]", "        for sign in (1, -1):\n            nspawn = 0\n            nodes = [halves[0][0] if halves else p0]", 1)
src = src.replace("stats = dict(no_anchor=0,", "stats = dict(spawned=0, spawns=0, no_anchor=0,", 1)
mod = {}; exec(compile(src, 'trace4_crawl', 'exec'), mod); walk_crawl = mod['walk']


class FineDots:
    """fine cubes of the model with their bend, computed per block of space on demand and kept"""
    def __init__(self, P, F, h, cube, flat, finer=3, area=0.087, unit=0.087, block=2.0, nbh=2.0):
        self.P, self.F, self.h, self.flat = P, F, h, flat; self.size = cube/finer; self.cells = h/self.size; self.area, self.unit, self.nbh = area, unit, nbh
        self.block = block; self.blocks = {}; self.lo = P.min(0); self.time = 0.0; self.n = 0
    def _make(self, key):
        t0 = time.time(); b0 = self.lo+np.array(key)*self.block; m = self.nbh*self.area+self.size
        T = self.P[self.F]; hit = np.all((T.max(1) >= b0-m) & (T.min(1) <= b0+self.block+m), axis=1)
        if not hit.any(): self.blocks[key] = (np.zeros((0, 3)), np.zeros(0), np.zeros((0, 3)), np.zeros(0, bool)); return
        try: r = voxel.run_direct(self.P, self.F, self.cells, None, height=self.h, flat=self.flat, box=(b0-m, b0+self.block+m), workers=1, blend=True)
        except (ValueError, IndexError): self.blocks[key] = (np.zeros((0, 3)), np.zeros(0), np.zeros((0, 3)), np.zeros(0, bool)); return   # triangles graze the box but no cube is hit
        C = r['lo']+(r['idx']+0.5)*r['size']
        if len(C) == 0: self.blocks[key] = (C, np.zeros(0), np.zeros((0, 3)), np.zeros(0, bool)); return
        inside = np.all((C >= b0) & (C < b0+self.block), axis=1); at = np.flatnonzero(inside)
        turn, _ = refine.bend(C, r['N'], self.area, self.unit, None, at=at, nbh=self.nbh) if len(at) else (np.zeros(0), None)
        self.blocks[key] = (C[at], turn, r['N'][at], np.zeros(len(at), bool)); self.time += time.time()-t0; self.n += 1
    def take(self, path, radius):
        """fine dots within radius of the polyline path not handed out before: positions (mm), bend, facings"""
        lo = (path.min(0)-radius-self.lo)//self.block; hi = (path.max(0)+radius-self.lo)//self.block; out = []
        for kx in range(int(lo[0]), int(hi[0])+1):
            for ky in range(int(lo[1]), int(hi[1])+1):
                for kz in range(int(lo[2]), int(hi[2])+1):
                    key = (kx, ky, kz)
                    if key not in self.blocks: self._make(key)
                    C, turn, N, used = self.blocks[key]
                    if len(C) == 0: continue
                    d = cKDTree(path).query(C)[0]; sel = (d <= radius) & ~used
                    if sel.any(): used[sel] = True; out.append((C[sel], turn[sel], N[sel]))
        if not out: return None
        return np.vstack([o[0] for o in out]), np.concatenate([o[1] for o in out]), np.vstack([o[2] for o in out])


if __name__ == '__main__':
    import pyvista as pv
    c = np.array([180., 150., 200.]); R = 22.0; Rv = 18.0
    P, F = linesio.load_mesh(os.path.join(root, 'out', 'orcfull_patch.stl'))
    inside = np.all(np.abs(P-c) < R, axis=1); keep = inside[F].all(1); F = F[keep]; used, inv = np.unique(F, return_inverse=True); P = P[used]; F = inv.reshape(-1, 3)
    h = float((P.max(0)-P.min(0)).max())
    S = json.load(open(os.path.join(root, 'out', 'debug', '20261023_114229', 'settings.json'))); cube, cut = S['cube_mm'], S['cut']
    wp = {k: v for k, v in S.items() if isinstance(v, (int, float)) and not isinstance(v, bool) and k not in ('faces', 'cube_mm', 'cut', 'dots', 'lines', 'refined', 'joined')}
    wp = dict(wp, reach=wp['step'])
    res = creases.find_dots(P, F, max(1, int(round(h/cube))), cut, 0.0, S['flat'], None, height=h, blend=S['blend'])
    line = res['line']; size = float(res['size']); lo = res['lo']
    X = (res['idx'][line]+0.5).astype(float); wgt = np.maximum(res['turn'][line]-cut, 0)+1e-3; N = res['N'][line]
    print(f'patch {len(F)} triangles, {len(X)} coarse dots')
    fine = FineDots(P, F, h, cube, S['flat'], finer=refine.DEFAULTS['finer'], area=refine.DEFAULTS['area'], unit=refine.DEFAULTS['unit']); ahead = 1.0; radius = 0.5
    def spawn(q, d, nodes):
        qm = lo+q*size; tail = lo+np.array(nodes[-12:])*size; path = np.vstack([tail, qm+np.outer(np.linspace(0, ahead, 12)[1:], d)])
        got = fine.take(path, radius)
        if got is None: return None
        C, turn, Nn = got; return (C-lo)/size, 2*turn+1e-3, Nn                     # no cut, no floor: the fine bend reads about half the coarse span on one crease
    import vtk
    from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray
    poly = vtk.vtkPolyData(); vp = vtk.vtkPoints(); vp.SetData(numpy_to_vtk(np.ascontiguousarray(P, np.float64), deep=True)); poly.SetPoints(vp)
    cells = vtk.vtkCellArray(); cells.SetData(3, numpy_to_vtkIdTypeArray(np.ascontiguousarray(F, np.int64).ravel(), deep=True)); poly.SetPolys(cells)
    locator = vtk.vtkStaticCellLocator(); locator.SetDataSet(poly); locator.BuildLocator()
    snap_stats = dict(ray=0, nearest=0, moved=[])
    tt = vtk.reference(0.0); xx = [0.0, 0.0, 0.0]; pc = [0.0, 0.0, 0.0]; sub = vtk.reference(0); cid = vtk.reference(0); d2 = vtk.reference(0.0); cp = [0.0, 0.0, 0.0]
    def snap(p, nrm, reach_mm=1.0):
        """the nearest surface point along the line through p (cube units) in the direction nrm, forward or back; the nearest surface point if the line misses"""
        pm = lo+p*size; n = np.asarray(nrm, float); nn = np.linalg.norm(n)
        best = None
        if nn > 1e-12:
            n = n/nn
            for sgn in (1.0, -1.0):
                if locator.IntersectWithLine(pm, pm+sgn*reach_mm*n, 1e-9, tt, xx, pc, sub, cid):
                    d = np.linalg.norm(np.array(xx)-pm)
                    if best is None or d < best[0]: best = (d, np.array(xx))
        if best is None:
            locator.FindClosestPoint(pm, cp, cid, sub, d2); best = (d2.get()**0.5, np.array(cp)); snap_stats['nearest'] += 1
        else: snap_stats['ray'] += 1
        snap_stats['moved'].append(best[0]); return (best[1]-lo)/size
    t0 = time.time(); lines1, owner1 = walk_crawl(X.copy(), wgt.copy(), N.copy(), size, log=print, spawn=spawn, snap=snap, **wp)
    mv = np.array(snap_stats['moved']); print(f'snap: {snap_stats["ray"]} nodes by ray, {snap_stats["nearest"]} by nearest point; moved p50 {np.median(mv):.3f}, p90 {np.percentile(mv, 90):.3f}, max {mv.max():.3f} mm')
    print(f'walk with the spawn: {len(lines1)} lines, {time.time()-t0:.0f} s ({fine.n} blocks of {fine.block} mm computed in {fine.time:.0f} s)')
    L1 = [lo+l*size for l in lines1]; diag = walk_crawl.enddiag
    for r_ in diag: r_['q'] = lo+r_['q']*size
    allp = np.vstack(L1); own = np.repeat(np.arange(len(L1)), [len(l) for l in L1]); tree = cKDTree(allp)
    ends = [(i, e, p) for i, L in enumerate(L1) if np.linalg.norm(L[0]-L[-1]) > 1e-9 for e, p in ((0, L[0]), (1, L[-1])) if all(own[j] == i for j in tree.query_ball_point(p, 1e-6)) and np.all(np.abs(p-c) < Rv)]
    print(f'{len(ends)} dead ends in view; {len(diag)} no-anchor records')
    # match records to the dead ends (the record's q is the last node)
    qs = np.array([r_['q'] for r_ in diag]); qt = cKDTree(qs); rows = []
    for i, e, p in ends:
        dd, j = qt.query(p)
        if dd > 1e-6: rows.append((i, e, None)); continue
        rows.append((i, e, diag[j]))
    with_rec = [r_ for r_ in rows if r_[2] is not None]; print(f'{len(with_rec)} dead ends ended with no anchor (the rest: short lines dropped elsewhere, collisions, patch border)')
    def cls(r_):
        if r_['ball'] == 0: return 'no dots at all within 3 reaches'
        if r_['ahead'] == 0: return 'dots around but none ahead within 3 reaches'
        if r_['along'] == 0: return 'dots ahead but beyond one reach'
        if r_['across'] == 0: return 'dots within reach but outside the band width + slack'
        if r_['cone'] == 0: return 'dots within band but outside the 45 deg cone'
        return 'anchor found?!'
    import collections; print('why the walker stopped:', dict(collections.Counter(cls(r_[2]) for r_ in with_rec)))
    for i, e, r_ in with_rec[:12]: print(f'  line {i} end {e}: ball {r_["ball"]}, ahead(3 reaches) {r_["ahead"]}, along {r_["along"]}, across {r_["across"]}, cone {r_["cone"]}; width {r_["width"]*size:.3f} mm, slack {r_["slack"]*size:.3f}, wmax {r_["wmax"]:.1f}, wmed {r_["wmed"]:.1f}, foreign {r_.get("foreign", 0)}, spawns {r_["spawns"]}')
    np.savez(os.path.join(root, 'out', 'walk_crawl2_hair.npz'), n1=[len(l) for l in L1], p1=np.vstack(L1), diag=np.array(diag, dtype=object), allow_pickle=True)
    mesh = pv.PolyData(P, np.hstack([np.full((len(F), 1), 3), F]).ravel()); box = [c[0]-Rv, c[0]+Rv, c[1]-Rv, c[1]+Rv, c[2]-Rv, c[2]+Rv]
    pl = pv.Plotter(off_screen=True, window_size=(2400, 1800)); pl.add_mesh(mesh.clip_box(box, invert=False), color='lightgrey', smooth_shading=False)
    for L in L1: pl.add_mesh(pv.PolyData(L, lines=np.r_[len(L), np.arange(len(L))]), color='blue', line_width=3)
    ee = np.array([p for _, _, p in ends]).reshape(-1, 3)
    if len(ee): pl.add_mesh(pv.PolyData(ee), color='red', point_size=9, render_points_as_spheres=True)
    pl.camera_position = [c+np.array([-1, 0.3, 2.2])*Rv*1.6, c, [0, 1, 0]]; pl.add_text(f'walk spawning cached fine dots: {len(L1)} lines, {len(ends)} dead ends', font_size=12)
    pl.screenshot(os.path.join(root, 'out', 'walk_crawl2_hair.png')); print('wrote out/walk_crawl2_hair.png')
