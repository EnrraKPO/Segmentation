"""The crawl with a moving cloud, on the hair patch: at each dead end, a fine cloud ahead, the fine walk through it, the new
end becomes the start of the next cloud, until the fine walk stops inside a cloud, runs into another coarse line, or the
round limit is hit. Renders out/crawl2_hair.png and out/crawl2_ends.png; prints how each crawl ended."""
import sys, os, json, time; root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); sys.path.insert(0, root)
import numpy as np
from scipy.spatial import cKDTree
import linesio, creases, trace4, refine

_G = {}
def _init(P, F, h, cube, flat): _G.update(P=P, F=F, h=h, cube=cube, flat=flat)

def crawl_end(args):
    """one dead end: returns (pieces of fine line past the end in order, how it ended, rounds)"""
    L, others_pts, others_own, own_i, prm = args; P, F, h, cube, flat = _G['P'], _G['F'], _G['h'], _G['cube'], _G['flat']
    ahead, radius, back, rounds, finer, area, unit, fw = prm
    tree = cKDTree(others_pts) if len(others_pts) else None
    cur = np.asarray(L, float); pieces = []; how = 'rounds'
    for r in range(rounds):
        s = np.r_[0, np.cumsum(np.linalg.norm(np.diff(cur, axis=0), axis=1))]; tail = cur[s >= s[-1]-back]
        d = tail[-1]-tail[0] if len(tail) > 1 else cur[-1]-cur[-2]; d /= np.linalg.norm(d)+1e-12; end = cur[-1]
        path = np.vstack([tail, end+np.outer(np.linspace(0, ahead, 12)[1:], d)])
        try: t = refine.tube(path, P, F, cube, h, finer=finer, radius=radius, cut=1.0, flat=flat, log=None, workers=1, area=area, unit=unit, rungs=False)
        except Exception: how = 'tube failed'; break
        if len(t['idx']) < 20: how = 'no cubes'; break
        fl = trace4.run(t, log=lambda *_: None, workers=1, **fw)
        best = None
        for f_ in fl:
            dn = np.linalg.norm(f_-end, axis=1)
            if dn.min() < 0.15:
                al = (f_-end) @ d; k = int(np.argmax(al))
                if best is None or al[k] > best[0]: best = (float(al[k]), f_, k, int(dn.argmin()))
        if best is None: how = 'no fine line through the end'; break
        al, f_, k, k0 = best
        fwd = f_[k0:k+1] if k >= k0 else f_[k:k0+1][::-1]                                 # the part past the end, oriented away from it
        fwd = fwd[np.r_[True, (fwd[1:]-end) @ d > 0.02]] if len(fwd) > 1 else fwd
        if len(fwd) < 3 or al < 0.15: how = 'stopped at the end'; break
        pieces.append(fwd); cur = np.vstack([cur, fwd[1:]])
        if tree is not None:
            dd, j = tree.query(fwd); hit = np.flatnonzero((dd < 0.12) & (others_own[j] != own_i))
            if len(hit): how = 'reached another line'; break
        if al < ahead-0.35: how = 'stopped inside the cloud'; break
    return pieces, how, r+1

if __name__ == '__main__':
    import pyvista as pv
    from concurrent.futures import ProcessPoolExecutor
    c = np.array([180., 150., 200.]); R = 22.0; Rv = 18.0
    P, F = linesio.load_mesh(os.path.join(root, 'out', 'orcfull_patch.stl'))
    inside = np.all(np.abs(P-c) < R, axis=1); keep = inside[F].all(1); F = F[keep]; used, inv = np.unique(F, return_inverse=True); P = P[used]; F = inv.reshape(-1, 3)
    h = float((P.max(0)-P.min(0)).max())
    S = json.load(open(os.path.join(root, 'out', 'debug', '20261023_114229', 'settings.json'))); cube, cut = S['cube_mm'], S['cut']
    wp = {k: v for k, v in S.items() if isinstance(v, (int, float)) and not isinstance(v, bool) and k not in ('faces', 'cube_mm', 'cut', 'dots', 'lines', 'refined', 'joined')}
    res = creases.find_dots(P, F, max(1, int(round(h/cube))), cut, 0.0, S['flat'], None, height=h, blend=S['blend'])
    lines = trace4.run(res, log=lambda *_: None, **dict(wp, reach=wp['step'])); print(f'patch {len(F)} triangles; walk with reach = one step: {len(lines)} lines')
    allp = np.vstack(lines); own = np.repeat(np.arange(len(lines)), [len(l) for l in lines]); tree = cKDTree(allp)
    ends = [(i, e) for i, L in enumerate(lines) if np.linalg.norm(L[0]-L[-1]) > 1e-9 for e, p in ((0, L[0]), (1, L[-1]))
            if all(own[j] == i for j in tree.query_ball_point(p, 1e-6)) and np.all(np.abs(p-c) < Rv)]
    prm = (1.0, 0.5, 1.0, 25, refine.DEFAULTS['finer'], refine.DEFAULTS['area'], refine.DEFAULTS['unit'], dict(refine.WALK, unit=refine.DEFAULTS['unit']))
    tasks = [(lines[i] if e == 1 else lines[i][::-1], allp, own, i, prm) for i, e in ends]
    t0 = time.time()
    with ProcessPoolExecutor(32, initializer=_init, initargs=(P, F, h, cube, S['flat'])) as ex: out = list(ex.map(crawl_end, tasks, chunksize=1))
    print(f'{len(ends)} dead ends crawled in {time.time()-t0:.0f} s')
    import collections; print('how the crawls ended:', dict(collections.Counter(o[1] for o in out)))
    lens = [sum(np.linalg.norm(np.diff(p, axis=0), axis=1).sum() for p in o[0]) for o in out]; print(f'length crawled: p10 {np.percentile(lens, 10):.1f}, p50 {np.median(lens):.1f}, p90 {np.percentile(lens, 90):.1f}, max {max(lens):.1f} mm; rounds p50 {np.median([o[2] for o in out]):.0f}')
    np.save(os.path.join(root, 'out', 'crawl2_hair.npy'), np.array([(ends[k], out[k][1], out[k][2], lens[k]) for k in range(len(out))], dtype=object), allow_pickle=True)
    mesh = pv.PolyData(P, np.hstack([np.full((len(F), 1), 3), F]).ravel()); box = [c[0]-Rv, c[0]+Rv, c[1]-Rv, c[1]+Rv, c[2]-Rv, c[2]+Rv]
    def polyline(L): return pv.PolyData(L, lines=np.r_[len(L), np.arange(len(L))])
    colour = {'stopped inside the cloud': 'magenta', 'reached another line': 'orange', 'rounds': 'cyan', 'no fine line through the end': 'black', 'stopped at the end': 'black', 'no cubes': 'black', 'tube failed': 'black'}
    pl = pv.Plotter(off_screen=True, window_size=(2400, 1800))
    pl.add_mesh(mesh.clip_box(box, invert=False), color='lightgrey', smooth_shading=False)
    for L in lines: pl.add_mesh(polyline(L), color='blue', line_width=3)
    for (i, e), (pieces, how, r) in zip(ends, out):
        for p in pieces: pl.add_mesh(polyline(p), color=colour[how], line_width=3)
    endpts = np.array([lines[i][-1] if e == 1 else lines[i][0] for i, e in ends]); pl.add_mesh(pv.PolyData(endpts), color='red', point_size=8, render_points_as_spheres=True)
    pl.camera_position = [c+np.array([-1, 0.3, 2.2])*Rv*1.6, c, [0, 1, 0]]
    pl.add_text('coarse lines blue; crawls: magenta = stopped inside a cloud, orange = reached another line, cyan = still going at 25 rounds, black = nothing to continue', font_size=11)
    pl.screenshot(os.path.join(root, 'out', 'crawl2_hair.png')); print('wrote out/crawl2_hair.png')
