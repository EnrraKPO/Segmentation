"""The crawl, tried on the hair patch of the orc.
  1. patch: triangles of the orc inside a box around the hair
  2. coarse dots (stage 1, settings of dump 20261023_114229) and the anchored walk with the reach set to one step: a line
     dies where the dots stop instead of grabbing a dot beside the trail
  3. crawl: at each dead end a cloud of fine cubes (refine's tube machinery, finer 3, no cut) around the end's last
     millimetre and 1.5 mm ahead; the fine walk (refine.WALK) runs on it; what it draws is shown
Renders out/crawl_hair.png (overview) and out/crawl_hair_ends.png (close-ups)."""
import sys, os, json, time; root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); sys.path.insert(0, root)
import numpy as np, pyvista as pv
from scipy.spatial import cKDTree
import linesio, creases, trace4, refine, voxel

if __name__ == '__main__':
    c = np.array([180., 150., 200.]); R = 22.0; Rv = 18.0
    P, F = linesio.load_mesh(os.path.join(root, 'out', 'orcfull_patch.stl'))
    inside = np.all(np.abs(P-c) < R, axis=1); keep = inside[F].all(1); F = F[keep]
    used, inv = np.unique(F, return_inverse=True); P = P[used]; F = inv.reshape(-1, 3)
    h = float((P.max(0)-P.min(0)).max()); print(f'patch: {len(F)} triangles, extent {h:.1f} mm')
    S = json.load(open(os.path.join(root, 'out', 'debug', '20261023_114229', 'settings.json')))
    cube, cut = S['cube_mm'], S['cut']; wp = {k: v for k, v in S.items() if isinstance(v, (int, float)) and not isinstance(v, bool) and k not in ('faces', 'cube_mm', 'cut', 'dots', 'lines', 'refined', 'joined')}
    # 2. coarse dots and the walk, reach = one step
    t0 = time.time(); res = creases.find_dots(P, F, max(1, int(round(h/cube))), cut, 0.0, S['flat'], print, height=h, blend=S['blend']); print(f'dots in {time.time()-t0:.1f} s')
    wp_die = dict(wp, reach=wp['step'])
    t0 = time.time(); lines = trace4.run(res, log=print, **wp_die); print(f'walk (reach = one step {wp["step"]} mm): {len(lines)} lines in {time.time()-t0:.1f} s')
    lines_old = trace4.run(res, log=lambda *_: None, **wp); print(f'walk as in the app (reach {wp["reach"]} mm): {len(lines_old)} lines')
    # dead ends: ends not on another line's node
    allp = np.vstack(lines); own = np.repeat(np.arange(len(lines)), [len(l) for l in lines]); tree = cKDTree(allp)
    ends = []
    for i, L in enumerate(lines):
        if np.linalg.norm(L[0]-L[-1]) < 1e-9: continue
        for e, p in ((0, L[0]), (1, L[-1])):
            if all(own[j] == i for j in tree.query_ball_point(p, 1e-6)) and np.all(np.abs(p-c) < Rv): ends.append((i, e))
    print(f'{len(ends)} dead ends inside the view box')
    # 3. the crawl
    fw = dict(refine.WALK, unit=refine.DEFAULTS['unit']); finer = refine.DEFAULTS['finer']; area = refine.DEFAULTS['area']; ahead = 1.5; back = 1.0; radius = 0.4
    crawls = []; t0 = time.time()
    for k, (i, e) in enumerate(ends):
        L = lines[i] if e == 1 else lines[i][::-1]
        s = np.r_[0, np.cumsum(np.linalg.norm(np.diff(L, axis=0), axis=1))]; tail = L[s >= s[-1]-back]
        d = tail[-1]-tail[0] if len(tail) > 1 else L[-1]-L[-2]; d /= np.linalg.norm(d)+1e-12
        path = np.vstack([tail, L[-1]+np.outer(np.linspace(0, ahead, 18)[1:], d)])            # the end's last mm and a straight guess ahead: only the cloud's extent
        t = refine.tube(path, P, F, cube, h, finer=finer, radius=radius, cut=1.0, flat=S['flat'], log=None, workers=1, area=area, unit=refine.DEFAULTS['unit'], rungs=False)
        fl = trace4.run(t, log=lambda *_: None, workers=1, **fw)
        crawls.append((i, e, path, t, fl))
        if k % 20 == 0: print(f'  crawl {k}/{len(ends)}: {len(t["idx"])} fine cubes, {len(fl)} fine lines, {time.time()-t0:.0f} s')
    # what did the crawl do at each end: the fine line through the end, how far past it along d, and whether it reaches another coarse line
    rows = []
    for i, e, path, t, fl in crawls:
        end = path[len(path)-18+1-1] if False else (lines[i][-1] if e == 1 else lines[i][0]); d = path[-1]-end; d /= np.linalg.norm(d)
        best = None
        for f_ in fl:
            dn = np.linalg.norm(f_-end, axis=1)
            if dn.min() < 0.15:
                reach_ = float(((f_-end) @ d).max())
                if best is None or reach_ > best[0]: best = (reach_, f_)
        if best is None: rows.append((i, e, 'no fine line through the end', 0.0, False)); continue
        f_ = best[1]; near_other = any(own[j] != i for j in tree.query_ball_point(f_[np.argmax((f_-end) @ d)], 0.15))
        rows.append((i, e, 'continued', best[0], near_other))
    cont = [r for r in rows if r[2] == 'continued']; print(f'crawl: {len(cont)} of {len(rows)} ends continued by a fine line through the end; of those {sum(r[4] for r in cont)} reach another coarse line; distance past the end p50 {np.median([r[3] for r in cont]) if cont else 0:.2f} mm, p90 {np.percentile([r[3] for r in cont], 90) if cont else 0:.2f} mm')
    # renders
    mesh = pv.PolyData(P, np.hstack([np.full((len(F), 1), 3), F]).ravel()); box = [c[0]-Rv, c[0]+Rv, c[1]-Rv, c[1]+Rv, c[2]-Rv, c[2]+Rv]
    def polyline(L): return pv.PolyData(L, lines=np.r_[len(L), np.arange(len(L))])
    pl = pv.Plotter(off_screen=True, shape=(1, 2), window_size=(3000, 1500), border=True)
    for k in range(2):
        pl.subplot(0, k); pl.add_mesh(mesh.clip_box(box, invert=False), color='lightgrey', smooth_shading=False)
        for L in (lines_old if k == 0 else lines): pl.add_mesh(polyline(L), color='blue', line_width=2)
        if k == 1:
            for i, e, path, t, fl in crawls:
                for f_ in fl: pl.add_mesh(polyline(f_), color='magenta', line_width=3)
            endpts = np.array([lines[i][-1] if e == 1 else lines[i][0] for i, e in ends]); pl.add_mesh(pv.PolyData(endpts), color='red', point_size=9, render_points_as_spheres=True)
        pl.camera_position = [c+np.array([-1, 0.3, 2.2])*Rv*1.6, c, [0, 1, 0]]; pl.add_text(['walk as in the app (reach 0.22 mm)', 'walk with reach = one step (blue), dead ends (red), crawl lines (magenta)'][k], font_size=12)
    pl.screenshot(os.path.join(root, 'out', 'crawl_hair.png')); print('wrote out/crawl_hair.png')
    pl = pv.Plotter(off_screen=True, shape=(2, 4), window_size=(3200, 1600), border=True)
    for k, (i, e, path, t, fl) in enumerate(crawls[:8]):
        pl.subplot(k//4, k % 4); end = lines[i][-1] if e == 1 else lines[i][0]; r = 1.2
        b = [end[0]-r, end[0]+r, end[1]-r, end[1]+r, end[2]-r, end[2]+r]; pl.add_mesh(mesh.clip_box(b, invert=False), color='lightgrey', smooth_shading=False)
        C = t['lo']+(t['idx']+0.5)*t['size']; cl = pv.PolyData(C); cl['bend'] = t['turn']; pl.add_mesh(cl, scalars='bend', cmap='viridis', point_size=3, render_points_as_spheres=True, show_scalar_bar=False, opacity=0.6)
        for L in lines:
            if np.any(np.all(np.abs(L-end) < r, axis=1)): pl.add_mesh(polyline(L), color='blue', line_width=3)
        for f_ in fl: pl.add_mesh(polyline(f_), color='magenta', line_width=4)
        pl.add_mesh(pv.PolyData(end.reshape(1, 3)), color='red', point_size=14, render_points_as_spheres=True)
        n = np.cross(P[F[:, 1]]-P[F[:, 0]], P[F[:, 2]]-P[F[:, 0]]); fi = cKDTree(P[F].mean(1)).query(end)[1]; n = n[fi]/np.linalg.norm(n[fi])
        pl.camera_position = [end+3*r*n, end, [0, 0, 1] if abs(n[2]) < 0.9 else [1, 0, 0]]; pl.add_text(f'end {k}: line {i}, {rows[k][2]}, {rows[k][3]:.2f} mm past the end', font_size=11)
    pl.screenshot(os.path.join(root, 'out', 'crawl_hair_ends.png')); print('wrote out/crawl_hair_ends.png')
