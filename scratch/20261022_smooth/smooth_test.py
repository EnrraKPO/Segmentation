"""Backward smoothing (layer C) test: the app's coarse dots in a 12 mm box around three spots, walked with the app's
defaults, smoothing off / free / bound to the band. Prints wobble and off-band numbers, renders out/trend_<spot>_<run>.png
and out/trend_<spot>_compare.png (off | free | bound). usage: python smooth_test.py [S]"""
import sys, os, time, numpy as np, trimesh, pyvista as pv
from scipy.spatial import cKDTree
sys.path.insert(0, r'D:\Stuff\PrintSuit\SegmentationTool'); os.chdir(r'D:\Stuff\PrintSuit\SegmentationTool'); import trace4, app
S = float(os.environ.get('S', 0.5))
DOTS = 'out/dots/orcfull_patch_0.087mm_cut17_20261020_164845.npz'
SPOTS = {'rim': ([219.51, 167.43, 165.39], 1.2), 'shield': ([162.13, 137.7, 107.92], 2.5), 'ring': ([194.72, 123.34, 169.49], 2.5)}
BOX = 6.0


def walk_params():
    t = app.TUNE_DEFAULTS; kw = {k: t[k] for k, g, *_ in app.TUNE if g in ('Walk', 'Walk experiments')}
    kw.update(step=0.087, slice=0.7, reach=0.22, pull=10.0, look=0.35, levp=6.0, minlen=8); return kw


def stats(lines, X):
    tree = cKDTree(X); tr = []; off = []
    for l in lines:
        v = np.diff(l, axis=0); v /= np.linalg.norm(v, axis=1, keepdims=True)+1e-12
        if len(v) > 1: tr.append(np.degrees(np.arccos(np.clip((v[1:]*v[:-1]).sum(1), -1, 1))))
        off.append(tree.query(l)[0])
    tr = np.concatenate(tr) if tr else np.zeros(1); off = np.concatenate(off)
    return dict(lines=len(lines), nodes=sum(len(l) for l in lines), t50=np.median(tr), t90=np.percentile(tr, 90), tmax=tr.max(), offn=int((off > 0.087).sum()), offmax=off.max())


def main():
    from PIL import Image, ImageDraw
    z = np.load(DOTS, allow_pickle=True); lo = z['lo']; size = float(z['size']); cut = float(z['cut']); idx, N, turn = z['idx'], z['N'], z['turn']
    C = lo+(idx+0.5)*size; m = trimesh.load('out/orcfull_patch.stl', process=False); P = np.asarray(m.vertices, float); F = np.asarray(m.faces); cen = P[F].mean(1)
    lift = lambda Q: (lambda r: r[0]+0.01*m.face_normals[r[2]])(trimesh.proximity.closest_point(m, Q))
    R = int(os.environ.get('R', 6)); runs = [('off', dict(smooth=0.0, curve=0, trend=6)), ('free', dict(smooth=0.0, curve=1, trend=R)), ('bound', dict(smooth=S, curve=1, trend=R))]
    if 'PULLS' in os.environ:   # pull handle sweep: chord trend, no smoothing, three pull values
        P = [float(x) for x in os.environ['PULLS'].split(',')]; runs = [('off', dict(smooth=0.0, curve=0, trend=6, pull=P[0])), ('free', dict(smooth=0.0, curve=0, trend=6, pull=P[1])), ('bound', dict(smooth=0.0, curve=0, trend=6, pull=P[2]))]
        labels = {n: f'chord trend, pull {v:g}' for n, v in zip(('off', 'free', 'bound'), P)}
    else: labels = None   # off = today's chord; free = arc trend; bound = arc trend + whole-chain smoothing
    print(f'{"spot":7} {"run":6} {"lines":>5} {"nodes":>6} {"turn p50":>8} {"p90":>6} {"max":>6} {"off-band nodes":>14} {"farthest":>8}')
    for spot, (focus, half) in SPOTS.items():
        focus = np.array(focus); sel = (np.abs(C-focus) < BOX).all(1); X = (idx[sel]+0.5).astype(float); wgt = np.maximum(turn[sel]-cut, 0)+1e-3; Nn = N[sel]
        res = {}
        for name, extra in runs:
            t0 = time.time(); lines, _ = trace4.walk(X, wgt, Nn, size, log=lambda *_: None, **{**walk_params(), **extra}); L = [lo+l*size for l in lines]; res[name] = L
            s = stats(L, C[sel]); print(f'{spot:7} {name:6} {s["lines"]:5} {s["nodes"]:6} {s["t50"]:8.1f} {s["t90"]:6.1f} {s["tmax"]:6.1f} {s["offn"]:14} {s["offmax"]:8.3f}  ({time.time()-t0:.1f} s)')
        if (sum(len(l) for l in res['off']) == sum(len(l) for l in res['free'])): pass
        np.savez(f'out/trend_{spot}.npz', **{f'{k}_n': [len(l) for l in v] for k, v in res.items()}, **{f'{k}_pts': np.vstack(v) for k, v in res.items()})
        if os.environ.get('NOREND'): continue
        near = (np.abs(cen-focus) < 2*half).all(1); Fp = F[near]; vv = m.face_normals[near].mean(0); vv /= np.linalg.norm(vv)
        for name, L in res.items():
            mesh = pv.PolyData(P, np.c_[np.full(len(Fp), 3), Fp].ravel()); pl = pv.Plotter(off_screen=True, window_size=(1600, 1200)); pl.add_mesh(mesh, color='lightgrey', smooth_shading=True)
            inv = sel & (np.abs(C-focus) < 2*half).all(1); pts = pv.PolyData(C[inv]+0.5*size*N[inv]); pts['t'] = turn[inv]
            pl.add_mesh(pts, scalars='t', cmap='plasma', clim=(cut, float(np.percentile(turn, 99))), point_size=4, render_points_as_spheres=True, show_scalar_bar=False, reset_camera=False)
            for l in L:
                if np.linalg.norm(l-focus, axis=1).min() < 2*half: ll = lift(l); pl.add_mesh(pv.PolyData(ll, lines=np.r_[len(ll), np.arange(len(ll))]), color='blue', line_width=3, reset_camera=False); pl.add_mesh(pv.PolyData(ll[[0, -1]]), color='black', point_size=8, render_points_as_spheres=True, reset_camera=False)
            pl.camera.focal_point = focus.tolist(); pl.camera.position = (focus+vv*30).tolist(); pl.camera.up = (0, 0, 1); pl.camera.parallel_projection = True; pl.camera.parallel_scale = half
            what = labels[name] if labels else {'off': 'chord trend, no smoothing (today)', 'free': f'arc trend reach {R}, no smoothing', 'bound': f'arc trend reach {R} + smoothing down the chain, strength {S:g}'}[name]
            pl.add_text(f'{what}; dots by strength, lines blue, ends black', font_size=11); pl.screenshot(f'out/trend_{spot}_{name}.png')
        ims = [Image.open(f'out/trend_{spot}_{n}.png') for n, _ in runs]; w, hh = ims[0].size; out = Image.new('RGB', (3*w+40, hh+60), 'white'); dr = ImageDraw.Draw(out)
        for i, (im, (n, _)) in enumerate(zip(ims, runs)): out.paste(im, (i*(w+20), 60)); dr.text((i*(w+20)+20, 15), (('LEFT: ', 'MIDDLE: ', 'RIGHT: ')[i]+labels[n] if labels else {'off': 'LEFT: chord trend (today)', 'free': f'MIDDLE: arc trend, reach {R}', 'bound': f'RIGHT: arc trend reach {R} + smoothing, strength {S:g}'}[n]), fill='black')
        out.resize((out.width//2, out.height//2)).save(f'out/trend_{spot}_'+('pull_' if labels else f'reach{R}_')+'compare.png')
    print('saved out/trend_<spot>_compare.png for', list(SPOTS))


if __name__ == '__main__': main()
