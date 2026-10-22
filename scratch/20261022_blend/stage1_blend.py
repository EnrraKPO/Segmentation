"""Stage 1 with the blended cube arrow against the nearest-triangle arrow, on the orc patch, app defaults (0.087 mm, cut 17),
then the stage 2 walk on both dot sets. Renders to out/blend1_<spot>_<off|on>.png: dots coloured by strength (plasma,
cut..99th pct), walk lines blue, legend in the picture. usage: python stage1_blend.py"""
import sys, os, time, numpy as np, trimesh, pyvista as pv
from scipy.spatial import cKDTree
sys.path.insert(0, r'D:\Stuff\PrintSuit\SegmentationTool'); os.chdir(r'D:\Stuff\PrintSuit\SegmentationTool'); import creases, trace4

WALK = dict(step=0.087, slice=0.7, reach=0.22, pull=10.0, look=0.35, levp=6.0, minlen=8)
SPOTS = {'rim': ([219.51, 167.43, 165.39], 1.2), 'shield': ([162.13, 137.7, 107.92], 2.5), 'ring': ([194.72, 123.34, 169.49], 2.5)}   # focus, half-view (mm)


def main():
    m = trimesh.load('out/orcfull_patch.stl', process=False); P = np.asarray(m.vertices, float); F = np.asarray(m.faces); h = float((P.max(0)-P.min(0)).max())
    cells = max(1, int(round(h/0.087))); out = {}
    for name, blend in (('off', False), ('on', True)):
        t0 = time.time(); res = creases.find_dots(P, F, cells, 17.0, 0.0, True, print, height=h, blend=blend); t1 = time.time()-t0
        t0 = time.time(); lines = trace4.run(res, log=lambda *_: None, **WALK); t2 = time.time()-t0
        print(f'blend {name}: cubes {len(res["idx"])}, dots {int(res["line"].sum())}, dots {t1:.0f} s; walk: {len(lines)} lines, {sum(len(l) for l in lines)} nodes, {t2:.0f} s')
        out[name] = (res, lines)
    r0, r1 = out['off'][0], out['on'][0]
    ang = np.degrees(np.arccos(np.clip((r0['N']*r1['N']).sum(1), -1, 1)))
    print(f'arrow change over all cubes: changed > 1 deg {int((ang > 1).sum())} of {len(ang)}, p90 {np.percentile(ang, 90):.1f}, max {ang.max():.1f} deg')
    print(f'dot strength (span) change on the common dots: ', end='')
    both = r0['line'] & r1['line']; d = r1['turn'][both]-r0['turn'][both]; print(f'{int(both.sum())} common, p50 {np.median(d):+.1f} p10 {np.percentile(d, 10):+.1f} p90 {np.percentile(d, 90):+.1f} deg; only off {int((r0["line"] & ~r1["line"]).sum())}, only on {int((r1["line"] & ~r0["line"]).sum())}')
    np.savez('out/blend1_dots.npz', idx=r0['idx'], N0=r0['N'], N1=r1['N'], turn0=r0['turn'], turn1=r1['turn'], line0=r0['line'], line1=r1['line'], lo=r0['lo'], size=r0['size'],
             n0=[len(l) for l in out['off'][1]], p0=np.vstack(out['off'][1]), n1=[len(l) for l in out['on'][1]], p1=np.vstack(out['on'][1]))
    lift = lambda Q: (lambda r: r[0]+0.01*m.face_normals[r[2]])(trimesh.proximity.closest_point(m, Q))
    cen = P[F].mean(1)
    for spot, (focus, half) in SPOTS.items():
        focus = np.array(focus); near = (np.abs(cen-focus) < 2*half).all(1); Fp = F[near]; vv = m.face_normals[near].mean(0); vv /= np.linalg.norm(vv)
        for name, (res, lines) in out.items():
            C = res['lo']+(res['idx'][res['line']]+0.5)*res['size']; sel = (np.abs(C-focus) < 2*half).all(1)
            mesh = pv.PolyData(P, np.c_[np.full(len(Fp), 3), Fp].ravel()); pl = pv.Plotter(off_screen=True, window_size=(1600, 1200)); pl.add_mesh(mesh, color='lightgrey', smooth_shading=True)
            pts = pv.PolyData(C[sel]+0.5*res['size']*res['N'][res['line']][sel]); pts['turn'] = res['turn'][res['line']][sel]
            pl.add_mesh(pts, scalars='turn', cmap='plasma', clim=(17, float(np.percentile(res['turn'][res['line']], 99))), point_size=5, render_points_as_spheres=True, show_scalar_bar=False, reset_camera=False)
            for l in lines:
                if np.linalg.norm(l-focus, axis=1).min() < 2*half: ll = lift(l); pl.add_mesh(pv.PolyData(ll, lines=np.r_[len(ll), np.arange(len(ll))]), color='blue', line_width=3, reset_camera=False)
            pl.camera.focal_point = focus.tolist(); pl.camera.position = (focus+vv*30).tolist(); pl.camera.up = (0, 0, 1); pl.camera.parallel_projection = True; pl.camera.parallel_scale = half
            pl.add_text(f'stage 1 blend {name}: dots = cubes above cut 17 coloured by strength (dark = 17, yellow = strong); blue = walk lines', font_size=11)
            pl.screenshot(f'out/blend1_{spot}_{name}.png')
    print('saved out/blend1_{rim,shield,ring}_{off,on}.png and out/blend1_dots.npz')


if __name__ == '__main__': main()
