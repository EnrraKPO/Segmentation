"""Summit dots: is the cube near the summit of a bend? Along the across direction (the neighbour direction of the largest
arrow angle, as fold_change uses), the rate of turn at a position x is the angle between the arrows r cubes behind and
r cubes ahead of x. The rise at the cube = rate(0) - mean(rate(-R), rate(+R)): how much higher the turn is here than at a
distance R on either side. A smoothly curving surface rises 0; a bend, soft or sharp, rises by its summit height.
Where a side reading leaves the surface (no cube there), that side counts 0, so the rise falls back to the rate itself.
Renders out/summit_<spot>_<name>.png for the current span (cut 17) and the rise at R = 2, 3 cubes (cuts 17 and 10).
usage: python summit_test.py"""
import sys, os, time, numpy as np, trimesh, pyvista as pv
sys.path.insert(0, r'D:\Stuff\PrintSuit\SegmentationTool'); os.chdir(r'D:\Stuff\PrintSuit\SegmentationTool'); import voxel, creases

OFFS = creases.OFFS; KB = creases.KB
SPOTS = {'rim': ([219.51, 167.43, 165.39], 1.2), 'shield': ([162.13, 137.7, 107.92], 2.5), 'ring': ([194.72, 123.34, 169.49], 2.5), 'rimwide': ([219.51, 167.43, 165.39], 6.0)}


def arrows_at(idx, N, key_sorted, srt, k_off):
    """the arrow of the cube at idx + k_off (per cube), or NaN where there is no cube"""
    kk = ((idx[:, 0]+k_off[:, 0])*KB+(idx[:, 1]+k_off[:, 1]))*KB+(idx[:, 2]+k_off[:, 2])
    pos = np.minimum(np.searchsorted(key_sorted, kk), len(key_sorted)-1); j = np.where(key_sorted[pos] == kk, srt[pos], -1)
    out = np.full((len(idx), 3), np.nan); h = j >= 0; out[h] = N[j[h]]; return out


def angle(a, b):
    d = (a*b).sum(1); return np.degrees(np.arccos(np.clip(d, -1, 1)))


def summit(idx, N, pdir, r, R):
    """rise per cube (deg): rate at the cube minus the mean rate R cubes either side along the across direction; rate at x =
    angle between the arrows x-r and x+r"""
    key = (idx[:, 0]*KB+idx[:, 1])*KB+idx[:, 2]; srt = np.argsort(key, kind='stable'); ks = key[srt]; o = OFFS[pdir]
    def rate(x):
        a = arrows_at(idx, N, ks, srt, (x-r)*o); b = arrows_at(idx, N, ks, srt, (x+r)*o); v = angle(a, b)
        return np.where(np.isnan(v), 0.0, v)
    return rate(0)-0.5*(rate(-R)+rate(R)), rate(0)


def main():
    m = trimesh.load('out/orcfull_patch.stl', process=False); P = np.asarray(m.vertices, float); F = np.asarray(m.faces); h = float((P.max(0)-P.min(0)).max())
    cells = max(1, int(round(h/0.087))); t0 = time.time()
    r = voxel.run_direct(P, F, cells, None, height=h, flat=True, blend=True); idx, N = r['idx'], r['N']; size = float(r['size']); lo = r['lo']
    print(f'cubes {len(idx)}, {time.time()-t0:.0f} s')
    t0 = time.time(); change, crest, span = creases.fold_change(idx, N, 17.0); print(f'fold change {time.time()-t0:.0f} s')
    # the across direction as fold_change picks it: the neighbour with the largest arrow angle
    key = (idx[:, 0]*KB+idx[:, 1])*KB+idx[:, 2]; srt = np.argsort(key, kind='stable'); ks = key[srt]; pmax = np.zeros(len(idx)); pdir = np.zeros(len(idx), np.int32)
    for q, o in enumerate(OFFS):
        a = arrows_at(idx, N, ks, srt, np.repeat(o[None], len(idx), 0)); d = (N*a).sum(1); ang = np.degrees(np.arccos(np.clip(d, -1, 1))); ang[(d < 0) | np.isnan(d)] = 0
        upd = ang > pmax; pmax[upd] = ang[upd]; pdir[upd] = q
    fields = {'span': span}
    for R in (2, 3):
        t0 = time.time(); rise, rate0 = summit(idx, N, pdir, 1, R); fields[f'rise{R}'] = rise; print(f'rise R={R}: {time.time()-t0:.0f} s; cubes >= 17: {int((rise >= 17).sum())}, >= 10: {int((rise >= 10).sum())} (span >= 17: {int((span >= 17).sum())})')
    np.savez('out/summit.npz', idx=idx, lo=lo, size=size, **fields)
    cen = P[F].mean(1)
    views = [('span17', span, 17), ('rise2_17', fields['rise2'], 17), ('rise2_10', fields['rise2'], 10), ('rise3_17', fields['rise3'], 17), ('rise3_10', fields['rise3'], 10)]
    for spot, (focus, half) in SPOTS.items():
        focus = np.array(focus); near = (np.abs(cen-focus) < 2*half).all(1); Fp = F[near]; vv = m.face_normals[near].mean(0); vv /= np.linalg.norm(vv)
        C = lo+(idx+0.5)*size; box = (np.abs(C-focus) < 2*half).all(1)
        for name, f, cut in views:
            sel = box & (f >= cut)
            mesh = pv.PolyData(P, np.c_[np.full(len(Fp), 3), Fp].ravel()); pl = pv.Plotter(off_screen=True, window_size=(1600, 1200)); pl.add_mesh(mesh, color='lightgrey', smooth_shading=True)
            if sel.any():
                pts = pv.PolyData(C[sel]+0.5*size*N[sel]); pts['v'] = f[sel]
                pl.add_mesh(pts, scalars='v', cmap='plasma', clim=(cut, max(cut+1, float(np.percentile(f[f >= cut], 99)))), point_size=5 if half < 4 else 3, render_points_as_spheres=True, show_scalar_bar=False, reset_camera=False)
            pl.camera.focal_point = focus.tolist(); pl.camera.position = (focus+vv*30).tolist(); pl.camera.up = (0, 0, 1); pl.camera.parallel_projection = True; pl.camera.parallel_scale = half
            what = 'span: turn between the cubes one behind and one ahead' if name.startswith('span') else f'rise: rate here minus the mean rate {name[4]} cubes either side'
            pl.add_text(f'{what}; dots = cubes >= {cut} deg, dark = {cut}, yellow = strong; {int(sel.sum())} dots in view', font_size=11)
            pl.screenshot(f'out/summit_{spot}_{name}.png')
    print('saved out/summit_<spot>_<name>.png, spots', list(SPOTS), 'names', [v[0] for v in views])


if __name__ == '__main__': main()
