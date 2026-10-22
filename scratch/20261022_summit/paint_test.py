"""Summit dots by painting outward (the user's rule, 2026-10-22). From every cube, walk outward over the surface in M
directions of its tangent plane, cube by cube. At each new cube note the change of arrow from the cube before; a change
of (nearly) zero is dismissed. The trend is the mean direction of the changes so far; a change pointing against the trend
(more than AGAINST degrees away from it) means the next slope has begun: the walk stops there, that cube excluded.
Region = the centre plus the end cubes of the M walks. Bend = the largest angle between any two arrows of the region.
The cube is a dot when bend >= cut and its own arrow sits about halfway between the two extremes
(|angle to one - angle to the other| <= HALF * bend). Handles: M, AGAINST (deg), ZERO (deg), MAXSTEP (cubes), HALF, cut.
Run on a box around each spot (whole orc too slow for a first look). Renders out/paint_<spot>_<name>.png and
out/paint_<spot>_compare.png (left: today's span dots at 17; right: paint dots at 17)."""
import sys, os, time, numpy as np, trimesh, pyvista as pv
from scipy.spatial import cKDTree
sys.path.insert(0, r'D:\Stuff\PrintSuit\SegmentationTool'); os.chdir(r'D:\Stuff\PrintSuit\SegmentationTool'); import voxel, creases

H = dict(M=8, AGAINST=45.0, ZERO=0.5, MAXSTEP=20, HALF=0.5)
for k in H: H[k] = type(H[k])(os.environ.get(k, H[k]))
SPOTS = {'rim': ([219.51, 167.43, 165.39], 1.2, 6.0), 'shield': ([162.13, 137.7, 107.92], 2.5, 6.0), 'ring': ([194.72, 123.34, 169.49], 2.5, 6.0)}   # focus, zoom half-view, box half (mm)


def paint(C, N, size, log=print):
    """per cube: bend (deg) and halfness (0 = the cube's arrow is one of the extremes, 1 = exactly halfway)"""
    n = len(C); tree = cKDTree(C); M = H['M']; cos_against = np.cos(np.radians(H['AGAINST'])); zero = 2*np.sin(np.radians(H['ZERO'])/2)
    a = np.where((np.abs(N[:, 0]) < 0.9)[:, None], np.array([1., 0, 0]), np.array([0, 1., 0]))
    t1 = np.cross(N, a); t1 /= np.linalg.norm(t1, axis=1, keepdims=True)+1e-12; t2 = np.cross(N, t1)
    th = np.arange(M)*2*np.pi/M
    rows = np.repeat(np.arange(n), M); d = np.cos(th)[None, :, None]*t1[:, None]+np.sin(th)[None, :, None]*t2[:, None]; d = d.reshape(-1, 3)
    cur = rows.copy(); p = C[cur].copy(); nprev = N[cur].copy(); trend = np.zeros((n*M, 3)); alive = np.ones(n*M, bool); end = rows.copy()   # end: the last accepted cube
    t0 = time.time()
    for step in range(H['MAXSTEP']):
        idx = np.nonzero(alive)[0]
        if len(idx) == 0: break
        q = p[idx]+size*d[idx]; dd, jj = tree.query(q, k=3)
        pick = jj[:, 0]; same = pick == cur[idx]; pick[same] = jj[same, 1]                  # never the cube we stand on
        far = dd[np.arange(len(idx)), np.where(same, 1, 0)] > 1.5*size                     # left the surface (box edge, hole)
        nn = N[pick]; delta = nn-nprev[idx]; mag = np.linalg.norm(delta, axis=1); changed = mag > zero
        u = delta/np.maximum(mag, 1e-12)[:, None]; tr = trend[idx]; has_tr = np.linalg.norm(tr, axis=1) > 1e-9
        tn = tr/np.maximum(np.linalg.norm(tr, axis=1), 1e-12)[:, None]; against = changed & has_tr & ((u*tn).sum(1) < cos_against)
        stop = far | against
        go = idx[~stop]; pk = pick[~stop]
        trend[go[changed[~stop]]] += u[~stop][changed[~stop]]
        cur[go] = pk; end[go] = pk; p[go] = C[pk]; nprev[go] = N[pk]
        dg = d[go]-((d[go]*N[pk]).sum(1))[:, None]*N[pk]; d[go] = dg/np.maximum(np.linalg.norm(dg, axis=1), 1e-12)[:, None]   # the direction follows the surface
        alive[idx[stop]] = False
    if log: log(f'paint: {n} cubes x {M} directions, {step+1} steps, {time.time()-t0:.0f} s; walks still going at the end: {int(alive.sum())}')
    E = N[end].reshape(n, M, 3); R = np.concatenate([N[:, None], E], 1)                       # the region's arrows: centre + M ends
    G = np.einsum('nik,njk->nij', R, R); ang = np.degrees(np.arccos(np.clip(G, -1, 1)))
    flat = ang.reshape(n, -1); best = flat.argmax(1); bend = flat[np.arange(n), best]; i, j = best // (M+1), best % (M+1)
    ai = ang[np.arange(n), 0, i]; aj = ang[np.arange(n), 0, j]; half = 1-np.abs(ai-aj)/np.maximum(bend, 1e-9)
    return bend, half


def render(P, F, m, C, N, size, f, sel, focus, half, title, path):
    cen = P[F].mean(1); near = (np.abs(cen-focus) < 2*half).all(1); Fp = F[near]; vv = m.face_normals[near].mean(0); vv /= np.linalg.norm(vv)
    mesh = pv.PolyData(P, np.c_[np.full(len(Fp), 3), Fp].ravel()); pl = pv.Plotter(off_screen=True, window_size=(1600, 1200)); pl.add_mesh(mesh, color='lightgrey', smooth_shading=True)
    box = (np.abs(C-focus) < 2*half).all(1); s = sel & box
    if s.any():
        pts = pv.PolyData(C[s]+0.5*size*N[s]); pts['v'] = f[s]; lo = float(f[s].min())
        pl.add_mesh(pts, scalars='v', cmap='plasma', clim=(lo, max(lo+1, float(np.percentile(f[s], 99)))), point_size=5 if half < 4 else 3, render_points_as_spheres=True, show_scalar_bar=False, reset_camera=False)
    pl.camera.focal_point = focus.tolist(); pl.camera.position = (focus+vv*30).tolist(); pl.camera.up = (0, 0, 1); pl.camera.parallel_projection = True; pl.camera.parallel_scale = half
    pl.add_text(f'{title}; {int(s.sum())} dots in view', font_size=11); pl.screenshot(path)


def main():
    from PIL import Image, ImageDraw
    m = trimesh.load('out/orcfull_patch.stl', process=False); P = np.asarray(m.vertices, float); F = np.asarray(m.faces); h = float((P.max(0)-P.min(0)).max())
    cube = 0.087; cells = max(1, int(round(h/cube))); cut = float(os.environ.get('CUT', 17))
    for spot, (focus, zoom, bh) in SPOTS.items():
        focus = np.array(focus); box = (focus-bh, focus+bh)
        r = voxel.run_direct(P, F, cells, None, height=h, flat=True, box=box, blend=True); idx, N = r['idx'], r['N']; size = float(r['size']); C = r['lo']+(idx+0.5)*size
        change, crest, span = creases.fold_change(idx, N, cut)
        bend, half = paint(C, N, size); dot = (bend >= cut) & (half >= H['HALF'])
        print(f'{spot}: {len(C)} cubes in the box; span >= {cut:g}: {int((span >= cut).sum())}; paint dots: {int(dot.sum())} (bend >= {cut:g}: {int((bend >= cut).sum())}, of which halfway: {int(dot.sum())})')
        np.savez(f'out/paint_{spot}.npz', C=C, N=N, size=size, span=span, bend=bend, half=half, handles=str(H))
        for name, zh in (('zoom', zoom), ('wide', bh)):
            render(P, F, m, C, N, size, span, span >= cut, focus, zh, f"today's dots: span >= {cut:g}", f'out/paint_{spot}_{name}_span.png')
            render(P, F, m, C, N, size, bend, dot, focus, zh, f'paint dots: bend >= {cut:g} and halfway (HALF {H["HALF"]}, AGAINST {H["AGAINST"]:g} deg, MAXSTEP {H["MAXSTEP"]})', f'out/paint_{spot}_{name}_paint.png')
            a = Image.open(f'out/paint_{spot}_{name}_span.png'); b = Image.open(f'out/paint_{spot}_{name}_paint.png'); w, hh = a.size
            out = Image.new('RGB', (2*w+20, hh+60), 'white'); out.paste(a, (0, 60)); out.paste(b, (w+20, 60)); dr = ImageDraw.Draw(out)
            dr.text((20, 15), "LEFT: today's dots (span, cut %g)" % cut, fill='black'); dr.text((w+40, 15), 'RIGHT: paint dots (bend >= %g, halfway)' % cut, fill='black')
            out.resize((out.width//2, out.height//2)).save(f'out/paint_{spot}_{name}_compare.png')
    print('saved out/paint_<spot>_<zoom|wide>_compare.png for', list(SPOTS))


if __name__ == '__main__': main()
