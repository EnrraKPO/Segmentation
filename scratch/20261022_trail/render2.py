"""render out/trail_<spot><tag>.npz as out/trail_<spot><tag>_compare.png: left today, right trail. usage: render2.py spot tag label_today label_trail"""
import sys, os, numpy as np, trimesh, pyvista as pv
sys.path.insert(0, r'D:\Stuff\PrintSuit\SegmentationTool'); os.chdir(r'D:\Stuff\PrintSuit\SegmentationTool'); sys.path.insert(0, 'scratch/20261022_smooth'); import smooth_test as T
from PIL import Image, ImageDraw
spot, tag = sys.argv[1], sys.argv[2]; labels = dict(zip(('today', 'trail'), sys.argv[3:5])); focus, half = T.SPOTS[spot]; focus = np.array(focus)
z = np.load(T.DOTS, allow_pickle=True); lo = z['lo']; size = float(z['size']); cut = float(z['cut']); idx, N, turn = z['idx'], z['N'], z['turn']; C = lo+(idx+0.5)*size
m = trimesh.load('out/orcfull_patch.stl', process=False); P = np.asarray(m.vertices, float); F = np.asarray(m.faces); cen = P[F].mean(1)
lift = lambda Q: (lambda r: r[0]+0.01*m.face_normals[r[2]])(trimesh.proximity.closest_point(m, Q))
r = np.load(f'out/trail_{spot}{tag}.npz'); sel = (np.abs(C-focus) < T.BOX).all(1)
near = (np.abs(cen-focus) < 2*half).all(1); Fp = F[near]; vv = m.face_normals[near].mean(0); vv /= np.linalg.norm(vv)
for name in ('today', 'trail'):
    L = np.split(r[f'{name}_pts'], np.cumsum(r[f'{name}_n'])[:-1])
    mesh = pv.PolyData(P, np.c_[np.full(len(Fp), 3), Fp].ravel()); pl = pv.Plotter(off_screen=True, window_size=(1600, 1200)); pl.add_mesh(mesh, color='lightgrey', smooth_shading=True)
    inv = sel & (np.abs(C-focus) < 2*half).all(1); pts = pv.PolyData(C[inv]+0.5*size*N[inv]); pts['t'] = turn[inv]
    pl.add_mesh(pts, scalars='t', cmap='plasma', clim=(cut, float(np.percentile(turn, 99))), point_size=4, render_points_as_spheres=True, show_scalar_bar=False, reset_camera=False)
    for l in L:
        if np.linalg.norm(l-focus, axis=1).min() < 2*half: ll = lift(l); pl.add_mesh(pv.PolyData(ll, lines=np.r_[len(ll), np.arange(len(ll))]), color='blue', line_width=3, reset_camera=False); pl.add_mesh(pv.PolyData(ll[[0, -1]]), color='black', point_size=8, render_points_as_spheres=True, reset_camera=False)
    pl.camera.focal_point = focus.tolist(); pl.camera.position = (focus+vv*30).tolist(); pl.camera.up = (0, 0, 1); pl.camera.parallel_projection = True; pl.camera.parallel_scale = half
    pl.add_text(f'{labels[name]}; dots by strength, lines blue, ends black', font_size=11); pl.screenshot(f'out/trail_{spot}{tag}_{name}.png'); pl.close()
ims = [Image.open(f'out/trail_{spot}{tag}_{n}.png') for n in ('today', 'trail')]; w, hh = ims[0].size; out = Image.new('RGB', (2*w+20, hh+60), 'white'); dr = ImageDraw.Draw(out)
for i, (im, n) in enumerate(zip(ims, ('today', 'trail'))): out.paste(im, (i*(w+20), 60)); dr.text((i*(w+20)+20, 15), ('LEFT: ', 'RIGHT: ')[i]+labels[n], fill='black')
out.resize((out.width//2, out.height//2)).save(f'out/trail_{spot}{tag}_compare.png'); print('saved', f'out/trail_{spot}{tag}_compare.png')
