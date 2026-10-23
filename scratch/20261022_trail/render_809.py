import sys, os, numpy as np, trimesh, pyvista as pv
os.chdir(r'D:\Stuff\PrintSuit\SegmentationTool'); from PIL import Image, ImageDraw
I = int(sys.argv[1]); z = np.load(f'out/refine{I}.npz'); coarse = z['coarse']
m = trimesh.load('out/orcfull_patch.stl', process=False); P = np.asarray(m.vertices, float); F = np.asarray(m.faces); cen = P[F].mean(1)
focus = coarse.mean(0); half = 0.6*np.linalg.norm(coarse-focus, axis=1).max()+0.5; near = (np.abs(cen-focus) < 2*half).all(1); Fp = F[near]; vv = m.face_normals[near].mean(0); vv /= np.linalg.norm(vv)
lift = lambda Q: (lambda r: r[0]+0.01*m.face_normals[r[2]])(trimesh.proximity.closest_point(m, Q))
for name in ('saved', 'walk'):
    L = np.split(z[f'{name}_pts'], np.cumsum(z[f'{name}_n'])[:-1])
    pl = pv.Plotter(off_screen=True, window_size=(1600, 1200)); pl.add_mesh(pv.PolyData(P, np.c_[np.full(len(Fp), 3), Fp].ravel()), color='lightgrey', smooth_shading=True)
    cl = lift(coarse); pl.add_mesh(pv.PolyData(cl, lines=np.r_[len(cl), np.arange(len(cl))]), color='orange', line_width=2, reset_camera=False)
    for l in L:
        ll = lift(l); pl.add_mesh(pv.PolyData(ll, lines=np.r_[len(ll), np.arange(len(ll))]), color='blue', line_width=3, reset_camera=False); pl.add_mesh(pv.PolyData(ll[[0, -1]]), color='black', point_size=8, render_points_as_spheres=True, reset_camera=False)
    pl.camera.focal_point = focus.tolist(); pl.camera.position = (focus+vv*30).tolist(); pl.camera.up = (0, 0, 1); pl.camera.parallel_projection = True; pl.camera.parallel_scale = half
    pl.add_text({'saved': 'saved run: fine walk = trail, arc trend', 'walk': 'fine walk = anchored walk, chord trend (committed)'}[name]+'; coarse line orange, fine lines blue', font_size=11); pl.screenshot(f'out/refine{I}_{name}.png'); pl.close()
ims = [Image.open(f'out/refine{I}_{n}.png') for n in ('saved', 'walk')]; w, hh = ims[0].size; out = Image.new('RGB', (2*w+20, hh+60), 'white'); dr = ImageDraw.Draw(out)
for i, (im, t) in enumerate(zip(ims, ('LEFT: saved run (trail, arc)', 'RIGHT: anchored walk, chord'))): out.paste(im, (i*(w+20), 60)); dr.text((i*(w+20)+20, 15), t, fill='black')
out.resize((out.width//2, out.height//2)).save(f'out/refine{I}_compare.png'); print('saved', f'out/refine{I}_compare.png')
