"""Definitive image-space check: shaded render; darkest pixel per row (red); each trail file drawn in its own colour in
the SAME plotter and camera as the mesh. Prints per-trail median offset from the darkest pixel, in um (2.25 um/px).
usage: python verdict.py mesh out.png chains1.npz [chains2.npz ...]
env: SNAP=1 moves each node to the nearest surface point (+0.01 mm) and reports how far nodes sat off the surface. FOCUS=x,y,z (default Harley waist 12,135,12) VIEW=zy|xz|yz|xy ZOOM=25 LIFT=0.03 (judge with LIFT=0)"""
import sys, os, numpy as np, pyvista as pv, trimesh, imageio.v2 as iio
from scipy.spatial import cKDTree
M, outp, files = sys.argv[1], sys.argv[2], sys.argv[3:]
m = trimesh.load(M, process=True); P, F = np.asarray(m.vertices), np.asarray(m.faces); VN = m.vertex_normals; kt = cKDTree(P)
SNAP = os.environ.get('SNAP', '0') == '1'
def snap(Q):
    """nodes onto the surface: nearest mesh point, then 0.01 mm out along the face normal"""
    C, dist, tri = trimesh.proximity.closest_point(m, Q)
    print(f'   snap: nodes sat {np.median(dist)*1000:.0f} um (median) / {np.percentile(dist, 90)*1000:.0f} um (p90) off the surface')
    return C+0.01*m.face_normals[tri]
mesh = pv.PolyData(P, np.hstack([np.full((len(F), 1), 3), F]).ravel()); W = 1500; UMPX = np.linalg.norm(P.max(0)-P.min(0))*1000/W/float(os.environ.get('ZOOM', 25)); focus = np.array([float(v) for v in os.environ.get('FOCUS', '12,135,12').split(',')]); zoom = float(os.environ.get('ZOOM', 25)); VIEW = os.environ.get('VIEW', 'zy')
def cam(pl):
    pl.camera_position = VIEW; pl.camera.focal_point = focus.tolist(); pl.camera.position = (focus+(np.array(pl.camera.position)-np.array(pl.camera.focal_point))).tolist(); pl.camera.zoom(zoom)
pl = pv.Plotter(off_screen=True, window_size=(W, W)); pl.add_mesh(mesh, color='lightgrey', smooth_shading=True); cam(pl); shade = pl.screenshot(return_img=True).astype(float)
dark = np.array([300+np.argmin(shade[r, 300:1300, 0]) for r in range(W)])
cols = [(0, 200, 255), (255, 140, 0), (0, 255, 0), (255, 0, 255)]; out = shade.copy()
for fi, f in enumerate(files):
    z = np.load(f); pts, n = z['pts'], z['n']; chains = np.split(np.arange(len(pts)), np.cumsum(n)[:-1])
    for ch in chains:
        if np.linalg.norm(pts[ch]-focus, axis=1).min() > 100/zoom or len(ch) < 40: continue
        LIFT = float(__import__("os").environ.get("LIFT", 0.03)); Q = snap(pts[ch]) if SNAP else pts[ch]+LIFT*VN[kt.query(pts[ch])[1]]
        p = pv.Plotter(off_screen=True, window_size=(W, W)); p.add_mesh(mesh, color='white', lighting=False); cam(p)
        p.add_mesh(pv.PolyData(Q, lines=np.r_[len(ch), np.arange(len(ch))]), color='black', line_width=3, lighting=False, reset_camera=False)
        mask = p.screenshot(return_img=True)[:, :, 0] < 128; out[mask] = cols[fi % 4]
        offs = [(np.nonzero(mask[r])[0].mean()-dark[r])*UMPX for r in range(100, 1400, 25) if mask[r].any()]
        if offs: print(f'{f}: trail {len(ch):4d} nodes, offset from darkest pixel: median {np.median(offs):+5.0f} um, spread {np.std(offs):4.0f} um  (colour {cols[fi % 4]})')
for r in range(0, W, 5): out[max(r-1, 0):r+2, dark[r]-2:dark[r]+3] = (255, 0, 0)
iio.imwrite(outp, out.astype(np.uint8)); print('saved', outp)
