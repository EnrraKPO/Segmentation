"""Render a view with a pixel grid, and/or turn pixels into 3D points on the surface.
usage: python pick.py mesh out.png [u,v u,v ...]     env: VIEW=xz AZ=0 EL=0 ZOOM=1 FOCUS=x,y,z (optional)
Prints for each pixel the 3D point and its distance to the nearest mesh vertex (sanity: should be well under 1 mm)."""
import sys, os, numpy as np, pyvista as pv, trimesh
from scipy.spatial import cKDTree
M, outp = sys.argv[1:3]; pix = [tuple(int(x) for x in a.split(',')) for a in sys.argv[3:]]
VIEW = os.environ.get('VIEW', 'xz'); AZ = float(os.environ.get('AZ', 0)); EL = float(os.environ.get('EL', 0)); zoom = float(os.environ.get('ZOOM', 1)); W = 1600
m = trimesh.load(M, process=True); P, F = np.asarray(m.vertices), np.asarray(m.faces)
mesh = pv.PolyData(P, np.hstack([np.full((len(F), 1), 3), F]).ravel())
pl = pv.Plotter(off_screen=True, window_size=(W, W)); pl.add_mesh(mesh, color='lightgrey', smooth_shading=True)
pl.camera_position = VIEW
if 'FOCUS' in os.environ:
    focus = np.array([float(v) for v in os.environ['FOCUS'].split(',')]); pl.camera.focal_point = focus.tolist(); pl.camera.position = (focus+(np.array(pl.camera.position)-np.array(pl.camera.focal_point))).tolist()
if AZ: pl.camera.Azimuth(AZ)
if EL: pl.camera.Elevation(EL)
pl.camera.zoom(zoom)
img = pl.screenshot(return_img=True); depth = pl.get_image_depth(fill_value=np.nan)
cam = pl.camera; pos = np.array(cam.position); fp = np.array(cam.focal_point); up = np.array(cam.up)
fwd = fp-pos; fwd /= np.linalg.norm(fwd); right = np.cross(fwd, up); right /= np.linalg.norm(right); up = np.cross(right, fwd)
half = np.tan(np.radians(cam.view_angle/2))
kt = cKDTree(P)
for (u, v) in pix:
    d = -depth[v, u]                                   # distance along the view direction
    x = (u-W/2)/(W/2)*half; y = (W/2-v)/(W/2)*half
    p = pos+(fwd+x*right+y*up)*d
    print(f'pixel ({u},{v}) -> {p.round(2)}  nearest vertex {kt.query(p)[0]*1000:.0f} um')
# grid every 100 px for reading coordinates off the image
img = img.copy(); img[::100, :] = (255, 0, 0); img[:, ::100] = (255, 0, 0)
import imageio.v2 as iio; iio.imwrite(outp, img); print('saved', outp)
