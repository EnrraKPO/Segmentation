"""Experiment: direction per dot from its neighbours, grouping by direction agreement.
For each line cell: neighbours within `reach` cells; the vectors to them are averaged as LINES (orientation tensor, so
left and right count the same); the main axis is the dot's direction and the anisotropy (how much stronger the main
axis is than the cross axis) says how definite it is. Dots with a definite direction group with neighbours whose
direction agrees within `agree` degrees; dots without one join the group of their directed neighbours (mass, not
direction). Output: close-up with dots coloured by group and their direction drawn.
usage: python exp_dirs.py mesh out.png [cells sens reach agree]"""
import sys, os, numpy as np, trimesh, pyvista as pv
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import voxel
path, out = sys.argv[1], sys.argv[2]
cells, sens = (int(sys.argv[3]) if len(sys.argv) > 3 else 1000), (float(sys.argv[4]) if len(sys.argv) > 4 else 6.0)
reach, agree = (float(sys.argv[5]) if len(sys.argv) > 5 else 2.0), (float(sys.argv[6]) if len(sys.argv) > 6 else 25.0)
m = trimesh.load(path, process=True); P, F = np.asarray(m.vertices), np.asarray(m.faces)
res = voxel.run(P, F, cells, sens, print)
li = np.nonzero(res['line'])[0]; C = (res['idx'][li]+0.5).astype(float); N = res['N'][li]; k = len(li)
tree = cKDTree(C); pairs = tree.query_pairs(reach, output_type='ndarray'); i, j = pairs[:, 0], pairs[:, 1]
v = C[j]-C[i]; v /= np.linalg.norm(v, axis=1, keepdims=True)
# orientation tensor per dot: sum of v v^T over its neighbours (both directions of each pair)
T = np.zeros((k, 3, 3))
for a, b in ((i, j), (j, i)):
    np.add.at(T, a, v[:, :, None]*v[:, None, :])
nn = np.bincount(i, minlength=k)+np.bincount(j, minlength=k)
# project into the tangent plane so the direction lies on the surface
w, e = np.linalg.eigh(T)                                           # ascending eigenvalues
d = e[:, :, -1]; d -= (d*N).sum(1, keepdims=True)*N; d /= np.linalg.norm(d, axis=1, keepdims=True)+1e-12
aniso = (w[:, -1]-w[:, -2])/(w[:, -1]+1e-12)                        # 1 = all neighbours on one axis, 0 = spread evenly
definite = (nn >= 2) & (aniso > 0.5)
print(f'{k} dots; with neighbours {(nn > 0).sum()}; definite direction {definite.sum()} ({definite.mean()*100:.0f}%)')
# grouping: neighbouring definite dots whose directions agree
cos_ok = np.abs((d[i]*d[j]).sum(1)) >= np.cos(np.radians(agree))
link = definite[i] & definite[j] & cos_ok
_, grp = connected_components(coo_matrix((np.ones(link.sum()), (i[link], j[link])), shape=(k, k)), directed=False)
grp = np.where(definite, grp, -1)
# undirected dots join the largest group among their neighbours
for _ in range(3):
    und = np.nonzero(grp < 0)[0]
    if len(und) == 0: break
    for a in und:
        nb = tree.query_ball_point(C[a], reach); g = grp[nb]; g = g[g >= 0]
        if len(g): grp[a] = np.bincount(g).argmax()
sizes = np.bincount(grp[grp >= 0]); print(f'groups {len(sizes)}; >= 20 dots: {(sizes >= 20).sum()}; largest {sizes.max()}; dots in groups >= 20: {sizes[sizes >= 20].sum()}')
np.savez(os.path.splitext(out)[0]+'.npz', C=res['lo']+C*res['size'], d=d, grp=grp, aniso=aniso, definite=definite, size=res['size'])
# picture: waist close-up, dots coloured by group, direction as a short segment
W = res['lo']+C*res['size']; rng = np.random.default_rng(0); cols = rng.uniform(0.1, 0.95, (grp.max()+2, 3)); col = cols[np.where(grp >= 0, grp, grp.max()+1)]; col[grp < 0] = [0, 0, 0]
mesh = pv.PolyData(P, np.hstack([np.full((len(F), 1), 3), F]).ravel())
seg = np.stack([W-0.6*res['size']*d, W+0.6*res['size']*d], 1)
pl = pv.Plotter(off_screen=True, shape=(1, 2), window_size=(3000, 1500))
for s in (0, 1):
    pl.subplot(0, s); pl.add_mesh(mesh, color='lightgrey', smooth_shading=True)
    if s == 0:
        pd = pv.PolyData(W); pd['c'] = col; pl.add_mesh(pd, scalars='c', rgb=True, point_size=7, render_points_as_spheres=True); pl.add_text('dots by group (black = no direction)', font_size=10)
    else:
        ln = pv.PolyData(seg.reshape(-1, 3), lines=np.c_[np.full(k, 2), np.arange(0, 2*k, 2), np.arange(1, 2*k, 2)].ravel()); ln['c'] = np.repeat(col, 2, 0)
        pl.add_mesh(ln, scalars='c', rgb=True, line_width=3); pl.add_text('direction per dot', font_size=10)
    focus = np.array([12, 135, 12.]); pl.camera_position = 'zy'; pl.camera.focal_point = focus.tolist(); pl.camera.position = (focus+(np.array(pl.camera.position)-np.array(pl.camera.focal_point))).tolist(); pl.camera.zoom(6)
pl.screenshot(out); print('saved', out)
