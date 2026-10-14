"""Ground truth off the mesh, no grid: walk a straight path across a feature, project each step onto the surface,
record the surface normal. Plot the normal's angle (vs the start) against distance, and its change per 0.35 mm.
env GRID=grid.npz adds the grid turn (max within 0.5 mm) for comparison.
usage: python mesh_profile.py mesh out.png "label x,y,z dx,dy,dz" ...   (dir = walking direction, +-2 mm, 0.05 mm steps)"""
import sys, os, numpy as np, trimesh, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
m = trimesh.load(sys.argv[1], process=True); outp = sys.argv[2]; specs = sys.argv[3:]
G = np.load(os.environ['GRID']) if 'GRID' in os.environ else None
if G is not None:
    from scipy.spatial import cKDTree; GX = G['lo']+(G['idx']+0.5)*G['size']; gt = cKDTree(GX)
fig, ax = plt.subplots(len(specs), 1, figsize=(8, 2.6*len(specs)), squeeze=False)
for i, s in enumerate(specs):
    label, p, d = s.split(); p = np.array([float(v) for v in p.split(',')]); d = np.array([float(v) for v in d.split(',')]); d /= np.linalg.norm(d)
    t = np.arange(-2, 2.0001, 0.05); Q = p+t[:, None]*d
    C, dist, tri = trimesh.proximity.closest_point(m, Q)
    b = trimesh.triangles.points_to_barycentric(m.triangles[tri], C); Nn = (m.vertex_normals[m.faces[tri]]*b[:, :, None]).sum(1); Nn /= np.linalg.norm(Nn, axis=1, keepdims=True)   # smooth normals, as the render shades
    # re-walk on the surface so the path does not leave it: use projected points, cumulative arc length
    s_ = np.r_[0, np.cumsum(np.linalg.norm(np.diff(C, axis=0), axis=1))]; s_ -= s_[len(t)//2]
    ang = np.degrees(np.arccos(np.clip(Nn @ Nn[len(t)//2], -1, 1)))
    step = np.degrees(np.arccos(np.clip((Nn[7:]*Nn[:-7]).sum(1), -1, 1)))     # turn over 7 steps = 0.35 mm = one cube
    ax[i, 0].plot(s_, ang, 'k', label='normal angle vs centre (deg)'); ax[i, 0].plot(s_[3:-4], step, 'r', label='turn per 0.35 mm (deg)')
    ax[i, 0].axhline(6, color='r', ls=':', lw=0.8); ax[i, 0].set_title(label); ax[i, 0].set_xlabel('mm along the surface'); ax[i, 0].legend(fontsize=8); ax[i, 0].grid(alpha=0.3)
    msg = f'{label}: TRUE max turn per cube {step.max():.1f} deg; total swing {ang.max():.1f} deg'
    if G is not None: nb = gt.query_ball_point(p, 0.5); msg += f'; GRID turn within 0.5 mm: max {G["turn"][nb].max():.1f}, median {np.median(G["turn"][nb]):.1f}, line cells {int(G["line"][nb].sum())}/{len(nb)}'
    print(msg)
plt.tight_layout(); plt.savefig(outp, dpi=110); print('saved', outp)
