import sys, os, pickle, time; root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); sys.path.insert(0, root)
import numpy as np
res = pickle.load(open(os.path.join(root, 'out', 'paint_walk.pkl'), 'rb')); emb = res['emb']
fid = np.unique(emb['fid']); Q = emb['nodes']; nf = emb['nface']; P = emb['P']; F = emb['F']
jump = np.linalg.norm(Q[fid+1]-Q[fid], axis=1); print('failed stretches', len(fid), 'surface-point distance p10/p50/p90/max', np.round(np.percentile(jump, [10, 50, 90, 100]), 4))
# edge lengths of the faces involved
e = np.linalg.norm(P[F[:, 1]]-P[F[:, 0]], axis=1); print('edge of fa p10/p50/p90', np.round(np.percentile(e[nf[fid]], [10, 50, 90]), 4), 'mesh p10', np.round(np.percentile(e, [1, 10]), 4))
# normals of fa and fb: facing away (two sides of a sheet)?
N = np.cross(P[F[:, 1]]-P[F[:, 0]], P[F[:, 2]]-P[F[:, 0]]); N /= np.linalg.norm(N, axis=1)[:, None]
dot = (N[nf[fid]]*N[nf[fid+1]]).sum(1); print('normal agreement fa.fb: fraction < 0 (opposite sides)', np.round((dot < 0).mean(), 3), 'fraction > 0.9', np.round((dot > 0.9).mean(), 3))
# VTK locator timing for projection
import pyvista as pv, vtk
mesh = pv.PolyData(P, np.hstack([np.full((len(F), 1), 3), F]).ravel())
t0 = time.time(); loc = vtk.vtkStaticCellLocator(); loc.SetDataSet(mesh); loc.BuildLocator(); print('locator built', round(time.time()-t0, 2), 's')
pts = Q[:50000]; out = np.zeros((len(pts), 3)); cid = vtk.reference(0); sub = vtk.reference(0); d2 = vtk.reference(0.0); cp = [0.0, 0.0, 0.0]; fids = np.zeros(len(pts), np.int64)
t0 = time.time()
for i, p in enumerate(pts):
    loc.FindClosestPoint(p, cp, cid, sub, d2); out[i] = cp; fids[i] = cid.get()
print('50k points projected by vtk in', round(time.time()-t0, 2), 's; agreement with earlier faces', (fids == nf[:50000]).mean(), 'max point diff', np.abs(out-Q[:50000]).max())
