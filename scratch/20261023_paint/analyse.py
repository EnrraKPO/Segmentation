import sys, os, pickle; root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); sys.path.insert(0, root)
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import shortest_path
res = pickle.load(open(os.path.join(root, 'out', 'paint_walk.pkl'), 'rb')); emb = res['emb']; con = res['con']
fid = emb['fid']; Q = emb['nodes']; nf = emb['nface']; P = emb['P']; F = emb['F']; nbr = emb['nbr']
jump = np.linalg.norm(Q[fid+1]-Q[fid], axis=1); print('failed stretches', len(fid), 'projected-point distance p50/p90/max', np.round(np.percentile(jump, [50, 90, 100]), 3))
# face-graph distance between fa and fb for a sample
m = len(F); I = np.repeat(np.arange(m), 3); J = nbr.ravel(); ok = J >= 0
G = coo_matrix((np.ones(ok.sum()), (I[ok], J[ok])), shape=(m, m)).tocsr()
samp = fid[:: max(1, len(fid)//60)]
d = []
for i in samp:
    dist = shortest_path(G, indices=nf[i], unweighted=True, limit=20); d.append(dist[nf[i+1]])
d = np.array(d); print('face-graph hops fa->fb (sample):', np.round(np.sort(d), 0))
# tiny regions: how many pieces, and are they at line ends?
A = con['areas']; tiny = np.where(A < 0.01)[0]; region = con['region']
cnt = np.bincount(region[region >= 0], minlength=con['nreg']); print('tiny regions', len(tiny), 'pieces per tiny region p50/p90/max', np.percentile(cnt[tiny], [50, 90, 100]))
# how far is a tiny region's piece from the nearest line END
start = emb['start']; ends = np.concatenate([Q[start[:-1]], Q[start[1:]-1]])
from scipy.spatial import cKDTree
t = cKDTree(ends)
DV, DF, pid = res['DV'], res['DF'], res['pid']; cen = DV[DF].mean(1); rp = region[pid]
istiny = np.isin(rp, tiny); dd, _ = t.query(cen[istiny]); print('tiny-region display triangles', istiny.sum(), 'distance to nearest line end p50/p90', np.round(np.percentile(dd, [50, 90]), 3))
