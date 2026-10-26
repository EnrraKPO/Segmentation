import sys, os, pickle, time; root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); sys.path.insert(0, root)
import numpy as np, regions
res = pickle.load(open(os.path.join(root, 'out', 'paint_walk.pkl'), 'rb')); emb = res['emb']
fid = np.unique(emb['fid']); Q = emb['nodes']; nf = emb['nface']; P = emb['P']; F = emb['F']; node2 = emb['node2']
gap = np.linalg.norm(Q[fid+1]-Q[fid], axis=1); print('failed', len(fid), 'gap p10/p50/p90/max', np.round(np.percentile(gap, [10, 50, 90, 100]), 4))
N = np.cross(P[F[:, 1]]-P[F[:, 0]], P[F[:, 2]]-P[F[:, 0]]); N /= np.linalg.norm(N, axis=1)[:, None]
dot = (N[nf[fid]]*N[nf[fid+1]]).sum(1); print('normal dot < 0:', (dot < 0).sum(), ' < 0.5:', (dot < 0.5).sum(), ' > 0.9:', (dot > 0.9).sum())
shared = np.array([len(set(F[nf[i]]) & set(F[nf[i+1]])) for i in fid]); print('faces sharing 2 verts (adjacent):', (shared == 2).sum(), '1 vert:', (shared == 1).sum(), 'none:', (shared == 0).sum())
# re-walk the failing pairs alone: with pruning, without, deeper
start = np.arange(0, 2*len(fid)+1, 2).astype(np.int64)
Qp = np.empty((2*len(fid), 3)); Qp[0::2] = Q[fid]; Qp[1::2] = Q[fid+1]; nfp = np.empty(2*len(fid), np.int64); nfp[0::2] = nf[fid]; nfp[1::2] = nf[fid+1]
n2 = np.empty((2*len(fid), 2)); n2[0::2] = node2[fid]; n2[1::2] = node2[fid+1]
for depth, prune in [(8, 1), (8, 0), (16, 0), (24, 1)]:
    t0 = time.time(); out = regions._walk(P, F, Qp, nfp, n2, start, emb['C2'], emb['nbr'], emb['fedge'], emb['E'], depth, 100000, prune); print(f'depth {depth} prune {prune}: {len(out[9])} still fail, {time.time()-t0:.1f} s')
