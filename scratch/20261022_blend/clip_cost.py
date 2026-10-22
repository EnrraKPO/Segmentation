"""how long the clipping maths itself takes per pair, on one core, on the tube box's multi-triangle pairs"""
import sys, os, time, numpy as np, trimesh
sys.path.insert(0, r'D:\Stuff\PrintSuit\SegmentationTool'); os.chdir(r'D:\Stuff\PrintSuit\SegmentationTool'); import voxel
m = trimesh.load('out/orcfull_patch.stl', process=False); P = np.asarray(m.vertices, float); F = np.asarray(m.faces); h = float((P.max(0)-P.min(0)).max())
d = np.load('out/debug/20261020_144432/lines.npz'); n = d['n']; pts = d['pts']; off = np.r_[0, np.cumsum(n)]; L = pts[off[318]:off[319]]
size = 0.087/3; box = (L.min(0)-3, L.max(0)+3)      # a bigger box: more pairs
r = voxel.run_direct(P, F, h/size, None, height=h, flat=True, box=box, workers=1, blend=False)
# rebuild the pairs the way run_direct does, single process
A, B, C = P[F[:, 0]], P[F[:, 1]], P[F[:, 2]]; fn = np.cross(B-A, C-A); fn /= np.linalg.norm(fn, axis=1, keepdims=True)
voxel._W.update(A=A, B=B, C=C, FN=fn, lo=r['lo'], size=float(r['size']))
# pairs: every cube-triangle hit in the box (reuse the task machinery)
lo = r['lo']; sz = float(r['size']); T = np.stack([A, B, C], 1); g = (T-lo)/sz; imin = np.floor(g.min(1)).astype(np.int64); imax = np.floor(g.max(1)).astype(np.int64)
bmin = np.floor((box[0]-lo)/sz).astype(np.int64); bmax = np.floor((box[1]-lo)/sz).astype(np.int64); imin = np.maximum(imin, bmin); imax = np.minimum(imax, bmax); keep = (imax >= imin).all(1)
tri = np.nonzero(keep)[0]; keys = []; reps = []
for t in tri:
    sp = imax[t]-imin[t]+1; offs = np.stack(np.meshgrid(*[np.arange(s) for s in sp], indexing='ij'), -1).reshape(-1, 3); cid = imin[t]+offs
    Q = lo+(cid+0.5)*sz; ok = voxel._tri_box_overlap(Q, np.repeat(A[t][None], len(Q), 0), np.repeat(B[t][None], len(Q), 0), np.repeat(C[t][None], len(Q), 0), sz*0.5)
    k = (cid[ok][:, 0]*voxel.KB+cid[ok][:, 1])*voxel.KB+cid[ok][:, 2]; keys.append(k); reps.append(np.full(len(k), t))
pk = np.concatenate(keys); pr = np.concatenate(reps); o = np.argsort(pk, kind='stable'); pk, pr = pk[o], pr[o]
first = np.r_[True, pk[1:] != pk[:-1]]; grp = np.cumsum(first)-1; multi = np.bincount(grp)[grp] >= 2; pk, pr = pk[multi], pr[multi]
t0 = time.time(); ub, tot = voxel._blend_task((pk, pr.astype(np.int32))); t = time.time()-t0
print(f'{len(pk)} pairs clipped on one core in {t:.2f} s = {1e6*t/len(pk):.2f} us per pair; whole orc 32.7M pairs would be {32.7e6*t/len(pk):.0f} s on one core, {32.7e6*t/len(pk)/32:.0f} s on 32')
