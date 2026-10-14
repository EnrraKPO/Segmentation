"""Shape features on a single grid fixed to the model: a shell of voxels over the surface.

Each surface cell holds the surface direction at the surface point nearest its centre (interpolated, as a texel does). A line is where the
direction turns fastest across the surface, measured over the whole neighbourhood at once (shape operator fit),
judged per cell against `sensitivity`
degrees per cell. Nothing else: no averaging, no viewpoint.

run(P, F, cells=300, sensitivity=10) -> dict(idx (k,3) cell indices, N (k,3) directions, turn (k,), line (k,) bool,
                                            face (k,) a mesh face in the cell, size (mm))
"""
import numpy as np
from scipy.spatial import cKDTree

OFFS = np.array([[i, j, k] for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1) if (i, j, k) != (0, 0, 0)])


def _cell_samples(P, F, VN, size, spp):
    """sample points on every face densely enough that every cell the surface passes through gets hit"""
    A, B, C = P[F[:, 0]], P[F[:, 1]], P[F[:, 2]]
    longest = np.maximum.reduce([np.linalg.norm(B-A, axis=1), np.linalg.norm(C-B, axis=1), np.linalg.norm(A-C, axis=1)])
    n_sub = np.maximum(1, np.ceil(spp*longest/size)).astype(int)         # subdivisions per face edge
    pts, nrm, fid = [], [], []
    for k in np.unique(n_sub):
        sel = np.nonzero(n_sub == k)[0]
        u = np.arange(k)+0.5
        bu, bv = np.meshgrid(u, u); m = (bu+bv) < k; bu, bv = bu[m]/k, bv[m]/k                # points in the triangle
        S = A[sel][:, None]+(B[sel]-A[sel])[:, None]*bu[None, :, None]+(C[sel]-A[sel])[:, None]*bv[None, :, None]
        Ns = VN[F[sel, 0]][:, None]*(1-bu-bv)[None, :, None]+VN[F[sel, 1]][:, None]*bu[None, :, None]+VN[F[sel, 2]][:, None]*bv[None, :, None]
        pts.append(S.reshape(-1, 3)); nrm.append(Ns.reshape(-1, 3)); fid.append(np.repeat(sel, len(bu)))
    return np.vstack(pts), np.vstack(nrm), np.concatenate(fid)


def run(P, F, cells=300, sensitivity=10.0, log=None):
    lo = P.min(0); ext = P.max(0)-lo; size = ext.max()/cells
    import trimesh
    VN = trimesh.Trimesh(P, F, process=False).vertex_normals
    pts, nrm, fid = _cell_samples(P, F, VN, size, 2.0)
    ijk = np.floor((pts-lo)/size).astype(np.int64)
    key = (ijk[:, 0]*4096+ijk[:, 1])*4096+ijk[:, 2]
    uk, inv = np.unique(key, return_inverse=True); k = len(uk)
    idx = np.stack([uk // (4096*4096), (uk // 4096) % 4096, uk % 4096], 1)
    # direction per cell = the interpolated surface direction at the surface point nearest the cell centre (a texel)
    centre = lo+(idx[inv]+0.5)*size; dc = np.linalg.norm(pts-centre, axis=1)
    order = np.lexsort((dc, inv)); first = np.r_[True, inv[order][1:] != inv[order][:-1]]; pick = order[first]
    N = nrm[pick]; N /= np.linalg.norm(N, axis=1, keepdims=True)+1e-12
    face = fid[pick]
    if log: log(f'grid {cells} cells/height, cell {size:.3f} mm: {k} surface cells from {len(pts)} samples')
    # neighbours: the 26 surrounding cells that are also surface cells
    lut = dict(zip(uk.tolist(), range(k)))
    nb = np.full((k, 26), -1, np.int64)
    for q, o in enumerate(OFFS):
        kk = ((idx[:, 0]+o[0])*4096+(idx[:, 1]+o[1]))*4096+(idx[:, 2]+o[2])
        nb[:, q] = [lut.get(int(x), -1) for x in kk]
    has = nb >= 0
    # turn of the direction across the whole neighbourhood at once: least-squares fit of dN = S . dx over all present
    # neighbours (dx in the cell's tangent plane, dN the change of direction). S is the shape operator; its largest
    # eigenvalue is the strongest rate of turn, its eigenvector the direction across which the surface turns.
    ref = cells/1000.0                                                  # degrees per reference cell (1/1000 height)
    a = np.where((np.abs(N[:, 0]) < 0.9)[:, None], np.array([1., 0, 0]), np.array([0, 1., 0]))
    t1 = np.cross(N, a); t1 /= np.linalg.norm(t1, axis=1, keepdims=True)+1e-12; t2 = np.cross(N, t1)
    A11 = np.zeros(k); A12 = np.zeros(k); A22 = np.zeros(k); B11 = np.zeros(k); B12 = np.zeros(k); B21 = np.zeros(k); B22 = np.zeros(k)
    for q, o in enumerate(OFFS):
        h = has[:, q]; j = nb[h, q]
        dx = (o.astype(float))[None, :]*np.ones((h.sum(), 1))               # offset in cells
        dx -= (dx*N[h]).sum(1, keepdims=True)*N[h]                         # in the tangent plane
        dn = N[j]-N[h]
        ex, ey = (dx*t1[h]).sum(1), (dx*t2[h]).sum(1); nx, ny = (dn*t1[h]).sum(1), (dn*t2[h]).sum(1)
        A11[h] += ex*ex; A12[h] += ex*ey; A22[h] += ey*ey
        B11[h] += nx*ex; B12[h] += nx*ey; B21[h] += ny*ex; B22[h] += ny*ey
    det = A11*A22-A12*A12+1e-12; I11, I12, I22 = A22/det, -A12/det, A11/det
    S11 = B11*I11+B12*I12; S12 = B11*I12+B12*I22; S21 = B21*I11+B22*I12; S22 = B21*I12+B22*I22
    s12 = 0.5*(S12+S21); tr = S11+S22; dd = np.sqrt(np.maximum(((S11-S22)/2)**2+s12**2, 0))
    l1, l2 = tr/2+dd, tr/2-dd
    big = np.abs(l1) >= np.abs(l2); k1 = np.where(big, l1, l2)
    turn = np.degrees(np.abs(k1))*ref                                      # rate of turn, degrees per reference cell
    vx = np.where(np.abs(s12) > 1e-12, s12, np.where(np.abs(S11-k1) < np.abs(S22-k1), 0.0, 1.0))
    vy = np.where(np.abs(s12) > 1e-12, k1-S11, np.where(np.abs(S11-k1) < np.abs(S22-k1), 1.0, 0.0))
    nv = np.sqrt(vx*vx+vy*vy)+1e-18; d1 = (vx/nv)[:, None]*t1+(vy/nv)[:, None]*t2
    nnb = has.sum(1); turn[nnb < 6] = 0                                    # a cell needs a neighbourhood to be judged
    # crest: the turn must peak across d1 (nearest cells along +-d1)
    def along(sign):
        step = np.rint(sign*d1).astype(int)
        kk = ((idx[:, 0]+step[:, 0])*4096+(idx[:, 1]+step[:, 1]))*4096+(idx[:, 2]+step[:, 2])
        j = np.array([lut.get(int(x), -1) for x in kk]); return np.where(j >= 0, turn[np.maximum(j, 0)], 0.0)
    tp, tm = along(1.0), along(-1.0)
    line = (turn >= sensitivity) & (turn >= tp) & (turn > tm)
    if log: log(f'turn per cell: p50 {np.percentile(turn, 50):.1f} p90 {np.percentile(turn, 90):.1f} p99 {np.percentile(turn, 99):.1f} deg; line cells {line.sum()}')
    return dict(idx=idx, N=N, turn=turn, line=line, face=face, size=size, lo=lo, d1=d1)
