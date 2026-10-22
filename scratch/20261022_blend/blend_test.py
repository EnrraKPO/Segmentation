"""Blended arrow test on line 318's tube: each fine cube's arrow = the facings of ALL triangles inside the cube, weighted by
the area each covers in the cube (area measured by dense sample points), against the current rule (the one triangle
nearest the cube centre). Same tube, same bend read on the same rungs, same walk, both ways. usage: python blend_test.py [cut]"""
import sys, os, time, numpy as np, trimesh
from scipy.spatial import cKDTree
sys.path.insert(0, r'D:\Stuff\PrintSuit\SegmentationTool'); os.chdir(r'D:\Stuff\PrintSuit\SegmentationTool'); import voxel, refine, trace4
cut = float(sys.argv[1]) if len(sys.argv) > 1 else 17.0
m = trimesh.load('out/orcfull_patch.stl', process=False); P = np.asarray(m.vertices, float); F = np.asarray(m.faces); h = float((P.max(0)-P.min(0)).max())
d = np.load('out/debug/20261020_144432/lines.npz'); n = d['n']; pts = d['pts']; off = np.r_[0, np.cumsum(n)]; L = pts[off[318]:off[319]]
cube, finer, radius, area, unit = 0.087, 3, 0.3, 0.087, 0.087; tn = refine.TUNE; size = cube/finer; cells = h/size
box = (L.min(0)-radius-size, L.max(0)+radius+size)
r = voxel.run_direct(P, F, cells, None, height=h, flat=True, box=box, workers=1)
centre = r['lo']+(r['idx']+0.5)*r['size']; dist = cKDTree(refine.resample(L, tn['tuberes']*size)).query(centre)[0]; keep = dist <= radius
idx, N1, C, size = r['idx'][keep], r['N'][keep], centre[keep], float(r['size']); lo = r['lo']
# blended arrow: sample points on every triangle whose box meets the tube's box, binned to cubes; a cube's arrow = mean facing of its samples
KB = voxel.KB; key = (idx[:, 0]*KB+idx[:, 1])*KB+idx[:, 2]
T = P[F]; inbox = (T.max(1) >= box[0]-size).all(1) & (T.min(1) <= box[1]+size).all(1); Fb = F[inbox]
spp = 8; S, Ns, fid = voxel._cell_samples(P, Fb, None, size, spp, flat=True)
ijk = np.floor((S-lo)/size).astype(np.int64); skey = (ijk[:, 0]*KB+ijk[:, 1])*KB+ijk[:, 2]
order = np.argsort(key); pos = np.searchsorted(key[order], skey); pos = np.minimum(pos, len(key)-1); hit = key[order][pos] == skey
cube_of = order[pos[hit]]; Nsum = np.zeros((len(idx), 3)); np.add.at(Nsum, cube_of, Ns[hit]); cnt = np.bincount(cube_of, minlength=len(idx))
ntri = np.zeros(len(idx), int); pair = np.unique(cube_of*10**7+fid[hit]); np.add.at(ntri, pair // 10**7, 1)
N2 = np.where(cnt[:, None] > 0, Nsum, N1); N2 /= np.linalg.norm(N2, axis=1, keepdims=True)+1e-12
ang = np.degrees(np.arccos(np.clip((N1*N2).sum(1), -1, 1)))
print(f'tube: {len(idx)} fine cubes of {size:.4f} mm, {inbox.sum()} triangles sampled, {len(S)} samples, samples per cube median {np.median(cnt[cnt>0]):.0f}')
print(f'cubes with 1 triangle {int((ntri==1).sum())}, 2 {int((ntri==2).sum())}, 3+ {int((ntri>=3).sum())}, none hit {int((cnt==0).sum())}')
print(f'arrow change nearest -> blended (deg): p50 {np.median(ang):.1f} p90 {np.percentile(ang, 90):.1f} max {ang.max():.1f}; cubes changed > 1 deg: {int((ang>1).sum())}')
# the same rungs, the same bend read, the same walk
Lr = refine.resample(L, tn['rungstep']*cube); Tg = np.gradient(Lr, axis=0); Tg /= np.linalg.norm(Tg, axis=1, keepdims=True)+1e-12
dr, jr = cKDTree(Lr).query(C); along = ((C-Lr[jr])*Tg[jr]).sum(1); at = np.flatnonzero(np.abs(along) <= tn['rungw']*size)
rung = jr[at]
out = {}
for name, N in (('nearest', N1), ('blended', N2)):
    turn, d1 = refine.bend(C, N, area, unit, None, at=at, nbh=tn['nbh'])
    band = turn >= cut
    # peak bend per rung along the line
    pk = np.full(len(Lr), 0.0); np.maximum.at(pk, rung, turn); pk = pk[np.bincount(rung, minlength=len(Lr)) > 0]
    jump = np.abs(np.diff(pk))
    t = dict(idx=idx[at], N=N[at], lo=lo, size=size, turn=turn, line=band, crest=band, d1=d1, cut=cut)
    lines = trace4.run(t, log=lambda *_: None, workers=1, **refine.WALK)
    print(f'{name}: dots {int(band.sum())} of {len(at)} rung cubes ({100*band.mean():.0f} %); peak bend per rung min {pk.min():.1f} p10 {np.percentile(pk,10):.1f} p50 {np.median(pk):.1f}; '
          f'rung-to-rung peak jump p50 {np.median(jump):.1f} p90 {np.percentile(jump,90):.1f} max {jump.max():.1f} deg; walk: {len(lines)} line(s), nodes {[len(l) for l in lines]}')
    out[name] = dict(turn=turn, lines=lines, pk=pk)
np.savez('out/blend_test.npz', C=C[at], N1=N1[at], N2=N2[at], turn1=out['nearest']['turn'], turn2=out['blended']['turn'], ang=ang[at], ntri=ntri[at],
         n1=[len(l) for l in out['nearest']['lines']], p1=np.vstack(out['nearest']['lines']) if out['nearest']['lines'] else np.zeros((0, 3)),
         n2=[len(l) for l in out['blended']['lines']], p2=np.vstack(out['blended']['lines']) if out['blended']['lines'] else np.zeros((0, 3)), coarse=L)
if out['nearest']['lines'] and out['blended']['lines']:
    A = np.vstack(out['nearest']['lines']); B = np.vstack(out['blended']['lines']); dd = cKDTree(B).query(A)[0]
    print(f'line difference nearest vs blended: median {np.median(dd):.3f} mm, max {dd.max():.3f} mm')
print('saved out/blend_test.npz')
