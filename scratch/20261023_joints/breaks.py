"""at every break between two pieces of one walk line, and at every shy end: are there fine tube dots on the walk line in the
gap (the walk stopped) or none (the clustering split the dots)? Prints the numbers; renders a few spots to out/joints_*.png"""
import sys, os, numpy as np; from scipy.spatial import cKDTree
d = sys.argv[1]; cube = 0.087; os.makedirs('out', exist_ok=True)
w = np.load(d+'/lines.npz'); W = [w['pts'][a:a+k] for a, k in zip(np.cumsum(np.r_[0, w['n'][:-1]]), w['n'])]
r = np.load(d+'/refined.npz'); R = [r['pts'][a:a+k] for a, k in zip(np.cumsum(np.r_[0, r['n'][:-1]]), r['n'])]; ri = r['index']
D = r['dots']; T = r['turn']; dt = cKDTree(D); print('tube dots', len(D), 'fine cube', float(r['size']))
pieces = np.bincount(ri)
def along(L, a, b, step=0.03):
    """points on walk line L between the nearest nodes to a and b"""
    t = cKDTree(L); i, j = t.query(a)[1], t.query(b)[1]; i, j = min(i, j), max(i, j); S = L[i:j+1]
    if len(S) < 2: return S
    s = np.r_[0, np.cumsum(np.linalg.norm(np.diff(S, axis=0), axis=1))]; u = np.arange(0, s[-1], step); return np.stack([np.interp(u, s, S[:, k]) for k in range(3)], 1)
# breaks
rows = []
for i in np.flatnonzero(pieces >= 2):
    P = [R[k] for k in np.flatnonzero(ri == i)]; L = W[i]; t = cKDTree(L)
    pos = sorted([(t.query(p[0])[1], t.query(p[-1])[1], k) for k, p in enumerate(P)])
    for (a0, a1, k), (b0, b1, m) in zip(pos[:-1], pos[1:]):
        e = P[k][-1] if a1 >= a0 else P[k][0]; s = P[m][0] if b1 >= b0 else P[m][-1]; S = along(L, e, s)
        if len(S) < 2: continue
        dd = dt.query(S)[0]; gap = np.linalg.norm(s-e); rows.append((i, gap, dd.max(), np.median(T[dt.query(S)[1]])))
rows = np.array(rows); print(f'breaks {len(rows)}: gap p50 {np.percentile(rows[:,1],50):.2f} mm; largest hole in the tube dots along the walk line across the break: p50 {np.percentile(rows[:,2],50):.3f} p90 {np.percentile(rows[:,2],90):.3f} mm; breaks with a hole over 3 units (0.26 mm, cluster split) {(rows[:,2]>0.26).mean()*100:.0f} %, over one step (0.087) {(rows[:,2]>0.087).mean()*100:.0f} %')
print(f'  median bend of the nearest tube dots along the break: p50 {np.percentile(rows[:,3],50):.1f} deg (whole tube p50 {np.percentile(T,50):.1f})')
# shy ends
wt = cKDTree(np.vstack(W)); wo = np.repeat(np.arange(len(W)), w['n']); ends = []
for i in np.flatnonzero(pieces > 0):
    L = W[i]
    for e in (L[0], L[-1]):
        dd, j = wt.query(e, k=50); m = wo[j] != i
        if not m.any() or dd[m][0] > cube: continue
        o = wo[j][m][0]; E = np.array([q for k in np.flatnonzero(ri == i) for q in (R[k][0], R[k][-1])]); re = E[np.argmin(np.linalg.norm(E-e, axis=1))]
        S = along(L, re, e); hole = dt.query(S)[0].max() if len(S) > 1 else 0.0
        ends.append((i, o, np.linalg.norm(re-e), hole))
ends = np.array(ends); print(f'shy ends {len(ends)}: refined end short of the walk end by p50 {np.percentile(ends[:,2],50):.2f} mm; largest hole in tube dots along that stretch p50 {np.percentile(ends[:,3],50):.3f} p90 {np.percentile(ends[:,3],90):.3f}; stretches with a hole over 0.087 {(ends[:,3]>0.087).mean()*100:.0f} %, over 0.174 (the end rule) {(ends[:,3]>0.174).mean()*100:.0f} %')
np.savez('out/joints_rows.npz', breaks=rows, ends=ends)
# renders: 2 breaks with no hole, 1 break with a hole, 2 shy ends
import pyvista as pv, trimesh
m = trimesh.load('out/orcfull_patch.stl', process=False); mesh = pv.PolyData(np.asarray(m.vertices), np.hstack([np.full((len(m.faces), 1), 3), m.faces]).ravel())
def snap(p): C, _, f = trimesh.proximity.closest_point(m, p); return C+0.01*m.face_normals[f]
def render(name, centre, i, half=1.2):
    pl = pv.Plotter(off_screen=True, window_size=(1400, 1000)); pl.set_background('white'); pl.add_mesh(mesh, color='lightgrey', smooth_shading=True)
    sel = (np.abs(D-centre) < half).all(1); pd = pv.PolyData(D[sel]); pd['bend'] = T[sel]; pl.add_mesh(pd, scalars='bend', cmap='plasma', point_size=3, render_points_as_spheres=True, show_scalar_bar=False)
    for j, L in enumerate(W):
        if (np.abs(L-centre) < half).any(1).any(): p = snap(L); pl.add_mesh(pv.PolyData(p, lines=np.r_[len(p), np.arange(len(p))]), color='blue' if j == i else 'cyan', line_width=3 if j == i else 2)
    for k in np.flatnonzero((np.abs(np.vstack(R)-centre) < half).all(1)[np.cumsum(np.r_[0, r['n'][:-1]])] | True):
        if (np.abs(R[k]-centre) < half).any(1).any(): p = snap(R[k]); pl.add_mesh(pv.PolyData(p, lines=np.r_[len(p), np.arange(len(p))]), color='red' if ri[k] == i else 'orange', line_width=4 if ri[k] == i else 2)
    n = m.face_normals[trimesh.proximity.closest_point(m, centre[None])[2][0]]; pl.camera.focal_point = centre; pl.camera.position = centre+n*8; pl.camera.up = (0, 0, 1); pl.camera.parallel_projection = True; pl.camera.parallel_scale = half
    pl.add_text(name, font_size=10); pl.screenshot(f'out/joints_{name}.png'); pl.close()
nohole = rows[rows[:, 2] <= 0.087]; hole = rows[rows[:, 2] > 0.26]; rng = np.random.default_rng(1)
for tag, src in (('break_nohole', nohole), ('break_hole', hole)):
    for q, row in enumerate(src[rng.choice(len(src), 2, replace=False)]):
        i = int(row[0]); P = [R[k] for k in np.flatnonzero(ri == i)]; E = np.vstack([np.array([p[0], p[-1]]) for p in P]); t = cKDTree(E); pr = t.query_pairs(row[1]*1.01); 
        c = np.mean([E[a] for a in list(pr)[0]], 0) if pr else E[0]; render(f'{tag}{q}_line{i}', c, i)
for q, row in enumerate(ends[rng.choice(len(ends), 2, replace=False)]):
    i = int(row[0]); L = W[i]; render(f'shy{q}_line{i}', L[0] if row[2] < 1e9 else L[-1], i)
print('renders written')
