"""at each no-anchor end of the fine walk on one tube: how many dots pass each part of the anchor test, and a close-up render.
usage: python noanchor.py dumpdir line"""
import sys, os, ast, numpy as np; from scipy.spatial import cKDTree
sys.path.insert(0, r'D:\Stuff\PrintSuit\SegmentationTool'); os.chdir(r'D:\Stuff\PrintSuit\SegmentationTool'); import refine, trace4, trimesh
def main():
    d = sys.argv[1]; i = int(sys.argv[2])
    w = np.load(d+'/lines.npz'); W = [w['pts'][a:a+k] for a, k in zip(np.cumsum(np.r_[0, w['n'][:-1]]), w['n'])]
    line12 = [l for l in open(d+'/log.txt', encoding='utf-8') if l.startswith('--- refine all')][0]
    walk = ast.literal_eval(line12.split('fine walk ')[1].split('; tune')[0]); tune = ast.literal_eval(line12.split('tune ')[-1].strip())
    m = trimesh.load('out/orcfull_patch.stl', process=True); P, F = np.asarray(m.vertices), np.asarray(m.faces); h = float((P.max(0)-P.min(0)).max()); cube = 0.087
    tree, owner = refine.nodes_tree(W); L = W[i]; others = [W[j] for j in refine.end_lines(W, i, tune['slack']*cube, tree, owner)]
    t = refine.tube(L, P, F, cube, h, 3, 0.3, 1.0, True, None, 1, 0.087, 0.087, tune=tune); t = refine.avoid(t, others, 2*0.087, None, res=tune['avoidres'])
    line = t['line']; size = float(t['size']); X = (t['idx'][line]+0.5).astype(float); wgt = np.maximum(t['turn'][line]-1.0, 0)+1e-3; N = t['N'][line]
    lines, ow = trace4.walk(X, wgt, N, size, log=lambda *_: None, **walk); diag = trace4.walk.diag
    step = walk['step']/size; reach = walk['reach']/size; slack = walk['slack']/size; cone = np.cos(np.radians(walk['cone'])); xt = cKDTree(X)
    print(f'line {i}: {len(lines)} pieces; fine cube {size:.4f} mm, step {step:.2f} cubes, reach {reach:.2f} cubes, slack {slack:.2f}')
    spots = []
    for tid, reasons, ends_at, widths in diag:
        Ln = lines[tid]
        for side, r in enumerate(reasons):   # reasons[0] = forward half end, reasons[1] = backward half end
            if r != 'no_anchor': continue
            nodes = Ln[::-1] if side == 0 else Ln                        # the chain ending at this end: last node = the end
            if side == 0: nodes = Ln; q = Ln[-1]; prev = Ln[-7:-1]
            else: q = Ln[0]; prev = Ln[1:7][::-1]
            dd = q-prev[0]; dd /= np.linalg.norm(dd)+1e-12
            c = np.array(xt.query_ball_point(q, reach+1.0/size), int); off = X[c]-q; al = off @ dd; ac = np.linalg.norm(off-al[:, None]*dd, axis=1)
            ahead = (al >= 0.5*step) & (al <= reach); inwid = ac <= 0.5*(0.6/size)+slack; incone = al/np.maximum(np.linalg.norm(off, axis=1), 1e-9) >= cone
            # the nearest dot ahead at all
            fwd = al > 0.5*step; nearest = al[fwd].min()*size if fwd.any() else np.inf
            own = ow[c]; print(f'  piece {tid} end {side}: {len(Ln)} nodes; dots in reach {len(c)}; ahead(0.5 step..reach) {ahead.sum()}, of which within width {(ahead&inwid).sum()}, within cone {(ahead&inwid&incone).sum()}; nearest dot ahead {nearest:.3f} mm; owners of the dots within width ahead: {dict(zip(*np.unique(own[ahead&inwid], return_counts=True)))}')
            spots.append((tid, side, q, dd))
    # close-ups: each no-anchor end, 0.4 mm half-window, looking along the surface normal
    import pyvista as pv
    lo = t['lo']; C = lo+X*size; T = t['turn'][line]
    mesh = pv.PolyData(P, np.hstack([np.full((len(F), 1), 3), F]).ravel())
    def snap(p): Cq, _, f = trimesh.proximity.closest_point(m, p); return Cq+0.005*m.face_normals[f]
    for k, (tid, side, q, dd) in enumerate(spots[:6]):
        c = lo+q*size; half = 0.45
        pl = pv.Plotter(off_screen=True, window_size=(1200, 1000)); pl.set_background('white'); pl.add_mesh(mesh, color='lightgrey', smooth_shading=True)
        sel = (np.abs(C-c) < half).all(1); pd = pv.PolyData(C[sel]+0.003*N[sel]); pd['bend'] = T[sel]; pl.add_mesh(pd, scalars='bend', cmap='plasma', point_size=6, render_points_as_spheres=True, show_scalar_bar=False)
        p = snap(L); pl.add_mesh(pv.PolyData(p, lines=np.r_[len(p), np.arange(len(p))]), color='blue', line_width=3)
        for j, l in enumerate(lines):
            p = snap(lo+l*size); pl.add_mesh(pv.PolyData(p, lines=np.r_[len(p), np.arange(len(p))]), color='red' if j == tid else 'orange', line_width=5)
        pl.add_mesh(pv.Sphere(radius=0.015, center=snap(c[None])[0]), color='black'); pl.add_mesh(pv.Arrow(start=snap(c[None])[0], direction=dd, scale=0.15), color='green')
        n = m.face_normals[trimesh.proximity.closest_point(m, c[None])[2][0]]; pl.camera.focal_point = c; pl.camera.position = c+n*5; pl.camera.up = (0, 0, 1); pl.camera.parallel_projection = True; pl.camera.parallel_scale = half
        pl.add_text(f'line {i} piece {tid}: no-anchor end (black), trend (green), walk blue, pieces red/orange, dots by bend', font_size=9); pl.screenshot(f'out/noanchor_line{i}_{k}.png'); pl.close()
    print('renders', len(spots[:6]))
if __name__ == '__main__': main()
