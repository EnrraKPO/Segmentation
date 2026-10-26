"""one-off: regions.py: exact projection through VTK's cell locator, consecutive nodes within TOL merged, funnel pruning in
the strip search (depth default 32)"""
import os
root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); path = os.path.join(root, 'regions.py')
src = open(path, encoding='utf-8').read()

# 1. project via VTK
a = src.index('def _closest_on_triangles('); b = src.index('# ------------------------------------------------------------------ 1. strokes across the surface')
src = src[:a]+'''def project(P, F, pts, log=None):
    """closest surface point and its triangle for every point (VTK's static cell locator: exact, about a microsecond a
    point after a build of under a second for a million triangles)"""
    import vtk
    from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray
    poly = vtk.vtkPolyData(); vp = vtk.vtkPoints(); vp.SetData(numpy_to_vtk(np.ascontiguousarray(P, np.float64), deep=True)); poly.SetPoints(vp)
    cells = vtk.vtkCellArray(); cells.SetData(3, numpy_to_vtkIdTypeArray(np.ascontiguousarray(F, np.int64).ravel(), deep=True)) if hasattr(cells, 'SetData') else None
    if not hasattr(cells, 'SetData'):
        cells.SetCells(len(F), numpy_to_vtkIdTypeArray(np.hstack([np.full((len(F), 1), 3), F]).ravel().astype(np.int64), deep=True))
    poly.SetPolys(cells)
    loc = vtk.vtkStaticCellLocator(); loc.SetDataSet(poly); loc.BuildLocator()
    n = len(pts); Q = np.zeros((n, 3)); fid = np.zeros(n, np.int64); dist = np.zeros(n)
    cid = vtk.reference(0); sub = vtk.reference(0); d2 = vtk.reference(0.0); cp = [0.0, 0.0, 0.0]
    for i in range(n):
        loc.FindClosestPoint(pts[i], cp, cid, sub, d2); Q[i] = cp; fid[i] = cid.get(); dist[i] = d2.get()**0.5
    return Q, fid, dist


'''+src[b:]

# 2. dedupe + 3. funnel pruning: new _walk DFS block
old = src[src.index('            # depth-first search over strips'):src.index('            if not found: fid[fails] = i; frs[fails] = 2; fails += 1; continue')]
new = '''            # depth-first search over strips; a funnel (the directions from a that pass through every edge so far) prunes it
            d = 0; cf[0] = fa; cen[0] = -1; cnext[0] = 0; flo[0] = -4.0; fhi[0] = 4.0
            cM[0, 0] = 1.0; cM[0, 1] = 0.0; cM[0, 2] = 0.0; cM[0, 3] = 1.0; cM[0, 4] = 0.0; cM[0, 5] = 0.0
            found = False; ref = 0.0
            while d >= 0:
                if cnext[d] >= 3: d -= 1; continue
                k = cnext[d]; cnext[d] += 1
                if k == cen[d]: continue
                f = cf[d]; g = nbr[f, k]
                if g < 0: continue
                k1 = (k+1) % 3
                px, py = _apply(cM[d], C2[f, k, 0], C2[f, k, 1]); qx, qy = _apply(cM[d], C2[f, k1, 0], C2[f, k1, 1])
                a1 = np.arctan2(py-av, px-au); a2 = np.arctan2(qy-av, qx-au)
                if d == 0 and k == cnext[0]-1 and flo[0] == -4.0 and fhi[0] == 4.0: pass
                if d == 0: ref = 0.5*(a1+a2) if abs(a1-a2) <= np.pi else 0.5*(a1+a2)+np.pi   # the funnel's reference direction: the first edge's middle
                r1 = a1-ref; r2 = a2-ref
                while r1 > np.pi: r1 -= 2*np.pi
                while r1 < -np.pi: r1 += 2*np.pi
                while r2 > np.pi: r2 -= 2*np.pi
                while r2 < -np.pi: r2 += 2*np.pi
                lo = min(r1, r2); hi = max(r1, r2)
                if hi-lo > np.pi: lo, hi = hi-2*np.pi, lo                                  # the edge spans the back of the circle (a sits almost on it)
                nlo = max(flo[d], lo-1e-9); nhi = min(fhi[d], hi+1e-9)
                if nlo > nhi: continue                                                    # no straight line from a crosses all these edges
                _unfold(cM[d], f, k, g, C2, F, Mg)
                ck[d] = k
                if g == fb:
                    bx, by = _apply(Mg, bu, bv); dx = bx-au; dy = by-av
                    ok = True; sprev = -1.0
                    for lv in range(d+1):
                        ff = cf[lv]; kk = ck[lv]; kk1 = (kk+1) % 3
                        px, py = _apply(cM[lv], C2[ff, kk, 0], C2[ff, kk, 1]); qx, qy = _apply(cM[lv], C2[ff, kk1, 0], C2[ff, kk1, 1])
                        ex = qx-px; ey = qy-py; den = dx*ey-dy*ex
                        if abs(den) < 1e-300: ok = False; break
                        rx = px-au; ry = py-av
                        s = (rx*ey-ry*ex)/den; t = (rx*dy-ry*dx)/den; ttol = TOL/(ex*ex+ey*ey)**0.5
                        if s < sprev-1e-12 or s < -1e-12 or s > 1+1e-12 or t < -ttol or t > 1+ttol: ok = False; break
                        sprev = s; ct[lv] = min(max(t, 0.0), 1.0)
                    if ok: found = True; break
                if d+1 < D:
                    d += 1; cf[d] = g; cnext[d] = 0; flo[d] = nlo; fhi[d] = nhi
                    for c in range(6): cM[d, c] = Mg[c]
                    for c in range(3):
                        if fedge[g, c] == fedge[f, k]: cen[d] = c
'''
src = src.replace(old, new)
src = src.replace("    cM = np.empty((D+1, 6)); ct = np.empty(D+1); Mg = np.empty(6)\n", "    cM = np.empty((D+1, 6)); ct = np.empty(D+1); Mg = np.empty(6); flo = np.empty(D+1); fhi = np.empty(D+1)\n")
src = src.replace("                if d == 0 and k == cnext[0]-1 and flo[0] == -4.0 and fhi[0] == 4.0: pass\n", "")

# dedupe in embed, depth default 32
old = "    def frame2(Q, nf): d = Q-v0[nf]; return np.stack([(d*e1[nf]).sum(1), (d*e2[nf]).sum(1)], 1)\n    node2 = frame2(Q, nf); added = 0\n"
new = '''    def frame2(Q, nf): d = Q-v0[nf]; return np.stack([(d*e1[nf]).sum(1), (d*e2[nf]).sum(1)], 1)
    # consecutive nodes that land on the same surface point (within TOL) are one node; lines left with one node are dropped
    same = np.linalg.norm(np.diff(Q, axis=0), axis=1) < TOL; same = np.r_[False, same]; same[start[1:-1]] = False
    keep = ~same; cum = np.cumsum(keep); start = np.r_[0, cum[start[1:]-1]].astype(np.int64); Q = Q[keep]; nf = nf[keep]; nodes = nodes[keep]
    n_ = np.diff(start); ok = n_ > 1
    if not ok.all():
        drop = np.zeros(len(Q), bool)
        for a, b in zip(start[:-1][~ok], start[1:][~ok]): drop[a:b] = True
        Q = Q[~drop]; nf = nf[~drop]; nodes = nodes[~drop]; start = np.r_[0, np.cumsum(n_[ok])].astype(np.int64)
    node2 = frame2(Q, nf); added = 0
'''
assert old in src; src = src.replace(old, new)
src = src.replace("def embed(P, F, lines, log=None, depth=8, rounds=8):", "def embed(P, F, lines, log=None, depth=32, rounds=8):")
src = src.replace("def run(P, F, lines, log=None, k=9, workers=None, min_area=0.0, depth=8):", "def run(P, F, lines, log=None, k=9, workers=None, min_area=0.0, depth=32):")
src = src.replace("        fid = np.unique(fid); mid = 0.5*(Q[fid]+Q[fid+1]); mq, mf, _ = project(P, F, mid)\n",
                  "        fid = np.unique(fid); fid = fid[np.linalg.norm(Q[fid+1]-Q[fid], axis=1) > 2*TOL]\n        if len(fid) == 0: break\n        mid = 0.5*(Q[fid]+Q[fid+1]); mq, mf, _ = project(P, F, mid)\n")
open(path, 'w', encoding='utf-8').write(src); print('patched')
