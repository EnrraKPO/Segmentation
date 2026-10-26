"""Regions of a surface bounded by lines, with the lines as their exact borders.

The mesh is only the canvas: nothing here changes it, and nothing of the triangle layout reaches the result.
  1. embed   : every node of every line is dropped onto its closest surface point; the stretch between two consecutive
               nodes is drawn across the surface as straight strokes, one per triangle it passes (a stroke that leaves a
               triangle is continued in the neighbour from the crossing point on the shared edge; the walk re-aims at the
               target node in each triangle's own plane). Result: per triangle, the strokes in that triangle's flat frame,
               and per crossing an 'edge point' known to both triangles.
  2. partition: inside a crossed triangle the strokes divide the area into pieces (shapely polygonize). Each piece knows
               the stretches of triangle edges it touches (between consecutive edge points). Pieces of neighbouring
               triangles are connected where they share a stretch; uncrossed triangles connect whole; a stroke lying on an
               edge makes that stretch a wall. No connection crosses a stroke, so the connected components are the
               regions and the lines are their exact borders. A line that does not close leaves both sides in one region.
  3. colour  : regions on opposite sides of a stroke are neighbours; greedy colouring with a palette of k colours.
A display mesh (pieces fanned into small triangles, display only) is built for the viewport and for picking.
"""
import os, time
import numpy as np
from numba import njit

TOL = 1e-6    # mm: a point this close to a corner or an edge is on it (far below any geometric scale, above float noise)


# ------------------------------------------------------------------ mesh bookkeeping
def edges_of(F):
    """unique edges E (lo, hi); fedge[f, k] = id of the edge (F[f,k], F[f,k+1]); nbr[f, k] = face across it (-1 at a border)"""
    m = len(F); a = F.ravel(); b = np.roll(F, -1, axis=1).ravel()
    pairs = np.stack([np.minimum(a, b), np.maximum(a, b)], 1)
    E, inv = np.unique(pairs, axis=0, return_inverse=True); inv = inv.ravel()
    order = np.argsort(inv, kind='stable'); es = inv[order]; same = es[1:] == es[:-1]
    i1 = order[:-1][same]; i2 = order[1:][same]
    nbr = -np.ones(3*m, np.int64); nbr[i1] = i2//3; nbr[i2] = i1//3
    return E.astype(np.int64), inv.reshape(m, 3).astype(np.int64), nbr.reshape(m, 3)


def frames(P, F):
    """each triangle's flat frame: origin v0, unit axes e1 (along the first edge) and e2 in the plane; C2 = the three corners in it"""
    v0 = P[F[:, 0]]; e1 = P[F[:, 1]]-v0; L1 = np.linalg.norm(e1, axis=1); e1 = e1/np.maximum(L1, 1e-300)[:, None]
    w = P[F[:, 2]]-v0; n = np.cross(e1, w); n /= np.maximum(np.linalg.norm(n, axis=1), 1e-300)[:, None]; e2 = np.cross(n, e1)
    C2 = np.zeros((len(F), 3, 2)); C2[:, 1, 0] = L1; C2[:, 2, 0] = (w*e1).sum(1); C2[:, 2, 1] = (w*e2).sum(1)
    return v0, e1, e2, C2


# ------------------------------------------------------------------ 1. nodes onto the surface
def project(P, F, pts, log=None):
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


# ------------------------------------------------------------------ 1. strokes across the surface
@njit(cache=True)
def _apply(M, x, y): return M[0]*x+M[1]*y+M[4], M[2]*x+M[3]*y+M[5]


@njit(cache=True)
def _unfold(M, f, k, g, C2, F, Mg):
    """Mg <- the rigid map from g's frame into the base plane that lays g flat next to f across edge k of f (f is already
    in the base plane through M). A face wound the other way is reflected so that it still lands on the far side."""
    k1 = (k+1) % 3; va = F[f, k]; vb = F[f, k1]
    ax, ay = _apply(M, C2[f, k, 0], C2[f, k, 1]); bx, by = _apply(M, C2[f, k1, 0], C2[f, k1, 1])
    ia = 0; ib = 0; ic = 0
    for c in range(3):
        if F[g, c] == va: ia = c
        elif F[g, c] == vb: ib = c
        else: ic = c
    gax = C2[g, ia, 0]; gay = C2[g, ia, 1]; gbx = C2[g, ib, 0]; gby = C2[g, ib, 1]
    ux = gbx-gax; uy = gby-gay; vx = bx-ax; vy = by-ay; L2 = ux*ux+uy*uy
    co = (ux*vx+uy*vy)/L2; si = (ux*vy-uy*vx)/L2
    Mg[0] = co; Mg[1] = -si; Mg[2] = si; Mg[3] = co; Mg[4] = ax-(co*gax-si*gay); Mg[5] = ay-(si*gax+co*gay)
    # the third corners must lie on opposite sides of the edge
    fx, fy = _apply(M, C2[f, (k+2) % 3, 0], C2[f, (k+2) % 3, 1]); gx, gy = _apply(Mg, C2[g, ic, 0], C2[g, ic, 1])
    sf = vx*(fy-ay)-vy*(fx-ax); sg = vx*(gy-ay)-vy*(gx-ax)
    if sf*sg > 0:                                                                     # reflect across the edge line: x -> a + (I-2nn^T)(x-a)
        n2 = vx*vx+vy*vy; nx = -vy/n2**0.5; ny = vx/n2**0.5
        r0 = 1-2*nx*nx; r1 = -2*nx*ny; r3 = 1-2*ny*ny
        m0 = r0*Mg[0]+r1*Mg[2]; m1 = r0*Mg[1]+r1*Mg[3]; m2 = r1*Mg[0]+r3*Mg[2]; m3 = r1*Mg[1]+r3*Mg[3]
        tx = Mg[4]-ax; ty = Mg[5]-ay; m4 = ax+r0*tx+r1*ty; m5 = ay+r1*tx+r3*ty
        Mg[0] = m0; Mg[1] = m1; Mg[2] = m2; Mg[3] = m3; Mg[4] = m4; Mg[5] = m5


@njit(cache=True)
def _walk(P, F, nodes3, nface, node2, start, C2, nbr, fedge, E, depth, cap, prune=0):
    """strokes between consecutive nodes of every line. The stretch from node a (in face fa) to node b (in face fb) is the
    straight segment in the strip of triangles unfolded flat between them: a depth-first search over the neighbours of fa
    unfolds each candidate strip into fa's plane; the strip whose every edge the segment crosses inside the edge, in
    order, is the path. The crossings are the edge points, shared by the two faces of the edge; a crossing within TOL of
    a corner is taken as the corner itself. Returns stroke arrays (face, p, q, edge point ids at p and q or -1, line),
    edge point arrays (edge id, canonical parameter along E[e,0]->E[e,1], 3D position), nodes whose stretch found no
    strip (reason 2) and overflow flag."""
    sf = np.empty(cap, np.int64); sp = np.empty((cap, 2)); sq = np.empty((cap, 2)); s0 = np.empty(cap, np.int64); s1 = np.empty(cap, np.int64); sl = np.empty(cap, np.int64)
    epe = np.empty(cap, np.int64); ept = np.empty(cap); ep3 = np.empty((cap, 3)); fid = np.empty(cap, np.int64); frs = np.empty(cap, np.int64)
    ns = 0; ne = 0; fails = 0; over = False
    D = depth
    cf = np.empty(D+1, np.int64); ck = np.empty(D+1, np.int64); cen = np.empty(D+1, np.int64); cnext = np.empty(D+1, np.int64)
    cM = np.empty((D+1, 6)); ct = np.empty(D+1); Mg = np.empty(6); flo = np.empty(D+1); fhi = np.empty(D+1)
    for li in range(len(start)-1):
        for i in range(start[li], start[li+1]-1):
            if ns+D+2 >= cap or ne+D+2 >= cap: over = True; break
            b = i+1; fa = nface[i]; fb = nface[b]
            au = node2[i, 0]; av = node2[i, 1]; bu = node2[b, 0]; bv = node2[b, 1]
            if fa == fb:
                if (bu-au)**2+(bv-av)**2 > 1e-24:
                    sf[ns] = fa; sp[ns, 0] = au; sp[ns, 1] = av; sq[ns, 0] = bu; sq[ns, 1] = bv; s0[ns] = -1; s1[ns] = -1; sl[ns] = li; ns += 1
                continue
            # depth-first search over strips; a funnel (the directions from a that pass through every edge so far) prunes it
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
                if d == 0: ref = 0.5*(a1+a2) if abs(a1-a2) <= np.pi else 0.5*(a1+a2)+np.pi   # the funnel's reference direction: the first edge's middle
                r1 = a1-ref; r2 = a2-ref
                while r1 > np.pi: r1 -= 2*np.pi
                while r1 < -np.pi: r1 += 2*np.pi
                while r2 > np.pi: r2 -= 2*np.pi
                while r2 < -np.pi: r2 += 2*np.pi
                lo = min(r1, r2); hi = max(r1, r2)
                if hi-lo > np.pi: lo, hi = hi-2*np.pi, lo                                  # the edge spans the back of the circle (a sits almost on it)
                dn1 = ((px-au)**2+(py-av)**2)**0.5; dn2 = ((qx-au)**2+(qy-av)**2)**0.5
                eps = TOL/max(min(dn1, dn2), TOL)                                          # the leaf test allows TOL past an end: the same slack in angle
                nlo = max(flo[d], lo-eps); nhi = min(fhi[d], hi+eps)
                if prune == 1 and nlo > nhi: continue                                     # no straight line from a crosses all these edges
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
            if not found: fid[fails] = i; frs[fails] = 2; fails += 1; continue
            # emit the strokes along the strip
            pu = au; pv = av; prev = -1; f = fa
            for lv in range(d+1):
                k = ck[lv]; k1 = (k+1) % 3; t = ct[lv]; g = nbr[f, k]
                ttol = TOL/((C2[f, k1, 0]-C2[f, k, 0])**2+(C2[f, k1, 1]-C2[f, k, 1])**2)**0.5
                if t < ttol or t > 1-ttol:                                            # through a corner: no edge point, the corner is a ring node of both faces
                    kc = k if t < ttol else k1; v = F[f, kc]; xu = C2[f, kc, 0]; xv = C2[f, kc, 1]
                    if (xu-pu)**2+(xv-pv)**2 > 1e-24:
                        sf[ns] = f; sp[ns, 0] = pu; sp[ns, 1] = pv; sq[ns, 0] = xu; sq[ns, 1] = xv; s0[ns] = prev; s1[ns] = -1; sl[ns] = li; ns += 1
                    for c in range(3):
                        if F[g, c] == v: pu = C2[g, c, 0]; pv = C2[g, c, 1]
                    prev = -1; f = g; continue
                e = fedge[f, k]; va = F[f, k]; vb = F[f, k1]
                tcan = t if va == E[e, 0] else 1.0-t
                xu = C2[f, k, 0]+t*(C2[f, k1, 0]-C2[f, k, 0]); xv = C2[f, k, 1]+t*(C2[f, k1, 1]-C2[f, k, 1])
                epe[ne] = e; ept[ne] = tcan
                for c in range(3): ep3[ne, c] = P[va, c]+t*(P[vb, c]-P[va, c])
                if (xu-pu)**2+(xv-pv)**2 > 1e-24:
                    sf[ns] = f; sp[ns, 0] = pu; sp[ns, 1] = pv; sq[ns, 0] = xu; sq[ns, 1] = xv; s0[ns] = prev; s1[ns] = ne; sl[ns] = li; ns += 1
                kk = 0
                for c in range(3):
                    if fedge[g, c] == e: kk = c
                kk1 = (kk+1) % 3
                tg = tcan if F[g, kk] == E[e, 0] else 1.0-tcan
                pu = C2[g, kk, 0]+tg*(C2[g, kk1, 0]-C2[g, kk, 0]); pv = C2[g, kk, 1]+tg*(C2[g, kk1, 1]-C2[g, kk, 1])
                prev = ne; ne += 1; f = g
            if (bu-pu)**2+(bv-pv)**2 > 1e-24:
                sf[ns] = f; sp[ns, 0] = pu; sp[ns, 1] = pv; sq[ns, 0] = bu; sq[ns, 1] = bv; s0[ns] = prev; s1[ns] = -1; sl[ns] = li; ns += 1
        if over: break
    return sf[:ns], sp[:ns], sq[:ns], s0[:ns], s1[:ns], sl[:ns], epe[:ne], ept[:ne], ep3[:ne], fid[:fails], frs[:fails], over


def embed(P, F, lines, log=None, depth=8, rounds=8):
    """lines (list of (k,3) arrays) onto the surface: see _walk. A stretch that finds no strip (two nodes projected far
    apart on the surface, or across a sharp fold) is bisected: the midpoint of its two surface points is dropped onto the
    surface as a new node, and the walk runs again, up to `rounds` times. Returns dict(strokes, edge points, nodes, ...)."""
    t0 = time.time(); P = np.ascontiguousarray(P, np.float64); F = np.ascontiguousarray(F, np.int64)
    lines = [np.asarray(L, float).reshape(-1, 3) for L in lines if len(L) > 1]
    nodes = np.vstack(lines); start = np.r_[0, np.cumsum([len(L) for L in lines])].astype(np.int64)
    E, fedge, nbr = edges_of(F); v0, e1, e2, C2 = frames(P, F)
    Q, nf, dist = project(P, F, nodes, log=log)
    if log: log(f'  {len(nodes)} nodes onto the surface in {time.time()-t0:.1f} s (distance median {np.median(dist):.3f}, max {dist.max():.3f} mm)')
    def frame2(Q, nf): d = Q-v0[nf]; return np.stack([(d*e1[nf]).sum(1), (d*e2[nf]).sum(1)], 1)
    # consecutive nodes that land on the same surface point (within TOL) are one node; lines left with one node are dropped
    same = np.linalg.norm(np.diff(Q, axis=0), axis=1) < TOL; same = np.r_[False, same]; same[start[1:-1]] = False
    keep = ~same; cum = np.cumsum(keep); start = np.r_[0, cum[start[1:]-1]].astype(np.int64); Q = Q[keep]; nf = nf[keep]; nodes = nodes[keep]
    n_ = np.diff(start); ok = n_ > 1
    if not ok.all():
        drop = np.zeros(len(Q), bool)
        for a, b in zip(start[:-1][~ok], start[1:][~ok]): drop[a:b] = True
        Q = Q[~drop]; nf = nf[~drop]; nodes = nodes[~drop]; start = np.r_[0, np.cumsum(n_[ok])].astype(np.int64)
    node2 = frame2(Q, nf); added = 0; via_vertex = 0; vtree = None
    for rnd in range(rounds+1):
        cap = 4*len(Q)+1000
        while True:
            sf, sp, sq, s0, s1, sl, epe, ept, ep3, fid, frs, over = _walk(P, F, Q, nf, node2, start, C2, nbr, fedge, E, depth, cap)
            if not over: break
            cap *= 2
        if len(fid) == 0 or rnd == rounds: break
        fid = np.unique(fid); fid = fid[np.linalg.norm(Q[fid+1]-Q[fid], axis=1) > 2*TOL]
        if len(fid) == 0: break
        # the new node: the mesh vertex the stretch passes by, when one lies within its span (a straight strip cannot pass
        # a vertex with a fold on the near side: the path goes through the vertex); otherwise the surface point under the middle
        if vtree is None:
            from scipy.spatial import cKDTree
            vtree = cKDTree(P); vface = np.full(len(P), -1, np.int64); vface[F.ravel()[::-1]] = np.repeat(np.arange(len(F)), 3)[::-1]
        mid = 0.5*(Q[fid]+Q[fid+1]); half = 0.5*np.linalg.norm(Q[fid+1]-Q[fid], axis=1)
        dv, iv = vtree.query(mid); mq, mf, md = project(P, F, mid)
        use = (dv <= half) & (np.linalg.norm(P[iv]-Q[fid], axis=1) > TOL) & (np.linalg.norm(P[iv]-Q[fid+1], axis=1) > TOL)
        mq[use] = P[iv[use]]; mf[use] = vface[iv[use]]; via_vertex += int(use.sum())
        near = use | (md <= half)                                                          # a middle that lands outside the stretch's span is a jump between sheets: left as a gap
        fid = fid[near]; mq = mq[near]; mf = mf[near]
        if len(fid) == 0: break
        Q = np.insert(Q, fid+1, mq, axis=0); nf = np.insert(nf, fid+1, mf); node2 = np.insert(node2, fid+1, frame2(mq, mf), axis=0)
        start = start+np.searchsorted(fid, start, side='left'); added += len(fid)                  # a line's offset grows by the insertions before it
    if log: log(f'  {len(sf)} strokes in {len(np.unique(sf))} crossed triangles, {len(epe)} edge crossings, {added} nodes added by bisection ({via_vertex} at a vertex), {len(fid)} stretches not drawn through, {time.time()-t0:.1f} s')
    return dict(P=P, F=F, E=E, fedge=fedge, nbr=nbr, v0=v0, e1=e1, e2=e2, C2=C2, nodes=Q, nface=nf, node2=node2, start=start,
                sf=sf, sp=sp, sq=sq, s0=s0, s1=s1, sl=sl, epe=epe, ept=ept, ep3=ep3, fails=len(fid), fid=fid, added=added)


def loose_ends(emb):
    """line ends that sit on no other line's node (within TOL), and the midpoints of the stretches that were not drawn through"""
    from scipy.spatial import cKDTree
    Q = emb['nodes']; start = emb['start']; line = np.repeat(np.arange(len(start)-1), np.diff(start))
    ends = np.concatenate([start[:-1], start[1:]-1]); t = cKDTree(Q); loose = []
    for i in ends:
        if np.linalg.norm(Q[start[line[i]]]-Q[start[line[i]+1]-1]) < TOL: continue      # a ring closed on itself
        if all(line[j] == line[i] for j in t.query_ball_point(Q[i], TOL)): loose.append(Q[i])
    fid = emb['fid']; gaps = 0.5*(Q[fid]+Q[fid+1]) if len(fid) else np.zeros((0, 3))
    return np.array(loose).reshape(-1, 3), gaps


# ------------------------------------------------------------------ 2. partition of a crossed triangle (worker)
GRID = 1e-7   # snapping grid for the noding, mm in the triangle's frame (far below any geometric scale, above float noise)


def _key(x, y): return (round(x, 7), round(y, 7))


def _partition_face(task):
    """task = (f, C2 (3,2), strokes [(pu,pv,qu,qv,ep0,ep1)], ring [(e, canonical, [(ep, t_local)...]) for k in 0..2]).
    Returns (f, pieces [(area, tris (T,3,2), poly_wkb, stretches [(e,j)...])], neighbours [(i,j)], walls [(e,j)])."""
    import shapely
    from shapely.geometry import LineString, LinearRing, Polygon
    from trimesh.creation import triangulate_polygon
    f, C, strokes, ring = task
    # ring nodes: corners with the edge points in order; coordinates of an edge point from its parameter, the same formula the strokes use
    pts = []; seg_stretch = {}; ep_xy = {}
    for k in range(3):
        e, canon, eps = ring[k]; k1 = (k+1) % 3; c0 = C[k]; c1 = C[k1]; n = len(eps)
        order = sorted(range(n), key=lambda i: eps[i][1])
        chain = [(c0[0], c0[1])]
        for j, i in enumerate(order):
            ep, t = eps[i]; xy = (c0[0]+t*(c1[0]-c0[0]), c0[1]+t*(c1[1]-c0[1])); ep_xy[ep] = xy; chain.append(xy)
        chain.append((c1[0], c1[1]))
        for jl in range(n+1):                                         # local stretch jl -> canonical stretch index
            j = jl if canon else n-jl
            a = _key(*chain[jl]); b = _key(*chain[jl+1]); seg_stretch[(a, b)] = (e, j); seg_stretch[(b, a)] = (e, j)
        pts.extend(chain[:-1])
    ring_geom = LinearRing(pts)
    segs = []; walls = []
    for pu, pv, qu, qv, ep0, ep1 in strokes:
        p = ep_xy[ep0] if ep0 >= 0 else (pu, pv); q = ep_xy[ep1] if ep1 >= 0 else (qu, qv)
        if p == q: continue
        st = seg_stretch.get((_key(*p), _key(*q)))
        if st is not None: walls.append(st); continue                 # a stroke along an edge: the stretch is a wall, no piece to make
        segs.append(LineString([p, q]))
    tri = Polygon(pts)
    if not segs:
        return f, [(tri.area, np.array([C]), shapely.to_wkb(tri), [seg_stretch[(_key(*pts[k]), _key(*pts[(k+1) % len(pts)]))] for k in range(len(pts))])], [], walls
    noded = shapely.unary_union([ring_geom]+segs, grid_size=GRID)
    polys = [g for g in shapely.get_parts(shapely.polygonize([noded])) if g.area > 0]
    pieces = []; seg_owner = {}
    for pi, poly in enumerate(polys):
        stretches = []
        for rg in [poly.exterior, *poly.interiors]:
            xy = np.asarray(rg.coords)
            for a, b in zip(xy[:-1], xy[1:]):
                ka = _key(*a); kb = _key(*b); st = seg_stretch.get((ka, kb))
                if st is not None: stretches.append(st)
                else: seg_owner.setdefault((min(ka, kb), max(ka, kb)), []).append(pi)
        try: v, t = triangulate_polygon(poly, engine='earcut'); tris = v[t]
        except Exception: tris = np.zeros((0, 3, 2))
        pieces.append((poly.area, tris, shapely.to_wkb(poly), stretches))
    nb = sorted({(min(o), max(o)) for o in seg_owner.values() if len(o) == 2 and o[0] != o[1]})
    return f, pieces, nb, walls


def partition(emb, log=None, workers=None):
    """pieces of every crossed triangle, on the cores. Returns dict(faces, pieces per face, neighbours, walls)."""
    from concurrent.futures import ProcessPoolExecutor
    t0 = time.time(); F = emb['F']; E = emb['E']; fedge = emb['fedge']; C2 = emb['C2']
    sf = emb['sf']; order = np.argsort(sf, kind='stable'); faces, first = np.unique(sf[order], return_index=True); bounds = np.r_[first, len(sf)]
    # edge points per (face, local edge): every edge point lies on the rings of both faces of its edge
    epe = emb['epe']; ept = emb['ept']; ep_by_edge = {}
    for i, e in enumerate(epe): ep_by_edge.setdefault(int(e), []).append(i)
    tasks = []
    for fi, f in enumerate(faces):
        idx = order[bounds[fi]:bounds[fi+1]]
        strokes = [(float(emb['sp'][i, 0]), float(emb['sp'][i, 1]), float(emb['sq'][i, 0]), float(emb['sq'][i, 1]), int(emb['s0'][i]), int(emb['s1'][i])) for i in idx]
        ring = []
        for k in range(3):
            e = int(fedge[f, k]); canon = bool(F[f, k] == E[e, 0])
            eps = [(i, float(ept[i]) if canon else 1.0-float(ept[i])) for i in ep_by_edge.get(e, [])]
            ring.append((e, canon, eps))
        tasks.append((int(f), C2[f], strokes, ring))
    workers = workers or os.cpu_count() or 1
    if workers > 1 and len(tasks) > 256:
        with ProcessPoolExecutor(workers) as ex: res = list(ex.map(_partition_face, tasks, chunksize=max(1, len(tasks)//(workers*8))))
    else: res = [_partition_face(t) for t in tasks]
    if log: log(f'  {len(tasks)} crossed triangles partitioned into {sum(len(r[1]) for r in res)} pieces in {time.time()-t0:.1f} s')
    return dict(faces=np.array([r[0] for r in res], np.int64), pieces=[r[1] for r in res], nb=[r[2] for r in res], walls=[r[3] for r in res], ep_by_edge=ep_by_edge)


# ------------------------------------------------------------------ 2. regions, 3. colours
def connect(emb, part, log=None):
    """pieces -> regions. Piece ids: 0..m-1 the (uncrossed) faces themselves, then the pieces of the crossed faces in order.
    Returns dict(region of every piece, number of regions, region neighbours, areas, first piece id of each crossed face)."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    t0 = time.time(); F = emb['F']; E = emb['E']; fedge = emb['fedge']; nbr = emb['nbr']; m = len(F)
    crossed = np.zeros(m, bool); crossed[part['faces']] = True
    pid0 = m+np.r_[0, np.cumsum([len(p) for p in part['pieces']])]            # first piece id per crossed face
    n_pieces = int(pid0[-1]); face_of = np.arange(n_pieces); face_of[m:] = np.repeat(part['faces'], [len(p) for p in part['pieces']])
    owner = {}                                                                     # (edge, stretch j) -> [piece ids on its two sides]
    walls = set()
    for fi, f in enumerate(part['faces']):
        for pi, (area, tris, wkb, stretches) in enumerate(part['pieces'][fi]):
            for st in stretches: owner.setdefault(st, []).append(int(pid0[fi]+pi))
        walls.update(part['walls'][fi])
    ne = {e: len(v) for e, v in part['ep_by_edge'].items()}
    # edges: between two uncrossed faces, connect the faces; otherwise connect the owners of every stretch (an uncrossed side owns all)
    fe = fedge.ravel(); nb = nbr.ravel(); fa = np.repeat(np.arange(m), 3); keep = (nb > fa)
    ea = fe[keep]; A = fa[keep]; B = nb[keep]
    both_plain = ~crossed[A] & ~crossed[B]
    I = list(A[both_plain]); J = list(B[both_plain])
    for e, a, b in zip(ea[~both_plain], A[~both_plain], B[~both_plain]):
        e = int(e)
        for j in range(ne.get(e, 0)+1):
            if (e, j) in walls: continue
            side = list(owner.get((e, j), []))
            for h in (a, b):
                if not crossed[h]: side.append(int(h))
            for x in side[1:]: I.append(side[0]); J.append(x)
    g = coo_matrix((np.ones(len(I)), (np.array(I, np.int64), np.array(J, np.int64))), shape=(n_pieces, n_pieces))
    _, lab = connected_components(g, directed=False)
    valid = np.ones(n_pieces, bool); valid[:m][crossed] = False                     # a crossed face's own id stands for nothing: its pieces do
    keep = np.unique(lab[valid]); remap = np.full(lab.max()+1, -1, np.int64); remap[keep] = np.arange(len(keep)); region = remap[lab]; nreg = len(keep)
    # areas
    T = emb['P'][F]; farea = 0.5*np.linalg.norm(np.cross(T[:, 1]-T[:, 0], T[:, 2]-T[:, 0]), axis=1)
    parea = np.zeros(n_pieces); parea[:m] = farea; parea[:m][crossed] = 0
    for fi in range(len(part['faces'])):
        for pi, pc in enumerate(part['pieces'][fi]): parea[pid0[fi]+pi] = pc[0]
    areas = np.bincount(region[valid], weights=parea[valid], minlength=nreg)
    # neighbours: pieces on the two sides of a stroke, and the owners of a wall stretch
    pairs = set()
    for fi in range(len(part['faces'])):
        for i, j in part['nb'][fi]:
            a, b = region[pid0[fi]+i], region[pid0[fi]+j]
            if a != b: pairs.add((min(a, b), max(a, b)))
    for st in walls:
        side = owner.get(st, [])
        if len(side) == 2 and region[side[0]] != region[side[1]]: a, b = sorted((region[side[0]], region[side[1]])); pairs.add((a, b))
    if log: log(f'  {nreg} regions from {n_pieces} pieces ({int(crossed.sum())} crossed triangles), {len(pairs)} neighbour pairs, {time.time()-t0:.1f} s')
    return dict(region=region, nreg=nreg, pairs=sorted(pairs), areas=areas, pid0=pid0, crossed=crossed, face_of=face_of)


def absorb(con, min_area, log=None):
    """regions below min_area (mm2) are merged into their largest neighbour, smallest first; regions without a neighbour
    stay. Rewrites con (region, nreg, pairs, areas) in place; returns the number absorbed."""
    t0 = time.time(); nreg = con['nreg']; areas = con['areas'].copy(); root = np.arange(nreg)
    def find(r):
        while root[r] != r: root[r] = root[root[r]]; r = root[r]
        return r
    nb = [set() for _ in range(nreg)]
    for a, b in con['pairs']: nb[a].add(b); nb[b].add(a)
    done = 0
    for r in sorted(range(nreg), key=lambda r: areas[r]):
        if areas[r] >= min_area: break
        if find(r) != r: continue
        cand = {find(q) for q in nb[r]}; cand.discard(r)
        if not cand: continue
        t = max(cand, key=lambda q: areas[q]); root[r] = t; areas[t] += areas[r]; nb[t] |= nb[r]; done += 1
    roots = np.array([find(r) for r in range(nreg)]); keep = np.unique(roots); remap = np.full(nreg, -1, np.int64); remap[keep] = np.arange(len(keep))
    reg = con['region']; m = reg >= 0; reg[m] = remap[roots[reg[m]]]
    con['nreg'] = len(keep); con['areas'] = np.bincount(remap[roots], weights=areas*(roots == np.arange(nreg)), minlength=len(keep))
    con['pairs'] = sorted({(min(remap[roots[a]], remap[roots[b]]), max(remap[roots[a]], remap[roots[b]])) for a, b in con['pairs'] if roots[a] != roots[b]})
    if log: log(f'  {done} regions below {min_area} mm2 absorbed into a neighbour: {con["nreg"]} regions, {time.time()-t0:.1f} s')
    return done


PALETTE = ['#E6194B', '#3CB44B', '#FFE119', '#4363D8', '#F58231', '#911EB4', '#42D4F4', '#F032E6', '#BFEF45']   # 9 well-separated colours


def colour(nreg, pairs, k=9, log=None):
    """greedy colouring, most-neighboured regions first; when no colour is free the one least used among the neighbours is
    taken and the clash counted. Returns (colour per region, clashes)."""
    nb = [[] for _ in range(nreg)]
    for a, b in pairs: nb[a].append(b); nb[b].append(a)
    col = -np.ones(nreg, np.int64); clashes = []
    for r in sorted(range(nreg), key=lambda r: -len(nb[r])):
        used = np.zeros(k, np.int64)
        for q in nb[r]:
            if col[q] >= 0: used[col[q]] += 1
        free = np.where(used == 0)[0]
        if len(free): col[r] = free[0]
        else: col[r] = int(np.argmin(used)); clashes.append(r)
    if log: log(f'  coloured with {k} colours, {len(clashes)} regions share a colour with a neighbour')
    return col, clashes


# ------------------------------------------------------------------ display mesh
def display_mesh(emb, part, con):
    """small triangles for the viewport: uncrossed faces as they are, crossed faces as their pieces fanned. Returns
    (vertices, faces, piece id per display triangle)."""
    P = emb['P']; F = emb['F']; m = len(F); crossed = con['crossed']; pid0 = con['pid0']
    V = [P]; T = [F[~crossed]]; pid = [np.where(~crossed)[0]]; nv = len(P)
    v0 = emb['v0']; e1 = emb['e1']; e2 = emb['e2']
    for fi, f in enumerate(part['faces']):
        for pi, (area, tris, wkb, stretches) in enumerate(part['pieces'][fi]):
            if len(tris) == 0: continue
            xy = tris.reshape(-1, 2); X = v0[f]+xy[:, :1]*e1[f]+xy[:, 1:]*e2[f]
            V.append(X); T.append(nv+np.arange(len(X)).reshape(-1, 3)); pid.append(np.full(len(tris), pid0[fi]+pi)); nv += len(X)
    return np.vstack(V), np.vstack(T), np.concatenate(pid)


# ------------------------------------------------------------------ all of it
def run(P, F, lines, log=None, k=9, workers=None, min_area=0.0, depth=8):
    t0 = time.time()
    if log: log(f'embedding {len(lines)} lines on {len(F)} triangles')
    emb = embed(P, F, lines, log=log, depth=depth)
    part = partition(emb, log=log, workers=workers)
    con = connect(emb, part, log=log)
    if min_area > 0: absorb(con, min_area, log=log)
    col, clashes = colour(con['nreg'], con['pairs'], k=k, log=log)
    DV, DF, pid = display_mesh(emb, part, con)
    if log: log(f'regions done in {time.time()-t0:.1f} s')
    return dict(emb=emb, part=part, con=con, colour=col, clashes=clashes, DV=DV, DF=DF, pid=pid, k=k)
