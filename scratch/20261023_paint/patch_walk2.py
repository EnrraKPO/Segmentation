"""one-off: replace the aim-based walk in regions.py by the unfolding walk (the stretch between two nodes is a straight
line in the unfolded strip of triangles between them)"""
import os
root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); path = os.path.join(root, 'regions.py')
src = open(path, encoding='utf-8').read()
a = src.index('@njit(cache=True)\ndef _face_at_vertex('); b = src.index('def embed(')
new = r'''@njit(cache=True)
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
def _walk(P, F, nodes3, nface, node2, start, C2, nbr, fedge, E, depth, cap):
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
    cM = np.empty((D+1, 6)); ct = np.empty(D+1); Mg = np.empty(6)
    for li in range(len(start)-1):
        for i in range(start[li], start[li+1]-1):
            if ns+D+2 >= cap or ne+D+2 >= cap: over = True; break
            b = i+1; fa = nface[i]; fb = nface[b]
            au = node2[i, 0]; av = node2[i, 1]; bu = node2[b, 0]; bv = node2[b, 1]
            if fa == fb:
                if (bu-au)**2+(bv-av)**2 > 1e-24:
                    sf[ns] = fa; sp[ns, 0] = au; sp[ns, 1] = av; sq[ns, 0] = bu; sq[ns, 1] = bv; s0[ns] = -1; s1[ns] = -1; sl[ns] = li; ns += 1
                continue
            # depth-first search over strips
            d = 0; cf[0] = fa; cen[0] = -1; cnext[0] = 0
            cM[0, 0] = 1.0; cM[0, 1] = 0.0; cM[0, 2] = 0.0; cM[0, 3] = 1.0; cM[0, 4] = 0.0; cM[0, 5] = 0.0
            found = False
            while d >= 0:
                if cnext[d] >= 3: d -= 1; continue
                k = cnext[d]; cnext[d] += 1
                if k == cen[d]: continue
                f = cf[d]; g = nbr[f, k]
                if g < 0: continue
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
                    d += 1; cf[d] = g; cnext[d] = 0
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


'''
src = src[:a]+new+src[b:]
old = src[src.index("    vf_ptr, vf_idx = vertex_faces(F, len(P)); cap = 4*len(nodes)+1000"):src.index("        if not over: break")]
src = src.replace(old, "    cap = 4*len(nodes)+1000\n    while True:\n        sf, sp, sq, s0, s1, sl, epe, ept, ep3, fid, frs, over = _walk(P, F, Q, nf, node2, start, C2, nbr, fedge, E, depth, cap)\n")
src = src.replace("def embed(P, F, lines, log=None, budget=200):", "def embed(P, F, lines, log=None, depth=8):")
src = src.replace("(reasons: no aim 1, no exit 2, no face at corner 3, border 4, budget 5: {list(frs)})", "(no strip of up to {depth} triangles holds a straight stretch)")
open(path, 'w', encoding='utf-8').write(src); print('patched')
