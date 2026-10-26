"""one-off: replace _walk in regions.py by the version with corner hops"""
import os
root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); path = os.path.join(root, 'regions.py')
src = open(path, encoding='utf-8').read()
a = src.index('@njit(cache=True)\ndef _walk('); b = src.index('def embed(')
new = '''@njit(cache=True)
def _face_at_vertex(v, B, P, F, v0, e1, e2, C2, vf_ptr, vf_idx, exclude):
    """the face around vertex v whose wedge at v contains the direction towards B (projected into that face's plane);
    the best-fitting one otherwise. exclude: a face not to take (-1 for none)."""
    best = -1; bests = -1e300
    for idx in range(vf_ptr[v], vf_ptr[v+1]):
        h = vf_idx[idx]
        if h == exclude: continue
        kv = 0
        for c in range(3):
            if F[h, c] == v: kv = c
        dx = B[0]-v0[h, 0]; dy = B[1]-v0[h, 1]; dz = B[2]-v0[h, 2]
        du = dx*e1[h, 0]+dy*e1[h, 1]+dz*e1[h, 2]-C2[h, kv, 0]; dv = dx*e2[h, 0]+dy*e2[h, 1]+dz*e2[h, 2]-C2[h, kv, 1]
        dn = (du*du+dv*dv)**0.5
        if dn < 1e-300: continue
        du /= dn; dv /= dn
        ax = C2[h, (kv+1) % 3, 0]-C2[h, kv, 0]; ay = C2[h, (kv+1) % 3, 1]-C2[h, kv, 1]; an = (ax*ax+ay*ay)**0.5; ax /= an; ay /= an
        cx = C2[h, (kv+2) % 3, 0]-C2[h, kv, 0]; cy = C2[h, (kv+2) % 3, 1]-C2[h, kv, 1]; cn = (cx*cx+cy*cy)**0.5; cx /= cn; cy /= cn
        s = 1.0 if ax*cy-ay*cx >= 0 else -1.0
        score = min((ax*dv-ay*du)*s, (du*cy-dv*cx)*s)                 # both positive: inside the wedge
        if score > bests: bests = score; best = h
    return best


@njit(cache=True)
def _walk(P, F, nodes3, nface, node2, start, v0, e1, e2, C2, nbr, fedge, E, vf_ptr, vf_idx, budget, cap):
    """strokes between consecutive nodes of every line. Returns stroke arrays (face, p, q, edge point ids at p and q or
    -1, line), edge point arrays (edge id, canonical parameter along E[e,0]->E[e,1], 3D position), the nodes whose
    stretch to the next node could not be drawn through, overflow flag. A stretch that reaches a triangle corner hops to
    the face around that vertex which holds the direction onward (the corner is a ring node of every face around it)."""
    sf = np.empty(cap, np.int64); sp = np.empty((cap, 2)); sq = np.empty((cap, 2)); s0 = np.empty(cap, np.int64); s1 = np.empty(cap, np.int64); sl = np.empty(cap, np.int64)
    epe = np.empty(cap, np.int64); ept = np.empty(cap); ep3 = np.empty((cap, 3)); fid = np.empty(cap, np.int64)
    ns = 0; ne = 0; fails = 0; over = False
    for li in range(len(start)-1):
        for i in range(start[li], start[li+1]-1):
            if ns+2*budget+2 >= cap or ne+budget+2 >= cap: over = True; break
            b = i+1; f = nface[i]; fb = nface[b]; B = nodes3[b]
            pu = node2[i, 0]; pv = node2[i, 1]; prev = -1; entry = -1; fresh = True; ok = False
            for step in range(budget):
                if f == fb:
                    qu = node2[b, 0]; qv = node2[b, 1]
                    if (qu-pu)**2+(qv-pv)**2 > 1e-24:
                        sf[ns] = f; sp[ns, 0] = pu; sp[ns, 1] = pv; sq[ns, 0] = qu; sq[ns, 1] = qv; s0[ns] = prev; s1[ns] = -1; sl[ns] = li; ns += 1
                    ok = True; break
                if fresh:                                                             # a node sitting on a corner: start from the right face around it
                    for k in range(3):
                        if (C2[f, k, 0]-pu)**2+(C2[f, k, 1]-pv)**2 < 1e-18:
                            v = F[f, k]; h = _face_at_vertex(v, B, P, F, v0, e1, e2, C2, vf_ptr, vf_idx, -1)
                            if h >= 0 and h != f:
                                for c in range(3):
                                    if F[h, c] == v: pu = C2[h, c, 0]; pv = C2[h, c, 1]
                                f = h
                            fresh = False; break
                    if f == fb: continue
                dx = B[0]-v0[f, 0]; dy = B[1]-v0[f, 1]; dz = B[2]-v0[f, 2]
                du = dx*e1[f, 0]+dy*e1[f, 1]+dz*e1[f, 2]-pu; dv = dx*e2[f, 0]+dy*e2[f, 1]+dz*e2[f, 2]-pv
                dn = (du*du+dv*dv)**0.5
                if dn < 1e-300: fid[fails] = i; fails += 1; break
                bs = 1e300; bk = -1; bt = 0.0
                for rnd in range(2):                                                  # second round: the entry edge is allowed (the aim turned back)
                    for k in range(3):
                        if rnd == 0 and k == entry: continue
                        k1 = (k+1) % 3
                        cx = C2[f, k, 0]; cy = C2[f, k, 1]; ex = C2[f, k1, 0]-cx; ey = C2[f, k1, 1]-cy
                        den = du*ey-dv*ex
                        if abs(den) < 1e-300: continue
                        rx = cx-pu; ry = cy-pv
                        s = (rx*ey-ry*ex)/den; t = (rx*dv-ry*du)/den
                        smin = -1e-9/dn if fresh else 1e-9/dn                           # a fresh node may sit on the edge it leaves through; afterwards the exit must be a step away
                        if s > smin and t >= -1e-9 and t <= 1+1e-9 and s < bs: bs = s; bk = k; bt = t
                    if bk >= 0: break
                fresh = False
                if bk < 0: fid[fails] = i; fails += 1; break
                k = bk; t = min(max(bt, 0.0), 1.0); k1 = (k+1) % 3
                if t < 1e-9 or t > 1-1e-9:                                            # through a corner: hop to the face that holds the direction onward
                    kc = k if t < 1e-9 else k1; v = F[f, kc]; xu = C2[f, kc, 0]; xv = C2[f, kc, 1]
                    if (xu-pu)**2+(xv-pv)**2 > 1e-24:
                        sf[ns] = f; sp[ns, 0] = pu; sp[ns, 1] = pv; sq[ns, 0] = xu; sq[ns, 1] = xv; s0[ns] = prev; s1[ns] = -1; sl[ns] = li; ns += 1
                    h = _face_at_vertex(v, B, P, F, v0, e1, e2, C2, vf_ptr, vf_idx, f)
                    if h < 0: fid[fails] = i; fails += 1; break
                    for c in range(3):
                        if F[h, c] == v: pu = C2[h, c, 0]; pv = C2[h, c, 1]
                    f = h; entry = -1; prev = -1; continue
                e = fedge[f, k]; va = F[f, k]; vb = F[f, k1]
                tcan = t if va == E[e, 0] else 1.0-t
                xu = C2[f, k, 0]+t*(C2[f, k1, 0]-C2[f, k, 0]); xv = C2[f, k, 1]+t*(C2[f, k1, 1]-C2[f, k, 1])
                epe[ne] = e; ept[ne] = tcan
                for c in range(3): ep3[ne, c] = P[va, c]+t*(P[vb, c]-P[va, c])
                if (xu-pu)**2+(xv-pv)**2 > 1e-24:
                    sf[ns] = f; sp[ns, 0] = pu; sp[ns, 1] = pv; sq[ns, 0] = xu; sq[ns, 1] = xv; s0[ns] = prev; s1[ns] = ne; sl[ns] = li; ns += 1
                g = nbr[f, k]
                if g < 0: ne += 1; fid[fails] = i; fails += 1; break
                kk = 0
                for c in range(3):
                    if fedge[g, c] == e: kk = c
                kk1 = (kk+1) % 3
                tg = tcan if F[g, kk] == E[e, 0] else 1.0-tcan
                pu = C2[g, kk, 0]+tg*(C2[g, kk1, 0]-C2[g, kk, 0]); pv = C2[g, kk, 1]+tg*(C2[g, kk1, 1]-C2[g, kk, 1])
                entry = kk; prev = ne; ne += 1; f = g
            else:
                fid[fails] = i; fails += 1
        if over: break
    return sf[:ns], sp[:ns], sq[:ns], s0[:ns], s1[:ns], sl[:ns], epe[:ne], ept[:ne], ep3[:ne], fid[:fails], over


def vertex_faces(F, nv):
    """CSR: faces around each vertex"""
    order = np.argsort(F.ravel(), kind='stable'); cnt = np.bincount(F.ravel(), minlength=nv)
    return np.r_[0, np.cumsum(cnt)].astype(np.int64), (order//3).astype(np.int64)


'''
src = src[:a]+new+src[b:]
old = "    cap = 4*len(nodes)+1000\n    while True:\n        sf, sp, sq, s0, s1, sl, epe, ept, ep3, fid, over = _walk(P, F, Q, nf, node2, start, v0, e1, e2, C2, nbr, fedge, E, budget, cap)"
assert old in src
src = src.replace(old, "    vf_ptr, vf_idx = vertex_faces(F, len(P)); cap = 4*len(nodes)+1000\n    while True:\n        sf, sp, sq, s0, s1, sl, epe, ept, ep3, fid, over = _walk(P, F, Q, nf, node2, start, v0, e1, e2, C2, nbr, fedge, E, vf_ptr, vf_idx, budget, cap)")
open(path, 'w', encoding='utf-8').write(src); print('patched')
