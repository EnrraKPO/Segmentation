"""Stage 4: bridge the pieces and seek the joints (brief agreed 2026-10-23).

The refined pieces sit on the crease but break inside one crease and stop short of junctions (the fine walk's anchor
window, not the measurement). No node is moved here. At each piece end the trajectory is read and continued:

  Trajectory: a curve of constant shape, a helix (one bend rate, one twist rate: a straight line, a circle and a curve
     round a curved body are all helices), fitted to the end's last segment. The segment is found by growing back from
     the end while the nodes' scatter about the fitted helix (rms) stays within their own wobble (WOBBLE, measured on the
     refined pieces); the first node that breaks that marks where the bend changed. The helix continued forward is the trajectory. (The
     parabola over a fixed fade, the first version, stays selectable: MODEL.)
  A  continuation: consecutive pieces of one walk line face each other across a gap. Both ends' trajectories are
     continued; where they meet within TOL_A and agree in direction within ANG_A the gap is bridged by the smooth blend
     of the two. Refused bridges leave both ends loose.
  B  joint: every loose end of every line (the walk's record of what the end sat on is NOT trusted, every line is a
     candidate) is continued up to MAX_B. The nearest thing its trajectory reaches within TOL_B wins: the body of another
     line (a tee: the bridge is pulled onto that line's nearest node, which both lines then share) or another loose
     end's trajectory (a corner: both bridges meet at one shared point). An end with nothing in reach stays loose.

Lengths in mm. run(pieces, index, curves, ...) -> dict(lines, measured, walk, shared, loose, report)."""
import numpy as np
from scipy.spatial import cKDTree

DEFAULTS = dict(model='helix', wobble=0.01, minfit=5, fade=0.5, step=0.087, tol_a=0.1, ang_a=30.0, max_a=3.0, tol_b=0.1, max_b=1.5)


def smoothstep(s):
    s = np.clip(s, 0, 1); return s*s*(3-2*s)


def helix_points(prm, t):
    """points of the helix prm at lengths t along it (0 = the end, negative back along the piece). prm = (r0, T, N, k, tau):
    the end point, the unit direction there, the unit normal (toward the bend's centre), the bend rate, the twist rate."""
    r0, T, N, k, tau = prm; B = np.cross(T, N)
    t = np.asarray(t, float); w = np.hypot(k, tau); wt = w*t
    small = w < 1e-6; wsafe = np.where(small, 1.0, w)
    s1 = np.where(small, t, np.sin(wt)/wsafe)                                                # sin(wt)/w -> t
    c1 = np.where(small, t*t/2, (1-np.cos(wt))/(wsafe*wsafe))                                # (1-cos wt)/w^2 -> t^2/2
    s3 = np.where(small, t**3/6, (t-s1)/(wsafe*wsafe))                                       # (t - sin(wt)/w)/w^2 -> t^3/6
    along = np.where(small, t, (tau*tau*t+k*k*s1)/(wsafe*wsafe))
    return r0+np.multiply.outer(along, T)+np.multiply.outer(k*c1, N)+np.multiply.outer(k*tau*s3, B)


def _frame(th, ph, ps):
    """unit direction T from two angles, unit normal N from one angle in the plane across T"""
    T = np.array([np.cos(th)*np.cos(ph), np.cos(th)*np.sin(ph), np.sin(th)])
    a = np.array([0.0, 0.0, 1.0]) if abs(T[2]) < 0.9 else np.array([1.0, 0.0, 0.0]); u = np.cross(T, a); u /= np.linalg.norm(u); v = np.cross(T, u)
    return T, np.cos(ps)*u+np.sin(ps)*v


def _unpack(x):
    T, N = _frame(x[3], x[4], x[5]); return x[0:3], T, N, float(x[6]), float(x[7])


def _angles(T):
    return np.arcsin(np.clip(T[2], -1, 1)), np.arctan2(T[1], T[0])


def _init(P, t):
    """starting helix: the end point, the direction of the last chord, the circle through the first, middle and last point"""
    d = P[-1]-P[max(len(P)-3, 0)]; d /= np.linalg.norm(d)+1e-12; th, ph = _angles(d); x = np.array([*P[-1], th, ph, 0.0, 0.0, 0.0])
    if len(P) >= 3:
        a, b, c = P[0], P[len(P)//2], P[-1]; ab, ac = b-a, c-a; n = np.cross(ab, ac); nn = np.linalg.norm(n)
        if nn > 1e-12:
            k = 2*nn/(np.linalg.norm(ab)*np.linalg.norm(ac)*np.linalg.norm(c-b)+1e-18)          # curvature of the circle through three points
            T, _ = _frame(th, ph, 0.0); Nc = np.cross(n/nn, T); Nc -= (Nc @ T)*T; Nc /= np.linalg.norm(Nc)+1e-12   # normal toward the centre
            _, u = _frame(th, ph, 0.0); v = np.cross(T, u); x[5] = np.arctan2(Nc @ v, Nc @ u); x[6] = k
    return x


def fit_helix(P, t, x0=None, iters=40):
    """least-squares helix through points P at chain lengths t (<= 0, 0 at the end): Levenberg-Marquardt with a finite-difference
    Jacobian over 8 parameters (end point, two direction angles, one normal angle, bend, twist). Returns the parameter
    vector and the largest distance of a point from the helix."""
    x = _init(P, t) if x0 is None else x0.copy(); lam = 1e-2
    def resid(x): return (helix_points(_unpack(x), t)-P).ravel()
    r = resid(x); cost = r @ r
    for _ in range(iters):
        J = np.empty((len(r), len(x)))
        for j in range(len(x)):
            h = 1e-6*max(1.0, abs(x[j])); xp = x.copy(); xp[j] += h; J[:, j] = (resid(xp)-r)/h
        A = J.T @ J; g = J.T @ r; D = np.diag(np.maximum(np.diag(A), 1e-3*np.diag(A).max()+1e-12)); moved = False
        for _ in range(10):
            dx = np.linalg.solve(A+lam*D, -g); xn = x+dx; rn = resid(xn); cn = rn @ rn
            if cn < cost: gain = (cost-cn)/(cost+1e-30); x, r, cost = xn, rn, cn; lam = max(lam/3, 1e-9); moved = True; break
            lam *= 10
        if not moved or gain < 1e-6: break
    far = float(np.sqrt((np.linalg.norm(helix_points(_unpack(x), t)-P, axis=1)**2).mean()))   # the nodes' scatter about the helix (rms)
    return x, far


def fit_end_helix(L, end, wobble, minfit):
    """the end's trajectory: a helix fitted to the last segment, found by growing back from the end while the nodes' rms scatter
    about the helix stays within `wobble`. Returns f(t) (t mm beyond the end), the segment length in nodes, the largest residual."""
    P = np.asarray(L if end else L[::-1], float); seg = np.linalg.norm(np.diff(P, axis=0), axis=1); s = np.concatenate([[0.0], np.cumsum(seg)]); t = s-s[-1]
    n = len(P); m0 = min(max(int(minfit), 4), n)
    def fit(m, x0=None): return fit_helix(P[n-m:], t[n-m:], x0)
    ok_x, far = fit(m0); ok_m = m0
    if far <= wobble and n > m0:
        m = m0; bad_m = None
        while True:                                                                  # grow by doubling, then bisect the boundary
            m2 = min(2*m, n); x2, far2 = fit(m2, ok_x)
            if far2 <= wobble: ok_m, ok_x, m = m2, x2, m2
            else: bad_m = m2; break
            if m2 == n: break
        if bad_m is not None:
            lo, hi = ok_m, bad_m
            while hi-lo > 1:
                mid = (lo+hi)//2; xm, fm = fit(mid, ok_x)
                if fm <= wobble: lo, ok_x = mid, xm
                else: hi = mid
            ok_m = lo
    prm = _unpack(ok_x); far = float(np.sqrt((np.linalg.norm(helix_points(prm, t[n-ok_m:])-P[n-ok_m:], axis=1)**2).mean()))
    fit_end_helix.prm = prm
    return (lambda tt: helix_points(prm, tt)), ok_m, far


def fit_end(L, end, fade):
    """curve through the last stretch of line L at `end` (0 = first node, 1 = last): a quadratic in the length along the
    chain, fitted by least squares with every node weighted exp(-(distance back / fade)^2). Returns f(t): the point t mm
    beyond the end (t >= 0 extends, t < 0 runs back along the fit)."""
    P = L if end else L[::-1]; P = np.asarray(P, float)
    seg = np.linalg.norm(np.diff(P, axis=0), axis=1); s = np.concatenate([[0.0], np.cumsum(seg)]); t = s-s[-1]   # 0 at the end, negative back
    if len(P) < 3 or s[-1] < 1e-9:
        d = (P[-1]-P[0]); d /= np.linalg.norm(d)+1e-12; e = P[-1]; return lambda tt: e+np.multiply.outer(np.asarray(tt, float), d)
    w = np.exp(-(t/fade)**2); A = np.stack([np.ones(len(P)), t, t*t], 1)*w[:, None]
    a0, b, c = np.linalg.lstsq(A, P*w[:, None], rcond=None)[0]
    return lambda tt: a0+np.multiply.outer(np.asarray(tt, float), b)+np.multiply.outer(np.asarray(tt, float)**2, c)


def trajectory_spec(L, end, p):
    """the trajectory at an end of line L under the settings p, as data: ('helix', prm) over the grown segment, or ('para', coef)
    over the fade"""
    if p.get('model', 'helix') == 'helix':
        P = np.asarray(L if end else L[::-1], float); seg = np.linalg.norm(np.diff(P, axis=0), axis=1); t = np.r_[0.0, np.cumsum(seg)]; t -= t[-1]
        _, m, _ = fit_end_helix(L, end, p['wobble'], p['minfit']); return ('helix', fit_end_helix.prm, m)
    P = np.asarray(L if end else L[::-1], float); seg = np.linalg.norm(np.diff(P, axis=0), axis=1); t = np.r_[0.0, np.cumsum(seg)]; t -= t[-1]
    if len(P) < 3 or t[0] > -1e-9: d = P[-1]-P[0]; d /= np.linalg.norm(d)+1e-12; return ('para', (P[-1], d, np.zeros(3)), len(P))
    w = np.exp(-(t/p['fade'])**2); A = np.stack([np.ones(len(P)), t, t*t], 1)*w[:, None]; return ('para', tuple(np.linalg.lstsq(A, P*w[:, None], rcond=None)[0]), len(P))


def make_f(spec):
    kind, prm = spec[0], spec[1]
    if kind == 'helix': return lambda tt: helix_points(prm, tt)
    a0, b, c = prm; return lambda tt: a0+np.multiply.outer(np.asarray(tt, float), b)+np.multiply.outer(np.asarray(tt, float)**2, c)


def trajectory(L, end, p): return make_f(trajectory_spec(L, end, p))


def _spec_task(args):
    L, end, p = args; return trajectory_spec(L, end, p)


def trajectories(items, p, workers=None):
    """trajectory specs for many (line, end) pairs, spread over the cores"""
    import os
    from concurrent.futures import ProcessPoolExecutor
    workers = workers or os.cpu_count() or 1; tasks = [(L, e, p) for L, e in items]
    if workers > 1 and len(tasks) > 64:
        with ProcessPoolExecutor(workers) as ex: return list(ex.map(_spec_task, tasks, chunksize=max(1, len(tasks)//(workers*8))))
    return [_spec_task(t) for t in tasks]


def _tangent(f, t, h=1e-3):
    d = f(t+h)-f(t-h); return d/(np.linalg.norm(d)+1e-12)


def _blend_bridge(fa, ta, fb, tb, step):
    """A: interior points of the bridge from end a (extension fa, meeting at ta) to end b (fb, meeting at tb): the two
    extensions blended by a smooth share of the way across, both reaching the meeting point at the same share"""
    T = ta+tb; n = max(int(np.ceil(T/step)), 1)
    s = np.arange(1, n)/n
    if len(s) == 0: return np.zeros((0, 3))
    w = smoothstep(s)[:, None]; return (1-w)*fa(s*T)+w*fb((1-s)*T)


def _pulled_bridge(f, t, target, step):
    """B: points of the bridge along the extension f up to length t, pulled smoothly onto `target` (the last point IS the target)"""
    n = max(int(np.ceil(t/step)), 1); tt = np.arange(1, n+1)/n*t
    pts = f(tt); off = target-f(t); return pts+smoothstep(tt/t)[:, None]*off


def continuation(pieces, index, curves, p, log=None, workers=None):
    """rule A. pieces: refined lines (mm), index[k]: the walk line each came from, curves: the walk lines. Returns the
    joined lines, per-node measured flags, each line's walk index, and the counts."""
    step, tol_a, ang_a, max_a = p['step'], p['tol_a'], p['ang_a'], p['max_a']
    cos_a = np.cos(np.radians(ang_a)); out, meas, widx = [], [], []; n_ok = n_ref = n_touch = 0
    by = {}
    for k, i in enumerate(index): by.setdefault(int(i), []).append(k)
    items = [(np.asarray(pieces[k], float), e) for k in range(len(pieces)) for e in (0, 1)]; specs = trajectories(items, p, workers)
    F = {}                                                                        # (piece, end) -> trajectory function; a piece reversed swaps its ends
    for q, (k, e) in enumerate((k, e) for k in range(len(pieces)) for e in (0, 1)): F[(k, e)] = make_f(specs[q])
    for i, ks in by.items():
        L = curves[i]; t = cKDTree(L); ordered = []
        for k in ks:                                                             # orient each piece along its walk line and order them
            q = np.asarray(pieces[k], float); a, b = t.query(q[0])[1], t.query(q[-1])[1]; fs, fe = F[(k, 0)], F[(k, 1)]
            if a > b: q = q[::-1]; a, b = b, a; fs, fe = fe, fs
            ordered.append((a, b, q, fs, fe))
        ordered.sort(key=lambda x: (x[0], x[1]))
        cur = ordered[0][2]; cm = np.ones(len(cur), bool); cur_fe = ordered[0][4]
        for a, b, q, fs, fe in ordered[1:]:
            gap = np.linalg.norm(q[0]-cur[-1])
            if gap <= step:                                                      # touching: just continue
                cur = np.vstack([cur, q]); cm = np.r_[cm, np.ones(len(q), bool)]; cur_fe = fe; n_touch += 1; continue
            if gap > max_a: out.append(cur); meas.append(cm); widx.append(i); cur, cm, cur_fe = q, np.ones(len(q), bool), fe; n_ref += 1; continue
            fa = cur_fe; fb = fs                                                    # the trajectory out of the joined line's end is its last piece's
            ts = np.arange(0, gap+2*step, step*0.5); A = fa(ts); B = fb(ts)        # where the two extensions pass closest
            ia, ib = np.unravel_index(np.argmin(((A[:, None, :]-B[None, :, :])**2).sum(2)), (len(ts), len(ts)))
            ta, tb = ts[ia], ts[ib]; dist = np.linalg.norm(A[ia]-B[ib]); agree = -(_tangent(fa, ta) @ _tangent(fb, tb))
            if dist <= tol_a and agree >= cos_a and ta+tb > 0:
                br = _blend_bridge(fa, ta, fb, tb, step); cur = np.vstack([cur, br, q]); cm = np.r_[cm, np.zeros(len(br), bool), np.ones(len(q), bool)]; cur_fe = fe; n_ok += 1
            else: out.append(cur); meas.append(cm); widx.append(i); cur, cm, cur_fe = q, np.ones(len(q), bool), fe; n_ref += 1
        out.append(cur); meas.append(cm); widx.append(i)
    if log: log(f'A continuation: {len(pieces)} pieces of {len(by)} walk lines -> {len(out)} lines; bridged {n_ok}, touching {n_touch}, refused {n_ref}')
    return out, meas, widx, dict(a_bridged=n_ok, a_touch=n_touch, a_refused=n_ref)


def joints(lines, meas, p, log=None, workers=None):
    """rule B on every end of every line. Returns the lines with bridges appended, measured flags, shared nodes
    [(line a, line b, point)], loose ends [(line, end, point)] and per-end records [(line, end, kind, target line, length)]."""
    step, tol_b, max_b = p['step'], p['tol_b'], p['max_b']
    lines = [np.asarray(l, float) for l in lines]; meas = [np.asarray(m, bool) for m in meas]; n = len(lines)
    allp = np.vstack(lines); owner = np.repeat(np.arange(n), [len(l) for l in lines]); tree = cKDTree(allp)
    ends = [(i, e) for i in range(n) for e in (0, 1) if len(lines[i]) >= 2]
    tt = np.arange(step, max_b+step*0.5, step); F = {}; S = []; sid = []
    specs = trajectories([(lines[i], e) for i, e in ends], p, workers)
    for q, (i, e) in enumerate(ends):
        f = make_f(specs[q]); F[(i, e)] = f; pts = f(tt); S.append(pts); sid.append(np.full(len(tt), q))
    S = np.vstack(S); sid = np.concatenate(sid); st = np.tile(tt, len(ends)); stree = cKDTree(S)
    cand = []                                                                     # (length, kind, end q, target)
    for q, (i, e) in enumerate(ends):
        pts = S[sid == q]
        dd, j = tree.query(pts); hit = np.flatnonzero((dd <= tol_b) & (owner[j] != i))     # body of another line
        if len(hit): h = hit[0]; cand.append((tt[h], 'tee', q, int(j[h])))
        near = stree.query_ball_point(pts, tol_b)                                          # another end's extension
        best = None
        for h, lst in enumerate(near):
            for m in lst:
                q2 = int(sid[m])
                if q2 == q or ends[q2][0] == i: continue
                tot = tt[h]+st[m]
                if best is None or tot < best[0]: best = (tot, q2, tt[h], st[m])
        if best is not None: cand.append((best[2], 'corner', q, best))
    cand.sort(key=lambda c: c[0]); used = set(); shared = []; rec = {}; add = {q: None for q in range(len(ends))}
    for length, kind, q, tgt in cand:
        if q in used: continue
        i, e = ends[q]
        if kind == 'tee':
            target = allp[tgt]; br = _pulled_bridge(F[(i, e)], length, target, step); add[q] = br; used.add(q)
            shared.append((i, int(owner[tgt]), target)); rec[q] = ('tee', int(owner[tgt]), length)
        else:
            tot, q2, ta, tb = tgt
            if q2 in used: continue
            i2, e2 = ends[q2]; M = 0.5*(F[(i, e)](ta)+F[(i2, e2)](tb))
            add[q] = _pulled_bridge(F[(i, e)], ta, M, step); add[q2] = _pulled_bridge(F[(i2, e2)], tb, M, step); used.update((q, q2))
            shared.append((i, i2, M)); rec[q] = ('corner', i2, ta); rec[q2] = ('corner', i, tb)
    out, om, loose, records = [], [], [], []; qof = {k: q for q, k in enumerate(ends)}
    for i in range(n):
        L, m = lines[i], meas[i]
        for e in (1, 0):
            q = qof.get((i, e))
            if q is None: continue
            br = add[q]
            if br is None: loose.append((i, e, L[-1] if e else L[0])); records.append((i, e, 'loose', -1, 0.0)); continue
            records.append((i, e)+rec[q])
            if e: L = np.vstack([L, br]); m = np.r_[m, np.zeros(len(br), bool)]
            else: L = np.vstack([br[::-1], L]); m = np.r_[np.zeros(len(br), bool), m]
        out.append(L); om.append(m)
    kinds = [r[2] for r in records]
    if log: log(f'B joints: {len(ends)} ends; tees {kinds.count("tee")}, corners {kinds.count("corner")}, loose {kinds.count("loose")}; shared nodes {len(shared)}')
    return out, om, shared, loose, records


def run(pieces, index, curves, log=None, workers=None, **kw):
    """pieces: refined lines (mm) with index[k] = walk line of each; curves: all walk lines. Walk lines with no piece
    take part as they are. Returns dict(lines, measured (per-node: True = a measured node, False = bridge), walk (walk
    index per line), shared, loose, report)."""
    p = {**DEFAULTS, **kw}
    lines, meas, widx, rep = continuation(pieces, index, curves, p, log, workers)
    done = set(int(i) for i in index)
    for i, c in enumerate(curves):
        if i not in done and len(c) >= 2: lines.append(np.asarray(c, float)); meas.append(np.ones(len(c), bool)); widx.append(i)
    lines, meas, shared, loose, records = joints(lines, meas, p, log, workers)
    rep.update(tees=sum(r[2] == 'tee' for r in records), corners=sum(r[2] == 'corner' for r in records), loose=len(loose), lines=len(lines), bridge_nodes=int(sum((~m).sum() for m in meas)))
    return dict(lines=lines, measured=meas, walk=widx, shared=shared, loose=loose, records=records, report=rep)


def walk_agreement(result, curves, slack):
    """how the B joints compare with the walk's own record: for each bridged end, did the walk's end sit on the line joined
    to (agree), on another line (disagree), or on none (new)?"""
    import refine
    tree, owner = refine.nodes_tree(curves); widx = result['walk']; agree = disagree = new = 0
    for i, e, kind, tgt, _ in result['records']:
        if kind == 'loose' or widx[i] < 0 or tgt < 0: continue
        near = set(refine.end_lines(curves, widx[i], slack, tree, owner))
        if not near: new += 1
        elif widx[tgt] in near: agree += 1
        else: disagree += 1
    return dict(agree=agree, disagree=disagree, new=new)
