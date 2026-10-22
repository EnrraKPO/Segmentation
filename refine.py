"""Stage 3: refine a walked line in a tube of fine cubes (brief agreed 2026-10-21, evening).

The walk's line says where to look. Around it a tube of cubes about ten times finer than the walk's is laid: the
same cube stage as stage 1 (voxel.run_direct) restricted to the box around the line, with each cube's facing the blend of
all the triangles inside it, weighted by the area each covers in the cube (blend=True, 2026-10-22: a crease cube is no
longer one side or the other by where its centre falls); of those cubes only the ones within RADIUS of the line are kept.
On the tube each fine cube's bend is read over the surrounding area (bend): how fast the facing turns across a
neighbourhood AREA mm wide around the cube, with smooth weights over it, as degrees per UNIT mm (the coarse cube,
so the cut keeps its meaning). The cube next door decides nothing. A dot where that bend is at least the cut. The same anchored walk as stage 2 (trace4.run)
then runs on the fine dots, with its lengths in mm unchanged. The result is the line at the fine cube's resolution.
The fine walk keeps the walk's lengths in mm, with the settings agreed for the fine band (WALK): reach one step so it
never leaves the band, each dot's say its bend squared, trend pull equal to the dots'. End rule (avoid): the fine dots
within the area's reach (2*AREA) of the other walked lines our own ends sit on are removed before the walk, since
there the bend read mixes both creases and the band curls toward the neighbour; no other line counts.
Nothing else: no triangles are chosen, no planes, no curve fit."""
import os, sys, time, numpy as np
from scipy.spatial import cKDTree
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import voxel, trace4

DEFAULTS = dict(finer=3, radius=0.3, area=0.087, unit=0.087)
TUNE = dict(nbh=2.0, rungstep=0.5, rungw=0.5, cresth=0.5, tuberes=0.5, avoidk=2.0, avoidres=0.02, slack=1.0)   # the stage's inner constants (see TUNE_TIPS in app.py)
WALK = dict(step=0.087, slice=0.7, reach=0.087, pull=1.0, look=0.45, levp=2.0, minlen=8)   # the fine walk's settings (2026-10-21)


def resample(L, step):
    d = np.r_[0, np.cumsum(np.linalg.norm(np.diff(L, axis=0), axis=1))]
    if d[-1] < 2*step: return L.copy()
    t = np.arange(0, d[-1], step); return np.stack([np.interp(t, d, L[:, i]) for i in range(3)], 1)


def bend(C, N, area=0.087, unit=0.087, log=None, chunk=1500, at=None, nbh=2.0):
    """C cube centres (mm), N facings. Per cube: the rate of turn of the facing over the surrounding area (shape operator
    fitted to all neighbours within nbh*area, weighted exp(-(d/area)^2)), in degrees per unit mm, and the across direction."""
    t0 = time.time(); tree = cKDTree(C); k = len(C); turn = np.zeros(k); d1 = np.zeros((k, 3)); at = np.arange(k) if at is None else np.asarray(at)   # at: the cubes to read (all neighbours still count)
    a = np.where((np.abs(N[:, 0]) < 0.9)[:, None], np.array([1., 0, 0]), np.array([0, 1., 0]))
    t1 = np.cross(N, a); t1 /= np.linalg.norm(t1, axis=1, keepdims=True)+1e-12; t2 = np.cross(N, t1)
    for s0 in range(0, len(at), chunk):
        ii = at[s0:s0+chunk]; nb = tree.query_ball_point(C[ii], nbh*area)
        cnt = np.array([len(b) for b in nb]); i = np.repeat(ii, cnt); j = np.concatenate([np.asarray(b, int) for b in nb])
        dx = C[j]-C[i]; w = np.exp(-(np.linalg.norm(dx, axis=1)/area)**2); dx -= (dx*N[i]).sum(1, keepdims=True)*N[i]; dn = N[j]-N[i]
        ex, ey = (dx*t1[i]).sum(1), (dx*t2[i]).sum(1); nx, ny = (dn*t1[i]).sum(1), (dn*t2[i]).sum(1)
        loc = np.repeat(np.arange(len(ii)), cnt); acc = lambda v: np.bincount(loc, v*w, minlength=len(ii))
        A11, A12, A22 = acc(ex*ex), acc(ex*ey), acc(ey*ey); B11, B12, B21, B22 = acc(nx*ex), acc(nx*ey), acc(ny*ex), acc(ny*ey)
        det = A11*A22-A12*A12+1e-18; I11, I12, I22 = A22/det, -A12/det, A11/det
        S11 = B11*I11+B12*I12; S12 = B11*I12+B12*I22; S21 = B21*I11+B22*I12; S22 = B21*I12+B22*I22
        s12 = 0.5*(S12+S21); tr = S11+S22; dd = np.sqrt(np.maximum(((S11-S22)/2)**2+s12**2, 0)); l1, l2 = tr/2+dd, tr/2-dd
        k1 = np.where(np.abs(l1) >= np.abs(l2), l1, l2); turn[ii] = np.degrees(np.abs(k1))*unit
        vx = np.where(np.abs(s12) > 1e-12, s12, np.where(np.abs(S11-k1) < np.abs(S22-k1), 0.0, 1.0)); vy = np.where(np.abs(s12) > 1e-12, k1-S11, np.where(np.abs(S11-k1) < np.abs(S22-k1), 1.0, 0.0))
        nv = np.sqrt(vx*vx+vy*vy)+1e-18; d1[ii] = (vx/nv)[:, None]*t1[ii]+(vy/nv)[:, None]*t2[ii]
    turn, d1 = turn[at], d1[at]
    if log: log(f'bend over {area} mm at {len(at)} cubes: p50 {np.percentile(turn, 50):.1f} p90 {np.percentile(turn, 90):.1f} p99 {np.percentile(turn, 99):.1f} deg per {unit} mm, {time.time()-t0:.1f} s')
    return turn, d1


def tube(L, P, F, cube, height, finer=10, radius=0.3, cut=17.0, flat=True, log=None, workers=None, area=0.087, unit=0.087, rungs=True, tune=None):
    """fine cubes within `radius` mm of line L (mm), with their fold change and dot flag; same fields as creases.find_dots"""
    tn = {**TUNE, **(tune or {})}; size = cube/finer; cells = height/size; L = np.asarray(L, float)
    box = (L.min(0)-radius-size, L.max(0)+radius+size); t0 = time.time()
    r = voxel.run_direct(P, F, cells, log, height=height, flat=flat, box=box, workers=workers, blend=True)
    centre = r['lo']+(r['idx']+0.5)*r['size']; dist = cKDTree(resample(L, tn['tuberes']*size)).query(centre)[0]; keep = dist <= radius
    idx, N = r['idx'][keep], r['N'][keep]
    if log: log(f'tube: {len(idx)} fine cubes of {r["size"]:.4f} mm within {radius} mm of the line ({len(r["idx"])} in its box), {time.time()-t0:.1f} s')
    C = centre[keep]
    if rungs:                                                                   # read the bend only on rungs across the line, one per coarse step
        Lr = resample(L, tn['rungstep']*cube); T = np.gradient(Lr, axis=0); T /= np.linalg.norm(T, axis=1, keepdims=True)+1e-12
        dr, jr = cKDTree(Lr).query(C); along = ((C-Lr[jr])*T[jr]).sum(1); at = np.flatnonzero(np.abs(along) <= tn['rungw']*size)
        turn, d1 = bend(C, N, area, unit, log, at=at, nbh=tn['nbh']); idx, N, C = idx[at], N[at], C[at]
    else: turn, d1 = bend(C, N, area, unit, log, nbh=tn['nbh'])
    tree = cKDTree(C); h = tn['cresth']*unit                                              # crest: the bend peaks across the line (read half a unit either side)
    tp = turn[tree.query(C+h*d1)[1]]; tm = turn[tree.query(C-h*d1)[1]]; crest = (turn >= cut) & (turn >= tp) & (turn >= tm)
    band = turn >= cut
    if log: log(f'tube: fine dots {int(band.sum())} (bend >= {cut} deg), crest {int(crest.sum())}')
    return dict(idx=idx, N=N, lo=r['lo'], size=r['size'], turn=turn, line=band, crest=crest, d1=d1, cut=cut)


def end_lines(curves, i, slack, tree=None, owner=None):
    """the other walked lines line i's ends sit on: those with a node within `slack` mm of either end point. tree, owner: a
    cKDTree over all nodes and each node's line (nodes_tree), built once when many lines are asked"""
    if tree is None: tree, owner = nodes_tree(curves)
    L = curves[i]; near = set(owner[k] for e in (L[0], L[-1]) for k in tree.query_ball_point(e, slack)); near.discard(i)
    return sorted(near)


def nodes_tree(curves):
    return cKDTree(np.vstack(curves)), np.repeat(np.arange(len(curves)), [len(c) for c in curves])


def slim(t):
    """the tube reduced to what is kept per line: its dots (centres, mm) and their bend"""
    C = t['lo']+(t['idx'][t['line']]+0.5)*t['size']; return dict(dots=C.astype(np.float32), turn=t['turn'][t['line']].astype(np.float32), size=t['size'])


def avoid(t, others, reach, log=None, res=0.02):
    """drop the tube's dots within `reach` mm of the lines in `others` (mm): there the bend read mixes two creases"""
    if not others: return t
    C = t['lo']+(t['idx']+0.5)*t['size']; near = cKDTree(np.vstack([resample(c, res) for c in others])).query(C)[0] <= reach
    gone = int((t['line'] & near).sum()); t = dict(t); t['line'] = t['line'] & ~near; t['crest'] = t['crest'] & ~near
    if log: log(f'end rule: {gone} fine dots within {reach:.3f} mm of {len(others)} end line(s) removed, {int(t["line"].sum())} dots left')
    return t


def refine_line(L, P, F, cube, height, finer=3, radius=0.3, cut=17.0, flat=True, log=None, workers=None, walk=None, area=0.087, unit=0.087, others=(), tune=None):
    """returns (fine lines in mm from the walk on the tube's dots, the tube). others: the walked lines (mm) this line's ends sit on"""
    tn = {**TUNE, **(tune or {})}; t = tube(L, P, F, cube, height, finer, radius, cut, flat, log, workers, area, unit, tune=tn)
    t = avoid(t, list(others), tn['avoidk']*area, log, res=tn['avoidres'])
    lines = trace4.run(t, log=log or (lambda *_: None), workers=workers, **{**WALK, **(walk or {})})
    return lines, t


_G = {}
def _init(P, F): _G['P'] = P; _G['F'] = F
def _task(args):
    """one line on one core: no process pool inside (starting one per line cost ~5 s, the work is under a second)"""
    L, others, kw = args; lines, t = refine_line(L, _G['P'], _G['F'], log=None, workers=1, others=others, **kw); return lines, slim(t)


def refine_all(curves, P, F, cube, height, log=None, only=None, finer=3, radius=0.3, cut=17.0, flat=True, walk=None, area=0.087, unit=0.087, slack=None, workers=None, tune=None):
    """curves: the walk's lines (mm). only: index of one line, or None for all. slack: how near (mm) an end must be to another
    line for it to count as the line the end sits on (default TUNE slack coarse cubes). One line is refined here, with its log; many
    lines are spread over the cores, one line per core at a time. Returns dict(lines, index, tubes) with tubes slimmed to
    each line's dots (slim)."""
    idx = [only] if only is not None else list(range(len(curves))); t0 = time.time(); lines, tubes, index = [], [], []
    tn = {**TUNE, **(tune or {})}; slack = tn['slack']*cube if slack is None else slack; tree, owner = nodes_tree(curves)
    kw = dict(cube=cube, height=height, finer=finer, radius=radius, cut=cut, flat=flat, walk=walk, area=area, unit=unit, tune=tn)
    jobs = [(curves[i], [curves[j] for j in end_lines(curves, i, slack, tree, owner)], kw) for i in idx]
    workers = workers or os.cpu_count() or 1
    if len(jobs) == 1:
        fl, t = refine_line(jobs[0][0], P, F, log=log, workers=1, others=jobs[0][1], **kw); results = [(fl, slim(t))]
    else:
        from concurrent.futures import ProcessPoolExecutor
        results = []
        with ProcessPoolExecutor(min(workers, len(jobs)), initializer=_init, initargs=(P, F)) as ex:
            for n, r in enumerate(ex.map(_task, jobs, chunksize=4)):
                results.append(r)
                if log and (n+1) % 100 == 0: log(f'  refined {n+1} of {len(jobs)}, {time.time()-t0:.0f} s')
    for i, (fl, t) in zip(idx, results): lines += fl; tubes.append(t); index += [i]*len(fl)
    if log: log(f'refine done: {len(idx)} lines -> {len(lines)} fine lines, {time.time()-t0:.1f} s on {1 if len(jobs) == 1 else min(workers, len(jobs))} cores')
    return dict(lines=lines, index=index, tubes=tubes)


def nearest_line(curves, p):
    """index of the line with a point nearest p"""
    p = np.asarray(p, float); return int(np.argmin([np.linalg.norm(np.asarray(c)-p, axis=1).min() for c in curves]))
