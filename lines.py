"""Line cells -> chains, by merging, aligning and closing gaps (the construction test-flown in exp_merge.py).

Every line cell is a dot with a position (cells) and an unsigned direction (the line direction at that dot, from all
dots within reach). One cycle: merge dots that share a direction and lie beside each other; align each dot towards
the midpoint of its neighbours ahead and behind; close gaps between ends that point at each other. The chains of
linked dots are the lines.

run(P, F, res, reach=4, agree=25, beside=35, merge_reach=2, pull=0.5, gap=4, cycles=6, min_len=10) -> curves (n,3)
"""
import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components


def directions(X, N, reach):
    k = len(X); tree = cKDTree(X); pairs = tree.query_pairs(reach, output_type="ndarray")
    i, j = pairs[:, 0], pairs[:, 1]; v = X[j]-X[i]; v /= np.linalg.norm(v, axis=1, keepdims=True)+1e-12
    T = np.zeros((k, 3, 3))
    for a in (i, j): np.add.at(T, a, v[:, :, None]*v[:, None, :])
    w, e = np.linalg.eigh(T)
    d = e[:, :, -1]; d -= (d*N).sum(1, keepdims=True)*N; d /= np.linalg.norm(d, axis=1, keepdims=True)+1e-12
    return d


def unsign_mean(d1, d2):
    d2 = d2 if d1 @ d2 >= 0 else -d2; v = d1+d2; return v/(np.linalg.norm(v)+1e-12)


class P_:
    pass


def make_ops(reach_m, agree, beside, pull, gap):
    prm = P_(); prm.reach = reach_m; prm.cos_agree = np.cos(np.radians(agree)); prm.beside = beside; prm.cos_beside = np.cos(np.radians(55)); prm.pull = pull; prm.gap = gap; prm.link = 1.5   # links look one step along the chain; gaps handle the rest
    return prm


def merge(X, D, alive, prm):
    idx_alive = np.nonzero(alive)[0]; tree = cKDTree(X[idx_alive]); pairs = tree.query_pairs(prm.reach, output_type='ndarray')
    if len(pairs) == 0: return 0
    p, q = idx_alive[pairs[:, 0]], idx_alive[pairs[:, 1]]
    v = X[q]-X[p]; dist = np.linalg.norm(v, axis=1); u = v/(dist[:, None]+1e-12)
    par = np.abs((D[p]*D[q]).sum(1)) >= prm.cos_agree                     # same direction
    dm = D[p]+np.where((D[p]*D[q]).sum(1)[:, None] < 0, -D[q], D[q]); dm /= np.linalg.norm(dm, axis=1, keepdims=True)+1e-12
    along = np.abs((v*dm).sum(1)); lateral = np.sqrt(np.maximum(dist**2-along**2, 0))
    bes = lateral > along                                               # beside: sideways offset exceeds along offset (their perpendiculars cross); in-line dots never merge
    ok = par & bes; p, q, dist = p[ok], q[ok], dist[ok]
    order = np.argsort(dist); taken = np.zeros(len(X), bool); n = 0
    for a_, b_ in zip(p[order], q[order]):
        if taken[a_] or taken[b_]: continue
        X[a_] = 0.5*(X[a_]+X[b_]); D[a_] = unsign_mean(D[a_], D[b_]); alive[b_] = False; taken[a_] = taken[b_] = True; n += 1
    return n




def links(X, D, alive, prm):
    """for each live dot: nearest live dot ahead and behind along its direction within prm.reach (and beside-test passed as NOT beside)"""
    idx_alive = np.nonzero(alive)[0]; tree = cKDTree(X[idx_alive]); ahead = np.full(len(X), -1); behind = np.full(len(X), -1)
    for k in idx_alive:
        nb = idx_alive[tree.query_ball_point(X[k], prm.link)]; nb = nb[nb != k]
        if len(nb) == 0: continue
        v = X[nb]-X[k]; dist = np.linalg.norm(v, axis=1); u = v/(dist[:, None]+1e-12); along = u @ D[k]
        par = np.abs((D[nb]*D[k]).sum(1)) >= prm.cos_agree; lateral = dist*np.sqrt(np.maximum(1-along**2, 0))
        on_line = (lateral <= 0.5) & par                                   # chain neighbours are nearly in line; sideways dots are merge candidates
        fw = nb[on_line & (along > 0)]; bw = nb[on_line & (along < 0)]
        if len(fw): ahead[k] = fw[np.argmin(dist[on_line & (along > 0)])]
        if len(bw): behind[k] = bw[np.argmin(dist[on_line & (along < 0)])]
    return ahead, behind




def align(X, D, alive, ahead, behind, prm):
    """position: a dot with neighbours on both sides moves prm.pull of the way to their midpoint.
    direction: every linked dot takes the direction of its own chain, the chord through its neighbours (both sides
    when it has them, one side at an end), blended with its current direction; ends thus turn to face the chain."""
    both = alive & (ahead >= 0) & (behind >= 0); k = np.nonzero(both)[0]
    mid = 0.5*(X[ahead[k]]+X[behind[k]]); Xn = X.copy(); Xn[k] = X[k]+prm.pull*(mid-X[k])
    for kk in np.nonzero(alive & ((ahead >= 0) | (behind >= 0)))[0]:
        a_ = X[ahead[kk]] if ahead[kk] >= 0 else X[kk]; b_ = X[behind[kk]] if behind[kk] >= 0 else X[kk]
        chord = a_-b_; nrm = np.linalg.norm(chord)
        if nrm < 1e-9: continue
        D[kk] = unsign_mean(D[kk], chord/nrm)
    return Xn, len(k)


def close_gaps(X, D, alive, ahead, behind, prm):
    idx_alive = np.nonzero(alive)[0]; tree = cKDTree(X[idx_alive]); n = 0; newX, newD = [], []
    ends = [k for k in idx_alive if ahead[k] < 0 or behind[k] < 0]
    for k in ends:
        for sgn, slot in ((1, ahead), (-1, behind)):
            if slot[k] >= 0: continue
            probe = X[k]+sgn*D[k]*np.arange(1, int(prm.gap)+1)[:, None]
            for pp in probe:
                nb = idx_alive[tree.query_ball_point(pp, 1.0)]; nb = nb[nb != k]
                if len(nb) == 0: continue
                par = np.abs((D[nb]*D[k]).sum(1)) >= prm.cos_agree
                if not par.any(): continue
                t = nb[par][0]; dist = np.linalg.norm(X[t]-X[k]); dm = unsign_mean(D[k], D[t])
                m_ = int(np.ceil(dist))-1                     # bridge with dots at most one cell apart, so the link pass joins them this cycle
                for q in range(1, m_+1): newX.append(X[k]+(X[t]-X[k])*q/(m_+1)); newD.append(dm)
                n += 1; break
    if newX:
        X = np.vstack([X, np.array(newX)]); D = np.vstack([D, np.array(newD)]); alive = np.r_[alive, np.ones(len(newX), bool)]
    return X, D, alive, n




def chains(X, alive, ahead, behind, min_len):
    k = len(X); a = np.nonzero(alive & (ahead >= 0))[0]; b = np.nonzero(alive & (behind >= 0))[0]
    src_ = np.r_[a, b]; dst = np.r_[ahead[a], behind[b]]
    _, comp = connected_components(coo_matrix((np.ones(len(src_)), (src_, dst)), shape=(k, k)), directed=False)
    comp = np.where(alive, comp, -1); out = []
    for c in np.unique(comp[comp >= 0]):
        mem = np.nonzero(comp == c)[0]
        if len(mem) < min_len: continue
        ms = set(mem.tolist()); nbrs = {m: set() for m in mem}
        for m in mem:
            for t in (ahead[m], behind[m]):
                if t in ms: nbrs[m].add(t); nbrs[t].add(m)
        ends = [m for m in mem if len(nbrs[m]) <= 1]; start = ends[0] if ends else mem[0]
        path = [start]; seen = {start}; cur = start
        while True:
            nxt = [t for t in nbrs[cur] if t not in seen]
            if not nxt: break
            cur = min(nxt, key=lambda t: np.linalg.norm(X[t]-X[path[-1]])); path.append(cur); seen.add(cur)
        if len(path) >= min_len: out.append(X[path])
    return out


class State:
    """the chain construction in progress: dots X (cells), directions D, alive mask, links; kept between stages"""
    def __init__(self, X, D):
        self.X, self.D, self.alive = X, D, np.ones(len(X), bool); self.ahead = np.full(len(X), -1); self.behind = np.full(len(X), -1); self.cycles_done = 0


def prepare(res, reach=4.0):
    idx, line, N_all = res["idx"], res["line"], res["N"]
    li = np.nonzero(line)[0]; X = (idx[li]+0.5).astype(float)
    return State(X, directions(X, N_all[li], reach))


def refine(st, agree=25.0, beside=1.0, merge_reach=2.0, pull=0.5, gap=4.0, cycles=6, log=None):
    prm = make_ops(merge_reach, agree, beside, pull, gap)
    st.ahead, st.behind = links(st.X, st.D, st.alive, prm)
    for c in range(int(cycles)):
        nm = merge(st.X, st.D, st.alive, prm); st.ahead, st.behind = links(st.X, st.D, st.alive, prm); st.X, na = align(st.X, st.D, st.alive, st.ahead, st.behind, prm)
        st.ahead, st.behind = links(st.X, st.D, st.alive, prm); st.X, st.D, st.alive, ng = close_gaps(st.X, st.D, st.alive, st.ahead, st.behind, prm); st.ahead, st.behind = links(st.X, st.D, st.alive, prm)
        st.cycles_done += 1
        if log: log(f"  cycle {st.cycles_done}: merged {nm}, aligned {na}, gaps closed {ng} -> {st.alive.sum()} dots")
    return st


def finish(st, P, F, res, min_len=10, log=None):
    lo, size = res["lo"], res["size"]
    paths = chains(st.X, st.alive, st.ahead, st.behind, int(min_len))
    tri_c = P[F].mean(1); allp = np.vstack([P, tri_c]); stree = cKDTree(allp)
    curves = []
    for pth in paths:
        pts = lo+pth*size; _, jj = stree.query(pts); curves.append(allp[jj])
    if log: log(f"{len(curves)} chains of at least {min_len} dots, longest {max((len(c) for c in curves), default=0)}")
    return curves


def run(P, F, res, reach=4.0, agree=25.0, beside=1.0, merge_reach=2.0, pull=0.5, gap=4.0, cycles=6, min_len=10, log=None):
    st = prepare(res, reach); refine(st, agree, beside, merge_reach, pull, gap, cycles, log)
    curves = finish(st, P, F, res, min_len, log); run.groups = curves; return curves
