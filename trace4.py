"""The anchored walk (docs/brief_anchored_walk.md). Lines through the dots, one node per slice across the band.

Per line: seed = strongest free dot. At an anchor dot: axis from the dots in a sphere (no preferred direction);
band = the hill of the across-profile the anchor sits on; node = weighted mean of the band's dots in a thin slice;
next anchor = strongest dot about one step ahead along the trend of the last nodes, inside the band's width.
Ends: no dot ahead; the new node grazes a node of another line (it becomes that node, shared by both); uncertain axis.
Dot ownership only blocks seeds; it never ends a line.

usage: python trace4.py grid.npz out.npz        env (mm unless noted): STEP 0.087, SLICE 0.7 (half-length of the
slice along the axis), SPHERE 0.7 (radius for the axis read, grows to 3x the band width), WMAX 1.0 (widest band
half-width looked at), ALIGN 1.5 (axis spread ratio long/wide), TREND 6 (nodes), CONE 45 (deg), MINLEN 8 (nodes),
PULL 10 (A: weight of the trend's predicted point against the dots' total, 0 = dots only; LOOK 0.35 probes ahead and scales it: dot support
along the trend vs along the local axis), LEVP 6 (dot say = (strength/line standard)^LEVP, LEVN 6 nodes of memory), AFALL 2.0 (B: cubes off the trend's
centre line at which a candidate's score has fallen to 1/e), ASTR 0 (B: exponent on relative strength; 0 = alignment only), OWNED 0.5 (fraction of the seed band owned by another line that blocks the seed), VALLEY 1.3 (the profile must rise
this much above its running minimum to count as a second hill)"""
import os, sys, time, numpy as np
from scipy.spatial import cKDTree


def env(k, d): return float(os.environ.get(k, d))


def line_power(wgt_claimed, n_nodes):
    """Power of a line: the strength (degrees above the cut) of every dot it gathered, divided by its node count.
    A dense row of strong dots scores high; a sparse row low; a scatter lowest. One place to change."""
    return float(wgt_claimed.sum())/max(n_nodes, 1)


def walk(X, wgt, N, size, log=print, **kw):
    """X dots (cubes), wgt strengths, N facings. Returns list of node arrays (cubes) and the owner array."""
    step = kw.get('step', 0.087)/size; slice_h = kw.get('slice', 0.7)/size; sph0 = kw.get('sphere', 0.7)/size
    wmax = kw.get('wmax', 1.0)/size; align = kw.get('align', 1.5); trend_n = int(kw.get('trend', 6)); cos_cone = np.cos(np.radians(kw.get('cone', 45)))
    reach = kw.get('reach', 0.22)/size   # how far ahead (cubes) the next anchor may lie: a gap in the dot row shorter than this does not end the line
    along = kw.get('along', 0.0)/size or reach   # ALONG: falloff (mm) of the along-score around one step; default = reach (the two roles coupled, as before)
    pull = kw.get('pull', 10.0); afall = kw.get('afall', 2.0); astr = kw.get('astr', 0.0)   # A: trend weight on the node; B: score falloff (cubes) off the centre line, strength exponent
    cos_bend = np.cos(np.radians(kw.get('bend', 45.0))); tan_bend = np.tan(np.radians(min(kw.get('bend', 45.0), 89.0)))   # BEND (deg): the dots' direction may steer the node only within this angle of the heading
    powercut = kw.get('powercut', 3.0)   # POWERCUT: a line whose power (line_power) is below this never existed; its dots are released
    graze = kw.get('graze', 1.0)   # GRAZE: a node within this many cubes of another line's node lands on that node (shared); the line ends there
    minlen = int(kw.get('minlen', 8)); owned_frac = kw.get('owned', 0.5); valley = kw.get('valley', 1.3); core_frac = kw.get('core', 0.0); grow = kw.get('grow', 99.0); lead_k = kw.get('lead', 0.0); lead_min = int(kw.get('leadmin', 8)); fringe = kw.get('fringe', 1.0) > 0; ridge = kw.get('ridge', 0.0) > 0; sparse = kw.get('sparse', 0.0); hill = kw.get('hill', 0.0) > 0
    walk.probe = []
    tree = cKDTree(X); n = len(X); owner = np.full(n, -1); lines = []; node_tree = None; node_owner = np.zeros(0, int); node_index = np.zeros(0, int); all_nodes = np.zeros((0, 3))
    look = kw.get('look', 0.35)/size; lookw = kw.get('lookw', 1.0); lookn = int(kw.get('lookn', 0)); levp = kw.get('levp', 6.0); levn = kw.get('levn', 6.0)   # LEVP: power of the strength/standard curve (0 = off); LEVN: nodes of memory of the standard  # LOOK: how far ahead (mm) the probe corridors run, 0 = off; LOOKW: corridor half-width (cubes)
    joins = []; powers = []   # joins: (line, node, other line, other node) shared nodes; powers: line_power per kept line
    stats = dict(no_anchor=0, collision=0, self_hit=0, uncertain=0, short=0, trend_used=0, clipped=0, look_trend=0, look_dots=0, look_mixed=0, bend_held=0, weak=0); kept_ends = dict(); diag = []

    def support(p, dd):
        """dot weight in a corridor from point p along heading dd: half a step to LOOK ahead, LOOKW across"""
        c = np.array(tree.query_ball_point(p+0.5*look*dd, 0.5*look+lookw), int)
        if len(c) == 0: return 0.0
        off = X[c]-p; al = off @ dd; ac = np.linalg.norm(off-al[:, None]*dd, axis=1); m = (al >= 0.5*step) & (al <= look) & (ac <= 2*lookw)
        return float((wgt[c][m]*np.exp(-(ac[m]/lookw)**2)).sum())                  # dots near the corridor's centre line count most

    def meets(p_new):
        """the new node grazes a node of another line (within GRAZE cubes): returns (line, node index) of that node, else None.
        Nodes sit one step apart along a line, so a crossing line always passes within a step of one of its nodes."""
        if node_tree is None: return None
        d, k = node_tree.query(p_new)
        return None if d > graze else (int(node_owner[k]), int(node_index[k]))

    def axis_at(a, radius):
        nb = np.array(tree.query_ball_point(X[a], radius), int)
        if len(nb) < 4: return None, None
        d = X[nb]-X[nb].mean(0); ev, vec = np.linalg.eigh(d.T @ d)
        if ev[-1] < align*align*max(ev[-2], 1e-9): return None, None
        return vec[:, -1], nb

    def band_at(a, h, expect=None):
        """slice across the axis at dot a: returns the dots of the band (hill of the across-profile around a)"""
        u = np.cross(N[a], h); nu = np.linalg.norm(u)
        if nu < 1e-6: return None
        u /= nu
        nb = np.array(tree.query_ball_point(X[a], np.hypot(wmax, slice_h)+1), int)
        off = X[nb]-X[a]; al = off @ h; ac = off @ u
        keep = (np.abs(al) <= slice_h) & (np.abs(ac) <= wmax); nb, ac = nb[keep], ac[keep]
        if len(nb) == 0: return None
        strong = np.ones(len(nb), bool)
        if lead_k > 0 and len(nb) >= 4:                                          # with a clear lead, only the strong dots carry connectivity: a weak
            w = wgt[nb]; top = w.max(); lead = 1-np.median(w)/top               # dot belongs to the band only if it touches a strong dot of it
            strong = w >= min(lead_k*lead, 0.95)*top; strong[nb == a] = True
            if strong.sum() < lead_min: strong[:] = True                      # too few strong dots to carry a band on their own: keep all
        # the band is the dots connected to the anchor through neighbouring dots (within one cube diagonal); a gap of empty
        # cubes separates bands, whatever lies behind or below in projection
        from scipy.sparse import coo_matrix
        from scipy.sparse.csgraph import connected_components
        if sparse > 0 and len(nb) >= 4:
            dnn = cKDTree(X[nb]).query(X[nb], k=2)[0][:, 1]; keep2 = dnn <= sparse*np.median(dnn); keep2[nb == a] = True
            nb, ac = nb[keep2], ac[keep2]
        pr = cKDTree(X[nb]).query_pairs(1.8, output_type='ndarray')
        if not strong.all():
            pr = pr[strong[pr[:, 0]] & strong[pr[:, 1]]]                        # connectivity runs through strong dots only
        lab = connected_components(coo_matrix((np.ones(len(pr), bool), (pr[:, 0], pr[:, 1])), shape=(len(nb), len(nb))), directed=False)[1]
        ia = int(np.nonzero(nb == a)[0][0]) if (nb == a).any() else int(np.argmin(np.linalg.norm(X[nb]-X[a], axis=1)))
        conn = lab == lab[ia]
        if not strong.all() and fringe:                                          # weak dots touching the strong core join it (one hop)
            core = np.nonzero(conn & strong)[0]
            if len(core): conn |= cKDTree(X[nb][core]).query(X[nb], k=1)[0] <= 1.8
        nb, ac = nb[conn], ac[conn]
        b = np.rint(ac).astype(int); k0 = int(round(wmax))+1; prof = np.zeros(2*k0+1); np.add.at(prof, b+k0, wgt[nb])
        # start at the highest bin within one cube of the anchor, then descend both ways
        sm = np.convolve(prof, [0.25, 0.5, 0.25], 'same')                      # one-cube jitter in the profile is not a valley
        c = k0+int(np.argmax(sm[k0-1:k0+2]))-1; lo = hi = c
        def descend(i, d):
            m = sm[i]
            while 0 <= i+d < len(prof) and prof[i+d] > 0 and sm[i+d] <= valley*m: i += d; m = min(m, sm[i])
            return i
        lo = descend(c, -1); hi = descend(c, 1)
        if not hill: lo, hi = int((b+k0).min()), int((b+k0).max())            # HILL=0: the band is the whole connected set across
        if expect is not None and (hi-lo+1) > grow*expect:                   # the hill is suddenly much wider than the band so far: it has
            half = int(np.ceil(expect/2)); lo = max(lo, c-half); hi = min(hi, c+half)   # caught something beside the band; keep the expected width
            stats['clipped'] += 1
        sel = (b+k0 >= lo) & (b+k0 <= hi)
        core = sel & (sm[b+k0] >= core_frac*sm[c])                           # the strong core of the hill: where the node is read
        if ridge:
            pk = int(np.argmax(np.where((np.arange(len(sm)) >= lo) & (np.arange(len(sm)) <= hi), sm, -1))); core = sel & (np.abs(b+k0-pk) <= 1)
        if not core.any(): core = sel
        return nb[sel], (hi-lo+1), u, nb[core]

    def node_w(bd, level=None):
        """weights of the band's dots for the node. Lead = how far the strongest dots stand above the typical dot of the band
        (0: all alike, ->1: a few tower over the rest). With a clear lead the others are noise and lose their say; with no lead
        every dot counts. LEAD scales how hard the cut follows the lead."""
        w = wgt[bd]; top = w.max()
        if levp > 0 and level is not None: return np.minimum(w/level, 1.0)**levp   # LEVP: dots below the line's standard lose their say steeply
        if lead_k <= 0 or len(bd) < 4: return w
        lead = 1-np.median(w)/top
        return np.maximum(w-min(lead_k*lead, 0.95)*top, 0)      # the strongest dot always keeps a say
    def finish(line, tid, claimed):
        nonlocal node_tree, node_owner, node_index, all_nodes
        if len(line) < minlen: owner[claimed] = -1; stats['short'] += 1; return False
        pw = line_power(wgt[np.array(claimed, int)], len(line))
        if pw < powercut: owner[claimed] = -1; stats['weak'] += 1; return False
        L = np.array(line); lines.append(L); powers.append(pw)
        all_nodes = np.vstack([all_nodes, L]); node_owner = np.r_[node_owner, np.full(len(L), len(lines)-1)]; node_index = np.r_[node_index, np.arange(len(L))]
        node_tree = cKDTree(all_nodes); return True

    for s in np.argsort(-wgt):
        if owner[s] != -1: continue
        tid = len(lines); h, _ = axis_at(s, sph0)
        if h is None: stats['uncertain'] += 1; continue
        band = band_at(s, h)
        if band is None: continue
        bd, bw, _, bc = band
        if ((owner[bd] != -1) & (owner[bd] != tid)).mean() > owned_frac: continue      # that line is drawn already
        level0 = float(wgt[bc].max()); w0 = node_w(bc, level0); p0 = (X[bc]*w0[:, None]).sum(0)/w0.sum(); p0 = X[s]+(p0-X[s])-((p0-X[s]) @ h)*h; claimed = list(bd); owner[bd] = tid
        halves = []; back = max(int(2*slice_h/step)+3, 15); halves_j = []
        reasons = []; ends_at = []; halves_w = []
        for sign in (1, -1):
            nodes = [p0]; d = sign*h; width = bw; widths = []; level = level0; joins_half = []
            while True:
                # anchor: strongest dot about one step ahead along the trend, within the band's width (own dots allowed:
                # slices overlap along the line; dots of another line mean a collision)
                q = nodes[-1]; cand = np.array(tree.query_ball_point(q, reach+width), int)
                if len(cand):
                    off = X[cand]-q; al = off @ d; ac = np.linalg.norm(off-al[:, None]*d, axis=1)
                    ok = (al >= 0.5*step) & (al <= reach) & (ac <= 0.5*width+1) & (al/np.maximum(np.linalg.norm(off, axis=1), 1e-9) >= cos_cone)
                    cand = cand[ok]; acd = ac[ok]
                if len(cand) == 0: stats['no_anchor'] += 1; reasons.append('no_anchor'); ends_at.append((X[a] if 'a' in dir() else q).copy()); break
                rel = wgt[cand]/wgt[cand].max()                                     # strength relative to the peers in reach
                ald = al[ok]; score = np.exp(-(acd/afall)**2)*rel**astr*np.exp(-((ald-step)/along)**2); a = cand[np.argmax(score)]
                ha, _ = axis_at(a, max(sph0, 3*width))
                ha_from_dots = ha is not None
                if ha is None: stats['trend_used'] += 1; ha = d                    # corner or crossing: the dots give no axis; the trend carries the line
                elif ha @ d < 0: ha = -ha
                held = False
                if ha_from_dots and ha @ d < cos_bend:                             # BEND: the dots run off at a sharp angle: that is another crease, not
                    stats['bend_held'] += 1; ha = d; held = True                    # this one bending; the heading wins outright (no referee below)
                band = band_at(a, ha, expect=width)
                if band is None: reasons.append('no_band'); ends_at.append((X[a] if 'a' in dir() else q).copy()); break
                bd, bw, _, bc = band
                ow = owner[bd]
                f = None; stp = None; pn_dots = None; pred = None
                wn = node_w(bc, level); level += ((wn*wgt[bc]).sum()/wn.sum()-level)/levn; pn = (X[bc]*wn[:, None]).sum(0)/wn.sum(); pn = X[a]+(pn-X[a])-((pn-X[a]) @ ha)*ha   # node: anchor moved across the band only
                if pull > 0 or look > 0:                                                # A: the trend's prediction votes too
                    pred = nodes[-1]+((X[a]-nodes[-1]) @ d)*d; pred = X[a]+(pred-X[a])-((pred-X[a]) @ ha)*ha
                    f = pull/(1+pull)                                                   # fixed share of the trend's say
                    stp = pn-nodes[-1]; s_al = stp @ d; s_ac = np.linalg.norm(stp-s_al*d)       # the turn the dots ask for, judged over a full step:
                    if not held and s_ac > max(s_al, step)*tan_bend: stats['bend_held'] += 1; held = True   # a short step with a small side offset is not a sharp turn
                    if held:                                                            # heading held: only the dots inside the allowed cone steer the node;
                        offc = X[bc]-nodes[-1]; alc = offc @ d; acc = np.linalg.norm(offc-alc[:, None]*d, axis=1)   # the rest belong to another crease
                        inc = acc <= np.maximum(alc, step)*tan_bend
                        if inc.any(): wc = wn[inc]; pc = (X[bc][inc]*wc[:, None]).sum(0)/wc.sum(); pn = X[a]+(pc-X[a])-((pc-X[a]) @ ha)*ha
                        else: pn = pred
                    if look > 0 and not held:                                           # LOOK: probe ahead; the candidate with dot support ahead of it wins
                        dl, pl = d, pred
                        if lookn > 0 and len(nodes) > lookn:                           # LOOKN: probe the trend of the nodes before the last few, which a drift under way has not yet bent
                            dl = nodes[-4]-nodes[-1-lookn]; dl /= np.linalg.norm(dl)+1e-12; pl = nodes[-4]+((X[a]-nodes[-4]) @ dl)*dl; pl = X[a]+(pl-X[a])-((pl-X[a]) @ ha)*ha
                        st = support(pl, dl); sa = support(pn, ha)
                        if st > sa: pred = pl                                           # the long trend wins: the node goes back onto it
                        if st+sa > 0: f = pull*st/(pull*st+sa); stats['look_trend' if f > 0.8 else 'look_dots' if f < 0.2 else 'look_mixed'] += 1
                        walk.probe.append((f, st, sa, float(np.linalg.norm(pn-pred))))
                    pn_dots = pn; pn = (1-f)*pn+f*pred                              # (held: pn is already the in-cone dots' node, no referee)
                hit = meets(pn)                                                         # the node grazes a node of another line: it becomes that node, shared
                if hit is not None:
                    stats['collision'] += 1; reasons.append('collision'); ends_at.append(X[a].copy()); nodes.append(lines[hit[0]][hit[1]].copy()); joins_half.append((len(nodes)-1, hit))
                    free = bd[ow == -1]; owner[free] = tid; claimed.extend(free); break
                old = np.array(nodes[:-back]) if len(nodes) > back else None
                if old is not None:
                    dist = np.linalg.norm(old-pn, axis=1); kmin = int(np.argmin(dist))
                    if dist[kmin] < 1.5: stats['self_hit'] += 1; reasons.append('loop'); ends_at.append((X[a] if 'a' in dir() else q).copy()); nodes.append(old[kmin]); break   # the ring has come round onto itself
                free = bd[ow == -1]; owner[free] = tid; claimed.extend(free)
                nodes.append(pn); width = 0.7*width+0.3*bw; widths.append((bw, float(np.linalg.norm(pn-X[a])), len(bd), ha_from_dots, held, f, float(np.degrees(np.arccos(np.clip(ha @ d, -1, 1)))), None if pred is None else float(np.linalg.norm(pn_dots-pred)), None if pred is None else float(np.linalg.norm(pred-X[a]))))
                k = min(trend_n, len(nodes)-1); d = nodes[-1]-nodes[-1-k]; d /= np.linalg.norm(d)+1e-12
                if len(nodes) > 20000: break
            halves.append(nodes); halves_w.append(widths); halves_j.append(joins_half)
        line = halves[1][::-1]+halves[0][1:]; nb_ = len(halves[1])
        line_joins = [(nb_-1-i, hit) for i, hit in halves_j[1]]+[(nb_-1+i, hit) for i, hit in halves_j[0]]
        if len(line) >= minlen:
            for r in reasons: kept_ends[r] = kept_ends.get(r, 0)+1
            diag.append((len(lines), reasons, ends_at, halves_w[1][::-1]+[None]+halves_w[0]))
        if finish(line, tid, claimed): joins.extend((len(lines)-1, i, j, k) for i, (j, k) in line_joins)
    if log: log(f'lines {len(lines)}, nodes {sum(len(l) for l in lines)}, dots unclaimed {int((owner == -1).sum())} of {n}; trials ended: {stats}; kept lines ended: {kept_ends}')
    walk.diag = diag; walk.joins = joins; walk.power = powers
    return lines, owner


def _task(t):
    X, wgt, N, size, kw = t
    lines, owner = walk(X, wgt, N, size, log=lambda *_: None, **kw)
    return lines, int((owner == -1).sum()), walk.joins, walk.power


def walk_all(X, wgt, N, size, log=print, workers=None, **kw):
    """The walk over a whole model: the dots are split into clusters (connected within a few cubes, so no band or
    reach crosses between them), each cluster walked on its own core, biggest first. Lines never cross clusters."""
    from concurrent.futures import ProcessPoolExecutor
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    t0 = time.time(); r = max(3.0, kw.get('reach', 0.22)/size+1)
    pairs = cKDTree(X).query_pairs(r, output_type='ndarray')
    lab = connected_components(coo_matrix((np.ones(len(pairs), bool), (pairs[:, 0], pairs[:, 1])), shape=(len(X), len(X))), directed=False)[1]
    order = np.argsort(lab, kind='stable'); groups = np.split(order, np.flatnonzero(np.diff(lab[order]))+1)
    cos_bend = np.cos(np.radians(kw.get('bend', 45.0))); tan_bend = np.tan(np.radians(min(kw.get('bend', 45.0), 89.0)))   # BEND (deg): the dots' direction may steer the node only within this angle of the heading
    powercut = kw.get('powercut', 3.0)   # POWERCUT: a line whose power (line_power) is below this never existed; its dots are released
    graze = kw.get('graze', 1.0)   # GRAZE: a node within this many cubes of another line's node lands on that node (shared); the line ends there
    minlen = int(kw.get('minlen', 8)); groups = sorted([g for g in groups if len(g) >= minlen], key=len, reverse=True)
    if log: log(f'{len(X)} dots in {len(groups)} clusters (largest {len(groups[0]) if groups else 0}), {time.time()-t0:.1f} s')
    tasks = [(X[g], wgt[g], N[g], size, kw) for g in groups]; workers = workers or os.cpu_count() or 1
    if workers > 1 and len(tasks) > 1:
        with ProcessPoolExecutor(workers) as ex: results = list(ex.map(_task, tasks, chunksize=1))
    else: results = [_task(t) for t in tasks]
    lines = []; joins = []; powers = []
    for ls, _, js, pw in results:                                                 # renumber the lines of each cluster into one list
        base = len(lines); lines.extend(ls); joins.extend((base+i, k, base+j, m) for i, k, j, m in js); powers.extend(pw)
    unclaimed = sum(r_[1] for r_ in results)+(len(X)-sum(len(g) for g in groups))
    if log: log(f'lines {len(lines)}, nodes {sum(len(l) for l in lines)}, junctions {len(joins)}, dots unclaimed {unclaimed} of {len(X)}, {time.time()-t0:.1f} s on {workers} cores')
    walk_all.joins = joins; walk_all.power = powers
    return lines


def run(res, log=print, **kw):
    """App entry: res is creases.find_dots output (idx, line, turn, N, lo, size, cut). Returns node arrays in mm.
    run.joins afterwards: (line, node, other line, other node) for every shared node (a line ended on another line's node)."""
    line = res['line']; size = float(res['size']); cut = float(res.get('cut', 0))
    X = (res['idx'][line]+0.5).astype(float); wgt = np.maximum(res['turn'][line]-cut, 0)+1e-3; N = res['N'][line]
    lines = walk_all(X, wgt, N, size, log=log, **kw); run.joins = walk_all.joins; run.power = walk_all.power
    return [res['lo']+l*size for l in lines]


def main():
    G, outp = sys.argv[1], sys.argv[2]; g = np.load(G); line = g['line']; size = float(g['size']); cut = float(g['cut']) if 'cut' in g.files else 17.0
    X = (g['idx'][line]+0.5).astype(float); wgt = (np.maximum(g['turn'][line]-cut, 0)+1e-3)**env('WPOW', 1.0); N = g['N'][line]   # WPOW: how aggressively strong dots outweigh weak ones
    kw = {k.lower(): env(k, d) for k, d in [('STEP', 0.087), ('SLICE', 0.7), ('SPHERE', 0.7), ('WMAX', 1.0), ('ALIGN', 1.5), ('TREND', 6), ('CONE', 45), ('MINLEN', 8), ('OWNED', 0.5), ('VALLEY', 1.3), ('CORE', 0.0), ('GROW', 99), ('PULL', 10.0), ('AFALL', 2.0), ('ASTR', 0.0), ('HILL', 0), ('REACH', 0.22), ('ALONG', 0.0), ('LOOK', 0.35), ('LOOKW', 1.0), ('LOOKN', 0), ('LEVP', 6.0), ('LEVN', 6), ('BEND', 45.0), ('GRAZE', 1.0), ('POWERCUT', 3.0), ('LEAD', 0.0), ('LEADMIN', 8), ('FRINGE', 1), ('RIDGE', 0), ('SPARSE', 0)]}
    t0 = time.time(); lines = walk_all(X, wgt, N, size, **kw) if env('PAR', 1) > 0 else walk(X, wgt, N, size, **kw)[0]; print(f'{time.time()-t0:.1f} s')
    cv = [g['lo']+l*size for l in lines]
    src = walk_all if env('PAR', 1) > 0 else walk; joins = src.joins
    np.savez(outp, n=np.array([len(c) for c in cv]), pts=np.vstack(cv) if cv else np.zeros((0, 3)), joins=np.array(joins, int).reshape(-1, 4), power=np.array(src.power)); print('saved', outp, f'({len(joins)} junctions)')


if __name__ == '__main__': main()
