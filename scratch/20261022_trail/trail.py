"""Trail (2026-10-22): the crease followed as the centre of a soft window of dots.
At every step the walker stands at a point p heading d. It looks at the dots around the point one step ahead through a
window that fades smoothly along the heading (ALONG) and across it (ACROSS): every dot counts by its strength times the
window's fade, nothing is cut. The window's weighted centre is the crease here; the window's weighted principal axis is
the crease's direction here. The node is placed one step ahead at the centre's across position; the new heading blends
the dots' axis with the trend of the chain (an arc fitted over the last TREACH nodes), by TSHARE. The line ends when the
window's weight falls below FADE of its running level (the dots run out), when most of the window's weight belongs to
another line (collision; the line ends on that line), or when it comes round onto itself.
All lengths in mm; the dot cloud X is in cubes (size = mm per cube) as in trace4.walk."""
import numpy as np
from scipy.spatial import cKDTree


def trail(X, wgt, N, size, log=print, **kw):
    step = kw.get('step', 0.087)/size; along = kw.get('along', 0.7)/size; across = kw.get('across', 0.35)/size
    tshare = kw.get('tshare', 0.5); treach = max(float(kw.get('treach', 6)), 1e-9); fade = kw.get('fade', 0.25); levn = kw.get('levn', 6)
    owned = kw.get('owned', 0.5); minlen = int(kw.get('minlen', 8)); align = kw.get('align', 1.5); unit = kw.get('unit', 0.087)/size
    selfhit = kw.get('selfhit', 1.5)*unit; back_min = int(kw.get('backmin', 15)); maxnodes = int(kw.get('maxnodes', 20000))
    n = len(X); tree = cKDTree(X); owner = np.full(n, -1); lines = []; stats = dict(faded=0, collision=0, self_hit=0, short=0, no_axis=0)
    R = 3*max(along, across)

    def window(q, d):
        """dots around q weighted by the window along d: returns indices, weights"""
        nb = np.array(tree.query_ball_point(q, R), int)
        if len(nb) == 0: return nb, np.zeros(0)
        off = X[nb]-q; al = off @ d; ac = np.linalg.norm(off-al[:, None]*d, axis=1)
        return nb, wgt[nb]*np.exp(-(al/along)**2)*np.exp(-(ac/across)**2)

    def axis(nb, w, c):
        """weighted principal axis of the dots about c, or None when the dots have no clear direction"""
        if w.sum() <= 0: return None
        off = (X[nb]-c)*np.sqrt(w)[:, None]; ev, vec = np.linalg.eigh(off.T @ off)
        if ev[-1] < align*align*max(ev[-2], 1e-9): return None
        return vec[:, -1]

    def trend_at(nodes, d_prev):
        P = np.array(nodes); k = len(P)-1
        if k < 1: return d_prev
        if k < 2: d = P[-1]-P[-2]; return d/(np.linalg.norm(d)+1e-12)
        seg = np.linalg.norm(np.diff(P, axis=0), axis=1); t = np.concatenate([[0.0], np.cumsum(seg)]); t = t-t[-1]
        w = np.exp(-(k-np.arange(k+1))/treach); A = np.stack([np.ones(k+1), t, t*t], 1)*w[:, None]
        b = np.linalg.lstsq(A, P*w[:, None], rcond=None)[0][1]; d = b/(np.linalg.norm(b)+1e-12)
        return d if d @ d_prev >= 0 else -d

    for s in np.argsort(-wgt):
        if owner[s] != -1: continue
        tid = len(lines)
        nb, w = window(X[s], np.zeros(3)); w = wgt[nb]*np.exp(-(np.linalg.norm(X[nb]-X[s], axis=1)/along)**2)   # the seed's window: a sphere
        c0 = (X[nb]*w[:, None]).sum(0)/w.sum(); h = axis(nb, w, c0)
        if h is None: stats['no_axis'] += 1; continue
        claimed = []; halves = []; reasons = []
        for sign in (1, -1):
            nodes = [c0]; d = sign*h; level = None; back = max(int(2*along/step)+3, back_min)
            while True:
                q = nodes[-1]+step*d; nb, w = window(q, d); tot = float(w.sum())
                if level is None: level = tot
                if tot <= 0 or tot < fade*level: stats['faded'] += 1; reasons.append('faded'); break
                level += (tot-level)/levn
                ow = owner[nb]; foreign = (ow != -1) & (ow != tid)
                if w[foreign].sum() > owned*tot:                                        # most of the weight belongs to another line: end on it
                    j = int(np.bincount(ow[foreign], weights=w[foreign]).argmax()); other = lines[j]; dn = np.linalg.norm(other-q, axis=1)
                    if dn.min() <= 2*across: nodes.append(other[int(dn.argmin())])
                    stats['collision'] += 1; reasons.append('collision'); break
                c = (X[nb]*w[:, None]).sum(0)/tot; e = axis(nb, w, c)
                if e is not None and e @ d < 0: e = -e
                dt = trend_at(nodes, d); dn_ = (1-tshare)*(e if e is not None else dt)+tshare*dt; dn_ /= np.linalg.norm(dn_)+1e-12
                pn = q+((c-q)-((c-q) @ dn_)*dn_)                                       # one step ahead, at the crease centre across
                old = np.array(nodes[:-back]) if len(nodes) > back else None
                if old is not None:
                    dist = np.linalg.norm(old-pn, axis=1); kmin = int(np.argmin(dist))
                    if dist[kmin] < selfhit: stats['self_hit'] += 1; reasons.append('loop'); nodes.append(old[kmin]); break
                mine = nb[(ow == -1) & (w >= np.exp(-1.0)*wgt[nb])]                     # claim the dots inside the window's core (within its fade lengths)
                owner[mine] = tid; claimed.extend(mine)
                nodes.append(pn); d = dn_
                if len(nodes) > maxnodes: break
            halves.append(nodes)
        line = halves[1][::-1]+halves[0][1:]
        if len(line) < minlen: owner[claimed] = -1; stats['short'] += 1; continue
        lines.append(np.array(line))
    if log: log(f'trail: lines {len(lines)}, nodes {sum(len(l) for l in lines)}, dots unclaimed {int((owner == -1).sum())} of {n}; ends: {stats}')
    return lines, owner
