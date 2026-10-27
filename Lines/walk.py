"""The painting walk (designed 2026-10-25 with the user, the node read redesigned the same day; junctions reworked 2026-10-27).

The paint. The dots are cubes on a grid; two dots touch when their cubes lie within the reach of each other (reach 1:
the 26 next-door cubes; reach 2 crosses one empty cube, so a missing dot, grid noise, does not split the band; a hole is
crossed by reaching "just a bit", in cubes, since noise holes are cube-sized). Paint starts at the seed. Each round
paints every unpainted dot touching the dots painted last round (the front); the new dots are grouped into clumps (new
dots touching each other are one clump). While the paint spreads evenly the ring is one clump; once the band's sides
have stopped, the new dots fall in two opposite clumps: the line's two arms (or, with the seed at a band's end, one clump
lying wholly past the seed: one arm, the seed its end; agreed 2026-10-27 after freed band pieces, blocked at one end by the
line that freed them, were painted by a seed whose ring never split). The first arm is painted to its halt, then the
second. Halts: no new dots (the band ends, or the paint has come round onto its own line), a joint with another line
(below), a faulty row (the fork, below). Finding halts too when the ring never splits within the round cap.

Rows (agreed 2026-10-27, replacing the rounds: a round's new dots are a ripple from the seed, so the same band walked from
another seed was sliced at other places and slants, and two lines meeting at a joint wobbled against each other; the
user: rows are slices one dot wide across the band, not layers of the ring). The paint still finds the band's dots by
touch, round by round, and the dots an arm paints go to its pool. Rows are cut from the pool: a row is the pool's dots
within one cube along the arm's slice direction from the row's origin; the first origin is the seed node, each next one
the centre of the row just cut moved along to the row's far edge, so the origin follows the band's bends; the slice
direction is the crease direction read from the arrows of the last window of rows (membership applied), so a row is
always across the crease where it is. A row is cut only once the paint's front lies past its far edge (a new dot is
always past the front's centre, so no dot can still fall into it), or when the arm halts. When the line is found every
painted dot (rings and clumps) goes to the arm on its side of the seed node along the crease; the rows behind an arm's
first (for its first nodes' windows) are the other side's dots sliced likewise. A row with no dot (a hole) gets no node.
The node read is unchanged; only what a row is.

The node (agreed 2026-10-25 after the row centres gave a sawtooth: the dots' positions are grid points, their arrows
carry the crease). For row k the region is rows k-window .. k+window of this arm, so the paint runs `window` rows ahead
of the node. Direction: every arrow on either side of a crease is at right angles to the crease line (the two sides are
the same arrow turned about that line), so the crease direction is the direction the arrows have least component along:
the eigenvector with the smallest eigenvalue of the strength-weighted sum of n n^T over the region's arrows; its sign
follows the arm. Sureness = smallest eigenvalue over the middle one (0 a clean crease, near 1 no direction). Across:
the strength-weighted centre of the region's dots, many rows, so the grid jitter of one row averages out. Along: the
row's own centre. node = region centre moved along the direction to the row's position, then pushed onto the surface
along the region's mean arrow. Handles: reach (cubes), window (rows), the finding round cap.

Plain lines (agreed 2026-10-27, replacing the junction tree, the split test and the loose-branch rule, after line 1 arm
6 of the chimp was ended at a "junction" where its paint had merely touched a parallel crease: two clumps are not a
meeting, and the walk never tested whether the directions meet). A line is one band and nothing else. When a round's
new dots fall in two or more clumps the arm keeps the clump that continues its band, the one whose centre lies nearest
straight ahead (least offset across from the front's centre), and the other clumps are unpainted again: they are other
bands, each free for a line of its own later, found from its own seed. A missing dot that makes two clumps for one
round costs the same as a trim: the paint grows back over the freed dots by touch where they are the band's. A seed
ring that falls in three or more clumps keeps the two most opposite each other (the strongest, and the one facing it
across the seed) as its arms and frees the rest. There are no branches, no tree of arms, no junction to correct.

Joints (agreed the same day): a joint is a relation between two lines, decided from their directions, not from the
paint. When a round's new dots touch dots of another line (the contact), the arm's direction is continued ahead from
its last node and crossed with that line: the closest approach between the ray and the other line's segments. It is a
joint when the two lines come within touch of each other there (the cubes of the two closest points lie within `reach`
of each other) and the crossing lies within the dots where the two bands touch (its cube within `reach` of a contact
dot): creases meet where their dots meet, a crossing out in bare surface is none. The joint is the point on the other
line; it is recorded on both lines and the arm ends at it. No joint (parallel bands, a graze, a crossing far ahead or
behind): the other line is ignored and the walk goes on over the free dots; the test is made again at every round of
contact. A ring in the finding stage ignores other lines likewise.

The gap (agreed 2026-10-25, for bands one dot wide that hit empty cubes): when a front finds no dot next door, the
unpainted dots within `gap` mm of the front that lie ahead along the arm's direction become the next row, and the
paint goes on from them. Nothing is invented: the found dots are the row. Dots behind or beside are left alone. gap 0
turns it off.

Membership (agreed 2026-10-25, night, after the node before a split drifted off its crest: the last rows before a
split hold dots of both bands): a dot of this band has its arrow at right angles to the band's direction; over the arm's
rows so far that holds within a small scatter, the band's own. A dot's say in the direction and the node is scaled by
1 / (1 + (n.t)^2 / s^2), s^2 the band's scatter from its history (the strength-weighted mean of (n.t)^2 over every dot
of the rows noded so far, not only the dots that passed, else the scatter shrinks on itself; the arm's heading is the
reference direction, so the window being judged does not judge itself), so a dot at the band's scatter keeps half its
say and an intruding dot from a band that crosses at an angle weighs little and the node stays on its own crest. Two
bands that run parallel are not told apart. No handle.

The fork (agreed 2026-10-25, last, replacing the trend break; reduced 2026-10-27): a row is faulty when its own arrows
no longer fit the arm's line (the rms angle of its arrows off right angles to the arm's heading is over `fork` degrees;
2-8 on a stem, 15 at the first fork row, 37-53 beyond). A node reads its window with the faulty rows ahead dropped, so
a faulty influence somewhere ahead does not pull it off the band. When the node's own row is faulty the line is cut
before it: the arm ends at its last node and the rows from the faulty one on are unpainted again, free for other lines.
`fork` is the handle, in degrees. Nodes are coloured by the misfit of their whole window (blue 5, red 20).

Sub-bands, the trim (agreed 2026-10-26, after L1 a6 n13: a second, stronger crease joined the row from one side and the
nodes, the row's centre of weight, slid across onto it). A band is decided by cohesion, not proximity: a row whose
strength profile across shows two crests with a valley between them is two rows, one per crest, whether or not the
dots still touch. The trim is done the moment a row is painted (_trim, in _spread): the row is read across in bins one
cube wide (each bin its strongest dot, every dot at full strength); a valley is a run of bins lower than the runs on
both sides, its rise past the valley over `lead` of the lower neighbouring crest; our piece is the one holding the
reference, the strength-weighted centre of the row before; the dots past a valley are unpainted again and out of
everything this line does. A dip that does not repeat costs nothing, the next round grows back over them by touch.
(History: node-time cuts, lanes and a window confirmation were built and replaced the same night; the leak near line
1 arm 7 of the chimp is let go with the paint-time trim.) The paint also grabs only dots past the front's centre along
the heading (agreed 2026-10-27), so a dot behind the front stays free.

Leads (agreed 2026-10-25, last; measured only, nothing acts on it): across a crease, from the crest outward, strength
only falls or holds; a rise after a fall is a second lead, a second crease in the same band. Per node: the window's dots
are ordered across the band (along a = n x t), binned one cube wide, each bin's strength its strongest dot; from the
strongest bin outward both ways the running minimum is kept and the climb above it noted; climb = the largest climb,
as a ratio to the crest's strength. The dots beyond the valley on the climb's side are the second lead. Shown in the
rig, logged per node.

Walk(dots, snap) holds the dots, the neighbour table and the owner of every dot. start(seed) -> step_read() /
step_node() / run(). state() describes what is on screen. Every read is recorded (records) for the log and the dump.
"""
import numpy as np

KB = 1 << 16
CLUMP_COLOURS = ['orange', 'darkorchid', 'saddlebrown', 'teal', 'gold', 'maroon', 'olive', 'steelblue']   # none of these is used for anything else in the rig


def offsets(reach):
    """every grid offset within `reach` cubes of a cube (all three axes), the cube itself left out"""
    r = int(reach); return np.array([[i, j, k] for i in range(-r, r+1) for j in range(-r, r+1) for k in range(-r, r+1) if (i, j, k) != (0, 0, 0)])


def neighbours(idx, reach=1):
    """(k, m) table: the dot in each of the m cubes within `reach` of the cube, -1 when that cube holds no dot"""
    idx = idx.astype(np.int64); k = len(idx); key = (idx[:, 0]*KB+idx[:, 1])*KB+idx[:, 2]; srt = np.argsort(key, kind='stable'); skey = key[srt]
    offs = offsets(reach); nb = np.full((k, len(offs)), -1, np.int32)
    for q, o in enumerate(offs):
        kk = ((idx[:, 0]+o[0])*KB+(idx[:, 1]+o[1]))*KB+(idx[:, 2]+o[2]); pos = np.minimum(np.searchsorted(skey, kk), k-1)
        nb[:, q] = np.where(skey[pos] == kk, srt[pos], -1)
    return nb


def crease_direction(N, w):
    """the direction the arrows N (weights w) have least component along: smallest eigenvector of sum w n n^T.
    Returns (direction, sureness = smallest eigenvalue / middle eigenvalue)."""
    M = (N*w[:, None]).T @ N; ev, vec = np.linalg.eigh(M)
    return vec[:, 0], float(ev[0]/max(ev[1], 1e-18))


def ray_segment(o, t, p0, p1):
    """closest approach between the ray o + s t (s >= 0) and the segment p0..p1: returns (s, point on the segment, distance)"""
    d = p1-p0; w0 = o-p0; a = t @ t; b = t @ d; c = d @ d; dd = w0 @ t; e = w0 @ d; den = a*c-b*b
    if den < 1e-18*max(a*c, 1e-30): u = 0.0                               # parallel: the segment's start
    else: u = (a*e-b*dd)/den
    u = min(1.0, max(0.0, u)); q = p0+u*d; s = max(0.0, (q-o) @ t/max(a, 1e-18))
    return s, q, float(np.linalg.norm(o+s*t-q))


class Walk:
    def __init__(self, dots, snap=None, log=None, cap=30, reach=1, window=3, gap=0.0, fork=20.0, lead=0.0):
        self.idx = dots['idx']; self.N = dots['N'].astype(np.float64); self.w = np.asarray(dots['turn'], np.float64); self.size = float(dots['size']); self.lo = np.asarray(dots['lo'], float)
        self.X = self.lo+(self.idx+0.5)*self.size; self.k = len(self.idx); self.reach = int(reach); self.window = int(window); self.nb = neighbours(self.idx, self.reach)
        self.owner = np.full(self.k, -1, np.int32); self.lines = []; self.snap = snap or (lambda p, n: p); self.log = log or (lambda s: None); self.cap = int(cap)
        self.gap = float(gap); self.fork = float(fork); self.lead = float(lead); self._tree = None   # gap radius (mm); the dot tree for it, built when first needed
        self.cur = None   # the walk under way
        self.stop_at = None; self.paused = False   # a breakpoint: the dot the walk pauses at, the moment a read paints it (once per line); `paused` is set by that read and cleared by the next step

    @property
    def tree(self):
        if self._tree is None:
            from scipy.spatial import cKDTree; self._tree = cKDTree(self.X)
        return self._tree

    # ---- the walk under way ----
    def start(self, seed):
        """a new line from dot `seed`; refused when the dot is painted already"""
        if self.owner[seed] != -1: self.log(f'dot {seed} belongs to line {self.owner[seed]} already'); return False
        tid = len(self.lines); self.owner[seed] = tid
        self.cur = dict(id=tid, seed=int(seed), stage='find', phase='spread', round=0, front=np.array([seed], int), prev=np.zeros(0, int), new=np.zeros(0, int), clumps=[], foreign=np.zeros(0, int),
                        ring=[np.array([seed], int)], arms=[], arm=-1, queue=[], joints=[], records=[], halt=None, seed_node=None, region=np.zeros(0, int), direction=None, sure=None, rownode=None, gapped=0, member=None, share=None, cut=False, climb=None, second=None, misfit=None, dropped=[], freed=[])
        self.log(f'--- line {tid}: seed dot {seed} at {self._fmt(self.X[seed])}, strength {self.w[seed]:.1f}')
        return True

    def strongest_free(self):
        free = np.flatnonzero(self.owner == -1)
        return None if len(free) == 0 else int(free[np.argmax(self.w[free])])

    def _fmt(self, p): return '(' + ', '.join(f'{v:.4f}' for v in p) + ')'

    def step_read(self):
        """one read. Finding: spread, clump, node (decide). Advancing: spread, clump, direction, node. Returns the read done."""
        c = self.cur
        if c is None or c['halt'] is not None: return None
        ph = c['phase']
        if ph == 'spread': self._spread(); c['phase'] = 'clump'
        elif ph == 'clump': self._clump(); c['phase'] = 'node' if c['stage'] == 'find' else 'direction'
        elif ph == 'direction': self._direction(); c['phase'] = 'node'
        elif ph == 'node': self._node(); c['phase'] = 'spread'
        if self.stop_at is not None and c.get('break_hit') != self.stop_at and self.owner[self.stop_at] == c['id']:
            c['break_hit'] = self.stop_at; self.paused = True; self.log(f'*** stopped: dot {self.stop_at} was just painted by line {c["id"]} ({ph} read, round {c["round"]}, arm {c["arm"]})')
        return ph

    def step_node(self):
        """reads until a node read is done (or the walk halts)"""
        c = self.cur
        if c is None or c['halt'] is not None: return
        self.paused = False
        while True:
            ph = self.step_read()
            if ph == 'node' or ph is None or c['halt'] is not None or self.paused: return

    def run(self, max_nodes=100000):
        """to the end of the walk"""
        c = self.cur; n = 0; self.paused = False
        while c is not None and c['halt'] is None and n < max_nodes and not self.paused: self.step_node(); n += 1

    def run_arm(self, max_nodes=100000):
        """to the next halt: the arm under way halts (its end, a joint, a cut)"""
        c = self.cur; n = 0; arm = c['arm']; self.paused = False
        while c['halt'] is None and n < max_nodes:
            self.step_node(); n += 1
            if self.paused or c['arm'] != arm: return

    # ---- the paint ----
    def _heading(self, arm, front):
        """the arm's direction for the paint: its heading once it has a node, else its slice direction"""
        return arm['heading'] if arm['heading'] is not None else arm['t']

    def _gap_seek(self, front, heading):
        """unpainted dots within the gap radius of the front that lie ahead along the heading"""
        if self.gap <= 0: return np.zeros(0, int)
        cand = np.unique(np.concatenate([np.asarray(b, int) for b in self.tree.query_ball_point(self.X[front], self.gap)] or [np.zeros(0, int)]))
        cand = cand[self.owner[cand] == -1]
        if len(cand) == 0: return cand
        centre = self.X[front].mean(0); ahead = (self.X[cand]-centre) @ heading > 0; return cand[ahead]

    def _paint_from(self, front, heading=None):
        """the next row from a front: the unpainted dots next door that lie past the front in the direction of the walk (agreed
        2026-10-27: a dot behind the front's strength-weighted centre along the heading is not grabbed, it stays free); at a dead
        end the gap is sought ahead. Returns new, foreign (the dots of other lines touched), gapped"""
        c = self.cur; tid = c['id']; cand = self.nb[front]; cand = cand[cand >= 0]; cand = np.unique(cand)
        foreign = cand[(self.owner[cand] != -1) & (self.owner[cand] != tid)]; new = cand[self.owner[cand] == -1]; gapped = False
        if heading is not None and len(new):
            ahead = (self.X[new]-self._centre(front)) @ heading > 0
            if not ahead.all(): self.log(f'  {int((~ahead).sum())} touching dot(s) behind the front left free'); new = new[ahead]
        if len(new) == 0 and heading is not None:
            new = self._gap_seek(front, heading); gapped = len(new) > 0
        self.owner[new] = tid; return new, foreign, gapped

    def _spread(self):
        c = self.cur; c['clumps'] = []; c['region'] = np.zeros(0, int); c['direction'] = None; c['sure'] = None; c['member'] = None; c['climb'] = None; c['second'] = None
        heading = self._heading(c['arms'][c['arm']], c['front']) if c['stage'] == 'advance' else None
        new, foreign, gapped = self._paint_from(c['front'], heading)
        if c['stage'] == 'advance': new = self._trim(c['arms'][c['arm']], new, c['front'])
        c['new'] = new; c['foreign'] = foreign
        if gapped: c['gapped'] += 1
        c['records'].append(dict(round=c['round'], stage=c['stage'], arm=c['arm'], read='spread', front=len(c['front']), new=len(new), foreign=len(foreign), gapped=gapped))
        self.log(f'round {c["round"]} spread ({c["stage"]}{"" if c["arm"] < 0 else " arm "+str(c["arm"])}): front {len(c["front"])} dots -> {len(new)} new' + (f' (dead end: found within {self.gap:g} mm ahead)' if gapped else '') + (f', {len(foreign)} touch another line' if len(foreign) else ''))

    def _clump(self):
        """group the new dots by touching; advancing: the joint test on contact, the end, which clump is ours, then pool and slice"""
        c = self.cur; new = c['new']
        if len(new) == 0: c['clumps'] = []; self.log('clump: nothing new')
        else:
            from scipy.sparse import coo_matrix
            from scipy.sparse.csgraph import connected_components
            local = np.full(self.k, -1, np.int64); local[new] = np.arange(len(new))
            nbs = self.nb[new]; ok = nbs >= 0; ii = np.repeat(np.arange(len(new)), nbs.shape[1]).reshape(nbs.shape)[ok]; jj = local[nbs[ok]]; keep = jj >= 0
            lab = connected_components(coo_matrix((np.ones(keep.sum(), bool), (ii[keep], jj[keep])), shape=(len(new), len(new))), directed=False)[1]
            clumps = [new[lab == l] for l in range(lab.max()+1)]; clumps.sort(key=lambda g: -self.w[g].sum()); c['clumps'] = clumps
            c['records'].append(dict(round=c['round'], stage=c['stage'], arm=c['arm'], read='clump', clumps=[len(g) for g in clumps]))
            self.log(f'clump: {len(clumps)} clump(s), sizes {[len(g) for g in clumps]}, strength {[round(float(self.w[g].sum()), 1) for g in clumps]}')
        if c['stage'] == 'find': return
        arm = c['arms'][c['arm']]; cl = c['clumps']; n = len(cl)
        if len(c['foreign']):
            p = self._joint(arm, c['foreign'])
            if p is not None:
                if n: self._free(np.concatenate(cl), 'past the joint')
                c['clumps'] = []; arm['end'] = p; self._arm_halt(f'a joint with line {c["joints"][-1]["line"]}'); return
        if n == 0: self._arm_halt('the dots end' if not self._touches_own(c) else 'the paint came round onto its own line'); return
        row = cl[0] if n == 1 else self._own_clump(arm, cl)
        arm['pool'] = np.concatenate([arm['pool'], row]); c['prev'] = c['front']; c['front'] = row; arm['front'] = row; c['round'] += 1
        self._slice(arm)

    def _free(self, dots, why):
        """dots this line painted are unpainted again: other bands, free for a line of their own"""
        dots = np.asarray(dots, int)
        if len(dots) == 0: return
        self.owner[dots] = -1; self.cur['freed'].append(dots); self.log(f'  {len(dots)} dot(s) freed: {why}')

    def _own_clump(self, arm, clumps):
        """the clump that continues the arm's band: its centre nearest straight ahead, least offset across from the front's centre.
        The others are freed."""
        c = self.cur; t = self._heading(arm, c['front']); a = self._across(np.concatenate(clumps), t); o = self._centre(c['front'])
        offs = [abs(float((self._centre(g)-o) @ a)) for g in clumps]; q = int(np.argmin(offs))
        for i, g in enumerate(clumps):
            if i != q: self._free(g, f'clump {i} ({len(g)} dots) lies {offs[i]/self.size:.1f} cubes across: another band')
        c['records'].append(dict(round=c['round'], stage='advance', arm=c['arm'], read='split', kept=q, across=[round(v/self.size, 2) for v in offs], sizes=[len(g) for g in clumps]))
        self.log(f'clump: {len(clumps)} clumps; the arm keeps clump {q} ({len(clumps[q])} dots, {offs[q]/self.size:.1f} cubes across), the rest are other bands')
        return clumps[q]

    # ---- rows: slices across the crease ----
    def _slabs(self, dots, o, t):
        """`dots` sliced one cube thick along t from the point o, nearest slice first (only the slices up to the farthest dot)"""
        if len(dots) == 0: return []
        b = np.floor((self.X[dots]-o) @ t/self.size).astype(int); b = np.maximum(b, 0)
        return [dots[b == i] for i in range(b.max()+1)]

    def _slice_direction(self, arm):
        """the slice direction for the next row: the crease direction from the arrows of the last window of rows (membership
        applied, the arm's current direction the reference), its sign along the arm"""
        rows = [r for r in arm['rows'][-(self.window+1):] if len(r)]
        if not rows: return arm['t']
        g = np.concatenate(rows); m, _ = self._membership(arm, g, arm['t']); t = crease_direction(self.N[g], self.w[g]*m)[0]
        return -t if t @ arm['t'] < 0 else t

    def _slice(self, arm, final=False):
        """cut rows from the arm's pool: each row the pool's dots within one cube along the slice direction from the origin. A row
        is cut once the front's centre lies past its far edge (every new dot lies past the front's centre, so none can still
        fall into it); final: the whole pool, the arm having halted. Then the origin moves to the row's far edge at the band's
        centre and the direction is re-read."""
        c = self.cur; cut = 0
        while len(arm['pool']):
            o, t = arm['origin'], arm['t']; al = (self.X[arm['pool']]-o) @ t
            if not final and (len(arm['front']) == 0 or (self._centre(arm['front'])-o) @ t < self.size): break
            inslab = al < self.size; row = arm['pool'][inslab]; arm['pool'] = arm['pool'][~inslab]; arm['rows'].append(row); cut += 1
            cr = self._centre(row) if len(row) else o
            arm['origin'] = cr+(self.size-(cr-o) @ t)*t; arm['t'] = self._slice_direction(arm)
        if cut: self.log(f'rows: {cut} row(s) cut across the crease, {len(arm["rows"])} so far, {len(arm["pool"])} painted dots ahead of them' + (' (the arm halted: every dot sliced)' if final else ''))

    # ---- joints ----
    def _joint(self, arm, foreign):
        """the contact with another line is a joint when the arm's direction, continued from its last node, crosses that line
        within the dots where the two bands touch (see the module text). Returns the joint point, or None."""
        c = self.cur; other = int(np.bincount(self.owner[foreign]).argmax()); fo = foreign[self.owner[foreign] == other]
        L = next((l for l in self.lines if l['id'] == other), None)
        if L is None or not L.get('polylines'): self.log(f'  contact with line {other}: it has no line drawn, ignored'); return None
        o = arm['nodes'][-1] if arm['nodes'] else arm['start']; t = self._heading(arm, c['front'])
        best = None
        for pi, pl in enumerate(L['polylines']):
            for q in range(len(pl)-1):
                s, p, dist = ray_segment(o, t, pl[q], pl[q+1])
                if best is None or dist < best[2]: best = (s, p, dist, pi, q)
        s, p, dist, pi, q = best; mine = o+s*t
        nbs = self.nb[fo]; nbs = nbs[nbs >= 0]; contact = np.unique(np.concatenate([fo, nbs[self.owner[nbs] == c['id']]]))   # the touching dots of both bands
        cube = np.floor((p-self.lo)/self.size).astype(int); mcube = np.floor((mine-self.lo)/self.size).astype(int)
        apart = int(np.abs(cube-mcube).max()); near = int(np.abs(self.idx[contact]-cube).max(1).min())
        joint = apart <= self.reach and near <= self.reach
        c['records'].append(dict(round=c['round'], stage=c['stage'], arm=c['arm'], read='joint', line=other, ahead=s/self.size, apart=apart, near=near, joint=joint, pos=[float(v) for v in p]))
        self.log(f'  contact with line {other} ({len(fo)} of its dots touched): our direction meets its line {s/self.size:.1f} cubes ahead, the two lines {apart} cube(s) apart there, {near} cube(s) from the touching dots: ' + ('a joint' if joint else 'no joint, ignored'))
        if not joint: return None
        c['joints'].append(dict(pos=p, line=other, arm=arm['id']))
        L.setdefault('joints', []).append(dict(pos=p, line=c['id'])); L.setdefault('junctions', []).append((p, 2, True))
        return p

    # ---- the node ----
    def row_misfit(self, g, t):
        """how far a row's arrows are from right angles to the line t: strength-weighted rms angle, degrees"""
        off = np.degrees(np.arcsin(np.clip(np.abs(self.N[g] @ t), 0, 1))); return float(np.sqrt((self.w[g]*off**2).sum()/max(self.w[g].sum(), 1e-12)))

    def _rows_region(self, arm, k, drop_faulty=False):
        """the dots of rows k-window .. k+window of the arm (rows before the arm's first are the finding rings, nearest first).
        drop_faulty: rows ahead whose own arrows do not fit the arm's heading are left out; returns (dots, dropped row numbers)"""
        rows = arm['rows']; before = arm['before']; W = self.window; parts = []; dropped = []
        for j in range(k-W, k+W+1):
            if 0 <= j < len(rows):
                rj = rows[j]
                if len(rj) == 0: continue
                if drop_faulty and j > k and arm['heading'] is not None and self.row_misfit(rj, arm['heading']) > self.fork: dropped.append(j); continue
                parts.append(rj)
            elif j < 0 and -j <= len(before): parts.append(before[len(before)+j])   # j = -1 is the row just before the arm's first
        g = np.concatenate(parts) if parts else np.zeros(0, int)
        return (g, dropped) if drop_faulty else g

    def _direction(self):
        """the crease direction from the arrows of the region around the row to be noded (nothing yet while the rows have not run `window` ahead)"""
        c = self.cur; arm = c['arms'][c['arm']]; k = arm['k']
        if len(arm['rows'])-1 < k+self.window and arm['halt'] is None: c['region'] = np.zeros(0, int); self.log(f'direction: {len(arm["rows"])-1-k} row(s) cut ahead of row {k}, needs {self.window}; nothing yet'); return
        self._direction_at(arm, k)

    def _membership(self, arm, g, t):
        """each dot's say in this arm's band: 1 / (1 + (n.t)^2 / s^2), s^2 the band's scatter (history, else the region's own weighted median)"""
        fit = (self.N[g] @ t)**2; w = self.w[g]
        if arm['fit_w'] > 0: s2 = arm['fit_sum']/arm['fit_w']                 # the band's own scatter, from its history
        else:
            o = np.argsort(fit); cw = np.cumsum(w[o]); s2 = float(fit[o][min(np.searchsorted(cw, 0.5*cw[-1]), len(o)-1)])   # no history yet: the region's weighted median
        s2 = max(s2, 1e-12); return 1.0/(1.0+fit/s2), s2                  # a dot at the band's scatter keeps half its say; far outside, little, never nothing (a Gaussian underflowed on a bent row and lost the node)

    def leads(self, g, m, t):
        """the strength profile across the window: returns (climb ratio, mask of the second lead's dots, profile bins, strengths)"""
        w = self.w[g]*m; wg = w/max(w.sum(), 1e-12); cg = (self.X[g]*wg[:, None]).sum(0); ng = (self.N[g]*wg[:, None]).sum(0); ng /= np.linalg.norm(ng)+1e-12
        a = np.cross(ng, t); a /= np.linalg.norm(a)+1e-12; sx = (self.X[g]-cg) @ a; b = np.rint(sx/self.size).astype(int); b -= b.min()
        prof = np.zeros(b.max()+1); np.maximum.at(prof, b, self.w[g]*m)                    # the strongest dot per bin, membership applied
        top = int(np.argmax(prof)); best = 0.0; second = np.zeros(len(g), bool)
        for step in (1, -1):
            runmin = prof[top]; i = top+step; climb = 0.0; v = top; valley = top
            while 0 <= i < len(prof):
                if prof[i] < runmin: runmin = prof[i]; v = i
                elif prof[i]-runmin > climb: climb = prof[i]-runmin; valley = v   # the valley is the dip before the best rise, not the last dip
                i += step
            if climb > best and prof[top] > 0:
                best = climb; second = (b > valley) if step == 1 else (b < valley)
        return (best/prof[top] if prof[top] > 0 else 0.0), second, prof

    def _direction_at(self, arm, k):
        c = self.cur; g_all = self._rows_region(arm, k); g, dropped = self._rows_region(arm, k, drop_faulty=True); c['region'] = g; c['dropped'] = dropped
        if len(g) == 0: g = g_all; c['region'] = g
        if len(g_all) == 0:                                                   # every row of the window is a hole: nothing to read
            c['direction'] = arm['heading'] if arm['heading'] is not None else arm['t']; c['sure'] = 1.0; c['member'] = np.zeros(0); c['climb'] = 0.0; c['second'] = np.zeros(0, bool); c['misfit'] = 0.0; self.log(f'direction at row {k}: no dot in its window'); return
        ref = arm['heading'] if arm['heading'] is not None else arm['t']   # sign: the way the arm goes
        t0 = arm['heading'] if arm['heading'] is not None else crease_direction(self.N[g], self.w[g])[0]   # the reference is the arm's history (its heading); an arm with none starts from every dot
        m, s2 = self._membership(arm, g, t0); t, sure = crease_direction(self.N[g], self.w[g]*m); m, s2 = self._membership(arm, g, t)
        if t @ ref < 0: t = -t
        c['direction'] = t; c['sure'] = sure; c['member'] = m; odd = int((m < 0.1).sum())
        climb, second, prof = self.leads(g, m, t); c['climb'] = climb; c['second'] = second; hist = arm['climb_sum']/arm['climb_n'] if arm['climb_n'] else None
        c['misfit'] = self.row_misfit(g_all, crease_direction(self.N[g_all], self.w[g_all])[0])   # the whole window's arrows against the one line that fits them best (rms, deg): a crease 3-7, a fork 18-22; red when high
        c['records'].append(dict(round=c['round'], stage='advance', arm=c['arm'], read='direction', row=k, region=len(g), direction=[float(v) for v in t], sureness=sure, scatter=s2, odd=odd, climb=climb, second=int(second.sum()), profile=[round(float(v), 1) for v in prof]))
        self.log(f'direction at row {k} from {len(g)} dots ({2*self.window+1-len(dropped)} rows' + (f', faulty rows {dropped} ahead dropped' if dropped else '') + f'): {self._fmt(t)}, the whole window fits one line within {c["misfit"]:.0f} deg, sureness {sure:.3f} (0 = a clean crease, 1 = none); band scatter {np.degrees(np.arcsin(min(1.0, np.sqrt(s2)))):.1f} deg, {odd} dot(s) of another band; across profile {[int(v) for v in prof]}: climb {climb:.2f} of the crest' + (f' (arm mean {hist:.2f})' if hist is not None else '') + (f', second lead {int(second.sum())} dots' if second.any() else ''))

    def _place(self, arm, k, check=True):
        """the node of row k: the region's centre across, the row's centre along, onto the surface. check: the fork is watched
        (not while an arm is flushed at its halt, where the rows end anyway). Returns True when the arm was cut."""
        c = self.cur; g = c['region']; t = c['direction']; row = arm['rows'][k]
        if len(row) == 0: arm['k'] = k+1; self.log(f'row {k} of arm {c["arm"]} holds no dot (a hole): no node'); return False
        rm = self.row_misfit(row, arm['heading']) if arm['heading'] is not None else 0.0; c['share'] = rm
        ahead = [j for j in range(k+1, min(len(arm['rows']), k+self.window+1))]; none_valid = len(ahead) > 0 and all(j in c['dropped'] for j in ahead)   # rows ahead exist and every one was faulty
        if check and arm['heading'] is not None and (rm > self.fork or none_valid): self._cut(arm, k, rm, none_valid); return True
        arm['climb_sum'] += c['climb']; arm['climb_n'] += 1
        wg = self.w[g]*c['member']; wg = wg/max(wg.sum(), 1e-12); cg = (self.X[g]*wg[:, None]).sum(0); ng = (self.N[g]*wg[:, None]).sum(0); ng /= np.linalg.norm(ng)+1e-12
        mr, _ = self._membership(arm, row, t); wr = self.w[row]*mr; wr = wr/max(wr.sum(), 1e-12); cr = (self.X[row]*wr[:, None]).sum(0)
        arm['fit_sum'] += float((self.w[row]*(self.N[row] @ t)**2).sum()); arm['fit_w'] += float(self.w[row].sum())   # the band's scatter: every dot of the row, strength-weighted (weighing by membership shrank it on itself: 0.7 deg, a third of the band cut, a staircase)
        p = cg+((cr-cg) @ t)*t; q = self.snap(p, ng); mv = float(np.linalg.norm(q-p))
        arm['nodes'].append(q); arm['misfits'].append(c['misfit']); arm['heading'] = t; arm['k'] = k+1; c['rownode'] = (cr, cg)
        c['records'].append(dict(round=c['round'], stage='advance', arm=c['arm'], read='node', row=k, node=[float(v) for v in q], moved=mv, dots=len(row), region=len(g), sureness=c['sure']))
        self.log(f'node {len(arm["nodes"])} of arm {c["arm"]} (row {k}, {len(row)} dots): {self._fmt(q)}; across from {len(g)} dots, along from the row, snap moved {mv*1000:.1f} um; own row off the arm by {rm:.0f} deg')
        return False

    def _cut(self, arm, k, share, none_valid=False):
        """the line is cut at row k: the arm ends at its last node; the rows from k on and the pool are unpainted again"""
        why = f'cut at row {k}: ' + (f'no valid row ahead of it (every row ahead is faulty)' if none_valid and share <= self.fork else f'its arrows are {share:.0f} deg off the arm')
        arm['halt'] = why; self.log(f'arm {arm["id"]} halts: {why}')
        rest = arm['rows'][k:]+[arm['pool']]; arm['rows'] = arm['rows'][:k]; arm['k'] = len(arm['rows']); arm['pool'] = np.zeros(0, int)
        self._free(np.concatenate(rest), f'{len(rest)-1} row(s) past the cut and the dots ahead of them, another band'); self.cur['cut'] = True

    def _node(self):
        c = self.cur; cl = c['clumps']; n = len(cl)
        if c['stage'] == 'find':
            if n >= 2 or (n == 1 and self._one_sided(cl[0])):
                q, mv = self._seed_node(); c['seed_node'] = q
                if n > 2:                                                       # the two most opposite clumps are the arms; the rest are other bands
                    d = [self._centre(g)-q for g in cl]; d = [v/(np.linalg.norm(v)+1e-12) for v in d]; j = int(np.argmin([d[0] @ v for v in d]))
                    for i, g in enumerate(cl):
                        if i not in (0, j): self._free(g, f'clump {i} at the seed is neither of the two opposite arms')
                    cl = [cl[0], cl[j]]
                allp = np.concatenate(c['ring']+cl); t0 = crease_direction(self.N[allp], self.w[allp])[0]   # the crease at the seed, from every arrow painted
                if (self._centre(cl[0])-q) @ t0 < 0: t0 = -t0                   # arm 0 goes the way of the strongest clump
                al = (self.X[allp]-q) @ t0; side = [allp[al >= 0], allp[al < 0]]
                fronts = [cl[0], cl[1]] if n >= 2 else [cl[0]]
                if n >= 2 and (self._centre(cl[1])-q) @ t0 > (self._centre(cl[0])-q) @ t0: fronts = [cl[1], cl[0]]
                for i, f in enumerate(fronts):
                    t = t0 if i == 0 else -t0; a = self._new_arm(f, q, self._slabs(side[1-i], q, -t), t, side[i]); self._slice(a, final=False)
                c['stage'] = 'advance'; self._next_arm()
                self.log(f'node: {"two clumps, the line found" if n == 2 else "one clump, all of it past the seed: the line found, the seed is its end" if n == 1 else str(n)+" clumps, the two most opposite kept as arms"}; seed node {self._fmt(q)} (snap moved {mv*1000:.1f} um); the rows run {self.window} ahead of the nodes')
            elif n == 0: self._halt('no dots around the seed')
            else:
                if c['round']+1 >= self.cap: self._halt(f'the ring never split within {self.cap} rounds')
                else: c['ring'].append(c['new']); c['prev'] = c['front']; c['front'] = c['new']; c['round'] += 1; self.log('node: one ring, still spreading')
            return
        arm = c['arms'][c['arm']]
        if c['direction'] is None: self.log('node: nothing yet (the rows run ahead first)'); return
        c['cut'] = False; first = True
        while arm['halt'] is None and len(arm['rows'])-1 >= arm['k']+self.window:   # every row whose window is cut gets its node (more than one when the paint ran ahead)
            if not first: self._direction_at(arm, arm['k'])
            first = False
            if self._place(arm, arm['k']): break
        if c['cut']: self._next_arm()

    # ---- sub-bands: the trim ----
    def _across(self, row, t):
        """the across axis of a row: a = n x t, n the row's strength-weighted arrow"""
        n = (self.N[row]*self.w[row, None]).sum(0); a = np.cross(n/(np.linalg.norm(n)+1e-12), t); return a/(np.linalg.norm(a)+1e-12)

    def _centre(self, g):
        wl = self.w[g]/max(self.w[g].sum(), 1e-12); return (self.X[g]*wl[:, None]).sum(0)

    def _trim(self, arm, row, prev):
        """a freshly painted row keeps only the dots up to the valley (agreed 2026-10-26, late): read across in bins one cube
        wide (each bin its strongest dot, every dot at full strength), a valley is a run of bins lower than the runs on both
        sides, its rise past the valley over `lead` of the lower neighbouring crest. Our piece is the one holding the reference,
        the strength-weighted centre of the row before (`prev`; the row before the arm, or its start, for the first row), the
        valley bins with it. The dots past a valley cease to exist for this line: unpainted again, out of the row, the front,
        the clumps and everything after. A dip that does not repeat costs nothing: the next round grows back over them by
        touch. Returns the trimmed round (its dots go to the pool, to be sliced into rows). `j` is the round, for the records."""
        t = self._heading(arm, prev); j = self.cur['round']
        if len(row) < 2: return row
        origin = self._centre(prev) if prev is not None and len(prev) else (self._centre(arm['before'][-1]) if arm['before'] and len(arm['before'][-1]) else arm['start'])
        a = self._across(row, t); sx = (self.X[row]-origin) @ a/self.size; b = np.rint(sx).astype(int); lo = b.min(); b -= lo
        prof = np.zeros(b.max()+1); np.maximum.at(prof, b, self.w[row])
        runs = []                                                             # (first bin, last bin, strength) of each run of equal bins
        for i, v in enumerate(prof):
            if runs and runs[-1][2] == v: runs[-1][1] = i
            else: runs.append([i, i, float(v)])
        valleys = []
        for r in range(1, len(runs)-1):
            if not (runs[r][2] < runs[r-1][2] and runs[r][2] < runs[r+1][2]): continue
            left = max(x[2] for x in runs[:r]); right = max(x[2] for x in runs[r+1:]); v = runs[r][2]; rise = min(left, right)-v
            if rise > self.lead*min(left, right): valleys.append((runs[r][0], runs[r][1], rise/max(min(left, right), 1e-12)))
        if not valleys: return row
        b0 = int(np.clip(-lo, 0, len(prof)-1))                               # the bin the reference falls in
        L = [v for v in valleys if v[1] < b0]; R = [v for v in valleys if v[0] > b0]
        keep = np.ones(len(row), bool); vmask = np.zeros(len(row), bool); gone = {}
        if L: v = L[-1]; keep &= b >= v[0]; vmask |= (b >= v[0]) & (b <= v[1]); gone[-1] = (int((b < v[0]).sum()), round(v[2], 2))
        if R: v = R[0]; keep &= b <= v[1]; vmask |= (b >= v[0]) & (b <= v[1]); gone[1] = (int((b > v[1]).sum()), round(v[2], 2))
        if not gone or keep.all() or not keep.any(): return row
        far = row[~keep]; self.owner[far] = -1; arm['aside'].append(far); arm['valleys'][j] = np.concatenate([arm['valleys'][j], row[vmask]]) if j in arm['valleys'] else row[vmask]
        self.cur['records'].append(dict(round=self.cur['round'], stage='advance', arm=arm['id'], read='trim', row=j, kept=int(keep.sum()), gone={int(k): v for k, v in gone.items()}, profile=[int(v) for v in prof]))
        self.log(f'trim: round {j} of arm {arm["id"]} keeps {int(keep.sum())} of {len(row)} dots: ' + ', '.join(f'{n} past a valley on the {"right" if sd == 1 else "left"} (rise {rs:g} of the crest)' for sd, (n, rs) in gone.items()) + ' let go')
        return row[keep]

    def bands(self, arms):
        """the trims of a line's arms and the dots it freed, for the view: valley dots, and the dots let go that no line has taken since"""
        val = [v for a in arms for v in a['valleys'].values() if len(v)]
        aside = [g for a in arms for g in a['aside'] if len(g)]+[g for g in self.cur['freed'] if len(g)] if self.cur is not None and self.cur['arms'] is arms else [g for a in arms for g in a['aside'] if len(g)]
        cat = lambda L: np.unique(np.concatenate(L)) if L else np.zeros(0, int)
        aside = cat(aside); aside = aside[self.owner[aside] == -1]
        return dict(valleys=cat(val), aside=aside)

    def _one_sided(self, g):
        """does the ring's one clump lie wholly past the seed along the crease (the crease direction read from the arrows of
        every dot painted so far, every dot of the clump strictly ahead of the seed along it)? Then the paint has stopped
        surrounding the seed and grows along the crease one way: the seed is the band's end. (Measured along the crease, not
        toward the clump's centre: a seed on the band's side edge has its first ring all on one side across.)"""
        c = self.cur; allp = np.concatenate(c['ring']+[g]); t = crease_direction(self.N[allp], self.w[allp])[0]
        d = self.X[g]-self.X[c['seed']]
        if (d.mean(0) @ t) < 0: t = -t
        return bool((d @ t > 0).all())

    def _seed_node(self):
        """the seed's node: the seed moved across onto the centre of its rings, pushed onto the surface"""
        g = np.concatenate(self.cur['ring']); w = self.w[g]/self.w[g].sum(); p = (self.X[g]*w[:, None]).sum(0); n = (self.N[g]*w[:, None]).sum(0); n /= np.linalg.norm(n)+1e-12
        q = self.snap(p, n); return q, float(np.linalg.norm(q-p))

    def _touches_own(self, c):
        """does the front touch dots of this line painted before the last two rounds (a loop), rather than nothing"""
        nbs = self.nb[c['front']]; nbs = nbs[nbs >= 0]; mine = nbs[self.owner[nbs] == c['id']]
        recent = set(c['front'].tolist()) | set(c['new'].tolist()) | set(c['prev'].tolist())
        return any(int(m) not in recent for m in mine) and len(c['arms'][c['arm']]['nodes']) > 2

    def _new_arm(self, front, start, before, t, pool):
        """an arm: its front (the clump it paints on from), the seed node it starts from, the rows before it (the other side's dots
        sliced from the seed node, nearest last), its slice direction, its pool (the dots painted on its side, to be sliced)"""
        c = self.cur; a = dict(id=len(c['arms']), rows=[], before=list(before[::-1]), front=np.asarray(front, int), pool=np.asarray(pool, int), origin=np.asarray(start, float), t=np.asarray(t, float), nodes=[], halt=None, heading=None, k=0, start=np.asarray(start, float), end=None, valleys={}, aside=[], fit_sum=0.0, fit_w=0.0, climb_sum=0.0, climb_n=0, misfits=[])
        c['arms'].append(a); c['queue'].append(a['id'])
        return a

    def _next_arm(self):
        """the next queued arm becomes the one painted"""
        c = self.cur
        if not c['queue']: self._finish(); return
        a = c['arms'][c['queue'].pop(0)]; c['arm'] = a['id']; c['prev'] = a['before'][-1] if a['before'] else np.zeros(0, int)
        c['new'] = np.zeros(0, int); c['clumps'] = []; c['region'] = np.zeros(0, int); c['direction'] = None; c['rownode'] = None; c['front'] = a['front']
        self.log(f'arm {a["id"]} goes'); c['cut'] = False; self._slice(a)
        while a['halt'] is None and len(a['rows'])-1 >= a['k']+self.window:
            self._direction_at(a, a['k'])
            if self._place(a, a['k']): break
        if c['cut']: self._next_arm()

    def _arm_halt(self, why):
        c = self.cur; arm = c['arms'][c['arm']]; arm['halt'] = why; self.log(f'arm {c["arm"]} halts: {why}'); c['cut'] = False; self._slice(arm, final=True)
        while arm['k'] < len(arm['rows']):                                   # the rows cut ahead still get their nodes, with the rows there are
            self._direction_at(arm, arm['k'])
            if self._place(arm, arm['k']): break                              # a faulty row among them: the arm is cut there instead, the joint (if any) dropped
        if c['cut'] and arm['end'] is not None:
            arm['end'] = None; j = c['joints'].pop(); L = next(l for l in self.lines if l['id'] == j['line']); L['joints'].pop(); L['junctions'].pop(); self.log(f'  the joint with line {j["line"]} is dropped: the arm was cut before it')
        self._next_arm()

    def _halt(self, why):
        c = self.cur; c['halt'] = why; self.log(f'line {c["id"]} halts in the finding stage: {why}')
        self.lines.append(dict(id=c['id'], seed=c['seed'], nodes=np.zeros((0, 3)), polylines=[], misfits=[], names=[], junctions=[], joints=[], halt=why, arms=[], records=c['records']))

    def _layout(self):
        """the polylines as drawn: the two arms of the seed as one polyline, from the end of one through the seed node to the end
        of the other; an arm that ended at a joint ends on it. Each entry: (points, names, misfits), the name of a node (arm id,
        node number from 0), (-1, -1) the seed node, (arm, -2) the joint an arm ends at; names are stable for a dot set and its
        settings, so a node can be named across sessions."""
        c = self.cur; arms = c['arms']; out = []
        if not arms: return out
        def chain(a):
            pts = list(a['nodes']); nm = [(a['id'], k) for k in range(len(a['nodes']))]; mf = list(a['misfits'])
            if a['end'] is not None: pts.append(a['end']); nm.append((a['id'], -2)); mf.append(mf[-1] if mf else 0.0)
            return pts, nm, mf
        if len(arms) == 2:
            a0, a1 = arms; p0, n0, m0 = chain(a0); p1, n1, m1 = chain(a1)
            out.append((np.array(p1[::-1]+[c['seed_node']]+p0), n1[::-1]+[(-1, -1)]+n0, np.array(m1[::-1]+[m0[0] if m0 else 0.0]+m0)))
        else:
            for a in arms:
                p, nm, mf = chain(a); pts = [a['start']]+p                 # a single arm: from the seed node
                if len(pts) > 1: out.append((np.array(pts), [(-1, -1)]+nm, np.array(([mf[0]] if mf else [0.0])+mf)))
        return out

    def node_misfits(self): return [mf for _, _, mf in self._layout()]
    def node_names(self): return [nm for _, nm, _ in self._layout()]
    def polylines(self): return [pts for pts, _, _ in self._layout()]

    def _finish(self):
        c = self.cur; arms = c['arms']; pl = self.polylines(); nj = len(c['joints'])
        c['halt'] = 'done: ' + '; '.join(f'arm {a["id"]} {a["halt"]}' for a in arms)
        self.lines.append(dict(id=c['id'], seed=c['seed'], nodes=np.vstack(pl) if pl else np.zeros((0, 3)), polylines=pl, misfits=self.node_misfits(), names=self.node_names(), junctions=[(j['pos'], 2, True) for j in c['joints']], joints=list(c['joints']), halt=c['halt'], arms=[a['halt'] for a in arms], records=c['records'], bands=self.bands(arms)))
        self.log(f'line {c["id"]} done: {len(arms)} arms, {sum(len(a["nodes"]) for a in arms)} nodes, {nj} joint(s); ' + '; '.join(f'arm {a["id"]} {a["halt"]}' for a in arms))

    def reset(self):
        """drop the walk under way; its dots are freed"""
        c = self.cur
        if c is None: return
        if c['halt'] is None or not any(l['id'] == c['id'] for l in self.lines): self.owner[self.owner == c['id']] = -1
        self.cur = None

    def clear(self):
        self.owner[:] = -1; self.lines = []; self.cur = None

    def state(self):
        """what is on screen: painted dots of the line under way, the front, the new dots by clump, the region and direction of the node being read, the polylines so far, the joints"""
        c = self.cur
        if c is None: return None
        mine = np.flatnonzero(self.owner == c['id']); newest = prev = None; pl = []; mf = []; nm = []
        if c['stage'] == 'advance':
            pl = self.polylines(); mf = self.node_misfits(); nm = self.node_names(); a = c['arms'][c['arm']]
            if a['nodes']: newest = a['nodes'][-1]; prev = a['nodes'][-2] if len(a['nodes']) > 1 else a['start']
        bands = self.bands(c['arms'])
        return dict(id=c['id'], seed=c['seed'], stage=c['stage'], bands=bands, phase=c['phase'], round=c['round'], arm=c['arm'], halt=c['halt'], painted=mine, front=c['front'], new=c['new'], clumps=c['clumps'], foreign=c['foreign'],
                    region=c['region'], member=c['member'], direction=c['direction'], sure=c['sure'], share=c['share'], misfit=c.get('misfit'), climb=c['climb'], second=c['second'], climb_mean=(c['arms'][c['arm']]['climb_sum']/c['arms'][c['arm']]['climb_n'] if c['arm'] >= 0 and c['arms'][c['arm']]['climb_n'] else None), rownode=c['rownode'], polylines=pl, misfits=mf, names=nm, junctions=[(j['pos'], True) for j in c['joints']], newest=newest, prev=prev, queued=len(c['queue']), gapped=c['gapped'])
