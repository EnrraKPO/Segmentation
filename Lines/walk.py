"""The painting walk (designed 2026-10-25 with the user, the node read redesigned the same day).

The paint. The dots are cubes on a grid; two dots touch when their cubes lie within the reach of each other (reach 1:
the 26 next-door cubes; reach 2 crosses one empty cube, so a missing dot, grid noise, does not split the band; a hole is
crossed by reaching "just a bit", in cubes, since noise holes are cube-sized). Paint starts at the seed. Each round
paints every unpainted dot touching the dots painted last round (the front); the new dots are grouped into clumps (new
dots touching each other are one clump). While the paint spreads evenly the ring is one clump; once the band's sides
have stopped, the new dots fall in two opposite clumps: the line's two arms. From there each round's new dots is one
row of the band. The first arm is painted to its halt, then the second. Halts: a round's new dots in more than one
clump (a branch: a junction, see below), no new dots (the band ends, or the paint has come round
onto its own line), a new dot already painted by another line (a collision). Finding halts too when the ring never
splits within the round cap.

The node (agreed 2026-10-25 after the row centres gave a sawtooth: the dots' positions are grid points, their arrows
carry the crease). For row k the region is rows k-window .. k+window of this arm, so the paint runs `window` rows ahead
of the node. Direction: every arrow on either side of a crease is at right angles to the crease line (the two sides are
the same arrow turned about that line), so the crease direction is the direction the arrows have least component along:
the eigenvector with the smallest eigenvalue of the strength-weighted sum of n n^T over the region's arrows; its sign
follows the arm. Sureness = smallest eigenvalue over the middle one (0 a clean crease, near 1 no direction). Across:
the strength-weighted centre of the region's dots, many rows, so the grid jitter of one row averages out. Along: the
row's own centre. node = region centre moved along the direction to the row's position, then pushed onto the surface
along the region's mean arrow. Handles: reach (cubes), window (rows), the finding round cap.

Junctions (agreed 2026-10-25, last): a split is a junction. The arm ends at a provisional junction node (the node of
its last row), every clump starts a new arm from that node, queued and painted like any other; a split inside a branch
does the same again, so a line is a tree of arms joined at junction nodes. Once every arm at a junction has its first
node and heading (or has halted), the junction is corrected: each arm is a line (its node nearest the junction, along
its heading) and lines that meet do so at one point, so the junction moves to the point closest to all of them at once
(least squares), onto the surface; the arms start from it, so the joint stays shared. The correction is bounded: the
junction lies between the provisional node and the arms' first nodes, so a point further away than those is refused
and the provisional node stays (with near-parallel arms the closest point slides along them; seen at 150 cubes).
A seed ring that splits in two is a plain line, its seed node kept; three or more clumps at the seed make the seed a
junction.

The split test (agreed 2026-10-25, after 199 of 457 arms on the chimp were stubs of one or two nodes: a missing dot in
a row makes two clumps for one round): a split is not decided at once. Both halves are painted on, each from its own
front, for up to `window` rounds. Halves whose paint touches again were one band with a hole: their rows are merged
round by round and the arm goes on as one. A half that finds nothing was a stub and is dropped. Halves still apart
after the window are a junction, each with the rows it painted.

The gap (agreed the same day, for bands one dot wide that hit empty cubes): when a front finds no dot next door, the
unpainted dots within `gap` mm of the front that lie ahead along the arm's direction become the next row, and the
paint goes on from them. Nothing is invented: the found dots are the row. Dots behind or beside are left alone. In two
clumps they are a split like any other. gap 0 turns it off.

Membership (agreed 2026-10-25, night, after the node before a junction drifted off its crest: the last rows before a
split hold dots of both bands): a dot of this band has its arrow at right angles to the band's direction; over the arm's
rows so far that holds within a small scatter, the band's own. A dot's say in the direction and the node is scaled by
1 / (1 + (n.t)^2 / s^2), s^2 the band's scatter from its history (the strength-weighted mean of (n.t)^2 over every dot
of the rows noded so far, not only the dots that passed, else the scatter shrinks on itself; the arm's heading is the
reference direction, so the window being judged does not judge itself), so a dot
at the band's scatter keeps half its say and an intruding dot from a band that
crosses at an angle weighs little and the node stays on its own crest. Two bands that run parallel are not told apart.
No handle.

The fork (agreed 2026-10-25, last, replacing the trend break): a row is faulty when its own arrows no longer fit the
stem's line (the rms angle of its arrows off right angles to the arm's heading is over `fork` degrees; 2-8 on a stem,
15 at the first fork row, 37-53 beyond). A node reads its window with the faulty rows ahead dropped, so a faulty
influence somewhere ahead does not pull it off the band (the stray node, a cube off, read a window of fork rows). When
the node's own row is faulty the stem is cut before it: its last node is the anchor; branches already found by the
paint open from it, else the rows painted ahead become a new arm from it. `fork` is the handle, in degrees. Nodes are
coloured by the misfit of their whole window (blue 5, red 20): a red node has seen a faulty row ahead.

The loose branch (agreed 2026-10-26, after a junction at L1 a3 n5 joined two bands that merely touched at one diagonal
pair of dots: the segment drawn between them ran 81 deg off the stem and 77 deg off the branch, across both creases).
A segment must point along the trail so far or along the node it reaches; one that points along neither follows no
crease. So at a branch's first node its segment from the junction is measured against the stem's heading and against
the branch's own direction: off both by more than `kink` degrees, the branch is cut loose from the junction (removed
from its arms, drawn from its first node) and walked on as a line of its own. Measured on lines 0 and 1 of the chimp:
real forks read 1-15 deg off their own line and never kink again; the false ones turn 30-90 deg once and then run
straight at 1-8 deg a step. `kink` is the handle, 30 to start. Not built: dropping the loose branch so a later line
captures it.

Sub-bands (agreed 2026-10-26, after L1 a6 n13: a second, stronger crease joined the row from one side and the nodes,
the row's centre of weight, slid across onto it; the original crease was left behind and only re-found as arm 10 when
the dots came apart seven rows later). A band is decided by cohesion, not proximity: a row whose strength profile
across shows two crests with a valley between them is two rows, one per crest, whether or not the dots still touch.
The paint stays one (by touch, as now); each row is then read across and divided into sub-bands by its crests, and
each sub-band carries its own nodes as an arm, the two sharing one paint until their dots part. Not two lines split
from one band: one line from each of two bands, from the row where the second crest appears. The leads reading is the
tool for it, taken per row with every dot at full strength: membership dimmed the newcomer (its arrows fit another
direction) and the window stacked a drifting crest over the valley. A crease crossing ours is a junction, not noise to
filter; a graze is caught after the fact by the loose-branch rule.
Built the same night: the paint stays one (arm rows are painted whole). When row k is noded the paint has run `window`
rows ahead of it, so the decision takes them in evidence (the user's point: no blur needed, the rows ahead are there):
row k and every row ahead read across, one row at a time (bins one cube wide from the centre of the arm's last own row,
each bin its strongest dot, from the strongest bin outward a rise after a dip is a second crest, `lead` the share of
the crest the rise must reach, 0 = any dip); all showing a second crest on the same side opens the lane (a junction on
the arm's node of row k-1, carrying the stem's heading there); all showing one crest closes it. While open, row k is
divided at its valley (at the lane's last valley when row k alone shows one crest), the arm keeps the crest its last
own row sits on, the lane takes the other; the rows ahead are divided the same way, provisionally, for the window read.
The lane becomes an arm from its junction when it closes, when the paint comes apart (the clump touching the lane's
dots goes with it; one clump each way and the stem simply goes on), when a split comes under test, or when the stem
halts. Its first segment is then judged by the loose-branch rule like any branch. First run on the chimp (lead 0):
line 1 arm 6 opens the lane at row 15 (climbs 0.07-0.46 over four rows), the lane is arm 10 with 7 nodes at the
parting, cut loose (78/52 deg); four one-row lanes remain on line 1 (a lane that opened on four rows closed on the next
four: the reading shifts with its origin, the centre of the last own row, as the bins move by a fraction of a cube).
Rebuilt later the same night (agreed with the user): the reference was the flaw. A row spanning two crests has its
centre in the dip, so "which crest is ours" was a coin flip row to row and line 1's arm onto the chimp's brow drifted
from a concave crease onto the convex one beside it although every row showed the dip (22, 14, 16, 20 across). Now every
row is cut at its valleys, always, with no starting point: a valley is wherever the across profile drops and rises
again (over `lead` of the lower crest). Our piece is the one holding the reference, the strength-weighted centre of
our piece of the row before; the node and the next reference come from our piece alone; the pieces beyond a valley are
the other side's. Cuts are made once, the first time a row is read with a heading in hand, and never move. The lane
opens when row k and the window rows ahead all have an other piece on the same side, closes when none has; a piece a
lane takes is consumed, so a later lane cannot take it again. A row's other pieces on a side with no lane stay painted
but noded by nobody. First run: the brow arm opens its lane at row 5 (rises 1.0, 1.0, 0.3, 1.0), the convex crease is
an arm of its own from there; line 1 has 4 lane arms (4-7 nodes) instead of four one-row ones.
Replaced the same night, later still (the user: "any row, at any time, that spots a valley counts only the dots up to the
point the tendency shifts; anything past it ceases to exist"): the node-time cut left the front whole, so the paint kept
growing from the far side of a valley; on line 1 it leaked sideways onto a neighbouring stripe near arm 7's end, crept
along it, and when it parted the parting rule made it a branch from a node 8 cubes away. Now the trim is done the
moment a row is painted (_trim, in _spread, for the main front and for every half of a split under test): the row keeps
the dots up to the valley, the rest are unpainted again and out of everything this line does. A dip that does not
repeat costs nothing, the next round grows back over them by touch. The lanes, the node-time cut and the window
confirmation are gone with it: the other band is simply free for another line. First run: line 1 goes 36 -> 177 arms,
749 -> 6011 nodes, 4.9k -> 38.9k dots painted, with no retracing (its long arms walk new ground: the line was being
choked by its own leaked paint before). The leak near arm 7 is let go (dots 41385, 41390 free again).

The arrow turn-back (agreed 2026-10-28, after L1 a10 slid from its crease onto an equally strong one three cubes beside it,
with no dip in strength between them): across one crease the arrows' lean about the heading changes in one direction only,
from flank to flank; a lean that turns back is a second crease. In the trim, the row's dots are ordered by their distance
from the heading line through the reference (height counted, so stacked dots do not collide; one side negative) and read
edge to edge; a turn-back over `swing` degrees (a handle: two dots on one flank differ by a few degrees from their
triangles) cuts the row after the dot of the extreme lean, the shared flank; our piece is the one the reference falls in.
The seed's direction (same day, after a seed between two parallel creases placed its first nodes between them: an arm had
no heading before its first node, and rows without a heading were never trimmed): the seed node is not a node like the
others, it has no direction of its own; it is given one from the arrows of its rings, read as every node reads its own,
and each arm leaves the seed with it as its heading, toward its clump, so its rows are trimmed from the first one (the
clump itself with the seed node as the reference). Where a crossing crease lies inside the rings the read is blurred;
the sureness says so.

Leads (agreed 2026-10-25, last; measured only, nothing acts on it yet): across a crease, from the crest outward, strength
only falls or holds; a rise after a fall is a second lead, a second crease in the same band (a split before the dots
come apart, or two parallel creases walked as one band). Per node: the window's dots are ordered across the band
(along a = n x t), binned one cube wide, each bin's strength its strongest dot; from the strongest bin outward both
ways the running minimum is kept and the climb above it noted; climb = the largest climb, as a ratio to the crest's
strength (0 on a normal crease, up to 1 for a second crest as strong as the first). The dots beyond the valley on the
climb's side are the second lead. The arm keeps the running mean of the climb (its history): a climb well beyond it is
a second lead. Shown in the rig, logged per node.

Walk(dots, snap) holds the dots, the neighbour table and the owner of every dot. start(seed) -> step_read() /
step_node() / run(). state() describes what is on screen. Every read is recorded (records) for the log and the dump.
"""
import numpy as np

KB = 1 << 16
CLUMP_COLOURS = ['orange', 'darkorchid', 'saddlebrown', 'teal', 'gold', 'maroon', 'olive', 'steelblue']   # none of these is used for anything else in the rig (sky blue and cyan were indistinguishable)


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


def meeting_point(points, dirs):
    """the point closest to all the lines (point, direction) at once: sum over lines of the squared distance to the line is
    least. None when the lines are all parallel (the normal matrix is singular)."""
    A = np.zeros((3, 3)); b = np.zeros(3)
    for p, t in zip(points, dirs):
        t = t/(np.linalg.norm(t)+1e-12); Q = np.eye(3)-np.outer(t, t); A += Q; b += Q @ p
    if np.linalg.cond(A) > 1e8: return None
    return np.linalg.solve(A, b)


class Walk:
    def __init__(self, dots, snap=None, log=None, cap=30, reach=1, window=3, gap=0.0, fork=20.0, kink=30.0, lead=0.0, swing=15.0):
        self.idx = dots['idx']; self.N = dots['N'].astype(np.float64); self.w = np.asarray(dots['turn'], np.float64); self.size = float(dots['size']); self.lo = np.asarray(dots['lo'], float)
        self.X = self.lo+(self.idx+0.5)*self.size; self.k = len(self.idx); self.reach = int(reach); self.window = int(window); self.nb = neighbours(self.idx, self.reach)
        self.owner = np.full(self.k, -1, np.int32); self.lines = []; self.snap = snap or (lambda p, n: p); self.log = log or (lambda s: None); self.cap = int(cap)
        self.gap = float(gap); self.fork = float(fork); self.kink = float(kink); self.lead = float(lead); self.swing = float(swing); self._tree = None; self._sub = np.full(self.k, -1, np.int32)   # gap radius (mm); the dot tree for it, built when first needed; the half each dot was painted by during a split test
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
                        ring=[np.array([seed], int)], arms=[], arm=-1, queue=[], junctions=[], records=[], halt=None, seed_node=None, region=np.zeros(0, int), direction=None, sure=None, rownode=None, probe=None, gapped=0, member=None, share=None, broken=False, climb=None, second=None, misfit=None, dropped=[])
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
        elif ph == 'clump': self._clump(); c['phase'] = 'node' if c['stage'] == 'find' or c['probe'] is not None else 'direction'
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
        """to the next halt: the arm under way halts (its end, a junction, a collision), or a split comes under test"""
        c = self.cur; n = 0; arm = c['arm']; probing = c['probe'] is not None; self.paused = False
        while c['halt'] is None and n < max_nodes:
            self.step_node(); n += 1
            if self.paused or c['arm'] != arm or (c['probe'] is not None and not probing): return
            probing = c['probe'] is not None

    # ---- the paint ----
    def _heading(self, arm, front):
        """the arm's direction: its heading once it has a node, else from where it started toward the front"""
        if arm['heading'] is not None: return arm['heading']
        d = self.X[front].mean(0)-arm['start']; return d/(np.linalg.norm(d)+1e-12)

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
        end the gap is sought ahead. Returns new, foreign, gapped"""
        c = self.cur; tid = c['id']; cand = self.nb[front]; cand = cand[cand >= 0]; cand = np.unique(cand)
        foreign = cand[(self.owner[cand] != -1) & (self.owner[cand] != tid)]; new = cand[self.owner[cand] == -1]; gapped = False
        if heading is not None and len(new):
            ahead = (self.X[new]-self._centre(front)) @ heading > 0
            if not ahead.all(): self.log(f'  {int((~ahead).sum())} touching dot(s) behind the front left free'); new = new[ahead]
        if len(new) == 0 and len(foreign) == 0 and heading is not None:
            new = self._gap_seek(front, heading); gapped = len(new) > 0
        self.owner[new] = tid; return new, foreign, gapped

    def _spread(self):
        c = self.cur; c['clumps'] = []; c['region'] = np.zeros(0, int); c['direction'] = None; c['sure'] = None; c['member'] = None; c['climb'] = None; c['second'] = None
        if c['probe'] is not None:                                            # a split under test: every live half paints from its own front
            pr = c['probe']; news = []; foreign = []
            for i, sub in enumerate(pr['subs']):
                if sub['dead']: continue
                new, fo, gapped = self._paint_from(sub['front'], self._heading(c['arms'][c['arm']], sub['front'])); foreign.append(fo)
                new = self._trim(c['arms'][c['arm']], new, sub['rows'][-1])
                if len(new) == 0: sub['dead'] = True; sub['new'] = np.zeros(0, int); self.log(f'  half {i}: nothing ahead, a stub'); continue
                self._sub[new] = i; sub['new'] = new; sub['rows'].append(new); news.append(new)
                if gapped: c['gapped'] += 1; self.log(f'  half {i}: dead end, {len(new)} dots found within {self.gap:g} mm ahead')
            c['new'] = np.concatenate(news) if news else np.zeros(0, int); c['foreign'] = np.concatenate(foreign) if foreign else np.zeros(0, int)
            c['records'].append(dict(round=c['round'], stage=c['stage'], arm=c['arm'], read='spread', probe=pr['rounds'], halves=[len(sub['new']) for sub in pr['subs']], foreign=len(c['foreign'])))
            self.log(f'round {c["round"]} spread (arm {c["arm"]}, split under test, round {pr["rounds"]+1} of {self.window}): halves -> {[len(sub["new"]) for sub in pr["subs"]]} new')
            return
        heading = self._heading(c['arms'][c['arm']], c['front']) if c['stage'] == 'advance' else None
        new, foreign, gapped = self._paint_from(c['front'], heading)
        if c['stage'] == 'advance': new = self._trim(c['arms'][c['arm']], new, c['front'])
        c['new'] = new; c['foreign'] = foreign
        if gapped: c['gapped'] += 1
        c['records'].append(dict(round=c['round'], stage=c['stage'], arm=c['arm'], read='spread', front=len(c['front']), new=len(new), foreign=len(foreign), gapped=gapped))
        self.log(f'round {c["round"]} spread ({c["stage"]}{"" if c["arm"] < 0 else " arm "+str(c["arm"])}): front {len(c["front"])} dots -> {len(new)} new' + (f' (dead end: found within {self.gap:g} mm ahead)' if gapped else '') + (f', {len(foreign)} touch another line' if len(foreign) else ''))

    def _clump(self):
        c = self.cur; new = c['new']
        if c['probe'] is not None: self._probe_touch(); return
        if len(new) == 0: c['clumps'] = []; self.log('clump: nothing new'); return
        from scipy.sparse import coo_matrix
        from scipy.sparse.csgraph import connected_components
        local = np.full(self.k, -1, np.int64); local[new] = np.arange(len(new))
        nbs = self.nb[new]; ok = nbs >= 0; ii = np.repeat(np.arange(len(new)), nbs.shape[1]).reshape(nbs.shape)[ok]; jj = local[nbs[ok]]; keep = jj >= 0
        lab = connected_components(coo_matrix((np.ones(keep.sum(), bool), (ii[keep], jj[keep])), shape=(len(new), len(new))), directed=False)[1]
        clumps = [new[lab == l] for l in range(lab.max()+1)]; clumps.sort(key=lambda g: -self.w[g].sum()); c['clumps'] = clumps
        c['records'].append(dict(round=c['round'], stage=c['stage'], arm=c['arm'], read='clump', clumps=[len(g) for g in clumps]))
        self.log(f'clump: {len(clumps)} clump(s), sizes {[len(g) for g in clumps]}, strength {[round(float(self.w[g].sum()), 1) for g in clumps]}')

    # ---- the split test ----
    def _start_probe(self, clumps):
        """a split: every clump becomes a half with its own front, painted on for up to `window` rounds"""
        c = self.cur
        c['probe'] = dict(subs=[dict(front=g, new=g, rows=[g], dead=False) for g in clumps], rounds=0)
        for i, g in enumerate(clumps): self._sub[g] = i
        self.log(f'split into {len(clumps)} halves {[len(g) for g in clumps]}: under test for {self.window} rounds (a hole rejoins, a stub dies, a junction stays apart)')

    def _probe_touch(self):
        """halves whose new dots touch another half's paint are one band: merged, their rows united round by round"""
        c = self.cur; pr = c['probe']; subs = pr['subs']; parent = list(range(len(subs)))
        def find(i):
            while parent[i] != i: i = parent[i]
            return i
        for i, sub in enumerate(subs):
            if sub['dead'] or len(sub['new']) == 0: continue
            nbs = self.nb[sub['new']]; nbs = nbs[nbs >= 0]; nbs = nbs[self.owner[nbs] == c['id']]; others = np.unique(self._sub[nbs]); others = others[(others >= 0) & (others != i)]
            for j in others: parent[find(int(j))] = find(i)
        groups = {}
        for i in range(len(subs)): groups.setdefault(find(i), []).append(i)
        merged = [g for g in groups.values() if len(g) > 1]
        if merged:
            keep = []; marks = self._sub.copy()
            for g in groups.values():
                g = sorted(g); base = subs[g[0]]
                for j in g[1:]:
                    o = subs[j]; n = max(len(base['rows']), len(o['rows']))
                    base['rows'] = [np.concatenate([r for r in (base['rows'][q] if q < len(base['rows']) else None, o['rows'][q] if q < len(o['rows']) else None) if r is not None]) for q in range(n)]
                    base['front'] = np.concatenate([base['front'], o['front']]) if not o['dead'] else base['front']; base['new'] = np.concatenate([base['new'], o['new']]); base['dead'] = base['dead'] and o['dead']
                for old_i in g: self._sub[marks == old_i] = len(keep)        # the halves are renumbered: every dot's mark follows
                keep.append(base)
            pr['subs'] = keep; self.log(f'clump: halves touched, {len(merged)} merge(s); {sum(1 for sub in keep if not sub["dead"])} live half/halves left')
        else: self.log(f'clump: halves apart ({sum(1 for sub in subs if not sub["dead"])} live)')
        for sub in pr['subs']:
            if not sub['dead']: sub['front'] = sub['new']
        c['clumps'] = [sub['new'] for sub in pr['subs'] if len(sub['new'])]

    def _probe_decide(self):
        """after a round of the test: one live half (or none) ends the test, apart after the window is a junction"""
        c = self.cur; pr = c['probe']; arm = c['arms'][c['arm']]; pr['rounds'] += 1; live = [sub for sub in pr['subs'] if not sub['dead']]
        if len(live) <= 1 or pr['rounds'] >= self.window:
            if len(live) >= 2:                                                  # a junction: each live half an arm with its rows
                c['probe'] = None; self._sub[self._sub >= 0] = -1
                self._arm_halt(f'{len(live)} halves still apart after {pr["rounds"]} rounds: a junction', branches=[sub['rows'] for sub in live]); return
            # one band: the rows of every half, united round by round, become the arm's rows
            n = max(len(sub['rows']) for sub in pr['subs']); rows = [np.concatenate([sub['rows'][q] for sub in pr['subs'] if q < len(sub['rows'])]) for q in range(n)]
            why = 'a hole, rejoined' if len(pr['subs']) < 2 or any(len(sub['rows']) == n for sub in pr['subs'] if not sub['dead']) and len(live) == 1 and len(pr['subs']) == 1 else 'stubs dropped'
            c['probe'] = None; self._sub[self._sub >= 0] = -1; c['broken'] = False
            for r in rows:
                if arm['halt'] is not None: break
                arm['rows'].append(r); c['prev'] = c['front']; c['front'] = r; c['round'] += 1
                while arm['halt'] is None and len(arm['rows'])-1 >= arm['k']+self.window:
                    self._direction_at(arm, arm['k'])
                    if self._place(arm, arm['k']): break
            if c['broken']: self._next_arm(); return
            if arm['halt'] is not None: return
            self.log(f'node: the split was not a junction ({why}): {len(rows)} rows taken as one band, the arm goes on' if live else f'node: every half died')
            if not live: self._arm_halt('the dots end'); return
            c['front'] = np.concatenate([sub['front'] for sub in live])
        else: self.log(f'node: still under test ({len(live)} halves apart, round {pr["rounds"]} of {self.window})')

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
                if drop_faulty and j > k and arm['heading'] is not None and self.row_misfit(rj, arm['heading']) > self.fork: dropped.append(j); continue
                parts.append(rj)
            elif j < 0 and -j <= len(before): parts.append(before[len(before)+j])   # j = -1 is the row just before the arm's first
        g = np.concatenate(parts) if parts else np.zeros(0, int)
        return (g, dropped) if drop_faulty else g

    def _direction(self):
        """the crease direction from the arrows of the region around the row to be noded (nothing yet while the paint has not run `window` rows ahead)"""
        c = self.cur; arm = c['arms'][c['arm']]; k = arm['k']
        if len(arm['rows'])-1 < k+self.window and arm['halt'] is None: c['region'] = np.zeros(0, int); self.log(f'direction: the paint is {len(arm["rows"])-1-k} row(s) ahead of row {k}, needs {self.window}; nothing yet'); return
        self._direction_at(arm, k)

    def _membership(self, arm, g, t):
        """each dot's say in this arm's band: exp(-(n.t)^2 / 2 s^2), s^2 the band's scatter (history, else the region's own weighted median)"""
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
        ref = arm['heading'] if arm['heading'] is not None else (self.X[arm['rows'][0]].mean(0)-arm['start'])   # sign: the way the arm goes, away from where it started
        t0 = arm['heading'] if arm['heading'] is not None else crease_direction(self.N[g], self.w[g])[0]   # the reference is the arm's history (its heading); an arm with none starts from every dot
        m, s2 = self._membership(arm, g, t0); t, sure = crease_direction(self.N[g], self.w[g]*m); m, s2 = self._membership(arm, g, t)
        if t @ ref < 0: t = -t
        c['direction'] = t; c['sure'] = sure; c['member'] = m; odd = int((m < 0.1).sum())
        climb, second, prof = self.leads(g, m, t); c['climb'] = climb; c['second'] = second; hist = arm['climb_sum']/arm['climb_n'] if arm['climb_n'] else None
        c['misfit'] = self.row_misfit(g_all, crease_direction(self.N[g_all], self.w[g_all])[0])   # the whole window's arrows against the one line that fits them best (rms, deg): a crease 3-7, a fork 18-22; red when high
        c['records'].append(dict(round=c['round'], stage='advance', arm=c['arm'], read='direction', row=k, region=len(g), direction=[float(v) for v in t], sureness=sure, scatter=s2, odd=odd, climb=climb, second=int(second.sum()), profile=[round(float(v), 1) for v in prof]))
        self.log(f'direction at row {k} from {len(g)} dots ({2*self.window+1-len(dropped)} rows' + (f', faulty rows {dropped} ahead dropped' if dropped else '') + f'): {self._fmt(t)}, the whole window fits one line within {c["misfit"]:.0f} deg, sureness {sure:.3f} (0 = a clean crease, 1 = none); band scatter {np.degrees(np.arcsin(min(1.0, np.sqrt(s2)))):.1f} deg, {odd} dot(s) of another band; across profile {[int(v) for v in prof]}: climb {climb:.2f} of the crest' + (f' (arm mean {hist:.2f})' if hist is not None else '') + (f', second lead {int(second.sum())} dots' if second.any() else ''))

    def _place(self, arm, k, check=True):
        """the node of row k: the region's centre across, the row's centre along, onto the surface. check: the trend break is
        watched (not while an arm is flushed at its halt, where the rows end anyway)"""
        c = self.cur; g = c['region']; t = c['direction']; row = arm['rows'][k]
        rm = self.row_misfit(row, arm['heading']) if arm['heading'] is not None else 0.0; c['share'] = rm
        ahead = [j for j in range(k+1, min(len(arm['rows']), k+self.window+1))]; none_valid = len(ahead) > 0 and all(j in c['dropped'] for j in ahead)   # rows ahead exist and every one was faulty
        if check and arm['heading'] is not None and (rm > self.fork or none_valid): self._trend_break(arm, k, rm, none_valid); return True
        arm['climb_sum'] += c['climb']; arm['climb_n'] += 1
        wg = self.w[g]*c['member']; wg = wg/max(wg.sum(), 1e-12); cg = (self.X[g]*wg[:, None]).sum(0); ng = (self.N[g]*wg[:, None]).sum(0); ng /= np.linalg.norm(ng)+1e-12
        mr, _ = self._membership(arm, row, t); wr = self.w[row]*mr; wr = wr/max(wr.sum(), 1e-12); cr = (self.X[row]*wr[:, None]).sum(0)
        arm['fit_sum'] += float((self.w[row]*(self.N[row] @ t)**2).sum()); arm['fit_w'] += float(self.w[row].sum())   # the band's scatter: every dot of the row, strength-weighted (weighing by membership shrank it on itself: 0.7 deg, a third of the band cut, a staircase)
        p = cg+((cr-cg) @ t)*t; q = self.snap(p, ng); mv = float(np.linalg.norm(q-p))
        arm['nodes'].append(q); arm['misfits'].append(c['misfit']); arm['heading'] = t; arm['k'] = k+1; c['rownode'] = (cr, cg)
        c['records'].append(dict(round=c['round'], stage='advance', arm=c['arm'], read='node', row=k, node=[float(v) for v in q], moved=mv, dots=len(row), region=len(g), sureness=c['sure']))
        self.log(f'node {len(arm["nodes"])} of arm {c["arm"]} (row {k}, {len(row)} dots): {self._fmt(q)}; across from {len(g)} dots, along from the row, snap moved {mv*1000:.1f} um; own row off the stem by {rm:.0f} deg')
        if len(arm['nodes']) == 1 and arm['junction'] is not None: self._check_loose(arm); self._try_fix(arm['junction'])

    def _check_loose(self, arm):
        """the branch's first segment, from the junction to its first node, must run along the stem's line or along the
        branch's own: off both by more than `kink` degrees it follows no crease, and the branch is cut loose from the junction"""
        c = self.cur; J = c['junctions'][arm['junction']]; d = arm['nodes'][0]-J['pos']; L = float(np.linalg.norm(d))
        if L < 1e-12: return
        d /= L; offs = {}
        stem = c['arms'][J['arm_in']]['heading'] if J['arm_in'] is not None else None
        if stem is not None: offs['stem'] = float(np.degrees(np.arccos(min(1.0, abs(d @ stem)))))
        offs['own'] = float(np.degrees(np.arccos(min(1.0, abs(d @ arm['heading'])))))
        off = min(offs.values()); loose = off > self.kink
        c['records'].append(dict(round=c['round'], stage='advance', arm=c['arm'], read='loose', junction=J['id'], length=L, offs=offs, loose=loose))
        text = ' and '.join(f'{v:.0f} deg off the {k}' for k, v in offs.items())
        if not loose: self.log(f'arm {arm["id"]}: its first segment ({L/self.size:.1f} cubes) is {text}: it leaves junction {J["id"]} along a crease'); return
        J['arms'].remove(arm['id']); arm['loose'] = True
        self.log(f'arm {arm["id"]}: its first segment ({L/self.size:.1f} cubes) is {text}, over {self.kink:g}: it follows no crease; cut loose from junction {J["id"]}, walked on as a line of its own ({len(J["arms"])} arm(s) left at the junction)')

    def _node(self):
        c = self.cur; cl = c['clumps']; n = len(cl)
        if len(c['foreign']): self._halt_arm_or_line(f'collision with line {int(np.bincount(self.owner[c["foreign"]]).argmax())}'); return
        if c['stage'] == 'find':
            if n >= 2:
                q, mv = self._seed_node(); c['seed_node'] = q; before = c['ring']
                jid = None
                if n > 2: jid = self._new_junction(q, None); self.log(f'node: {n} clumps at the seed: the seed is a junction')
                # the seed's direction (agreed 2026-10-28): read from the arrows of the rings, as every node reads its own; each arm
                # leaves the seed with it as its heading, pointed toward its clump, so its rows are trimmed from the first one
                rings = np.concatenate(before); t, sure = crease_direction(self.N[rings], self.w[rings]); c['direction'] = t; c['sure'] = sure
                self.log(f'seed direction from the arrows of {len(rings)} ring dots: {self._fmt(t)}, sureness {sure:.3f} (0 = a clean crease, 1 = none)')
                for g in cl:
                    a = self._new_arm(g, q, before, jid); d = self._centre(g)-q; a['heading'] = t if d @ t >= 0 else -t
                    a['rows'][0] = self._trim(a, g, None, origin=q)
                c['records'].append(dict(round=c['round'], stage='find', read='seed', node=[float(v) for v in q], direction=[float(v) for v in t], sureness=sure, rings=len(rings)))
                c['stage'] = 'advance'; self._next_arm()
                self.log(f'node: {"two clumps, the line found" if n == 2 else str(n)+" arms"}; seed node {self._fmt(q)} (snap moved {mv*1000:.1f} um); the paint runs {self.window} rows ahead of the nodes')
            elif n == 0: self._halt('no dots around the seed')
            else:
                if c['round']+1 >= self.cap: self._halt(f'the ring never split within {self.cap} rounds')
                else: c['ring'].append(c['new']); c['prev'] = c['front']; c['front'] = c['new']; c['round'] += 1; self.log('node: one ring, still spreading')
            return
        if c['probe'] is not None: self._probe_decide(); return
        arm = c['arms'][c['arm']]
        if n == 1:
            arm['rows'].append(cl[0]); c['prev'] = c['front']; c['front'] = cl[0]; c['round'] += 1
            if c['direction'] is not None:
                c['broken'] = False
                if self._place(arm, arm['k']): self._next_arm()
            else: self.log('node: nothing yet (the paint runs ahead first)')
        elif n == 0: self._arm_halt('the dots end' if not self._touches_own(c) else 'the paint came round onto its own line')
        else: self._start_probe(cl)

    # ---- sub-bands: the trim ----
    def _across(self, row, t):
        """the across axis of a row: a = n x t, n the row's strength-weighted arrow"""
        n = (self.N[row]*self.w[row, None]).sum(0); a = np.cross(n/(np.linalg.norm(n)+1e-12), t); return a/(np.linalg.norm(a)+1e-12)

    def _centre(self, g):
        wl = self.w[g]/max(self.w[g].sum(), 1e-12); return (self.X[g]*wl[:, None]).sum(0)

    def _trim(self, arm, row, prev, origin=None):
        """a freshly painted row keeps only the dots up to the valley (agreed 2026-10-26, late): read across in bins one cube
        wide (each bin its strongest dot, every dot at full strength), a valley is a run of bins lower than the runs on both
        sides, its rise past the valley over `lead` of the lower neighbouring crest. Our piece is the one holding the reference,
        the strength-weighted centre of the row before (`prev`; the row before the arm, or its start, for the first row), the
        valley bins with it. The dots past a valley cease to exist for this line: unpainted again, out of the row, the front,
        the clumps and everything after. A dip that does not repeat costs nothing: the next round grows back over them by
        touch. Returns the trimmed row. No heading yet (an arm's first rows): the row stays whole."""
        t = arm['heading']; j = len(arm['rows'])
        if t is None or len(row) < 2: return row
        if origin is None: origin = self._centre(prev) if prev is not None and len(prev) else (self._centre(arm['before'][-1]) if arm['before'] and len(arm['before'][-1]) else arm['start'])
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
        b0 = int(np.clip(-lo, 0, len(prof)-1))                               # the bin the reference falls in
        L = [v for v in valleys if v[1] < b0]; R = [v for v in valleys if v[0] > b0]
        keep = np.ones(len(row), bool); vmask = np.zeros(len(row), bool); gone = {}
        if L: v = L[-1]; keep &= b >= v[0]; vmask |= (b >= v[0]) & (b <= v[1]); gone[-1] = (int((b < v[0]).sum()), round(v[2], 2))
        if R: v = R[0]; keep &= b <= v[1]; vmask |= (b >= v[0]) & (b <= v[1]); gone[1] = (int((b > v[1]).sum()), round(v[2], 2))
        # the arrow swing (agreed 2026-10-28): the dots ordered by their distance from the heading line through the reference,
        # one side negative; read from one edge of the row to the other, the arrows' lean about the heading must keep turning
        # the same way (one crease is one swing); a lean that turns back by more than `swing` degrees is a second crease, and the
        # row is cut after the dot of the extreme lean (the shared flank). Our piece is the one the reference falls in.
        off = self.X[row]-origin; off -= (off @ t)[:, None]*t; s = np.sign(off @ a)*np.linalg.norm(off, axis=1)
        u = np.cross(t, a); lean = np.degrees(np.arctan2(self.N[row] @ a, self.N[row] @ u)); o = np.argsort(s)
        cuts = []; d = 0; ext = lean[o[0]]; ei = 0                            # cuts: (position of the cut along s, size of the turn-back)
        for q in range(1, len(o)):
            lq = lean[o[q]]
            if d == 0:
                if abs(lq-ext) > self.swing: d = 1 if lq > ext else -1; ext = lq; ei = q
            elif (lq-ext)*d > 0: ext = lq; ei = q
            elif (ext-lq)*d > self.swing: cuts.append((0.5*(s[o[ei]]+s[o[ei+1]]), abs(ext-lq))); vmask[o[ei]] = True; d = 0; ext = lq; ei = q
        for pos, size in cuts:
            if pos < 0: m = s < pos; gone[-2] = (int(m.sum()), round(float(size), 0)); keep &= ~m
            else: m = s > pos; gone[2] = (int(m.sum()), round(float(size), 0)); keep &= ~m
        if not gone or keep.all() or not keep.any(): return row
        far = row[~keep]; self.owner[far] = -1; arm['aside'].append(far); arm['valleys'][j] = row[vmask]
        self.cur['records'].append(dict(round=self.cur['round'], stage='advance', arm=arm['id'], read='trim', row=j, kept=int(keep.sum()), gone={int(k): v for k, v in gone.items()}, profile=[int(v) for v in prof]))
        self.log(f'trim: row {j} of arm {arm["id"]} keeps {int(keep.sum())} of {len(row)} dots: ' + ', '.join((f'{n} past a valley on the {"right" if sd == 1 else "left"} (rise {rs:g} of the crest)' if abs(sd) == 1 else f'{n} past an arrow turn-back of {rs:g} deg on the {"right" if sd > 0 else "left"}') for sd, (n, rs) in gone.items()) + ' let go')
        return row[keep]

    def bands(self, arms):
        """the trims of a line's arms, for the view: valley dots, and the dots let go that no line has taken since"""
        val = [v for a in arms for v in a['valleys'].values() if len(v)]
        aside = [g for a in arms for g in a['aside'] if len(g)]
        cat = lambda L: np.unique(np.concatenate(L)) if L else np.zeros(0, int)
        aside = cat(aside); aside = aside[self.owner[aside] == -1]
        return dict(valleys=cat(val), aside=aside)

    def _seed_node(self):
        """the seed's node: the seed moved across onto the centre of its rings, pushed onto the surface"""
        g = np.concatenate(self.cur['ring']); w = self.w[g]/self.w[g].sum(); p = (self.X[g]*w[:, None]).sum(0); n = (self.N[g]*w[:, None]).sum(0); n /= np.linalg.norm(n)+1e-12
        q = self.snap(p, n); return q, float(np.linalg.norm(q-p))

    def _touches_own(self, c):
        """does the front touch dots of this line painted before the last two rounds (a loop), rather than nothing"""
        nbs = self.nb[c['front']]; nbs = nbs[nbs >= 0]; mine = nbs[self.owner[nbs] == c['id']]
        recent = set(c['front'].tolist()) | set(c['new'].tolist()) | set(c['prev'].tolist())
        return any(int(m) not in recent for m in mine) and len(c['arms'][c['arm']]['nodes']) > 2

    def _trend_break(self, arm, k, share, none_valid=False):
        """the stem is cut at row k: the arm ends at its last node (the anchor); the branches open from it"""
        c = self.cur; why = f'the stem is cut at row {k}: ' + (f'no valid row ahead of it (every row ahead is faulty)' if none_valid and share <= self.fork else f'its arrows are {share:.0f} deg off the stem')
        arm['halt'] = why; self.log(f'arm {arm["id"]} halts: {why}')
        q = arm['nodes'][-1] if arm['nodes'] else arm['start']; jid = self._new_junction(q, arm['id']); arm['end_junction'] = jid
        rest = arm['rows'][k:]; arm['rows'] = arm['rows'][:k]; arm['k'] = len(arm['rows']); pending = arm.get('pending'); arm['pending'] = None
        if pending is not None:                                                # the paint has already found the branches: they open from the anchor; the fork rows between stay painted, unnoded
            for g in pending: self._new_arm(g, q, [], jid)
            self.log(f'junction {jid} at {self._fmt(q)}: {len(pending)} arms open from it, {len(rest)} fork row(s) left unnoded')
        else:
            a = self._new_arm(rest, q, [], jid); self.log(f'anchor {jid} at {self._fmt(q)}: arm {a["id"]} re-aligns from it with the {len(rest)} rows painted ahead')
        if arm['junction'] is not None: self._try_fix(arm['junction'])
        c['probe'] = None; c['broken'] = True

    def _halt_arm_or_line(self, why):
        if self.cur['stage'] == 'find': self._halt(why)
        else: self._arm_halt(why)

    def _new_arm(self, rows, start, before, junction):
        """an arm: its first row (or rows), the point it starts from, the rows before it (nearest last), the junction it leaves (or None)"""
        rows = list(rows) if isinstance(rows, list) else [rows]
        c = self.cur; a = dict(id=len(c['arms']), rows=rows, before=list(before), nodes=[], halt=None, heading=None, k=0, start=np.asarray(start, float), junction=junction, loose=False, valleys={}, aside=[], end_junction=None, fit_sum=0.0, fit_w=0.0, share_sum=0.0, share_n=0, climb_sum=0.0, climb_n=0, misfits=[])   # before: the rows before the arm's first (the seed rings for the seed arms, nothing for a branch)
        c['arms'].append(a); c['queue'].append(a['id'])
        if junction is not None: c['junctions'][junction]['arms'].append(a['id'])
        return a

    def _new_junction(self, pos, arm_in):
        c = self.cur; j = dict(id=len(c['junctions']), pos=np.asarray(pos, float), provisional=np.asarray(pos, float), arm_in=arm_in, arms=[], fixed=False); c['junctions'].append(j); return j['id']

    def _next_arm(self):
        """the next queued arm becomes the one painted"""
        c = self.cur
        if not c['queue']: self._finish(); return
        a = c['arms'][c['queue'].pop(0)]; c['arm'] = a['id']; c['prev'] = a['before'][-1] if a['before'] else np.zeros(0, int); c['front'] = a['rows'][0]
        c['new'] = np.zeros(0, int); c['clumps'] = []; c['region'] = np.zeros(0, int); c['direction'] = None; c['rownode'] = None; c['front'] = a['rows'][-1]; c['prev'] = a['rows'][-2] if len(a['rows']) > 1 else c['prev']
        self.log(f'arm {a["id"]} goes' + (f' (from junction {a["junction"]}, {len(a["rows"])} rows painted already)' if a['junction'] is not None else ''))
        c['broken'] = False
        while a['halt'] is None and len(a['rows'])-1 >= a['k']+self.window:
            self._direction_at(a, a['k'])
            if self._place(a, a['k']): break
        if c['broken']: self._next_arm(); return
        if a.get('pending') is not None: self._arm_halt(f'{len(a["pending"])} halves apart: a junction (found before the re-alignment)', branches=a['pending'])

    def _arm_halt(self, why, branches=None):
        c = self.cur; arm = c['arms'][c['arm']]; arm['halt'] = why; self.log(f'arm {c["arm"]} halts: {why}')
        arm['pending'] = branches; c['broken'] = False
        while arm['k'] < len(arm['rows']):                                   # the rows the paint ran ahead over still get their nodes, with the rows there are
            self._direction_at(arm, arm['k'])
            if self._place(arm, arm['k']): break                              # the trend broke in these rows: the rest, and the branches, went to a new arm
        if c['broken']: self._next_arm(); return
        branches = arm['pending']; arm['pending'] = None
        if branches is not None:                                             # a junction: provisional node = the arm's last node; every clump an arm from it
            q = arm['nodes'][-1] if arm['nodes'] else arm['start']; jid = self._new_junction(q, arm['id']); arm['end_junction'] = jid
            for g in branches: self._new_arm(g, q, [], jid)                 # no rows before: the parent's rows are another crease and would drag the first nodes
            self.log(f'junction {jid} (provisional) at {self._fmt(q)}: {len(branches)} arms queued, {[sum(len(r) for r in g) if isinstance(g, list) else len(g) for g in branches]} dots')
        if arm['junction'] is not None: self._try_fix(arm['junction'])     # a halted arm counts as ready, with what it has
        self._next_arm()

    def _try_fix(self, jid):
        """correct the junction once every arm leaving it has a heading (its first node) or has halted"""
        c = self.cur; J = c['junctions'][jid]
        if J['fixed']: return
        arms = [c['arms'][i] for i in J['arms']]
        if any(a['heading'] is None and a['halt'] is None for a in arms): return
        pts, dirs = [], []
        if J['arm_in'] is not None:
            ai = c['arms'][J['arm_in']]
            if ai['heading'] is not None and ai['nodes']: pts.append(ai['nodes'][-1]); dirs.append(ai['heading'])
        for a in arms:
            if a['heading'] is not None and a['nodes']: pts.append(a['nodes'][0]); dirs.append(a['heading'])
        J['fixed'] = True
        if len(pts) < 2: self.log(f'junction {jid}: fewer than two arms with a direction, the provisional node stays'); return
        p = meeting_point(pts, dirs)
        if p is None: self.log(f'junction {jid}: the arms are parallel, the provisional node stays'); return
        bound = max(np.linalg.norm(np.asarray(x)-J['provisional']) for x in pts)   # the junction lies between the provisional node and the arms' first nodes
        if np.linalg.norm(p-J['provisional']) > bound: self.log(f'junction {jid}: the closest point lies {np.linalg.norm(p-J["provisional"])*1000:.1f} um away, beyond the arms\' first nodes ({bound*1000:.1f} um): refused, the provisional node stays'); return
        n = np.zeros(3)
        for a in arms: g = a['rows'][0]; n += self.N[g].mean(0)
        q = self.snap(p, n/(np.linalg.norm(n)+1e-12)); J['pos'] = q
        for a in arms: a['start'] = q
        if J['arm_in'] is not None and c['arms'][J['arm_in']]['nodes']: c['arms'][J['arm_in']]['nodes'][-1] = q   # the incoming arm ends at the junction: its last node is the junction node
        mv = float(np.linalg.norm(q-J['provisional']))
        c['records'].append(dict(round=c['round'], stage='advance', arm=c['arm'], read='junction', junction=jid, arms=len(pts), pos=[float(v) for v in q], moved=mv))
        self.log(f'junction {jid} corrected from {len(pts)} lines: {self._fmt(q)}, moved {mv*1000:.1f} um from the provisional node')

    def _halt(self, why):
        c = self.cur; c['halt'] = why; self.log(f'line {c["id"]} halts in the finding stage: {why}')
        self.lines.append(dict(id=c['id'], seed=c['seed'], nodes=np.zeros((0, 3)), polylines=[], misfits=[], names=[], junctions=[], halt=why, arms=[], records=c['records']))

    def _layout(self):
        """the polylines as drawn: one per arm, from the point it starts at through its nodes; the two arms of a plain seed are one
        polyline; a loose arm (cut from its junction) starts at its first node. Each entry: (points, names, misfits), the name of
        a node (arm id, node number from 0), (-1, -1) the seed node, (arm, -1) the point an arm starts from; names are stable
        for a dot set and its settings, so a node can be named across sessions."""
        c = self.cur; arms = c['arms']; out = []; seed_arms = [a for a in arms if a['junction'] is None]
        if len(seed_arms) == 2:
            a0, a1 = seed_arms
            out.append((np.array(a1['nodes'][::-1]+[c['seed_node']]+a0['nodes']), [(a1['id'], k) for k in range(len(a1['nodes'])-1, -1, -1)]+[(-1, -1)]+[(a0['id'], k) for k in range(len(a0['nodes']))],
                        np.array(a1['misfits'][::-1]+[a0['misfits'][0] if a0['misfits'] else 0.0]+a0['misfits'])))
            rest = [a for a in arms if a['junction'] is not None]
        else: rest = arms
        for a in rest:
            n = len(a['nodes'])
            if a['loose']: pts, nm, mf = list(a['nodes']), [(a['id'], k) for k in range(n)], list(a['misfits'])
            else: pts, nm, mf = [a['start']]+a['nodes'], [(a['id'], -1)]+[(a['id'], k) for k in range(n)], ([a['misfits'][0]] if a['misfits'] else [])+a['misfits']
            if len(pts) > 1: out.append((np.array(pts), nm, np.array(mf)))
        return out

    def node_misfits(self): return [mf for _, _, mf in self._layout()]
    def node_names(self): return [nm for _, nm, _ in self._layout()]
    def polylines(self): return [pts for pts, _, _ in self._layout()]

    def _finish(self):
        c = self.cur; arms = c['arms']; pl = self.polylines(); nj = len(c['junctions'])
        c['halt'] = 'done: ' + '; '.join(f'arm {a["id"]} {a["halt"]}' for a in arms)
        self.lines.append(dict(id=c['id'], seed=c['seed'], nodes=np.vstack(pl) if pl else np.zeros((0, 3)), polylines=pl, misfits=self.node_misfits(), names=self.node_names(), junctions=[(J['pos'], len(J['arms'])+(J['arm_in'] is not None), J['fixed']) for J in c['junctions']], halt=c['halt'], arms=[a['halt'] for a in arms], records=c['records'], bands=self.bands(arms)))
        self.log(f'line {c["id"]} done: {len(arms)} arms, {sum(len(a["nodes"]) for a in arms)} nodes, {nj} junction(s); ' + '; '.join(f'arm {a["id"]} {a["halt"]}' for a in arms))

    def reset(self):
        """drop the walk under way; its dots are freed"""
        c = self.cur
        if c is None: return
        if c['halt'] is None or not any(l['id'] == c['id'] for l in self.lines): self.owner[self.owner == c['id']] = -1
        self._sub[:] = -1; self.cur = None

    def clear(self):
        self.owner[:] = -1; self.lines = []; self.cur = None

    def state(self):
        """what is on screen: painted dots of the line under way, the front, the new dots by clump, the region and direction of the node being read, the polylines so far, the junctions"""
        c = self.cur
        if c is None: return None
        mine = np.flatnonzero(self.owner == c['id']); newest = prev = None; pl = []; mf = []; nm = []
        if c['stage'] == 'advance':
            pl = self.polylines(); mf = self.node_misfits(); nm = self.node_names(); a = c['arms'][c['arm']]
            if a['nodes']: newest = a['nodes'][-1]; prev = a['nodes'][-2] if len(a['nodes']) > 1 else a['start']
        bands = self.bands(c['arms'])
        return dict(id=c['id'], seed=c['seed'], stage=c['stage'], bands=bands, phase=c['phase'], round=c['round'], arm=c['arm'], halt=c['halt'], painted=mine, front=c['front'], new=c['new'], clumps=c['clumps'],
                    region=c['region'], member=c['member'], direction=c['direction'], sure=c['sure'], share=c['share'], misfit=c.get('misfit'), climb=c['climb'], second=c['second'], climb_mean=(c['arms'][c['arm']]['climb_sum']/c['arms'][c['arm']]['climb_n'] if c['arm'] >= 0 and c['arms'][c['arm']]['climb_n'] else None), rownode=c['rownode'], polylines=pl, misfits=mf, names=nm, junctions=[(J['pos'], J['fixed']) for J in c['junctions']], newest=newest, prev=prev, queued=len(c['queue']),
                    probe=None if c['probe'] is None else dict(rounds=c['probe']['rounds'], halves=[np.concatenate(sub['rows']) for sub in c['probe']['subs']], live=[not sub['dead'] for sub in c['probe']['subs']]), gapped=c['gapped'])
