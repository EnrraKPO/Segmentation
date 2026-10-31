# Lines notes

The lines stage. Started as a clean slate the night of 2026-10-24, after the Dots tool (../Dots) reached a usable state.
It consumes the dot sets the Dots tool saves in its out/dots (cubes above the cuts, with arrows, strength, spread,
direction, line, concentration, crest shares, and the settings they were made with) and the mesh the set was made from.

Nothing finds lines yet. The rule for this stage, agreed with the user: nothing is built that does not follow from a
definition, reasoned out together first (the box drawing method), so that a failure can only be a mistake in the doing.

## What the first app's walk taught (the failure list to design against, not to port)

- dies at holes in the dots
- hooks and curls at the end of a line
- floating nodes off the surface
- wobble from hard band edges
- coasting: the line keeps going past the crease on momentum
- bad landings on other lines
- duplicates of the same crease

The anchored walk, trail and arc trend of the first app were judged a tangle and are not carried over.

## 2026-10-25: the painting walk, and the rig that steps it

Designed with the user, one step at a time. The ball of the old axis read was dropped (a hard edge, blind to gaps, a
blob at crossings). In its place the paint: from the seed, every round paints the unpainted dots touching the dots
painted last round. The new dots of a round are grouped into clumps (touching each other). While the paint spreads
evenly the ring is one clump; once the band's sides have stopped, the new dots fall in two opposite clumps, the line's
two directions. From there the walk IS the painting: each round's new dots is a slice across the band, its
strength-weighted centre pushed onto the surface along the slice's facing is a node. No step, no slice, no anchor, no
reach, no trend, no handle but a round cap on the finding stage. The two arms are painted one after the other.

Halts, by design: more than one clump in a round (a branch, T, crossing: to be worked through in the next iteration),
no new dots (the band ends, or the paint has come round onto its own line), a new dot already painted by another line.
The user's rule for this stage: halt at any confusing signal and work it through next time; for now the line runs as
far as the trail is consistent.

Rig (app.py): dots coloured by strength; seed = strongest free dot or the dot clicked; Next read (spread / clump /
node), Next node, Run to the halt, Reset this walk, Clear all lines; the viewport shows the line's paint (silver),
the front (white), the new dots by clump (one colour each), the seed (red), the nodes (blue, newest large), finished
lines (navy), the picked dot (lime). Every read is logged with its numbers; Dump writes lines.npz, walk.json (every
read of every line), settings, log, screenshot.

First run, toes patch set (0.087 mm, reach 1, cut 17, lift 1), five strongest seeds: the ring split in two at round 1
every time; lines of 113-212 nodes, every arm ending at "the dots end" (none hit a branch or another line on this
set); snap moved nodes by 0.00-0.03 mm; a whole line runs in 0.1 s. Not yet looked at in pictures: where the lines
stop against the dots, and what the slices do where the band widens (clumps of 7-11 dots) or thins to 1.

## 2026-10-25 (later): the reach, and the node read from the arrows

The row-centre nodes gave a sawtooth on the chimp (1 um cubes): the band is 3-5 ragged dots wide, a row is 4-9 scattered
dots, and the centre of a row jumps about a cube from row to row. The snap was checked and is clean (under a cube).
Flattened onto the band's plane the rows do follow the band; the wobble is the dots' grid scatter, not the crease.

Reach (Settings, cubes, default 1): two dots touch when their cubes lie within this many cubes; 2 crosses one empty
cube, so grid-noise holes do not split the band. Agreed as "reach just a bit, in cubes, to capture noise". At reach 2
each row is two cubes deep, half the nodes, half the wobble (turn per node 20 -> 9 deg median).

Node read, agreed after the theory session (arrows, dot product, the door): the dots' positions are grid points; their
arrows carry the crease. Every arrow on either side of a crease is at right angles to the crease line, so the crease
direction is the direction the region's arrows have least component along (smallest eigenvector of the strength-
weighted sum of n n^T), with a sureness (smallest over middle eigenvalue). The node is placed across at the
strength-weighted centre of the rows around it (window rows before and after, Settings "Node window (rows)", the paint
runs that far ahead), along at its own row's centre, then pushed onto the surface. Direction and across split the job;
no blend number. Measured on the chimp seed, reach 1: turn per node median 20 deg (window 0, the row alone) -> 3 deg
(window 3) -> 2 deg (window 6); sureness median 0.01. Both arms still halt at a branch, as designed.
Rig: the window's dots light green, the crease direction a dark green arrow with its sureness in the legend, the row's
centre yellow, the window's centre light green. A direction read sits between clump and node.

## 2026-10-25 (last): junctions

A split is a junction (agreed): the arm ends at a provisional junction node (the node of its last row), every clump
starts a new arm from that node, queued and painted like any other, so a line is a tree of arms. Once every arm at a
junction has its first node and heading (or has halted), the junction is corrected: each arm is a line (its node
nearest the junction along its heading), lines that meet do so at one point, the junction moves to the point closest
to all of them (least squares, walk.meeting_point), onto the surface; the arms start from it, so the joint stays
shared. Drawn magenta while provisional, purple once corrected. A branch's node window holds its own rows only: with
the parent's last rows in it, the branch's first nodes were dragged toward the parent band (junction 2 moved 14 um
under that pull; 1.2 um without). Chimp seed 53873, reach 1, window 3: 16 arms, 562 nodes, 7 junctions, corrections
0.8-2.2 um, 0.3 s. Seen in the flattened picture, open: a cluster of three junctions within a few um at the top left
(arms of 1-2 dots: a ragged band splitting on noise), and a split of a wide band into two parallel arms (red, lower
left) where the band is two crests or one wide one.

Speed (2026-10-25, last): the snap cast a ray over all 1.5M triangles per node (70 ms); now it looks only at the
triangles of the 16 nearest vertices (app.make_snap, Moller-Trumbore both ways, else the closest point on them):
0.09 ms per node, identical result on 200 chimp nodes. The redraw recomputed the mesh's shading normals each time
(0.4 s on the chimp); the mesh is built once with its point normals and reused: a redraw with a walk is 0.03 s.

## 2026-10-25 (night): the split test, the gap, the bounded correction

The user's dump (chimp, seed 23829): lines through the air, landing on other surfaces. Measured: 457 arms of which
199 stubs of one or two nodes; 50 of 222 junction corrections moved over 20 um, the largest 151 um (150 cubes); 37 %
of splits had their two arms within 15 deg of each other. Two faults: a hole in a ragged row makes two clumps for one
round and was taken for a junction; and "the point closest to all the lines" slides along near-parallel lines.
Built, as agreed:
- The split test: a split is painted on from each half for up to `window` rounds; halves that touch again were a hole
  (rows merged round by round, the arm goes on as one), a half that finds nothing was a stub (dropped), halves still
  apart after the window are a junction with the rows they painted.
- The gap (Settings "Gap radius (mm)", 0 = off): at a dead end, the unpainted dots within the radius that lie ahead of
  the arm become the next row. Two clumps there are a split like any other.
- The correction is bounded: the junction lies between the provisional node and the arms' first nodes; a closest
  point further away is refused and the provisional node stays.
Both seeds, reach 1, window 3, gap 0.003 mm: stubs 199 -> 0; line 1 arms 457 -> 227, junctions 222 -> 110, of which
78 corrected (moved 2 um median, 6 max) and 32 refused; 75 splits proved holes or stubs; 16 gaps crossed; no drawn
segment over 10 um. 3 s for both trees. Open: the remaining junctions on wide bands that split into two parallel arms
(two crests, or one wide one) stay junctions, since their halves never touch.
Drawing (same night): the user saw the sawtooth back. The nodes were unchanged (turn per node 3 deg at window 3 without
the snap, 7 with it); the drawn copy, lifted to the dots' height along the nearest vertex's normal, wobbled 33 deg
because that normal flips sides at a crease from node to node. The lift now follows the mean normal of the 64 nearest
vertices (7.7 deg, the nodes' own). The nodes never moved; only the picture did.

## 2026-10-25 (night): membership, the dots of the other band at a junction

The user: the node before a junction drifts off the crest. Cause: the last rows before a split hold dots of both
bands, and their centre is the centre of a Y-shaped blob. Agreed: a dot of this band has its arrow at right angles to
the band's direction, within the band's own scatter, known from the arm's history; dots outside it are the other
band. Built: a dot's say in the direction and the node is exp(-(n.t)^2 / 2 s^2), s^2 the band's scatter (the
strength-weighted mean of (n.t)^2 over the rows noded so far), floored at the region's own weighted median. The floor
was found the hard way: a clean crease's history scatter is near zero, a slightly bent row then scored zero
everywhere, the normalised weights blew up and nodes landed 75-500 um away; with the floor, at least half the region
always keeps a full say. Rig: the window's dots that fit are light green, those of another band salmon; the log
prints the band's scatter in degrees and the count of other-band dots. Both seeds: corrections 2.5 um median, 7.9 max.
Parallel bands are not told apart by this; the Y point itself is still read where the rows separate, past the meeting
of the creases (joints.png / joints2.png in the session's scratch).
Membership, second build (same night): the user's fork (seed 23829, arm 1, row 22) still drifted: the median floor
let the window judge itself, and at the fork the other band is the majority, so the scatter rose 4 -> 19 deg and
every intruder counted as own; the direction flipped and the last nodes went down the slope. Now the reference is
the arm's own heading (its history) and the fall-off is 1/(1 + fit/s^2): half a say at the band's scatter, little far
out, never nothing, so no floor is needed. At that fork the scatter stays 2.8-3.1 deg through the rows and 5-12 dots
per window are called the other band; the direction turns smoothly. Settings "Line lift (cubes)", default 0: lines
drawn on the nodes (the earlier one-cube lift read as floating).
Membership, third build (same night): the user's dump showed a staircase (turn per node 9 deg median, 29 at the 90th;
a third of each window called "another band", scatter 0.7 deg). The band's scatter was learnt from the dots weighted
by their own membership, so it shrank on itself until only the best-fitting dots counted and the node jumped between
them. Now the scatter is the strength-weighted mean over every dot of each row. Seed 53873: turn per node back to 7.3
deg median, scatter 4 deg, 1 dot of 45 called the other band; the fork of seed 23829 turns smoothly (fork_b2.png).

## 2026-10-25 (last): the trend break, built and not yet sound

The stray node at the ear-hat-head fork (stray_node.png in Lines/in): the user's reading, agreed: the stem crease ends
where a new plane (the head) emerges; the arrows there no longer share the stem's relationship, which is legitimate,
not noise. Agreed design: the reference (direction, scatter) refined every node with its history kept, so a bend is
followed; the line ends where the trend's continuity breaks; a new arm re-aligns from that anchor.
Built: at each node the share of the window's strength clearly outside the band (membership under 0.1, an arrow more
than three scatters off the trend) is measured; above `brk` (Settings "Trend break (share outside)") the arm ends at
its last node, the rows painted ahead become a new arm from that anchor (a split already found follows it), and the
junction correction treats the anchor as a corner. The check runs in the flush too.
Measured on the stem of seed 23829 (arm 1): outside share 0.00-0.20 on most rows, 0.29-0.37 at rows 14-15 (a bend or
something else, not looked at), 0.37-0.43 at the fork rows 22-25. At brk 0.25 the fork is caught but so is row 14,
and the big seed breaks 137 times. The first statistic tried (mean membership) was worse: stem 0.42-0.72, fork 0.38.
Default set to 1.0 (off) until the measure separates a bend from a new plane. Open with the user.

## 2026-10-25 (night, last): the fork read from the arrows, the faulty row

Settled with the user on the ear-hat-head fork (seed 23829, arm 1). Measured: the arrows of a window fit one line
within 3-7 deg along the stem and 18-22 deg at the fork (no single crease direction exists there: three creases'
dots in one window). Per row, the row's own arrows against the stem's line: 2-8 deg on the stem (one bump of 13),
15 at row 23, 37 and 53 beyond. The strength across the rows shows one crest throughout (the user checked the
profiles): the arrows see the fork, the strength does not. The stray node (row 23, a cube off) was read from a
window of fork rows with the arrow test throwing out one side.
Rule (the user's): a red node has seen a faulty row somewhere ahead; it drops the faulty rows from its window and
reads the rest. A node whose own row is faulty cuts the stem: the node before it is the junction, and the branches
the paint has found open from it (the fork rows between stay painted, unnoded); with no branches known yet the rows
ahead become a new arm re-aligning from the anchor. Faulty = the row's arrows more than `fork` degrees off right
angles to the stem's line, Settings "Faulty row (deg)", default 20 (replaces the trend-break share). Nodes are
coloured by their whole window's fit (blue 5 deg, red 20): red = a faulty row ahead.
Result at the fork: the stray's deviation 1.08 -> 0.08 um, the stem cut at row 24, the junction at row 23 on the
band, the two branches from it. Side effect to look at: seed 53873 now cuts 16 times (0 before) where rows exceed
20 deg, its junction count 5 -> 16; the big seed 110 -> 145 junctions, 77 cuts. Turn per node unchanged.
Speed (last thing): with both chimp trees on screen (279 arms) every arm was drawn as its own actors, twice over with
the faint copy: over a thousand actors remade at every step, the app sluggish. All arms are now one line actor and
one node actor (coloured by misfit), plus their faint copies: 17 actors, a redraw 0.06 s, ten node steps 0.5 s. A
missing state key (misfit) that crashed the redraw during a direction read was fixed on the way.

## 2026-10-26: rig only, no walk change

Run all lines: seeds the strongest free dot, runs the line to its end, repeats until no free dot is left. The reads go to
the session log file only (a six-line run wrote a 10 MB log); the box gets one line per finished line and a summary.
Names: a dot is "dot N", its number in the loaded set, as the log already calls it (the cube index is shown too, the same
across sets of one cube size). A node is "L2 a5 n7": line 2, arm 5, node 7, numbered as the log does; "L2 seed" and
"L2 a5 start" are the seed node and the point an arm starts from. Nodes are pickable like dots (click), the picked one
is named in the panel, the legend and the log, and the Find box selects a name typed in and moves the camera to it.
Show painted dots: hides the khaki (earlier lines) and tan (this line) paint; the current read's markers stay.
Toes set: 91 lines, 1302 nodes, all 5537 dots painted, 0.8 s.

## 2026-10-26: the loose branch (the first walk change since the leads)

Found with the user at L1 a3 n5 → L1 a7 n1 on the chimp set (reach 1, window 3, gap 0, fork 20): the junction joined
two bands that touched at one diagonal pair of dots (37163 of arm 3, 37720 the one-dot clump that split the ring). The
segment drawn between them ran 81 deg off the stem and 77 deg off the branch: across both creases, a line where no
crease is. Rule agreed: a segment must point along the trail so far or along the node it reaches; off both it follows
no crease. Measured on every junction of lines 0 and 1 (56 segments): real forks 1-15 deg off their own line and
straight after; the false ones turn 30-90 deg once, then 1-8 deg a step. So at a branch's first node the segment from
the (provisional) junction is checked against the stem's heading and the branch's own direction; over `kink` degrees
off both (handle "Loose branch", 30) the branch is cut loose: removed from the junction, drawn from its first node,
walked on as a loose arm of the same line (its names stay L1 a7 n…). Line 0: arms 5, 6, 11, 12, 14 cut loose; line 1:
arms 7, 11, 15, 17, 21, 23, 27. Not built: dropping the loose branch so a later line captures it (the other option
the user named). Caveat: the check runs before the junction correction; two borderline arms (L0 a12 at 33 deg, L1 a21
at 32 deg off their own line) read 13 and 11 deg against the corrected junction and would have been kept.

Rig additions of the same day not noted above: "Copy the picked name" button, "Show lines and nodes" toggle.

## 2026-10-26, later: sub-bands (cohesion, not proximity)

Found with the user at L1 a6 n13: a second, stronger crease joined the row from one side; the node, the row's centre
of weight, slid onto it, the original crease was left behind and re-found as arm 10 seven rows later. Agreed: a band
is decided by cohesion, not proximity. A row with two crests and a valley between is two rows, one per crest; not two
lines split from one band but one line from each of two bands. Built in walk.py (see the docstring, "Sub-bands"): rows
stay whole in the paint; at node time row k and the `window` rows ahead are each read across; all two-crested on one
side opens a lane (a junction on the stem's node), all one-crested closes it; while open the rows are divided at the
valley and the lane's dots become an arm from the junction when it closes, the paint parts, or the stem halts. The
loose-branch rule then judges its first segment. Handle "Second crest (share of the crest)", `lead`, 0 = any dip (the
user's choice: try it at full weight first). Bugs found on the way: the leads reading's valley was the last dip, not
the dip before the best rise (fixed in both readings); membership had dimmed the newcomer in the leads reading (not
used in the per-row read).
Results, chimp set, lines 0 and 1, lead 0: line 0 unchanged (28 arms); line 1 32 -> 44 arms, 9 lane arms, of them
four of one row (opened on four agreeing rows, closed on the next four: the per-row read moves with its origin, the
last own row's centre, by a fraction of a cube and the bins shift). At arm 6 the lane opens at row 15 and is arm 10
(7 nodes), cut loose by the kink rule (78 deg off the stem, 52 off its own first direction); it halts "came round
onto its own line" after 7 nodes, not yet looked at. The rig draws the open lane in orchid. Toes set: 91 lines,
1329 nodes, no errors. The user to judge the rest on the model.

User check on the model: "much better", the solve works; some lines got cut (the lane arm at arm 6 halts "came round onto its own line"), parked, not minded for now.

## 2026-10-26, night: every row cut at its valleys (the reference was the flaw)

The brow arm of line 1 (then L1 a21) drifted from a concave crease onto the convex one beside it, fully engaged by its
node 12, although the user could see the dip between them in the dots (row 6 across: 22, 14, 16, 20). The reading
found the dip in every row; what failed was "which crest is ours": the reference was the centre of the previous row,
and a row spanning two crests has its centre in the dip, so the side of the dip flipped row to row and the lane never
opened.

Agreed mechanism (the user's, in their words: no "start from", a dip is wherever two consecutive bins do not continue
the trend; the only direction that matters runs along the trail):
1. Every row is cut at its valleys, always. A valley is a drop followed by a rise anywhere in the across profile
   (over `lead` of the lower neighbouring crest; 0 = any dip).
2. Our piece is the one that continues ours from the previous row: the piece holding the strength-weighted centre of
   our piece of the row before. The node and the next reference come from our piece alone.
3. The other piece becomes a lane only once the rows ahead confirm it, as before (row k and the window rows ahead all
   showing an other piece on the same side). A dip no neighbouring row continues is a cut nobody follows.

Built in walk.py as `_cut_row` (cuts once per row, stored in `arm['cut']`, never moved), `_reference`, `_decide`,
`_divide_ahead`; `_crests`, `_own_part`, `_cut_at` are gone. A piece a lane takes is consumed (the first run gave the
same dots to two lanes and to the same lane twice). Rows' other pieces on a side with no lane stay painted but are
noded by nobody (parked).

Result on the chimp: the brow arm opens its lane at row 5 (rises 1.0, 1.0, 0.3, 1.0), the convex crease is an arm of
its own from there; the user's dots 46655/47190 stay with the concave crease, 46658/46653 go with the convex one.
Line 1: 36 arms, 749 nodes, 4 lane arms of 4-7 nodes (before: four one-row lane arms). Line 0 unchanged but for 7 nodes.

## 2026-10-26, later still: the trim at paint time (the far side of a valley ceases to exist)

The leap at L1 a14 n1: the paint of the line ending at L1 a7 n26 leaked sideways onto a neighbouring stripe of dots
near n24, crept along it by touch for two rounds, and when it parted from the crest the parting rule made it a branch
from n26, 8 cubes away. The node-time cut had kept those dots out of the nodes but never out of the front.

The user's statement, after a long misunderstanding on my side: any row, at any time, that spots a valley counts only the
dots up to the point the tendency shifts; anything past it ceases to exist for this line. Not its row, not its front,
not its clumps, not its parting. A second band is another line and means nothing to this one. The front shrinking on a
noise dip is harmless: the paint grows by touch in every direction and grows back the next round if there is where to.

Built: `_trim` in walk.py, applied in `_spread` to every freshly painted row (the main front, and each half of a split
under test). The dots past a valley are unpainted again (owner -1), free for another line. The lanes, `_cut_row`,
`_decide` and the window confirmation are removed. The view keeps the toggle "Show valleys and dots let go": orange-red
valley dots, cyan dots let go that no line has taken since.

Rig: side panel rearranged (small buttons three to a row, toggles two to a row, the log gets the height); Find pans
the camera to the target without changing its angle; the tool remembers the last model and dot set (Lines/last.json).

First run on the chimp (lead 0): line 0 21 arms, 469 nodes (was 28 / 556); line 1 177 arms, 6011 nodes, 38.9k dots
painted (was 36 / 749 / 4.9k), none of its long arms retracing earlier ones: the line was being choked by its own leaked
paint before. 544 trims over both lines; 1328 dots let go, 354 of them more than once. 49 arms of line 1 start with a
first segment over 4 cubes (to look at). The brow dots 46655/47190 and 46658/46653 are all on arm 9 now, rows 13-19.

## 2026-10-27: the paint grabs only dots past the front

L1 a10 climbed from n12 onto a crease above and behind, then n20 leapt 16 cubes ahead. Cause: at row 12 the front touched free dots
behind it (a little higher up), and from then on every row had a part ahead on our crease and a part behind on the other, 16 cubes
apart by row 18; the trim reads across the heading only and both parts fell on the same bins. The nodes, centred on the rows, drifted
back and up; row 19 held one dot ahead, hence the leap.

Rule agreed: a dot is grabbed only if it lies past the front (the strength-weighted centre of the row painted the round before) in
the direction of the walk. Dots behind stay free. Built in `_paint_from`. Line 1 went from 177 arms / 6011 nodes to 15 / 583 with no
segment over 3 cubes; it no longer spreads backward from every arm onto neighbouring creases. Not perfect yet; the model needs a look.

Rig: a "Stop at" box next to Find: a dot number; every step or run pauses right after the read that paints that dot (once per line),
the state line says where, and the next step or run goes on from there.

## 2026-10-28: the arrow turn-back, and the seed's direction

Resumed from commit 7220219 (the morning walk of 2026-10-27 with the rig rework); the plain-lines rework of the 27th stays
reverted. Looked at with the user: L1 a10 left its crease from n12 and sat on another by n20. In the rows there a second
crease ran three cubes beside ours, equally strong, with no dip in strength between them, so the trim saw one band and the
node (the row's centre) slid across as our dots faded. The arrows told them apart: across our crease their lean swung one
way, three cubes on it swung back.

Built (the user's rule): the trim also cuts a row where the arrows' lean, read from one edge of the row to the other,
turns back by more than the handle "Arrow turn-back (deg)" (`swing`, default 15; the two creases here swing about 50, two
neighbouring dots on one flank differ by about 7 from their triangles). The dots are ordered by their distance from the
heading line, height counted: the sideways-only distance of the strength trim was my shortcut and puts a dot and the one
above it at the same place. On the user's run the climb at that spot is gone.

Then L1838 a0 (a seed between two parallel creases) placed its first nodes between them: an arm had no heading before its
first node, and a row without a heading was never trimmed. Agreed: the seed node is not a node like the others, it has no
direction; it is given one from the arrows of its rings, and each arm leaves the seed with it as heading, so its rows are
trimmed from the first. Where a crossing crease lies in the rings the read is blurred (the sureness says so). The user:
climbing improved significantly, not fully resolved; on that line the nodes still move from the upper crease to the lower
between n6 and n10. Line numbers change between runs now (the trims change the earlier lines' courses); dot numbers do not.
Headless check of a line: run all lines with `walk.log` swapped in at the wanted seed (scratch of this session).

## 2026-10-30/31: other lines' dots, the join, the seed rule, the trim reference, the seed's dip

Started from the small lines "crawling around a crease" (L965 seed). Found: the behind-the-front rule leaves about one dot
in seven of a band unpainted (a dot one cube behind, one and a half sideways); of 4505 lines, 2756 seeded on such leftovers.

Agreed with the user over the two days, and built (walk.py, app.py):
- Other lines' dots are not invisible. The paint spreads over them (mask `_in`: the dots in this walk's paint, own or borrowed;
  `owner` keeps the line a dot was found on). A row made entirely of another line's dots halts the arm: that crease is walked.
  The old halt on touching another line's dot is gone.
- The join (handle "Join: owned share of the row", 0.5; "Merge within (share of node spacing)", 0.5): a node whose row's
  strength is mostly another line's merges into that line's nearest node when within `near` of its node spacing, else it is
  inserted as a new node on the nearest segment (named "L{line} joined by L{k}"). Ownership, not distance, tells parallel
  creases apart. A seed ring that is all another line's joins its seed node into that line and the line halts.
- The seed rule: a free dot whose neighbours' strength is mostly (> join) one line's is that line's dot, not a seed (lines
  4505 -> 2040 on the chimp at that point).
- The trim reference: the strongest trail dot touching the row (before: the centre of the row before, an average, which
  could sit in the dip between two crests; L403's seed node did). A reference that falls inside a valley takes the nearest
  crest: "the lowest of a dip is as sharp a cut as the highest of a crest".
- The seed's dip (31st, after L798 seeded on a weak bump beside line 396): rings grow in every direction with no heading, so
  the paint walked over the dip around the seed onto line 396's band, and the clump past the dip, handed to the arm as its
  first row, put the first node on line 396; the trim reads only the row, never the ring the dip lay in. Built: at the split,
  before the seed node, each clump is read along the direction from the seed to it with the valley read the trim uses across
  (now one shared piece, `_dip_cut`), the seed as the reference. A seed that is a peak keeps a clump with no dip on its path
  and lets go of one beyond a dip; a seed in a dip goes with the nearest crest and lets go of the rest. The dots let go are
  unpainted again; the dip dots themselves stay painted but make no row. One clump left is a plain line with one arm; none
  left halts the line at its seed node, "a bump of its own". Both cases looked at (seeds 47704 and 42843) now halt there.

Whole chimp set, seed cut off / on: 2277 / 3447 lines, 31370 / 26957 nodes, 2843 / 1459 joins, 1225 / 2187 lines without a
node, 769 seeds halted as a bump. The user: "expected a bit more of an improvement; hard to tell, so much is still broken".
Open: lead sits at 0, so a dip of one unit cuts; whether the 769 bumps are real or dips read too finely needs the eye. The
join's insert half makes zigzag "joined by" nodes the user judges broken (L1009/L1346, L1418); no decision yet. Dots left
behind the front are to be monitored. Headless runs and replays live in the session scratchpad (run_seedcut.py, baseline.py).
Line numbers differ between the app's run and the headless replay (e.g. 403 vs 570 for seed 42843); dot numbers do not.

## 2026-10-31 afternoon: the seed's dip read sideways

Looked at with the user: L52 a1 n1 (seed 27751) halted after one node, and seed 29660 beside it (strength 66) was declared "a
bump of its own"; between them a patch of painted dots no line walked, and line 115 that seeded in the patch halted on
"the ring is all another line's" with its node lost into the bump's empty record. Cause, both: the seed's dip read ran along
the path from the seed to each clump. A seed is the strongest dot of its patch, so along the crease the strength can only
fall from it, and the wobble along a crest (one to three points on 29660's two sides, twelve on 27751's) read as a rise after
a dip. The user's framing: the dip read at the seed exists to tell a clump in another crease from the two ends of our own
or a branch; across a band a dip is a second crest, along it a dip is noise; the seed has no heading for the row read.
Agreed and built: the seed's direction (the rings' arrows) is read before the cut, and the cut is one sideways read over
the painted dots along the across axis of that heading, the seed as the reference, the trim's valley read. A clump ahead
or behind sits at the seed's own across position and is never cut; a clump to the side past a valley is let go. Whole set:
3447 -> 2308 lines, 26957 -> 31325 nodes, 769 -> 14 bump halts, 2187 -> 1244 lines without a node, 1459 -> 2757 joins.
Seed 29660 now walks a line of 143 nodes that also takes 27751's patch. Open, raised by the user: a seed that halts without a
node keeps a line record and its painted dots (no later seed can start there, a join into it is lost); not changed yet.
Log fix: the spread's "belong to line" count now reads the row as the trim left it.

## 2026-10-31 evening: the join is the nearest node; the first node is the row through the seed

L50 a0 n9 (seed on leftover flank dots beside line 50's walked band): rows 100% line 50's were inserted and merged into
line 50 as "joined by" nodes, branches landing on its nodes, the whole system collapsing onto the joined node. The user's
reading: a trimmed row landing on a walked band is on a crease that was read whole and already has a node for that spot;
our row sees part of the band and cannot do better. Agreed and built: a row that, after the trim, holds any dot of another
line takes that line's nearest node as its node; nothing is inserted; the two handles that served the join (join, near) are
gone, `seedshare` keeps the seed rule. The walk goes on past the joint along its heading. Whole set: 2302 lines, 29661
nodes, 3324 joins, 14 bump halts. Found on the way: bands touch before trails meet, so consecutive rows can merge onto the
same node (452 merged nodes repeat their neighbour exactly, 97 turn the line back on itself).

L22 a6 n16 (seed 37196, line 368): a seed on a leftover flank dot a cube and a half beside line 22's crest; its first node
off the crest, every later node one of line 22's, drawn as a loop out of line 22 and back. The user rejected the rings as
anything but the gathering of dots until the heading is known: once it is, the dots around the seed are paint on the band,
and the row through the seed is a row like any other. Built: the first node is the slice of the ring dots one cube thick
along the heading at the seed, read as a row (across from all the ring dots, along from the slice, onto the surface), join
included: a dot of another line in it makes the node that line's nearest. The seed dot has no standing of its own. The
`_seed_node` ring centre stays only for the halt on a ring entirely another line's, where no heading is read. Seed 37196's
first node is now L22 a6 n13. Whole set: 2301 lines, 29645 nodes, 3661 joins; 351 first nodes joined. Found on the way, not
changed: a seed whose first node joins still walks its arms, and where every row lands on the same band the line is a
string of that line's nodes (263 lines have no node of their own, 228 before); the seed cut on that seed let go of line 22's
crest dots on a dip of 51, 43, 50 across, the crest drifting over neighbouring one-cube bins (the "ripple", open).

## 2026-10-31 evening, later: a stripe reads its sum per cube of length

L279 (seed 46985 beside line 18) dented like L368: the seed cut let go of line 18's crest on 52, 47, 54 across, so the row
through the seed held no line 18 dot and there was nothing to join. Walked through with the user: the cut reads the whole
patch (rings and clumps, several rows' worth) as one cross-section; a stripe across runs the whole length of the patch and
took its strongest dot as its value, so one dot of one row set a whole stripe and the waver along the crest read as a dip
across (the "ripple"). The user: "this is noise; a bit of dust on one edge reads as a bump all across"; and a stripe with
few dots should carry its share of nothing. Agreed and built: a stripe reads the sum of its dots' strengths divided by its
length; at the seed the patch is the disc the rings make (radius = ring count, the clumps the ring beyond), so a stripe at
distance d across is a chord of it, never under one cube; in the trim a row's stripe is one cube long, so its sum. The lead
stays at zero for now; the tolerance from the variation is the next step. Whole set: 2301 -> 2360 lines, 29645 -> 30122
nodes, 3661 -> 4062 joins, 14 -> 7 bumps; seed 46985 is now painted by line 18's own walk (no line starts there); seed
37196 still cuts on a 1% dip (82, 81, 105 across) but its first node joins line 22. The dips that cut are deeper than
before (median rise 0.26 of the crest, 0.07 before) yet with the lead at zero 54 of 455 seed cuts fall on a rise under 5%.
