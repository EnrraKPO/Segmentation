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
