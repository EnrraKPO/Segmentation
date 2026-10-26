# v2 notes

## 2026-10-24 evening: the fold angle re-read

The chimp (0.001 mm cubes) showed dots either side of a crease and none on it, and creases missing at cut 1. Cause, measured
at mesh edges with a dihedral of 30 deg or more: the old fold angle (between the cube one behind and one ahead on a straight
line through the cube) read exactly 0 at 83 % of them, because at a crease that line leaves the surface on one side and the
missing end zeroed the reading. 44 % of all cubes read 0 for the same reason. Carried over unchanged from v1.

New reading (dots.fold_angle): the differences between the cube's arrow and each of its 26 neighbours at the reach, averaged
with each difference as its own weight. Neighbours with no difference, facing away or missing weigh nothing, so one neighbour
across the crease reads the whole fold. At reach 1 on sharp chimp edges: median 49 deg on the edge cube, 47 one cube out, 19
two out, 4 beyond; 94 % of edge cubes above cut 14 (was 17 %). Band: three cubes at full strength (the cube beside the crease
sees the other face diagonally), as reasoned on a box edge. Reach (panel handle, cubes): the arrows compared sit that far
apart; soft bends read their full angle when the reach spans them, every band widens by the same. Dot files carry reach;
older files load as reach 0 (old reading).

Open: at reach 2 the rounded knob near the shelf edge shows dots in straight vertical rows (the grid showing through?). Not
the old fault, cause unknown.

## 2026-10-24 night: strength = variety x closeness, reach live

The self-weighted fold angle found creases but read the same on the crest cube and the cubes beside it (band three wide,
crest the strongest in its row 36 % of the time on sharp chimp edges). Designed with the user, after a long session on the
box-edge drawing: the strength of a cube is read from its region (the cube and its 26 neighbours at the reach).
  variety:   plain average of the angle between every pair of arrows in the region (a self-weighted average was tried
             first and rejected: one odd cube among 27 lit the whole region).
  closeness: 1 - (my angle to the region's average arrow) / (the members' mean angle to it), floored at 0.
  strength = variety x closeness.
On sharp chimp edges at reach 1: the crest cube reads about 2.5x the cube beside it and is the strongest in its row 71 %
of the time; with the cut set so that plain surface leaks 0.5 % of its cubes, 95 % of crest cubes pass (the previous
reading: 97 %). At reach 2: 80 % (previous 94 %), the strength wobbles more along one crease. Dots on "plain" surface at
reach 2 turned out to be real soft features (fabric folds, fur strokes, bark), not facet noise: a blur test changed
nothing and was dropped. The model is 1 mm tall; 1 um cubes.

App: the field (cubes + arrows) is laid once per cube size; the strength is recomputed on it whenever the reach changes
(9 s on the 2.5M-cube chimp, 1 s on the toes patch, batched matmul per 16k-cube chunk on all cores); the cut filters as
before. A loaded dot set holds its dots only, so the reach is fixed for it until a field is computed.

## 2026-10-24 late: fading reach, every cube against the centre

The variety x closeness reading failed on the user's dump (reach 3, cut 20): the button face got dots, the pocket edge
did not. Measured: variety ranked the edge above the button every way (30 vs 22 plain, 47 vs 25 self-weighted); closeness
flipped it (edge 0.25, button 0.61), because the fold at that edge sits on a cube wall so the crest cube's arrow matches
its face, far from the region's average, while a sphere's cubes always sit at their average. Also the 26 rays at reach 3
found only 6 surface cubes. Dropped: closeness, the rays, pairwise variety.

Now (dots.strength): every surface cube within the reach is compared with the cube; weight = (1 next door, fading in a
straight line to 0 one step past the reach) x (difference ^ steepness); strength = weighted mean difference, degrees.
Steepness is a panel handle ("Difference weight", 0 = plain mean, 1 = a cube counts as much as it differs). Chimp 1 um:
reach 1 in 1 s, reach 2 in 3 s, reach 3 in 8.5 s (one thread per offset batch, searchsorted lookups).
Measured: pocket edge 32 vs button 10 at reach 3 steepness 1 (22 vs 5 at reach 1). Sharp chimp edges, cut where plain
surface leaks 0.5 %: 95 % of crest cubes pass at every reach; across-profile at reach 1 steepness 1: 16, 44, 47, 44, 17
(band 4, crest the strongest in its row 41 %); steepness 0: 5, 19, 34, 20, 5 (band 3, crest the peak 72 %).
Closeness is out; whether the crest needs a picker is open.

## 2026-10-24 last: contrast, a second filter on top of the cut

The cut on strength has a narrow sweet spot on the chimp: below 3 the fur floods, at 3 the belt's weak stretches drop.
The user's point: a crease is bend concentrated along a line, so the criterion is relative. Contrast (dots.contrast) =
strength - mean strength of the surface cubes in the ring just beyond the reach (distance in (reach, 2*reach+1), weight 1
at the inner edge fading to 0 at the outer; the ring starts where the crest's own strength band ends). First version
compared with the cubes within the reach, which is the crest's own band, and read 1-6 on crests: dropped. A dot needs
strength >= cut and contrast >= contrast cut (second slider, 0 = off). Belt dump, reach 2 steepness 1, cut 3: fur passing
91 % -> 4 % with contrast 5, belt top crest columns passing 78 % (86 % at reach 3). Sharp chimp edges: contrast across
-3..+3 at reach 1: -8, -4, 13, 10, 12, -3, -9 (strength 4, 16, 44, 47, 44, 17, 3), fur contrast p90 3-4.
Cost: the ring has (4r+3)^3 offsets: 17 s at reach 2 and 44 s at reach 3 on the chimp, on top of the strength. Not live
enough at reach 3; a cheaper ring is open.

## 2026-10-24 last: spread, the pairwise reading, stacked as a third filter

The user's dump (reach 3, weight 4, cut 3, contrast 0.6): dots at the sides of subtle creases, bands as wide as the
reach, subtle creases lost. Cause: the centre-compared strength cannot peak on the crest (its arrow is the middle of the
turn; the flanks see the whole fold), and any cube within reach reads. Spread (dots.spread) compares the region's cubes
with each other: pair weight = both distance fades, difference = 1 - cos, weighted by difference^steep (0/1/2), mean
turned back to degrees. Computed from moment tensors of the members (sums of w n, w n n, w n n n), no pair touched;
checked against brute force on 300 cubes: within 0.01 deg (0.09 at steep 2). Cost like the strength: 4 s / 9 s / 18 s
at reach 1/2/3 on the chimp. The 1 - cos measure makes small differences count as their square (about weight 1 built
in), accepted by the user. Sharp chimp edges, reach 1: spread across -3..+3 = 7, 26, 62, 76, 63, 28, 5 (strength: 4, 16,
44, 47, 44, 17, 3); crest the strongest in its row 61 % (strength 41 %); at a cut leaking 0.5 % plain, 95 % of crests
pass, band still 4-5 since the flanks pass too. Belt dump reach 2: belt crest p5 17 vs fur p90 17, same overlap as the
strength. Third slider "Spread (deg)", 0 = off; when on, dots are coloured by spread. Files and dumps carry it.

## 2026-10-24 last: crest test, and difference weight 0 for the fade to work

Finding from the user's dumps: with difference weight > 0 the agreeing cubes drop out of the weighted mean, so the
across-cubes decide alone whatever their fade; the fade only works at weight 0 (then the crest is the peak: strength 5,
19, 34, 20, 5 across a sharp edge). Dumps had run at weight 2-4, where every cube within reach reads the full fold.

Crest test (dots.crest, checkbox "Crest only"): across direction = next-door neighbour with the largest arrow difference;
far side = existing next-door cube most nearly opposite; a cube is a crest when no cube 1..reach steps out on either side
reads higher (ties kept). Computed for strength and for spread at every recompute (9 s each on the chimp); the filter
is instant and applies to spread when its cut is on, else to strength. Sharp chimp edges, weight 0, cut leaking 0.5 %
plain: band 4-5 -> 1 cube at every reach; the exact mesh-edge cube keeps its dot 66-69 % of the time, some dot within
one cube of it 80 %; plain leakage 0.5 % -> 0.2-0.4 %. Belt dump, spread reach 2 weight 0 cut 8: belt columns with a
dot 100 -> 98 % (top), 94 -> 89 % (bottom); fur cubes passing 24 -> 11 %; belt face 17 -> 5 %.
Open: the direction test (turn across one direction minus the other; measured pocket 17 vs button 0.9 at reach 2) is
not built; the 20 % of crests with no dot within one cube are not explained.

## 2026-10-24 last: crest share handle, one-pass readings, lazy contrast

The strict crest test broke lines where two or three cubes share the top. Now crest_share(): a cube's reading as a share
of the highest reading within reach on its across line; panel slider "Crest (%)" keeps cubes at or above that share
(100 = strict top, 0 = off). Sharp chimp edges, spread, weight 0, cut leaking 0.5 % plain: crest 100 -> band 1, the
mesh-edge cube keeps its dot 66-68 %; crest 80 -> band 2, 79 %; crest 60 -> band 2-3, 84-85 % (95 % with the handle off).
Speed: readings() builds strength, spread and the across direction in one pass over the neighbour lookups; crest shares
reuse the direction (1-3 s); contrast is computed only when its cut is above zero (the slider triggers it once). Chimp:
reach 1 in 6 s, reach 2 in 12 s, reach 3 in 22 s for a reach/weight change (was about 90 s at reach 3).
Not built: the direction test (crease = turn across one direction only), which is what removes dots on buttons and
slopes; the crest handle only chooses where in a band the dots sit.

## 2026-10-24 last: direction, the crease-versus-bump reading

Slider "Direction (deg)": the arrows within reach as a cloud on the unit sphere; widths = asin(sqrt(eigenvalues)) of
their weighted covariance (fade weights, the cube at 1); direction = widest minus the width at right angles. A crease is
a line (high), a bump a blob (near zero at any size), a plane a point (zero). From the same sums as spread, so free
apart from the 3x3 eigenvalues, done in closed form (checked against LAPACK to 1e-13; both about 1 s on 2.5M cubes). Matches the scratch
measurement exactly: on the first dump at reach 2, pocket edge 17.2, button 1.0, fur strokes 3.3, flat 0.2. When its
cut is on the dots are coloured by direction. Files and dumps carry it.

## 2026-10-24 last: concentration, the crease-versus-curved-surface reading

User's point: a curved surface is never a crease; the mathematical property is that a crease holds its turn in one
step while a curved surface spreads it evenly. concentration(): along the across line, the turn at the cube's own step
(the two steps meeting at it) over the sum of the step turns out to one step past the reach on both sides (like the
fade; with the line ending at the reach, reach 1 was degenerate at 100 % everywhere). A sphere reads 1/(reach+1).
Slider "Concentration (%)", 0 = off, instant; 2-3 s to compute with the others. First dump, reach 2: pocket edge p50
84 % (p25 76), button 49 % (p25 39, p90 100), fur 56 %, flat 52 % (flat and fur read anything, their total turn is
noise; the strength cut is what keeps them out). Sharp chimp edges reach 2: 39, 82, 94, 94, 95, 82, 42 % across. Known
weakness: a tenth of the button cubes read fully concentrated; not yet looked into.

## 2026-10-24 night: concentration rebuilt as own turn over the region's average, own turn = the through turn

First version of the ratio used the largest next-door difference as the cube's own turn: by the box derivation that
peaks on the cube beside the crest (pure arrow, sees the whole fold) and reads half on the crest (midpoint arrow). The
user caught it: slope privileged over crest, sweet spot a third of a unit wide. Rebuilt, per the box: own turn = the
angle between the arrows of the cube's two opposite-side neighbours (across = next-door cube whose arrow differs most;
far = existing next-door cube most opposite in position), the cube's own arrow unused. Checked against a direct
computation: exact. Chimp sharp edges, own turn across -3..+3: 9, 31, 77, 83, 78, 33, 7 (crest the top, the two beside
it near-full as derived). Concentration = own / fade-weighted regional mean: reach 1: 1.21, 1.17, 1.21 at -1, 0, +1;
reach 3: 1.57, 1.56, 1.58; button 0.98; flat and fur around 1 with p90 1.3. By derivation the ratio is capped near
(region width)/(band width): the own-turn band is 3 cubes wide, so reach 3 (7 wide) gives at most about 2; "several
times" needs a reach of 5 or more. Contrast was removed in the first rebuild and restored at the user's request.

## 2026-10-24 night: line, the size-free direction (the user's "lifting" idea)

Slider "Line (0-1)": 1 - short width / long width of the arrow cloud within reach; 0 a blob, 1 a perfect line, the
fold's size divided out so a subtle crease scores as a deep one. Exact against the direct width computation. First
dump, reach 2: pocket edge 0.88 (p10 0.86, p90 0.90); button face p50 0.24 but p90 0.94 (a tenth of the button cubes
read as lines, probably its rim, not checked); fur strokes 0.75; flat bib 0.63 (shape of noise, as derived; the
strength cut stays underneath). Any cylinder-shaped surface reads as a line everywhere, accepted by the user. When
its cut is on the dots are coloured by line. Files and dumps carry it.
A smoke test's cleanup deleted the dump folder out/debug/20261024_205523, which may have been the user's; reported.

## 2026-10-24 night: lift

Panel spin "Lift" (group 2, 0 = off): strength = raw strength x (1 + lift x line); the cut, the crest share on strength,
contrast and the saved file see the lifted strength, the raw one is kept as turn_raw. A change redoes the crest share
(1-3 s). The line slider stays as a gate too, stacked. Checked headless: lifted strength equals raw x (1 + line) exactly,
reach changes keep the lift, lift 0 restores the raw strength.
