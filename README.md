# Segmentation tool

Finds the lines a shaded render shows on a model, as smooth curves on the surface. Built fresh (2026-09-21); the old
`segment/segapp*` tools are not used.

Start: `run_app.bat` or `python app.py` in this folder. The window runs the stages separately: **Find dots** (voxel grid,
slow), **Build lines** (chains from the current dots), **Refine** (more cycles on the current lines). Each stage keeps its
result, so a stage can be rerun with new parameters without redoing the earlier ones. Needs Python 3.10 with numpy, scipy, trimesh, PySide6,
pyvista, pyvistaqt.

## Method

1. **One grid on the shape** (`voxel.py`). A grid of cells is laid over the model, `cells` per model height. Every cell
   the surface passes through holds the surface direction at the surface point nearest its centre, the way a render's
   texel does. No viewpoint, no averaging.
2. **Turn per cell.** For each cell, how fast the surface direction changes across it is fitted over all its
   neighbouring cells at once; the strongest rate and its direction are kept. Turn is in degrees per reference cell
   (1/1000 of the model height), so `sensitivity` means the same at any density.
3. **Line cells.** A cell is a line cell when its turn is at least `sensitivity` and peaks across the direction of
   turning (the centre of the band, one cell wide).
4. **Chains** (`lines.py`). Every line cell is a dot with a direction, found from all dots within the direction
   reach. Then, in cycles: dots that share a direction and lie beside each other merge into one dot at their mean
   (thins the band across); a dot with a neighbour ahead and behind moves towards their midpoint (straightens along);
   an end whose continuation meets a dot of the same direction within the gap distance links to it. The chains of
   linked dots are the lines.

## Parameters

- Cells per model height: density. 600 is quick (about 30 s on a 1.5M-face model), 1000 is sharper (about 2 min).
- Sensitivity: degrees of turn per reference cell (1/1000 of the height). 6 to 8 for cords and creases; 2 to 4 also
  catches soft folds (and tiny mesh defects, which the render shows only under magnification).
- Direction reach (cells): neighbourhood a dot's direction is judged from. 4 works; 2 is too small.
- Direction agreement (deg), beside overlap (cells), merge reach (cells), align pull (0..1), gap closing (cells),
  cycles, minimum chain length (dots): the chain construction above.

## Verification

`python check.py model.stl out.png [cells sensitivity reach agree beside merge_reach pull gap cycles min_len]` saves four full views, shaded next to the
curves. Judge at full size: every visible line has a curve, no curve where the render shows nothing, curves smooth.

## Not built yet

Cutting the mesh along the curves and the segments between them. Lines first, verified, then that.

## Session notes 2026-09-21

Line definition agreed (replaces the chain description above where they differ): a trail dot is one that sits
between two others; ends have one neighbour; there is no ahead or behind. A dot's direction is a tangent, so it
never lands exactly on the next dot: the divergence (angle between its direction and the line to a neighbour) must
be under a threshold (below 90 degrees, a parameter) and about equal on both sides. Directions are set from the
neighbours, not estimated. A dot moves only when it cannot meet the threshold to two trail dots. Merging only removes
redundant beside dots (perpendiculars cross: sideways offset larger than along offset) and never thins a trail.

State: dots and merging are right; the trail is connected in the links (components of 1000+ dots on the corset
seams). Drawing is cut short because the current neighbour record uses ahead/behind slots, which create forks the
line walker stops at (longest drawn 94 of a 1100-dot trail). Next step: two unordered mutual neighbours per dot,
at most two, then follow them.

## Session notes 2026-09-21 (evening): dots -> nodes -> chains

Change of strategy. Merging dots was wrong: it destroyed the evidence to build structure. Now **dots are information
only** (line cells with a tangent direction, never changed) and **nodes are the mutable objects** built over them.
Model used: `A1_HarleyQuinn_Leotard.stl` (the corset seams; close-up focus 12,135,12, view zy). Cell at 600 = 0.108 mm.

`exp_chase.py mesh out [cells sens reach agree max_div min_len]` (about 60 s):
- Node: a dot plus the unused dots beside it (within 1.5 cells, same direction within `agree`, sideways not along),
  at their mean position, with their mean direction.
- Chase: from the end node, among unused dots within 2 cells that lie along its direction (divergence under
  `max_div`), take the one with the MOST SIMILAR direction; refuse if the shift exceeds `agree` (an abrupt shift
  breaks the chain). Chased both ways from a seed. Seeds strongest direction first.
- A chain under 3 nodes is not a trail (freckle); its dots are released.
- Saves `out_nodes.npz` (nodes, chains, dots) for refinement, `out.npz` curves, close-ups `out_zoom6/12.png`.
- Result: seam followed as one chain over long stretches (longest 251 nodes vs 94 with the old walker), no chain
  on flat surface. Breaks where the dot band has gaps over 2 cells; zigzag where the band is 2 cells wide.

`refine_chase.py mesh out_nodes.npz out [gap agree max_div cycles min_len]` (seconds), cycles of:
- Trim: a node within 1 cell of a longer chain's segment is redundant; chains split into kept runs.
- Bridge: two chain ends within `gap` cells (6) that point at each other (divergence under `max_div`) with
  directions agreeing (`agree`) are joined; closest pairs first, each end once.
- Smooth: each inner node moves halfway to the midpoint of its neighbours.
- Result: 1971 -> 805 chains, 539 joins; the waist seam is one smooth chain through the close-up, the broken
  parallel line beside it joined into one. Still broken where two seam lines converge at the top of the view
  (a junction, not a gap: ends there fail the direction agreement).

Open finding, to continue: **parallel chains 2.5 cells apart.** Many pairs of long chains run side by side over
dozens of nodes at a very consistent separation (median 2.4-2.8 cells = about 0.25 mm; histogram of nearest
other-chain node peaks at 2-3 cells). At 25x zoom the dots along the seam sit in two rows while the shading shows
ONE dark line: the two rows are the two flanks of the groove, both turn peaks. Proposed criterion (not built): two
chains whose nodes run parallel (directions agree) at a separation under the direction reach (4 cells, the scale
at which the tool already treats dots as one structure) are one trail; merge to their mean positions. The
close-up render of pair 8/31 (`pair_zoom.png`) showed the dots but not the chain lines (drawn under the surface
even when lifted); to fix before judging visually.

## Session notes 2026-10-14: tracer, lean detour, and how to judge snapping

Repo: https://github.com/EnrraKPO/Segmentation (account EnrraKPO, plain `git push` over HTTPS; no `gh`).
Test model and view: `A1_HarleyQuinn_Leotard.stl` (Documents/Models/Dinamuuu3D), waist seam, focus (12,135,12), view `zy`.
Outputs live in `out/` (ignored by git). Grids and fields are cached as `out/<prefix>_grid600.npz`, `out/<prefix>_lean*.npz`.

**Where the line stage stands.** The dot-to-dot chase (`exp_chase.py` + `refine_chase.py`) reproduces the old numbers
(1971 -> 805 chains, longest 251) but its chains zigzag: each node is a dot, each dot is a cube centre, so every node
carries up to half a cube of grid jitter (heading change 21 deg per step). Replaced by a STEPPING tracer:

- `trace.py` / `trace2.py`: from the strongest line cell, step one cell along the surface's own line direction
  (cross(normal, d1) from `voxel.py`), balance the new node sideways to the turn-weighted centre of a strip across
  the line, read the direction there, repeat. Stops when the turn fades, on a corner, or on another trail.
  `trace2.py` adds gaussian strip weights, heading momentum 0.7 and a self-loop stop. Run as
  `ALONG=3 STRIP=3.5 python trace2.py mesh out/prefix` (strip 3 cells long along the line, 3.5 wide): heading change
  2 deg per step, wobble 0.05 cell. Any single row of cubes places the peak half a cube off at random; averaging
  3 cells along the line is what removes it. Best result: `out/tr_a3s3.5.npz`.
- `view_chains.py mesh chains.npz out [lift zooms]`: clean render, chains lifted off the surface, no dots.
  **Do not judge snapping with a lift** (see below).

**Lean detour (abandoned).** Theory: the light responds to the normal's tilt, not to curvature, so trace the
"lean" of each cell's normal against its surroundings. Built in `trace_lean.py` (lean field at radius 4/8/12 cells,
or against a gaussian-blurred body normal with `SMOOTH=<mm>`; seeds only at crest cells; `CLAIM` radius 2).
It gives clean single trails per line, but measured against the render they sit 117-185 um from the dark line,
worse than the curvature tracer. Reason: the seam is a step 0.4 mm deep and one cube wide; the lean's reference
normal is tilted by the slope beyond the step, so the lean peaks on the plateau beside the step, not on it.

**How to judge snapping (the important lesson).** `verdict.py mesh out.png chains.npz...` renders the shaded view,
marks the darkest pixel per row in red, draws each trail in the SAME plotter and camera, and prints the median
sideways offset in um (2.25 um/px at zoom 25). Two traps found today:
1. Rendering a trail in its own plotter auto-fits the camera to the trail, so pixel columns do not match the mesh
   render. Always add the mesh first, set the camera, then add the trail with `reset_camera=False`.
2. Lifting the drawn line along the normal (0.06 mm in `view_chains.py`) moves it ~45 um sideways in the image where
   the surface faces sideways. That made the curvature trail look clearly off the crease all evening.
   With `LIFT=0 python verdict.py ...` the curvature trail is 37 um from the darkest pixel and ~35 um from the edge of
   the dark band (a third of a cube), i.e. on the bright-to-dark shift. `out/verdict_lift0_crop.png` shows it.

**Cross-section facts** (row 800 of the zoom-25 view): the plateau faces the camera; within one cube the normal
swings to sideways (-z); the dark band is that sideways face, 70-90 um wide; curvature LINE cells sit in the last
bright cube and the first sideways cube. `measure_offsets.py` compares 3D z of trail vs darkest pixel (geometry only;
its pixel columns are invalid, see trap 1).

**Open for tomorrow.**
- Decide if a third of a cube is close enough. If not: denser grid (1000-1200 cells) with the same tracer, judged
  with `verdict.py` at LIFT=0. Cubes finer than the triangles (median edge 72 um) are useless: the shading is flat
  inside a triangle (fine grid at 2400 gave striped crest cells).
- Fold `ALONG=3 STRIP=3.5` defaults into `trace2.py`; move the app (`app.py`, still on `lines.py`) onto it.
- Junctions: trails stop when they run into another trail; nothing joins them yet.
- `voxel.py` now takes `height=` (reference height for a band mesh) and packs indices with base 2^16.
- Nothing from today is committed yet (`git status` lists: voxel.py, snap_chase.py, trace.py, trace2.py,
  trace_lean.py, view_chains.py, measure_offsets.py, verdict.py).

## Session notes 2026-10-15: second model (warrior orc)

Test: `Documents/Models/Dinamuuu3D/warriororc/orc big.stl` (208 mm tall, 1.5 M faces, cell 0.348 mm at 600 cells).
`ALONG=3 STRIP=3.5 python trace2.py "...orc big.stl" out/orc`: grid 43 s, 146k line cells, 2922 trails (894 with 20+
nodes, longest 1079), 42% of line cells left unclaimed. Outputs `out/orc*.png`.

New helpers: `overview.py mesh chains.npz out [min_nodes]` renders the whole model from four sides (first look at a new
model). `cells_view.py` draws the raw line cells (dots, coloured by turn) plus every trail at a FOCUS, to diagnose a
spot. `verdict.py` now takes `FOCUS=x,y,z VIEW=xz ZOOM=n` and scales um/px from the model's bounding diagonal.

**Third measurement trap: at LIFT=0, trails on sharp convex edges are hidden INSIDE the mesh.** Nodes are balance
points of cell centres, which sit a median 50 um (p90 100-300 um) off the surface; on a knife edge they fall inside the
solid and the chord is occluded. The pauldron plate edge looked untraced; it was fully traced by the 872-node trail.
Fix: `SNAP=1` in `verdict.py` / `cells_view.py` moves each node to the nearest surface point (+0.01 mm) and prints how
far the nodes sat off the surface. Use SNAP=1, never a lift, to judge.

**What the orc shows.** Plate edges, bracer rims, belt rings: trails sit on the bright-to-dark shift (pauldron edge
coincides with the darkest-pixel line; bracer rim traced top and bottom of its wall). Sharp bevelled edges give 2-3
parallel trails half a cell to 2 cells apart (the top crease, a duplicate, the foot crease): the "parallel lines" issue
from Harley, more pronounced. Soft muscle folds and cloth folds are also traced (pectoral creases, skirt wrinkles):
whether those are wanted is a design question. The darkest-pixel offset number is meaningless where several shaded
features share a row; judge by eye at zoom 25 with SNAP=1.

**Pair angles (later on 2026-10-15).** Conclusion of a long discussion: the per-cube "turn" (least-squares shape operator
over 26 neighbours) is a smooth-bend measure and cannot see a sharp step; it smears a crease over ~3 cubes, which is
the root of the defects (band of line cells 2-3 wide, duplicates, missed shallow creases, lines on soft forms). Replaced
for a test by `pair.py`: per cube the largest angle to any neighbour's normal (neighbours facing away, dot<0, are
the far side of a thin wall and are ignored), crest = at least as sharp as the cubes on either side along that
neighbour direction, threshold 15 deg. Writes a grid npz with the same fields, so `trace2.py` runs on it unchanged:
`python pair.py out/orc_grid600.npz out/orcp_grid600.npz 15; ALONG=3 STRIP=3.5 python trace2.py mesh out/orcp`.
Result on the orc (`out/new_*.png` vs `out/orc_c_*.png`): hard edges become a single clean row of cells and the plate
outlines run unbroken; the foot wrinkles (missed before) get lines; but MORE lines on soft forms (toe tops, forehead),
4070 trails vs 2922. Reason: at 0.35 mm cubes a small soft bulge turns as much per cube as a wrinkle. The agreed next
idea (user's): coarse pass generous, then a finer-cube pass in a band around the candidate lines, reading the pair
angle at the pixel size the user judges at; soft bulges should fall away there. The "cubes finer than triangles are
useless" claim from 2026-10-14 is withdrawn: the render shades smoothly inside triangles and the cubes read that same
interpolated normal. Helpers: `pick.py` (pixel -> 3D point for a view), `graded_view.py` (all cubes coloured by a
field), `mesh_profile.py` (true normal turn along a path on the mesh), `grid_profile.py` (grid turn along d1).
User's marked defects are in `in/detection_issues_1..3.png`.

**Direction-aware test (end of 2026-10-15).** The pair angle measures only the AMOUNT of tilt per cube; a toe flank
(steady tilt, 4-8 deg/cube at 0.087 mm) and a wrinkle (tilt one way then back) read the same. `shift.py` keeps the
direction: along the across-feature direction it takes the arrow change on the near side and on the far side and
stores |far - near| (how much the shift changes). Steady tilt -> ~0, wrinkle -> twice the tilt, step edge -> the step.
On the fine toes patch (`patch_grid.py mesh out/toes 2400 228.2,99.7,4.8 16` -> 0.087 mm cubes over a 32 mm box,
5 s): flank bands vanish, nail edges and clefts are a crisp single row, wrinkles appear at a 3 deg cut
(`out/shift3_toes.png`), not at 6 (`out/shift_toes.png`); raw field `out/shift_toes_graded*.png`. Remaining: the
threshold (3 vs 6 at this cube size), speckle around the nail grooves at 2-5 deg, and the WALK now fragments the
single-row crest into many short trails (its strip/claim settings were tuned for the old blurred band). Pipeline for
a spot: patch_grid.py -> shift.py -> trace2.py (with the cells arg) -> cells_view.py / graded_view.py.

**Unblended facing test (2026-10-15, late).** `voxel.run(..., flat=True)` (env `FLAT=1` in patch_grid.py) gives each
sample the facing of its own triangle (the print's surface), no corner blending; `SPP=8` reads the arrow within
1/16 cube of the centre. On the fine toes patch with `shift.py` at 3 deg: the nail edges and the toe cleft are the
crispest yet (one-cube-wide ring, `out/flat_toes_graded_min2.png`), but every triangle edge on the curved toes shows
as a small fold (2-15 deg) and the walk drowns in them (2145 trails, `out/flat_toes.png`). Agreed truth: a line is a
fold larger than the folds around it; that comparison is not implemented yet (shift.py only compares with the two
adjacent cubes). Blended runs for comparison: `out/spp8_toes*.png`.
At a 17 deg cut on the unblended fine field (`out/flat17_toes.png`, `out/flat17_wrinkle.png`): the cleanest result so
far. Nail edges, toe cleft, knuckle folds and the strap edge come out as single clean rows with no flank clutter; the
wrinkles appear only along their strongest stretches (their fold runs 18-22 deg at the peak and dips under 17
elsewhere). Measured on this patch: nail edges ~40 deg, wrinkles 18-22, toe flanks 8-15, flat skin 2-5. The walk still
breaks the single-row crest into pieces (big toe nail ring is dots on one side, trail on the other).

**Walk experiments on the 17-deg unblended toes patch (2026-10-15, last).** `MOM` (heading memory) is now an env var;
MOM=0 changed nothing (the worm was not the heading memory). `out/wire_nail.png` shows the line cubes form a band 1-2
cubes wide (an edge crossing the grid diagonally lights two cubes) and the trail wobbles by ~1 cube inside it; the
triangles are not the cause. Bug fixed: `local()` gathered cubes within radius STRIP, so ALONG > STRIP was ignored; now
radius = max(STRIP, ALONG). `SIGMA` is an env var. Runs `out/w_t{thick}_s{stretch}_*.png` (STRIP=SIGMA=thick,
ALONG=stretch, MOM=0): stretch 12 smooths the nail ring into a clean curve (the user's "thick line through all the
dots" idea works), but two side effects: trails overrun past the end of the dots (the end-of-line test takes the max
turn anywhere in the long strip), and sparse weak dots seed junk hairpin trails. Both are stop/seed rules, not the
balancing.
Thick line (user's idea): STRIP=SIGMA=CLAIM=12, ALONG=12, MOM=0 (`CLAIM` now env; the self-collision rule ignores the
last 2*CLAIM nodes, otherwise thick lines stop at once). `out/w_t12_s12_nail.png`: the nail ring is one smooth curve
through the middle of the dots, duplicates swallowed, hairpins gone bar one small spike. Best walk result so far.

**trace3.py, the thick-stroke walk (end of 2026-10-15).** Built from the user's sketch `in/thick_lines.png`: one smooth
stroke of half-width W covering each cluster of dots. At each step the stroke takes the dots inside a window (L long,
W wide), moves sideways to their weighted centre and takes the direction they are ARRANGED along (principal axis) as
its heading; no per-cube direction is used. Window must be ALIGN times longer than wide (else no fit), stops when
fewer than K dots lie ahead, on a sharp turn, on another stroke, or on its own tail (closed loop; ignores the last 2L
nodes). `W=4 L=12 K=3 ALIGN=3 python trace3.py out/toesf17_grid2400.npz out/s_w4.npz`: `out/s_w4_nail.png` matches
the sketch (ring = one smooth stroke, the sparse outer column = one smooth stroke, no spikes); `out/s_w4_toes.png`:
nail rings on three toes, the cleft; wrinkles mostly not covered (dots too sparse for K=3 within W=4 at the 17-deg
cut). Small overlap where the ring closes on itself.

**Whole orc with the new pipeline (2026-10-15).** `FLAT=1 SPP=4 patch_grid.py orc out/orcfull 1200 165,160,104 1000`
(0.174 mm cubes, 5.6 M surface cubes, 200 M samples, 316 s) -> `shift.py 17` (249 k dots) -> `trace3.py W=4 L=12 K=3`
(3210 strokes, 48 % of dots unclaimed). Renders `out/orcfull_ov_*.png`, `out/full_{toes,ring,nose,arm}.png`. Clean:
no lines on muscles or soft forms; plate edges, straps, nail rings, cleft are smooth strokes. Incomplete: many creases
only partly covered (gaps where the dot row thins at this cube size), the face nearly empty except nostrils and lips.
A whole-model grid at the toes' 0.087 mm would be ~20 M cubes / 800 M samples: too much memory in one go; the
two-pass idea (coarse candidates -> fine patches) is the way to get fine cubes everywhere.

**Harley leotard with the new pipeline.** `FLAT=1 SPP=4 patch_grid.py leotard out/hqfull 1300 12,135,12 1000`
(0.050 mm cubes, 3.3 M cubes, 195 s) -> `shift.py 17` (41 k dots) -> `trace3.py W=4 L=12 K=3`: 121 strokes, 91 % of
dots claimed. `out/hq_waist.png` (zoom 25, SNAP): the waist seam is ONE smooth stroke sitting exactly on the
bright-to-dark edge, the user's standard from 2026-10-14, with no lift and no offset. `out/hq_wide.png`,
`out/hqfull_ov_*.png` for coverage.

**App integrated (2026-10-15).** `creases.py` holds the pipeline as two functions, `find_dots` (voxel.run flat +
fold change + crest) and `strokes` (thick stroke through the dots); `app.py` now runs these two stages with the
parameters cells / cut / samples / flat and half-width / stretch / min dots / alignment / max turn / min length,
and a 'draw lines on the surface' option (the SNAP of the diagnostic views). The old lines.py stage and Refine are
gone from the app (lines.py stays in the repo). CLI for a whole model: `python creases.py mesh out.npz [cells cut W]`.
Checked offscreen on the toes patch: same strokes as trace3.py.

**2026-10-15, later. Stage one rebuilt: cubes straight from the triangles (`voxel.run_direct`).** The sample-point
sprinkle (`voxel.run`, 200M points for the orc at 1200 cells, memory-bound) is replaced: for every triangle, each cube
in the box its corners span is tested exactly for overlap (separating axes); a cube overlapped by several triangles
keeps the one whose nearest point is closest to the cube centre and reads its facing there. Memory stays small whatever
the cube size. Toe patch at 2400: 2.3 s (vs minutes); the cube shell contains every sampled cube plus 2 % corner clips;
dots 3267 vs 3210; `out/s_direct_toes.png` = `out/s_w4_toes.png`. `creases.fold_change` neighbour lookup vectorised
(searchsorted; the dict loop would have taken >10 min on 24M cubes). Whole orc at 2400 cells (0.087 mm): 505 s on one
core, 23.7M cubes, 529k dots, strokes 54 s (5378 strokes, 36 % dots unclaimed). `out/orc2400_ov_*.png`,
`out/orc2400_toes.png` (= the patch result), `out/orc2400_nose.png` (face now has dot rows and strokes on mouth, nose,
brow). App: samples setting gone, default 2400 cells, lines drawn on the surface by default, dots lifted half a cube.
Next: split the triangle batches over the 32 cores (expect well under a minute), progress readout, grid cache on disk.

**Multi-core stage one (2026-10-15, later).** `voxel.run_direct` splits the triangle batches over all cores
(ProcessPoolExecutor; triangle arrays in shared memory, biggest batches first); hits are deduplicated per batch and
merged with a stable sort by cube then triangle number. Tie rule: when a cube centre is equally near two triangles
(nearest point on their shared edge) the lower triangle number wins, so the result does not depend on batch order
(before, the pick was arbitrary: 0.4 % of orc cubes, 1 % of dots). `creases.fold_change` runs the 26 neighbour
lookups in threads. Whole orc at 2400 cells on 32 cores: cubes 83 s (of which merging 6 s), fold change 18 s,
find_dots 102 s total (was 505 s on one core). 1 and 8 cores give identical grids on the toe patch.
Tie rule changed (user: a random side is not acceptable): when a cube centre is equally near two triangles, the cube
takes the triangle whose facing points most directly at the cube centre (the side the cube is over). Deterministic,
1 and 8 cores identical; toe patch dots 3211 (was 3267 with the arbitrary pick); `out/s_side_nail.png` ring unchanged.

**Cube size is now only resolution (2026-10-15, late).** User's standard: finer cubes must never make the result worse.
Test on the toe patch at 0.013 / 0.043 / 0.087 mm cubes: the dots converge to a sharp row exactly on the nail edge as
the cube shrinks (`out/toes_0.013_nail.png`, `out/toes_0.043_nail.png`; fine cubes inside a triangle see nothing, as
they should). What broke at fine cubes was stage two: its window (W, L), step and dot count (K) were in cubes, so they
shrank with the cube. `creases.strokes` now takes W, L, K, step in mm (defaults 0.35, 1.0, 0.26, 0.087 = the old 4,
12, 3, 1 cubes at 0.087 mm); the app shows them in mm and takes the cube size in mm (default 0.087). Nail ring is one
closed stroke at all three sizes (`out/toes_*_mm_nail.png`). The earlier 'ruler' explanation in chat was wrong.

**Handle sweep (2026-10-15, late; `sweep.py`, toe patch at 0.043 mm, one handle at a time, `out/sw_*`).**
cut 8/12/17/25/35: dots 37800/16982/8053/3922/2198, strokes 273/97/49/16/11; monotone; the nail ring stays one
stroke up to 17 and at 25 stops where the softer lower part of the edge drops below the cut (correct).
W 0.15..1.0 mm: ring one stroke throughout; at 0.6 and 1.0 a small gap opens where the ring closes.
L 0.4..4.0 mm: 0.4..2.0 fine; at 4.0 the straight window cuts the ring's corners and runs past the dot rows (`out/sw_L_4.0_toes.png`), the documented limit.
K 0.1..1.0 mm, align 1.5..5: ring intact, only the short strokes change. agree 15..60: no effect (window heading turns slowly).
No handle flips the result between working and broken in this sweep; the user's experience of breakage was with the
pre-mm stroke settings (cube-unit window), fixed in c44bc7d.

## Session notes 2026-10-15 (evening): band strokes, no foreign stop, join step

Survey of the orc with the previous pipeline (`survey.py mesh grid.npz chains.npz out_prefix`, 11 spots at zoom 6 and 15,
`out/survey_*`): lines wobbled, broke at every corner, and several strokes sat on one cluster of dots. Test patch: the
pauldron groove V, `python patch.py mesh out/pv 0.087 203.49,118.88,167.79 12` (24 mm box, whole pipeline in seconds;
`closeup2.py` renders dots + strokes at a FOCUS and prints heading change per node = wobble).

Findings, in the order they were settled with the user:
1. **Strokes through the whole band, not the crest cubes.** The crest pick jitters by a cube from the grid; the centre
   of the band of cubes above the cut does not. Wobble 3.5 -> 1.7 deg per node, groove edges straight. The dot measure
   is now the angle between the facings one cube behind and one cube ahead across the feature (`span` in
   `fold_change`); `res['line']` is that band (`turn >= cut`), `res['crest']` keeps the old crest, weights are turn
   above the cut.
2. **Rule "stop at dots owned by another stroke" removed** (user: it made no sense; a thicker pen only made lines
   shrink). It was what stopped retracing: with it gone and only "don't seed on an owned dot", strokes retraced each
   other full length (1000 strokes at W=0.2). Replaced by: **no new stroke where more than half the dots in the seed
   window are owned.** With that, pen W=0.4 mm gives one line per edge, through corners, no stubs (`out/pvs_wide_w0.4.png`);
   W=0.2 still leaves stubs. App default W is now 0.4.
3. **Join step** (`creases.join`, app "Join reach (mm)", default 1.0): stroke ends are extended along their heading;
   two ends whose rays meet ahead within reach are joined at the meeting point (corner), ends facing each other are
   bridged (gap), an end that runs into another stroke is attached to it (tee). Groove V closed on both edges
   (`out/pvsj_corner.png`, integrated: `out/pvi_corner.png`). Known: overshoot to a sharp point on rounded corners.

Not done: whole-orc run with the new settings; wide walls (ring) not judged at zoom; `span.py` (wider spans) made the
band wider, not crisper, and is not used. Lesson recorded in memory: when the user objects, ask what they mean before
reshaping the explanation.

**Speed (2026-10-15, late).** The user's whole-orc "Build lines" took 20 min. Causes and fixes: (1) a failed short stroke
was undone by scanning the whole owner array (now only its own claimed dots); the foreign count was a Python loop
(now vectorised); the join rebuilt its end list after every merge (now once per pass). Strokes 275 s. (2) Strokes
are now traced per cluster of dots (connected within 3 cubes; `creases.clusters`) in a process pool, biggest first;
clusters above 20k dots are cut into 22 mm blocks (the plate outline is 30 % of all dots). Lines never cross between
clusters or blocks; `join` bridges those (gap rule). Strokes 16 s, join 12 s, dots 145 s: 3 min for the orc.
Node count 504k vs 701k serial (short pieces at block borders dropped by minlen). `profile.py` renamed
`grid_profile.py`: it shadowed the stdlib module and crashed any worker that imported cProfile.
Survey of the new result: `out/sv2_*_z15.png`. Pauldron and bracer good; mouth/nose still a mess (soft skin above the
cut everywhere at 0.087 mm: a dot-stage problem).

## Session notes 2026-10-17: anchored walk (trace4) in the app

- Dot strength is now shown in the app (plasma colour, cut..99th percentile, colour bar).
- New line stage `trace4.py` (anchored walk), concept in `docs/brief_anchored_walk.md`. Seed = strongest free dot;
  axis from PCA of a sphere of dots; band = dots connected (1.8 cubes) to the anchor inside a thin slice across the
  axis; node = weighted mean moved across only, pulled toward the trend (PULL); next anchor scored by closeness to the
  trend line (AFALL) within REACH; lines end exactly on other lines (collision), at uncertain spots, or when no anchor.
  Both directions from the seed; loop closure; short lines dropped.
- Settled defaults: SLICE 0.7, REACH 0.22, PULL 1, AFALL 2, ASTR 0, HILL 0, LEAD 0. Pauldron patch: 45 lines,
  6 nodes more than a cube off their local line (`out/t4_h0_r0.22_wide.png`).
- Tried and backed out: hill rule on the across-profile (red herring on rounded edges), CORE node, GROW clip,
  weak-dot filtering inside the band (LEAD): changed where nodes land but did not fix the rivet, and broke bands.
- App: "Method" switch in Build lines. Anchored walk (step, slice, reach, trend pull, min length) is the default;
  thick strokes + join stay selectable. `trace4.run(res, log, **kw)` is the entry point.
- Open: rivet (lines do not follow the geometry; no agreed rule yet), hooks at line ends, junction region C,
  backward smoothing (layer C), speed on full models (see timing below).
- Timing, full orc grid (`out/orcfull17_grid1200.npz`, 0.17 mm cubes, 249k dots): 111 s single core, 4394 lines,
  46k nodes, 163k dots (65%) unclaimed, most trials ending `no_anchor`/`short`. Usable, but the walk is a Python
  loop per node; the per-cluster multi-core split used by the strokes would cut it. The high unclaimed share at this
  coarse cube size is unjudged.
- Multi-core: `trace4.walk_all` splits the dots into clusters (connected within max(3 cubes, reach+1)) and walks each
  on its own core, as the strokes do; `run` uses it. Full orc: 111 s -> 32 s on 32 cores, bounded by the largest
  cluster (67k dots, 26 s alone). Going lower means cutting big clusters, which cuts lines (no join step now).
  Patch result nearly unchanged (51 lines vs 45; one messy blob near 192.5,122.8,161 reads "uncertain" on its own
  where the whole-model run borrowed a false axis from neighbouring dots). Render `out/t4_par_wide.png`.

## Session notes 2026-10-20: junction landing, bend hold, power cut; all backed out

Committed first: the 18/19 October work (trend pull 10, look-ahead probe LOOK, strength curve LEVP) as f16e15a.
End state: **the walk is back at f16e15a** (commit d9479a9); only the app's dot-set saving stays (`out/dots/`,
last set reloaded with the model, "Load saved dots…"). The walk changes below live in the history
(036921d, 0968ea8, aa6885a) and in the notes here; none is in the code.

**Found and agreed.** Lines stopped about a millimetre short of junctions: a line's slice claims every connected dot up
to 1 mm beside it, so the first line through a junction claimed the start of the crossing band, and the next line
ended at the first owned dot (measured on the pauldron patch: 21 collision ends, stopping dot 0.6-1.3 mm from the
owning line). Agreed rules: dot ownership must never end a line; a line ends when its new node grazes a node of
another line (within one cube) and that node becomes shared by both (the junction record; later the region stage
walks loops of pieces between shared nodes, turning the same way at each junction).

**Built and judged.** Graze rule alone: junction ends exact, but on the whole orc the lines looked sloppy: hooks at
line ends (the look-ahead referee hands the node to the dots at junctions, where the crossing crease has more dots
ahead than the own crease), side-by-side lines in weak clouds (second seeds that the ownership stop used to kill),
and stubs/curls on scatter. Then: BEND (dots may steer only within 45 deg of the heading; held = in-cone dots only)
removed hooks on the shield (dents 235 -> 7 in the test area) but its first version straightened curves into chords
(ring edges broke); line power = gathered strength per node with a cut (default 3) removed the stubs while keeping
soft folds. The user judged the whole-model result still worse than the morning's: stupid lines on no visible
crease, wobble, breaks; the notch rim (a strong trail) untraced because a neighbouring line's slice claimed its dots
from 0.6 mm away so no seed could start there. Two further attempts (band width/axis from strong dots only; the
LEVP say curve applied to band, width, axis, anchor, support and claims) fragmented lines and bred siblings on the
shield; reverted unrun.

**Principles stated by the user (to design from, not yet implemented):** a strong dot always earns a line, however
short; a weak row earns one only if dense along its trail; weak and sparse is not a line and must never be born,
not killed later. A line that reaches another line shares a node with it. Lines must trail the geometry faithfully:
a steep curve is fine, a sharp corner is two lines meeting. The curve of influence (LEVP) is applied today only to
the node average, not to the band, width, axis, anchor choice or claims; the user sees today's drift and wobble as
that inconsistency. Test spots: shield facet (dump 20261020_144432, view centre 162.13,137.7,107.92), ring
(20261020_163631, 194.72,123.34,169.49), notch rim (20261020_165259, rim point 219.51,167.43,165.39).
Helpers in the scratchpad of this session: siblings.py, kinks.py, follow.py (node-by-node record), rewalk_view.py,
strong_trail.py, held_check.py, bend_sweep.py, powerview.py.

## Session notes 2026-10-21: fine refinement of one line (scratch only, no project code changed)

Goal agreed: after the walk, refine each line at a much finer resolution; the eye sees the math and nothing else.
Spot: the notch rim line (dump 20261020_144432, line nearest 219.51,167.43,165.39), resampled to 0.02 mm.

What held up (scratch `refine_planes5.py`, this session's scratchpad):
- Per node, a slice across the line walked along the surface at 0.0087 mm steps; the facing read as the render blends
  it (vertex normals interpolated; corners must be merged first or the blend is a no-op).
- Each side's plane = weighted average of position and facing on that side, every weight a smooth function of
  distance (rises in from the node, fades out) and of distance along the line; the node = where the two planes cross.
  Repeated from the new point until it settles (3-6 passes). Back onto the surface by a push along the planes' mean
  facing, never nearest-point (a coin toss inside a corner).
- Result: a continuous line on the crease, 1.5 s per 600 nodes, no jumps between neighbouring points; the step at
  the wall dip shows as a soft S. Renders `rim23b_*`.

What failed, and the lesson: every version with a hard edge (hump cut-off, strongest-bend anchor, peak of the
blurred bend, plane patches of N steps, nearest-point drop) produced jumps of 0.05 mm between points 0.02 mm apart.
The surface is continuous; a jump is always the measurement's fault. No cut-offs, no argmax, no boxes, anywhere.

Open: the clean curve through the measurements (`fit_joints.py`: joints where the ironed path bends, one spline).
Fitted tight it keeps the triangle-scale ripples; loose, it loses the step. On this line the step and the ripples
have the same size and sharpness along the path; the eye tells them apart by the wall's dip above, which the line
does not carry. No rule for that yet. The user ended the session confused by the dots; next time start by showing
only the measured red line on the smooth-shaded surface (`view_smooth.py`), nothing else.
