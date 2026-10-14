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
