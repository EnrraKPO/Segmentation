# Brief: the anchored walk (2026-10-17)

Replacement for the line stage (`creases.strokes` + `creases.join`), agreed in conversation. Dots stay as they are:
one point per cube where the surface folds sharply, with a strength (degrees of fold).

## What was wrong with the current stroke

- The cursor is a point in free space, one step ahead along the previous heading. The window is laid along that
  heading, so the previous heading chooses which dots are seen, and those dots give the next heading. The old
  direction leaks into the new one. Lines go straight where the band bends, and overrun past the end of the dots.
- The window width is fixed. A plate edge makes a band several cubes wide, a wrinkle one cube. One width either
  swallows the neighbouring crease or misses the flanks.
- Lines do not stop at other lines; junctions are patched afterwards by extending ends (the join step).

## The new process

Every line is a chain of nodes. A node is the weighted average of the dots in one slice across the band.

1. **Seed.** The strongest dot not yet claimed by any line.
2. **Axis at a dot.** Take the dots in a sphere around it, no preferred direction. Fit a straight line through them
   (axis of greatest spread). If the dots are not clearly longer than wide along that axis, the reading is
   uncertain: drop the trial. The sphere is a few times wider than the band, so a wide band still reads as a line.
3. **Band across the axis.** Take a thin slice: a couple of cubes long along the axis, as wide as needed across.
   The band is the set of slice dots connected to the anchor: dots are joined when they touch (within about 1.8
   cubes), so a neighbouring band is told apart by the empty gap between them, not by a width limit. (A 'hill' rule
   on the across-profile was tried and dropped: on rounded edges the profile has no clear hill.)
4. **Node.** The weighted average of the dots in the band slice, moved only across the band (not along it), and
   pulled toward the line's trend (layer A, handle PULL). Those dots are claimed by this line.
5. **Next anchor.** From the newest node, extrapolate the line's trend (the last few nodes, not just the last
   heading). Among the dots roughly one step ahead along the trend and within the band's width, take the strongest.
   It becomes the anchor: steps 2 to 4 run on it, with no memory of the previous axis except to decide which end of
   the new axis points onward. The line grows both ways from the seed; after the first step one way leads into
   claimed dots and dies, so no notion of forward is needed.
6. **Ends.** No dot ahead: the band has ended. The anchor or its band belongs to another line (or to this line's
   own earlier part): collision. The line is ended exactly on that line, at its nearest node, so the junction is
   made at trace time. Where the band confuses the reading (uncertain axis, two hills with the anchor between),
   the line stops; the free dots beyond will seed a later line that grows back to the same spot. Crossings thus
   become two lines meeting at a junction, which is what regions need.
7. Lines shorter than a minimum are dropped and their dots released.

## Settled state (2026-10-17)

Defaults in `trace4.py`: SLICE 0.7, REACH 0.22, PULL 1, AFALL 2, ASTR 0, HILL 0, LEAD 0. On the pauldron patch:
45 lines, 6 nodes more than a cube off their local line, 471 dots unclaimed. Render `out/t4_h0_r0.22_wide.png`.

Not solved: the rivet. Its lines do not follow the geometry. Weak-dot filtering inside the band (LEAD) was tried
and backed out.
Also open: hooks at line ends, junction region C, backward smoothing (layer C).

## Open questions to test

- The band rule on a bevelled edge: two creases whose bands touch read as one hill. One line down the middle, or
  two? The eye decides.
- How many nodes the trend should use; how wide the search cone for the next anchor.
- Whether ending at every uncertain spot leaves too many short pieces.

## Test

`trace4.py grid.npz out.npz` on the pauldron groove patch (`out/pvi4_grid.npz`, 0.087 mm cubes) and the toes patch,
rendered with `closeup2.py` against the current `out/pvi4.npz`. Judge: one line per crease, through corners, no
overrun, junctions met.
