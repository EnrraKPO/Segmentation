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
