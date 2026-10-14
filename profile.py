"""Cross-section of the grid's turn across a point: walk +-n cells along d1 (the across-the-line direction) from the
surface cell nearest POINT, print per cell: turn (deg/cell), angle of the normal vs the centre normal, line flag.
usage: python profile.py grid.npz x,y,z [n]"""
import sys, numpy as np
from scipy.spatial import cKDTree
g = np.load(sys.argv[1]); p = np.array([float(v) for v in sys.argv[2].split(',')]); n = int(sys.argv[3]) if len(sys.argv) > 3 else 8
X = g['lo']+(g['idx']+0.5)*g['size']; tree = cKDTree(X); c = tree.query(p)[1]
d1 = g['d1'][c]; N0 = g['N'][c]; size = g['size']
print(f'centre cell {X[c].round(2)} turn {g["turn"][c]:.1f} line {bool(g["line"][c])}  (cell {size:.3f} mm)')
for k in range(-n, n+1):
    q = X[c]+k*size*d1; j = tree.query(q)[1]
    ang = np.degrees(np.arccos(np.clip(g['N'][j] @ N0, -1, 1)))
    print(f'{k:+3d}  turn {g["turn"][j]:5.1f}  normal vs centre {ang:5.1f} deg  {"LINE" if g["line"][j] else ""}')
