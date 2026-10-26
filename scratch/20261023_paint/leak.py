"""one ring alone on the orc: where do its two sides connect?"""
import sys, os, pickle, time; root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); sys.path.insert(0, root)
import numpy as np, shapely, regions, linesio
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import shortest_path
P, F = linesio.load_mesh(os.path.join(root, 'out', 'orcfull_patch.stl')); st = linesio.read_lines(os.path.join(root, 'out', 'debug', '20261023_114229', 'lines.npz'))
li = int(sys.argv[1]) if len(sys.argv) > 1 else 1260
L = st['walk'][li]; print('line', li, len(L), 'nodes, closed:', np.linalg.norm(L[0]-L[-1]) < 1e-9)
emb = regions.embed(P, F, [L], log=print); part = regions.partition(emb, log=print, workers=1); con = regions.connect(emb, part, log=print)
print('regions', con['nreg'], 'areas', np.round(np.sort(con['areas'])[::-1][:5], 3))
if con['nreg'] >= 2: sys.exit()
# rebuild the connection graph as connect() does and find a path between the two sides of the first stroke
region = con['region']; pid0 = con['pid0']; fi_of = {int(f): i for i, f in enumerate(part['faces'])}; m = len(F)
f0 = int(emb['sf'][0]); fi = fi_of[f0]; mid = shapely.Point(0.5*(emb['sp'][0]+emb['sq'][0]))
sides = [int(pid0[fi]+k) for k, pc in enumerate(part['pieces'][fi]) if shapely.from_wkb(pc[2]).distance(mid) < 1e-7]; print('first stroke in face', f0, 'pieces on its sides', sides)
# graph
crossed = con['crossed']; owner = {}; walls = set()
for fi_, f in enumerate(part['faces']):
    for pi, pc in enumerate(part['pieces'][fi_]):
        for stt in pc[3]: owner.setdefault(stt, []).append(int(pid0[fi_]+pi))
    walls.update(part['walls'][fi_])
ne = {e: len(v) for e, v in part['ep_by_edge'].items()}
fe = emb['fedge'].ravel(); nb = emb['nbr'].ravel(); fa = np.repeat(np.arange(m), 3); keep = nb > fa; ea = fe[keep]; A = fa[keep]; B = nb[keep]
both = ~crossed[A] & ~crossed[B]; I = list(A[both]); J = list(B[both]); why = {}
for e, a, b in zip(ea[~both], A[~both], B[~both]):
    e = int(e)
    for j in range(ne.get(e, 0)+1):
        if (e, j) in walls: continue
        side = list(owner.get((e, j), []))
        for h in (a, b):
            if not crossed[h]: side.append(int(h))
        for x in side[1:]: I.append(side[0]); J.append(x); why[(side[0], x)] = (e, j, len(side))
n = int(pid0[-1]); g = coo_matrix((np.ones(len(I)), (I, J)), shape=(n, n)).tocsr()
dist, pred = shortest_path(g, directed=False, unweighted=True, indices=sides[0], return_predecessors=True)
path = []; x = sides[1]
while x != sides[0] and x >= 0: path.append(x); x = pred[x]
path.append(sides[0]); path = path[::-1]; print('path length', len(path))
face_of = con['face_of']
for a, b in zip(path[:-1], path[1:]):
    if a >= m or b >= m:
        k = why.get((a, b)) or why.get((b, a)); print('  step', a, '->', b, 'faces', face_of[a], face_of[b], 'via stretch', k)
