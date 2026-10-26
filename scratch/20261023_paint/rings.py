"""closed walk lines (rings): do they separate two regions? For each ring, the pieces on both sides of its strokes."""
import sys, os, pickle; root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); sys.path.insert(0, root)
import numpy as np, shapely
res = pickle.load(open(os.path.join(root, 'out', 'paint_walk.pkl'), 'rb')); emb = res['emb']; part = res['part']; con = res['con']
Q = emb['nodes']; start = emb['start']; region = con['region']; pid0 = con['pid0']; m = len(emb['F'])
fi_of = {int(f): i for i, f in enumerate(part['faces'])}
polys = {}
def pieces(f):
    if f not in polys: polys[f] = [(shapely.from_wkb(pc[2]), int(pid0[fi_of[f]]+k)) for k, pc in enumerate(part['pieces'][fi_of[f]])]
    return polys[f]
fails = set(int(i) for i in emb['fid'])
rings = [li for li in range(len(start)-1) if np.linalg.norm(Q[start[li]]-Q[start[li+1]-1]) < 1e-9]
print(len(rings), 'closed rings of', len(start)-1, 'lines')
sl = emb['sl']; sf = emb['sf']; sp = emb['sp']; sq = emb['sq']
sep = 0; nosep = 0; nosep_fail = 0; nosep_ok = []
for li in rings:
    has_fail = any(i in fails for i in range(start[li], start[li+1]-1))
    idx = np.where(sl == li)[0]; sides = set()
    for i in idx[::max(1, len(idx)//40)]:
        f = int(sf[i]); mid = shapely.Point(0.5*(sp[i]+sq[i]))
        regs = {region[p] for poly, p in pieces(f) if poly.distance(mid) < 1e-7}
        sides |= regs
        if len(regs) == 2: break
    if len(regs) == 2: sep += 1
    else:
        nosep += 1
        if has_fail: nosep_fail += 1
        else: nosep_ok.append((li, start[li+1]-start[li], len(idx)))
print(f'separate two regions: {sep}; do not: {nosep} (of which with an undrawn stretch: {nosep_fail})')
print('rings without any undrawn stretch that still do not separate (line, nodes, strokes):', nosep_ok[:15])
if nosep_ok:
    li = nosep_ok[0][0]; print('nodes of line', li, Q[start[li]:start[li]+3], '...', Q[start[li+1]-3:start[li+1]])
