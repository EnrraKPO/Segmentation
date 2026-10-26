import sys, os, pickle; root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); sys.path.insert(0, root)
import numpy as np, shapely
res = pickle.load(open(os.path.join(root, 'out', 'paint_dump.pkl'), 'rb')); emb = res['emb']; part = res['part']; con = res['con']
region = con['region']; pid0 = con['pid0']; fi_of = {int(f): i for i, f in enumerate(part['faces'])}
li = 1260; sl = emb['sl']; sf = emb['sf']; sp = emb['sp']; sq = emb['sq']; idx = np.where(sl == li)[0]
same = []
for i in idx:
    f = int(sf[i]); mid = shapely.Point(0.5*(sp[i]+sq[i])); fi = fi_of[f]
    regs = [region[pid0[fi]+k] for k, pc in enumerate(part['pieces'][fi]) if shapely.from_wkb(pc[2]).distance(mid) < 1e-7]
    if len(set(regs)) < 2: same.append((i, f, len(regs), regs))
print(len(idx), 'strokes;', len(same), 'with one region on both sides (or one side only)'); print(same[:12])
# look at one such face
if same:
    i, f, n, regs = same[0]; fi = fi_of[f]; print('face', f, 'pieces', len(part['pieces'][fi]), 'nb', part['nb'][fi], 'walls', part['walls'][fi])
    for k, pc in enumerate(part['pieces'][fi]): print('  piece', k, 'area', round(pc[0], 8), 'region', region[pid0[fi]+k], 'stretches', pc[3][:6])
    print('  strokes in face:', [(j, np.round(sp[j], 5).tolist(), np.round(sq[j], 5).tolist(), int(emb['s0'][j]), int(emb['s1'][j]), int(sl[j])) for j in np.where(sf == f)[0]])
    print('  corners', np.round(emb['C2'][f], 5).tolist())
