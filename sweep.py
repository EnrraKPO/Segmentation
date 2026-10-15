"""Sweep one stroke/cut handle at a time on the toe patch grid (0.043 mm) and save strokes per value."""
import sys, numpy as np
sys.path.insert(0, 'D:/Stuff/PrintSuit/SegmentationTool'); import creases
g = np.load('D:/Stuff/PrintSuit/SegmentationTool/out/toes_0.043_grid.npz'); base = {k: g[k] for k in g.files}
sweeps = dict(cut=[8, 12, 17, 25, 35], W=[0.15, 0.25, 0.35, 0.6, 1.0], L=[0.4, 0.7, 1.0, 2.0, 4.0], K=[0.1, 0.26, 0.5, 1.0], align=[1.5, 2, 3, 5], agree=[15, 30, 60])
for name, vals in sweeps.items():
    for v in vals:
        res = dict(base)
        if name == 'cut':
            change, crest = creases.fold_change(res['idx'], res['N'], float(v)); res['turn'] = change; res['line'] = crest; kw = {}
        else: kw = {name: v}
        cv = creases.strokes(res, **kw); n = np.array([len(c) for c in cv])
        print(f'{name}={v}: dots {int(res["line"].sum())}, strokes {len(cv)}, nodes {int(n.sum())}, longest {int(n.max()) if len(n) else 0}', flush=True)
        np.savez(f'D:/Stuff/PrintSuit/SegmentationTool/out/sw_{name}_{v}.npz', n=n, pts=np.vstack(cv) if len(cv) else np.zeros((0, 3)))
        np.savez(f'D:/Stuff/PrintSuit/SegmentationTool/out/sw_{name}_{v}_grid.npz', **res)
