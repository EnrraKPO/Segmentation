"""export the orc painting (out/paint_walk.pkl) as a Bambu Studio project and check the notes decode"""
import sys, os, time, pickle, re, zipfile, collections; root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); sys.path.insert(0, root)
import numpy as np
import bambu3mf, regions

if __name__ == '__main__':
    res = pickle.load(open(os.path.join(root, 'out', 'paint_walk.pkl'), 'rb')); emb = res['emb']
    fine = float(sys.argv[1]) if len(sys.argv) > 1 else 0.05; depth = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    notes = bambu3mf.paint_notes(res, res['colour']+1, fine=fine, depth=depth, log=print)
    out = os.path.join(root, 'out', 'orc_painted.3mf'); bambu3mf.write_project(out, emb['P'], emb['F'], notes, regions.PALETTE, name='orc', log=print)
    L = np.array([len(s) for s in notes]); print('note length: mean', L.mean().round(2), 'p99', np.percentile(L, 99), 'max', L.max(), 'total nibbles', L.sum())
    sys.path.insert(0, os.path.dirname(__file__)); import importlib; bc = importlib.import_module('bambu_check')
    states = collections.Counter(); leaves = 0
    for s in notes[:200000]:
        t = bc.decode(bambu3mf.from_string(s))
        def walk(t):
            global leaves
            if isinstance(t, int): states[t] += 1
            else: [walk(c) for c in t[2]]
        walk(t)
    print('decoded 200k notes; leaf states', dict(sorted(states.items())))
    with zipfile.ZipFile(out) as z: print('package:', z.namelist())
