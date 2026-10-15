"""Strokes through the whole turn band (every cube with turn >= cut, weighted by turn-cut) instead of crest dots only.
usage: python bandstrokes.py grid.npz out.npz cut W [L K]"""
import sys, numpy as np
sys.path.insert(0, 'D:/Stuff/PrintSuit/SegmentationTool'); import creases
g = np.load(sys.argv[1]); outp = sys.argv[2]; cut = float(sys.argv[3]); W = float(sys.argv[4]); L = float(sys.argv[5]) if len(sys.argv) > 5 else 1.0; K = float(sys.argv[6]) if len(sys.argv) > 6 else 0.26
res = {k: g[k] for k in g.files}; res['line'] = res['turn'] >= cut; res['turn'] = np.maximum(res['turn']-cut, 0)+1e-3
cv = creases.strokes(res, W=W, L=L, K=K, log=print)
np.savez(outp, n=np.array([len(c) for c in cv]), pts=np.vstack(cv) if cv else np.zeros((0, 3))); np.savez(outp.replace('.npz', '_grid.npz'), **res); print('saved', outp)
