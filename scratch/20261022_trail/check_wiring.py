"""trace4.trail vs the scratch trail, and the app's path (trace4.run over clusters on several cores) with the trail on and off"""
import sys, os, numpy as np
sys.path.insert(0, r'D:\Stuff\PrintSuit\SegmentationTool'); os.chdir(r'D:\Stuff\PrintSuit\SegmentationTool'); sys.path.insert(0, 'scratch/20261022_smooth'); sys.path.insert(0, 'scratch/20261022_trail')
import smooth_test as T, trace4, trail as TR


def main():
    z = np.load(T.DOTS, allow_pickle=True); lo = z['lo']; size = float(z['size']); cut = float(z['cut']); idx, N, turn = z['idx'], z['N'], z['turn']; C = lo+(idx+0.5)*size
    for spot in T.SPOTS:
        focus = np.array(T.SPOTS[spot][0]); sel = (np.abs(C-focus) < T.BOX).all(1); X = (idx[sel]+0.5).astype(float); wgt = np.maximum(turn[sel]-cut, 0)+1e-3
        a, _ = TR.trail(X, wgt, N[sel], size, log=lambda *_: None); b, _ = trace4.trail(X, wgt, N[sel], size, log=lambda *_: None)
        same = len(a) == len(b) and all(x.shape == y.shape and np.allclose(x, y) for x, y in zip(a, b)); print(spot, 'trace4.trail == scratch trail:', same)
        res = dict(idx=idx[sel], line=np.ones(sel.sum(), bool), turn=turn[sel], N=N[sel], lo=lo, size=size, cut=cut)
        for tr in (1, 0):
            L = trace4.run(res, log=lambda *_: None, **{**T.walk_params(), 'trail': tr, 'curve': 1, 'trend': 6, 'talong': 0.7, 'tacross': 0.35, 'tshare': 0.5, 'tfade': 0.25})
            print(f'  run() trail={tr}: {len(L)} lines, {sum(len(l) for l in L)} nodes')


if __name__ == '__main__': main()
