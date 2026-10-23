"""today's walk vs the trail on the three spots: numbers, out/trail_<spot>.npz, and (in a separate process per spot, the
renderer crashes when run after the walk) out/trail_<spot>_compare.png: left today's walk, right the trail.
env: ALONG, ACROSS, TSHARE, TREACH, FADE override the trail's handles (mm / share / nodes)."""
import sys, os, time, subprocess, numpy as np
sys.path.insert(0, r'D:\Stuff\PrintSuit\SegmentationTool'); os.chdir(r'D:\Stuff\PrintSuit\SegmentationTool'); sys.path.insert(0, 'scratch/20261022_smooth'); sys.path.insert(0, 'scratch/20261022_trail')
import smooth_test as T, trace4, trail as TR
H = {k: float(os.environ[k]) for k in ('ALONG', 'ACROSS', 'TSHARE', 'TREACH', 'FADE') if k in os.environ}; hk = {k.lower(): v for k, v in H.items()}
tag = os.environ.get('TAG', '')


def main():
    z = np.load(T.DOTS, allow_pickle=True); lo = z['lo']; size = float(z['size']); cut = float(z['cut']); idx, N, turn = z['idx'], z['N'], z['turn']; C = lo+(idx+0.5)*size
    print(f'{"spot":7} {"run":6} {"lines":>5} {"nodes":>6} {"turn p50":>8} {"p90":>6} {"max":>6} {"off-band nodes":>14} {"farthest":>8}')
    for spot, (focus, half) in T.SPOTS.items():
        focus = np.array(focus); sel = (np.abs(C-focus) < T.BOX).all(1); X = (idx[sel]+0.5).astype(float); wgt = np.maximum(turn[sel]-cut, 0)+1e-3; Nn = N[sel]
        res = {}
        for name in ('today', 'trail'):
            t0 = time.time()
            if name == 'today': lines, _ = trace4.walk(X, wgt, Nn, size, log=lambda *_: None, **{**T.walk_params(), 'curve': 0, 'trend': 6})
            else: lines, _ = TR.trail(X, wgt, Nn, size, log=print, **hk)
            L = [lo+l*size for l in lines]; res[name] = L; s = T.stats(L, C[sel])
            print(f'{spot:7} {name:6} {s["lines"]:5} {s["nodes"]:6} {s["t50"]:8.1f} {s["t90"]:6.1f} {s["tmax"]:6.1f} {s["offn"]:14} {s["offmax"]:8.3f}  ({time.time()-t0:.1f} s)')
        np.savez(f'out/trail_{spot}{tag}.npz', **{f'{k}_n': [len(l) for l in v] for k, v in res.items()}, **{f'{k}_pts': np.vstack(v) for k, v in res.items()})
        lab = 'trail' + (': '+', '.join(f'{k.lower()} {v:g}' for k, v in H.items()) if H else ' (defaults)')
        subprocess.run([sys.executable, 'scratch/20261022_trail/render2.py', spot, tag, "today's walk", lab])


if __name__ == '__main__': main()
