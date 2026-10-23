"""stage 4 on a dump: counts, walk agreement, bridge error against the tube dots; saves out/bridge_<dump>.npz"""
import sys, os, time, json, numpy as np; from scipy.spatial import cKDTree
sys.path.insert(0, r'D:\Stuff\PrintSuit\SegmentationTool'); os.chdir(r'D:\Stuff\PrintSuit\SegmentationTool'); import bridge
def main():
    d = sys.argv[1]; kw = {k: (v if k == 'model' else float(v)) for k, v in (a.split('=') for a in sys.argv[2:])}; tag = os.environ.get('TAG', '')
    w = np.load(d+'/lines.npz'); W = [w['pts'][a:a+k] for a, k in zip(np.cumsum(np.r_[0, w['n'][:-1]]), w['n'])]
    r = np.load(d+'/refined.npz'); R = [r['pts'][a:a+k] for a, k in zip(np.cumsum(np.r_[0, r['n'][:-1]]), r['n'])]; ri = r['index']
    t0 = time.time(); res = bridge.run(R, ri, W, log=print, **kw); print(f'{time.time()-t0:.1f} s', res['report'])
    print('walk record vs B joints:', bridge.walk_agreement(res, W, 0.087))
    # ends: distance to the nearest other line, before and after
    def endstat(S):
        pts = np.vstack(S); t = cKDTree(pts); o = np.repeat(np.arange(len(S)), [len(c) for c in S]); dd = []
        for i, c in enumerate(S):
            for e in (c[0], c[-1]):
                d_, j = t.query(e, k=300); m = o[j] != i; dd.append(d_[m][0] if m.any() else np.inf)
        dd = np.array(dd); return f'ends {len(dd)}, within one cube of another line {(dd <= 0.087).mean()*100:.0f} %, median {np.median(dd):.3f} mm'
    print('before (refined pieces + unrefined walk lines):', endstat(R+[W[i] for i in range(len(W)) if i not in set(ri.tolist())]))
    print('after :', endstat(res['lines']))
    # bridge nodes against the tube dots: distance from each bridge node to the nearest tube dot (the measurement)
    B = np.vstack([l[~m] for l, m in zip(res['lines'], res['measured']) if (~m).any()]); D = r['dots']
    sel = (np.abs(D[:, None, :]-B[None, :5, :]) < 0).any(1) if False else None
    dt = cKDTree(D[::3]); dd = dt.query(B)[0]; print(f'bridge nodes {len(B)}: distance to the nearest tube dot p50 {np.percentile(dd,50):.3f} p90 {np.percentile(dd,90):.3f} max {dd.max():.2f} mm')
    np.savez(f'out/bridge_{os.path.basename(d)}{tag}.npz', n=np.array([len(l) for l in res['lines']]), pts=np.vstack(res['lines']), measured=np.concatenate(res['measured']), walk=np.array(res['walk']), shared=np.array([s[2] for s in res['shared']]), loose=np.array([l[2] for l in res['loose']]))


if __name__ == '__main__': main()
