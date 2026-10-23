"""stage 3b on the dump's refined pieces, then the join stage on the merged pieces: counts, and the doubles check after"""
import sys, os, time, numpy as np; from scipy.spatial import cKDTree
sys.path.insert(0, r'D:\Stuff\PrintSuit\SegmentationTool'); os.chdir(r'D:\Stuff\PrintSuit\SegmentationTool'); import merge, bridge
def doubles(R, ri):
    by = {}
    for k, i in enumerate(ri): by.setdefault(int(i), []).append(k)
    n = 0; pairs = 0
    for i, ks in by.items():
        for a in range(len(ks)):
            for b in range(a+1, len(ks)):
                pairs += 1; dd = cKDTree(R[ks[a]]).query(R[ks[b]])[0]
                if (dd < 0.3).sum() >= 3: n += 1
    return n, pairs
def main():
    d = 'out/debug/20261023_114229'; margin = float(sys.argv[1]) if len(sys.argv) > 1 else 0.05
    w = np.load(d+'/lines.npz'); W = [w['pts'][a:a+k] for a, k in zip(np.cumsum(np.r_[0, w['n'][:-1]]), w['n'])]
    r = np.load(d+'/refined.npz'); R = [r['pts'][a:a+k] for a, k in zip(np.cumsum(np.r_[0, r['n'][:-1]]), r['n'])]; ri = list(r['index'])
    print('before: overlapping piece pairs within one walk line', doubles(R, ri))
    t0 = time.time(); res = merge.run(R, margin, log=print, labels=ri); print(f'{time.time()-t0:.1f} s', res['report'])
    M, lab = res['lines'], res['labels']; print('after : overlapping piece pairs within one walk line', doubles(M, lab))
    np.savez(f'out/pieces_merged_{os.path.basename(d)}.npz', n=np.array([len(l) for l in M]), pts=np.vstack(M), index=np.array(lab))
    t0 = time.time(); j = bridge.run(M, lab, W, log=print); print(f'join on the merged pieces: {time.time()-t0:.0f} s', j['report'])
if __name__ == '__main__': main()
