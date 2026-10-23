"""end reasons of the fine walk on one line's tube: runs the tube and the walk directly (one core) and prints trace4.walk's stats"""
import sys, os, ast, numpy as np
sys.path.insert(0, r'D:\Stuff\PrintSuit\SegmentationTool'); os.chdir(r'D:\Stuff\PrintSuit\SegmentationTool'); import refine, trace4, trimesh
def main():
    d = sys.argv[1]; ids = [int(a) for a in sys.argv[2:]]
    w = np.load(d+'/lines.npz'); W = [w['pts'][a:a+k] for a, k in zip(np.cumsum(np.r_[0, w['n'][:-1]]), w['n'])]
    line12 = [l for l in open(d+'/log.txt', encoding='utf-8') if l.startswith('--- refine all')][0]
    walk = ast.literal_eval(line12.split('fine walk ')[1].split('; tune')[0]); tune = ast.literal_eval(line12.split('tune ')[-1].strip())
    m = trimesh.load('out/orcfull_patch.stl', process=True); P, F = np.asarray(m.vertices), np.asarray(m.faces); h = float((P.max(0)-P.min(0)).max()); cube = 0.087
    tree, owner = refine.nodes_tree(W)
    for i in ids:
        L = W[i]; others = [W[j] for j in refine.end_lines(W, i, tune['slack']*cube, tree, owner)]
        t = refine.tube(L, P, F, cube, h, 3, 0.3, 1.0, True, None, 1, 0.087, 0.087, tune=tune); t = refine.avoid(t, others, 2*0.087, None, res=tune['avoidres'])
        line = t['line']; size = float(t['size']); X = (t['idx'][line]+0.5).astype(float); wgt = np.maximum(t['turn'][line]-1.0, 0)+1e-3; N = t['N'][line]
        print(f'=== line {i}'); lines, ow = trace4.walk(X, wgt, N, size, log=print, **walk)
        for k, (tid, reasons, ends_at, _) in enumerate(trace4.walk.diag): print(f'   piece {tid}: {len(lines[tid])} nodes, ends: {reasons}')
if __name__ == '__main__': main()
