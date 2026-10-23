"""pieces per tube under different fine-walk / rung settings, no code change: usage python handles.py dumpdir line [line...]"""
import sys, os, ast, time, numpy as np
sys.path.insert(0, r'D:\Stuff\PrintSuit\SegmentationTool'); os.chdir(r'D:\Stuff\PrintSuit\SegmentationTool'); import refine, trace4, trimesh
def main():
    d = sys.argv[1]; ids = [int(a) for a in sys.argv[2:]]
    w = np.load(d+'/lines.npz'); W = [w['pts'][a:a+k] for a, k in zip(np.cumsum(np.r_[0, w['n'][:-1]]), w['n'])]
    line12 = [l for l in open(d+'/log.txt', encoding='utf-8') if l.startswith('--- refine all')][0]
    walk = ast.literal_eval(line12.split('fine walk ')[1].split('; tune')[0]); tune = ast.literal_eval(line12.split('tune ')[-1].strip())
    m = trimesh.load('out/orcfull_patch.stl', process=True); P, F = np.asarray(m.vertices), np.asarray(m.faces); h = float((P.max(0)-P.min(0)).max()); cube = 0.087
    tree, owner = refine.nodes_tree(W)
    runs = [('as run (reach 0.087, rungs every 0.5 cube, cone 45)', {}, {}), ('reach 0.174', {'reach': 0.174}, {}), ('rungs every 0.25 cube', {}, {'rungstep': 0.25}), ('cone 70', {'cone': 70}, {}), ('reach 0.174 + rungs 0.25', {'reach': 0.174}, {'rungstep': 0.25})]
    for i in ids:
        L = W[i]; others = [W[j] for j in refine.end_lines(W, i, tune['slack']*cube, tree, owner)]; print(f'=== line {i} ({len(L)} walk nodes)')
        for name, wk, tk in runs:
            t0 = time.time(); t = refine.tube(L, P, F, cube, h, 3, 0.3, 1.0, True, None, 1, 0.087, 0.087, tune={**tune, **tk}); t = refine.avoid(t, others, 2*0.087, None, res=tune['avoidres'])
            line = t['line']; size = float(t['size']); X = (t['idx'][line]+0.5).astype(float); wgt = np.maximum(t['turn'][line]-1.0, 0)+1e-3; N = t['N'][line]
            lines, ow = trace4.walk(X, wgt, N, size, log=lambda *_: None, **{**walk, **wk}); ends = {}
            for _, reasons, *_r in trace4.walk.diag:
                for r in reasons: ends[r] = ends.get(r, 0)+1
            # how far the fine nodes sit from the walk line (the measurement's accuracy is not to change)
            from scipy.spatial import cKDTree; dv = cKDTree(L).query(t['lo']+np.vstack(lines)*size)[0] if lines else np.zeros(1)
            print(f'  {name:48} pieces {len(lines):2}, nodes {[len(l) for l in lines]}, ends {ends}, node-to-walk-line p50 {np.percentile(dv,50):.3f} p90 {np.percentile(dv,90):.3f} mm, {time.time()-t0:.0f} s')
if __name__ == '__main__': main()
