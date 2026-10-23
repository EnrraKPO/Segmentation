"""reproduce the app's refine of one line with the trail / curve toggles, single core"""
import sys, os, traceback, numpy as np, trimesh
os.environ['QT_QPA_PLATFORM'] = 'offscreen'; sys.path.insert(0, r'D:\Stuff\PrintSuit\SegmentationTool'); os.chdir(r'D:\Stuff\PrintSuit\SegmentationTool'); sys.path.insert(0, 'scratch/20261022_smooth')
sys.argv = ['app.py']; import app, trace4, refine, smooth_test as T
from PySide6 import QtWidgets


def main():
    q = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv); w = app.App()
    z = np.load(T.DOTS, allow_pickle=True); lo = z['lo']; size = float(z['size']); cut = float(z['cut']); idx, N, turn = z['idx'], z['N'], z['turn']; C = lo+(idx+0.5)*size
    focus = np.array(T.SPOTS['ring'][0]); sel = (np.abs(C-focus) < T.BOX).all(1)
    res = dict(idx=idx[sel], line=np.ones(sel.sum(), bool), turn=turn[sel], N=N[sel], lo=lo, size=size, cut=cut)
    m = trimesh.load('out/orcfull_patch.stl', process=False); P = np.asarray(m.vertices, float); F = np.asarray(m.faces); h = float((P.max(0)-P.min(0)).max())
    for tr, cv in ((0, 0), (0, 1), (1, 1), (1, 0)):
        w.chk_trail.setChecked(bool(tr)); w.chk_curve.setChecked(bool(cv)); wp = w.walk_params(); rp = w.refine_params()
        curves = trace4.run(res, log=lambda *_: None, **wp); L = max(curves, key=len)
        for i in range(len(curves)):                                            # every line on its own, single core, to find the one that crashes
            try: r = refine.refine_all(curves, P, F, size, h, log=lambda *_: None, only=i, flat=True, **rp)
            except Exception: print(f'trail={tr} curve={cv}: line {i} ({len(curves[i])} nodes) CRASH'); traceback.print_exc(); break
        else: print(f'trail={tr} curve={cv}: all {len(curves)} lines refined one by one')
        if os.environ.get('ALL'):
            try: r = refine.refine_all(curves, P, F, size, h, log=print, only=None, flat=True, **rp); print(f'  all at once on the cores: {len(r["lines"])} fine lines')
            except Exception: print('  all at once: CRASH'); traceback.print_exc()


if __name__ == '__main__': main()
