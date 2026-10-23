"""offscreen smoke test of app group 2b on a subset of the dump's walk lines"""
import os, sys, time, numpy as np
os.environ['QT_QPA_PLATFORM'] = 'offscreen'; sys.path.insert(0, r'D:\Stuff\PrintSuit\SegmentationTool'); os.chdir(r'D:\Stuff\PrintSuit\SegmentationTool')
from PySide6 import QtWidgets; import app as A, merge
def main():
    qa = QtWidgets.QApplication(sys.argv); w = A.App(); w.P, w.F = A.load_mesh('out/orcfull_patch.stl'); w.path = 'out/orcfull_patch_smoke.stl'
    z = np.load('out/debug/20261023_114229/lines.npz'); W = [z['pts'][a:a+k] for a, k in zip(np.cumsum(np.r_[0, z['n'][:-1]]), z['n'])]
    c = np.array([173.08, 143.54, 96.98]); w.curves = [L for L in W if (np.abs(L-c) < 6).all(1).any()]; print('lines in the box', len(w.curves))
    w.chk_snap.setChecked(False); w.chk_dots.setChecked(False); w.busy(False); print('merge button enabled:', w.btn_merge.isEnabled())
    t0 = time.time(); r = merge.run(w.curves, w.sp_margin.value(), log=w.logmsg); w._merge_done(r, t0); print(w.txt.toPlainText().splitlines()[-1])
    w.redraw(); print('actors', len(w.plotter.renderer.actors), 'shared drawn', len(w.walk_shared)); print('all_params merge_margin:', w.all_params()['merge_margin'])
if __name__ == '__main__': main()
