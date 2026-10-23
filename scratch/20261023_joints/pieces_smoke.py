"""offscreen smoke test of app group 3b on a subset of the dump's refined pieces"""
import os, sys, time, numpy as np
os.environ['QT_QPA_PLATFORM'] = 'offscreen'; sys.path.insert(0, r'D:\Stuff\PrintSuit\SegmentationTool'); os.chdir(r'D:\Stuff\PrintSuit\SegmentationTool')
from PySide6 import QtWidgets; import app as A, merge
def main():
    qa = QtWidgets.QApplication(sys.argv); w = A.App(); w.P, w.F = A.load_mesh('out/orcfull_patch.stl'); w.path = 'out/orcfull_patch_smoke.stl'
    d = 'out/debug/20261023_114229'; z = np.load(d+'/lines.npz'); w.curves = [z['pts'][a:a+k] for a, k in zip(np.cumsum(np.r_[0, z['n'][:-1]]), z['n'])]
    r = np.load(d+'/refined.npz'); R = [r['pts'][a:a+k] for a, k in zip(np.cumsum(np.r_[0, r['n'][:-1]]), r['n'])]; ri = list(r['index'])
    c = np.array([173.08, 143.54, 96.98]); keep = [k for k, L in enumerate(R) if (np.abs(L-c) < 6).all(1).any()]
    w.refined = dict(lines=[R[k] for k in keep], index=[ri[k] for k in keep], tubes=[]); w.chk_snap.setChecked(False); w.chk_dots.setChecked(False); w.busy(False)
    print('pieces in the box', len(keep), '| merge pieces button enabled:', w.btn_rmerge.isEnabled())
    t0 = time.time(); res = merge.run(w.refined['lines'], w.sp_rmargin.value(), log=w.logmsg, labels=w.refined['index']); w._pieces_merged(res, t0); print(w.txt.toPlainText().splitlines()[-1])
    print('index carried:', len(w.refined['index']) == len(w.refined['lines']), '| walk lines referenced:', len(set(w.refined['index'])))
    w.redraw(); print('actors', len(w.plotter.renderer.actors)); w.save(); print(w.txt.toPlainText().splitlines()[-1])
if __name__ == '__main__': main()
