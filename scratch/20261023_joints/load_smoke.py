"""offscreen smoke test of Load lines…: dump folder -> App via read_lines, a fake joined set, save with stages, load the save back, dump without tubes"""
import os, sys, time, shutil, numpy as np
os.environ['QT_QPA_PLATFORM'] = 'offscreen'; sys.path.insert(0, r'D:\Stuff\PrintSuit\SegmentationTool'); os.chdir(r'D:\Stuff\PrintSuit\SegmentationTool')
from PySide6 import QtWidgets; import app as A
def main():
    d = 'out/debug/20261023_114229'; before = set(os.listdir('out/log')); dbg = set(os.listdir('out/debug'))
    qa = QtWidgets.QApplication(sys.argv); w = A.App()
    w.P, w.F = A.load_mesh('out/orcfull_patch.stl'); w.path = 'out/orcfull_patch_smoke.stl'; w.busy(False); print('load button enabled:', w.btn_load_lines.isEnabled())
    st = A.read_lines(d+'/lines.npz'); w.curves = st['walk']; w.walk_shared = st['walk_shared']; w.refined = st['refined']; w.joined = st['joined']
    # a stand-in joined set: the pieces themselves, all measured, so the round trip of that stage is exercised without a 10-minute join
    L = w.refined['lines'][:50]; w.joined = dict(lines=L, measured=[np.ones(len(c), bool) for c in L], walk=w.refined['index'][:50], shared=[(0, 1, L[0][0])], loose=[(2, 'end', L[1][-1])], records=[], report={}, agreement=None)
    w.chk_snap.setChecked(False); w.chk_dots.setChecked(False); t0 = time.time(); w.redraw(); print(f'redraw {time.time()-t0:.1f} s, actors {len(w.plotter.renderer.actors)}')
    w.save(); print(w.txt.toPlainText().splitlines()[-1])
    p = 'out/orcfull_patch_smoke_lines.npz'; z = A._npz(p); print('keys:', list(z))
    back = A.read_lines(p); print('back:', back['note'], len(back['walk']), len(back['refined']['lines']), len(set(back['refined']['index'])), len(back['joined']['lines']), back['joined']['measured'][0].dtype, len(back['joined']['shared']), len(back['joined']['loose']))
    assert all(np.array_equal(a, b) for a, b in zip(back['walk'], w.curves)) and all(np.array_equal(a, b) for a, b in zip(back['refined']['lines'], w.refined['lines'])) and back['refined']['index'] == [int(i) for i in w.refined['index']]
    w.curves, w.walk_shared, w.refined, w.joined = back['walk'], back['walk_shared'], back['refined'], back['joined']; w.busy(False); w.redraw(); print('redraw after load ok; rmerge enabled:', w.btn_rmerge.isEnabled(), 'join enabled:', w.btn_join.isEnabled())
    w.dump(); print(w.txt.toPlainText().splitlines()[-1][:160])
    os.remove(p)
    for f in set(os.listdir('out/log'))-before: os.remove(os.path.join('out/log', f))
    for f in set(os.listdir('out/debug'))-dbg: shutil.rmtree(os.path.join('out/debug', f))
    print('cleaned')
if __name__ == '__main__': main()
