"""offscreen smoke test of app group 4: load the dump into an App, run join_lines' job synchronously, redraw, save, dump"""
import os, sys, time, numpy as np
os.environ['QT_QPA_PLATFORM'] = 'offscreen'; sys.path.insert(0, r'D:\Stuff\PrintSuit\SegmentationTool'); os.chdir(r'D:\Stuff\PrintSuit\SegmentationTool')
import pyvista as pv
from PySide6 import QtWidgets; import app as A
def main():
    d = 'out/debug/20261023_114229'
    qa = QtWidgets.QApplication(sys.argv); w = A.App(); print('session log line:', w.txt.toPlainText().splitlines()[-1])
    w.P, w.F = A.load_mesh('out/orcfull_patch.stl'); w.path = 'out/orcfull_patch_smoke.stl'
    z = np.load(d+'/lines.npz'); w.curves = [z['pts'][a:a+k] for a, k in zip(np.cumsum(np.r_[0, z['n'][:-1]]), z['n'])]
    r = np.load(d+'/refined.npz'); w.refined = dict(lines=[r['pts'][a:a+k] for a, k in zip(np.cumsum(np.r_[0, r['n'][:-1]]), r['n'])], index=list(r['index']), tubes=[dict(dots=r['dots'][:10], turn=r['turn'][:10], size=float(r['size']))])
    w.chk_snap.setChecked(False); w.chk_dots.setChecked(False); w.busy(False); print('join button enabled:', w.btn_join.isEnabled())
    # the button's job, run here instead of in the worker thread
    jp = w.join_params(); t0 = time.time()
    import bridge
    res = bridge.run(w.refined['lines'], w.refined['index'], w.curves, log=w.logmsg, **jp); res['agreement'] = bridge.walk_agreement(res, w.curves, w.tune['eslack']*w.sp_cube.value()); w._join_done(res, t0)
    print(w.txt.toPlainText().splitlines()[-1][:300])
    t0 = time.time(); w.redraw(); print(f'redraw {time.time()-t0:.1f} s, actors {len(w.plotter.renderer.actors)}')
    w.save(); print(w.txt.toPlainText().splitlines()[-1]); zz = np.load('out/orcfull_patch_smoke_lines.npz', allow_pickle=True); print('saved file: lines', len(zz['n']), 'measured flags', zz['measured'].shape, 'bridge nodes', int((~zz['measured']).sum()))
    w.dump(); print(w.txt.toPlainText().splitlines()[-1][:200])
    print('log files:', sorted(os.listdir('out/log'))[-2:])
if __name__ == '__main__': main()
