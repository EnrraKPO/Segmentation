"""line 318 refined with the app's current parameters, chord vs arc trend, against the saved reference out/app_refine_test.npz"""
import sys, os, numpy as np, trimesh
os.environ['QT_QPA_PLATFORM'] = 'offscreen'; sys.path.insert(0, r'D:\Stuff\PrintSuit\SegmentationTool'); os.chdir(r'D:\Stuff\PrintSuit\SegmentationTool')
sys.argv = ['app.py']; import app, refine
from PySide6 import QtWidgets
from scipy.spatial import cKDTree


def main():
    q = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv); w = app.App()
    d = np.load('out/debug/20261020_144432/lines.npz'); n = d['n']; pts = d['pts']; off = np.r_[0, np.cumsum(n)]; curves = [pts[off[i]:off[i+1]] for i in range(len(n))]
    m = trimesh.load('out/orcfull_patch.stl', process=False); P = np.asarray(m.vertices, float); F = np.asarray(m.faces); h = float((P.max(0)-P.min(0)).max())
    old = np.load('out/app_refine_test.npz'); op = old['pts']
    for tr, cv in ((0, 0), (0, 1), (1, 0), (1, 1)):
        w.chk_trail.setChecked(bool(tr)); w.chk_curve.setChecked(bool(cv)); rp = w.refine_params()
        r = refine.refine_all(curves, P, F, 0.087, h, log=lambda *_: None, only=318, flat=True, **rp); L = r['lines']
        allp = np.vstack(L); dd = cKDTree(op).query(allp)[0]; v = [np.diff(l, axis=0) for l in L if len(l) > 2]
        tr_ = np.concatenate([np.degrees(np.arccos(np.clip(((x[1:]*x[:-1]).sum(1))/(np.linalg.norm(x[1:], axis=1)*np.linalg.norm(x[:-1], axis=1)+1e-12), -1, 1))) for x in v])
        print(f'trail={tr} curve={cv}: {len(L)} fine line(s), nodes {[len(l) for l in L]}; vs reference (156 nodes): median {np.median(dd):.3f} mm, max {dd.max():.3f} mm; turn p90 {np.percentile(tr_, 90):.1f} max {tr_.max():.1f}')
        np.savez(f'out/check318_trail{tr}_curve{cv}.npz', n=[len(l) for l in L], pts=allp)


if __name__ == '__main__': main()
