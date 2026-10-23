"""line 809 of the saved run out/debug/20261023_021229: the trail's fine lines (as saved) next to the refinement with the
anchored walk + chord (the committed fine walk), same tube, cut 17. Renders out/refine809_compare.png (left saved trail, right walk)."""
import sys, os, numpy as np, trimesh, pyvista as pv, subprocess
os.environ['QT_QPA_PLATFORM'] = 'offscreen'; sys.path.insert(0, r'D:\Stuff\PrintSuit\SegmentationTool'); os.chdir(r'D:\Stuff\PrintSuit\SegmentationTool')
D = 'out/debug/20261023_021229'; I = int(os.environ.get('LINE', 809))


def main():
    sys.argv = ['app.py']; import app, refine
    from PySide6 import QtWidgets
    q = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv); w = app.App()
    d = np.load(f'{D}/lines.npz'); n = d['n']; pts = d['pts']; off = np.r_[0, np.cumsum(n)]; curves = [pts[off[i]:off[i+1]] for i in range(len(n))]
    r = np.load(f'{D}/refined.npz'); rn = r['n']; rp_ = r['pts']; ro = np.r_[0, np.cumsum(rn)]; saved = [rp_[ro[k]:ro[k+1]] for k in range(len(rn)) if r['index'][k] == I]
    m = trimesh.load('out/orcfull_patch.stl', process=False); P = np.asarray(m.vertices, float); F = np.asarray(m.faces); h = float((P.max(0)-P.min(0)).max())
    w.chk_trail.setChecked(False); w.chk_curve.setChecked(False); w.sp_rcut.setValue(float(os.environ.get('CUT', 17))); rp = w.refine_params()
    res = refine.refine_all(curves, P, F, 0.087, h, log=print, only=I, flat=True, **rp); walked = res['lines']
    print(f'line {I}: coarse {len(curves[I])} nodes; saved (trail, arc) {len(saved)} fine lines {sum(len(l) for l in saved)} nodes; walk+chord {len(walked)} fine lines {sum(len(l) for l in walked)} nodes')
    np.savez(f'out/refine{I}.npz', coarse=curves[I], saved_n=[len(l) for l in saved], saved_pts=np.vstack(saved), walk_n=[len(l) for l in walked], walk_pts=np.vstack(walked))
    subprocess.run([sys.executable, 'scratch/20261022_trail/render_809.py', str(I)])


if __name__ == '__main__': main()
