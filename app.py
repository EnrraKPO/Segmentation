"""Segmentation tool. Three stages, each run on its own and kept:
  1. Find dots   : the voxel grid over the shape; line cells (voxel.py). Parameters: cells, sensitivity.
  2. Build lines : dots -> directions -> chains by merging, aligning and closing gaps (lines.py).
  3. Refine      : more cycles on the current chains without starting over.
Run: python app.py (from this folder)."""
import os, sys, time, traceback
import numpy as np, trimesh
from PySide6 import QtCore, QtWidgets
import pyvista as pv
from pyvistaqt import QtInteractor
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import voxel, lines


def load_mesh(path):
    m = trimesh.load(path, process=True)
    if isinstance(m, trimesh.Scene): m = trimesh.util.concatenate([g for g in m.geometry.values()])
    return np.asarray(m.vertices, np.float64), np.asarray(m.faces, np.int64)


class Worker(QtCore.QThread):
    log = QtCore.Signal(str); done = QtCore.Signal(object); failed = QtCore.Signal(str)
    def __init__(self, fn): super().__init__(); self.fn = fn
    def run(self):
        try: self.done.emit(self.fn(self.log.emit))
        except Exception: self.failed.emit(traceback.format_exc())


class App(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__(); self.setWindowTitle('Segmentation tool'); self.resize(1500, 950)
        self.P = self.F = None; self.path = None; self.res = None; self.state = None; self.curves = None
        root = QtWidgets.QWidget(); self.setCentralWidget(root); lay = QtWidgets.QHBoxLayout(root)
        self.plotter = QtInteractor(root); lay.addWidget(self.plotter.interactor, 1)
        side = QtWidgets.QVBoxLayout(); lay.addLayout(side, 0)
        self.btn_load = QtWidgets.QPushButton('Load STL…'); self.btn_load.clicked.connect(self.load); side.addWidget(self.btn_load)
        self.lbl_file = QtWidgets.QLabel('no model'); self.lbl_file.setWordWrap(True); side.addWidget(self.lbl_file)
        def spin(v, lo, hi, step, dec=0):
            s = QtWidgets.QDoubleSpinBox(); s.setRange(lo, hi); s.setSingleStep(step); s.setDecimals(dec); s.setValue(v); return s
        g1 = QtWidgets.QGroupBox('1. Find dots'); f1 = QtWidgets.QFormLayout(g1)
        self.sp_cells = spin(600, 100, 4000, 100); f1.addRow('Cells per model height', self.sp_cells); self.sp_cells.setToolTip("Density of the grid laid over the model: how many cells span the model's height. Every cell the surface passes through becomes a sample, like a texel. Higher = finer lines and more detail, time grows with the square. 600 is quick (~30 s on 1.5M faces), 1000 is sharper (~2 min).")
        self.sp_sens = spin(6, 0.5, 90, 0.5, 1); f1.addRow('Sensitivity (deg per reference cell)', self.sp_sens); self.sp_sens.setToolTip('How much the surface must turn, in degrees per reference cell (1/1000 of the model height), for a cell to count as a line dot. The same value means the same on any model and density. 6-8 catches cords and creases; 2-4 also catches soft folds and tiny mesh defects.')
        self.btn_dots = QtWidgets.QPushButton('Find dots'); self.btn_dots.clicked.connect(self.find_dots); self.btn_dots.setEnabled(False); f1.addRow(self.btn_dots)
        for w_ in (self.sp_cells, self.sp_sens): f1.labelForField(w_).setToolTip(w_.toolTip())
        side.addWidget(g1)
        g2 = QtWidgets.QGroupBox('2. Build lines'); f2 = QtWidgets.QFormLayout(g2)
        self.sp_reach = spin(4, 1, 20, 1); f2.addRow('Direction reach (cells)', self.sp_reach); self.sp_reach.setToolTip("Neighbourhood, in cells, from which each dot's line direction is judged (all dots within this distance, left and right counted alike). Too small (2) and directions flicker; 4 works; larger smooths the direction over a longer stretch.")
        self.sp_agree = spin(25, 1, 90, 5); f2.addRow('Direction agreement (deg)', self.sp_agree); self.sp_agree.setToolTip('Two dots count as sharing a direction when their directions differ by at most this many degrees. Used for merging across and linking along.')
        self.sp_beside = spin(1.0, 0.1, 10, 0.5, 1); f2.addRow('Beside overlap (cells)', self.sp_beside); self.sp_beside.setToolTip('Two dots merge only when they lie BESIDE each other: their offset along the shared direction is at most this many cells (their strokes overlap, i.e. their perpendicular projections intersect), and their sideways offset is within the merge reach. Dots one behind the other on the line do not merge, so the line is not collapsed along itself.')
        self.sp_mreach = spin(2, 1, 10, 0.5, 1); f2.addRow('Merge reach (cells)', self.sp_mreach); self.sp_mreach.setToolTip('Distance, in cells, within which two dots that share a direction and lie beside each other merge into one dot at their mean. Sets how wide a band is collapsed to one line: 2-3 for cords, larger merges parallel flanks (and may swallow nearby separate lines).')
        self.sp_pull = spin(0.5, 0, 1, 0.1, 1); f2.addRow('Align pull (0..1)', self.sp_pull); self.sp_pull.setToolTip('How far each dot moves towards the midpoint of its neighbours ahead and behind, per cycle (0 = no straightening, 1 = all the way). Removes the cell-by-cell wobble along a chain.')
        self.sp_gap = spin(4, 0, 30, 1); f2.addRow('Gap closing (cells)', self.sp_gap); self.sp_gap.setToolTip('How far, in cells, a chain end may reach straight ahead along its direction to meet a dot of the same direction and link to it. 0 = off. Closes small breaks in a line.')
        self.sp_cycles = spin(6, 1, 30, 1); f2.addRow('Cycles', self.sp_cycles); self.sp_cycles.setToolTip('How many times the merge / align / close-gaps cycle is run. Merging converges in about 6; Refine adds more cycles on the current lines.')
        self.sp_minlen = spin(10, 1, 500, 5); f2.addRow('Minimum chain length (dots)', self.sp_minlen); self.sp_minlen.setToolTip('Chains with fewer dots than this are not drawn. Removes short strays; a dot is about one cell, so 10 is about 1% of the model height at 1000 cells.')
        self.btn_lines = QtWidgets.QPushButton('Build lines (from dots)'); self.btn_lines.clicked.connect(self.build_lines); self.btn_lines.setEnabled(False); f2.addRow(self.btn_lines)
        self.btn_refine = QtWidgets.QPushButton('Refine (more cycles on current lines)'); self.btn_refine.clicked.connect(self.refine); self.btn_refine.setEnabled(False); f2.addRow(self.btn_refine)
        for w_ in (self.sp_reach, self.sp_agree, self.sp_beside, self.sp_mreach, self.sp_pull, self.sp_gap, self.sp_cycles, self.sp_minlen): f2.labelForField(w_).setToolTip(w_.toolTip())
        side.addWidget(g2)
        self.chk_dots = QtWidgets.QCheckBox('Show dots'); self.chk_dots.setChecked(True); self.chk_dots.stateChanged.connect(self.redraw); side.addWidget(self.chk_dots)
        self.chk_curves = QtWidgets.QCheckBox('Show lines'); self.chk_curves.setChecked(True); self.chk_curves.stateChanged.connect(self.redraw); side.addWidget(self.chk_curves)
        self.chk_dirs = QtWidgets.QCheckBox('Show dot directions (short strokes)'); self.chk_dirs.setChecked(False); self.chk_dirs.stateChanged.connect(self.redraw); side.addWidget(self.chk_dirs)
        self.chk_links = QtWidgets.QCheckBox('Show links (ahead/behind)'); self.chk_links.setChecked(False); self.chk_links.stateChanged.connect(self.redraw); side.addWidget(self.chk_links)
        self.btn_save = QtWidgets.QPushButton('Save lines…'); self.btn_save.clicked.connect(self.save); self.btn_save.setEnabled(False); side.addWidget(self.btn_save)
        self.txt = QtWidgets.QPlainTextEdit(); self.txt.setReadOnly(True); self.txt.setMaximumWidth(380); side.addWidget(self.txt, 1)
        self.plotter.set_background('white'); self.plotter.add_axes()

    def logmsg(self, s): self.txt.appendPlainText(s)
    def busy(self, b):
        self.btn_dots.setEnabled(not b and self.P is not None); self.btn_lines.setEnabled(not b and self.res is not None); self.btn_refine.setEnabled(not b and self.state is not None)

    def load(self):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, 'Open STL', '', 'Meshes (*.stl *.obj *.ply *.3mf)')
        if not path: return
        try: self.P, self.F = load_mesh(path)
        except Exception as ex: QtWidgets.QMessageBox.critical(self, 'Load failed', str(ex)); return
        self.path = path; self.lbl_file.setText(os.path.basename(path)); self.res = self.state = self.curves = None; self.btn_save.setEnabled(False)
        h = float((self.P.max(0)-self.P.min(0)).max()); self.logmsg(f'loaded {os.path.basename(path)}: {len(self.F)} faces, height {h:.1f} mm')
        self.busy(False); self.redraw()

    def _start(self, fn, on_done):
        self.busy(True); self.worker = Worker(fn); self.worker.log.connect(self.logmsg); self.worker.done.connect(on_done); self.worker.failed.connect(self.failed); self.worker.start()
    def failed(self, msg): self.logmsg(msg); self.busy(False)

    def find_dots(self):
        cells, sens = int(self.sp_cells.value()), self.sp_sens.value(); self.logmsg(f'--- find dots: cells={cells}, sensitivity={sens}'); t0 = time.time()
        self._start(lambda log: voxel.run(self.P, self.F, cells, sens, log), lambda r: self._dots_done(r, t0))
    def _dots_done(self, r, t0):
        self.res = r; self.state = None; self.curves = None; self.logmsg(f'dots done in {time.time()-t0:.1f} s'); self.busy(False); self.redraw()

    def chain_params(self):
        return dict(agree=self.sp_agree.value(), beside=self.sp_beside.value(), merge_reach=self.sp_mreach.value(), pull=self.sp_pull.value(), gap=self.sp_gap.value(), cycles=int(self.sp_cycles.value()))

    def build_lines(self):
        reach, cp, ml = self.sp_reach.value(), self.chain_params(), int(self.sp_minlen.value()); self.logmsg('--- build lines: reach=%g, ' % reach + ', '.join(f'{k}={v}' for k, v in cp.items())); t0 = time.time()
        def fn(log):
            st = lines.prepare(self.res, reach); lines.refine(st, log=log, **cp); return st, lines.finish(st, self.P, self.F, self.res, ml, log)
        self._start(fn, lambda r: self._lines_done(r, t0))
    def refine(self):
        cp, ml = self.chain_params(), int(self.sp_minlen.value()); self.logmsg('--- refine: ' + ', '.join(f'{k}={v}' for k, v in cp.items())); t0 = time.time(); st = self.state
        def fn(log):
            lines.refine(st, log=log, **cp); return st, lines.finish(st, self.P, self.F, self.res, ml, log)
        self._start(fn, lambda r: self._lines_done(r, t0))
    def _lines_done(self, r, t0):
        self.state, self.curves = r; self.logmsg(f'lines done in {time.time()-t0:.1f} s'); self.busy(False); self.btn_save.setEnabled(True); self.redraw()

    def redraw(self):
        if self.P is None: return
        cam = self.plotter.camera_position; had = self.res is not None; self.plotter.clear()
        mesh = pv.PolyData(self.P, np.hstack([np.full((len(self.F), 1), 3), self.F]).ravel())
        self.plotter.add_mesh(mesh, color='lightgrey', smooth_shading=True)
        if self.res is not None and self.state is None and self.chk_dirs.isChecked():
            lo, size = self.res['lo'], self.res['size']; li = np.nonzero(self.res['line'])[0]; X = (self.res['idx'][li]+0.5).astype(float)
            if getattr(self, '_dirs_cache', (None,))[0] is not self.res: self._dirs_cache = (self.res, lines.directions(X, self.res['N'][li], self.sp_reach.value()))
            D = self._dirs_cache[1]; W = lo+X*size; seg = np.stack([W-0.6*size*D, W+0.6*size*D], 1).reshape(-1, 3); n = len(X)
            self.plotter.add_mesh(pv.PolyData(seg, lines=np.c_[np.full(n, 2), np.arange(0, 2*n, 2), np.arange(1, 2*n, 2)].ravel()), color='green', line_width=2)
        if self.res is not None and self.chk_dots.isChecked():
            if self.state is not None: C = self.res['lo']+self.state.X[self.state.alive]*self.res['size']
            else: C = self.res['lo']+(self.res['idx'][self.res['line']]+0.5)*self.res['size']
            self.plotter.add_points(C, color='red', point_size=3, render_points_as_spheres=False)
        if self.state is not None and (self.chk_dirs.isChecked() or self.chk_links.isChecked()):
            st, lo, size = self.state, self.res['lo'], self.res['size']; k = np.nonzero(st.alive)[0]; W = lo+st.X[k]*size
            if self.chk_dirs.isChecked():
                seg = np.stack([W-0.6*size*st.D[k], W+0.6*size*st.D[k]], 1).reshape(-1, 3); n = len(k)
                self.plotter.add_mesh(pv.PolyData(seg, lines=np.c_[np.full(n, 2), np.arange(0, 2*n, 2), np.arange(1, 2*n, 2)].ravel()), color='green', line_width=2)
            if self.chk_links.isChecked():
                lk = k[st.ahead[k] >= 0]
                if len(lk):
                    seg = np.stack([lo+st.X[lk]*size, lo+st.X[st.ahead[lk]]*size], 1).reshape(-1, 3); n = len(lk)
                    self.plotter.add_mesh(pv.PolyData(seg, lines=np.c_[np.full(n, 2), np.arange(0, 2*n, 2), np.arange(1, 2*n, 2)].ravel()), color='orange', line_width=2)
        if self.curves and self.chk_curves.isChecked():
            pts = np.vstack(self.curves); n = [len(c) for c in self.curves]; off = np.cumsum([0]+n[:-1])
            cells = np.concatenate([np.r_[k, np.arange(o, o+k)] for k, o in zip(n, off)])
            self.plotter.add_mesh(pv.PolyData(pts, lines=cells), color='blue', line_width=3)
        if cam is not None and had: self.plotter.camera_position = cam
        else: self.plotter.reset_camera()

    def save(self):
        if not self.curves: return
        out = os.path.splitext(self.path)[0]+'_lines.npz'
        prm = dict(cells=int(self.sp_cells.value()), sensitivity=self.sp_sens.value(), reach=self.sp_reach.value(), min_len=int(self.sp_minlen.value()), cycles_done=self.state.cycles_done, **self.chain_params())
        np.savez(out, n=np.array([len(c) for c in self.curves]), pts=np.vstack(self.curves), params=str(prm)); self.logmsg(f'saved {out}')


if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv); w = App(); w.show(); sys.exit(app.exec())
