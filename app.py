"""Segmentation tool. Two stages, each run on its own and kept:
  1. Find dots   : cubes over the surface, each reading its own triangle's facing; a dot where the facing shift
                   changes sharply from one cube to the next (creases.find_dots). Parameters: cube size (mm), cut.
  2. Build lines : one smooth thick stroke through each run of dots (creases.strokes), then stroke ends joined at
                   corners, gaps and tees (creases.join). Parameters: thickness, stretch, min row ahead, join reach, all in mm.
Run: python app.py (from this folder)."""
import os, sys, time, traceback
import numpy as np, trimesh
from PySide6 import QtCore, QtWidgets
import pyvista as pv
from pyvistaqt import QtInteractor
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import creases


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


TIPS = dict(
    cells="Edge of the cubes laid over the surface, in mm. Should be a few times smaller than the model's triangles but not far below them (orc: 0.41 mm triangles, 0.087 mm cubes; leotard: 0.072 mm triangles, 0.050 mm cubes). Halving it gives four times the cubes and the time.",
    cut="How much the facing shift must change from one cube to the next for the cube to be a dot. Measured on the orc: plate edges ~80, nail edges ~40, skin wrinkles 18-22, soft toe flanks 8-15, flat skin 2-5. 17 keeps creases and drops soft forms.",
    flat="On: the cube reads the facing of the flat triangle under it, the surface the print has. Off: the render's blended facing, which spreads every fold over the neighbouring triangles and hides faint creases.",
    W="Half the thickness of the stroke laid over the dots, in mm. Dots within this distance sideways belong to the stroke. 0.4 covers the band a crease makes at 0.087 mm cubes; much larger merges neighbouring parallel creases into one.",
    G="How far a stroke end is extended to meet another stroke: at a corner (two ends whose headings cross), across a gap (two ends facing each other) or into the side of another stroke. 0 = no joining.",
    L="How far ahead and behind the stroke looks, in mm, when deciding its position and heading. Longer = smoother line, but a long straight stretch cuts corners on tight curves. 1.0 worked on nail rings and seams.",
    K="The stroke needs at least this length of dot row ahead of it to continue, and in its window to start. Shorter follows sparse dot rows (faint creases) but may start strokes on scattered dots.",
    align="The dots in the window must be arranged along a line: the cloud must be at least this many times longer than it is wide. Keeps strokes from starting on or shooting toward scattered dots.",
    agree="The stroke stops where its heading would turn more than this in one step: a corner.",
    minlen="Strokes with fewer nodes than this are dropped. One node per cube.",
    snap="Moves each node onto the nearest surface point before drawing, so lines on sharp convex edges are not hidden inside the solid. Judge placement with this on; lifting lines off the surface shifts them sideways.",
)


class App(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__(); self.setWindowTitle('Segmentation tool'); self.resize(1500, 950)
        self.P = self.F = None; self.path = None; self.res = None; self.curves = None; self.mesh_t = None
        root = QtWidgets.QWidget(); self.setCentralWidget(root); lay = QtWidgets.QHBoxLayout(root)
        self.plotter = QtInteractor(root); lay.addWidget(self.plotter.interactor, 1)
        side = QtWidgets.QVBoxLayout(); lay.addLayout(side, 0)
        self.btn_load = QtWidgets.QPushButton('Load STL…'); self.btn_load.clicked.connect(self.load); side.addWidget(self.btn_load)
        self.lbl_file = QtWidgets.QLabel('no model'); self.lbl_file.setWordWrap(True); side.addWidget(self.lbl_file)

        def spin(v, lo, hi, step, dec=0):
            s = QtWidgets.QDoubleSpinBox(); s.setRange(lo, hi); s.setSingleStep(step); s.setDecimals(dec); s.setValue(v); return s

        def row(form, label, widget, tip):
            widget.setToolTip(tip); form.addRow(label, widget); form.labelForField(widget).setToolTip(tip)

        g1 = QtWidgets.QGroupBox('1. Find dots'); f1 = QtWidgets.QFormLayout(g1)
        self.sp_cube = spin(0.087, 0.01, 5, 0.01, 3); row(f1, 'Cube size (mm)', self.sp_cube, TIPS['cells'])
        self.sp_cut = spin(17, 1, 90, 1); row(f1, 'Cut (deg)', self.sp_cut, TIPS['cut'])
        self.chk_flat = QtWidgets.QCheckBox("Read each triangle's own facing (no corner blending)"); self.chk_flat.setChecked(True); self.chk_flat.setToolTip(TIPS['flat']); f1.addRow(self.chk_flat)
        self.btn_dots = QtWidgets.QPushButton('Find dots'); self.btn_dots.clicked.connect(self.find_dots); self.btn_dots.setEnabled(False); f1.addRow(self.btn_dots)
        side.addWidget(g1)

        g2 = QtWidgets.QGroupBox('2. Build lines'); f2 = QtWidgets.QFormLayout(g2)
        self.sp_W = spin(0.4, 0.02, 5, 0.05, 2); row(f2, 'Stroke half-width (mm)', self.sp_W, TIPS['W'])
        self.sp_L = spin(1.0, 0.1, 10, 0.1, 2); row(f2, 'Stroke stretch (mm)', self.sp_L, TIPS['L'])
        self.sp_K = spin(0.26, 0.02, 5, 0.05, 2); row(f2, 'Min dot row ahead (mm)', self.sp_K, TIPS['K'])
        self.sp_align = spin(3, 1, 10, 0.5, 1); row(f2, 'Alignment (long/wide)', self.sp_align, TIPS['align'])
        self.sp_agree = spin(30, 1, 90, 5); row(f2, 'Max turn per step (deg)', self.sp_agree, TIPS['agree'])
        self.sp_minlen = spin(10, 1, 500, 5); row(f2, 'Minimum stroke length (nodes)', self.sp_minlen, TIPS['minlen'])
        self.sp_G = spin(1.0, 0, 10, 0.1, 2); row(f2, 'Join reach (mm)', self.sp_G, TIPS['G'])
        self.btn_lines = QtWidgets.QPushButton('Build lines (from dots)'); self.btn_lines.clicked.connect(self.build_lines); self.btn_lines.setEnabled(False); f2.addRow(self.btn_lines)
        side.addWidget(g2)

        self.chk_dots = QtWidgets.QCheckBox('Show dots'); self.chk_dots.setChecked(True); self.chk_dots.stateChanged.connect(self.redraw); side.addWidget(self.chk_dots)
        self.chk_curves = QtWidgets.QCheckBox('Show lines'); self.chk_curves.setChecked(True); self.chk_curves.stateChanged.connect(self.redraw); side.addWidget(self.chk_curves)
        self.chk_snap = QtWidgets.QCheckBox('Draw lines on the surface (a few seconds)'); self.chk_snap.setChecked(True); self.chk_snap.setToolTip(TIPS['snap']); self.chk_snap.stateChanged.connect(self.redraw); side.addWidget(self.chk_snap)
        self.btn_save = QtWidgets.QPushButton('Save lines…'); self.btn_save.clicked.connect(self.save); self.btn_save.setEnabled(False); side.addWidget(self.btn_save)
        self.btn_dump = QtWidgets.QPushButton('Dump for debugging'); self.btn_dump.clicked.connect(self.dump); self.btn_dump.setToolTip('Writes dots, lines, log, settings and a screenshot of the viewport to out/debug/<time>/'); side.addWidget(self.btn_dump)
        self.txt = QtWidgets.QPlainTextEdit(); self.txt.setReadOnly(True); self.txt.setMaximumWidth(380); side.addWidget(self.txt, 1)
        self.plotter.set_background('white'); self.plotter.add_axes()

    def logmsg(self, s): self.txt.appendPlainText(s)
    def busy(self, b):
        self.btn_dots.setEnabled(not b and self.P is not None); self.btn_lines.setEnabled(not b and self.res is not None)

    def load(self):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, 'Open STL', '', 'Meshes (*.stl *.obj *.ply *.3mf)')
        if not path: return
        try: self.P, self.F = load_mesh(path)
        except Exception as ex: QtWidgets.QMessageBox.critical(self, 'Load failed', str(ex)); return
        self.path = path; self.lbl_file.setText(os.path.basename(path)); self.res = self.curves = None; self.mesh_t = None; self.btn_save.setEnabled(False)
        h = float((self.P.max(0)-self.P.min(0)).max())
        e = np.linalg.norm(self.P[self.F[:, 1]]-self.P[self.F[:, 0]], axis=1)
        self.logmsg(f'loaded {os.path.basename(path)}: {len(self.F)} faces, height {h:.1f} mm, typical triangle edge {np.median(e):.3f} mm')
        self.busy(False); self.redraw()

    def _start(self, fn, on_done):
        self.busy(True); self.worker = Worker(fn); self.worker.log.connect(self.logmsg); self.worker.done.connect(on_done); self.worker.failed.connect(self.failed); self.worker.start()
    def failed(self, msg): self.logmsg(msg); self.busy(False)

    def find_dots(self):
        size, cut, flat = self.sp_cube.value(), self.sp_cut.value(), self.chk_flat.isChecked()
        h = float((self.P.max(0)-self.P.min(0)).max()); cells = max(1, int(round(h/size)))
        self.logmsg(f'--- find dots: cube={size} mm ({cells} cells over {h:.1f} mm), cut={cut}, flat={flat}'); t0 = time.time()
        self._start(lambda log: creases.find_dots(self.P, self.F, cells, cut, 0.0, flat, log, height=h), lambda r: self._dots_done(r, t0))
    def _dots_done(self, r, t0):
        self.res = r; self.curves = None; self.logmsg(f'dots done in {time.time()-t0:.1f} s'); self.busy(False); self.redraw()

    def stroke_params(self):
        return dict(W=self.sp_W.value(), L=self.sp_L.value(), K=self.sp_K.value(), align=self.sp_align.value(), agree=self.sp_agree.value(), minlen=int(self.sp_minlen.value()))
    def all_params(self): return dict(**self.stroke_params(), join_mm=self.sp_G.value())

    def build_lines(self):
        sp = self.stroke_params(); self.logmsg('--- build lines: ' + ', '.join(f'{k}={v}' for k, v in sp.items())); t0 = time.time()
        G = self.sp_G.value()
        self._start(lambda log: (lambda cv: creases.join(cv, G, log=log) if G > 0 else cv)(creases.strokes(self.res, log=log, **sp)), lambda r: self._lines_done(r, t0))
    def _lines_done(self, r, t0):
        self.curves = r; self.logmsg(f'lines done in {time.time()-t0:.1f} s'); self.busy(False); self.btn_save.setEnabled(True); self.redraw()

    def redraw(self):
        if self.P is None: return
        cam = self.plotter.camera_position; had = self.res is not None; self.plotter.clear()
        mesh = pv.PolyData(self.P, np.hstack([np.full((len(self.F), 1), 3), self.F]).ravel())
        self.plotter.add_mesh(mesh, color='lightgrey', smooth_shading=True)
        if self.res is not None and self.chk_dots.isChecked():
            C = self.res['lo']+(self.res['idx'][self.res['line']]+0.5)*self.res['size']
            C = C+0.5*self.res['size']*self.res['N'][self.res['line']]   # half a cube out along the facing, so dots on convex edges are not buried
            self.plotter.add_points(C, color='red', point_size=5, render_points_as_spheres=True)
        if self.curves and self.chk_curves.isChecked():
            pts = np.vstack(self.curves)
            if self.chk_snap.isChecked():
                if self.mesh_t is None: self.mesh_t = trimesh.Trimesh(self.P, self.F, process=False)
                Cq, dq, tq = trimesh.proximity.closest_point(self.mesh_t, pts); pts = Cq+0.01*self.mesh_t.face_normals[tq]
            n = [len(c) for c in self.curves]; off = np.cumsum([0]+n[:-1])
            cells = np.concatenate([np.r_[k, np.arange(o, o+k)] for k, o in zip(n, off)])
            self.plotter.add_mesh(pv.PolyData(pts, lines=cells), color='blue', line_width=3)
        if cam is not None and had: self.plotter.camera_position = cam
        else: self.plotter.reset_camera()

    def save(self):
        if not self.curves: return
        out = os.path.splitext(self.path)[0]+'_lines.npz'
        prm = dict(cube_mm=self.sp_cube.value(), cut=self.sp_cut.value(), flat=self.chk_flat.isChecked(), **self.all_params())
        np.savez(out, n=np.array([len(c) for c in self.curves]), pts=np.vstack(self.curves), params=str(prm)); self.logmsg(f'saved {out}')

    def dump(self):
        """Everything on screen, to out/debug/<time>/: dots.npz (cube centres, drawn positions, grid), lines.npz (points,
        run lengths), log.txt, settings.json (parameters, model, camera) and viewport.png. Works with any stage done."""
        d = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out', 'debug', time.strftime('%Y%m%d_%H%M%S')); os.makedirs(d, exist_ok=True)
        try:
            prm = dict(model=self.path, faces=None if self.F is None else int(len(self.F)), cube_mm=self.sp_cube.value(), cut=self.sp_cut.value(), flat=self.chk_flat.isChecked(),
                       **self.all_params(), show_dots=self.chk_dots.isChecked(), show_lines=self.chk_curves.isChecked(), snap=self.chk_snap.isChecked(),
                       camera=[list(map(float, v)) for v in self.plotter.camera_position], dots=0, lines=0)
            if self.res is not None:
                r = self.res; C = r['lo']+(r['idx'][r['line']]+0.5)*r['size']; D = C+0.5*r['size']*r['N'][r['line']]
                np.savez(os.path.join(d, 'dots.npz'), centres=C, drawn=D, idx=r['idx'][r['line']], N=r['N'][r['line']], turn=r['turn'][r['line']], lo=r['lo'], size=r['size'])
                prm['dots'] = int(len(C))
            if self.curves:
                np.savez(os.path.join(d, 'lines.npz'), n=np.array([len(c) for c in self.curves]), pts=np.vstack(self.curves)); prm['lines'] = len(self.curves)
            with open(os.path.join(d, 'log.txt'), 'w', encoding='utf-8') as f: f.write(self.txt.toPlainText())
            import json
            with open(os.path.join(d, 'settings.json'), 'w', encoding='utf-8') as f: json.dump(prm, f, indent=1)
            self.plotter.screenshot(os.path.join(d, 'viewport.png'))
            self.logmsg(f'dumped {prm["dots"]} dots, {prm["lines"]} lines, log, settings, screenshot to {d}')
        except Exception: self.logmsg(traceback.format_exc())


if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv); w = App(); w.show(); sys.exit(app.exec())
