"""Segmentation tool, second version (2026-10-24): built stage by stage from what the first app and its experiments
settled. Stage 1, the field (dots.py): cubes over the surface, each with the facing of the flat triangles inside it and
the angle the facing folds by from one cube to the next. Stage 2, the cut: a cube is a dot where that angle is at least
the cut. The cut is a filter on the field held in memory, so sliding it redraws the dots at once.

The panel holds only the handles expected to change; everything else (display, cores) is in Settings. Dot sets are
saved to out/dots/ whenever Find dots finishes and the last set of a model is loaded with it. Run: python app.py."""
import os, sys, time, traceback, json
import numpy as np, trimesh
from PySide6 import QtCore, QtWidgets
import pyvista as pv
from pyvistaqt import QtInteractor
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dots

ROOT = os.path.dirname(os.path.abspath(__file__))
DOTS_DIR = os.path.join(ROOT, 'out', 'dots')   # saved dot sets; last.json maps model name -> last file
LOG_DIR = os.path.join(ROOT, 'out', 'log')     # one log file per session, written as it happens
DEBUG_DIR = os.path.join(ROOT, 'out', 'debug')   # dumps: everything on screen, one folder per dump

TIPS = dict(
    reach='How far the fold angle looks, in cubes: the arrows compared sit this many cubes apart. 1 reads crisp creases only; a bend spread over several cubes reads its full angle when the reach spans it, at the cost of a band that many cubes wider on every crease.',
    cube="Edge of the cubes laid over the surface, in mm. A few times smaller than the model's triangles but not far below them (orc: 0.41 mm triangles, 0.087 mm cubes). Halving it gives four times the cubes and the time. The only limit is the grid: fewer than 65 536 cubes along the model's height.",
    cut="How far the facing must fold from one cube to the next, in degrees, for the cube to be a dot. Measured on the orc: plate edges ~80, nail edges ~40, skin wrinkles 18-22, soft toe flanks 8-15, flat skin 2-5. 17 keeps creases and drops soft forms. A filter on the field: slide it and the dots redraw. A loaded dot set holds only the cubes above its own cut, so below that nothing more appears.",
    blend="On: a cube holding several triangles takes their facings blended by the area each covers inside the cube, so the arrow changes gradually as a crease crosses from one cube to the next (orc: 2.9M of 23.6M cubes, about 40 s more). Off: the cube takes the one triangle nearest its centre.",
)

# Settings: everything that is not expected to change, as (key, label, default, lo, hi, step, decimals, meaning).
SETTINGS = [
    ('dotlift', 'Dot lift (cubes)', 0.5, 0, 5, 0.1, 1, 'Dots are drawn this many cubes out along their facing, so dots on convex edges are not buried in the surface.'),
    ('dotsize', 'Dot size (px)', 5, 1, 20, 1, 0, 'Point size of the drawn dots.'),
    ('workers', 'Cores (0 = all)', 0, 0, 256, 1, 0, 'Processes the cube stage is spread over.'),
]
SETTINGS_DEFAULTS = {s[0]: s[2] for s in SETTINGS}


def open_session_log():
    """app_<start>.txt gets every line the text box shows, flushed at once; crash_<start>.txt gets the Python stack if the
    process dies below Python (faulthandler); uncaught exceptions on the main thread go to the app log too."""
    import faulthandler
    os.makedirs(LOG_DIR, exist_ok=True); stamp = time.strftime('%Y%m%d_%H%M%S')
    path = os.path.join(LOG_DIR, f'app_{stamp}.txt'); f = open(path, 'a', encoding='utf-8')
    crash = open(os.path.join(LOG_DIR, f'crash_{stamp}.txt'), 'w'); faulthandler.enable(file=crash, all_threads=True)
    def write(s):
        try: f.write(time.strftime('%H:%M:%S ')+s+'\n'); f.flush()
        except Exception: pass
    def hook(t, v, tb): write('uncaught exception on the main thread:\n'+''.join(traceback.format_exception(t, v, tb))); sys.__excepthook__(t, v, tb)
    sys.excepthook = hook
    write(f'session started, python {sys.version.split()[0]}, pid {os.getpid()}')
    return write, path


def load_mesh(path):
    m = trimesh.load(path, process=True)
    if isinstance(m, trimesh.Scene): m = trimesh.util.concatenate([g for g in m.geometry.values()])
    return np.asarray(m.vertices, np.float64), np.asarray(m.faces, np.int64)


def save_dots(fld, cut, model, cube, blend, reach=1, path=None):
    """The dots alone (cubes at or above the cut): idx, N, turn, plus lo, size, cut and the settings. Same file layout as the
    first app's, so its sets load here. Default path: out/dots/<model>_<cube>mm_cut<cut>[_nearest]_<time>.npz, noted as the
    model's last set."""
    m = dots.dots(fld, cut); base = os.path.splitext(os.path.basename(model))[0]
    if path is None:
        os.makedirs(DOTS_DIR, exist_ok=True)
        name = f'{base}_{cube:g}mm_reach{reach}_cut{cut:g}{"" if blend else "_nearest"}_{time.strftime("%Y%m%d_%H%M%S")}.npz'; path = os.path.join(DOTS_DIR, name)
        lj = os.path.join(DOTS_DIR, 'last.json'); last = json.load(open(lj)) if os.path.exists(lj) else {}
        last[os.path.basename(model)] = name; json.dump(last, open(lj, 'w'), indent=1)
    np.savez(path, idx=fld['idx'][m], N=fld['N'][m], turn=fld['turn'][m], lo=fld['lo'], size=fld['size'], cut=cut, cube=cube, flat=True, blend=blend, reach=reach, model=model)
    return path


def load_dots(path):
    """A saved dot set back as a field holding only its dots. Returns field, its cut, cube (mm), blend, reach, model path.
    Sets from before 2026-10-24 evening have no reach: their fold angle was the straight-line reading, flagged as reach 0."""
    with np.load(path, allow_pickle=True) as z:
        fld = dict(idx=z['idx'].astype(np.int32), N=z['N'].astype(np.float32), turn=z['turn'].astype(np.float32), lo=z['lo'], size=float(z['size']))
        return fld, float(z['cut']), float(z['cube']), bool(z['blend']) if 'blend' in z.files else False, int(z['reach']) if 'reach' in z.files else 0, str(z['model'])   # sets before 2026-10-22 took the nearest triangle


def last_dots(model):
    lj = os.path.join(DOTS_DIR, 'last.json')
    if not os.path.exists(lj): return None
    name = json.load(open(lj)).get(os.path.basename(model)); p = name and os.path.join(DOTS_DIR, name)
    return p if p and os.path.exists(p) else None


class Worker(QtCore.QThread):
    log = QtCore.Signal(str); done = QtCore.Signal(object); failed = QtCore.Signal(str)
    def __init__(self, fn): super().__init__(); self.fn = fn
    def run(self):
        try: self.done.emit(self.fn(self.log.emit))
        except Exception: self.failed.emit(traceback.format_exc())


class SettingsDialog(QtWidgets.QDialog):
    """the values not expected to change; edits land in app.settings at once"""
    def __init__(self, app):
        super().__init__(app); self.app = app; self.setWindowTitle('Settings'); self.resize(420, 300)
        lay = QtWidgets.QVBoxLayout(self); form = QtWidgets.QFormLayout(); lay.addLayout(form); self.spins = {}
        for key, label, default, lo, hi, step, dec, tip in SETTINGS:
            sp = QtWidgets.QDoubleSpinBox(); sp.setRange(lo, hi); sp.setSingleStep(step); sp.setDecimals(dec); sp.setValue(app.settings[key]); sp.setToolTip(tip)
            sp.valueChanged.connect(lambda v, k=key: self._set(k, v)); form.addRow(label, sp); form.labelForField(sp).setToolTip(tip); self.spins[key] = sp
        self.chk_flat = QtWidgets.QCheckBox('Flat shading (each triangle lit by its own facing, as the print has it)'); self.chk_flat.setChecked(app.flat_shading); self.chk_flat.stateChanged.connect(lambda s: self._shade(bool(s))); lay.addWidget(self.chk_flat)
        lay.addStretch(1); btns = QtWidgets.QHBoxLayout(); lay.addLayout(btns)
        b = QtWidgets.QPushButton('Reset to defaults'); b.clicked.connect(self.reset); btns.addWidget(b)
        c = QtWidgets.QPushButton('Close'); c.clicked.connect(self.close); btns.addWidget(c)
    def _set(self, k, v): self.app.settings[k] = v; self.app.redraw() if k in ('dotlift', 'dotsize') else None
    def _shade(self, on): self.app.flat_shading = on; self.app.redraw()
    def reset(self):
        for k, sp in self.spins.items(): sp.setValue(SETTINGS_DEFAULTS[k])
        self.chk_flat.setChecked(False)


class App(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__(); self.setWindowTitle('Segmentation tool v2'); self.resize(1500, 950)
        self.P = self.F = None; self.path = None; self.field = None; self.order = None; self.cube = None; self.settings = dict(SETTINGS_DEFAULTS); self.flat_shading = False; self.dlg = None
        try: self.file_log, log_path = open_session_log()
        except Exception: self.file_log, log_path = None, None
        root = QtWidgets.QWidget(); self.setCentralWidget(root); lay = QtWidgets.QHBoxLayout(root)
        self.plotter = QtInteractor(root); lay.addWidget(self.plotter.interactor, 1)
        side = QtWidgets.QVBoxLayout(); lay.addLayout(side, 0)
        self.btn_load = QtWidgets.QPushButton('Load mesh…'); self.btn_load.clicked.connect(self.load); side.addWidget(self.btn_load)
        self.lbl_file = QtWidgets.QLabel('no model'); self.lbl_file.setWordWrap(True); side.addWidget(self.lbl_file)

        def spin(v, lo, hi, step, dec=0):
            s = QtWidgets.QDoubleSpinBox(); s.setRange(lo, hi); s.setSingleStep(step); s.setDecimals(dec); s.setValue(v); return s
        def row(form, label, widget, tip):
            widget.setToolTip(tip); form.addRow(label, widget); form.labelForField(widget).setToolTip(tip)

        g1 = QtWidgets.QGroupBox('1. Field'); f1 = QtWidgets.QFormLayout(g1)
        self.sp_cube = spin(0.087, 0.0001, 5, 0.001, 4); row(f1, 'Cube size (mm)', self.sp_cube, TIPS['cube'])
        self.sp_reach = QtWidgets.QSpinBox(); self.sp_reach.setRange(1, 50); self.sp_reach.setValue(1); row(f1, 'Reach (cubes)', self.sp_reach, TIPS['reach'])
        self.chk_blend = QtWidgets.QCheckBox('Blend the triangles inside a cube by area'); self.chk_blend.setChecked(True); self.chk_blend.setToolTip(TIPS['blend']); f1.addRow(self.chk_blend)
        self.btn_field = QtWidgets.QPushButton('Compute the field'); self.btn_field.clicked.connect(self.compute_field); self.btn_field.setToolTip('Cubes over the surface with their arrows and fold angles. The cut below does not enter here.'); f1.addRow(self.btn_field)
        self.btn_load_dots = QtWidgets.QPushButton('Load dots…'); self.btn_load_dots.clicked.connect(self.load_dots); self.btn_load_dots.setToolTip('A saved set, from this app or the first one. Opening a mesh loads its last set by itself.'); f1.addRow(self.btn_load_dots)
        side.addWidget(g1)

        g2 = QtWidgets.QGroupBox('2. Cut'); f2 = QtWidgets.QFormLayout(g2)
        cw = QtWidgets.QWidget(); ch = QtWidgets.QHBoxLayout(cw); ch.setContentsMargins(0, 0, 0, 0)
        self.sl_cut = QtWidgets.QSlider(QtCore.Qt.Horizontal); self.sl_cut.setRange(10, 1800); self.sl_cut.setValue(170); ch.addWidget(self.sl_cut, 1)   # tenths of a degree
        self.sp_cut = spin(17, 1, 180, 1, 1); self.sp_cut.setMaximumWidth(70); ch.addWidget(self.sp_cut); row(f2, 'Cut (deg)', cw, TIPS['cut']); self.sp_cut.setToolTip(TIPS['cut'])
        self.sl_cut.valueChanged.connect(lambda v: self.sp_cut.setValue(v/10)); self.sp_cut.valueChanged.connect(self._cut_changed)
        self.lbl_dots = QtWidgets.QLabel('no field'); f2.addRow(self.lbl_dots)
        self.btn_save = QtWidgets.QPushButton('Save dots…'); self.btn_save.clicked.connect(self.save_dots); self.btn_save.setToolTip('Writes the dots at this cut where you choose. Every computed field is also saved at its cut to out/dots/ by itself.'); f2.addRow(self.btn_save)
        side.addWidget(g2)
        self.cut_timer = QtCore.QTimer(self); self.cut_timer.setSingleShot(True); self.cut_timer.setInterval(60); self.cut_timer.timeout.connect(self.redraw)   # one redraw per slider pause, not per tick

        self.chk_dots = QtWidgets.QCheckBox('Show dots'); self.chk_dots.setChecked(True); self.chk_dots.stateChanged.connect(self.redraw); side.addWidget(self.chk_dots)
        self.btn_dump = QtWidgets.QPushButton('Dump for debugging'); self.btn_dump.clicked.connect(self.dump); self.btn_dump.setToolTip('Writes the whole field (every cube, its arrow and fold angle), the settings, the camera, the log and a screenshot to out/debug/<time>/.'); side.addWidget(self.btn_dump)
        self.btn_settings = QtWidgets.QPushButton('Settings…'); self.btn_settings.clicked.connect(self.open_settings); self.btn_settings.setToolTip('Display and machine settings: nothing here changes a dot.'); side.addWidget(self.btn_settings)
        self.txt = QtWidgets.QPlainTextEdit(); self.txt.setReadOnly(True); self.txt.setMaximumWidth(380); side.addWidget(self.txt, 1)
        self.plotter.set_background('white'); self.plotter.add_axes()
        self.logmsg(f'session log: {log_path}' if log_path else 'session log could not be opened (out/log)'); self.busy(False)

    def logmsg(self, s):
        self.txt.appendPlainText(s)
        if self.file_log: self.file_log(s)
    def busy(self, b):
        self.btn_load.setEnabled(not b); self.btn_field.setEnabled(not b and self.P is not None); self.btn_load_dots.setEnabled(not b and self.P is not None); self.btn_save.setEnabled(not b and self.field is not None)
    def open_settings(self):
        if self.dlg is None: self.dlg = SettingsDialog(self)
        self.dlg.show(); self.dlg.raise_()

    def load(self):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, 'Open mesh', '', 'Meshes (*.stl *.obj *.ply *.3mf)')
        if not path: return
        try: self.P, self.F = load_mesh(path)
        except Exception as ex: QtWidgets.QMessageBox.critical(self, 'Load failed', str(ex)); return
        self.path = path; self.lbl_file.setText(os.path.basename(path)); self._set_field(None)
        h = float((self.P.max(0)-self.P.min(0)).max()); e = np.linalg.norm(self.P[self.F[:, 1]]-self.P[self.F[:, 0]], axis=1)
        self.logmsg(f'loaded {os.path.basename(path)}: {len(self.F)} faces, height {h:.1f} mm, typical triangle edge {np.median(e):.3f} mm')
        p = last_dots(path)
        if p: self._use_dots(p)
        self.busy(False); self.redraw(reset=True)

    def _set_field(self, fld):
        """the field held, with its cubes ordered by fold angle so a cut is a prefix of the order (searchsorted, instant)"""
        self.field = fld; self.order = None if fld is None else np.argsort(-fld['turn'], kind='stable')
        self.sorted_turn = None if fld is None else -fld['turn'][self.order]

    def _use_dots(self, p):
        try: fld, cut, cube, blend, reach, model = load_dots(p)
        except Exception as ex: self.logmsg(f'could not load dots from {p}: {ex}'); return
        self._set_field(fld); self.cube = cube; self.sp_cube.setValue(cube); self.sp_cut.setValue(cut); self.chk_blend.setChecked(blend); self.sp_reach.setValue(max(reach, 1))
        self.logmsg(f'loaded {len(fld["idx"])} dots from {os.path.basename(p)} (cube {cube:g} mm, cut {cut:g}, {f"reach {reach}" if reach else "old straight-line fold angle"}, {"triangles blended in the cube" if blend else "nearest triangle"}): the field holds these cubes only, so the cut filters upward from {cut:g}')

    def _cut_changed(self, v):
        if abs(self.sl_cut.value()-v*10) >= 1: self.sl_cut.blockSignals(True); self.sl_cut.setValue(int(round(v*10))); self.sl_cut.blockSignals(False)
        self.cut_timer.start()

    def load_dots(self):
        os.makedirs(DOTS_DIR, exist_ok=True)
        p, _ = QtWidgets.QFileDialog.getOpenFileName(self, 'Load dots', DOTS_DIR, 'Dot sets (*.npz)')
        if not p: return
        self._use_dots(p); self.busy(False); self.redraw()

    def save_dots(self):
        if self.field is None: return
        os.makedirs(DOTS_DIR, exist_ok=True)
        p, _ = QtWidgets.QFileDialog.getSaveFileName(self, 'Save dots', DOTS_DIR, 'Dot sets (*.npz)')
        if not p: return
        try: self.logmsg(f'saved {self._ndots()} dots to {save_dots(self.field, self.sp_cut.value(), self.path, self.cube, self.chk_blend.isChecked(), self.sp_reach.value(), p)}')
        except Exception: self.logmsg(traceback.format_exc())

    def compute_field(self):
        cube, blend, reach = self.sp_cube.value(), self.chk_blend.isChecked(), self.sp_reach.value(); workers = int(self.settings['workers']) or None
        self.logmsg(f'--- field: cube={cube} mm, reach={reach}, blend={blend}'); t0 = time.time(); self.busy(True)
        self.worker = Worker(lambda log: dots.field(self.P, self.F, cube, blend, log, workers=workers, reach=reach))
        self.worker.log.connect(self.logmsg); self.worker.done.connect(lambda r: self._field_done(r, cube, blend, reach, t0)); self.worker.failed.connect(self.failed); self.worker.start()
    def failed(self, msg): self.logmsg(msg); self.busy(False)
    def _field_done(self, r, cube, blend, reach, t0):
        self._set_field(r); self.cube = cube; self.logmsg(f'field done in {time.time()-t0:.1f} s: {len(r["idx"])} cubes')
        try: self.logmsg(f'saved the dots at cut {self.sp_cut.value():g} to {save_dots(r, self.sp_cut.value(), self.path, cube, blend, reach)}')
        except Exception: self.logmsg(traceback.format_exc())
        self.busy(False); self.redraw()

    def _dot_index(self):
        """indices of the field's cubes at or above the cut: the prefix of the order"""
        k = int(np.searchsorted(self.sorted_turn, -self.sp_cut.value(), side='right')); return self.order[:k]
    def _ndots(self): return 0 if self.field is None else len(self._dot_index())

    def redraw(self, reset=False):
        if self.P is None: return
        cam = self.plotter.camera_position; self.plotter.renderer.clear_actors()
        if self.plotter.scalar_bars: self.plotter.remove_scalar_bar()
        mesh = pv.PolyData(self.P, np.hstack([np.full((len(self.F), 1), 3), self.F]).ravel())
        self.plotter.add_mesh(mesh, color='lightgrey', smooth_shading=not self.flat_shading)
        if self.field is not None:
            r = self.field; m = self._dot_index(); cut = self.sp_cut.value(); self.lbl_dots.setText(f'{len(m)} dots of {len(r["idx"])} cubes at cut {cut:g}')
            if len(m) and self.chk_dots.isChecked():
                C = r['lo']+(r['idx'][m]+0.5)*r['size']+self.settings['dotlift']*r['size']*r['N'][m]
                d = pv.PolyData(C); d['fold (deg)'] = r['turn'][m]; hi = max(float(np.percentile(d['fold (deg)'], 99)), cut+1e-3)
                self.plotter.add_mesh(d, scalars='fold (deg)', cmap='plasma', clim=(cut, hi), point_size=int(self.settings['dotsize']), render_points_as_spheres=True,
                                      scalar_bar_args=dict(title='dot strength (deg)', vertical=True, position_x=0.88, position_y=0.3, height=0.4, width=0.06))
        else: self.lbl_dots.setText('no field')
        if reset: self.plotter.reset_camera()
        else: self.plotter.camera_position = cam
        self.plotter.render()


    def dump(self):
        """out/debug/<time>/: field.npz (every cube: idx, N, turn, lo, size), settings.json (model, cube, blend, cut, camera, display), log.txt, viewport.png"""
        d = os.path.join(DEBUG_DIR, time.strftime('%Y%m%d_%H%M%S')); os.makedirs(d, exist_ok=True)
        try:
            prm = dict(model=self.path, faces=None if self.F is None else int(len(self.F)), cube_mm=self.cube, reach=self.sp_reach.value(), blend=self.chk_blend.isChecked(), cut=self.sp_cut.value(), show_dots=self.chk_dots.isChecked(),
                       flat_shading=self.flat_shading, settings=dict(self.settings), camera=[list(map(float, v)) for v in self.plotter.camera_position], cubes=0, dots=0)
            if self.field is not None:
                f = self.field; np.savez(os.path.join(d, 'field.npz'), idx=f['idx'], N=f['N'], turn=f['turn'], lo=f['lo'], size=f['size']); prm['cubes'] = int(len(f['idx'])); prm['dots'] = self._ndots()
            with open(os.path.join(d, 'log.txt'), 'w', encoding='utf-8') as fh: fh.write(self.txt.toPlainText())
            with open(os.path.join(d, 'settings.json'), 'w', encoding='utf-8') as fh: json.dump(prm, fh, indent=1)
            self.plotter.screenshot(os.path.join(d, 'viewport.png'))
            self.logmsg(f'dumped {prm["cubes"]} cubes ({prm["dots"]} dots at cut {prm["cut"]:g}), settings, log, screenshot to {d}')
        except Exception: self.logmsg(traceback.format_exc())


if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv); w = App(); w.show(); sys.exit(app.exec())
