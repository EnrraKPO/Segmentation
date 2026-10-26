"""Dots: the dots stage of the segmentation tool (built 2026-10-24 as "v2", renamed Dots the same night). Built stage by stage from what the first app and its experiments
settled. Stage 1, the field (dots.py): cubes over the surface, each with the facing of the flat triangles inside it,
laid once per cube size. Then the strength per cube: the angles between its arrow and the arrows of the surface cubes
within reach, weighted by nearness and by size (dots.py header), recomputed on the field held whenever the reach or the
difference weight changes. Stage 2, the cut:
a cube is a dot where its strength is at least the cut, a filter on the field held in memory, so sliding it redraws
the dots at once. Further filters on the same field: spread, direction, concentration (the cube's own turn against the
average turn around it, so uniform bend never floods the view) and the crest share.

The panel holds only the handles expected to change; everything else (display, cores) is in Settings. Dot sets are
saved to out/dots/ whenever Find dots finishes and the last set of a model is loaded with it. Run: python app.py."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')   # the strength runs its own threads; BLAS must not add more per call
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
    reach='How far the strength looks, in cubes: every surface cube within this distance is compared with the cube, with a weight of 1 next door fading in a straight line to nothing one step past the reach. 1 reads crisp creases only; a bend spread over several cubes reads its full angle when the reach spans it. Changing it recomputes the strength on the field held (seconds), not the cubes.',
    steep='How much a cube that differs more counts more: each compared cube is weighted by its difference to this power. 0: every cube within reach counts the same (a smooth bump reads its mean difference). 1: a cube counts as much as it differs, so the cubes across a crease decide and the ones that agree weigh nothing. Recomputes the strength on the field held.',
    scut='The spread a cube needs to be a dot: the pairwise reading, in degrees, where the region\'s cubes are compared with each other (weighted by nearness of both and by their difference), so the crest of a crease reads highest and its flanks fall off. 0 turns this filter off. Instant. When above 0 the dots are coloured by spread instead of strength.',
    crest='Crest share: keep a cube only if its reading (spread when its cut is on, else strength) is at least this share of the highest reading within reach on its across line. 100 keeps the top cube only, 80 keeps everything within 80 %% of the top, 0 turns it off. Bands become lines of a chosen width at any reach. Instant.',
    dcut='The direction a cube needs to be a dot, in degrees: how much wider the cloud of arrows within reach is along its longest direction than at right angles. A crease turns across one direction (high); a bump turns every way (near zero at any size); a plane is zero. 0 turns this filter off. Instant. When above 0 the dots are coloured by direction.',
    conc='Concentration: the turn at this cube divided by the average turn of the surface cubes within reach (fade-weighted, the cube included). 1 means it turns like its surroundings (a button, fur, a plane); a crest whose region is mostly flank reads several times. Silent on a crease whose turn is wider than the reach; dumbed down when the reach spans several features. 0 turns it off. Instant.',
    ccut='The contrast a cube needs to be a dot: its strength minus the mean strength of the surface cubes within reach, in degrees. Uniform bend (fur, a sphere) has zero contrast at any strength, so it cannot flood; a crest stands above its flanks by most of its strength. 0 turns this filter off. Instant, like the cut.',
    lcut='Line: how far the arrows within reach lie along one line on the sphere, 0 a blob (a button) to 1 a perfect line (a crease, however subtle: the size of the fold is divided out). A cylinder-shaped surface scores as a line everywhere. 0 turns it off. Instant. When above 0 the dots are coloured by line.',
    lift='Lift: strength, spread, direction and concentration each become raw x (1 + lift x line). A perfect line at lift 1 doubles its readings, a blob keeps its own, so subtle creases are helped past every cut instead of gated. 0 turns it off. Takes a second or two (the crest shares are redone).',
    cube="Edge of the cubes laid over the surface, in mm. A few times smaller than the model's triangles but not far below them (orc: 0.41 mm triangles, 0.087 mm cubes). Halving it gives four times the cubes and the time. The only limit is the grid: fewer than 65 536 cubes along the model's height.",
    cut="The strength a cube needs to be a dot: the weighted mean angle (degrees) between its arrow and the arrows of the surface cubes within reach. A cube on a crisp 90-degree edge reads near 90 at steepness 1, a smooth curve a few degrees. A filter on the field: slide it and the dots redraw. A loaded dot set holds only the cubes above its own cut, so below that nothing more appears.",
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


def save_dots(fld, cut, model, cube, blend, reach=1, steep=1.0, scut=0.0, crest_pct=0.0, dcut=0.0, conc=0.0, ccut=0.0, lcut=0.0, lift=0.0, path=None):
    """The dots alone (cubes at or above the cut): idx, N, turn, plus lo, size, cut and the settings. Same file layout as the
    first app's, so its sets load here. Default path: out/dots/<model>_<cube>mm_cut<cut>[_nearest]_<time>.npz, noted as the
    model's last set."""
    m = dots.dots(fld, cut, scut, crest_pct, dcut, conc, ccut, lcut); base = os.path.splitext(os.path.basename(model))[0]
    if path is None:
        os.makedirs(DOTS_DIR, exist_ok=True)
        name = f'{base}_{cube:g}mm_reach{reach}_cut{cut:g}{"" if blend else "_nearest"}_{time.strftime("%Y%m%d_%H%M%S")}.npz'; path = os.path.join(DOTS_DIR, name)
        lj = os.path.join(DOTS_DIR, 'last.json'); last = json.load(open(lj)) if os.path.exists(lj) else {}
        last[os.path.basename(model)] = name; json.dump(last, open(lj, 'w'), indent=1)
    np.savez(path, idx=fld['idx'][m], N=fld['N'][m], turn=fld['turn'][m], spread=fld.get('spread', fld['turn'])[m], direction=fld.get('direction', fld['turn'])[m], concentration=fld.get('concentration', np.ones(len(m), np.float32))[m], own=fld.get('own', fld['turn'])[m], line=fld.get('line', np.ones(len(m), np.float32))[m], **{k+'_raw': fld.get(k+'_raw', fld.get(k, fld['turn']))[m] for k in LIFTED}, contrast=(fld['contrast'] if fld.get('contrast') is not None else fld['turn'])[m], crest_turn=fld.get('crest_turn', np.ones(len(m), np.float32))[m], crest_spread=fld.get('crest_spread', np.ones(len(m), np.float32))[m], lo=fld['lo'], size=fld['size'], cut=cut, scut=scut, crest_pct=crest_pct, dcut=dcut, conc=conc, ccut=ccut, lcut=lcut, lift=lift, cube=cube, flat=True, blend=blend, reach=reach, steep=steep, model=model)
    return path


def load_dots(path):
    """A saved dot set back as a field holding only its dots. Returns field, its cut, cube (mm), blend, reach, steepness, model path.
    Sets from before 2026-10-24 evening have no reach: their fold angle was the straight-line reading, flagged as reach 0."""
    with np.load(path, allow_pickle=True) as z:
        fld = dict(idx=z['idx'].astype(np.int32), N=z['N'].astype(np.float32), turn=z['turn'].astype(np.float32), spread=(z['spread'] if 'spread' in z.files else z['turn']).astype(np.float32), contrast=(z['contrast'] if 'contrast' in z.files else z['turn']).astype(np.float32), line=(z['line'] if 'line' in z.files else np.ones(len(z['idx']))).astype(np.float32), direction=(z['direction'] if 'direction' in z.files else z['turn']).astype(np.float32), concentration=(z['concentration'] if 'concentration' in z.files else np.ones(len(z['idx']))).astype(np.float32), crest_turn=(z['crest_turn'] if 'crest_turn' in z.files else np.ones(len(z['idx']))).astype(np.float32), crest_spread=(z['crest_spread'] if 'crest_spread' in z.files else np.ones(len(z['idx']))).astype(np.float32), lo=z['lo'], size=float(z['size']))
        return fld, float(z['cut']), float(z['cube']), bool(z['blend']) if 'blend' in z.files else False, int(z['reach']) if 'reach' in z.files else 0, float(z['steep']) if 'steep' in z.files else 1.0, str(z['model'])   # sets before 2026-10-22 took the nearest triangle


def last_dots(model):
    lj = os.path.join(DOTS_DIR, 'last.json')
    if not os.path.exists(lj): return None
    name = json.load(open(lj)).get(os.path.basename(model)); p = name and os.path.join(DOTS_DIR, name)
    return p if p and os.path.exists(p) else None


LIFTED = ('turn', 'spread', 'direction', 'concentration')   # the readings the lift multiplies; line itself and the crest shares are relative and stay
def lifted(raw, line, lift): return (raw*(1+lift*line)).astype(np.float32) if lift > 0 else raw


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
        super().__init__(); self.setWindowTitle('Dots'); self.resize(1500, 950)
        self.P = self.F = None; self.path = None; self.field = None; self.order = None; self.cube = None; self.full = False; self.worker = None; self.reach_pending = False; self.settings = dict(SETTINGS_DEFAULTS); self.flat_shading = False; self.dlg = None
        try: self.file_log, log_path = open_session_log()
        except Exception: self.file_log, log_path = None, None
        root = QtWidgets.QWidget(); self.setCentralWidget(root); lay = QtWidgets.QHBoxLayout(root)
        self.plotter = QtInteractor(root); lay.addWidget(self.plotter.interactor, 1)
        panel = QtWidgets.QWidget(); panel.setMinimumWidth(460); panel.setMaximumWidth(560); lay.addWidget(panel, 0); side = QtWidgets.QVBoxLayout(panel); side.setContentsMargins(0, 0, 0, 0)   # the panel keeps its width whatever the window does
        self.btn_load = QtWidgets.QPushButton('Load mesh…'); self.btn_load.clicked.connect(self.load); side.addWidget(self.btn_load)
        self.lbl_file = QtWidgets.QLabel('no model'); self.lbl_file.setWordWrap(True); side.addWidget(self.lbl_file)

        def spin(v, lo, hi, step, dec=0):
            s = QtWidgets.QDoubleSpinBox(); s.setRange(lo, hi); s.setSingleStep(step); s.setDecimals(dec); s.setValue(v); s.setMinimumWidth(110); return s
        def row(form, label, widget, tip):
            widget.setToolTip(tip); form.addRow(label, widget); form.labelForField(widget).setToolTip(tip)

        g1 = QtWidgets.QGroupBox('1. Field'); f1 = QtWidgets.QFormLayout(g1)
        self.sp_cube = spin(0.087, 0.0001, 5, 0.001, 4); row(f1, 'Cube size (mm)', self.sp_cube, TIPS['cube'])
        self.sp_reach = QtWidgets.QSpinBox(); self.sp_reach.setRange(1, 50); self.sp_reach.setValue(1); row(f1, 'Reach (cubes)', self.sp_reach, TIPS['reach'])
        self.sp_steep = QtWidgets.QDoubleSpinBox(); self.sp_steep.setRange(0, 4); self.sp_steep.setSingleStep(0.5); self.sp_steep.setDecimals(1); self.sp_steep.setValue(1); row(f1, 'Difference weight', self.sp_steep, TIPS['steep'])
        self.reach_timer = QtCore.QTimer(self); self.reach_timer.setSingleShot(True); self.reach_timer.setInterval(300); self.reach_timer.timeout.connect(self.recompute_strength); self.sp_reach.valueChanged.connect(lambda v: self.reach_timer.start()); self.sp_steep.valueChanged.connect(lambda v: self.reach_timer.start())
        self.chk_blend = QtWidgets.QCheckBox('Blend the triangles inside a cube by area'); self.chk_blend.setChecked(True); self.chk_blend.setToolTip(TIPS['blend']); f1.addRow(self.chk_blend)
        self.btn_field = QtWidgets.QPushButton('Compute the field'); self.btn_field.clicked.connect(self.compute_field); self.btn_field.setToolTip('Cubes over the surface with their arrows, then the strength at the reach. Reach and cut can be changed afterwards without laying the cubes again.'); f1.addRow(self.btn_field)
        self.btn_load_dots = QtWidgets.QPushButton('Load dots…'); self.btn_load_dots.clicked.connect(self.load_dots); self.btn_load_dots.setToolTip('A saved set, from this app or the first one. Opening a mesh loads its last set by itself.'); f1.addRow(self.btn_load_dots)
        side.addWidget(g1)

        g2 = QtWidgets.QGroupBox('2. Cut'); f2 = QtWidgets.QFormLayout(g2)
        cw = QtWidgets.QWidget(); ch = QtWidgets.QHBoxLayout(cw); ch.setContentsMargins(0, 0, 0, 0)
        self.sl_cut = QtWidgets.QSlider(QtCore.Qt.Horizontal); self.sl_cut.setRange(10, 1800); self.sl_cut.setValue(170); ch.addWidget(self.sl_cut, 1)   # tenths of a degree
        self.sp_cut = spin(17, 1, 180, 1, 1); self.sp_cut.setMinimumWidth(96); ch.addWidget(self.sp_cut); row(f2, 'Cut (deg)', cw, TIPS['cut']); self.sp_cut.setToolTip(TIPS['cut'])
        self.sl_cut.valueChanged.connect(lambda v: self.sp_cut.setValue(v/10)); self.sp_cut.valueChanged.connect(self._cut_changed)
        kw = QtWidgets.QWidget(); kh = QtWidgets.QHBoxLayout(kw); kh.setContentsMargins(0, 0, 0, 0)
        self.sl_ccut = QtWidgets.QSlider(QtCore.Qt.Horizontal); self.sl_ccut.setRange(0, 900); self.sl_ccut.setValue(0); kh.addWidget(self.sl_ccut, 1)   # tenths of a degree
        self.sp_ccut = spin(0, 0, 90, 1, 1); kh.addWidget(self.sp_ccut); row(f2, 'Contrast (deg)', kw, TIPS['ccut']); self.sp_ccut.setToolTip(TIPS['ccut'])
        self.sl_ccut.valueChanged.connect(lambda v: self.sp_ccut.setValue(v/10)); self.sp_ccut.valueChanged.connect(self._ccut_changed)
        pw = QtWidgets.QWidget(); ph = QtWidgets.QHBoxLayout(pw); ph.setContentsMargins(0, 0, 0, 0)
        self.sl_scut = QtWidgets.QSlider(QtCore.Qt.Horizontal); self.sl_scut.setRange(0, 1800); self.sl_scut.setValue(0); ph.addWidget(self.sl_scut, 1)   # tenths of a degree
        self.sp_scut = spin(0, 0, 180, 1, 1); self.sp_scut.setMinimumWidth(96); ph.addWidget(self.sp_scut); row(f2, 'Spread (deg)', pw, TIPS['scut']); self.sp_scut.setToolTip(TIPS['scut'])
        self.sl_scut.valueChanged.connect(lambda v: self.sp_scut.setValue(v/10)); self.sp_scut.valueChanged.connect(self._scut_changed)
        dw = QtWidgets.QWidget(); dh = QtWidgets.QHBoxLayout(dw); dh.setContentsMargins(0, 0, 0, 0)
        self.sl_dcut = QtWidgets.QSlider(QtCore.Qt.Horizontal); self.sl_dcut.setRange(0, 900); self.sl_dcut.setValue(0); dh.addWidget(self.sl_dcut, 1)   # tenths of a degree
        self.sp_dcut = spin(0, 0, 90, 1, 1); self.sp_dcut.setMinimumWidth(96); dh.addWidget(self.sp_dcut); row(f2, 'Direction (deg)', dw, TIPS['dcut']); self.sp_dcut.setToolTip(TIPS['dcut'])
        self.sl_dcut.valueChanged.connect(lambda v: self.sp_dcut.setValue(v/10)); self.sp_dcut.valueChanged.connect(self._dcut_changed)
        lw = QtWidgets.QWidget(); lh = QtWidgets.QHBoxLayout(lw); lh.setContentsMargins(0, 0, 0, 0)
        self.sl_lcut = QtWidgets.QSlider(QtCore.Qt.Horizontal); self.sl_lcut.setRange(0, 100); self.sl_lcut.setValue(0); lh.addWidget(self.sl_lcut, 1)   # hundredths
        self.sp_lcut = spin(0, 0, 1, 0.01, 2); lh.addWidget(self.sp_lcut); row(f2, 'Line (0-1)', lw, TIPS['lcut']); self.sp_lcut.setToolTip(TIPS['lcut'])
        self.sl_lcut.valueChanged.connect(lambda v: self.sp_lcut.setValue(v/100)); self.sp_lcut.valueChanged.connect(self._lcut_changed)
        self.sp_lift = spin(0, 0, 10, 0.1, 1); row(f2, 'Lift', self.sp_lift, TIPS['lift'])
        self.lift_timer = QtCore.QTimer(self); self.lift_timer.setSingleShot(True); self.lift_timer.setInterval(300); self.lift_timer.timeout.connect(self.apply_lift); self.sp_lift.valueChanged.connect(lambda v: self.lift_timer.start())
        qw = QtWidgets.QWidget(); qh = QtWidgets.QHBoxLayout(qw); qh.setContentsMargins(0, 0, 0, 0)
        self.sl_conc = QtWidgets.QSlider(QtCore.Qt.Horizontal); self.sl_conc.setRange(0, 200); self.sl_conc.setValue(0); qh.addWidget(self.sl_conc, 1)   # tenths of the ratio
        self.sp_conc = spin(0, 0, 20, 0.1, 1); qh.addWidget(self.sp_conc); row(f2, 'Concentration (x)', qw, TIPS['conc']); self.sp_conc.setToolTip(TIPS['conc'])
        self.sl_conc.valueChanged.connect(lambda v: self.sp_conc.setValue(v/10)); self.sp_conc.valueChanged.connect(self._conc_changed)
        xw = QtWidgets.QWidget(); xh = QtWidgets.QHBoxLayout(xw); xh.setContentsMargins(0, 0, 0, 0)
        self.sl_crest = QtWidgets.QSlider(QtCore.Qt.Horizontal); self.sl_crest.setRange(0, 100); self.sl_crest.setValue(0); xh.addWidget(self.sl_crest, 1)
        self.sp_crest = spin(0, 0, 100, 1, 0); self.sp_crest.setMinimumWidth(96); xh.addWidget(self.sp_crest); row(f2, 'Crest (%)', xw, TIPS['crest']); self.sp_crest.setToolTip(TIPS['crest'])
        self.sl_crest.valueChanged.connect(lambda v: self.sp_crest.setValue(v)); self.sp_crest.valueChanged.connect(self._crest_changed)
        self.lbl_dots = QtWidgets.QLabel('no field'); f2.addRow(self.lbl_dots)
        self.btn_save = QtWidgets.QPushButton('Save dots…'); self.btn_save.clicked.connect(self.save_dots); self.btn_save.setToolTip('Writes the dots at this cut where you choose. Every computed field is also saved at its cut to out/dots/ by itself.'); f2.addRow(self.btn_save)
        side.addWidget(g2)
        self.cut_timer = QtCore.QTimer(self); self.cut_timer.setSingleShot(True); self.cut_timer.setInterval(60); self.cut_timer.timeout.connect(self.redraw)   # one redraw per slider pause, not per tick

        self.chk_dots = QtWidgets.QCheckBox('Show dots'); self.chk_dots.setChecked(True); self.chk_dots.stateChanged.connect(self.redraw); side.addWidget(self.chk_dots)
        self.btn_dump = QtWidgets.QPushButton('Dump for debugging'); self.btn_dump.clicked.connect(self.dump); self.btn_dump.setToolTip('Writes the whole field (every cube, its arrow and strength), the settings, the camera, the log and a screenshot to out/debug/<time>/.'); side.addWidget(self.btn_dump)
        self.btn_settings = QtWidgets.QPushButton('Settings…'); self.btn_settings.clicked.connect(self.open_settings); self.btn_settings.setToolTip('Display and machine settings: nothing here changes a dot.'); side.addWidget(self.btn_settings)
        self.txt = QtWidgets.QPlainTextEdit(); self.txt.setReadOnly(True); self.txt.setMinimumHeight(160); side.addWidget(self.txt, 1)
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
        """the field held (idx, N and the readings); the cuts are masks over it, instant"""
        self.field = fld

    def _use_dots(self, p):
        try: fld, cut, cube, blend, reach, steep, model = load_dots(p)
        except Exception as ex: self.logmsg(f'could not load dots from {p}: {ex}'); return
        self._set_field(fld); self.full = False; self.cube = cube; self.sp_cube.setValue(cube); self.sp_cut.setValue(cut); self.chk_blend.setChecked(blend); self.sp_reach.blockSignals(True); self.sp_reach.setValue(max(reach, 1)); self.sp_reach.blockSignals(False); self.sp_steep.blockSignals(True); self.sp_steep.setValue(steep); self.sp_steep.blockSignals(False)
        self.logmsg(f'loaded {len(fld["idx"])} dots from {os.path.basename(p)} (cube {cube:g} mm, cut {cut:g}, {f"reach {reach}, steepness {steep:g}" if reach else "old fold angle"}, {"triangles blended in the cube" if blend else "nearest triangle"}): the field holds these cubes only, so the cut filters upward from {cut:g}')

    def _cut_changed(self, v):
        if abs(self.sl_cut.value()-v*10) >= 1: self.sl_cut.blockSignals(True); self.sl_cut.setValue(int(round(v*10))); self.sl_cut.blockSignals(False)
        self.cut_timer.start()
    def _ccut_changed(self, v):
        if abs(self.sl_ccut.value()-v*10) >= 1: self.sl_ccut.blockSignals(True); self.sl_ccut.setValue(int(round(v*10))); self.sl_ccut.blockSignals(False)
        if v > 0 and self.field is not None and self.full and self.field.get('contrast') is None: self.compute_contrast(); return
        self.cut_timer.start()
    def compute_contrast(self):
        """contrast is costly (the ring beyond the reach), so it is computed only once its cut is raised above zero"""
        if self.worker is not None and self.worker.isRunning(): self.reach_pending = True; return
        fld = self.field; reach = self.sp_reach.value(); self.busy(True)
        self.worker = Worker(lambda log: dots.contrast(fld, fld['turn'], reach, log))
        def done(c): fld['contrast'] = c; self.busy(False); self.redraw()
        self.worker.log.connect(self.logmsg); self.worker.done.connect(done); self.worker.failed.connect(self.failed); self.worker.start()
    def _scut_changed(self, v):
        if abs(self.sl_scut.value()-v*10) >= 1: self.sl_scut.blockSignals(True); self.sl_scut.setValue(int(round(v*10))); self.sl_scut.blockSignals(False)
        self.cut_timer.start()
    def _dcut_changed(self, v):
        if abs(self.sl_dcut.value()-v*10) >= 1: self.sl_dcut.blockSignals(True); self.sl_dcut.setValue(int(round(v*10))); self.sl_dcut.blockSignals(False)
        self.cut_timer.start()
    def _lcut_changed(self, v):
        if abs(self.sl_lcut.value()-v*100) >= 1: self.sl_lcut.blockSignals(True); self.sl_lcut.setValue(int(round(v*100))); self.sl_lcut.blockSignals(False)
        self.cut_timer.start()
    def _conc_changed(self, v):
        if abs(self.sl_conc.value()-v*10) >= 1: self.sl_conc.blockSignals(True); self.sl_conc.setValue(int(round(v*10))); self.sl_conc.blockSignals(False)
        self.cut_timer.start()
    def _crest_changed(self, v):
        if abs(self.sl_crest.value()-v) >= 1: self.sl_crest.blockSignals(True); self.sl_crest.setValue(int(round(v))); self.sl_crest.blockSignals(False)
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
        try: self.logmsg(f'saved {self._ndots()} dots to {save_dots(self.field, self.sp_cut.value(), self.path, self.cube, self.chk_blend.isChecked(), self.sp_reach.value(), self.sp_steep.value(), self.sp_scut.value(), self.sp_crest.value(), self.sp_dcut.value(), self.sp_conc.value(), self.sp_ccut.value(), self.sp_lcut.value(), self.sp_lift.value(), p)}')
        except Exception: self.logmsg(traceback.format_exc())

    def compute_field(self):
        cube, blend, reach, steep = self.sp_cube.value(), self.chk_blend.isChecked(), self.sp_reach.value(), self.sp_steep.value(); workers = int(self.settings['workers']) or None
        self.logmsg(f'--- field: cube={cube} mm, blend={blend}; then the strength at reach {reach}, steepness {steep:g}'); t0 = time.time(); self.busy(True)
        def job(log):
            r = dots.field(self.P, self.F, cube, blend, log, workers=workers); r.update(self._readings(r, reach, steep, log)); return r
        self.worker = Worker(job)
        self.worker.log.connect(self.logmsg); self.worker.done.connect(lambda r: self._field_done(r, cube, blend, reach, steep, t0)); self.worker.failed.connect(self.failed); self.worker.start()
    def failed(self, msg): self.logmsg(msg); self.busy(False)
    def _field_done(self, r, cube, blend, reach, steep, t0):
        self._set_field(r); self.full = True; self.cube = cube; self.logmsg(f'field done in {time.time()-t0:.1f} s: {len(r["idx"])} cubes')
        try: self.logmsg(f'saved the dots at cut {self.sp_cut.value():g} to {save_dots(r, self.sp_cut.value(), self.path, cube, blend, reach, steep, self.sp_scut.value(), self.sp_crest.value(), self.sp_dcut.value(), self.sp_conc.value(), self.sp_ccut.value(), self.sp_lcut.value(), self.sp_lift.value())}')
        except Exception: self.logmsg(traceback.format_exc())
        self.busy(False); self.redraw()
        if self.reach_pending: self.reach_pending = False; self.recompute_strength()

    def _readings(self, fld, reach, steep, log):
        """everything that depends on reach and weight: strength, spread, direction, line, concentration, across direction, crest shares; the lift applied to strength"""
        r = dots.readings(fld, reach, steep, log); lift = self.sp_lift.value()
        for key in LIFTED: r[key+'_raw'] = r[key]; r[key] = lifted(r[key], r['line'], lift)
        if lift > 0: log(f'lift {lift:g}: strength, spread, direction and concentration x (1 + {lift:g} x line)')
        r['crest_turn'] = dots.crest_share(fld, r['turn'], r['across'], r['far'], reach, log); r['crest_spread'] = dots.crest_share(fld, r['spread'], r['across'], r['far'], reach, log)
        r['contrast'] = dots.contrast(fld, r['turn'], reach, log) if self.sp_ccut.value() > 0 else None
        return r

    def apply_lift(self):
        """the lift changed: every size-carrying reading is rebuilt from its raw one and the line, the crest shares redone; the cuts then filter as before"""
        fld = self.field
        if fld is None or fld.get('turn_raw') is None or fld.get('line') is None: return
        if not self.full: self.logmsg('the lift applies to a computed field; the loaded set holds its dots only'); return
        if self.worker is not None and self.worker.isRunning(): self.reach_pending = True; return
        lift = self.sp_lift.value(); reach = self.sp_reach.value(); self.busy(True)
        def job(log):
            out = {key: lifted(fld[key+'_raw'], fld['line'], lift) for key in LIFTED if fld.get(key+'_raw') is not None}
            log(f'lift {lift:g}: strength, spread, direction and concentration x (1 + {lift:g} x line)')
            out['crest_turn'] = dots.crest_share(fld, out['turn'], fld['across'], fld['far'], reach, log); out['crest_spread'] = dots.crest_share(fld, out['spread'], fld['across'], fld['far'], reach, log)
            if self.sp_ccut.value() > 0: out['contrast'] = dots.contrast(fld, out['turn'], reach, log)
            return out
        self.worker = Worker(job)
        def done(res):
            fld.update(res); self._set_field(fld); self.busy(False); self.redraw()
            if self.reach_pending: self.reach_pending = False; self.recompute_strength()
        self.worker.log.connect(self.logmsg); self.worker.done.connect(done); self.worker.failed.connect(self.failed); self.worker.start()

    def recompute_strength(self):
        """the reach or the steepness changed: the strength is read again on the field held (not the cubes); the cut then filters as before"""
        if self.field is None: return
        if not self.full: self.logmsg('the reach applies to a computed field; the loaded set holds its dots only, compute the field to change the reach'); return
        if self.worker is not None and self.worker.isRunning(): self.reach_pending = True; return
        reach, steep = self.sp_reach.value(), self.sp_steep.value(); fld = self.field; t0 = time.time(); self.busy(True)
        self.worker = Worker(lambda log: self._readings(fld, reach, steep, log))
        def done(res):
            fld.update(res); self._set_field(fld); self.busy(False); self.redraw()
            if self.reach_pending: self.reach_pending = False; self.recompute_strength()
        self.worker.log.connect(self.logmsg); self.worker.done.connect(done); self.worker.failed.connect(self.failed); self.worker.start()

    def _dot_index(self):
        """indices of the field's cubes passing every cut"""
        return np.nonzero(dots.dots(self.field, self.sp_cut.value(), self.sp_scut.value(), self.sp_crest.value(), self.sp_dcut.value(), self.sp_conc.value(), self.sp_ccut.value(), self.sp_lcut.value()))[0]
    def _ndots(self): return 0 if self.field is None else len(self._dot_index())

    def redraw(self, reset=False):
        if self.P is None: return
        cam = self.plotter.camera_position; self.plotter.renderer.clear_actors()
        if self.plotter.scalar_bars: self.plotter.remove_scalar_bar()
        mesh = pv.PolyData(self.P, np.hstack([np.full((len(self.F), 1), 3), self.F]).ravel())
        self.plotter.add_mesh(mesh, color='lightgrey', smooth_shading=not self.flat_shading)
        if self.field is not None:
            r = self.field; m = self._dot_index(); cut = self.sp_cut.value(); scut = self.sp_scut.value(); self.lbl_dots.setText(f'{len(m)} dots of {len(r["idx"])} cubes at cut {cut:g}' + (f', contrast {self.sp_ccut.value():g}' if self.sp_ccut.value() > 0 else '') + (f', spread {scut:g}' if scut > 0 else '') + (f', direction {self.sp_dcut.value():g}' if self.sp_dcut.value() > 0 else '') + (f', concentration x{self.sp_conc.value():g}' if self.sp_conc.value() > 0 else '') + (f', line {self.sp_lcut.value():g}' if self.sp_lcut.value() > 0 else '') + (f', crest {self.sp_crest.value():g} %' if self.sp_crest.value() > 0 else ''))
            if len(m) and self.chk_dots.isChecked():
                C = r['lo']+(r['idx'][m]+0.5)*r['size']+self.settings['dotlift']*r['size']*r['N'][m]
                by, lowv = ('line', self.sp_lcut.value()) if self.sp_lcut.value() > 0 and r.get('line') is not None else ('direction', self.sp_dcut.value()) if self.sp_dcut.value() > 0 and r.get('direction') is not None else ('spread', scut) if scut > 0 and 'spread' in r else ('turn', cut)
                d = pv.PolyData(C); d['strength'] = r[by][m]; hi = max(float(np.percentile(d['strength'], 99)), lowv+1e-3)
                self.plotter.add_mesh(d, scalars='strength', cmap='plasma', clim=(lowv, hi), point_size=int(self.settings['dotsize']), render_points_as_spheres=True,
                                      scalar_bar_args=dict(title='strength', vertical=True, position_x=0.88, position_y=0.3, height=0.4, width=0.06))
        else: self.lbl_dots.setText('no field')
        if reset: self.plotter.reset_camera()
        else: self.plotter.camera_position = cam
        self.plotter.render()


    def dump(self):
        """out/debug/<time>/: field.npz (every cube: idx, N, turn, lo, size), settings.json (model, cube, blend, cut, camera, display), log.txt, viewport.png"""
        d = os.path.join(DEBUG_DIR, time.strftime('%Y%m%d_%H%M%S')); os.makedirs(d, exist_ok=True)
        try:
            prm = dict(model=self.path, faces=None if self.F is None else int(len(self.F)), cube_mm=self.cube, reach=self.sp_reach.value(), steep=self.sp_steep.value(), lift=self.sp_lift.value(), blend=self.chk_blend.isChecked(), cut=self.sp_cut.value(), scut=self.sp_scut.value(), crest_pct=self.sp_crest.value(), dcut=self.sp_dcut.value(), conc=self.sp_conc.value(), ccut=self.sp_ccut.value(), lcut=self.sp_lcut.value(), show_dots=self.chk_dots.isChecked(),
                       flat_shading=self.flat_shading, settings=dict(self.settings), camera=[list(map(float, v)) for v in self.plotter.camera_position], cubes=0, dots=0)
            if self.field is not None:
                f = self.field; np.savez(os.path.join(d, 'field.npz'), idx=f['idx'], N=f['N'], turn=f['turn'], spread=f.get('spread', f['turn']), direction=f.get('direction', f['turn']), concentration=f.get('concentration', np.ones(len(f['idx']), np.float32)), crest_turn=f.get('crest_turn', np.ones(len(f['idx']), np.float32)), crest_spread=f.get('crest_spread', np.ones(len(f['idx']), np.float32)), own=f.get('own', f['turn']), **{k+'_raw': f.get(k+'_raw', f.get(k, f['turn'])) for k in LIFTED}, line=f.get('line', np.ones(len(f['idx']), np.float32)), contrast=(f['contrast'] if f.get('contrast') is not None else f['turn']), lo=f['lo'], size=f['size']); prm['cubes'] = int(len(f['idx'])); prm['dots'] = self._ndots()
            with open(os.path.join(d, 'log.txt'), 'w', encoding='utf-8') as fh: fh.write(self.txt.toPlainText())
            with open(os.path.join(d, 'settings.json'), 'w', encoding='utf-8') as fh: json.dump(prm, fh, indent=1)
            self.plotter.screenshot(os.path.join(d, 'viewport.png'))
            self.logmsg(f'dumped {prm["cubes"]} cubes ({prm["dots"]} dots at cut {prm["cut"]:g}), settings, log, screenshot to {d}')
        except Exception: self.logmsg(traceback.format_exc())


if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv); w = App(); w.show(); sys.exit(app.exec())
