"""Lines: the lines stage of the segmentation tool, a clean slate started the night of 2026-10-24.

It consumes what the Dots tool (../Dots) saves: a dot set in out/dots (the cubes above the cuts, with their arrows and
readings), and the mesh the set was made from. The walk is the painting walk (walk.py, designed 2026-10-25): paint
outward from a seed over touching dots, one ring per round; the ring splits in two where the band's sides stop, and
from there each round's new dots is a slice across the band whose strength-weighted centre is a node. It halts at
anything confusing (a branch, no dots, another line).

The node read (walk.py header): the crease direction from the arrows of the rows around the node, the node across at
those rows' centre, along at its own row, then onto the surface.

The rig: every read is stepped and watched. Dots are coloured by strength. A seed is the strongest free dot or the dot
clicked in the viewport. Next read = one read (spread, clump, node); Next node = reads until a node; Run = to the halt.
Dots are a quiet greyscale by strength (toggleable). The viewport shows the line's paint (tan), the front (white),
the new dots by clump (one colour each), the seed (red), the direction as an arrow, the nodes (blue, newest large),
finished lines (navy), with a legend. The log prints every read's numbers; Dump writes it all.

Run: python app.py (or run.bat)."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import sys, time, traceback, json
import numpy as np, trimesh
from PySide6 import QtCore, QtWidgets
import pyvista as pv
from pyvistaqt import QtInteractor

ROOT = os.path.dirname(os.path.abspath(__file__))
DOTS_ROOT = os.path.join(os.path.dirname(ROOT), 'Dots')
import importlib.util
sys.path.insert(0, DOTS_ROOT)   # so the Dots app finds its dots.py
sys.path.insert(0, ROOT)
_spec = importlib.util.spec_from_file_location('dots_app', os.path.join(DOTS_ROOT, 'app.py')); _dots_app = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(_dots_app)
load_mesh, load_dots, last_dots, DOTS_DIR = _dots_app.load_mesh, _dots_app.load_dots, _dots_app.last_dots, _dots_app.DOTS_DIR   # the Dots tool's loaders: one definition of the file format
import walk as walkmod

LOG_DIR = os.path.join(ROOT, 'out', 'log')
DEBUG_DIR = os.path.join(ROOT, 'out', 'debug')
SETTINGS = [   # key, label, default, lo, hi, step, decimals, tip
    ('dotlift', 'Dot lift (cubes)', 0.5, 0, 5, 0.1, 1, 'Dots are drawn this many cubes out along their arrow.'),
    ('dotsize', 'Dot size (px)', 5, 1, 20, 1, 0, 'Point size of the drawn dots.'),
    ('linelift', 'Line lift (cubes)', 0, 0, 5, 0.1, 1, 'Lines and nodes are drawn this many cubes above the surface. 0: on the surface, where the nodes are (the dots, drawn lifted, may then hide them). About dot lift + 0.5 puts them above the dots.'),
    ('picktol', 'Pick tolerance (window fraction)', 0.01, 0.001, 0.1, 0.001, 3, 'A dot within this fraction of the window of the cursor is picked.'),
    ('reach', 'Paint reach (cubes)', 1, 1, 10, 1, 0, 'Two dots touch when their cubes lie within this many cubes of each other. 1: next door only, one missing dot splits the band. 2: one empty cube is crossed, so grid noise does not split the band; each round then paints two cubes deep. Changing it rebuilds the neighbour table on the next Start walk.'),
    ('window', 'Node window (rows)', 3, 0, 50, 1, 0, 'A node reads the crease direction and its across position from the dots of this many rows before and after its own row (2 x window + 1 rows). The paint runs this many rows ahead of the nodes. 0: the row alone.'),
    ('gap', 'Gap radius (mm)', 0.0, 0, 10, 0.001, 4, 'At a dead end, the unpainted dots within this distance of the front that lie ahead along the arm become the next row, so a hole in a narrow band is crossed. 0 turns it off.'),
    ('fork', 'Faulty row (deg)', 20, 1, 90, 1, 0, "A row is faulty when its own arrows are more than this many degrees off right angles to the stem's line (a stem reads 2-8, a fork 15-50). Faulty rows ahead are dropped from a node's window; a node whose own row is faulty cuts the stem: the branches open from the node before."),
    ('kink', 'Loose branch (deg)', 30, 1, 90, 1, 0, "A branch's first segment, from the junction to its first node, must run along the stem's line or along the branch's own direction. Off both by more than this many degrees it follows no crease: the branch is cut loose from the junction and walked on as a line of its own (a real fork reads 1-15 off its own line, a false one 40-90)."),
    ('lead', 'Second crest (share of the crest)', 0.0, 0, 1, 0.05, 2, "A row is read across, every dot at full strength: from its strongest crest outward, a rise after a dip is a second crest. A rise over this share of the crest's strength makes the row two rows: the arm keeps the crest its last row sits on, the other crest becomes a second band under the same paint, an arm of its own from a junction on the stem. 0: any dip counts."),
    ('cap', 'Finding cap (rounds)', 30, 2, 1000, 1, 0, 'The finding stage halts when the ring has not split in two within this many rounds. A cap for Run, not a rule.'),
    ('maxnodes', 'Run cap (nodes)', 100000, 10, 10000000, 1000, 0, 'Run stops after this many nodes.'),
]
SETTINGS_DEFAULTS = {k: v for k, _, v, *_ in SETTINGS}
TIPS = dict(
    seed='Strongest free dot: the walk starts at the strongest dot no line has painted. Picked dot: the dot you clicked in the viewport (click a dot without dragging).',
    start='Begins a line at the seed: its dot is painted and the finding stage is set to spread. Nothing is read yet.',
    read='One read. Spread: paint every unpainted dot touching the front. Clump: group the new dots by touching. Direction: the crease direction from the arrows of the rows around the node. Node: decide (ring, line found, branch, end) and place the node: across at the centre of those rows, along at its own row.',
    node='Reads until a node is placed or the walk halts.',
    run='Reads until the arm under way halts (its end, a junction, a collision) or a split comes under test, then stops so you can look.',
    runall='Reads until the whole walk is done: every arm, every junction.',
    runlines='Walks every line in the set: the strongest free dot seeds a line, it is run to its end, and so on until no free dot is left. Every read goes to the session log file; the box gets one line per finished line.',
    find='Type a name and press Enter to select it and move the camera there. A dot: "dot 53873" (or just 53873), its number in the loaded set, as the log names it. A node: "L2 a5 n7" (line 2, arm 5, node 7, numbered as the log does), "L2 seed" or "L2 a5 start".',
    paint='Draws the dots the lines have painted: khaki for earlier lines, tan for the line under way. Off: only the current read (front, clumps, window) stays on the dots.',
    reset='Drops the walk under way and frees its dots. Finished lines stay.',
    clear='Frees every dot and drops every line.',
)


LAST_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'last.json')   # the model and dot set last in use

def open_session_log():
    os.makedirs(LOG_DIR, exist_ok=True); path = os.path.join(LOG_DIR, time.strftime('lines_%Y%m%d_%H%M%S.txt')); fh = open(path, 'a', encoding='utf-8')
    def write(s): fh.write(time.strftime('%H:%M:%S ')+s+'\n'); fh.flush()
    return write, path


def _closest_on_triangles(Q, A, B, C):
    """closest point on triangle ABC to Q, per row (Ericson 5.1.5)"""
    ab = B-A; ac = C-A; ap = Q-A; d1 = (ab*ap).sum(1); d2 = (ac*ap).sum(1); bp = Q-B; d3 = (ab*bp).sum(1); d4 = (ac*bp).sum(1); cp = Q-C; d5 = (ab*cp).sum(1); d6 = (ac*cp).sum(1)
    vc = d1*d4-d3*d2; vb = d5*d2-d1*d6; va = d3*d6-d5*d4; out = np.empty_like(Q); done = np.zeros(len(Q), bool)
    def put(m, pts): m = m & ~done; out[m] = pts[m] if pts.ndim == 2 else pts; done[m] = True
    put((d1 <= 0) & (d2 <= 0), A); put((d3 >= 0) & (d4 <= d3), B); put((d6 >= 0) & (d5 <= d6), C)
    v = d1/np.where(d1-d3 != 0, d1-d3, 1); put((vc <= 0) & (d1 >= 0) & (d3 <= 0), A+v[:, None]*ab)
    w = d2/np.where(d2-d6 != 0, d2-d6, 1); put((vb <= 0) & (d2 >= 0) & (d6 <= 0), A+w[:, None]*ac)
    den = (d4-d3)+(d5-d6); w = (d4-d3)/np.where(den != 0, den, 1); put((va <= 0) & (d4-d3 >= 0) & (d5-d6 >= 0), B+w[:, None]*(C-B))
    s_ = va+vb+vc; s_ = np.where(s_ != 0, s_, 1); v = vb/s_; w = vc/s_; put(~done, A+v[:, None]*ab+w[:, None]*ac); return out


def make_snap(P, F, k=16):
    """(point, facing) -> the surface point along the facing (the nearer hit of the ray both ways), else the nearest surface
    point. Local: only the triangles at the k nearest vertices are looked at (a point a cube off the surface is always
    nearest to the vertices of the triangles under it), so a node costs microseconds, not a cast over the whole mesh."""
    from scipy.spatial import cKDTree
    P = np.asarray(P, float); F = np.asarray(F); tree = cKDTree(P)
    fn = np.cross(P[F[:, 1]]-P[F[:, 0]], P[F[:, 2]]-P[F[:, 0]]); vn = np.zeros_like(P)
    for c_ in range(3): np.add.at(vn, F[:, c_], fn)
    vn /= np.linalg.norm(vn, axis=1, keepdims=True)+1e-12
    order = np.argsort(F.ravel(), kind='stable'); vf = F.ravel()[order]; ff = np.repeat(np.arange(len(F)), 3)[order]; starts = np.searchsorted(vf, np.arange(len(P)+1))   # faces per vertex
    def snap(p, n):
        p = np.asarray(p, float); n = np.asarray(n, float); n = n/(np.linalg.norm(n)+1e-12)
        vs = tree.query(p, k)[1]; fs = np.unique(np.concatenate([ff[starts[v]:starts[v+1]] for v in vs]))
        A, B, C = P[F[fs, 0]], P[F[fs, 1]], P[F[fs, 2]]
        e1 = B-A; e2 = C-A; h = np.cross(n, e2); a = (e1*h).sum(1); ok = np.abs(a) > 1e-12   # Moller-Trumbore, the ray along n in both senses
        f = 1/np.where(ok, a, 1); sv = p-A; u = f*(sv*h).sum(1); q = np.cross(sv, e1); v = f*(n*q).sum(1); t = f*(e2*q).sum(1)
        hit = ok & (u >= -1e-9) & (v >= -1e-9) & (u+v <= 1+1e-9)
        if hit.any(): th = t[hit]; return p+th[np.argmin(np.abs(th))]*n
        X = _closest_on_triangles(np.repeat(p[None], len(fs), 0), A, B, C); return X[np.argmin(np.linalg.norm(X-p, axis=1))]
    def normals(pts, k=64):
        """a smooth surface normal at points, for lifting drawn lines: the mean normal of the k nearest vertices. The nearest
        vertex alone flips from one side of a crease to the other from node to node and drew a sawtooth (turn per node 7 -> 33 deg)"""
        n = vn[tree.query(np.asarray(pts, float).reshape(-1, 3), k)[1]].mean(1); return n/(np.linalg.norm(n, axis=1, keepdims=True)+1e-12)
    snap.normals = normals
    return snap


class SettingsDialog(QtWidgets.QDialog):
    def __init__(self, app):
        super().__init__(app); self.app = app; self.setWindowTitle('Settings'); self.resize(460, 300)
        lay = QtWidgets.QVBoxLayout(self); form = QtWidgets.QFormLayout(); lay.addLayout(form); self.spins = {}
        for key, label, default, lo, hi, step, dec, tip in SETTINGS:
            sp = QtWidgets.QDoubleSpinBox(); sp.setRange(lo, hi); sp.setSingleStep(step); sp.setDecimals(dec); sp.setValue(app.settings[key]); sp.setToolTip(tip); sp.setMinimumWidth(110)
            sp.valueChanged.connect(lambda v, k=key: self._set(k, v)); form.addRow(label, sp); form.labelForField(sp).setToolTip(tip); self.spins[key] = sp
        self.chk_flat = QtWidgets.QCheckBox('Flat shading'); self.chk_flat.setChecked(app.flat_shading); self.chk_flat.stateChanged.connect(lambda s: self._shade(bool(s))); lay.addWidget(self.chk_flat)
        lay.addStretch(1); c = QtWidgets.QPushButton('Close'); c.clicked.connect(self.close); lay.addWidget(c)
    def _set(self, k, v):
        self.app.settings[k] = v
        if k == 'cap' and self.app.walk is not None: self.app.walk.cap = int(v)
        if k == 'reach' and self.app.walk is not None: self.app.logmsg(f'reach {int(v)}: takes effect at the next Start walk')
        if k == 'window' and self.app.walk is not None: self.app.walk.window = int(v); self.app.logmsg(f'window {int(v)}: applies from the next node')
        if k == 'gap' and self.app.walk is not None: self.app.walk.gap = float(v); self.app.logmsg(f'gap radius {v:g} mm: applies from the next dead end')
        if k == 'fork' and self.app.walk is not None: self.app.walk.fork = float(v); self.app.logmsg(f'faulty row {v:g} deg: applies from the next node')
        if k == 'kink' and self.app.walk is not None: self.app.walk.kink = float(v); self.app.logmsg(f'loose branch {v:g} deg: applies from the next branch')
        if k == 'lead' and self.app.walk is not None: self.app.walk.lead = float(v); self.app.logmsg(f'second crest {v:g} of the crest: applies from the next row')
        if k in ('dotlift', 'dotsize', 'linelift'): self.app.redraw()
    def _shade(self, on): self.app.flat_shading = on; self.app.redraw()


class App(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__(); self.setWindowTitle('Lines'); self.resize(1500, 950)
        self.P = self.F = None; self.path = None; self.dots = None; self.dots_path = None; self.dots_meta = {}; self.settings = dict(SETTINGS_DEFAULTS); self.flat_shading = False; self.dlg = None
        self.walk = None; self.picked = None; self.pnode = None; self.act_dots = None; self.pickables = []; self._press = None; self.snap = None; self.settings_mesh_colour = 'lightsteelblue'; self.mesh_pv = None
        try: self.file_log, log_path = open_session_log()
        except Exception: self.file_log, log_path = None, None
        root = QtWidgets.QWidget(); self.setCentralWidget(root); lay = QtWidgets.QHBoxLayout(root)
        self.plotter = QtInteractor(root); lay.addWidget(self.plotter.interactor, 1)
        self.plotter.iren.add_observer('LeftButtonPressEvent', self._on_press); self.plotter.iren.add_observer('LeftButtonReleaseEvent', self._on_release)
        from vtkmodules.vtkRenderingCore import vtkRenderer
        self.ren_top = vtkRenderer(); self.ren_top.SetLayer(1); self.ren_top.InteractiveOff(); self.plotter.ren_win.SetNumberOfLayers(2); self.plotter.ren_win.AddRenderer(self.ren_top)
        self.ren_top.SetActiveCamera(self.plotter.renderer.GetActiveCamera())   # the overlay: lines, nodes and arrows are drawn here, over the mesh and the dots, never hidden by them
        panel = QtWidgets.QWidget(); panel.setMinimumWidth(460); panel.setMaximumWidth(560); lay.addWidget(panel, 0); side = QtWidgets.QVBoxLayout(panel); side.setContentsMargins(0, 0, 0, 0)
        self.btn_load = QtWidgets.QPushButton('Load mesh…'); self.btn_load.clicked.connect(self.load); side.addWidget(self.btn_load)
        self.lbl_file = QtWidgets.QLabel('no model'); self.lbl_file.setWordWrap(True); side.addWidget(self.lbl_file)
        g1 = QtWidgets.QGroupBox('1. Dots (from the Dots tool)'); f1 = QtWidgets.QFormLayout(g1)
        self.btn_load_dots = QtWidgets.QPushButton('Load dots…'); self.btn_load_dots.clicked.connect(self.load_dots); self.btn_load_dots.setToolTip('A dot set saved by the Dots tool (its out/dots). Opening a mesh loads its last set by itself.'); f1.addRow(self.btn_load_dots)
        self.lbl_dots = QtWidgets.QLabel('no dots'); self.lbl_dots.setWordWrap(True); f1.addRow(self.lbl_dots)
        side.addWidget(g1)

        def row(form, label, widget, tip):
            widget.setToolTip(tip); form.addRow(label, widget); form.labelForField(widget).setToolTip(tip)
        def button(form, text, fn, tip):
            b = QtWidgets.QPushButton(text); b.clicked.connect(fn); b.setToolTip(tip); form.addRow(b); return b

        g2 = QtWidgets.QGroupBox('2. Walk (the painting walk, one read at a time)'); f2 = QtWidgets.QFormLayout(g2)
        self.cb_seed = QtWidgets.QComboBox(); self.cb_seed.addItems(['Strongest free dot', 'Picked dot']); row(f2, 'Seed', self.cb_seed, TIPS['seed'])
        self.lbl_pick = QtWidgets.QLabel('no dot picked (click a dot in the viewport)'); self.lbl_pick.setWordWrap(True); f2.addRow(self.lbl_pick)
        self.btn_copy = button(f2, 'Copy the picked name', self.copy_name, 'Puts the name of the picked dot ("dot 53873") or node ("L2 a5 n7") on the clipboard.')
        self.btn_start = button(f2, 'Start walk', self.start_walk, TIPS['start'])
        self.btn_read = button(f2, 'Next read', self.next_read, TIPS['read'])
        self.btn_node = button(f2, 'Next node', self.next_node, TIPS['node'])
        self.btn_run = button(f2, 'Run to the next halt', self.run_walk, TIPS['run'])
        self.btn_runall = button(f2, 'Run to the end', self.run_all, TIPS['runall'])
        self.btn_runlines = button(f2, 'Run all lines', self.run_all_lines, TIPS['runlines'])
        self.btn_reset = button(f2, 'Reset this walk', self.reset_walk, TIPS['reset'])
        self.btn_clear = button(f2, 'Clear all lines', self.clear_lines, TIPS['clear'])
        self.lbl_state = QtWidgets.QLabel('no walk'); self.lbl_state.setWordWrap(True); f2.addRow(self.lbl_state)
        self.ed_find = QtWidgets.QLineEdit(); self.ed_find.setPlaceholderText('dot 53873  /  L2 a5 n7'); self.ed_find.returnPressed.connect(self.find_named); row(f2, 'Find by name', self.ed_find, TIPS['find'])
        side.addWidget(g2)
        self.btn_look = button(f2, 'Look at the seed', self.look_at_seed, 'Moves the camera to the seed of the walk under way (or the picked dot), at a few millimetres.')
        self.chk_dots = QtWidgets.QCheckBox('Show dots'); self.chk_dots.setChecked(True); self.chk_dots.stateChanged.connect(lambda s: self.redraw()); side.addWidget(self.chk_dots)
        self.chk_lines = QtWidgets.QCheckBox('Show lines and nodes'); self.chk_lines.setChecked(True); self.chk_lines.setToolTip('Off: the lines, their nodes and junctions are not drawn, finished or under way, so the dots behind them can be seen. The current read (front, clumps, window, direction) stays.'); self.chk_lines.stateChanged.connect(lambda s: self.redraw()); side.addWidget(self.chk_lines)
        self.chk_paint = QtWidgets.QCheckBox('Show painted dots (khaki: earlier lines, tan: this line)'); self.chk_paint.setChecked(True); self.chk_paint.setToolTip(TIPS['paint']); self.chk_paint.stateChanged.connect(lambda s: self.redraw()); side.addWidget(self.chk_paint)
        self.chk_shade = QtWidgets.QCheckBox('Shade dots by strength (light = weak, black = strong)'); self.chk_shade.setChecked(True); self.chk_shade.stateChanged.connect(lambda s: self.redraw()); side.addWidget(self.chk_shade)
        self.chk_top = QtWidgets.QCheckBox('Hidden parts of lines shown faint through the mesh; arrows always on top'); self.chk_top.setChecked(True); self.chk_top.stateChanged.connect(lambda s: self.redraw()); side.addWidget(self.chk_top)
        fw = QtWidgets.QWidget(); fh = QtWidgets.QHBoxLayout(fw); fh.setContentsMargins(0, 0, 0, 0); fh.addWidget(QtWidgets.QLabel('Hidden parts'))
        self.sl_faint = QtWidgets.QSlider(QtCore.Qt.Horizontal); self.sl_faint.setRange(0, 100); self.sl_faint.setValue(10); self.sl_faint.setToolTip('How visible the hidden parts of lines and nodes are: 0 invisible, 1 as solid as the visible parts.'); fh.addWidget(self.sl_faint, 1)
        self.lbl_faint = QtWidgets.QLabel('0.10'); fh.addWidget(self.lbl_faint); side.addWidget(fw)
        self.faint_timer = QtCore.QTimer(self); self.faint_timer.setSingleShot(True); self.faint_timer.setInterval(80); self.faint_timer.timeout.connect(self.redraw)
        self.sl_faint.valueChanged.connect(lambda v: (self.lbl_faint.setText(f'{v/100:.2f}'), self.faint_timer.start()))
        self.chk_legend = QtWidgets.QCheckBox('Legend'); self.chk_legend.setChecked(True); self.chk_legend.stateChanged.connect(lambda s: self.redraw()); side.addWidget(self.chk_legend)
        self.btn_dump = QtWidgets.QPushButton('Dump for debugging'); self.btn_dump.clicked.connect(self.dump); side.addWidget(self.btn_dump)
        self.btn_settings = QtWidgets.QPushButton('Settings…'); self.btn_settings.clicked.connect(self.open_settings); side.addWidget(self.btn_settings)
        self.txt = QtWidgets.QPlainTextEdit(); self.txt.setReadOnly(True); self.txt.setMinimumHeight(160); side.addWidget(self.txt, 1)
        self.plotter.set_background('white'); self.plotter.add_axes()
        self.logmsg(f'session log: {log_path}' if log_path else 'no session log')
        self.busy(False)

    def logmsg(self, s):
        self.txt.appendPlainText(s)
        if self.file_log: self.file_log(s)
    def busy(self, b):
        self.btn_load.setEnabled(not b); self.btn_load_dots.setEnabled(not b and self.P is not None)
        w = self.walk; live = w is not None and w.cur is not None and w.cur['halt'] is None
        self.btn_start.setEnabled(not b and self.dots is not None and not live)
        for x in (self.btn_read, self.btn_node, self.btn_run, self.btn_runall): x.setEnabled(not b and live)
        self.btn_runlines.setEnabled(not b and self.dots is not None)
        self.btn_reset.setEnabled(not b and w is not None and w.cur is not None); self.btn_clear.setEnabled(not b and w is not None)
    def open_settings(self):
        if self.dlg is None: self.dlg = SettingsDialog(self)
        self.dlg.show(); self.dlg.raise_()

    def load(self):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, 'Load mesh', '', 'Meshes (*.stl *.obj *.ply)')
        if path: self._use_mesh(path)

    def _use_mesh(self, path, dots=None):
        """load a mesh, then its dots: the set given, else the newest set saved for it"""
        try: self.P, self.F = load_mesh(path)
        except Exception as ex: QtWidgets.QMessageBox.critical(self, 'Load failed', str(ex)); return
        self.path = path; self.lbl_file.setText(os.path.basename(path)); self.dots = None; self.dots_path = None; self.walk = None; self.picked = None; self.pnode = None; self.snap = None; self.mesh_pv = None
        self.logmsg(f'loaded {os.path.basename(path)}: {len(self.F)} faces')
        p = dots if dots and os.path.exists(dots) else last_dots(path)
        if p: self._use_dots(p)
        self._remember(); self.busy(False); self.redraw(reset=True)

    def _remember(self):
        """the model and dot set in use, kept in last.json beside the app, reloaded at the next start"""
        try:
            with open(LAST_FILE, 'w', encoding='utf-8') as fh: json.dump(dict(model=self.path, dots=self.dots_path), fh, indent=1)
        except Exception as ex: self.logmsg(f'could not remember the files in use: {ex}')

    def reload_last(self):
        try:
            with open(LAST_FILE, encoding='utf-8') as fh: last = json.load(fh)
        except Exception: return
        if last.get('model') and os.path.exists(last['model']): self._use_mesh(last['model'], last.get('dots'))

    def _use_dots(self, p):
        try: fld, cut, cube, blend, reach, steep, model = load_dots(p)
        except Exception as ex: self.logmsg(f'could not load dots from {p}: {ex}'); return
        self.dots = fld; self.dots_path = p; self.dots_meta = dict(cut=cut, cube=cube, reach=reach, steep=steep); self.walk = None; self.picked = None; self.pnode = None
        self.lbl_dots.setText(f'{os.path.basename(p)}: {len(fld["idx"])} dots, cube {cube:g} mm, reach {reach}, cut {cut:g}')
        self.logmsg(f'loaded {len(fld["idx"])} dots from {os.path.basename(p)} (cube {cube:g} mm, reach {reach}, cut {cut:g})')
        self.lbl_pick.setText('no dot picked (click a dot in the viewport)'); self.lbl_state.setText('no walk')

    def load_dots(self):
        p, _ = QtWidgets.QFileDialog.getOpenFileName(self, 'Load dots', DOTS_DIR, 'Dot sets (*.npz)')
        if not p: return
        self._use_dots(p); self._remember(); self.busy(False); self.redraw()

    # ---- the walk ----
    def _ensure_walk(self):
        if self.walk is None:
            t0 = time.time()
            if self.snap is None: self.snap = make_snap(self.P, self.F)
            self.walk = walkmod.Walk(self.dots, snap=self.snap, log=self.logmsg, cap=int(self.settings['cap']), reach=int(self.settings['reach']), window=int(self.settings['window']), gap=float(self.settings['gap']), fork=float(self.settings['fork']), kink=float(self.settings['kink']), lead=float(self.settings['lead']))
            self.logmsg(f'walk ready: {self.walk.k} dots, neighbour table at reach {self.walk.reach} in {time.time()-t0:.1f} s')
        return self.walk

    def start_walk(self):
        if self.dots is None: return
        if self.walk is not None and self.walk.reach != int(self.settings['reach']):   # the reach changed: a new neighbour table, the finished lines kept
            old = self.walk; self.walk = None; w = self._ensure_walk(); w.owner[:] = old.owner; w.lines = old.lines
        w = self._ensure_walk()
        if w.cur is not None and w.cur['halt'] is None: self.logmsg('a walk is under way: step it, run it or reset it'); return
        if self.cb_seed.currentIndex() == 1:
            if self.picked is None: self.logmsg('pick a dot first: click it in the viewport'); return
            seed = self.picked
        else:
            seed = w.strongest_free()
            if seed is None: self.logmsg('no free dot left'); return
        if w.start(seed): self._after()

    def _after(self):
        w = self.walk; c = w.cur
        if c is None: self.lbl_state.setText('no walk')
        else:
            s = f'line {c["id"]}: {c["stage"]}, round {c["round"]}' + (f', arm {c["arm"]} of {len(c["arms"])} ({len(c["queue"])} queued)' if c['arm'] >= 0 else '') + f'; next read: {c["phase"]}'
            if c['halt'] is not None: s = f'line {c["id"]} halted: {c["halt"]}'
            self.lbl_state.setText(s)
        self.busy(False); self.redraw()
    def next_read(self):
        if self.walk is None or self.walk.cur is None: return
        self.walk.step_read(); self._after()
    def next_node(self):
        if self.walk is None or self.walk.cur is None: return
        self.walk.step_node(); self._after()
    def run_walk(self):
        if self.walk is None or self.walk.cur is None: return
        t0 = time.time(); self.walk.run_arm(int(self.settings['maxnodes'])); self.logmsg(f'run to the next halt: {time.time()-t0:.1f} s'); self._after()
    def run_all(self):
        if self.walk is None or self.walk.cur is None: return
        t0 = time.time(); self.walk.run(int(self.settings['maxnodes'])); self.logmsg(f'run to the end: {time.time()-t0:.1f} s'); self._after()
    def run_all_lines(self):
        """every line of the set, one after the other: the walk under way is finished first, then the strongest free dot seeds the
        next line until none is left. The reads go to the session log file only; the box gets one line per finished line."""
        if self.dots is None: return
        w = self._ensure_walk(); box_log = w.log; w.log = self.file_log or (lambda s: None); t0 = time.time(); n0 = len(w.lines); cap = int(self.settings['maxnodes']); self.busy(True)
        try:
            if w.cur is not None and w.cur['halt'] is None: w.run(cap)
            while True:
                seed = w.strongest_free()
                if seed is None: break
                if not w.start(seed): break
                w.run(cap); L = w.lines[-1]
                self.logmsg(f'line {L["id"]}: seed dot {L["seed"]}, {len(L.get("polylines", []))} arm(s), {len(L["nodes"])} nodes, {len(L.get("junctions", []))} junction(s); {L["halt"][:80]}')
                if (len(w.lines)-n0) % 20 == 0: QtWidgets.QApplication.processEvents()
        finally: w.log = box_log
        done = w.lines[n0:]; self.logmsg(f'run all lines: {len(done)} lines, {sum(len(L["nodes"]) for L in done)} nodes, {sum(len(L.get("junctions", [])) for L in done)} junctions, {int((w.owner != -1).sum())} of {w.k} dots painted, {time.time()-t0:.1f} s')
        self._after()

    def reset_walk(self):
        if self.walk is None: return
        self.walk.reset(); self.logmsg('walk reset'); self._after()
    def clear_lines(self):
        if self.walk is None: return
        self.walk.clear(); self.logmsg('all lines cleared'); self._after()

    def look_at_seed(self):
        i = self.walk.cur['seed'] if self.walk is not None and self.walk.cur is not None else self.picked
        if i is None or self.dots is None: return
        self._look_at(self._pts(np.array([i]))[0], self.dots['N'][i])
    def _look_at(self, c, n=None):
        r = self.dots; d = 40*r['size']; n = n if n is not None else self.snap.normals(c)[0] if self.snap is not None else np.array([0, 0, 1.0])
        self.plotter.camera_position = [tuple(c+d*n), tuple(c), (0, 0, 1)]; self.plotter.render()

    # ---- names: a dot by its number in the set, a node by line, arm and number as the log has them ----
    def _node_name(self, line, arm, k):
        return f'L{line} seed' if arm < 0 else f'L{line} a{arm} start' if k < 0 else f'L{line} a{arm} n{k+1}'
    def _select_dot(self, i):
        r = self.dots; self.picked = int(i); self.pnode = None; c = r['lo']+(r['idx'][i]+0.5)*r['size']
        extra = ''.join(f', {k} {float(r[k][i]):.2f}' for k in ('spread', 'direction', 'line', 'concentration') if k in r)
        own = '' if self.walk is None or self.walk.owner[i] == -1 else f', painted by line {self.walk.owner[i]}'
        s = f'dot {i} (cube {tuple(int(v) for v in r["idx"][i])}): strength {float(r["turn"][i]):.1f}{extra}{own}'
        self.lbl_pick.setText(s); self.logmsg(f'picked {s}, at {np.round(c, 4).tolist()}'); return c
    def _select_node(self, name, pos):
        self.picked = None; self.pnode = (name, np.asarray(pos, float)); self.lbl_pick.setText(f'node {name} at {np.round(pos, 4).tolist()}'); self.logmsg(f'picked node {name} at {np.round(pos, 4).tolist()}'); return self.pnode[1]
    def copy_name(self):
        name = f'dot {self.picked}' if self.picked is not None else self.pnode[0] if self.pnode is not None else None
        if name is None: self.logmsg('nothing picked to copy'); return
        QtWidgets.QApplication.clipboard().setText(name); self.logmsg(f'copied "{name}" to the clipboard')
    def find_named(self):
        """the Find box: 'dot N' or 'N' selects dot N; 'L2 a5 n7', 'L2 seed', 'L2 a5 start' a node of the lines on screen"""
        if self.dots is None: return
        t = self.ed_find.text().strip().replace('.', ' ').lower(); parts = t.split()
        if not parts: return
        if parts[0] == 'dot' or parts[0].isdigit():
            try: i = int(parts[-1])
            except ValueError: self.logmsg(f'find: "{t}" is not a dot number'); return
            if not 0 <= i < len(self.dots['idx']): self.logmsg(f'find: dot {i} is not in this set ({len(self.dots["idx"])} dots)'); return
            c = self._select_dot(i); self.redraw(); self._look_at(c, self.dots['N'][i]); return
        want = ' '.join(parts)
        for actor, kind, names, pts in self.pickables:
            if kind != 'node': continue
            for q, nm in enumerate(names):
                if nm.lower() == want: c = self._select_node(nm, pts[q]); self.redraw(); self._look_at(c); return
        self.logmsg(f'find: no node named "{t}" on screen (names are like "L2 a5 n7", "L2 seed", "L2 a5 start")')

    # ---- picking ----
    def _on_press(self, *_): self._press = self.plotter.iren.get_event_position()
    def _on_release(self, *_):
        pos = self.plotter.iren.get_event_position()
        if self._press is None or abs(pos[0]-self._press[0]) > 3 or abs(pos[1]-self._press[1]) > 3: return   # a drag rotated the view
        if not self.pickables or self.dots is None: return
        from vtkmodules.vtkRenderingCore import vtkPointPicker
        pk = vtkPointPicker(); pk.SetTolerance(self.settings['picktol']); pk.PickFromListOn()
        for actor, *_ in self.pickables: pk.AddPickList(actor)
        if not pk.Pick(pos[0], pos[1], 0, self.plotter.renderer) or pk.GetPointId() < 0: return
        q = int(pk.GetPointId()); hit = pk.GetActor()
        for actor, kind, names, pts in self.pickables:
            if actor is not hit: continue
            if kind == 'dot': self._select_dot(q)
            else: self._select_node(names[q], pts[q])
            break
        self.redraw()

    # ---- drawing ----
    def _pts(self, sel=None):
        r = self.dots; idx = r['idx'] if sel is None else r['idx'][sel]; N = r['N'] if sel is None else r['N'][sel]
        return r['lo']+(idx+0.5)*r['size']+self.settings['dotlift']*r['size']*N
    def _top(self, actor):
        """moves an actor to the overlay layer when the overlay is on"""
        if actor is not None and self.chk_top.isChecked(): self.plotter.renderer.RemoveActor(actor); self.ren_top.AddActor(actor)
        return actor
    def _lift(self, pts):
        """drawn lines and nodes sit at the dots' height above the surface (dot lift plus half a cube), so the dots do not hide them"""
        pts = np.asarray(pts, float).reshape(-1, 3)
        if self.snap is None or self.dots is None or len(pts) == 0: return pts
        if self.settings['linelift'] <= 0: return pts
        return pts+self.settings['linelift']*float(self.dots['size'])*self.snap.normals(pts)
    @property
    def FAINT(self): return self.sl_faint.value()/100   # opacity of the hidden copy, the slider
    def _points(self, pts, color, size, top=False):
        if not len(pts): return
        pts = self._lift(pts) if top else np.asarray(pts, float).reshape(-1, 3)
        self.plotter.add_mesh(pv.PolyData(pts), color=color, point_size=size, render_points_as_spheres=True, reset_camera=False)
        if top and self.chk_top.isChecked() and self.FAINT > 0: self._top(self.plotter.add_mesh(pv.PolyData(pts), color=color, point_size=max(1, size//2), render_points_as_spheres=True, opacity=self.FAINT, reset_camera=False))
    def _polylines_batch(self, pls, mfs, color, width, size, names=None, line=None):
        """many polylines as one actor, their nodes as one actor coloured by misfit (two more for the faint copies): a tree of
        hundreds of arms drew thousands of actors and made every step sluggish"""
        pls = [np.asarray(pl, float).reshape(-1, 3) for pl in pls if len(pl)]
        if not pls: return
        lifted = [self._lift(pl) for pl in pls]; pts = np.vstack(lifted); n = [len(pl) for pl in lifted]; off = np.cumsum([0]+n[:-1])
        cells = np.concatenate([np.r_[k, np.arange(o, o+k)] for k, o in zip(n, off) if k > 1]) if any(k > 1 for k in n) else None
        mf = np.concatenate([np.asarray(m, float) if m is not None and len(m) == k else np.full(k, self.FIT_LO) for m, k in zip(mfs, n)])
        if cells is not None:
            self.plotter.add_mesh(pv.PolyData(pts, lines=cells), color=color, line_width=width, reset_camera=False)
            if self.chk_top.isChecked() and self.FAINT > 0: self._top(self.plotter.add_mesh(pv.PolyData(pts, lines=cells), color=color, line_width=max(1, width//2), opacity=self.FAINT, reset_camera=False))
        d = pv.PolyData(pts); d['misfit'] = mf
        act = self.plotter.add_mesh(d, scalars='misfit', cmap='coolwarm', clim=(self.FIT_LO, self.FIT_HI), point_size=size, render_points_as_spheres=True, show_scalar_bar=False, reset_camera=False)
        if names is not None:   # the nodes are pickable, each with its name (line, arm, node number); positions unlifted, as the nodes are
            nm = [self._node_name(li, a, k) for nl, li in zip(names, line) for a, k in nl]
            if len(nm) == len(pts): self.pickables.append((act, 'node', nm, np.vstack(pls)))
        if self.chk_top.isChecked() and self.FAINT > 0: self._top(self.plotter.add_mesh(pv.PolyData(pts), color=color, point_size=max(1, size//2), render_points_as_spheres=True, opacity=self.FAINT, reset_camera=False))

    FIT_LO, FIT_HI = 5.0, 20.0   # misfit (deg) drawn fully blue, fully red
    def _fit_points(self, pts, mf, size):
        """nodes coloured by their arrows' misfit to one line: blue 5 deg or less, red 20 deg or more"""
        pts = self._lift(pts); mf = np.asarray(mf, float)
        if len(pts) == 0 or len(mf) != len(pts): self._points(pts, 'blue', size, top=True); return
        d = pv.PolyData(pts); d['misfit'] = mf
        self.plotter.add_mesh(d, scalars='misfit', cmap='coolwarm', clim=(self.FIT_LO, self.FIT_HI), point_size=size, render_points_as_spheres=True, show_scalar_bar=False, reset_camera=False)
        if self.chk_top.isChecked() and self.FAINT > 0: self._top(self.plotter.add_mesh(pv.PolyData(pts), color='blue', point_size=max(1, size//2), render_points_as_spheres=True, opacity=self.FAINT, reset_camera=False))
    def _polyline(self, pts, color, width):
        pts = self._lift(pts)
        if len(pts) < 2: return
        self.plotter.add_mesh(pv.PolyData(pts, lines=np.r_[len(pts), np.arange(len(pts))]), color=color, line_width=width, reset_camera=False)
        if self.chk_top.isChecked() and self.FAINT > 0: self._top(self.plotter.add_mesh(pv.PolyData(pts, lines=np.r_[len(pts), np.arange(len(pts))]), color=color, line_width=max(1, width//2), opacity=self.FAINT, reset_camera=False))

    def _arrow(self, p, d, color, length):
        d = np.asarray(d, float); n = np.linalg.norm(d)
        if n < 1e-12: return
        self._top(self.plotter.add_mesh(pv.Arrow(start=np.asarray(p, float), direction=d/n, scale=length, tip_length=0.3, tip_radius=0.12, shaft_radius=0.05), color=color, reset_camera=False))

    def redraw(self, reset=False):
        if self.P is None: return
        cam = self.plotter.camera_position; self.plotter.renderer.clear_actors(); self.ren_top.RemoveAllViewProps(); self.act_dots = None; self.pickables = []
        if self.plotter.scalar_bars: self.plotter.remove_scalar_bar()
        try: self.plotter.remove_legend()
        except Exception: pass
        if self.mesh_pv is None:   # built once: the shading normals cost 0.4 s per redraw on 1.5M faces when left to add_mesh
            m = pv.PolyData(self.P, np.hstack([np.full((len(self.F), 1), 3), self.F]).ravel()); self.mesh_pv = (m.compute_normals(cell_normals=False, point_normals=True, split_vertices=False), m.copy())
        self.plotter.add_mesh(self.mesh_pv[1] if self.flat_shading else self.mesh_pv[0], color=self.settings_mesh_colour, smooth_shading=False)   # a tint, so the grey dots stand out on it; point normals give the smooth look
        r = self.dots; ds = int(self.settings['dotsize']); legend = []
        if r is not None and self.chk_dots.isChecked() and len(r['idx']):
            d = pv.PolyData(self._pts())
            if self.chk_shade.isChecked():
                d['strength'] = r['turn']; lo_ = float(self.dots_meta.get('cut', 0)); hi = max(float(np.percentile(r['turn'], 99)), lo_+1e-3)
                self.act_dots = self.plotter.add_mesh(d, scalars='strength', cmap='Greys', clim=(lo_, hi), point_size=ds, render_points_as_spheres=True,
                                                      scalar_bar_args=dict(title='strength', vertical=True, position_x=0.88, position_y=0.3, height=0.4, width=0.06))
            else: self.act_dots = self.plotter.add_mesh(d, color='dimgrey', point_size=ds, render_points_as_spheres=True)
            self.pickables.append((self.act_dots, 'dot', None, None))
        w = self.walk; big = ds+6; huge = ds+12; cube = 0 if r is None else float(r['size'])
        if w is not None and r is not None:
            owned = np.flatnonzero(w.owner != -1)
            if len(owned) and self.chk_paint.isChecked(): self._points(self._pts(owned), 'khaki', ds); legend.append(('painted by an earlier line', 'khaki'))
            pls = [pl for L in w.lines for pl in L.get('polylines', [])]; mfs = [mf for L in w.lines for mf in L.get('misfits', [None]*len(L.get('polylines', [])))]
            nms = [nl for L in w.lines for nl in L.get('names', [])]; lids = [L['id'] for L in w.lines for nl in L.get('names', [])]
            show_lines = self.chk_lines.isChecked()
            if show_lines:
                self._polylines_batch(pls, mfs, 'navy', 4, 6, names=nms if len(nms) == len(pls) else None, line=lids)
                jf = [pos for L in w.lines for pos, deg, fixed in L.get('junctions', []) if fixed]; jp = [pos for L in w.lines for pos, deg, fixed in L.get('junctions', []) if not fixed]
                if jf: self._points(jf, 'purple', huge, top=True)
                if jp: self._points(jp, 'magenta', huge, top=True)
                if w.lines: legend.append(('finished lines', 'navy'))
            st = w.state()
            if st is not None:
                if self.chk_paint.isChecked(): self._points(self._pts(st['painted']), 'tan', ds+1); legend.append(("this line's paint", 'tan'))
                self._points(self._pts(st['front']), 'white', big); legend.append(('the front (painted last round)', 'white'))
                if st.get('lane') is not None: self._points(self._pts(st['lane']), 'orchid', big); legend.append((f'second band under the same paint ({len(st["lane"])} dots): an arm of its own once the paint parts', 'orchid'))
                if st['probe'] is not None:
                    for q, (g, live) in enumerate(zip(st['probe']['halves'], st['probe']['live'])):
                        col = walkmod.CLUMP_COLOURS[q % len(walkmod.CLUMP_COLOURS)]; self._points(self._pts(g), col, big); legend.append((f'split under test, half {q} ({len(g)} dots{", dead" if not live else ""}), round {st["probe"]["rounds"]} of {self.walk.window}', col))
                elif st['clumps']:
                    for q, g in enumerate(st['clumps']):
                        col = walkmod.CLUMP_COLOURS[q % len(walkmod.CLUMP_COLOURS)]; self._points(self._pts(g), col, big); legend.append((f'new dots, clump {q} ({len(g)})', col))
                elif len(st['new']): self._points(self._pts(st['new']), 'yellow', big); legend.append((f'new dots, not yet clumped ({len(st["new"])})', 'yellow'))
                seedp = self._pts(np.array([st['seed']]))[0]; self._points([seedp], 'red', huge, top=True); legend.append(('seed', 'red'))
                if show_lines:
                    self._polylines_batch(st['polylines'], st['misfits'], 'blue', 5, big, names=st['names'] if len(st['names']) == len(st['polylines']) else None, line=[st['id']]*len(st['polylines']))
                    legend.append(('nodes: blue where the arrows fit one line (a crease), red where they do not (a fork)', 'blue'))
                    jf = [pos for pos, fixed in st['junctions'] if fixed]; jp = [pos for pos, fixed in st['junctions'] if not fixed]
                    if jf: self._points(jf, 'purple', huge, top=True)
                    if jp: self._points(jp, 'magenta', huge, top=True)
                    if st['junctions']: legend.append(('junction, provisional', 'magenta')); legend.append(('junction, corrected', 'purple'))
                    if st['newest'] is not None: self._points([st['newest']], 'blue', huge+4, top=True)
                # the direction: at the finding stage from the seed to each clump; while advancing, from the newest node (or the seed) toward the front
                if len(st['region']):
                    m = st['member'] if st['member'] is not None else np.ones(len(st['region'])); sec = st['second'] if st['second'] is not None else np.zeros(len(st['region']), bool)
                    own = st['region'][(m >= 0.1) & ~sec]; odd = st['region'][m < 0.1]; lead2 = st['region'][(m >= 0.1) & sec]
                    self._points(self._pts(own), 'lightgreen', ds+2); legend.append((f'node window: the rows read ({len(own)} dots)', 'lightgreen'))
                    if len(lead2): self._points(self._pts(lead2), 'hotpink', ds+4); legend.append((f'second lead in the band ({len(lead2)} dots): climb {st["climb"]:.2f} of the crest' + (f', arm mean {st["climb_mean"]:.2f}' if st['climb_mean'] is not None else ''), 'hotpink'))
                    elif st['climb'] is not None: legend.append((f'one lead: climb {st["climb"]:.2f} of the crest' + (f', arm mean {st["climb_mean"]:.2f}' if st['climb_mean'] is not None else ''), 'lightgreen'))
                    if len(odd): self._points(self._pts(odd), 'salmon', ds+2); legend.append((f'dots of another band, arrows off this direction ({len(odd)})', 'salmon'))
                if st['direction'] is not None:
                    tail = st['rownode'][0] if st['rownode'] is not None else (st['newest'] if st['newest'] is not None else seedp)
                    self._arrow(tail, st['direction'], 'darkgreen', 8*cube); legend.append((f'crease direction: the arrows fit one line within {st["misfit"] if st["misfit"] is not None else 0:.0f} deg' + (f'; own row off the stem by {st["share"]:.0f} deg' if st['share'] is not None else ''), 'darkgreen'))
                if st['rownode'] is not None:
                    cr, cg = st['rownode']; self._points([cr], 'yellow', big, top=True); self._points([cg], 'lightgreen', big+2, top=True); legend.append(("the row's centre (along)", 'yellow')); legend.append(("the window's centre (across)", 'lightgreen'))
        if self.picked is not None and r is not None: self._points(self._pts(np.array([self.picked])), 'lime', huge, top=True); legend.append((f'picked dot {self.picked}', 'lime'))
        if self.pnode is not None: self._points([self.pnode[1]], 'lime', huge+4, top=True); legend.append((f'picked node {self.pnode[0]}', 'lime'))
        if w is not None and r is not None and self.chk_top.isChecked(): legend.append(('faint: hidden behind geometry', 'grey'))
        if legend and self.chk_legend.isChecked(): self.plotter.add_legend(legend, bcolor='white', border=True, size=(0.28, 0.02*len(legend)+0.02), loc='upper left', face='circle')
        if reset: self.plotter.reset_camera()
        else: self.plotter.camera_position = cam
        self.plotter.render()

    def dump(self):
        d = os.path.join(DEBUG_DIR, time.strftime('%Y%m%d_%H%M%S')); os.makedirs(d, exist_ok=True)
        try:
            prm = dict(model=self.path, dots=self.dots_path, dots_meta=self.dots_meta, settings=dict(self.settings), seed_mode=self.cb_seed.currentText(), picked=self.picked, camera=[list(map(float, v)) for v in self.plotter.camera_position])
            w = self.walk
            if w is not None:
                lines = w.lines+([dict(id=w.cur['id'], seed=w.cur['seed'], polylines=w.state()['polylines'], halt=w.cur['halt'] or 'under way', arms=[], records=w.cur['records'])] if w.cur is not None and not any(l['id'] == w.cur['id'] for l in w.lines) else [])
                pls = [(l['id'], pl) for l in lines for pl in l.get('polylines', []) if len(pl)]
                np.savez(os.path.join(d, 'lines.npz'), n=np.array([len(pl) for _, pl in pls]), line=np.array([i for i, _ in pls]), pts=np.vstack([np.asarray(pl).reshape(-1, 3) for _, pl in pls]) if pls else np.zeros((0, 3)), owner=w.owner, seeds=np.array([l['seed'] for l in lines]))
                with open(os.path.join(d, 'walk.json'), 'w', encoding='utf-8') as fh: json.dump([dict(id=l['id'], seed=int(l['seed']), polylines=[len(pl) for pl in l.get('polylines', [])], halt=l['halt'], arms=l['arms'], records=l['records']) for l in lines], fh, indent=1)
                prm['lines'] = len(lines)
            with open(os.path.join(d, 'settings.json'), 'w', encoding='utf-8') as fh: json.dump(prm, fh, indent=1)
            with open(os.path.join(d, 'log.txt'), 'w', encoding='utf-8') as fh: fh.write(self.txt.toPlainText())
            self.plotter.screenshot(os.path.join(d, 'viewport.png')); self.logmsg(f'dumped lines, walk records, settings, log, screenshot to {d}')
        except Exception: self.logmsg(traceback.format_exc())


if __name__ == '__main__':
    qa = QtWidgets.QApplication(sys.argv); w = App(); w.show(); QtCore.QTimer.singleShot(0, w.reload_last); sys.exit(qa.exec())
