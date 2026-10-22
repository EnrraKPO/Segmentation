"""Segmentation tool. Two stages, each run on its own and kept:
  1. Find dots   : cubes over the surface, each reading its own triangle's facing; a dot where the facing shift
                   changes sharply from one cube to the next (creases.find_dots). Parameters: cube size (mm), cut.
  3. Refine lines: a tube of fine cubes is laid around each line (refine.py), each fine cube's bend read over the
                   surrounding area, the fine dots are the crease, and the same walk draws the line through them.
                   Parameters: finer (times smaller than the coarse cube), tube radius (mm), area (mm), cut. The line
                   to refine is picked in the viewport: hover highlights the line under the cursor, a click selects it.
  2. Build lines : anchored walk (trace4.walk): from the strongest free dot, read the band of connected dots in a thin
                   slice across it, place a node, step along the line's trend to the next anchor; lines end on other
                   lines. Parameters in mm: step, slice, reach; trend pull; min length. The old thick-stroke method
                   (creases.strokes + creases.join) stays selectable for comparison.
Run: python app.py (from this folder)."""
import os, sys, time, traceback
import numpy as np, trimesh
from PySide6 import QtCore, QtWidgets
import pyvista as pv
from pyvistaqt import QtInteractor
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import creases, trace4, refine


def pick_line(plotter, act_mesh, drawn, x, y, tol=0.006, depth_mm=0.5):
    """The drawn line under window position (x, y), or None. drawn: {kind: (actor, pts, n, off, index)}. Lines are picked
    within tol (fraction of the window diagonal) of the cursor; the surface is picked at the cursor too, and a line lying more than depth_mm behind it
    (hidden inside or behind the model) is not taken. Returns (kind, k, point): the actor's kind, the line's position in
    its curve list, the picked point."""
    from vtkmodules.vtkRenderingCore import vtkCellPicker
    ren = plotter.renderer; cam = ren.GetActiveCamera(); cp = np.array(cam.GetPosition()); vd = np.array(cam.GetDirectionOfProjection())
    depth = lambda p: float((np.asarray(p)-cp) @ vd)
    lp = vtkCellPicker(); lp.SetTolerance(tol); lp.PickFromListOn()
    for a, *_ in drawn.values(): lp.AddPickList(a)
    if not lp.Pick(x, y, 0, ren) or lp.GetCellId() < 0: return None
    hit = lp.GetActor(); cell = lp.GetCellId(); pt = np.array(lp.GetPickPosition())
    if act_mesh is not None:
        mp = vtkCellPicker(); mp.SetTolerance(0.0005); mp.PickFromListOn(); mp.AddPickList(act_mesh)
        if mp.Pick(x, y, 0, ren) and mp.GetCellId() >= 0 and depth(pt)-depth(mp.GetPickPosition()) > depth_mm: return None
    for kind, (a, pts, n, off, index) in drawn.items():
        if a is hit: return kind, index[cell], pt
    return None


DOTS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out', 'dots')   # saved dot sets; last.json maps model name -> last file


def save_dots(res, model, cube, cut, flat):
    """Keep only the dots (cubes above the cut): idx, N, turn, plus lo, size, cut and the settings. Written to out/dots/
    as <model>_<cube>mm_cut<cut>[_blend]_<time>.npz and noted as the last set for this model."""
    import json
    os.makedirs(DOTS_DIR, exist_ok=True); m = res['line']; base = os.path.splitext(os.path.basename(model))[0]
    name = f'{base}_{cube:g}mm_cut{cut:g}{"" if flat else "_blend"}_{time.strftime("%Y%m%d_%H%M%S")}.npz'; path = os.path.join(DOTS_DIR, name)
    np.savez(path, idx=res['idx'][m], N=res['N'][m], turn=res['turn'][m], lo=res['lo'], size=res['size'], cut=res['cut'], cube=cube, flat=flat, model=model)
    lj = os.path.join(DOTS_DIR, 'last.json'); last = json.load(open(lj)) if os.path.exists(lj) else {}
    last[os.path.basename(model)] = name; json.dump(last, open(lj, 'w'), indent=1)
    return path


def load_dots(path):
    """A saved dot set back into the res dict the line stage reads (every stored cube is a dot)."""
    z = np.load(path, allow_pickle=True); k = len(z['idx'])
    res = dict(idx=z['idx'], N=z['N'], turn=z['turn'], line=np.ones(k, bool), lo=z['lo'], size=float(z['size']), cut=float(z['cut']))
    return res, float(z['cube']), bool(z['flat']), str(z['model'])


def last_dots(model):
    import json
    lj = os.path.join(DOTS_DIR, 'last.json')
    if not os.path.exists(lj): return None
    name = json.load(open(lj)).get(os.path.basename(model)); p = name and os.path.join(DOTS_DIR, name)
    return p if p and os.path.exists(p) else None


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
    method="Anchored walk: the new line stage (see docs/brief_anchored_walk.md). Thick strokes: the previous method, kept for comparison.",
    step="Distance between nodes along the line, in mm. About one cube.",
    slice="Half-length, in mm, of the slice across the band that each node averages. Longer = smoother nodes, but a long slice blurs tight corners.",
    look="Probe length in mm, 0 = off. Before each node, two corridors this long are probed for dots: one continuing the trend, one along the local dot axis. Dot support ahead scales the trend pull: the trend wins where the dots continue along it, the dots win where the trend leads onto bare surface.",
    levp="0 = off. Each line keeps a running standard, the typical strength of the dots it follows. A dot's say is (strength / standard) to this power, capped at 1: with power 6 a dot 10 % below the standard has half a say, one 50 % below almost none. Weak dots lead only where nothing stronger is there.",
    reach="How far ahead, in mm, the next anchor may lie. A gap in the dot row shorter than this does not end the line; longer reaches may jump to another band.",
    pull="How much the line's trend pulls each node toward where the line was heading (0 = node sits on the band alone, 1 = equal say). Smooths dents.",
    wminlen="Lines with fewer nodes than this are dropped and their dots released.",
    pick="Move the mouse over a line in the viewport: it lights up yellow. Click it (without dragging) to select it, green. Click bare surface to clear. 'Refine one line' refines the selected line; a selected refined line refines its walk line again.",
    finer="How many times smaller than the coarse cube the tube's cubes are. 3 gives the same line as 10 within a few hundredths of a millimetre at a fraction of the time.",
    radius="Half-width of the tube around the walked line, in mm. Must hold the crease's band with room for the line to move.",
    area="Width in mm of the neighbourhood over which each fine cube's bend is read, smooth weights. The coarse cube: the cut keeps its meaning.",
    rcut="Bend, in degrees per coarse cube, a fine cube needs to be a dot. Same meaning as the stage 1 cut.",
    snap="Moves each node onto the nearest surface point before drawing, so lines on sharp convex edges are not hidden inside the solid. Judge placement with this on; lifting lines off the surface shifts them sideways.",
)


# Every constant of the pipeline that is not on the main panel, as a handle: (key, group, label, default, lo, hi, step, decimals, meaning).
# Defaults are the values the code ran with before the handles existed (2026-10-22); nothing here changes a rule.
TUNE = [
    # walk (trace4.walk), lengths in mm unless said; 'unit' is the length the walk's cube-sized rules were written for
    ('wmax', 'Walk', 'Slab half-width across (mm)', 1.0, 0.05, 10, 0.05, 2, 'Widest band half-width the slice looks at across the axis. Dots further out across are never part of a band.'),
    ('sphere', 'Walk', 'Axis sphere radius (mm)', 0.7, 0.05, 10, 0.05, 2, 'Radius of the sphere of dots the axis (the direction the line runs in) is read from at each anchor.'),
    ('sphgrow', 'Walk', 'Axis sphere grows to (x band width)', 3.0, 0, 20, 0.5, 1, 'On a wide band the axis sphere grows to this many times the band width, so a wide band still reads as a line.'),
    ('align', 'Walk', 'Alignment ratio (long/wide)', 1.5, 1, 10, 0.1, 2, 'The dots in the sphere must spread at least this many times more along the axis than across it, else the axis is uncertain and the trend carries the line (or the seed is dropped).'),
    ('axismin', 'Walk', 'Dots needed for an axis', 4, 2, 100, 1, 0, 'Fewer dots than this in the sphere: no axis.'),
    ('cone', 'Walk', 'Anchor cone (deg)', 45, 1, 90, 1, 0, 'The next anchor must lie within this angle of the trend direction, seen from the last node.'),
    ('trend', 'Walk', 'Trend length (nodes)', 6, 1, 100, 1, 0, 'The trend direction is the chord over this many latest nodes.'),
    ('afall', 'Walk', 'Anchor score falloff (cubes)', 2.0, 0.1, 50, 0.1, 1, 'A candidate anchor this many cubes off the trend line scores 1/e; the best score wins.'),
    ('along', 'Walk', 'Along-score falloff (mm, 0 = reach)', 0.0, 0, 5, 0.01, 3, 'A candidate anchor this far from one step ahead scores 1/e along; 0 couples it to the reach.'),
    ('slack', 'Walk', 'Anchor slack beyond band half-width (mm)', 0.087, 0, 2, 0.01, 3, "How far beyond the band's half-width the next anchor may sit off the trend line."),
    ('touch', 'Walk', 'Touch distance (units)', 1.8, 0.5, 10, 0.1, 2, 'Dots within this many units of each other are joined; the band is the dots joined to the anchor. 1.73 is the cube diagonal (neighbours including diagonals).'),
    ('unit', 'Walk', 'Unit (mm)', 0.087, 0.001, 5, 0.001, 3, 'The length the cube-sized rules (touch, self-hit, profile bins) were written for, so the walk reads the same at any cube size.'),
    ('selfhit', 'Walk', 'Self-hit distance (units)', 1.5, 0.1, 20, 0.1, 2, 'A new node this close to an earlier node of its own line closes the loop.'),
    ('backmin', 'Walk', 'Self-hit ignores latest nodes (min)', 15, 1, 1000, 1, 0, 'At least this many latest nodes are ignored for the self-hit (more when the slice is long against the step).'),
    ('wmem', 'Walk', 'Band width memory (share kept)', 0.7, 0, 1, 0.05, 2, 'The remembered band width keeps this share of its old value per node and takes the rest from the new slice.'),
    ('owned', 'Walk', 'Owned share that ends a line', 0.5, 0, 1, 0.05, 2, "A band with more than this share of its dots owned by another line is that line's: the walk ends on it (and no seed starts there)."),
    ('levn', 'Walk', 'Strength memory (nodes)', 6, 1, 100, 1, 0, "The line's running strength standard moves toward the band's strength by 1/this per node."),
    ('lookw', 'Walk', 'Look-ahead corridor half-width (cubes)', 1.0, 0.1, 20, 0.1, 1, 'Half-width of the probe corridors; dots near the corridor centre line count most.'),
    ('cluster', 'Walk', 'Cluster radius (units)', 3.0, 1, 50, 0.5, 1, 'Dots within this many units (or reach + one unit, whichever is larger) are one cluster; each cluster is walked on its own core and lines never cross clusters.'),
    ('maxnodes', 'Walk', 'Node cap per half line', 20000, 10, 10000000, 1000, 0, 'A half line stops after this many nodes.'),
    # walk experiments that were tried and backed out; all off at these values
    ('hill', 'Walk experiments', 'Hill rule (1 = on)', 0, 0, 1, 1, 0, 'Band = the hill of the across-profile around the anchor instead of the whole connected set. Backed out 2026-10-17.'),
    ('valley', 'Walk experiments', 'Hill: valley rise', 1.3, 1, 5, 0.05, 2, 'The profile must rise this much above its running minimum to count as a second hill.'),
    ('core', 'Walk experiments', 'Hill: core share', 0.0, 0, 1, 0.05, 2, 'The node is read from the dots whose profile bin is at least this share of the peak.'),
    ('grow', 'Walk experiments', 'Hill: width clip (x expected)', 99, 1, 99, 1, 0, 'A hill wider than this many times the band so far is clipped to the expected width.'),
    ('ridge', 'Walk experiments', 'Hill: ridge node (1 = on)', 0, 0, 1, 1, 0, 'The node is read from the peak bin and its neighbours only.'),
    ('lead', 'Walk experiments', 'Lead cut', 0.0, 0, 5, 0.1, 1, 'Weak-dot filtering inside the band: dots below lead x (1 - median/top) of the top strength carry no connectivity. 0 = off. Backed out.'),
    ('leadmin', 'Walk experiments', 'Lead: min strong dots', 8, 1, 100, 1, 0, 'Too few strong dots to carry a band on their own: keep all.'),
    ('fringe', 'Walk experiments', 'Lead: fringe (1 = on)', 1, 0, 1, 1, 0, 'Weak dots touching the strong core join it.'),
    ('sparse', 'Walk experiments', 'Sparse cut (x median spacing)', 0.0, 0, 10, 0.1, 1, 'Drop slice dots whose nearest neighbour is further than this times the median spacing. 0 = off.'),
    ('lookn', 'Walk experiments', 'Look-ahead: older trend (nodes)', 0, 0, 100, 1, 0, 'Probe the trend of the nodes before the last few. 0 = off.'),
    ('astr', 'Walk experiments', 'Anchor strength exponent', 0.0, 0, 5, 0.1, 1, 'Exponent on relative strength in the anchor score. 0 = alignment only.'),
    # refine stage (refine.py)
    ('nbh', 'Refine', 'Bend neighbourhood (x area)', 2.0, 0.5, 10, 0.1, 1, 'The bend at a fine cube is fitted to all cubes within this many times the area; weights exp(-(d/area)^2).'),
    ('rungstep', 'Refine', 'Rung spacing (coarse cubes)', 0.5, 0.05, 5, 0.05, 2, 'The bend is read only on rungs across the line, one every this many coarse cubes. One per step broke the walk; half is the spacing found.'),
    ('rungw', 'Refine', 'Rung thickness (fine cubes)', 0.5, 0.1, 10, 0.1, 1, 'A fine cube belongs to a rung when it lies within this many fine cubes of it along the line.'),
    ('cresth', 'Refine', 'Crest read distance (x unit)', 0.5, 0.05, 5, 0.05, 2, 'The crest flag reads the bend this far either side across the line.'),
    ('tuberes', 'Refine', 'Tube distance resample (fine cubes)', 0.5, 0.05, 5, 0.05, 2, 'The line is resampled at this spacing before the distance of each cube to it is measured.'),
    ('avoidk', 'Refine', 'End rule reach (x area)', 2.0, 0, 10, 0.1, 1, 'Fine dots within this many times the area of the lines our ends sit on are removed before the fine walk.'),
    ('avoidres', 'Refine', 'End rule resample (mm)', 0.02, 0.001, 1, 0.005, 3, 'Those lines are resampled at this spacing for the distance measure.'),
    ('eslack', 'Refine', 'End slack (coarse cubes)', 1.0, 0, 10, 0.1, 1, 'Another line counts as the line an end sits on when it has a node within this many coarse cubes of the end point.'),
    ('funit', 'Refine', 'Fine walk unit (mm)', 0.087, 0.001, 5, 0.001, 3, 'The unit the fine walk reads its cube-sized rules in (the coarse cube, so they keep their meaning).'),
    # the fine walk's own settings (refine.WALK); the hidden walk constants above apply to it too
    ('fstep', 'Fine walk', 'Step (mm)', 0.087, 0.001, 2, 0.01, 3, 'Node spacing of the fine walk.'),
    ('fslice', 'Fine walk', 'Slice half-length (mm)', 0.7, 0.05, 5, 0.05, 2, 'Half-length of the fine slice across the band.'),
    ('freach', 'Fine walk', 'Reach (mm)', 0.087, 0.001, 3, 0.01, 3, 'How far ahead the next fine anchor may lie. One step: the walk never leaves the band.'),
    ('fpull', 'Fine walk', 'Trend pull', 1.0, 0, 20, 0.1, 1, "Trend's say on each fine node; 1 = equal to the dots."),
    ('flook', 'Fine walk', 'Look ahead (mm)', 0.45, 0, 2, 0.05, 2, 'Probe length of the fine walk.'),
    ('flevp', 'Fine walk', 'Strength power', 2.0, 0, 20, 0.5, 1, "A fine dot's say is (bend / standard) to this power: 2 = bend squared, every dot counts."),
    ('fminlen', 'Fine walk', 'Minimum line length (nodes)', 8, 1, 500, 1, 0, 'Fine lines shorter than this are dropped.'),
    # display and picking (app.py)
    ('dotlift', 'Display', 'Dot lift (cubes)', 0.5, 0, 5, 0.1, 1, 'Dots are drawn this many cubes out along the facing so dots on convex edges are not buried.'),
    ('snaplift', 'Display', 'Line lift on the surface (mm)', 0.01, 0, 1, 0.005, 3, 'Lines drawn on the surface are lifted this much along the face normal so they are not hidden inside it.'),
    ('picktol', 'Display', 'Pick tolerance (window fraction)', 0.006, 0.0005, 0.05, 0.001, 4, 'A line within this fraction of the window diagonal of the cursor is picked (about 7 px at 1200 px).'),
    ('pickdepth', 'Display', 'Pick depth (mm)', 0.5, 0, 10, 0.1, 1, 'A line this far behind the surface under the cursor is hidden by it and not picked.'),
]
TUNE_DEFAULTS = {t[0]: t[3] for t in TUNE}


class TuneDialog(QtWidgets.QDialog):
    """All the pipeline's inner constants as handles, grouped by stage. Edits land in app.tune at once; the next stage run reads them."""
    def __init__(self, app):
        super().__init__(app); self.app = app; self.setWindowTitle('Tuning: every constant of the pipeline'); self.resize(560, 900)
        lay = QtWidgets.QVBoxLayout(self); scroll = QtWidgets.QScrollArea(); scroll.setWidgetResizable(True); lay.addWidget(scroll, 1)
        body = QtWidgets.QWidget(); scroll.setWidget(body); vb = QtWidgets.QVBoxLayout(body); self.spins = {}
        note = QtWidgets.QLabel('Defaults are the values the code ran with before these handles existed. A handle that breaks the lines somewhere in its range is a missing mechanism, not a setting.'); note.setWordWrap(True); vb.addWidget(note)
        groups = {}
        for key, grp, label, default, lo, hi, step, dec, tip in TUNE:
            if grp not in groups: g = QtWidgets.QGroupBox(grp); groups[grp] = QtWidgets.QFormLayout(g); vb.addWidget(g)
            sp = QtWidgets.QDoubleSpinBox(); sp.setRange(lo, hi); sp.setSingleStep(step); sp.setDecimals(dec); sp.setValue(app.tune[key]); sp.setToolTip(tip)
            sp.valueChanged.connect(lambda v, k=key: app.tune.__setitem__(k, v)); groups[grp].addRow(label, sp); groups[grp].labelForField(sp).setToolTip(tip); self.spins[key] = sp
        vb.addStretch(1)
        btns = QtWidgets.QHBoxLayout(); lay.addLayout(btns)
        b = QtWidgets.QPushButton('Reset all to defaults'); b.clicked.connect(self.reset); btns.addWidget(b)
        c = QtWidgets.QPushButton('Close'); c.clicked.connect(self.close); btns.addWidget(c)
    def reset(self):
        for k, sp in self.spins.items(): sp.setValue(TUNE_DEFAULTS[k])


class App(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__(); self.setWindowTitle('Segmentation tool'); self.resize(1500, 950)
        self.P = self.F = None; self.path = None; self.res = None; self.curves = None; self.refined = None; self.mesh_t = None
        self.drawn = {}; self.act_mesh = None; self.hover = None; self.sel = None; self._press = None   # picking: drawn line actors, hovered and selected (kind, k)
        self.tune = dict(TUNE_DEFAULTS); self.tune_dlg = None   # every inner constant of the pipeline (TUNE), edited in the Tuning window
        root = QtWidgets.QWidget(); self.setCentralWidget(root); lay = QtWidgets.QHBoxLayout(root)
        self.plotter = QtInteractor(root); lay.addWidget(self.plotter.interactor, 1)
        self.plotter.iren.add_observer('MouseMoveEvent', self._on_move); self.plotter.iren.add_observer('LeftButtonPressEvent', self._on_press); self.plotter.iren.add_observer('LeftButtonReleaseEvent', self._on_release)
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
        self.btn_load_dots = QtWidgets.QPushButton('Load saved dots…'); self.btn_load_dots.clicked.connect(self.load_dots); self.btn_load_dots.setEnabled(False); self.btn_load_dots.setToolTip('Dot sets are saved to out/dots/ every time Find dots finishes; opening a model loads its last set by itself.'); f1.addRow(self.btn_load_dots)
        side.addWidget(g1)

        g2 = QtWidgets.QGroupBox('2. Build lines'); f2 = QtWidgets.QFormLayout(g2)
        self.cb_method = QtWidgets.QComboBox(); self.cb_method.addItems(['Anchored walk', 'Thick strokes (old)']); row(f2, 'Method', self.cb_method, TIPS['method'])
        self.w_walk = QtWidgets.QWidget(); fw = QtWidgets.QFormLayout(self.w_walk); fw.setContentsMargins(0, 0, 0, 0); f2.addRow(self.w_walk)
        self.sp_step = spin(0.087, 0.01, 2, 0.01, 3); row(fw, 'Step (mm)', self.sp_step, TIPS['step'])
        self.sp_slice = spin(0.7, 0.05, 5, 0.05, 2); row(fw, 'Slice half-length (mm)', self.sp_slice, TIPS['slice'])
        self.sp_reach = spin(0.22, 0.02, 3, 0.02, 2); row(fw, 'Reach (mm)', self.sp_reach, TIPS['reach'])
        self.sp_pull = spin(10.0, 0, 20, 0.1, 1); row(fw, 'Trend pull', self.sp_pull, TIPS['pull'])
        self.sp_look = spin(0.35, 0, 2, 0.05, 2); row(fw, 'Look ahead (mm)', self.sp_look, TIPS['look'])
        self.sp_levp = spin(6.0, 0, 20, 0.5, 1); row(fw, 'Strength power', self.sp_levp, TIPS['levp'])
        self.sp_wminlen = spin(8, 1, 500, 1); row(fw, 'Minimum line length (nodes)', self.sp_wminlen, TIPS['wminlen'])
        self.w_old = QtWidgets.QWidget(); f2o = QtWidgets.QFormLayout(self.w_old); f2o.setContentsMargins(0, 0, 0, 0); f2.addRow(self.w_old); f2, f2w = f2o, f2
        self.sp_W = spin(0.4, 0.02, 5, 0.05, 2); row(f2, 'Stroke half-width (mm)', self.sp_W, TIPS['W'])
        self.sp_L = spin(1.0, 0.1, 10, 0.1, 2); row(f2, 'Stroke stretch (mm)', self.sp_L, TIPS['L'])
        self.sp_K = spin(0.26, 0.02, 5, 0.05, 2); row(f2, 'Min dot row ahead (mm)', self.sp_K, TIPS['K'])
        self.sp_align = spin(3, 1, 10, 0.5, 1); row(f2, 'Alignment (long/wide)', self.sp_align, TIPS['align'])
        self.sp_agree = spin(30, 1, 90, 5); row(f2, 'Max turn per step (deg)', self.sp_agree, TIPS['agree'])
        self.sp_minlen = spin(10, 1, 500, 5); row(f2, 'Minimum stroke length (nodes)', self.sp_minlen, TIPS['minlen'])
        self.sp_G = spin(1.0, 0, 10, 0.1, 2); row(f2, 'Join reach (mm)', self.sp_G, TIPS['G'])
        self.btn_lines = QtWidgets.QPushButton('Build lines (from dots)'); self.btn_lines.clicked.connect(self.build_lines); self.btn_lines.setEnabled(False); f2w.addRow(self.btn_lines)
        self.cb_method.currentIndexChanged.connect(self._method_changed); self._method_changed()
        side.addWidget(g2)

        g3 = QtWidgets.QGroupBox('3. Refine lines'); f3 = QtWidgets.QFormLayout(g3)
        self.lbl_pick = QtWidgets.QLabel('no line selected (hover a line, click to select)'); self.lbl_pick.setWordWrap(True); self.lbl_pick.setToolTip(TIPS['pick']); f3.addRow(self.lbl_pick)
        self.sp_finer = spin(3, 1, 20, 1); row(f3, 'Finer than the coarse cube (times)', self.sp_finer, TIPS['finer'])
        self.sp_radius = spin(0.3, 0.05, 2, 0.05, 2); row(f3, 'Tube radius (mm)', self.sp_radius, TIPS['radius'])
        self.sp_area = spin(0.087, 0.01, 1, 0.01, 3); row(f3, 'Bend area (mm)', self.sp_area, TIPS['area'])
        self.sp_rcut = spin(17, 1, 180, 1); row(f3, 'Cut (deg per coarse cube)', self.sp_rcut, TIPS['rcut'])
        self.btn_refine1 = QtWidgets.QPushButton('Refine one line (the selected one)'); self.btn_refine1.clicked.connect(lambda: self.refine_lines(True)); self.btn_refine1.setEnabled(False); f3.addRow(self.btn_refine1)
        self.btn_refine = QtWidgets.QPushButton('Refine all lines (a few seconds per line)'); self.btn_refine.clicked.connect(lambda: self.refine_lines(False)); self.btn_refine.setEnabled(False); f3.addRow(self.btn_refine)
        side.addWidget(g3)

        self.btn_tune = QtWidgets.QPushButton('Tuning: all constants…'); self.btn_tune.clicked.connect(self.open_tune); self.btn_tune.setToolTip('Every constant of the pipeline not on this panel, as a handle. The main panel keeps the handles found to be sensible.'); side.addWidget(self.btn_tune)
        self.chk_flatshade = QtWidgets.QCheckBox('Flat shading (creases visible)'); self.chk_flatshade.setChecked(False); self.chk_flatshade.setToolTip('Each triangle lit by its own facing, as the print has it. Off: smooth shading, which blends every fold over the neighbouring triangles and hides faint creases.'); self.chk_flatshade.stateChanged.connect(self.redraw); side.addWidget(self.chk_flatshade)
        self.chk_dots = QtWidgets.QCheckBox('Show dots'); self.chk_dots.setChecked(True); self.chk_dots.stateChanged.connect(self.redraw); side.addWidget(self.chk_dots)
        self.chk_curves = QtWidgets.QCheckBox('Show walk lines (blue)'); self.chk_curves.setChecked(True); self.chk_curves.stateChanged.connect(self.redraw); side.addWidget(self.chk_curves)
        self.chk_refined = QtWidgets.QCheckBox('Show refined lines (red)'); self.chk_refined.setChecked(True); self.chk_refined.stateChanged.connect(self.redraw); side.addWidget(self.chk_refined)
        self.chk_snap = QtWidgets.QCheckBox('Draw lines on the surface (a few seconds)'); self.chk_snap.setChecked(True); self.chk_snap.setToolTip(TIPS['snap']); self.chk_snap.stateChanged.connect(self.redraw); side.addWidget(self.chk_snap)
        self.btn_save = QtWidgets.QPushButton('Save lines…'); self.btn_save.clicked.connect(self.save); self.btn_save.setEnabled(False); side.addWidget(self.btn_save)
        self.btn_dump = QtWidgets.QPushButton('Dump for debugging'); self.btn_dump.clicked.connect(self.dump); self.btn_dump.setToolTip('Writes dots, lines, log, settings and a screenshot of the viewport to out/debug/<time>/'); side.addWidget(self.btn_dump)
        self.txt = QtWidgets.QPlainTextEdit(); self.txt.setReadOnly(True); self.txt.setMaximumWidth(380); side.addWidget(self.txt, 1)
        self.plotter.set_background('white'); self.plotter.add_axes()

    def logmsg(self, s): self.txt.appendPlainText(s)
    def _method_changed(self, *_):
        walk = self.cb_method.currentIndex() == 0; self.w_walk.setVisible(walk); self.w_old.setVisible(not walk)
    def busy(self, b):
        self.btn_dots.setEnabled(not b and self.P is not None); self.btn_load_dots.setEnabled(not b and self.P is not None); self.btn_lines.setEnabled(not b and self.res is not None); self.btn_refine1.setEnabled(not b and bool(self.curves)); self.btn_refine.setEnabled(not b and bool(self.curves))

    def load(self):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, 'Open STL', '', 'Meshes (*.stl *.obj *.ply *.3mf)')
        if not path: return
        try: self.P, self.F = load_mesh(path)
        except Exception as ex: QtWidgets.QMessageBox.critical(self, 'Load failed', str(ex)); return
        self.path = path; self.lbl_file.setText(os.path.basename(path)); self.res = self.curves = self.refined = None; self.mesh_t = None; self.btn_save.setEnabled(False)
        h = float((self.P.max(0)-self.P.min(0)).max())
        e = np.linalg.norm(self.P[self.F[:, 1]]-self.P[self.F[:, 0]], axis=1)
        self.logmsg(f'loaded {os.path.basename(path)}: {len(self.F)} faces, height {h:.1f} mm, typical triangle edge {np.median(e):.3f} mm')
        p = last_dots(path)
        if p: self._use_dots(p)
        self.busy(False); self.redraw()

    def _use_dots(self, p):
        try: res, cube, flat, model = load_dots(p)
        except Exception as ex: self.logmsg(f'could not load dots from {p}: {ex}'); return
        self.res = res; self.curves = self.refined = None; self.btn_save.setEnabled(False); self.sp_cube.setValue(cube); self.sp_cut.setValue(res['cut']); self.chk_flat.setChecked(flat)
        self.logmsg(f'loaded {len(res["idx"])} dots from {os.path.basename(p)} (cube {cube:g} mm, cut {res["cut"]:g}, {"flat" if flat else "blended"})')

    def load_dots(self):
        os.makedirs(DOTS_DIR, exist_ok=True)
        p, _ = QtWidgets.QFileDialog.getOpenFileName(self, 'Load saved dots', DOTS_DIR, 'Dot sets (*.npz)')
        if not p: return
        self._use_dots(p); self.busy(False); self.redraw()

    def _start(self, fn, on_done):
        self.busy(True); self.worker = Worker(fn); self.worker.log.connect(self.logmsg); self.worker.done.connect(on_done); self.worker.failed.connect(self.failed); self.worker.start()
    def failed(self, msg): self.logmsg(msg); self.busy(False)

    def find_dots(self):
        size, cut, flat = self.sp_cube.value(), self.sp_cut.value(), self.chk_flat.isChecked()
        h = float((self.P.max(0)-self.P.min(0)).max()); cells = max(1, int(round(h/size)))
        self.logmsg(f'--- find dots: cube={size} mm ({cells} cells over {h:.1f} mm), cut={cut}, flat={flat}'); t0 = time.time()
        self._start(lambda log: creases.find_dots(self.P, self.F, cells, cut, 0.0, flat, log, height=h), lambda r: self._dots_done(r, t0))
    def _dots_done(self, r, t0):
        self.res = r; self.curves = self.refined = None; self.logmsg(f'dots done in {time.time()-t0:.1f} s')
        try: self.logmsg(f'saved dots to {save_dots(r, self.path, self.sp_cube.value(), self.sp_cut.value(), self.chk_flat.isChecked())}')
        except Exception: self.logmsg(traceback.format_exc())
        self.busy(False); self.redraw()

    def stroke_params(self):
        return dict(W=self.sp_W.value(), L=self.sp_L.value(), K=self.sp_K.value(), align=self.sp_align.value(), agree=self.sp_agree.value(), minlen=int(self.sp_minlen.value()))
    def walk_tune(self):
        """the walk's hidden constants (TUNE groups Walk and Walk experiments) as trace4.walk keywords"""
        return {k: self.tune[k] for k, g, *_ in TUNE if g in ('Walk', 'Walk experiments')}
    def walk_params(self):
        return dict(**self.walk_tune(), step=self.sp_step.value(), slice=self.sp_slice.value(), reach=self.sp_reach.value(), pull=self.sp_pull.value(), look=self.sp_look.value(), levp=self.sp_levp.value(), minlen=int(self.sp_wminlen.value()))
    def refine_params(self):
        """main handles plus the stage's inner constants (tune) and the fine walk's settings (walk: hidden walk constants + the Fine walk group)"""
        t = self.tune; walk = {**self.walk_tune(), 'step': t['fstep'], 'slice': t['fslice'], 'reach': t['freach'], 'pull': t['fpull'], 'look': t['flook'], 'levp': t['flevp'], 'minlen': int(t['fminlen']), 'unit': t['funit']}   # the fine walk's own unit replaces the walk group's
        tune = dict(nbh=t['nbh'], rungstep=t['rungstep'], rungw=t['rungw'], cresth=t['cresth'], tuberes=t['tuberes'], avoidk=t['avoidk'], avoidres=t['avoidres'], slack=t['eslack'])
        return dict(finer=int(self.sp_finer.value()), radius=self.sp_radius.value(), area=self.sp_area.value(), cut=self.sp_rcut.value(), unit=t['funit'], walk=walk, tune=tune)
    def all_params(self):
        return dict(method=self.cb_method.currentText(), **({} if self.cb_method.currentIndex() else self.walk_params()), **({'join_mm': self.sp_G.value(), **self.stroke_params()} if self.cb_method.currentIndex() else {}), tune=dict(self.tune))
    def open_tune(self):
        if self.tune_dlg is None: self.tune_dlg = TuneDialog(self)
        self.tune_dlg.show(); self.tune_dlg.raise_()

    def build_lines(self):
        t0 = time.time()
        if self.cb_method.currentIndex() == 0:
            wp = self.walk_params(); self.logmsg('--- build lines (anchored walk): ' + ', '.join(f'{k}={v}' for k, v in wp.items()))
            self._start(lambda log: trace4.run(self.res, log=log, **wp), lambda r: self._lines_done(r, t0)); return
        sp = self.stroke_params(); self.logmsg('--- build lines (thick strokes): ' + ', '.join(f'{k}={v}' for k, v in sp.items()))
        G = self.sp_G.value()
        self._start(lambda log: (lambda cv: creases.join(cv, G, log=log) if G > 0 else cv)(creases.strokes(self.res, log=log, **sp)), lambda r: self._lines_done(r, t0))
    def _lines_done(self, r, t0):
        self.curves = r; self.refined = None; self.logmsg(f'lines done in {time.time()-t0:.1f} s'); self.busy(False); self.btn_save.setEnabled(True); self.redraw()

    def refine_lines(self, one):
        """stage 3: refine.refine_all on one line (nearest the point) or all; result kept in self.refined (dict of lists)"""
        if not self.curves: return
        rp = self.refine_params(); only = None
        if one:
            if self.sel is None: self.logmsg('select a line first: hover it in the viewport and click'); return
            only = self._walk_index(self.sel)
        h = float((self.P.max(0)-self.P.min(0)).max()); cube = self.sp_cube.value(); flat = self.chk_flat.isChecked()
        self.logmsg(f'--- refine {"line "+str(only) if one else "all lines"}: ' + ', '.join(f'{k}={v}' for k, v in rp.items() if k not in ('walk', 'tune')) + f'; fine walk {rp["walk"]}; tune {rp["tune"]}'); t0 = time.time()
        self._start(lambda log: refine.refine_all(self.curves, self.P, self.F, cube, h, log=log, only=only, flat=flat, **rp), lambda r: self._refine_done(r, t0))
    def _refine_done(self, r, t0):
        if self.sel is not None and self.sel[0] == 'refined': self.sel = ('walk', self._walk_index(self.sel))   # the refined list is about to change
        if self.refined is not None and len(set(r['index'])) == 1:                  # one refined line joins those already refined
            old = [(i, c) for i, c in zip(self.refined['index'], self.refined['lines']) if i not in r['index']]
            self.refined = dict(index=[i for i, _ in old]+r['index'], lines=[c for _, c in old]+r['lines'], tubes=self.refined['tubes']+r['tubes'])
        else: self.refined = r
        self.logmsg(f'refined in {time.time()-t0:.1f} s; {len(self.refined["index"])} refined lines held'); self.busy(False); self.redraw()

    def _draw_lines(self, curves, color, width, kind=None, index=None):
        """draws the curves as one actor, one polyline cell per curve. kind: registers the actor for picking, with index[cell] the
        curve's position in its list (walk: self.curves, refined: self.refined['lines'])"""
        pts = np.vstack(curves)
        if self.chk_snap.isChecked():
            if self.mesh_t is None: self.mesh_t = trimesh.Trimesh(self.P, self.F, process=False)
            Cq, dq, tq = trimesh.proximity.closest_point(self.mesh_t, pts); pts = Cq+self.tune['snaplift']*self.mesh_t.face_normals[tq]
        n = [len(c) for c in curves]; off = np.cumsum([0]+n[:-1])
        cells = np.concatenate([np.r_[k, np.arange(o, o+k)] for k, o in zip(n, off)])
        a = self.plotter.add_mesh(pv.PolyData(pts, lines=cells), color=color, line_width=width)
        if kind: self.drawn[kind] = (a, pts, n, off, list(range(len(curves))) if index is None else index)

    def _line_pts(self, key):
        """the drawn (snapped) points of line key=(kind, k)"""
        a, pts, n, off, index = self.drawn[key[0]]; c = index.index(key[1]); return pts[off[c]:off[c]+n[c]]
    def _walk_index(self, key):
        """the walk line behind key: itself, or the one a refined line came from"""
        return key[1] if key[0] == 'walk' else self.refined['index'][key[1]]
    def _mark(self, name, key, color, width):
        """draws (or removes, key=None) the highlight actor `name` over line key"""
        self.plotter.remove_actor(name, render=False)
        if key is not None and key[0] in self.drawn:
            p = self._line_pts(key); self.plotter.add_mesh(pv.PolyData(p, lines=np.r_[len(p), np.arange(len(p))]), color=color, line_width=width, name=name, reset_camera=False)
        self.plotter.render()
    def _describe(self, key):
        if key is None: return 'no line selected (hover a line, click to select)'
        i = self._walk_index(key); n = len(self.curves[i]); kind = 'refined ' if key[0] == 'refined' else ''
        return f'selected: {kind}line {i} ({n} walk nodes{", "+str(len(self.refined["lines"][key[1]]))+" fine nodes" if key[0] == "refined" else ""})'
    def _pick(self):
        if not self.drawn: return None
        x, y = self.plotter.iren.get_event_position(); r = pick_line(self.plotter, self.act_mesh, self.drawn, x, y, self.tune['picktol'], self.tune['pickdepth'])
        return None if r is None else (r[0], r[1])
    def _on_move(self, *_):
        key = self._pick()
        if key != self.hover: self.hover = key; self._mark('hover', None if key == self.sel else key, 'yellow', 7)
    def _on_press(self, *_): self._press = self.plotter.iren.get_event_position()
    def _on_release(self, *_):
        pos = self.plotter.iren.get_event_position()
        if self._press is None or abs(pos[0]-self._press[0]) > 3 or abs(pos[1]-self._press[1]) > 3: return   # a drag rotated the view
        key = self._pick(); self.sel = key; self._mark('hover', None, 'yellow', 7); self._mark('sel', key, 'lime', 7); self.lbl_pick.setText(self._describe(key))

    def redraw(self):
        if self.P is None: return
        cam = self.plotter.camera_position; had = self.res is not None; self.plotter.renderer.clear_actors(); self.plotter.remove_scalar_bar() if self.plotter.scalar_bars else None; self.drawn = {}; self.hover = None   # not clear(): that drops the lights too
        mesh = pv.PolyData(self.P, np.hstack([np.full((len(self.F), 1), 3), self.F]).ravel())
        self.act_mesh = self.plotter.add_mesh(mesh, color='lightgrey', smooth_shading=not self.chk_flatshade.isChecked())
        if self.res is not None and self.chk_dots.isChecked():
            C = self.res['lo']+(self.res['idx'][self.res['line']]+0.5)*self.res['size']
            C = C+self.tune['dotlift']*self.res['size']*self.res['N'][self.res['line']]   # out along the facing, so dots on convex edges are not buried
            d = pv.PolyData(C); d['fold (deg)'] = self.res['turn'][self.res['line']]   # dot strength: how sharply the surface folds there
            hi = float(np.percentile(d['fold (deg)'], 99))
            self.plotter.add_mesh(d, scalars='fold (deg)', cmap='plasma', clim=(float(self.res['cut']), hi), point_size=5, render_points_as_spheres=True,
                                  scalar_bar_args=dict(title='dot strength (deg)', vertical=True, position_x=0.88, position_y=0.3, height=0.4, width=0.06))
        if self.curves and self.chk_curves.isChecked(): self._draw_lines(self.curves, 'blue', 3, 'walk')
        if self.refined is not None and self.refined['lines'] and self.chk_refined.isChecked():
            keep = [k for k, c in enumerate(self.refined['lines']) if len(c) > 1]
            if keep: self._draw_lines([self.refined['lines'][k] for k in keep], 'red', 4, 'refined', keep)   # refined lines over the walk's lines
        if self.sel is not None and (self.sel[0] not in self.drawn or self.sel[1] not in self.drawn[self.sel[0]][4]): self.sel = None   # the selected line is gone
        if cam is not None and had: self.plotter.camera_position = cam
        else: self.plotter.reset_camera()
        self._mark('sel', self.sel, 'lime', 7); self.lbl_pick.setText(self._describe(self.sel))

    def save(self):
        if not self.curves: return
        out = os.path.splitext(self.path)[0]+'_lines.npz'
        prm = dict(cube_mm=self.sp_cube.value(), cut=self.sp_cut.value(), flat=self.chk_flat.isChecked(), **self.all_params())
        cv = list(self.curves)
        if self.refined is not None:                                                  # refined lines replace their walk lines (a line refined into
            done = sorted(set(self.refined['index']))                                 # several pieces gives several lines; one into none is dropped)
            cv = [c for i, c in enumerate(cv) if i not in done]+[c for c in self.refined['lines'] if len(c) > 1]
            prm['refined'] = len(done); prm.update({'refine_'+k: v for k, v in self.refine_params().items() if k not in ('walk', 'tune')})
        np.savez(out, n=np.array([len(c) for c in cv]), pts=np.vstack(cv), params=str(prm)); self.logmsg(f'saved {out} ({prm.get("refined", 0)} refined lines)')

    def dump(self):
        """Everything on screen, to out/debug/<time>/: dots.npz (cube centres, drawn positions, grid), lines.npz (points,
        run lengths), log.txt, settings.json (parameters, model, camera) and viewport.png. Works with any stage done."""
        d = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out', 'debug', time.strftime('%Y%m%d_%H%M%S')); os.makedirs(d, exist_ok=True)
        try:
            prm = dict(model=self.path, faces=None if self.F is None else int(len(self.F)), cube_mm=self.sp_cube.value(), cut=self.sp_cut.value(), flat=self.chk_flat.isChecked(),
                       **self.all_params(), show_dots=self.chk_dots.isChecked(), flat_shading=self.chk_flatshade.isChecked(), show_lines=self.chk_curves.isChecked(), show_refined=self.chk_refined.isChecked(), snap=self.chk_snap.isChecked(),
                       camera=[list(map(float, v)) for v in self.plotter.camera_position], dots=0, lines=0)
            if self.res is not None:
                r = self.res; C = r['lo']+(r['idx'][r['line']]+0.5)*r['size']; D = C+self.tune['dotlift']*r['size']*r['N'][r['line']]
                np.savez(os.path.join(d, 'dots.npz'), centres=C, drawn=D, idx=r['idx'][r['line']], N=r['N'][r['line']], turn=r['turn'][r['line']], lo=r['lo'], size=r['size'])
                prm['dots'] = int(len(C))
            if self.curves:
                np.savez(os.path.join(d, 'lines.npz'), n=np.array([len(c) for c in self.curves]), pts=np.vstack(self.curves)); prm['lines'] = len(self.curves)
            if self.refined is not None and self.refined['lines']:
                r_ = self.refined; T = r_['tubes']
                np.savez(os.path.join(d, 'refined.npz'), n=np.array([len(c) for c in r_['lines']]), pts=np.vstack(r_['lines']), index=np.array(r_['index']),
                         dots=np.vstack([t['dots'] for t in T]), turn=np.concatenate([t['turn'] for t in T]), size=T[0]['size'])
                prm['refined'] = len(set(r_['index'])); prm.update({'refine_'+k: v for k, v in self.refine_params().items() if k not in ('walk', 'tune')})
            with open(os.path.join(d, 'log.txt'), 'w', encoding='utf-8') as f: f.write(self.txt.toPlainText())
            import json
            with open(os.path.join(d, 'settings.json'), 'w', encoding='utf-8') as f: json.dump(prm, f, indent=1)
            self.plotter.screenshot(os.path.join(d, 'viewport.png'))
            self.logmsg(f'dumped {prm["dots"]} dots, {prm["lines"]} lines, log, settings, screenshot to {d}')
        except Exception: self.logmsg(traceback.format_exc())


if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv); w = App(); w.show(); sys.exit(app.exec())
