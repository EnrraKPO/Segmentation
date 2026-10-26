"""Painting tool. Loads a model and the lines the line app saved, finds the regions the lines enclose (regions.py: the
lines are the exact borders, the mesh is only the canvas and is never changed), pre-colours them from a 9-colour palette
and lets one bucket stroke recolour a region. Export: a Bambu Studio project (bambu3mf.py) with the untouched mesh and a
paint note per triangle, the way Bambu Studio itself stores a brush stroke.
Run: python paint.py (from this folder)."""
import os, sys, time, traceback, pickle
import numpy as np
from PySide6 import QtCore, QtWidgets
import pyvista as pv
from pyvistaqt import QtInteractor
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import regions, bambu3mf
from linesio import read_lines, load_mesh

LOG_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out', 'log')
TIPS = dict(
    set="Which line set bounds the regions: the walk lines (with the merge's shared nodes), the refined pieces, or the joined set. Only the sets the file holds are offered.",
    min_area="Regions smaller than this are absorbed into their largest neighbour (smallest first). Below a nozzle width squared (0.4 mm -> 0.16 mm2) a region cannot be printed as its own colour anyway; 0 keeps every region the lines make, slivers included.",
    depth="How many triangles a straight stretch between two consecutive nodes may cross (the search unfolds that many). Nodes are closer than a triangle edge nearly everywhere, so a few suffice; a stretch that finds no strip is bisected.",
    fine="Export: a triangle crossed by a border is split in four, again and again where the border still runs through a child, until the child is smaller than this (mm). Bambu Studio redoes the same splitting in memory; the mesh in the file is untouched.",
    edepth="Export: at most this many splittings of one triangle (4^depth children at the worst).",
    lines="Draw the lines (through their surface points) on top of the colours.",
    bucket="Click a region to give it the active colour. Dragging rotates the view as usual.",
    loose="A loose end is a line end that sits on no other line's node: a loop through it is open, and the regions on its two sides are one. An undrawn stretch is a gap the region finder could not draw across (two consecutive nodes projected onto different sheets of the surface). Both are where a colour leaks.",
)


def open_session_log():
    os.makedirs(LOG_DIR, exist_ok=True); path = os.path.join(LOG_DIR, 'paint_'+time.strftime('%Y%m%d_%H%M%S')+'.txt'); f = open(path, 'a', encoding='utf-8')
    def write(s): f.write(s+'\n'); f.flush()
    return write, path


class Worker(QtCore.QThread):
    log = QtCore.Signal(str); done = QtCore.Signal(object); failed = QtCore.Signal(str)
    def __init__(self, fn): super().__init__(); self.fn = fn
    def run(self):
        try: self.done.emit(self.fn(self.log.emit))
        except Exception: self.failed.emit(traceback.format_exc())


class App(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__(); self.setWindowTitle('Painting tool'); self.resize(1500, 950)
        self.P = self.F = None; self.path = None; self.st = None; self.res = None; self.mesh = None; self.act = None; self.locator = None; self._press = None; self.hover = None; self.worker = None
        try: self.file_log, _ = open_session_log()
        except Exception: self.file_log = None
        root = QtWidgets.QWidget(); self.setCentralWidget(root); lay = QtWidgets.QHBoxLayout(root)
        self.plotter = QtInteractor(root); lay.addWidget(self.plotter.interactor, 1)
        self.plotter.iren.add_observer('LeftButtonPressEvent', self._on_press); self.plotter.iren.add_observer('LeftButtonReleaseEvent', self._on_release)
        panel = QtWidgets.QWidget(); panel.setFixedWidth(360); side = QtWidgets.QVBoxLayout(panel); lay.addWidget(panel)
        def spin(v, lo, hi, step, dec=0):
            s = QtWidgets.QDoubleSpinBox() if dec else QtWidgets.QSpinBox(); s.setRange(lo, hi); s.setSingleStep(step)
            if dec: s.setDecimals(dec)
            s.setValue(v); return s
        def row(form, label, widget, tip):
            lab = QtWidgets.QLabel(label); lab.setToolTip(tip); widget.setToolTip(tip); form.addRow(lab, widget)
        g0 = QtWidgets.QGroupBox('Input'); f0 = QtWidgets.QFormLayout(g0)
        self.btn_model = QtWidgets.QPushButton('Load model…'); self.btn_model.clicked.connect(self.load_model); f0.addRow(self.btn_model)
        self.lbl_model = QtWidgets.QLabel('no model'); self.lbl_model.setWordWrap(True); f0.addRow(self.lbl_model)
        self.btn_lines = QtWidgets.QPushButton('Load lines…'); self.btn_lines.clicked.connect(self.load_lines); self.btn_lines.setEnabled(False); f0.addRow(self.btn_lines)
        self.lbl_lines = QtWidgets.QLabel('no lines'); self.lbl_lines.setWordWrap(True); f0.addRow(self.lbl_lines)
        self.cb_set = QtWidgets.QComboBox(); row(f0, 'Line set', self.cb_set, TIPS['set'])
        side.addWidget(g0)
        g1 = QtWidgets.QGroupBox('Regions'); f1 = QtWidgets.QFormLayout(g1)
        self.sp_min = spin(0.16, 0, 100, 0.01, 3); row(f1, 'Absorb regions below (mm²)', self.sp_min, TIPS['min_area'])
        self.sp_depth = spin(8, 2, 64, 1); row(f1, 'Strip search depth (triangles)', self.sp_depth, TIPS['depth'])
        self.btn_find = QtWidgets.QPushButton('Find regions'); self.btn_find.clicked.connect(self.find_regions); self.btn_find.setEnabled(False); f1.addRow(self.btn_find)
        self.lbl_regions = QtWidgets.QLabel(''); self.lbl_regions.setWordWrap(True); f1.addRow(self.lbl_regions)
        side.addWidget(g1)
        g2 = QtWidgets.QGroupBox('Paint (bucket: click a region)'); g2.setToolTip(TIPS['bucket']); f2 = QtWidgets.QVBoxLayout(g2)
        sw = QtWidgets.QHBoxLayout(); self.swatches = QtWidgets.QButtonGroup(self); self.swatches.setExclusive(True)
        for i, c in enumerate(regions.PALETTE):
            b = QtWidgets.QPushButton(); b.setCheckable(True); b.setFixedSize(30, 30); b.setToolTip(f'colour {i+1}: {c}')
            b.setStyleSheet(f'QPushButton {{ background: {c}; border: 2px solid #888; }} QPushButton:checked {{ border: 3px solid black; }}'); self.swatches.addButton(b, i); sw.addWidget(b)
        self.swatches.button(0).setChecked(True); f2.addLayout(sw)
        self.lbl_hover = QtWidgets.QLabel('click a region'); self.lbl_hover.setWordWrap(True); f2.addWidget(self.lbl_hover)
        self.chk_lines = QtWidgets.QCheckBox('Show lines'); self.chk_lines.setChecked(True); self.chk_lines.setToolTip(TIPS['lines']); self.chk_lines.stateChanged.connect(self.redraw); f2.addWidget(self.chk_lines)
        self.chk_loose = QtWidgets.QCheckBox('Show loose ends (red) and undrawn stretches (magenta)'); self.chk_loose.setChecked(True); self.chk_loose.setToolTip(TIPS['loose']); self.chk_loose.stateChanged.connect(self.redraw); f2.addWidget(self.chk_loose)
        side.addWidget(g2)
        g3 = QtWidgets.QGroupBox('Export'); f3 = QtWidgets.QFormLayout(g3)
        self.sp_fine = spin(0.05, 0.005, 5, 0.005, 3); row(f3, 'Finest split (mm)', self.sp_fine, TIPS['fine'])
        self.sp_edepth = spin(8, 1, 12, 1); row(f3, 'Most splits per triangle', self.sp_edepth, TIPS['edepth'])
        self.btn_export = QtWidgets.QPushButton('Export to Bambu Studio (3MF)…'); self.btn_export.clicked.connect(self.export); self.btn_export.setEnabled(False); f3.addRow(self.btn_export)
        side.addWidget(g3)
        self.btn_save = QtWidgets.QPushButton('Save painting…'); self.btn_save.clicked.connect(self.save); self.btn_save.setEnabled(False); side.addWidget(self.btn_save)
        self.btn_load = QtWidgets.QPushButton('Load painting…'); self.btn_load.clicked.connect(self.load_painting); self.btn_load.setEnabled(False); side.addWidget(self.btn_load)
        self.txt = QtWidgets.QPlainTextEdit(); self.txt.setReadOnly(True); side.addWidget(self.txt, 1)
        self.plotter.add_axes(); self.plotter.set_background('white')

    # ---------------------------------------------------------------- plumbing
    def logmsg(self, s):
        self.txt.appendPlainText(s)
        if self.file_log: self.file_log(s)
    def busy(self, b):
        has = self.P is not None; self.btn_model.setEnabled(not b); self.btn_lines.setEnabled(not b and has); self.btn_load.setEnabled(not b and has)
        self.btn_find.setEnabled(not b and has and self.st is not None); self.btn_export.setEnabled(not b and self.res is not None); self.btn_save.setEnabled(not b and self.res is not None)
    def _start(self, fn, on_done):
        self.busy(True); self.worker = Worker(fn); self.worker.log.connect(self.logmsg); self.worker.done.connect(on_done); self.worker.failed.connect(self.failed); self.worker.start()
    def failed(self, msg): self.logmsg(msg); self.busy(False)

    # ---------------------------------------------------------------- input
    def load_model(self):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, 'Open model', '', 'Meshes (*.stl *.obj *.ply *.3mf)')
        if not path: return
        try: self.P, self.F = load_mesh(path)
        except Exception as ex: QtWidgets.QMessageBox.critical(self, 'Load failed', str(ex)); return
        self.path = path; self.res = None; self.st = None; self.lbl_model.setText(os.path.basename(path)); self.lbl_lines.setText('no lines'); self.cb_set.clear()
        e = np.linalg.norm(self.P[self.F[:, 1]]-self.P[self.F[:, 0]], axis=1)
        self.logmsg(f'loaded {os.path.basename(path)}: {len(self.F)} triangles, typical edge {np.median(e):.3f} mm'); self.busy(False); self.redraw()
    def load_lines(self):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, 'Open lines', os.path.dirname(self.path) if self.path else '', 'Lines (*.npz)')
        if not path: return
        try: st = read_lines(path)
        except Exception as ex: QtWidgets.QMessageBox.critical(self, 'Load failed', str(ex)); return
        self.st = st; self.cb_set.clear(); self.cb_set.addItem(f'{st.get("flat", "walk lines")} ({len(st["walk"])})', 'walk')
        if st['refined'] is not None: self.cb_set.addItem(f'refined pieces ({len(st["refined"]["lines"])})', 'refined')
        if st['joined'] is not None: self.cb_set.addItem(f'joined lines ({len(st["joined"]["lines"])})', 'joined')
        self.lbl_lines.setText(os.path.basename(path)+f' ({st["note"]})'); self.logmsg(f'loaded lines from {path} ({st["note"]})'); self.busy(False)
    def current_lines(self):
        key = self.cb_set.currentData()
        return self.st['walk'] if key == 'walk' else self.st[key]['lines']

    # ---------------------------------------------------------------- regions
    def find_regions(self):
        lines = [L for L in self.current_lines() if len(L) > 1]; min_area = self.sp_min.value(); depth = int(self.sp_depth.value())
        self.logmsg(f'--- find regions: {self.cb_set.currentText()}, absorb below {min_area} mm2, depth {depth}'); t0 = time.time(); P, F = self.P, self.F
        self._start(lambda log: regions.run(P, F, lines, log=log, min_area=min_area, depth=depth), lambda r: self._regions_done(r, t0))
    def _regions_done(self, r, t0):
        self.res = r; con = r['con']; A = np.sort(con['areas'])[::-1]
        loose, gaps = regions.loose_ends(r['emb']); n_ends = 2*(len(r['emb']['start'])-1)
        self.lbl_regions.setText(f'{con["nreg"]} regions; largest {A[0]:.0f} mm², median {np.median(A):.2f} mm²; {len(r["clashes"])} colour clashes; {len(loose)} of {n_ends} line ends loose; {len(gaps)} stretches not drawn')
        self.logmsg(f'regions found in {time.time()-t0:.1f} s: {con["nreg"]} regions, {len(loose)} of {n_ends} line ends loose, {len(gaps)} stretches not drawn through'); self.busy(False); self.redraw()

    # ---------------------------------------------------------------- view
    def redraw(self):
        if self.P is None: return
        cam = self.plotter.camera_position; had = self.act is not None; self.plotter.renderer.clear_actors(); self.act = None; self.locator = None
        if self.res is None:
            mesh = pv.PolyData(self.P, np.hstack([np.full((len(self.F), 1), 3), self.F]).ravel()); self.act = self.plotter.add_mesh(mesh, color='lightgrey')
        else:
            r = self.res; DV, DF, pid = r['DV'], r['DF'], r['pid']
            self.mesh = pv.PolyData(DV, np.hstack([np.full((len(DF), 1), 3), DF]).ravel()); self.mesh.cell_data['c'] = r['colour'][r['con']['region'][pid]].astype(float)
            self.act = self.plotter.add_mesh(self.mesh, scalars='c', cmap=list(regions.PALETTE), clim=(-0.5, len(regions.PALETTE)-0.5), show_scalar_bar=False)
            if self.chk_lines.isChecked():
                emb = r['emb']; Q = emb['nodes']; start = emb['start']; nf = emb['nface']
                N = np.cross(emb['e1'], emb['e2']); Q = Q+0.01*N[nf]                   # a hair off the surface, else the lines hide in the colours
                pts = []; cells = []; n = 0
                for a, b in zip(start[:-1], start[1:]): pts.append(Q[a:b]); cells.append(np.r_[b-a, n+np.arange(b-a)]); n += b-a
                self.plotter.add_mesh(pv.PolyData(np.vstack(pts), lines=np.concatenate(cells)), color='black', line_width=1.5)
            if self.chk_loose.isChecked():
                loose, gaps = regions.loose_ends(r['emb'])
                if len(loose): self.plotter.add_mesh(pv.PolyData(loose), color='red', point_size=9, render_points_as_spheres=True)
                if len(gaps): self.plotter.add_mesh(pv.PolyData(gaps), color='magenta', point_size=9, render_points_as_spheres=True)
        if cam is not None and had: self.plotter.camera_position = cam
        else: self.plotter.reset_camera()
        self.plotter.render()
    def _pick_region(self):
        """the region under the cursor, or None"""
        if self.res is None or self.act is None: return None
        import vtk
        if self.locator is None:
            self.locator = vtk.vtkStaticCellLocator(); self.locator.SetDataSet(self.mesh); self.locator.BuildLocator()
            self.picker = vtk.vtkCellPicker(); self.picker.SetTolerance(0.0005); self.picker.PickFromListOn(); self.picker.AddPickList(self.act); self.picker.AddLocator(self.locator)
        x, y = self.plotter.iren.get_event_position()
        if not self.picker.Pick(x, y, 0, self.plotter.renderer) or self.picker.GetCellId() < 0: return None
        r = int(self.res['con']['region'][self.res['pid'][self.picker.GetCellId()]]); return r if r >= 0 else None
    def _describe(self, r):
        if r is None: return 'click a region'
        con = self.res['con']; deg = sum(1 for a, b in con['pairs'] if a == r or b == r)
        return f'region {r}: {con["areas"][r]:.2f} mm², {deg} neighbours, colour {int(self.res["colour"][r])+1}'
    def _on_press(self, *_): self._press = self.plotter.iren.get_event_position()
    def _on_release(self, *_):
        pos = self.plotter.iren.get_event_position()
        if self._press is None or abs(pos[0]-self._press[0]) > 3 or abs(pos[1]-self._press[1]) > 3: return   # a drag rotated the view
        r = self._pick_region()
        if r is None: return
        c = self.swatches.checkedId(); self.res['colour'][r] = c
        self.mesh.cell_data['c'] = self.res['colour'][self.res['con']['region'][self.res['pid']]].astype(float); self.plotter.render()
        self.lbl_hover.setText(self._describe(r)); self.logmsg(f'region {r} -> colour {c+1}')

    # ---------------------------------------------------------------- files
    def save(self):
        if self.res is None: return
        out = os.path.splitext(self.path)[0]+'_paint.pkl'
        path, _ = QtWidgets.QFileDialog.getSaveFileName(self, 'Save painting', out, 'Painting (*.pkl)')
        if not path: return
        t0 = time.time()
        with open(path, 'wb') as f: pickle.dump(dict(model=self.path, lines=self.lbl_lines.text(), set=self.cb_set.currentText(), res=self.res), f, protocol=pickle.HIGHEST_PROTOCOL)
        self.logmsg(f'saved {path} ({os.path.getsize(path)/1e6:.0f} MB) in {time.time()-t0:.1f} s')
    def load_painting(self):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, 'Open painting', os.path.dirname(self.path) if self.path else '', 'Painting (*.pkl)')
        if not path: return
        t0 = time.time()
        try:
            with open(path, 'rb') as f: d = pickle.load(f)
        except Exception as ex: QtWidgets.QMessageBox.critical(self, 'Load failed', str(ex)); return
        if 'res' not in d: d = dict(res=d, set='?', model='?')                      # a raw regions.run result (the smoke scripts write those)
        if len(d['res']['emb']['F']) != len(self.F): QtWidgets.QMessageBox.critical(self, 'Load failed', 'this painting belongs to a model with a different triangle count'); return
        self.res = d['res']; con = self.res['con']; self.lbl_regions.setText(f'{con["nreg"]} regions (from {os.path.basename(path)}, {d["set"]})')
        self.logmsg(f'loaded painting {path} ({d["set"]}, {con["nreg"]} regions) in {time.time()-t0:.1f} s'); self.busy(False); self.redraw()
    def export(self):
        if self.res is None: return
        out = os.path.splitext(self.path)[0]+'_painted.3mf'
        path, _ = QtWidgets.QFileDialog.getSaveFileName(self, 'Export to Bambu Studio', out, 'Bambu Studio project (*.3mf)')
        if not path: return
        fine = self.sp_fine.value(); depth = int(self.sp_edepth.value()); res = self.res; P, F = self.P, self.F; name = os.path.splitext(os.path.basename(self.path))[0]
        self.logmsg(f'--- export: finest split {fine} mm, depth {depth}'); t0 = time.time()
        def job(log):
            notes = bambu3mf.paint_notes(res, res['colour']+1, fine=fine, depth=depth, log=log)
            return bambu3mf.write_project(path, P, F, notes, regions.PALETTE, name=name, log=log)
        self._start(job, lambda p: (self.logmsg(f'exported {p} in {time.time()-t0:.1f} s'), self.busy(False)))


if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv); w = App(); w.show(); sys.exit(app.exec())
