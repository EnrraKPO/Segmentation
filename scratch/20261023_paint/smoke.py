"""regions.py on the orc: walk lines of dump 20261023_114229 (or a Save lines file given as argument), counts, timings, renders."""
import sys, os, time, pickle; root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); sys.path.insert(0, root)
import numpy as np
import regions, linesio

if __name__ == '__main__':
    path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(root, 'out', 'debug', '20261023_114229', 'lines.npz')
    which = sys.argv[2] if len(sys.argv) > 2 else 'walk'
    t0 = time.time(); P, F = linesio.load_mesh(os.path.join(root, 'out', 'orcfull_patch.stl')); st = linesio.read_lines(path)
    lines = st['walk'] if which == 'walk' else st[which]['lines']
    print(f'mesh {len(F)} faces, {len(lines)} {which} lines ({st["note"]}), loaded in {time.time()-t0:.1f} s')
    res = regions.run(P, F, lines, log=print, min_area=float(sys.argv[3]) if len(sys.argv) > 3 else 0.16)
    con = res['con']; A = np.sort(con['areas'])[::-1]
    print(f'regions {con["nreg"]}: largest areas {np.round(A[:8], 1)} mm2; regions under 1 mm2: {(A < 1).sum()}, under 0.01 mm2: {(A < 0.01).sum()}; clashes {len(res["clashes"])}')
    with open(os.path.join(root, 'out', f'paint_{which}.pkl'), 'wb') as f: pickle.dump(res, f)
    import pyvista as pv
    DV, DF, pid = res['DV'], res['DF'], res['pid']; col = res['colour'][con['region'][pid]]
    mesh = pv.PolyData(DV, np.hstack([np.full((len(DF), 1), 3), DF]).ravel()); mesh.cell_data['c'] = col
    cmap = [regions.PALETTE[i] for i in range(9)]
    for name, cam in [('front', None), ('rivet', [(173.18+6, 143.4+6, 98.66+6), (173.18, 143.4, 98.66), (0, 0, 1)])]:
        pl = pv.Plotter(off_screen=True, window_size=(1600, 1200)); pl.add_mesh(mesh, scalars='c', cmap=cmap, clim=(0, 8), show_scalar_bar=False)
        Q = res['emb']['nodes']; start = res['emb']['start']
        pts = []; cells = []; n = 0
        for a, b in zip(start[:-1], start[1:]):
            pts.append(Q[a:b]); cells.append(np.r_[b-a, n+np.arange(b-a)]); n += b-a
        pl.add_mesh(pv.PolyData(np.vstack(pts), lines=np.concatenate(cells)), color='black', line_width=1)
        if cam is not None: pl.camera_position = cam; pl.camera.zoom(1.0)
        out = os.path.join(root, 'out', f'paint_{which}_{name}.png'); pl.screenshot(out); print('wrote', out)
