"""rerun refine on a few walk lines of the dump with the log on (end reasons of the fine tracer), render each: walk line blue,
other walk lines cyan, rerun pieces red, tube dots by bend. usage: python rerun.py dumpdir line [line...]"""
import sys, os, ast, time, numpy as np; from scipy.spatial import cKDTree
sys.path.insert(0, r'D:\Stuff\PrintSuit\SegmentationTool'); os.chdir(r'D:\Stuff\PrintSuit\SegmentationTool'); import refine, trimesh
def main():
  d = sys.argv[1]; ids = [int(a) for a in sys.argv[2:]]
  w = np.load(d+'/lines.npz'); W = [w['pts'][a:a+k] for a, k in zip(np.cumsum(np.r_[0, w['n'][:-1]]), w['n'])]
  line12 = [l for l in open(d+'/log.txt', encoding='utf-8') if l.startswith('--- refine all')][0]
  walk = ast.literal_eval(line12.split('fine walk ')[1].split('; tune')[0]); tune = ast.literal_eval(line12.split('tune ')[-1].strip())
  head = dict(x.split('=') for x in line12.split(': ')[1].split(';')[0].split(', ')); finer, radius, area, cut, unit = int(head['finer']), float(head['radius']), float(head['area']), float(head['cut']), float(head['unit'])
  print('fine walk method: trail', walk.get('trail'), 'curve', walk.get('curve'), '| finer', finer, 'radius', radius, 'area', area, 'cut', cut, flush=True)
  m = trimesh.load('out/orcfull_patch.stl', process=True); P, F = np.asarray(m.vertices), np.asarray(m.faces); h = float((P.max(0)-P.min(0)).max()); cube = 0.087
  tree, owner = refine.nodes_tree(W)
  import pyvista as pv
  mesh = pv.PolyData(P, np.hstack([np.full((len(F), 1), 3), F]).ravel()); bb = np.array([(L.min(0), L.max(0)) for L in W])
  def snap(p): C, _, f = trimesh.proximity.closest_point(m, p); return C+0.01*m.face_normals[f]
  for i in ids:
      L = W[i]; others = refine.end_lines(W, i, tune['slack']*cube, tree, owner); print(f'\n=== line {i}: {len(L)} walk nodes, ends sit on lines {others}', flush=True); t0 = time.time()
      lines, t = refine.refine_line(L, P, F, cube, h, finer, radius, cut, True, print, 1, walk, area, unit, [W[j] for j in others], tune)
      print(f'-> {len(lines)} fine pieces, nodes {[len(l) for l in lines]}, {time.time()-t0:.0f} s', flush=True)
      E = np.array([q for l in lines for q in (l[0], l[-1])]); wt = cKDTree(np.vstack([W[j] for j in range(len(W)) if j != i]))
      print('   piece ends: distance to the nearest other walk line', np.round(wt.query(E)[0], 2).tolist(), flush=True)
      C = t['lo']+(t['idx'][t['line']]+0.5)*t['size']; T = t['turn'][t['line']]
      c = 0.5*(L.min(0)+L.max(0)); half = 0.6*max((L.max(0)-L.min(0)).max(), 1.0)+0.3
      pl = pv.Plotter(off_screen=True, window_size=(1600, 1200)); pl.set_background('white'); pl.add_mesh(mesh, color='lightgrey', smooth_shading=True)
      pd = pv.PolyData(C); pd['bend'] = T; pl.add_mesh(pd, scalars='bend', cmap='plasma', point_size=3, render_points_as_spheres=True, show_scalar_bar=False)
      near = np.flatnonzero((bb[:, 1] >= c-half).all(1) & (bb[:, 0] <= c+half).all(1))
      for j in near:
          p = snap(W[j]); pl.add_mesh(pv.PolyData(p, lines=np.r_[len(p), np.arange(len(p))]), color='blue' if j == i else 'cyan', line_width=4 if j == i else 2)
      for l in lines:
          p = snap(l); pl.add_mesh(pv.PolyData(p, lines=np.r_[len(p), np.arange(len(p))]), color='red', line_width=5)
          for e in (p[0], p[-1]): pl.add_mesh(pv.Sphere(0.03, e), color='black')
      n = m.face_normals[trimesh.proximity.closest_point(m, c[None])[2][0]]; pl.camera.focal_point = c; pl.camera.position = c+n*20; pl.camera.up = (0, 0, 1); pl.camera.parallel_projection = True; pl.camera.parallel_scale = half
      pl.add_text(f'line {i}: walk blue, rerun fine pieces red (ends black), other walk lines cyan, tube dots by bend', font_size=10); pl.screenshot(f'out/joints_line{i}.png'); pl.close()
      print(f'   out/joints_line{i}.png', flush=True)


if __name__ == '__main__': main()
