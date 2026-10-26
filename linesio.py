"""Reading the line app's files: Save lines files, debugging dumps, meshes. Shared by app.py and paint.py."""
import os
import numpy as np, trimesh


def _split(n, pts, dtype=float):
    off = np.r_[0, np.cumsum(n)]; return [np.asarray(pts[off[i]:off[i+1]], dtype) for i in range(len(n))]


def _npz(path):
    """an .npz read whole and closed (np.load keeps the file open otherwise)"""
    with np.load(path, allow_pickle=True) as z: return {k: z[k] for k in z.files}


def read_lines(path):
    """Lines of an earlier session, from a Save lines file (the stages are stored in it since 2026-10-23) or from a debugging
    dump folder (pick any of its lines.npz / refined.npz / joined.npz). An older save file holds one flat set: it comes back
    as walk lines. Returns dict(walk, walk_shared, refined or None, joined or None, note)."""
    z = _npz(path); d, base = os.path.split(path); out = dict(walk=[], walk_shared=None, refined=None, joined=None, note='')
    def joined_from(j):
        n = j['n']; return dict(lines=_split(n, j['pts']), measured=_split(n, j['measured'], bool), walk=list(j['walk']), shared=[(q,) for q in np.asarray(j['shared']).reshape(-1, 3)],
                                loose=[(q,) for q in np.asarray(j['loose']).reshape(-1, 3)], records=[], report={}, agreement=None)
    if 'walk_n' in z:                                                                 # a save file with the stages
        out['walk'] = _split(z['walk_n'], z['walk_pts']); ws = z['walk_shared']; out['walk_shared'] = ws if len(ws) else None
        if 'ref_n' in z: out['refined'] = dict(lines=_split(z['ref_n'], z['ref_pts']), index=[int(i) for i in z['ref_index']], tubes=[])
        if 'joined_n' in z: out['joined'] = joined_from({k[7:]: v for k, v in z.items() if k.startswith('joined_')})
        out['note'] = 'save file'
    elif base in ('lines.npz', 'refined.npz', 'joined.npz') and os.path.exists(os.path.join(d, 'lines.npz')):   # a dump folder
        w = _npz(os.path.join(d, 'lines.npz')); out['walk'] = _split(w['n'], w['pts'])
        if os.path.exists(os.path.join(d, 'refined.npz')):
            r = _npz(os.path.join(d, 'refined.npz')); out['refined'] = dict(lines=_split(r['n'], r['pts']), index=[int(i) for i in r['index']], tubes=[])
        if os.path.exists(os.path.join(d, 'joined.npz')): out['joined'] = joined_from(_npz(os.path.join(d, 'joined.npz')))
        out['note'] = 'dump folder'
    else:                                                                             # an older save: one flat set
        out['walk'] = _split(z['n'], z['pts']); prm = str(z['params']) if 'params' in z else ''
        what = 'joined lines' if "'joined':" in prm else 'walk lines with refined pieces in place' if "'refined':" in prm else 'walk lines'
        out['note'] = f'older save file, one flat set: {what}'; out['flat'] = what
    return out


def load_mesh(path):
    m = trimesh.load(path, process=True)
    if isinstance(m, trimesh.Scene): m = trimesh.util.concatenate([g for g in m.geometry.values()])
    return np.asarray(m.vertices, np.float64), np.asarray(m.faces, np.int64)
