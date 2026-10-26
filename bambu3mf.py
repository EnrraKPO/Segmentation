"""Bambu Studio project (3MF) with the painting as per-triangle paint notes. The mesh is written exactly as it came.

Bambu Studio (like PrusaSlicer) keeps colour painting as a note per triangle, attribute `paint_color` of <triangle>, read
from its source (libslic3r/TriangleSelector.cpp serialize/perform_split, Model.cpp get_triangle_as_string):
  - a bit stream, 4 bits a nibble; the hex string lists the nibbles LAST FIRST (get_triangle_as_string inserts at the front)
  - nibble, bits low to high: b0 b1 = number of split sides; if 0 (leaf): b2 b3 = state 0..2, or 0b11 then the next nibble
    holds state-3 (states up to 17; a run of 1111 nibbles adds 15 each); if split: b2 b3 = special side, then the children
    in REVERSE order (child index split_sides .. 0)
  - states: 0 none (default filament), 1 = filament 1, 2 = filament 2, 3.. = filament 3..
  - splitting with all three sides (our choice; special side 0): midpoints A of (V0,V1), B of (V1,V2), C of (V2,V0);
    children 0:(V0,A,C) 1:(A,V1,B) 2:(B,V2,C) 3:(A,B,C). Bambu rebuilds the same split in memory when it loads.
Our note for a crossed triangle: split while a child still has strokes inside it, until its longest edge is below `fine`
(mm) or the depth reaches `depth`; a leaf takes the colour of the region under its centre. Uncrossed triangles: a leaf.
"""
import os, io, time, zipfile, uuid, datetime
import numpy as np


def leaf_bits(state):
    bits = [0, 0]
    if state < 3: bits += [state & 1, (state >> 1) & 1]
    else:
        bits += [1, 1]; n = state-3
        while n >= 15: bits += [1, 1, 1, 1]; n -= 15
        bits += [(n >> i) & 1 for i in range(4)]
    return bits


def to_string(bits):
    """nibbles (4 bits each, low bit first) to Bambu's hex string, last nibble first"""
    assert len(bits) % 4 == 0
    out = []
    for i in range(0, len(bits), 4):
        v = bits[i] | bits[i+1] << 1 | bits[i+2] << 2 | bits[i+3] << 3; out.append('0123456789ABCDEF'[v])
    return ''.join(reversed(out))


def from_string(s):
    """the reverse, for checking: bits low first"""
    bits = []
    for ch in reversed(s):
        v = int(ch, 16); bits += [(v >> i) & 1 for i in range(4)]
    return bits


def _note_bits(tri, strokes, piece_of, fine, depth, prepared):
    """tri: (3,2) corners in the triangle frame; strokes: shapely geometry of the noded strokes; piece_of(x, y) -> state"""
    import shapely
    from shapely.geometry import Polygon
    V0, V1, V2 = tri; poly = Polygon([V0, V1, V2])
    inside = prepared.intersects(poly.buffer(-1e-9)) if poly.area > 0 else False     # strokes through the interior, not just along the boundary
    longest = max(np.linalg.norm(V1-V0), np.linalg.norm(V2-V1), np.linalg.norm(V0-V2))
    if not inside or longest < fine or depth <= 0:
        c = (V0+V1+V2)/3; return leaf_bits(piece_of(c[0], c[1]))
    A = (V0+V1)/2; B = (V1+V2)/2; C = (V2+V0)/2
    children = [(V0, A, C), (A, V1, B), (B, V2, C), (A, B, C)]
    bits = [1, 1, 0, 0]                                                             # 3 split sides, special side 0
    for ch in reversed(children): bits += _note_bits(np.array(ch), strokes, piece_of, fine, depth-1, prepared)
    return bits


def _face_note(task):
    """task = (f, C2 (3,2), piece polygons wkb, piece states, fine, depth). Returns (f, note string, leaves)"""
    import shapely
    f, C, wkbs, states, fine, depth = task
    polys = [shapely.from_wkb(w) for w in wkbs]
    strokes = shapely.unary_union([p.boundary for p in polys]).difference(shapely.Polygon([C[0], C[1], C[2]]).boundary.buffer(1e-9))
    from shapely.prepared import prep
    prepared = prep(strokes) if not strokes.is_empty else None
    if prepared is None: return f, to_string(leaf_bits(states[0] if len(states) else 0)), 1
    tree = shapely.STRtree(polys)
    def piece_of(x, y):
        pt = shapely.Point(x, y); hits = tree.query(pt, predicate='intersects')
        if len(hits) == 0: hits = [int(np.argmin([p.distance(pt) for p in polys]))]
        return states[hits[0]]
    bits = _note_bits(np.asarray(C, float), strokes, piece_of, fine, depth, prepared)
    return f, to_string(bits), bits.count(0)  # leaves not needed precisely; count is informative only


def paint_notes(res, state_of_region, fine=0.05, depth=8, log=None, workers=None):
    """a note per triangle. res: regions.run result; state_of_region: int array, Bambu state (filament index 1..) per region.
    Returns list of strings (empty = no note) indexed by triangle."""
    from concurrent.futures import ProcessPoolExecutor
    t0 = time.time(); emb = res['emb']; part = res['part']; con = res['con']; F = emb['F']; m = len(F); region = con['region']
    notes = [''] * m
    plain = np.where(~con['crossed'])[0]
    for f in plain: notes[f] = to_string(leaf_bits(int(state_of_region[region[f]])))
    tasks = []
    for fi, f in enumerate(part['faces']):
        pcs = part['pieces'][fi]; p0 = con['pid0'][fi]
        tasks.append((int(f), emb['C2'][f], [pc[2] for pc in pcs], [int(state_of_region[region[p0+i]]) for i in range(len(pcs))], fine, depth))
    workers = workers or os.cpu_count() or 1
    if workers > 1 and len(tasks) > 256:
        with ProcessPoolExecutor(workers) as ex: out = list(ex.map(_face_note, tasks, chunksize=max(1, len(tasks)//(workers*8))))
    else: out = [_face_note(t) for t in tasks]
    for f, s, _ in out: notes[f] = s
    if log: log(f'  paint notes for {m} triangles ({len(tasks)} split, fine {fine} mm, depth {depth}) in {time.time()-t0:.1f} s; longest note {max(len(s) for s in notes)} nibbles')
    return notes


# ------------------------------------------------------------------ the 3MF package
def _xml_escape(s): return s.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;').replace('"', '&quot;')


def write_project(path, P, F, notes, palette, name='painted', log=None):
    """Bambu Studio project: the mesh in 3D/Objects/object_1.model with a paint_color per triangle, the usual package
    files around it, and a project_settings.config naming one filament per palette colour."""
    t0 = time.time(); P = np.asarray(P, float); F = np.asarray(F)
    uu = lambda: str(uuid.uuid4())
    # the object file
    buf = io.StringIO(); w = buf.write
    w('<?xml version="1.0" encoding="UTF-8"?>\n<model unit="millimeter" xml:lang="en-US" xmlns="http://schemas.microsoft.com/3dmanufacturing/core/2015/02" xmlns:BambuStudio="http://schemas.bambulab.com/package/2021" xmlns:p="http://schemas.microsoft.com/3dmanufacturing/production/2015/06" requiredextensions="p">\n')
    w(' <metadata name="BambuStudio:3mfVersion">1</metadata>\n <resources>\n')
    w(f'  <object id="1" p:UUID="{uu()}" type="model">\n   <mesh>\n    <vertices>\n')
    w(''.join(f'     <vertex x="{x:.9g}" y="{y:.9g}" z="{z:.9g}"/>\n' for x, y, z in P))
    w('    </vertices>\n    <triangles>\n')
    w(''.join(f'     <triangle v1="{a}" v2="{b}" v3="{c}"' + (f' paint_color="{s}"' if s else '') + '/>\n' for (a, b, c), s in zip(F.tolist(), notes)))
    w('    </triangles>\n   </mesh>\n  </object>\n </resources>\n <build/>\n</model>\n')
    obj_xml = buf.getvalue()
    today = datetime.date.today().isoformat(); nm = _xml_escape(name)
    main_xml = ('<?xml version="1.0" encoding="UTF-8"?>\n<model unit="millimeter" xml:lang="en-US" xmlns="http://schemas.microsoft.com/3dmanufacturing/core/2015/02" xmlns:BambuStudio="http://schemas.bambulab.com/package/2021" xmlns:p="http://schemas.microsoft.com/3dmanufacturing/production/2015/06" requiredextensions="p">\n'
                f' <metadata name="Application">BambuStudio-02.00.00.00</metadata>\n <metadata name="BambuStudio:3mfVersion">1</metadata>\n <metadata name="CreationDate">{today}</metadata>\n <metadata name="ModificationDate">{today}</metadata>\n <metadata name="Title">{nm}</metadata>\n'
                ' <resources>\n'
                f'  <object id="2" p:UUID="{uu()}" type="model">\n   <components>\n    <component p:path="/3D/Objects/object_1.model" objectid="1" p:UUID="{uu()}" transform="1 0 0 0 1 0 0 0 1 0 0 0"/>\n   </components>\n  </object>\n'
                ' </resources>\n'
                f' <build p:UUID="{uu()}">\n  <item objectid="2" p:UUID="{uu()}" transform="1 0 0 0 1 0 0 0 1 0 0 0" printable="1"/>\n </build>\n</model>\n')
    content_types = ('<?xml version="1.0" encoding="UTF-8"?>\n<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">\n'
                     ' <Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>\n <Default Extension="model" ContentType="application/vnd.ms-package.3dmanufacturing-3dmodel+xml"/>\n <Default Extension="png" ContentType="image/png"/>\n</Types>')
    rels = ('<?xml version="1.0" encoding="UTF-8"?>\n<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">\n'
            ' <Relationship Target="/3D/3dmodel.model" Id="rel-1" Type="http://schemas.microsoft.com/3dmanufacturing/2013/01/3dmodel"/>\n</Relationships>')
    model_rels = ('<?xml version="1.0" encoding="UTF-8"?>\n<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">\n'
                  ' <Relationship Target="/3D/Objects/object_1.model" Id="rel-1" Type="http://schemas.microsoft.com/3dmanufacturing/2013/01/3dmodel"/>\n</Relationships>')
    model_settings = ('<?xml version="1.0" encoding="UTF-8"?>\n<config>\n  <object id="2">\n'
                      f'    <metadata key="name" value="{nm}"/>\n    <metadata key="extruder" value="1"/>\n    <metadata face_count="{len(F)}"/>\n'
                      f'    <part id="1" subtype="normal_part">\n      <metadata key="name" value="{nm}"/>\n      <metadata key="matrix" value="1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1"/>\n'
                      f'      <mesh_stat face_count="{len(F)}" edges_fixed="0" degenerate_facets="0" facets_removed="0" facets_reversed="0" backwards_edges="0"/>\n    </part>\n  </object>\n'
                      '  <plate>\n    <metadata key="plater_id" value="1"/>\n    <metadata key="plater_name" value=""/>\n    <metadata key="locked" value="false"/>\n'
                      '    <model_instance>\n      <metadata key="object_id" value="2"/>\n      <metadata key="instance_id" value="0"/>\n      <metadata key="identify_id" value="1"/>\n    </model_instance>\n  </plate>\n</config>\n')
    import json
    colours = [c if len(c) == 7 else c[:7] for c in palette]
    project_settings = json.dumps({'filament_colour': colours, 'filament_type': ['PLA']*len(colours), 'filament_settings_id': ['Generic PLA']*len(colours)}, indent=4)
    with zipfile.ZipFile(path, 'w', zipfile.ZIP_DEFLATED) as z:
        z.writestr('[Content_Types].xml', content_types); z.writestr('_rels/.rels', rels)
        z.writestr('3D/3dmodel.model', main_xml); z.writestr('3D/_rels/3dmodel.model.rels', model_rels); z.writestr('3D/Objects/object_1.model', obj_xml)
        z.writestr('Metadata/model_settings.config', model_settings); z.writestr('Metadata/project_settings.config', project_settings)
    if log: log(f'  wrote {path} ({os.path.getsize(path)/1e6:.1f} MB, {len(F)} triangles) in {time.time()-t0:.1f} s')
    return path
