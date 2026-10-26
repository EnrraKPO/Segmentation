"""Decode every paint_color string of a Bambu-saved object file into its split tree and encode it back: the strings must
come out identical (checks nibble layout, child order and string direction of bambu3mf.py)."""
import sys, os, re, collections; root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); sys.path.insert(0, root)
import bambu3mf

def decode(bits):
    """bits -> nested tree: int state for a leaf, (special, [children in index order]) for a split"""
    pos = [0]
    def nib():
        v = bits[pos[0]] | bits[pos[0]+1] << 1 | bits[pos[0]+2] << 2 | bits[pos[0]+3] << 3; pos[0] += 4; return v
    def node():
        code = nib(); split = code & 3
        if split == 0:
            if code >> 2 == 3:
                nxt = nib(); num = 0
                while nxt == 15: num += 1; nxt = nib()
                return nxt+15*num+3
            return code >> 2
        children = [node() for _ in range(split+1)][::-1]           # serialized in reverse
        return (code >> 2, split, children)
    t = node(); assert pos[0] == len(bits), (pos[0], len(bits)); return t

def encode(t):
    if isinstance(t, int): return bambu3mf.leaf_bits(t)
    special, split, children = t; bits = [split & 1, (split >> 1) & 1, special & 1, (special >> 1) & 1]
    for ch in reversed(children): bits += encode(ch)
    return bits

path = r"C:\Users\artur\AppData\Local\Temp\claude\D--Stuff-PrintSuit-SegmentationTool\67e9db24-c45b-46f8-9786-7ee80add74d8\scratchpad\bambu_3mf\bts head test 2 focus w highlight + lashes\3D\Objects\object_18.model"
strings = re.findall(r'paint_color="([0-9A-F]+)"', open(path, encoding='utf-8').read())
cnt = collections.Counter(strings); bad = 0; states = collections.Counter(); splits = collections.Counter()
for s in cnt:
    t = decode(bambu3mf.from_string(s)); back = bambu3mf.to_string(encode(t))
    if back != s: bad += 1; print('MISMATCH', s, back)
    def walk(t, d=0):
        if isinstance(t, int): states[t] += 1
        else: splits[t[1]] += 1; [walk(c, d+1) for c in t[2]]
    walk(t)
print(f'{len(strings)} strings, {len(cnt)} distinct, {bad} mismatches; leaf states {dict(states)}; split kinds {dict(splits)}')
print('our leaves:', {s: bambu3mf.to_string(bambu3mf.leaf_bits(s)) for s in range(1, 10)})
