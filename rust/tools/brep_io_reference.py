"""Independent reader of OCCT's .brep text (T2 of TOPOLOGY_MODEL.md), for the
data/occ corpus, without Rust or OCCT.

It follows dox/specification/brep_format.md and the pinned
BRepTools_ShapeSet / TopTools_LocationSet readers: matrix location records are
always numbered and an empty composite chain is not; UV points only in
version 2; B-spline records carry rational and periodic flags; a closed
surface's continuity may be glued to the second pcurve number (`6CN`).

For every file it reports the geometry records the kernel cannot represent
(by the kernel's names: every curve, 2D curve and surface other than a line,
a circle, a plane or a cylinder; a trimmed line or circle counts as its
basis), and for every solid reached from the root through compounds:

* whether its structure is representable: plane and cylinder faces (direct
  or indirect), edges that are not degenerate, have a forward and a reversed
  vertex and a line or circle 3D curve, a pcurve on every non-plane face they
  bound, forward or reversed orientations only, and rigid locations;
* OCCT's own distinct subshape counts of the original, seamed solid
  (vertices, edges, wires, faces, shells, solids by record and placement):
  the counts native `nbshapes` reports, and the kernel's synthesized counts
  must reproduce them for every solid it imports.
"""
from pathlib import Path
import re


class Reader:
    def __init__(self, text):
        self.words = text.split()
        self.at = 0

    def word(self):
        w = self.words[self.at]
        self.at += 1
        return w

    def real(self):
        return float(self.word())

    def int(self):
        return int(self.word())

    def reals(self, n):
        return [self.real() for _ in range(n)]


def curve3(r):
    kind = r.int()
    if kind == 1:
        r.reals(6)
        return 'line'
    if kind == 2:
        r.reals(13)
        return 'circle'
    if kind in (3, 5):
        r.reals(14)
        return {3: 'Ellipse', 5: 'Hyperbola'}[kind]
    if kind == 4:
        r.reals(13)
        return 'Parabola'
    if kind == 6:
        rational, degree = r.int(), r.int()
        r.reals((degree+1)*(4 if rational else 3))
        return 'BezierCurve'
    if kind == 7:
        rational, _periodic, _degree, poles, knots = (r.int() for _ in range(5))
        r.reals(poles*(4 if rational else 3)+2*knots)
        return 'BSplineCurve'
    if kind == 8:
        r.reals(2)
        basis = curve3(r)
        return basis if basis in ('line', 'circle') else 'TrimmedCurve'
    if kind == 9:
        r.reals(4)
        curve3(r)
        return 'OffsetCurve'
    raise ValueError(f'3D curve type {kind}')


def curve2(r):
    kind = r.int()
    sizes = {1: 4, 2: 7, 3: 8, 4: 7, 5: 8}
    names = {1: 'line', 2: 'circle', 3: 'Ellipse2d', 4: 'Parabola2d', 5: 'Hyperbola2d'}
    if kind in sizes:
        r.reals(sizes[kind])
        return names[kind]
    if kind == 6:
        rational, degree = r.int(), r.int()
        r.reals((degree+1)*(3 if rational else 2))
        return 'BezierCurve2d'
    if kind == 7:
        rational, _periodic, _degree, poles, knots = (r.int() for _ in range(5))
        r.reals(poles*(3 if rational else 2)+2*knots)
        return 'BSplineCurve2d'
    if kind == 8:
        r.reals(2)
        basis = curve2(r)
        return basis if basis in ('line', 'circle') else 'TrimmedCurve2d'
    if kind == 9:
        r.reals(1)
        curve2(r)
        return 'OffsetCurve2d'
    raise ValueError(f'2D curve type {kind}')


def surface(r):
    kind = r.int()
    if kind == 1:
        r.reals(12)
        return 'plane'
    if kind == 2:
        r.reals(13)
        return 'cylinder'
    names = {3: 'ConicalSurface', 4: 'SphericalSurface', 5: 'ToroidalSurface'}
    if kind in (3, 5):
        r.reals(14)
        return names[kind]
    if kind == 4:
        r.reals(13)
        return names[kind]
    if kind == 6:
        r.reals(3)
        curve3(r)
        return 'SurfaceOfLinearExtrusion'
    if kind == 7:
        r.reals(6)
        curve3(r)
        return 'SurfaceOfRevolution'
    if kind == 8:
        ru, rv, du, dv = (r.int() for _ in range(4))
        r.reals((du+1)*(dv+1)*(4 if ru or rv else 3))
        return 'BezierSurface'
    if kind == 9:
        ru, rv, _pu, _pv, _du, _dv, nu, nv, ku, kv = (r.int() for _ in range(10))
        r.reals(nu*nv*(4 if ru or rv else 3)+2*(ku+kv))
        return 'BSplineSurface'
    if kind == 10:
        r.reals(4)
        surface(r)
        return 'RectangularTrimmedSurface'
    if kind == 11:
        r.reals(1)
        surface(r)
        return 'OffsetSurface'
    raise ValueError(f'surface type {kind}')


def matmul(a, b):
    """3x4 affine a * b (b applied first)."""
    out = []
    for i in range(3):
        row = []
        for j in range(4):
            v = sum(a[4*i+k]*b[4*k+j] for k in range(3))+(a[4*i+3] if j == 3 else 0.0)
            row.append(v)
        out += row
    return out


IDENTITY = [1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0]


def inverse(m):
    a = [[m[4*i+j] for j in range(3)] for i in range(3)]
    det = (a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])-a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])
           + a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]))
    inv = [[(a[(j+1) % 3][(i+1) % 3]*a[(j+2) % 3][(i+2) % 3]-a[(j+1) % 3][(i+2) % 3]*a[(j+2) % 3][(i+1) % 3])/det
            for j in range(3)] for i in range(3)]
    t = [-(inv[i][0]*m[3]+inv[i][1]*m[7]+inv[i][2]*m[11]) for i in range(3)]
    return [inv[0][0], inv[0][1], inv[0][2], t[0], inv[1][0], inv[1][1], inv[1][2], t[1],
            inv[2][0], inv[2][1], inv[2][2], t[2]]


def rigid(m, tol=1e-12):
    cols = [[m[4*i+j] for i in range(3)] for j in range(3)]
    dot = lambda p, q: sum(x*y for x, y in zip(p, q))
    x, y, z = cols
    cross = (x[1]*y[2]-x[2]*y[1], x[2]*y[0]-x[0]*y[2], x[0]*y[1]-x[1]*y[0])
    return (abs(dot(x, x)-1) <= tol and abs(dot(y, y)-1) <= tol and abs(dot(z, z)-1) <= tol
            and abs(dot(x, y)) <= tol and abs(dot(y, z)) <= tol and abs(dot(x, z)) <= tol
            and abs(dot(cross, z)-1) <= tol)


def near(a, b):
    return all(abs(x-y) <= 1e-12*(1+max(abs(x), abs(y))) for x, y in zip(a, b))


def read(text):
    r = Reader(text)
    while r.word() != 'Topology':
        pass
    version = int(r.word()[1])
    r.word()
    r.word()
    assert r.word() == 'Locations'
    locations = []
    for _ in range(r.int()):
        if r.int() == 1:
            locations.append(r.reals(12))
        else:
            l, empty = IDENTITY, True
            while True:
                index = r.int()
                if index == 0:
                    break
                empty = False
                power = r.int()
                base = locations[index-1] if power >= 0 else inverse(locations[index-1])
                for _ in range(abs(power)):
                    l = matmul(base, l)
            if not empty:
                locations.append(l)
    tables = {}
    for name, parse in (('Curve2ds', curve2), ('Curves', curve3)):
        assert r.word() == name
        tables[name] = [parse(r) for _ in range(r.int())]
    assert r.word() == 'Polygon3D'
    for _ in range(r.int()):
        nodes, params = r.int(), r.int()
        r.reals(1+3*nodes+(nodes if params else 0))
    assert r.word() == 'PolygonOnTriangulations'
    for _ in range(r.int()):
        nodes = r.int()
        r.reals(nodes)
        assert r.word() == 'p'
        r.real()
        if r.int():
            r.reals(nodes)
    assert r.word() == 'Surfaces'
    tables['Surfaces'] = [surface(r) for _ in range(r.int())]
    assert r.word() == 'Triangulations'
    for _ in range(r.int()):
        nodes, triangles, uv = r.int(), r.int(), r.int()
        normals = r.int() if version >= 3 else 0
        r.real()
        r.reals(3*nodes+(2*nodes if uv else 0)+3*triangles+(3*nodes if normals else 0))
    assert r.word() == 'TShapes'
    n = r.int()
    shapes = []
    for _ in range(n):
        kind = r.word()
        data = None
        if kind == 'Ve':
            r.reals(4)
            while True:
                r.real()
                t = r.int()
                if t == 0:
                    break
                r.reals({1: 1, 2: 2, 3: 2}[t]+1)
        elif kind == 'Ed':
            r.real()
            r.int()
            r.int()
            degenerated = r.int() == 1
            reps = []
            while True:
                t = r.int()
                if t == 0:
                    break
                if t == 1:
                    c, loc = r.int(), r.int()
                    r.reals(2)
                    reps.append(('curve', c, loc))
                elif t in (2, 3):
                    pcs = [r.int()]
                    if t == 3:
                        w = r.word()
                        digits = len(w)-len(w.lstrip('0123456789'))
                        pcs.append(int(w[:digits]))
                        if digits == len(w):
                            r.word()
                    s, loc = r.int(), r.int()
                    r.reals(2+(4 if version == 2 else 0))
                    reps.append(('pcurve', pcs, s, loc))
                elif t == 4:
                    r.word()
                    r.reals(4)
                elif t == 5:
                    r.reals(2)
                else:
                    r.reals(4 if t == 7 else 3)
            data = (degenerated, reps)
        elif kind == 'Fa':
            flag = r.int()
            if flag == 2:
                r.int()
                data = None
            else:
                r.real()
                data = (r.int(), r.int())
                if r.words[r.at] == '2':
                    r.at += 2
        flags = r.word()
        assert len(flags) == 7
        subs = []
        while True:
            w = r.word()
            if w == '*':
                break
            subs.append((w[0], n-int(w[1:]), r.int()))
        shapes.append((kind, data, subs))
    w = r.word()
    # operator>> takes the location's leading digits and nothing follows.
    location = r.word()
    root = (w[0], n-int(w[1:]), int(re.match(r'\d+', location).group()))
    return locations, tables, shapes, root


def summary(text):
    """(unsupported geometry {name: count}, [(record, representable, counts)])."""
    locations, tables, shapes, root = read(text)
    loc = lambda i: IDENTITY if i == 0 else locations[i-1]
    unsupported = {}
    for name in ('Curves', 'Curve2ds', 'Surfaces'):
        for kind in tables[name]:
            if kind not in ('line', 'circle', 'plane', 'cylinder'):
                unsupported[kind] = unsupported.get(kind, 0)+1
    solids = []
    stack = [(root, IDENTITY, '+')]
    while stack:
        (o, index, l), parent, orient = stack.pop()
        t = matmul(parent, loc(l))
        kind, _, subs = shapes[index]
        if o not in '+-' or not rigid(t):
            continue
        if kind == 'Co':
            for s in reversed(subs):
                stack.append((s, t, o))
        elif kind == 'So':
            solids.append((index,)+solid(shapes, tables, loc, index, t))
    return unsupported, solids


def solid(shapes, tables, loc, record, t):
    ok = True
    seen = {k: [] for k in ('Ve', 'Ed', 'Wi', 'Fa', 'Sh', 'So')}

    def mark(kind, index, m):
        if any(i == index and near(u, m) for i, u in seen[kind]):
            return False
        seen[kind].append((index, m))
        return True
    mark('So', record, t)
    shells = shapes[record][2]
    if not shells:
        ok = False
    for o, s, l in shells:
        st = matmul(t, loc(l))
        if o not in '+-' or shapes[s][0] != 'Sh':
            ok = False
            continue
        mark('Sh', s, st)
        for fo, f, fl in shapes[s][2]:
            ft = matmul(st, loc(fl))
            kind, data, wires = shapes[f]
            if fo not in '+-' or kind != 'Fa' or data is None:
                ok = False
                continue
            mark('Fa', f, ft)
            surf = tables['Surfaces'][data[0]-1] if data[0] else None
            if surf not in ('plane', 'cylinder'):
                ok = False
            for wo, w, wl in wires:
                wt = matmul(ft, loc(wl))
                if wo not in '+-' or shapes[w][0] != 'Wi':
                    ok = False
                    continue
                mark('Wi', w, wt)
                for eo, e, el in shapes[w][2]:
                    et = matmul(wt, loc(el))
                    if eo not in '+-':
                        ok = False
                    mark('Ed', e, et)
                    degenerated, reps = shapes[e][1]
                    curves = [x[1] for x in reps if x[0] == 'curve']
                    if degenerated or not curves or tables['Curves'][curves[0]-1] not in ('line', 'circle'):
                        ok = False
                    on = [x for x in reps if x[0] == 'pcurve' and x[2] == data[0]
                          and near(matmul(et, loc(x[3])), matmul(ft, loc(data[1])))]
                    if on:
                        if any(tables['Curve2ds'][p-1] not in ('line', 'circle') for p in on[0][1]):
                            ok = False
                    elif surf != 'plane':
                        ok = False
                    orients = sorted(vo for vo, _, _ in shapes[e][2])
                    if orients != ['+', '-']:
                        ok = False
                    for vo, v, vl in shapes[e][2]:
                        mark('Ve', v, matmul(et, loc(vl)))
    counts = tuple(len(seen[k]) for k in ('Ve', 'Ed', 'Wi', 'Fa', 'Sh', 'So'))
    return ok, counts


def files(root):
    return sorted(Path(root).glob('*.brep'))
