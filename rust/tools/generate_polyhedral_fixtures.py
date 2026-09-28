#!/usr/bin/env python3
"""Fixtures for S9b of REVIEW_NOTES.md: Booleans of polyhedral prisms in any
relative position.

`boolean-polyhedra-cases.txt` lists each case in the Boolean protocol
(`identity_reference.encode_boolean_case`: the object's prism, op 91, a
`boolean OP 93` row, the tool's prism, op 92); the two frames' axes differ
(S9b's domain): `XY`, the split fixtures' `TILT` (normal `(0, 3, 4)`), `ROT`
(the `XY` normal, x along `(3, 4, 0)`), `SIDE` (normal along x, x along y:
the world's axes permuted, so its faces lie on the object's planes) and
`LEAN` (normal `(3, 0, 4)`, x along y).

`boolean-polyhedra-expected.tsv` gives per case, from
`polyhedral_reference.py`:

* `expect KIND S9b`: the declared outcome (`solid`: a result of one or more
  solids; `empty`: no solids; `degenerate`: the result touches itself along
  an edge or at a point, `Degenerate` until the kernel holds non-manifold
  bodies, or holds a piece thinner than the resolution);
* `result N volume area cx cy cz` (totals over the N solids, world
  coordinates) or `empty`.

The measures are the exact models' (the stored axes as rationals, not
exactly orthonormal: a prism's volume is its profile's area times its height
times `det(x, y, n)`, within `1e-16` of the orthonormal one).

Before writing, `reference_checks` compares the reference with independent
computations: every polygon case of S9a's fixtures against
`boolean_reference.py`'s slicing (solid counts equal, measures within
`1e-15` relative: the stored axes' departure from orthonormal); every case's
convex cells against inclusion and exclusion (built in, exact); `area(A u
B) + area(A n B) = area(A) + area(B)` for every pair whose boundaries meet
nowhere face to face with opposite orientations; and the `SIDE` cases, whose
tool is a box with the world's axes, against the closed forms of two boxes.
"""
import argparse
from fractions import Fraction as F
from pathlib import Path

from identity_reference import Boundary, Case, encode_boolean_case
import polyhedral_reference as ref

ROOT = Path(__file__).resolve().parents[1]
XY = (0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0)
BOOLEAN_OPERATION = 93


def at(frame, origin):
    return tuple(origin)+tuple(frame)


TILT = (0.0, 3.0, 4.0, 1.0, 0.0, 0.0)
ROT = (0.0, 0.0, 1.0, 3.0, 4.0, 0.0)
SIDE = (1.0, 0.0, 0.0, 0.0, 1.0, 0.0)
LEAN = (3.0, 0.0, 4.0, 0.0, 1.0, 0.0)
UP = (0.0, 0.0, 1.0, 1.0, 0.0, 0.0)


def square(x0, y0, x1, y1):
    return Boundary(points=[(x0, y0), (x1, y0), (x1, y1), (x0, y1)])


def polygon(*points):
    return Boundary(points=list(points))


def prism(boundaries, frame, start, end, op):
    return Case('', 1e-7, op, frame, start, end, boundaries)


class Boolean:
    def __init__(self, name, operation, obj, tool, kind='solid', opposite=False):
        self.name, self.operation, self.kind, self.opposite = name, operation, kind, opposite
        self.obj = prism(*obj, op=91)
        self.obj.name = name
        self.tool = prism(*tool, op=92)
        self.tool.name = name

    def encode(self):
        return encode_boolean_case(self.obj, self.operation, self.tool, BOOLEAN_OPERATION)


def trio(name, obj, tool, kinds=('solid', 'solid', 'solid'), opposite=False):
    """The three operations of one pair (`kinds`: fuse, cut, common)."""
    return [Boolean(f'{name}_{op}', op, obj, tool, kind, opposite)
            for op, kind in zip(('fuse', 'cut', 'common'), kinds)]


def cases():
    """Every class of the S9b decisions: a box and a box turned about its
    axis, a tilted bar through a box (a cut in two), caps coplanar with the
    same and the opposite orientation, walls coplanar with the object's
    faces (the world's axes permuted) either way, an edge on a face, a vertex
    on a face, a tilted corner through an edge, a box inside another (a
    cavity), apart, an L profile against a tilted box, a bar through a hole
    without touching and a tool across the hole's walls, and both inputs
    turned."""
    box = ([square(0.0, 0.0, 10.0, 10.0)], XY, 0.0, 5.0)
    holed = ([square(0.0, 0.0, 10.0, 10.0), square(3.0, 3.0, 7.0, 7.0)], XY, 0.0, 5.0)
    ell = ([polygon((0.0, 0.0), (8.0, 0.0), (8.0, 3.0), (3.0, 3.0), (3.0, 8.0), (0.0, 8.0))],
           at(ROT, (1.0, 1.0, 0.0)), 0.0, 4.0)
    out = []
    out += trio('turned', box, ([square(0.0, 0.0, 6.0, 6.0)], at(ROT, (4.0, 2.0, 0.0)), -1.0, 4.0))
    out += trio('bar_tilted', box, ([square(0.0, -1.0, 12.0, 1.0)], at(TILT, (-1.0, 5.0, 2.5)), -5.0, 5.0))
    out += trio('caps_same', box, ([square(0.0, 0.0, 6.0, 6.0)], at(ROT, (4.0, 2.0, 0.0)), 0.0, 5.0))
    out += trio('caps_opposite', box, ([square(0.0, 0.0, 6.0, 6.0)], at(ROT, (4.0, 2.0, 0.0)), 5.0, 8.0),
                ('solid', 'solid', 'empty'), opposite=True)
    out += trio('walls_same', box, ([square(2.0, 0.0, 8.0, 3.0)], at(SIDE, (0.0, 0.0, 0.0)), -2.0, 5.0))
    out += trio('walls_opposite', box, ([square(2.0, 5.0, 8.0, 8.0)], at(SIDE, (0.0, 0.0, 0.0)), 2.0, 8.0),
                ('solid', 'solid', 'empty'), opposite=True)
    out += trio('edge_on_face', box, ([square(0.0, -3.0, 4.0, 0.0)], at(TILT, (3.0, 5.0, 5.0)), 0.0, 3.0),
                ('degenerate', 'solid', 'empty'))
    out += trio('vertex_on_face', box,
                ([polygon((0.0, 0.0), (-2.0, -4.0), (2.0, -4.0))], at(TILT, (5.0, 5.0, 5.0)), 0.0, 3.0),
                ('degenerate', 'solid', 'empty'))
    out += trio('corner_tilted', box, ([square(-2.0, -2.0, 2.0, 2.0)], at(LEAN, (10.0, 5.0, 5.0)), -2.0, 2.0))
    out += trio('inside', box, ([square(0.0, 0.0, 2.0, 2.0)], at(TILT, (4.0, 4.0, 1.5)), 0.0, 2.0))
    out += trio('apart', box, ([square(0.0, 0.0, 2.0, 2.0)], at(TILT, (14.0, 4.0, 1.5)), 0.0, 2.0),
                ('solid', 'solid', 'empty'))
    # The tool's far top edge passes within rounding of the L's corner edge
    # (`y = 4 + 1.6 + 1.8`): the fuse's and the cut's wall there is two
    # triangles joined by a neck thinner than the resolution.
    out += trio('ell_tilted', ell, ([square(0.0, -2.0, 6.0, 2.0)], at(TILT, (2.0, 4.0, 1.0)), -2.0, 3.0),
                ('degenerate', 'degenerate', 'solid'))
    wide = ([square(0.0, 0.0, 10.0, 10.0), square(2.0, 2.0, 8.0, 8.0)], XY, 0.0, 5.0)
    out += trio('through_hole', wide, ([square(4.0, -0.25, 6.0, 0.25)], at(TILT, (0.0, 2.5, 0.0)), -2.0, 8.5),
                ('solid', 'solid', 'empty'))
    out += trio('across_hole', holed, ([square(-2.0, -1.0, 14.0, 1.0)], at(ROT, (-1.0, 0.5, 1.0)), 0.0, 2.0))
    out += trio('both_turned', ([square(0.0, 0.0, 6.0, 6.0)], at(ROT, (0.0, 0.0, 0.0)), 0.0, 5.0),
                ([square(-3.0, -3.0, 3.0, 3.0)], at(LEAN, (2.0, 3.25, 2.5)), -2.0, 2.0))
    return out


# ------------------------------------------------------------------ checks

def box_measures(lo, hi):
    """Volume, area and centre of an axis-aligned box (Fractions)."""
    d = [hi[i]-lo[i] for i in range(3)]
    if min(d) <= 0:
        return F(0), F(0), None
    return d[0]*d[1]*d[2], 2*(d[0]*d[1]+d[1]*d[2]+d[0]*d[2]), tuple((lo[i]+hi[i])/2 for i in range(3))


def side_closed_form(case):
    """The SIDE cases: two axis-aligned boxes. Common exactly; fuse and cut
    volumes by inclusion and exclusion, areas by faces' pieces."""
    a, b = ref.Prism(case.obj), ref.Prism(case.tool)

    def bounds(p):
        pts = [pt for cell in p.cells for f in cell for pt in f]
        return [min(q[i] for q in pts) for i in range(3)], [max(q[i] for q in pts) for i in range(3)]

    (alo, ahi), (blo, bhi) = bounds(a), bounds(b)
    clo = [max(alo[i], blo[i]) for i in range(3)]
    chi = [min(ahi[i], bhi[i]) for i in range(3)]
    va, sa, ca = box_measures(alo, ahi)
    vb, sb, cb = box_measures(blo, bhi)
    vc, sc, cc = box_measures(clo, chi)
    return (va, sa, ca), (vb, sb, cb), (vc, sc, cc)


def reference_checks():
    """Largest deviations of each check (relative)."""
    import boolean_reference as slicing
    import generate_boolean_fixtures as s9a
    mp = ref.mp
    worst = {'s9a_slicing': mp.mpf(0), 'area_identity': mp.mpf(0), 'side_boxes': mp.mpf(0)}
    covered = {'s9a_slicing': 0, 'area_identity': 0, 'side_boxes': 0}
    rel = lambda x, y: abs(ref.M(x)-ref.M(y))/max(abs(ref.M(y)), 1)
    for case in s9a.cases():
        try:
            ref.Prism(case.obj)
            ref.Prism(case.tool)
        except AssertionError:
            continue
        a = slicing.rows(case.obj, case.operation, case.tool)[0][0].split()
        b = ref.rows(case.obj, case.operation, case.tool)[0][0].split()
        assert a[0] == b[0], (case.name, a, b)
        if a[0] == 'empty':
            covered['s9a_slicing'] += 1
            continue
        assert a[1] == b[1], ('solids', case.name, a, b)
        for x, y in zip(a[2:], b[2:]):
            worst['s9a_slicing'] = max(worst['s9a_slicing'], rel(mp.mpf(y), mp.mpf(x)))
        covered['s9a_slicing'] += 1
    assert worst['s9a_slicing'] < 1e-15, worst
    pairs = {}
    for case in cases():
        pairs.setdefault(case.name.rsplit('_', 1)[0], []).append(case)
    for name, trio_ in pairs.items():
        c = trio_[0]
        pair = ref.Pair(c.obj, c.tool)
        own = {}
        for label, prism_ in (('A', pair.A), ('B', pair.B)):
            own[label] = sum((ref.area(f) for f in prism_.faces), mp.mpf(0))
        areas = {op: pair.result(op)[2] for op in ('fuse', 'common')}
        if not c.opposite:
            dev = rel(areas['fuse']+areas['common'], own['A']+own['B'])
            worst['area_identity'] = max(worst['area_identity'], dev)
            covered['area_identity'] += 1
        if name.startswith('walls'):
            (va, sa, _), (vb, sb, _), (vc, sc, cc) = side_closed_form(c)
            for case in trio_:
                n, vol, surface, centre = pair.result(case.operation)
                want = {'common': vc, 'cut': va-vc, 'fuse': va+vb-vc}[case.operation]
                assert vol == want, (case.name, vol, want)
                if case.operation == 'common' and vc:
                    assert centre == cc and ref.M(surface) == ref.M(sc), case.name
                covered['side_boxes'] += 1
    assert worst['area_identity'] < 1e-35, worst
    return worst, covered


def generate():
    blocks = []
    out = ['# case\trow (S9b, polyhedral_reference.py: expect KIND S9b, then result N volume area cx cy cz '
           'or empty)']
    for case in cases():
        blocks.append(case.encode())
        rows, pair = ref.rows(case.obj, case.operation, case.tool)
        n = 0 if rows[0] == 'empty' else int(rows[0].split()[1])
        want = {'solid': n > 0, 'empty': n == 0, 'degenerate': n > 0}[case.kind]
        assert want, f'{case.name}: declared {case.kind}, the reference gives {rows[0]}'
        out.append(f'{case.name}\texpect {case.kind} S9b')
        for row in rows:
            out.append(f'{case.name}\t{row}')
    return {'boolean-polyhedra-cases.txt': '\n'.join(blocks)+'\n',
            'boolean-polyhedra-expected.tsv': '\n'.join(out)+'\n'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    names = [c.name for c in cases()]
    assert len(names) == len(set(names)), 'duplicate case names'
    worst, covered = reference_checks()
    files = generate()
    for name, contents in files.items():
        target = ROOT/'fixtures'/name
        if args.check:
            if target.read_text() != contents:
                raise SystemExit(f'{target} is stale')
        else:
            target.write_text(contents)
    listed = cases()
    print(len(listed), 'cases:', ', '.join(f'{sum(1 for c in listed if c.operation == op)} {op}'
                                           for op in ('fuse', 'cut', 'common')))
    print('reference checks (largest relative deviation):',
          ', '.join(f'{k} {ref.mp.nstr(v, 3)}' for k, v in sorted(worst.items())),
          '- cases:', ', '.join(f'{k} {v}' for k, v in sorted(covered.items())))


if __name__ == '__main__':
    main()
