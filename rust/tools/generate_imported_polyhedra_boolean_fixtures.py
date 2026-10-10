#!/usr/bin/env python3
"""Fixtures for S9e.4b.2 of REVIEW_NOTES.md: imported polyhedra other than
prisms (bodies of plane faces and line edges OCCT wrote to `.brep` files)
given to Booleans, decided on their stored vertices (S9b.2's stored model).

The bodies are OCCT's own output, as S9e.4a's
(`generate_imported_boolean_fixtures.py`): `boolean-imported-polyhedra-
bodies.txt` holds one `write NAME imported/NAME.brep` block per body for
`occt_boolean_oracle.cpp` (a `wedge` row, `BRepPrimAPI_MakeWedge`; a
`polyhedron` row, faces of given points sewn into a solid; a Boolean of two
`box` rows, its result's one solid; or `MakePrism` of a profile), which
`compare_imported_polyhedra_boolean.py --write-bodies` runs to write
`rust/fixtures/imported/NAME.brep` (format version 1, no triangulations).
Nothing here reads the kernel; the files are read only to check that each
is the body it claims to be and that it is this step's (below).

`boolean-imported-polyhedra-cases.txt` lists each case in the Boolean
protocol, an imported input as its one `brep imported/NAME.brep` row, a
chain's further `then` rows as S9e.1's; `boolean-imported-polyhedra-
expected.tsv` gives per case `expect KIND S9e.4b.2 CLASS` (`solid`,
`empty`, `degenerate`: a `Degenerate` of S9's rules; `unsupported`:
`OutOfDomain`, a later sub-step's; the class `polyhedron`: an imported
polyhedron and a construction, `both`: two imported inputs, `chain`: an
imported input's result given to another Boolean), a degenerate or
unsupported case's `reason TEXT`, then `result N volume area cx cy cz`
(totals over the N solids, world coordinates) or `empty`;
`boolean-imported-polyhedra-frames.tsv` the stored axes of every solid the
kernel builds from its rows.

The bodies: `tetra`, a tetrahedron of the local points `(0, 0, 0)`, `(8, 0,
0)`, `(2, 7, 0)`, `(3, 2, 6)` in the `TURN30` frame at `(2, 1, -1/2)`, its
world points rounded to binary64 (a `polyhedron`); `octa`, an octahedron
about `(5, 5, 4)` of half diagonals 4 along the axes (exact points, every
vertex of four faces); `pyramid`, `MakeWedge` in the `TILT` frame at `(1, 1,
0)` of base 8 by 8 and height 6, its top the point over the base's middle
(an apex of four faces); `truncated`, `MakeWedge` in the `TURN30` frame at
`(1, 1, 0)` of base 8 by 8, height 3 and top `[2, 6] x [2, 6]`; `wedge`,
`MakeWedge` in the `R125` frame at `(2, 2, 1)` of base 8 by 6, height 4 and
top `[1, 5] x [0, 4]`; `notched`, the box `[0, 10]^2 x [0, 4]` less a box of
6 cubed in the `SKEW` frame (normal `(1, 2, 2)`, `x` along `(2, -2, 1)`) at
`(8.2, 7.9, 2.4)` (a result OCCT computed);
`ell`, the box `[0, 6]^2 x [0, 3]` fused with a box of 6 by 4 by 3 in the
`TURN30` frame at `(3, 2, 2)`; `hollow`, the box `[0, 10]^3` less the box
`[3, 7]^3` (a cavity); and the DRAW survey's shapes from exact points:
`steps_low` and `steps_high`, two frustums, the second's base the first's
top (`buc60803a`, `b`); `pedestal`, a frustum whose base is a kernel box's
top, its coordinates decimals of 15 digits (`pro9481b` on `pro9481a`);
`draft`, a prism along `x` of a profile with a reflex corner, two of its
walls drafted (`CTO900_pro12559a`), and `ridge`, `MakePrism` of the lower
part of its far cap along `x` (`CTO900_pro12559b`, S9e.4a's prism); `vane_up`
and `vane_down`, two frustums by `MakeWedge` in frames turned by right
angles about `x`, their bases one square (`OCC578_w1`, `w2`).

Cases (each the three operations): `tetra_box`, `box_tetra`, the
tetrahedron and a box across one of its edges, either way; `octa_slab`, the
octahedron and a `TILT` slab; `pyramid_box`, the pyramid and a box across
its apex; `truncated_box`, the truncated pyramid and a box; `wedge_slab`;
`notched_rod`, the notched box and a rod of square section in the `TILT`
frame through its notch; `ell_slab`;
`tetra_octa`, both imported; `steps`, both imported, sharing a face;
`pedestal_box`, the kernel box and the frustum on its top; `draft_ridge`,
both imported, sharing a cap's part; `chain_vanes`, the vanes fused, then
the box `[0, 10]^2 x [0, 5]` with that result (swapped: the box the
object); `chain_notched`, the notched box less a vertical rod, then with a
`TILT` slab; `pyramid_flush`, a slab in the `TILTX` frame on the pyramid's
base's plane (on it in the construction, within the resolution of it in the
file, whose base corners are roundings), their fuse declared `degenerate`
(S9's rules: their faces meet within the resolution, a face using an edge
both ways), the cut and the common the construction's (two stored base
corners on the slab's plane, two off it on the pyramid's side by 7e-16); `hollow_slab`, the box
with a cavity and a `TILT` slab through it (declared `unsupported` until
S9e.4b.4c.1 took polyhedra with a cavity).

The reference is the construction OCCT was given
(`imported_polyhedra_boolean_reference.py`: every solid as convex cells in
exact Fractions, a Boolean's result as S9b's cells): volumes three ways
(inclusion and exclusion, the result's cells, the divergence over its
boundary), areas two ways, the pair identities on the cells' volumes and
the area identity where no faces share a plane, each solid's closed form
(the prismatoid formula, a tetrahedron's determinant, an octahedron's
diagonals) exactly, Monte Carlo against the inputs' half-spaces, solid
counts by cells meeting in area, and outside the declared group each
input's vertices and edges at least 1e-6 of the size from the other's faces
and edges where not on them exactly (an exact contact only between bodies
whose stored points are their construction's exactly) and crossing faces
at sines of at least 0.05. With `--check`, where the bodies' files exist, each is read
(`stored_records`): its faces' kinds planes, every stored vertex within
1e-12 of the size of the construction's boundary (an exact body's its
construction's points exactly), and the construction no prism (no pair of
opposite parallel faces with every other face along their normal): the
reason each body is this step's (but `ridge`, S9e.4a's prism).
"""
import argparse
from fractions import Fraction as F
import os
import struct
import zlib

import mpmath as mp

from identity_reference import Boundary, Case
from curve_surface_reference import stored_axes
import generate_imported_boolean_fixtures as base
import imported_polyhedra_boolean_reference as ref
import polyhedral_reference as pr
from generate_curved_boolean_fixtures import square

ROOT = base.ROOT
OPS = base.OPS
OPERATIONS, BOOLEANS = base.OPERATIONS, base.BOOLEANS
CLASSES = ('polyhedron', 'both', 'chain')
BODIES_FILE = 'boolean-imported-polyhedra-bodies.txt'
STEP = 'S9e.4b.2'
FRAMES = dict(base.FRAMES, SKEW=(1.0, 2.0, 2.0, 2.0, -2.0, 1.0))
prism = base.prism


def at(name, origin):
    return tuple(float(c) for c in origin)+FRAMES[name]


def numbers(*values):
    return ' '.join(repr(float(v)) for v in values)


def case_of(spec, op):
    """The identity case of a prism spec."""
    return base.construction(spec, op)


# ------------------------------------------------------------------ bodies

class Body:
    """An imported body: the rows OCCT builds it from, its construction as
    the reference's body (`build`), whether its stored points are its
    construction's exactly (`exact`), its closed-form volume (or None) and
    whether it is S9e.4a's prism (`prism`)."""

    def __init__(self, name, rows, build, exact=False, closed=None, is_prism=False):
        self.name, self.rows, self.build, self.exact = name, rows, build, exact
        self.closed, self.is_prism = closed, is_prism
        self.path = f'imported/{name}.brep'
        self._body = None

    def body(self):
        if self._body is None:
            self._body = self.build()
        return self._body

    def native_rows(self):
        return self.rows


def wedge(name, frame, dims, exact=False):
    o, n, x = frame[:3], frame[3:6], frame[6:9]
    return Body(name, ['wedge '+numbers(*o, *n, *x, *dims)],
                lambda: ref.Body([ref.hull(ref.wedge_points(frame, *dims))]), exact=exact,
                closed=lambda: ref.wedge_closed(frame, *dims))


def polyhedron(name, points, faces, cells=None, exact=False, closed=None):
    """A solid of binary64 points and faces (index cycles); its convex cells
    (point lists; all its points one cell by default)."""
    pts = [tuple(float(c) for c in p) for p in points]
    row = f'polyhedron {len(pts)} '+' '.join(numbers(*p) for p in pts)
    row += f' {len(faces)} '+' '.join(f'{len(f)} '+' '.join(str(i) for i in f) for f in faces)
    cells = cells or [pts]
    return Body(name, [row], lambda: ref.Body([ref.hull(c) for c in cells]), exact=exact, closed=closed)


def box_row(frame, size):
    o, n, x = frame[:3], frame[3:6], frame[6:9]
    return 'box '+numbers(*o, *n, *x, *size)


def box_spec(frame, size):
    """`MakeBox(gp_Ax2(origin, normal, x), DX, DY, DZ)`: the prism of `[0,
    DX] x [0, DY]` on the frame over `[0, DZ]`."""
    return prism([square(0.0, 0.0, size[0], size[1])], frame, 0.0, size[2])


def boxes(name, a, op, b):
    """A Boolean of two boxes OCCT computed (`(frame, size)` each)."""
    rows = [box_row(*a), f'boolean {op}', box_row(*b)]
    return Body(name, rows, lambda: ref.boolean(ref.prism_body(case_of(box_spec(*a), 0)),
                                                ref.prism_body(case_of(box_spec(*b), 0)), op)[0])


def local_points(frame, pts):
    """Local points on the frame's stored axes, rounded to binary64."""
    o, x, y, n = (tuple(F(c) for c in v) for v in stored_axes(frame))
    return [tuple(float(o[i]+F(u)*x[i]+F(v)*y[i]+F(w)*n[i]) for i in range(3)) for u, v, w in pts]


def frustum_faces(k):
    """A frustum's faces over `k` bottom points `0..k-1` and their top
    points `k..2k-1`."""
    return [list(range(k)), list(range(k, 2*k))]+[[i, (i+1) % k, k+(i+1) % k, k+i] for i in range(k)]


def tetra_det(p):
    a, b, c, d = (tuple(F(x) for x in q) for q in p)
    return abs(pr.det3(pr.sub(b, a), pr.sub(c, a), pr.sub(d, a)))/6


TETRA = local_points(at('TURN30', (2, 1, -0.5)), [(0, 0, 0), (8, 0, 0), (2, 7, 0), (3, 2, 6)])
OCTA = [(9, 5, 4), (1, 5, 4), (5, 9, 4), (5, 1, 4), (5, 5, 8), (5, 5, 0)]
OCTA_FACES = [[0, 2, 4], [2, 1, 4], [1, 3, 4], [3, 0, 4], [2, 0, 5], [1, 2, 5], [3, 1, 5], [0, 3, 5]]


def rect(x0, y0, x1, y1, z):
    return [(x0, y0, z), (x1, y0, z), (x1, y1, z), (x0, y1, z)]


SHIFT = (100.0, -40.0, 6.0)


def shifted(pts):
    return [tuple(p[i]+SHIFT[i] for i in range(3)) for p in pts]


STEPS_LOW = shifted(rect(0, 0, 16, 8, 0)+rect(1, 1, 15, 7, 2))
STEPS_HIGH = shifted(rect(1, 1, 15, 7, 2)+rect(3, 2, 13, 6, 3))
# `pro9481a`'s top and `pro9481b`'s, decimals of 15 digits.
PED = (-12.7239959266648, -9.76227793225116, 9.94168505699231, 12.6109683775700)
PED_TOP = (-2.93153733476602, 0.0301806596475927, 0.149226465093557, 2.81850978567122)
PEDESTAL = rect(*PED, 0.0)+rect(*PED_TOP, 21.0)
PEDESTAL_BOX = prism([square(*PED)], at('XY', (0, 0, -31)), 0.0, 31.0)
# The drafted prism's caps at x = 0 and x = 4 in (y, z): a profile with a
# reflex corner at (7, 7), its walls on lines of slopes 0, infinite and 1 and
# -1, two moved in at the far cap (y + z = 18 to 17.5, z = 9 to 8.5).
DRAFT0 = [(0, 0), (12, 0), (12, 6), (9, 9), (7, 9), (7, 7), (3, 7), (0, 4)]
DRAFT1 = [(0, 0), (12, 0), (12, 5.5), (9, 8.5), (7, 8.5), (7, 7), (3, 7), (0, 4)]
DRAFT = [(0.0, y, z) for y, z in DRAFT0]+[(4.0, y, z) for y, z in DRAFT1]
# Its convex cells: below and above z = 7 (the profiles split there at
# (11, 7) and (10.5, 7)).
DRAFT_CELLS = [(0.0, y, z) for y, z in [(0, 0), (12, 0), (12, 6), (11, 7), (3, 7), (0, 4)]] \
    + [(4.0, y, z) for y, z in [(0, 0), (12, 0), (12, 5.5), (10.5, 7), (3, 7), (0, 4)]]
DRAFT_TOP = [(0.0, y, z) for y, z in [(11, 7), (9, 9), (7, 9), (7, 7)]] \
    + [(4.0, y, z) for y, z in [(10.5, 7), (9, 8.5), (7, 8.5), (7, 7)]]
RIDGE = Boundary(points=[(0.0, 0.0), (12.0, 0.0), (12.0, 5.5), (10.5, 7.0), (3.0, 7.0), (0.0, 4.0)])


def draft_volume():
    # Two convex pieces, each with its caps' corresponding points.
    lo = ref.prismatoid(DRAFT_CELLS[:6], DRAFT_CELLS[6:], 0)
    hi = ref.prismatoid(DRAFT_TOP[:4], DRAFT_TOP[4:], 0)
    return lo+hi


def make_bodies():
    out = [
        polyhedron('tetra', TETRA, [[0, 2, 1], [0, 1, 3], [1, 2, 3], [2, 0, 3]], closed=lambda: tetra_det(TETRA)),
        polyhedron('octa', OCTA, OCTA_FACES, exact=True, closed=lambda: F(8*8*8, 6)),
        wedge('pyramid', at('TILT', (1, 1, 0)), (8.0, 6.0, 8.0, 4.0, 4.0, 4.0, 4.0)),
        wedge('truncated', at('TURN30', (1, 1, 0)), (8.0, 3.0, 8.0, 2.0, 2.0, 6.0, 6.0)),
        wedge('wedge', at('R125', (2, 2, 1)), (8.0, 4.0, 6.0, 1.0, 0.0, 5.0, 4.0)),
        boxes('notched', (at('XY', (0, 0, 0)), (10.0, 10.0, 4.0)), 'cut', (at('SKEW', (8.2, 7.9, 2.4)), (6.0, 6.0, 6.0))),
        boxes('ell', (at('XY', (0, 0, 0)), (6.0, 6.0, 3.0)), 'fuse', (at('TURN30', (3, 2, 2)), (6.0, 4.0, 3.0))),
        boxes('hollow', (at('XY', (0, 0, 0)), (10.0, 10.0, 10.0)), 'cut', (at('XY', (3, 3, 3)), (4.0, 4.0, 4.0))),
        polyhedron('steps_low', STEPS_LOW, frustum_faces(4), exact=True,
                   closed=lambda: ref.prismatoid(STEPS_LOW[:4], STEPS_LOW[4:], 2)),
        polyhedron('steps_high', STEPS_HIGH, frustum_faces(4), exact=True,
                   closed=lambda: ref.prismatoid(STEPS_HIGH[:4], STEPS_HIGH[4:], 2)),
        polyhedron('pedestal', PEDESTAL, frustum_faces(4), exact=True,
                   closed=lambda: ref.prismatoid(PEDESTAL[:4], PEDESTAL[4:], 2)),
        # Not convex: its cells its two convex pieces.
        polyhedron('draft', DRAFT, frustum_faces(8), cells=[DRAFT_CELLS, DRAFT_TOP], exact=True,
                   closed=draft_volume),
        Body('ridge', base.Body('ridge', prism([RIDGE], at('SIDE', (4, 0, 0)), 0.0, 6.0), 'prism').native_rows(),
             lambda: ref.prism_body(case_of(prism([RIDGE], at('SIDE', (4, 0, 0)), 0.0, 6.0), 0)), exact=True,
             is_prism=True),
        wedge('vane_up', (2.5, 7.5, 5.0, 0.0, -1.0, 0.0, 1.0, 0.0, 0.0), (5.0, 0.5, 5.0, 1.0, 1.0, 4.0, 4.0),
              exact=True),
        wedge('vane_down', (2.5, 2.5, 5.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0), (5.0, 3.0, 5.0, 1.0, 1.0, 4.0, 4.0),
              exact=True),
    ]
    return {b.name: b for b in out}


BODIES = make_bodies()


def imported(name):
    return ('imported', name)


BOX = prism([square(0.0, 0.0, 10.0, 10.0)], at('XY', (0, 0, 0)), 0.0, 5.0)
SLAB = prism([square(-4.0, -15.0, 14.0, 15.0)], at('TILT', (0, 0, 2.25)), 0.0, 1.5)
ROD = prism([square(4.0, 3.5, 9.5, 6.5)], at('XY', (0, 0, -1)), 0.0, 7.0)


def solid_of(spec):
    """A spec's reference body: an imported body's construction, or a
    prism's."""
    if spec[0] == 'imported':
        return BODIES[spec[1]].body()
    return ref.prism_body(case_of(spec, 0))


def spec_frame(spec):
    return spec[2]


class Imported:
    """One case: its solids (`imported(NAME)` or a prism), the first
    Boolean, the fixed stages (`(op, swapped)`) and the last stage's
    operation."""

    def __init__(self, name, klass, specs, op1, stages, kind, solids, reason=None):
        assert klass in CLASSES, klass
        self.name, self.klass, self.op1, self.stages = name, klass, op1, stages
        self.kind, self.solids, self.reason = kind, solids, reason
        self.group = name.rsplit('_', 1)[0]
        self.specs = specs
        self.solid_cases = []
        for k, s in enumerate(specs):
            if s[0] == 'imported':
                c = Case(name, 1e-7, OPERATIONS[k], at('XY', (0, 0, 0)), 0.0, 0.0, [], brep=BODIES[s[1]].path)
            else:
                c = case_of(s, OPERATIONS[k])
                c.name = name
            self.solid_cases.append(c)

    @property
    def last(self):
        return self.stages[-1][0] if self.stages else self.op1

    def encode(self):
        from identity_reference import encode_boolean_case, encode_chained_case
        c = self.solid_cases
        if not self.stages:
            return encode_boolean_case(c[0], self.op1, c[1], BOOLEANS[0])
        (op2, sw2) = self.stages[0]
        return encode_chained_case(c[0], self.op1, c[1], BOOLEANS[0], op2, c[2], BOOLEANS[1], sw2, None)

    def native(self):
        from identity_reference import native_boolean_case, native_chained_case
        c = self.solid_cases
        if not self.stages:
            return native_boolean_case(c[0], self.op1, c[1])
        (op2, sw2) = self.stages[0]
        return native_chained_case(c[0], self.op1, c[1], op2, c[2], sw2, None)


def group(name, klass, specs, outcomes, first=None, swapped=False, reason=None, kind=None, refused=None):
    """Cases of one group: the last operation varied over `outcomes` (`{op:
    solids}`, 0 for empty); with a `reason`, the operations in `refused`
    (all by default) declared `degenerate`."""
    out = []
    for op, n in outcomes.items():
        degenerate = reason and (refused is None or op in refused)
        k = kind or ('degenerate' if degenerate else ('solid' if n else 'empty'))
        stages = [(op, swapped)] if first else []
        out.append(Imported(f'{name}_{op}', klass, specs, first or op, stages, k, n, reason))
    return out


FLUSH = ('a slab on the plane of the imported pyramid\'s base (on it in the construction, its stored corners '
         'within the resolution of it): their fuse')


def cases():
    out = []
    three = {'fuse': 1, 'cut': 1, 'common': 1}
    tetra_box = prism([square(5.0, -2.0, 12.0, 5.0)], at('XY', (0, 0, 1)), 0.0, 3.0)
    out += group('tetra_box', 'polyhedron', [imported('tetra'), tetra_box], three)
    out += group('box_tetra', 'polyhedron', [tetra_box, imported('tetra')], three)
    out += group('octa_slab', 'polyhedron', [imported('octa'), prism([square(-4.0, -15.0, 14.0, 15.0)],
                                                                     at('TILT', (0, 0, 0)), 4.5, 5.5)],
                 {'fuse': 1, 'cut': 2, 'common': 1})
    out += group('pyramid_box', 'polyhedron', [imported('pyramid'),
                                               prism([square(2.0, 4.5, 8.0, 10.5)], at('XY', (0, 0, 3)), 0.0, 4.0)],
                 three)
    out += group('truncated_box', 'polyhedron', [imported('truncated'),
                                               prism([square(-1.5, 3.5, 2.5, 7.5)], at('XY', (0, 0, 2.5)), 0.0, 3.0)],
                 three)
    out += group('wedge_slab', 'polyhedron', [imported('wedge'), SLAB], three)
    out += group('notched_rod', 'polyhedron', [imported('notched'),
                                               prism([square(7.0, -9.0, 9.0, 13.0)], at('TILT', (0, 4, 4)),
                                                     -3.0, 1.0)], three)
    out += group('ell_slab', 'polyhedron', [imported('ell'), SLAB], {'fuse': 1, 'cut': 2, 'common': 1})
    out += group('tetra_octa', 'both', [imported('tetra'), imported('octa')], three)
    out += group('steps', 'both', [imported('steps_low'), imported('steps_high')],
                 {'fuse': 1, 'cut': 1, 'common': 0})
    out += group('pedestal_box', 'polyhedron', [PEDESTAL_BOX, imported('pedestal')],
                 {'fuse': 1, 'cut': 1, 'common': 0})
    out += group('draft_ridge', 'both', [imported('draft'), imported('ridge')], {'fuse': 1, 'cut': 1, 'common': 0})
    out += group('chain_vanes', 'chain', [imported('vane_up'), imported('vane_down'), BOX], three, first='fuse',
                 swapped=True)
    out += group('chain_notched', 'chain', [imported('notched'), ROD, SLAB], {'fuse': 1, 'cut': 2, 'common': 1},
                 first='cut')
    out += group('pyramid_flush', 'polyhedron', [imported('pyramid'),
                                                 prism([square(-2.0, -2.0, 10.0, 10.0)], at('TILTX', (1, 1, 0)),
                                                       0.0, 2.0)],
                 {'fuse': 2, 'cut': 1, 'common': 0}, reason=FLUSH, refused=('fuse',))
    out += group('hollow_slab', 'polyhedron', [imported('hollow'), SLAB], {'fuse': 1, 'cut': 2, 'common': 1})
    return out


def all_cases():
    return cases()


def validate(listed):
    names = [c.name for c in listed]
    assert len(names) == len(set(names)), 'duplicate case names'
    for c in listed:
        assert len(c.specs) == 2+len(c.stages), c.name
        assert c.last in OPS, c.name
        assert any(s[0] == 'imported' for s in c.specs), c.name
        if c.klass == 'both':
            assert sum(s[0] == 'imported' for s in c.specs) == 2, c.name
        if c.klass == 'chain':
            assert c.stages, c.name
        for s in c.specs:
            if s[0] != 'imported':
                assert tuple(spec_frame(s)[3:]) in FRAMES.values(), c.name


# ------------------------------------------------------------------ the reference

def exact_frame(spec):
    """Whether a prism's stored axes are exact (every coordinate 0 or +-1)."""
    return all(c in (0.0, 1.0, -1.0) for v in stored_axes(spec_frame(spec))[1:] for c in v)


def evaluate(job):
    """One group: its rows, results and checks."""
    name, specs, op1, swapped, chained, ops, degenerate, mc_n = job
    inputs = [solid_of(s) for s in specs]
    size = max(b.size() for b in inputs)
    rows, res, checks = {}, {}, {}

    def note(key, value):
        checks[key] = max(checks.get(key, mp.mpf(0)), pr.M(value) if isinstance(value, F) else value)
    pairs = {}
    if chained:
        n1, v1, _, _, first, c1 = ref.result(inputs[0], inputs[1], op1)
        assert n1 == 1, f'{name}: the first result is {n1} solids'
        for k, v in c1.items():
            note(k, v)
        x, y = (inputs[2], first) if swapped else (first, inputs[2])
        exprs = {op: (op, 2, (op1, 0, 1)) if swapped else (op, (op1, 0, 1), 2) for op in ops}
    else:
        x, y = inputs
        exprs = {op: (op, 0, 1) for op in ops}
    for op, solids in ops.items():
        n, vol, area, centre, body, c = ref.result(x, y, op)
        for k, v in c.items():
            note(k, v)
        assert n == solids, f'{name} {op}: {n} solids, {solids} declared'
        rows[op] = ref.rows(n, vol, area, centre)
        res[op] = (n, vol, area, centre)
        pairs[op] = body
    # The pair identities on the results' own cells.
    if set(ops) == set(OPS):
        vx, mx = x.measure()
        vy, my = y.measure()
        vf, mf = pairs['fuse'].measure()
        vn, mn = pairs['common'].measure()
        vt, mt = pairs['cut'].measure()
        dev = max([abs(vf+vn-vx-vy)/size**3, abs(vt-(vx-vn))/size**3]
                  + [abs(mf[i]+mn[i]-mx[i]-my[i])/size**4 for i in range(3)]
                  + [abs(mt[i]-(mx[i]-mn[i]))/size**4 for i in range(3)])
        note('pair_identities', dev)
        if not ref.shares_plane(x, y):
            ax = sum((pr.area(f) for f in x.faces), mp.mpf(0))
            ay = sum((pr.area(f) for f in y.faces), mp.mpf(0))
            note('area_identity', abs(res['fuse'][2]+res['common'][2]-ax-ay)/size**2)
    # Each imported input's closed form.
    for s, b in zip(specs, inputs):
        if s[0] == 'imported' and BODIES[s[1]].closed is not None:
            note('solid_closed_forms', abs(b.measure()[0]-BODIES[s[1]].closed())/size**3)
    # Monte Carlo against the inputs' half-spaces.
    z = 0.0
    for op in ops:
        vol, c, sv, sc = ref.monte_carlo(inputs, exprs[op], mc_n, zlib.crc32(f'{name} {op}'.encode()))
        n, v, _, centre = res[op]
        z = max(z, abs(vol-float(v))/max(sv, 1e-12))
        if c is not None and n:
            z = max(z, max(abs(c[i]-float(centre[i]))/max(sc[i], 1e-12) for i in range(3)))
    note('monte_carlo_sigma', mp.mpf(z))
    clearance, contacts, sine = ref.margins(x, y)
    if chained:
        c0, k0, s0 = ref.margins(inputs[0], inputs[1])
        clearance, contacts, sine = min(clearance, c0), contacts+k0, min(sine, s0)
    exact = all((s[0] == 'imported' and BODIES[s[1]].exact) or (s[0] != 'imported' and exact_frame(s))
                for s in specs)
    if contacts and not degenerate:
        assert exact, f'{name}: an exact contact between bodies not stored exactly'
    margins = {'clearance': clearance, 'sine': sine, 'contacts': contacts}
    return name, rows, res, checks, margins, degenerate


def jobs(mc_n):
    groups = {}
    for c in all_cases():
        e = groups.setdefault(c.group, [c.specs, c.op1, bool(c.stages and c.stages[0][1]), bool(c.stages), {},
                                        c.kind in ('degenerate', 'unsupported')])
        e[4][c.last] = c.solids
    return [(name, specs, op1, swapped, chained, ops, degenerate, mc_n)
            for name, (specs, op1, swapped, chained, ops, degenerate) in groups.items()]


def cached(job):
    """`evaluate`, its result kept in `--cache` (a development aid)."""
    cache = os.environ.get('IMPORTED_POLYHEDRA_CACHE')
    if not cache:
        return evaluate(job)
    import hashlib
    import inspect
    import pickle
    from pathlib import Path
    key = hashlib.sha256((repr(job)+inspect.getsource(evaluate)+Path(ref.__file__).read_text()
                          + Path(__file__).read_text()).encode()).hexdigest()[:24]
    path = Path(cache)/f'{job[0]}-{key}.pickle'
    if path.exists():
        return pickle.loads(path.read_bytes())
    try:
        out = evaluate(job)
    except Exception as e:
        raise RuntimeError(f'{job[0]}: {e!r}') from e
    path.write_bytes(pickle.dumps(out))
    return out


LIMITS = {'volume_three_ways': 0, 'area_two_ways': 1e-30, 'pair_identities': 0, 'area_identity': 1e-30,
          'solid_closed_forms': 0, 'monte_carlo_sigma': 5, 'least_sine': 20, 'clearance': 1e6}


def reference_checks(results):
    worst, covered = {}, {}

    def note(key, value):
        worst[key] = max(worst.get(key, mp.mpf(0)), value)
        covered[key] = covered.get(key, 0)+1

    contacts = 0
    for name, rows, _, checks, margins, degenerate in results:
        for key, value in checks.items():
            note(key, value)
        if not degenerate:
            note('least_sine', 1/max(margins['sine'], mp.mpf(10)**-30))
            note('clearance', 1/max(margins['clearance'], mp.mpf(10)**-40))
            contacts += margins['contacts']
    return worst, covered, contacts


# ------------------------------------------------------------------ the files

def bodies_text():
    """The writer's blocks: each body's rows under `write NAME PATH`."""
    out = []
    for b in BODIES.values():
        out += [f'write {b.name} {b.path}', *b.native_rows(), 'end']
    return '\n'.join(out)+'\n'


def generate(results):
    by_group = {r[0]: r for r in results}
    blocks = []
    out = [f'# case\trow ({STEP}; imported_polyhedra_boolean_reference.py: expect KIND {STEP} CLASS, reason TEXT '
           'for a degenerate or unsupported case, then result N volume area cx cy cz or empty)']
    for case in all_cases():
        blocks.append(case.encode())
        row = by_group[case.group][1][case.last]
        n = 0 if row[0] == 'empty' else int(row[0].split()[1])
        want = {'solid': n > 0, 'empty': n == 0, 'degenerate': True, 'unsupported': True}[case.kind]
        assert want, f'{case.name}: declared {case.kind}, the reference gives {row[0]}'
        out.append(f'{case.name}\texpect {case.kind} {STEP} {case.klass}')
        if case.kind in ('degenerate', 'unsupported'):
            out.append(f'{case.name}\treason {case.reason}')
        for r in row:
            out.append(f'{case.name}\t{r}')
    frames = ['# case\tframe (solid k) axis (n, x, y)\tstored unit vector (stored_axes, as hex bits; the solids '
              'the kernel builds from rows, not the imported ones)']
    for case in all_cases():
        for k, spec in enumerate(case.specs):
            if spec[0] == 'imported':
                continue
            _, x, y, n = stored_axes(spec_frame(spec))
            for key, v in (('n', n), ('x', x), ('y', y)):
                frames.append(f'{case.name}\tsolid{k} {key}\t'+' '.join(struct.pack('>d', q).hex() for q in v))
    return {'boolean-imported-polyhedra-cases.txt': '\n'.join(blocks)+'\n',
            'boolean-imported-polyhedra-expected.tsv': '\n'.join(out)+'\n',
            'boolean-imported-polyhedra-frames.tsv': '\n'.join(frames)+'\n',
            BODIES_FILE: bodies_text()}


def native_input():
    """Every case's native rows (`compare_imported_polyhedra_boolean.py`)."""
    return '\n'.join(c.native() for c in all_cases())+'\n'


def case_size(case):
    """The case's size: its solids' reach from the origin."""
    return float(max(solid_of(s).size() for s in case.specs))


# ------------------------------------------------------------------ the files' check

def is_prism(body):
    """Whether a body's boundary is a prism's: a pair of planes of its
    faces with opposite normals, every other plane's normal perpendicular to
    theirs."""
    planes = []
    for f in body.faces:
        n = pr.normal_of(f)
        k = max(abs(x) for x in n)
        key = (tuple(x/k for x in n), pr.dot(n, f[0])/k)
        if key not in planes:
            planes.append(key)
    for a in planes:
        for b in planes:
            if pr.cross(a[0], b[0]) == (0, 0, 0) and pr.dot(a[0], b[0]) < 0:
                if all(pr.dot(c[0], a[0]) == 0 for c in planes if c not in (a, b)):
                    return True
    return False


def check_bodies():
    """Each existing body file: its faces planes, each stored vertex within
    1e-12 of the size of the construction's boundary (an exact body's its
    construction's points exactly), the construction no prism (but S9e.4a's
    `ridge`). Returns the largest deviation and the files read."""
    worst, read = mp.mpf(0), 0
    for b in BODIES.values():
        path = ROOT/'fixtures'/b.path
        if not path.exists():
            continue
        read += 1
        kinds, vertices = base.stored_records(path.read_text())
        assert set(kinds) == {'plane'}, (b.name, kinds)
        body = b.body()
        assert is_prism(body) == b.is_prism, (b.name, 'a prism' if b.is_prism else 'not a prism')
        size = pr.M(body.size())
        points = {tuple(float(c) for c in p) for p in body.points()}
        for v in vertices:
            if b.exact:
                assert v in points, (b.name, v)
            d = min(ref.point_face_distance(v, f) for f in body.faces)/size
            worst = max(worst, d)
    assert worst <= 1e-12, ('a stored vertex off its construction', mp.nstr(worst, 3))
    return worst, read


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--workers', type=int, default=max(1, min(4, (os.cpu_count() or 2)-2)))
    parser.add_argument('--samples', type=int, default=200000, help='Monte-Carlo points per operation')
    parser.add_argument('--only', help='one group, printed (no files)')
    parser.add_argument('--cache', help='a directory keeping each group\'s result (development)')
    args = parser.parse_args()
    if args.cache:
        os.environ['IMPORTED_POLYHEDRA_CACHE'] = args.cache
    listed = all_cases()
    validate(listed)
    todo = jobs(args.samples)
    if args.only:
        todo = [j for j in todo if j[0] == args.only]
        for r in base.run(todo, cached, 1):
            print(r)
        return
    results = base.run(todo, cached, args.workers)
    worst, covered, contacts = reference_checks(results)
    for key, value in worst.items():
        assert value <= LIMITS[key], (key, mp.nstr(value, 3))
    files = generate(results)
    for name, contents in files.items():
        target = ROOT/'fixtures'/name
        if args.check:
            if target.read_text() != contents:
                raise SystemExit(f'{target} is stale')
        else:
            target.write_text(contents)
    deviation, read = check_bodies()
    kinds = {}
    for c in listed:
        kinds[c.kind] = kinds.get(c.kind, 0)+1
    print(len(listed), 'cases:', ', '.join(f'{sum(1 for c in listed if c.last == op)} {op}' for op in OPS),
          '-', ', '.join(f'{v} {k}' for k, v in sorted(kinds.items())),
          '-', ', '.join(f'{sum(1 for c in listed if c.klass == k)} {k}' for k in CLASSES))
    print('reference checks (largest deviation, relative to the case size; least_sine and clearance as '
          'reciprocals):', ', '.join(f'{k} {mp.nstr(v, 3)}' for k, v in sorted(worst.items())))
    print('cases per check:', ', '.join(f'{k} {v}' for k, v in sorted(covered.items())),
          f'- exact contacts {contacts}')
    print(f'bodies: {read} of {len(BODIES)} files read, every stored vertex within {mp.nstr(deviation, 3)} of the '
          'size of its construction\'s boundary')


if __name__ == '__main__':
    main()
