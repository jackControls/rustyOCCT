#!/usr/bin/env python3
"""Fixtures for S9e.4b.3c.2 of REVIEW_NOTES.md: exact incidences of two
pieces of one sphere given to a Boolean (a vertex, a circle or a line of
both inputs, plane faces on one plane with overlapping edges: the DRAW
survey's `so1` and `so2`, `so2` and `so3`, `so5` and `so2`).

The bodies are OCCT's own output, as S9e.4b.3c.1's
(`generate_one_sphere_boolean_fixtures.py`): `boolean-one-sphere-incidence-
bodies.txt` holds one `write NAME imported/NAME.brep` block per body for
`occt_boolean_oracle.cpp` (a sphere's row, a `boolean common` row and a
`box` row, `BRepPrimAPI_MakeBox` on a frame; the hemisphere's rim divided by
a `divide` row as `so1`'s), which `compare_one_sphere_incidence_boolean.py
--write-bodies` runs to write `rust/fixtures/imported/NAME.brep`. Nothing
here reads the kernel; the files are read only to check that each is the
body it claims to be and that its stored data hold the incidences exactly.

`boolean-one-sphere-incidence-cases.txt` lists each case in the Boolean
protocol (an imported input its one `brep imported/NAME.brep` row, a
chain's further `then` row as S9e.1's); `boolean-one-sphere-incidence-
expected.tsv` gives per case `expect KIND S9e.4b.3c.2 CLASS` (`solid`,
`empty`, `degenerate`: a `Degenerate` of S9's rules; the class `pieces`: two
imported pieces of one sphere, `sphere`: a constructed cap of the sphere
and a piece or another cap, `chain`: two pieces' result given to a
Boolean with a zone of their sphere), a degenerate case's `reason TEXT`,
then `result N volume area cx cy cz` or `empty`;
`boolean-one-sphere-incidence-frames.tsv` the stored axes of every solid the
kernel builds from its rows.

Every body is a piece of the ball of radius 5 about `(5, 5, 4)` cut by
planes holding a world axis through its centre or normal to it, as the
survey's files are: an exact incidence of two inputs' planes along a line
holds in the stored data only where their normals' zero components survive
OCCT's normalization and the kernel's (`gp_Dir`, `Frame3::new`), so on the
world's axes; in turned rational frames two inputs' planes are one plane
only where stored bit for bit alike (S9e.4b.3c.1's `hemi_octant`, the
hemisphere and the octant on one frame). Each body is the ball (its
sphere's frame on the axis) common `MakeBox` on a frame whose normal is the
axis and whose origin lies on it: the box's faces through its origin are
the axis's half-planes and its parallels' planes (their stored normals
with exact zeros), the far faces apart from the ball. About the world's `z`:
`incidence_hemi`, above the centre's parallel (the box of side 16 about the
axis), its rim divided as `so1`'s; `incidence_wedge45`, between the
half-planes at 45 and 135 degrees (`so2`'s); `incidence_wedge23`, between
those at `atan2(5, 12)` and 90 degrees more (`so3`'s); `incidence_high23`,
the same above the parallel's plane at `5/2` (`so5`'s); `incidence_half`,
on one side of the meridian plane along `x` (two vertices of the rim, the
semicircle through the pole and the diameter its edges);
`incidence_octant`, between the half-planes at 0 and 90 degrees, and
`incidence_back`, at 180 and 270 degrees; the near ones, declared
`degenerate`: `incidence_turned45`, `incidence_wedge45`'s with its frame's
`x` turned by about `2^-36` (each half-plane within the resolution of the
other's over the ball), and `incidence_lifted45`, its box `2^-40` above the
centre (its base within the resolution of the hemisphere's). About the
world's `x`: `incidence_x45` and `incidence_x23`, the wedges of the same
angles about it. The constructed solids: caps above the parallels at
latitudes `1/4` and `1/2` and a zone between them, on the world's axes at
the centre (S9d.1's spheres).

Cases (each the three operations, the cut either way where it differs):
`hemi_wedge`, `wedge_hemi` (`so1` and `so2`: the wedge's equator arc on the
hemisphere's rim circle, its base on the hemisphere's base, its corner and
pole on the hemisphere's faces; the wedge's cut by the hemisphere empty);
`wedge_wedge`, `wedge_back` (`so2` and `so3`: the corner, the pole and the
axis edge of both, the base faces on one plane with overlapping edges, each
equator arc's end inside the other's, the meridian half-planes crossing
along the axis); `high_wedge`, `wedge_high` (`so5` and `so2`: the higher
wedge's axis edge inside the other's, its corner on it); `half_octant`,
`octant_half` (the octant inside the half: a vertex of both on the rim,
the octant's corner inside the half's diameter, its pole inside the half's
semicircle, an arc and a radius of the octant on the half's); `x_wedges`
(the wedges about the world's `x`); `same_wedge`, the wedge and itself
(every vertex, edge and face of both; its cut empty); `cap_cap`, `cap_back` (two caps on one
frame, their stored pole one vertex); `cap_wedge`, `wedge_cap` (a cap and
the wedge, the pole of both); `chain_zone`, the wedges' fuse then with the
zone (its cut two solids). Declared `degenerate`: `quadrants`, the octant
and the opposite octant touching along the axis only (their fuse a result
touching itself along it; S9's near-plane guard refuses their faces meeting
only there, their boxes' overlap thinner than the resolution: the cut the
octant and the common empty by the reference); `half_wedge`, the half and
the wedge (the wedge's axis edge inside the half's meridian face, off its
edges: an edge of one input on a face of the other, S9's; the cut also a
result touching itself along it); `near_turned`, the wedge and
its turned copy, and `near_lifted`, the hemisphere and the lifted wedge (two
faces within the resolution of one plane).

The reference is exact and of its own (independent of the kernel and of
the chained reference, whose sweeps follow the world's meridians and
parallels, which these bodies' faces lie on): every input of a case is the
ball, its heights along the common axis in an interval and its angles about
the axis in the half-turns of its half-planes through the axis (each read
off the rows, asserted to hold the axis or to be normal to it or to miss the
ball); the ball is cut into cells by every input's heights and every
half-plane's angle, each cell kept by the case's set function of the
inputs' memberships at its middle; each kept cell's volume, first moments
and sphere area in closed form (mpmath, rounding once), the faces between
kept and dropped cells (the parallels' sectors, the half-planes' pieces)
added to the area. Checks: the volume again by the divergence theorem over
those faces (the sphere's zones and the parallels' sectors; the half-planes
through the centre add none), the pair identities, every Monte-Carlo point's
membership by the rows' own inequalities against its cell's, the volume and
centre by Monte Carlo (100,000 points a group), the solid counts by the
kept cells' adjacency (a result whose kept cells about the axis are apart
touches itself along it: declared `degenerate`), and with `--check` each
body's file: its faces one sphere and planes, every stored vertex on the
construction's surfaces within 1e-12 of the size, every stored plane
holding the axis or normal to it exactly or apart from the ball (the
incidences exact in the stored data), the divided rim two arcs.
"""
import argparse
from fractions import Fraction as F
import os
import random
import struct
import zlib

import mpmath as mp

from identity_reference import Case
from curve_surface_reference import stored_axes
from torus_boolean_reference import number
import brep_io_reference as bio
import generate_imported_boolean_fixtures as base

mp.mp.dps = 40
ROOT = base.ROOT
OPS = base.OPS
OPERATIONS, BOOLEANS = base.OPERATIONS, base.BOOLEANS
HP = base.HP
STEP = 'S9e.4b.3c.2'
CLASSES = ('pieces', 'sphere', 'chain')
BODIES_FILE = 'boolean-one-sphere-incidence-bodies.txt'
PREFIX = 'boolean-one-sphere-incidence'
CENTRE = (5.0, 5.0, 4.0)
RADIUS = 5.0
# Exact coincidences of the cells' breakpoints (angles of one half-plane
# read twice, the sphere's own heights): closer than this are one.
MERGE = mp.mpf(10)**-30
TWO_PI = 2*mp.pi


# ------------------------------------------------------------------ axes

# A world axis: (e, u, v), `e` the axis, `u` the angles' zero, `v = e x u`.
AXES = {
    'z': ((0, 0, 1), (1, 0, 0), (0, 1, 0)),
    'x': ((1, 0, 0), (0, 1, 0), (0, 0, 1)),
}


def on_axis(axis, h, u=0.0, v=0.0):
    """The point at height `h` along the axis through the centre, `u` and
    `v` across it (binary64, exact for these values)."""
    e, uu, vv = AXES[axis]
    return tuple(CENTRE[i]+h*e[i]+u*uu[i]+v*vv[i] for i in range(3))


def along(axis, a, b):
    """The direction `a u + b v` across the axis (a frame's `x`)."""
    _, uu, vv = AXES[axis]
    return tuple(a*uu[i]+b*vv[i] for i in range(3))


# ------------------------------------------------------------------ bodies

class Body:
    """An imported piece: the ball on the axis common `MakeBox(gp_Ax2(o,
    e, x), d)`, its rim divided where `divide`."""

    def __init__(self, name, axis, o, x, d, divide=False):
        self.name, self.axis, self.o, self.x, self.d, self.divide = name, axis, o, x, d, divide
        self.path = f'imported/{name}.brep'

    def sphere_row(self):
        e, u, _ = AXES[self.axis]
        values = (*CENTRE, *e, *u, RADIUS, -HP, HP)
        return 'sphere '+' '.join(repr(float(v)) for v in values)

    def box_row(self):
        e, _, _ = AXES[self.axis]
        values = (*self.o, *e, *self.x, *self.d)
        return 'box '+' '.join(repr(float(v)) for v in values)

    def native_rows(self):
        return [self.sphere_row(), 'boolean common', self.box_row()]+(['divide'] if self.divide else [])


def make_bodies():
    out = [
        Body('incidence_hemi', 'z', on_axis('z', 0.0, -8.0, -8.0), along('z', 1, 0), (16.0, 16.0, 8.0),
             divide=True),
        Body('incidence_wedge45', 'z', on_axis('z', 0.0), along('z', 1, 1), (8.0, 8.0, 8.0)),
        Body('incidence_wedge23', 'z', on_axis('z', 0.0), along('z', 12, 5), (8.0, 8.0, 8.0)),
        Body('incidence_high23', 'z', on_axis('z', 2.5), along('z', 12, 5), (8.0, 8.0, 8.0)),
        Body('incidence_half', 'z', on_axis('z', 0.0, -8.0), along('z', 1, 0), (16.0, 8.0, 8.0)),
        Body('incidence_octant', 'z', on_axis('z', 0.0), along('z', 1, 0), (8.0, 8.0, 8.0)),
        Body('incidence_back', 'z', on_axis('z', 0.0), along('z', -1, 0), (8.0, 8.0, 8.0)),
        Body('incidence_turned45', 'z', on_axis('z', 0.0), along('z', 1, 1+2.0**-36), (8.0, 8.0, 8.0)),
        Body('incidence_lifted45', 'z', on_axis('z', 2.0**-40), along('z', 1, 1), (8.0, 8.0, 8.0)),
        Body('incidence_x45', 'x', on_axis('x', 0.0), along('x', 1, 1), (8.0, 8.0, 8.0)),
        Body('incidence_x23', 'x', on_axis('x', 0.0), along('x', 12, 5), (8.0, 8.0, 8.0)),
    ]
    return {b.name: b for b in out}


BODIES = make_bodies()


def imported(name):
    return ('imported', name)


def cap(low, high=HP):
    """S9d.1's sphere on the world's axes at the centre between two
    latitudes."""
    return base.sphere(RADIUS, base.at('XY', CENTRE), low, high)


CAP_LOW, CAP_HIGH = cap(0.25), cap(0.5)
ZONE = cap(0.25, 0.5)


# ------------------------------------------------------------------ cases

class Incident:
    """One case: its inputs (`imported(NAME)` or a construction), the first
    Boolean, the fixed stages (`(op, swapped)`) and the last stage's
    operation."""

    def __init__(self, name, klass, items, op1, stages, kind, solids, reason=None):
        assert klass in CLASSES, klass
        self.name, self.klass, self.op1, self.stages = name, klass, op1, stages
        self.kind, self.solids, self.reason = kind, solids, reason
        self.group = name.rsplit('_', 1)[0]
        self.items = items
        self.solid_cases = []
        for k, s in enumerate(items):
            if s[0] == 'imported':
                c = Case(name, 1e-7, OPERATIONS[k], base.at('XY', (0, 0, 0)), 0.0, 0.0, [],
                         brep=BODIES[s[1]].path)
            else:
                c = base.construction(s, OPERATIONS[k])
                c.name = name
            self.solid_cases.append(c)

    @property
    def last(self):
        return self.stages[-1][0] if self.stages else self.op1

    def expr(self):
        if not self.stages:
            return (self.op1, 0, 1)
        (op2, sw2) = self.stages[0]
        first = (self.op1, 0, 1)
        return (op2, 2, first) if sw2 else (op2, first, 2)

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


def group(name, klass, items, outcomes, first=None, reason=None):
    """Cases of one group: the last operation varied over `outcomes` (`{op:
    solids}`: 0 for empty, `None` for declared `degenerate` with `reason`)."""
    out = []
    for op, n in outcomes.items():
        kind = 'degenerate' if n is None else ('solid' if n else 'empty')
        stages = [(op, False)] if first else []
        out.append(Incident(f'{name}_{op}', klass, items, first or op, stages, kind, n,
                            reason if n is None else None))
    return out


INSIDE = ('an edge of one input on a face of the other (the wedge\'s axis edge inside the half\'s meridian face, '
          'off its edges)')
NEAR = 'two faces within the resolution of one plane (a half-plane turned or a base lifted by rounding)'
TOUCH = ('two octants about one axis touching along it only (their fuse a result touching itself there, their faces '
         'meeting only there within the resolution of one plane over their overlap)')


def cases():
    out = []
    three = {'fuse': 1, 'cut': 1, 'common': 1}
    wedge45, wedge23 = imported('incidence_wedge45'), imported('incidence_wedge23')
    out += group('hemi_wedge', 'pieces', [imported('incidence_hemi'), wedge45], three)
    out += group('wedge_hemi', 'pieces', [wedge45, imported('incidence_hemi')], {'cut': 0})
    out += group('wedge_wedge', 'pieces', [wedge45, wedge23], three)
    out += group('wedge_back', 'pieces', [wedge23, wedge45], {'cut': 1})
    out += group('high_wedge', 'pieces', [imported('incidence_high23'), wedge45], three)
    out += group('wedge_high', 'pieces', [wedge45, imported('incidence_high23')], {'cut': 1})
    out += group('half_octant', 'pieces', [imported('incidence_half'), imported('incidence_octant')], three)
    out += group('octant_half', 'pieces', [imported('incidence_octant'), imported('incidence_half')], {'cut': 0})
    out += group('x_wedges', 'pieces', [imported('incidence_x45'), imported('incidence_x23')], three)
    out += group('same_wedge', 'pieces', [wedge45, wedge45], {'fuse': 1, 'cut': 0, 'common': 1})
    out += group('cap_cap', 'sphere', [CAP_LOW, CAP_HIGH], three)
    out += group('cap_back', 'sphere', [CAP_HIGH, CAP_LOW], {'cut': 0})
    out += group('cap_wedge', 'sphere', [CAP_HIGH, wedge45], three)
    out += group('wedge_cap', 'sphere', [wedge45, CAP_HIGH], {'cut': 1})
    out += group('chain_zone', 'chain', [wedge45, wedge23, ZONE], {'fuse': 1, 'cut': 2, 'common': 1},
                 first='fuse')
    out += group('quadrants', 'pieces', [imported('incidence_octant'), imported('incidence_back')],
                 {'fuse': None, 'cut': None, 'common': None}, reason=TOUCH)
    out += group('half_wedge', 'pieces', [imported('incidence_half'), wedge45],
                 {'fuse': None, 'cut': None, 'common': None}, reason=INSIDE)
    out += group('near_turned', 'pieces', [wedge45, imported('incidence_turned45')],
                 {'fuse': None, 'cut': None, 'common': None}, reason=NEAR)
    out += group('near_lifted', 'pieces', [imported('incidence_hemi'), imported('incidence_lifted45')],
                 {'fuse': None, 'cut': None, 'common': None}, reason=NEAR)
    return out


def all_cases():
    return cases()


def validate(listed):
    names = [c.name for c in listed]
    assert len(names) == len(set(names)), 'duplicate case names'
    for c in listed:
        assert len(c.items) == 2+len(c.stages), c.name
        assert c.last in OPS, c.name
        assert any(s[0] == 'imported' for s in c.items) or c.klass == 'sphere', c.name
        if c.klass == 'pieces':
            assert all(s[0] == 'imported' for s in c.items), c.name
        if c.klass == 'chain':
            assert c.stages, c.name
        if c.kind == 'degenerate':
            assert c.reason, c.name


# ------------------------------------------------------------------ the reference

def M(x):
    """An exact binary64 or rational as mpf."""
    q = F(x)
    return mp.mpf(q.numerator)/q.denominator


class Axial:
    """An input as the reference takes it: the ball cut to heights `[lo,
    hi]` along the axis (relative to the centre) and to the half-turns of
    its half-planes through the axis (each by its outward normal's angle
    `phi`: inside where `cos(theta - phi) <= 0`); `inside(X)` its rows' own
    inequalities (binary64, the Monte-Carlo check's)."""

    def __init__(self, axis, lo, hi, phis, inside):
        self.axis, self.lo, self.hi, self.phis, self.inside = axis, lo, hi, phis, inside

    def holds(self, h, t):
        return self.lo < h < self.hi and all(mp.cos(t-phi) < 0 for phi in self.phis)


def unit_axes(axis):
    e, u, v = AXES[axis]
    return [tuple(mp.mpf(c) for c in w) for w in (e, u, v)]


def angle_of(axis, m):
    """The angle of a direction across the axis from its `u` towards `v`."""
    _, u, v = unit_axes(axis)
    return mp.atan2(sum(a*b for a, b in zip(m, v)), sum(a*b for a, b in zip(m, u)))


def body_input(b):
    """A body's axial description, read off its rows: the box's frame (`e`
    its normal, `x` across it normalized, `y = e x x`), each of its faces
    either a half-plane holding the axis (its plane through the centre's
    axis exactly, in the rows' rationals), a parallel's plane, or apart from
    the ball."""
    e = [F(t) for t in AXES[b.axis][0]]
    x = [F(t) for t in b.x]
    assert sum(x[i]*e[i] for i in range(3)) == 0, b.name
    y = [e[1]*x[2]-e[2]*x[1], e[2]*x[0]-e[0]*x[2], e[0]*x[1]-e[1]*x[0]]
    rel = [F(CENTRE[i])-F(b.o[i]) for i in range(3)]
    r = M(RADIUS)
    E = [mp.mpf(t) for t in AXES[b.axis][0]]
    X = [M(t) for t in x]
    n = mp.sqrt(sum(t*t for t in X))
    X = [t/n for t in X]
    Y = [E[1]*X[2]-E[2]*X[1], E[2]*X[0]-E[0]*X[2], E[0]*X[1]-E[1]*X[0]]
    phis = []
    for w_rat, w, d in ((x, X, b.d[0]), (y, Y, b.d[1])):
        # The near face `w . (P - O) = 0` (outward normal -w) and the far
        # one at `d`: the centre at `w . (C - O)` between them.
        s = sum(w_rat[i]*rel[i] for i in range(3))
        at = sum(w[i]*M(rel[i]) for i in range(3))
        if s == 0:
            phis.append(angle_of(b.axis, [-t for t in w]))
        else:
            assert at > r, (b.name, 'a near face across the ball off the axis')
        assert M(d)-at > r, (b.name, 'a far face across the ball')
    h0 = -sum(rel[i]*e[i] for i in range(3))
    assert M(h0+F(b.d[2])) > r, (b.name, 'a top across the ball')
    lo = max(M(h0), -r)
    O = [float(t) for t in b.o]
    D = [float(t) for t in b.d]
    Xf, Yf, Ef = [float(t) for t in X], [float(t) for t in Y], [float(t) for t in E]

    def inside(P):
        if sum((P[i]-CENTRE[i])**2 for i in range(3)) >= RADIUS**2:
            return False
        q = [P[i]-O[i] for i in range(3)]
        coords = [sum(q[i]*w[i] for i in range(3)) for w in (Xf, Yf, Ef)]
        return all(0 < coords[k] < D[k] for k in range(3))
    return Axial(b.axis, lo, r, phis, inside)


def sphere_input(spec):
    """S9d.1's sphere on the world's axes at the centre between two
    latitudes."""
    _, r, frame, low, high = spec
    assert frame == base.at('XY', CENTRE) and r == RADIUS, spec
    R = M(r)
    lo = -R if low == -HP else R*mp.sin(M(low))
    hi = R if high == HP else R*mp.sin(M(high))
    flo, fhi = float(lo), float(hi)

    def inside(P):
        h = P[2]-CENTRE[2]
        return sum((P[i]-CENTRE[i])**2 for i in range(3)) < RADIUS**2 and flo < h < fhi
    return Axial('z', lo, hi, [], inside)


def make_input(item):
    if item[0] == 'imported':
        return body_input(BODIES[item[1]])
    return sphere_input(item)


SET_OPS = {
    'fuse': lambda a, b: a or b,
    'cut': lambda a, b: a and not b,
    'common': lambda a, b: a and b,
}


def evaluate(expr, mem):
    if isinstance(expr, int):
        return mem[expr]
    op, x, y = expr
    return SET_OPS[op](evaluate(x, mem), evaluate(y, mem))


def leaves(expr):
    if isinstance(expr, int):
        return {expr}
    return leaves(expr[1]) | leaves(expr[2])


def merged(points, period=None):
    """Sorted points with those within `MERGE` of another one (on a circle
    of `period`, across its end too)."""
    out = []
    for p in sorted(points):
        if out and p-out[-1] <= MERGE:
            continue
        out.append(p)
    if period is not None and len(out) > 1 and out[0]+period-out[-1] <= MERGE:
        out.pop()
    return out


def G(r, h):
    """The volume integrand's antiderivative: `int (r^2 - h^2) dh`."""
    return r*r*h-h**3/3


def K(r, h):
    """`int (r^2 - h^2)^(3/2) / 3 dh`."""
    s = mp.sqrt(max(r*r-h*h, mp.mpf(0)))
    return ((h/8)*(5*r*r-2*h*h)*s+(3*r**4/8)*mp.asin(h/r))/3


def L(r, h):
    """`int sqrt(r^2 - h^2) dh`: a half-plane's area between heights."""
    s = mp.sqrt(max(r*r-h*h, mp.mpf(0)))
    return (h*s+r*r*mp.asin(h/r))/2


class Cells:
    """The ball cut by every input's heights and half-planes, each cell
    kept by `expr` at its middle: the result's measures, faces and
    components."""

    def __init__(self, inputs, expr):
        axes = {I.axis for I in inputs}
        assert len(axes) == 1, 'one axis'
        self.axis = axes.pop()
        r = self.r = M(RADIUS)
        used = sorted(leaves(expr))
        hs = [-r, r]+[x for k in used for x in (inputs[k].lo, inputs[k].hi) if -r < x < r]
        self.hs = merged(hs)
        ts = []
        for k in used:
            for phi in inputs[k].phis:
                for t in (phi+mp.pi/2, phi+3*mp.pi/2):
                    ts.append(t % TWO_PI)
        ts = merged(ts, TWO_PI)
        if len(ts) < 2:
            # A whole turn (or one half-plane's two ends, never alone here).
            assert not ts, 'a lone half-plane'
            self.arcs = [(mp.mpf(0), TWO_PI)]
        else:
            self.arcs = [(ts[k], ts[k+1]) for k in range(len(ts)-1)]+[(ts[-1], ts[0]+TWO_PI)]
        self.kept = {}
        for i in range(len(self.hs)-1):
            hm = (self.hs[i]+self.hs[i+1])/2
            for j, (t0, t1) in enumerate(self.arcs):
                tm = (t0+t1)/2
                mem = [inputs[k].holds(hm, tm) if k in used else False for k in range(len(inputs))]
                self.kept[(i, j)] = evaluate(expr, mem)

    def measures(self):
        """(volume, moments about the world's origin, area, volume by the
        divergence theorem)."""
        r, hs, arcs = self.r, self.hs, self.arcs
        e, u, v = unit_axes(self.axis)
        c = [M(x) for x in CENTRE]
        V, Mu, Mv, Me, A_sph = (mp.mpf(0),)*5
        A_flat, div = mp.mpf(0), mp.mpf(0)
        whole = len(arcs) == 1
        for (i, j), k in self.kept.items():
            if not k:
                continue
            h0, h1 = hs[i], hs[i+1]
            t0, t1 = arcs[j]
            dt = t1-t0
            V += dt/2*(G(r, h1)-G(r, h0))
            Mu += (mp.sin(t1)-mp.sin(t0))*(K(r, h1)-K(r, h0))
            Mv += (mp.cos(t0)-mp.cos(t1))*(K(r, h1)-K(r, h0))
            Me += dt/2*(r*r*(h1**2-h0**2)/2-(h1**4-h0**4)/4)
            A_sph += r*dt*(h1-h0)
        # Faces between kept and dropped cells: parallels' sectors and
        # half-planes' pieces.
        for i in range(1, len(hs)-1):
            for j, (t0, t1) in enumerate(arcs):
                below, above = self.kept[(i-1, j)], self.kept[(i, j)]
                if below != above:
                    a = (t1-t0)/2*(r*r-hs[i]**2)
                    A_flat += a
                    div += (hs[i] if below else -hs[i])*a
        if not whole:
            for i in range(len(hs)-1):
                for j in range(len(arcs)):
                    if self.kept[(i, j)] != self.kept[(i, (j+1) % len(arcs))]:
                        A_flat += L(r, hs[i+1])-L(r, hs[i])
        div = (r*A_sph+div)/3
        moments = [c[q]*V+Mu*u[q]+Mv*v[q]+Me*e[q] for q in range(3)]
        return V, moments, A_sph+A_flat, div

    def components(self):
        """(the kept cells' components by faces of positive area, whether
        kept cells meet only along the axis or at a point of it)."""
        parent = {key: key for key, k in self.kept.items() if k}

        def find(a):
            while parent[a] != a:
                parent[a] = parent[parent[a]]
                a = parent[a]
            return a

        def join(a, b):
            if a in parent and b in parent:
                ra, rb = find(a), find(b)
                if ra != rb:
                    parent[ra] = rb
        n, m = len(self.hs)-1, len(self.arcs)
        for i in range(n):
            for j in range(m):
                if i+1 < n:
                    join((i, j), (i+1, j))
                if m > 1:
                    join((i, j), (i, (j+1) % m))
        count = len({find(a) for a in parent})
        touching = False
        for i in range(n):
            row = [self.kept[(i, j)] for j in range(m)]
            if m > 1 and sum(1 for j in range(m) if row[j] and not row[j-1]) > 1:
                touching = True
        for i in range(n-1):
            a = [self.kept[(i, j)] for j in range(m)]
            b = [self.kept[(i+1, j)] for j in range(m)]
            if any(a) and any(b) and not any(x and y for x, y in zip(a, b)):
                touching = True
        return count, touching

    def cell_of(self, P):
        """The kept flag of the cell holding a world point (binary64)."""
        e, u, v = unit_axes(self.axis)
        d = [P[q]-CENTRE[q] for q in range(3)]
        if sum(t*t for t in d) >= RADIUS**2:
            return False
        h = sum(d[q]*float(e[q]) for q in range(3))
        t = mp.atan2(sum(d[q]*float(v[q]) for q in range(3)), sum(d[q]*float(u[q]) for q in range(3))) % TWO_PI
        i = max(k for k in range(len(self.hs)-1) if self.hs[k] <= h)
        if len(self.arcs) == 1:
            return self.kept[(i, 0)]
        for j, (t0, t1) in enumerate(self.arcs):
            if t0 <= t < t1 or t0 <= t+TWO_PI < t1:
                return self.kept[(i, j)]
        raise AssertionError('a point in no cell')


def result(inputs, expr, solids):
    """`(N, volume, area, centre, checks)` of an expression."""
    cells = Cells(inputs, expr)
    V, moments, A, div = cells.measures()
    count, touching = cells.components()
    size = M(RADIUS)
    checks = {'divergence': abs(V-div)/size**3}
    n = 0 if V == 0 else count
    centre = [m/V for m in moments] if V > 0 else None
    return n, V, A, centre, moments, touching, checks, cells


def rows(n, V, A, centre):
    if n == 0:
        return ['empty']
    return ['result {} {} {} {} {} {}'.format(n, number(V), number(A), *(number(x) for x in centre))]


def monte_carlo(inputs, exprs, cells, samples, seed):
    """Monte-Carlo volumes and centres over the ball's box, every point's
    membership by the rows against its cell's: {op: (volume, centre,
    errors)}, and the points whose two memberships differ."""
    rnd = random.Random(seed)
    lo = [c-RADIUS for c in CENTRE]
    box = (2*RADIUS)**3
    hits = {op: [] for op in exprs}
    differ = 0
    for _ in range(samples):
        P = [lo[q]+2*RADIUS*rnd.random() for q in range(3)]
        mem = [I.inside(P) for I in inputs]
        for op, expr in exprs.items():
            inside = evaluate(expr, mem)
            if inside != cells[op].cell_of(P):
                differ += 1
            if inside:
                hits[op].append(P)
    out = {}
    for op, pts in hits.items():
        p = len(pts)/samples
        vol = box*p
        # At least one point's weight: a sliver no point fell in.
        sv = box*mp.sqrt(max(p*(1-p), mp.mpf(1)/samples)/samples)
        if pts:
            cent = [mp.fsum(P[q] for P in pts)/len(pts) for q in range(3)]
            sc = [mp.sqrt(mp.fsum((P[q]-cent[q])**2 for P in pts)/len(pts)/len(pts)) for q in range(3)]
        else:
            cent, sc = None, None
        out[op] = (vol, cent, sv, sc)
    return out, differ


def evaluate_group(job):
    """One group: rows, results and checks."""
    name, items, op1, chained, ops, degenerate, mc_n = job
    inputs = [make_input(it) for it in items]

    def expr_of(op):
        return (op, (op1, 0, 1), 2) if chained else (op, 0, 1)
    exprs = {op: expr_of(op) for op in ops}
    rows_of, res, checks, cells = {}, {}, {}, {}
    size = M(RADIUS)
    for op, expr in exprs.items():
        n, V, A, centre, moments, touching, ch, cl = result(inputs, expr, ops[op])
        rows_of[op] = rows(n, V, A, centre)
        res[op] = (n, V, A, centre, moments, touching)
        cells[op] = cl
        for k, x in ch.items():
            checks[k] = max(checks.get(k, mp.mpf(0)), x)
        if not degenerate[op]:
            assert not touching, f'{name} {op}: touching along the axis, not declared'
            assert n == ops[op], f'{name} {op}: {n} solids, {ops[op]} declared'
    # Each input alone and the pair identities.
    if not chained and set(exprs) == set(OPS):
        a = Cells(inputs, 0).measures()
        b = Cells(inputs, 1).measures()
        vf, mf = res['fuse'][1], res['fuse'][4]
        vn, mn = res['common'][1], res['common'][4]
        vt, mt = res['cut'][1], res['cut'][4]
        dev = abs(vf+vn-a[0]-b[0])/size**3
        for q in range(3):
            dev = max(dev, abs(mf[q]+mn[q]-a[1][q]-b[1][q])/size**4, abs(mt[q]-(a[1][q]-mn[q]))/size**4)
        dev = max(dev, abs(vt-(a[0]-vn))/size**3)
        checks['pair_identities'] = dev
    mc, differ = monte_carlo(inputs, exprs, cells, mc_n, zlib.crc32(name.encode()))
    z = 0.0
    for op in exprs:
        vol, cent, sv, sc = mc[op]
        _, V, _, centre, _, _ = res[op]
        z = max(z, abs(vol-V)/max(sv, mp.mpf(10)**-12))
        if cent is not None and V > 0:
            z = max(z, max(abs(cent[q]-centre[q])/max(sc[q], mp.mpf(10)**-12) for q in range(3)))
    checks['monte_carlo_sigma'] = mp.mpf(z)
    checks['cell_memberships_differing'] = mp.mpf(differ)
    return name, rows_of, checks


def jobs(mc_n):
    """One job per group."""
    groups = {}
    for c in all_cases():
        e = groups.setdefault(c.group, [c.items, c.op1 if c.stages else None, bool(c.stages), {}, {}])
        e[3][c.last] = c.solids
        e[4][c.last] = c.kind == 'degenerate'
    return [(name, items, op1, chained, ops, degenerate, mc_n)
            for name, (items, op1, chained, ops, degenerate) in groups.items()]


LIMITS = {'divergence': 1e-30, 'pair_identities': 1e-30, 'monte_carlo_sigma': 5,
          'cell_memberships_differing': 0}


def reference_checks(results):
    worst, covered = {}, {}
    for _, _, checks in results:
        for key, value in checks.items():
            worst[key] = max(worst.get(key, mp.mpf(0)), value)
            covered[key] = covered.get(key, 0)+1
    return worst, covered


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
    out = [f'# case\trow ({STEP}; generate_one_sphere_incidence_boolean_fixtures.py: expect KIND {STEP} CLASS, '
           'reason TEXT for a degenerate case, then result N volume area cx cy cz or empty)']
    for case in all_cases():
        blocks.append(case.encode())
        row = by_group[case.group][1][case.last]
        out.append(f'{case.name}\texpect {case.kind} {STEP} {case.klass}')
        if case.kind == 'degenerate':
            out.append(f'{case.name}\treason {case.reason}')
        for r in row:
            out.append(f'{case.name}\t{r}')
    frames = ['# case\tframe (solid k) axis (n, x, y)\tstored unit vector (stored_axes, as hex bits; the solids '
              'the kernel builds from rows, not the imported ones)']
    for case in all_cases():
        for k, (item, c) in enumerate(zip(case.items, case.solid_cases)):
            if item[0] == 'imported':
                continue
            _, x, y, n = stored_axes(c.frame)
            for key, w in (('n', n), ('x', x), ('y', y)):
                frames.append(f'{case.name}\tsolid{k} {key}\t'+' '.join(struct.pack('>d', q).hex() for q in w))
    return {f'{PREFIX}-cases.txt': '\n'.join(blocks)+'\n',
            f'{PREFIX}-expected.tsv': '\n'.join(out)+'\n',
            f'{PREFIX}-frames.tsv': '\n'.join(frames)+'\n',
            BODIES_FILE: bodies_text()}


def native_input():
    """Every case's native rows (`compare_one_sphere_incidence_boolean.py`)."""
    return '\n'.join(c.native() for c in all_cases())+'\n'


def case_size(case):
    """The case's size: the ball's reach from the origin."""
    return max(abs(c) for c in CENTRE)+RADIUS


# ------------------------------------------------------------------ the files' check

def stored_planes(text):
    """Every face's stored plane `(point, x, y)` (rationals of its binary64
    record), its face without a location."""
    words = text.split()
    at = words.index('Surfaces')
    count = int(words[at+1])
    k = at+2
    surfaces = []
    for _ in range(count):
        kind = int(words[k])
        n = {1: 12, 4: 13}[kind]
        vals = [F(float(w)) for w in words[k+1:k+1+n]]
        surfaces.append((kind, vals))
        k += 1+n
    _, _, shapes, _ = bio.read(text)
    out = []
    for kind, data, _ in shapes:
        if kind != 'Fa':
            continue
        s, loc = data
        assert loc == 0, 'a face with a location'
        sk, vals = surfaces[s-1]
        if sk == 1:
            out.append((vals[0:3], vals[6:9], vals[9:12]))
    return out


def check_bodies():
    """Each existing body file: its faces one sphere and planes, each
    stored vertex on the construction's surfaces within 1e-12 of the size,
    every stored plane (its normal the cross product of its stored axes)
    holding the axis through the centre, normal to it, or apart from the
    ball, exactly; a divided body's rim two arcs (two stored vertices on its
    base's circle). Returns the largest deviation and the files read."""
    worst, read = mp.mpf(0), 0
    c = [F(x) for x in CENTRE]
    for b in BODIES.values():
        path = ROOT/'fixtures'/b.path
        if not path.exists():
            continue
        read += 1
        text = path.read_text()
        kinds, vertices = base.stored_records(text)
        assert sorted(k for k in kinds if k != 'plane') == ['sphere'], (b.name, kinds)
        e = [F(t) for t in AXES[b.axis][0]]
        for p, x, y in stored_planes(text):
            n = [x[1]*y[2]-x[2]*y[1], x[2]*y[0]-x[0]*y[2], x[0]*y[1]-x[1]*y[0]]
            rel = [p[i]-c[i] for i in range(3)]
            across = [n[1]*e[2]-n[2]*e[1], n[2]*e[0]-n[0]*e[2], n[0]*e[1]-n[1]*e[0]]
            holds = sum(n[i]*e[i] for i in range(3)) == 0 and sum(n[i]*rel[i] for i in range(3)) == 0
            normal = all(t == 0 for t in across)
            off = sum(n[i]*rel[i] for i in range(3))
            apart = off*off > F(RADIUS)**2*sum(t*t for t in n)
            assert holds or normal or apart, (b.name, 'a stored plane off the axis', p, n)
        inp = body_input(b)
        size = max([1.0]+[abs(t) for v in vertices for t in v])
        for v in vertices:
            X = [mp.mpf(t) for t in v]
            d_sph = abs(mp.sqrt(sum((X[i]-c[i])**2 for i in range(3)))-M(RADIUS))
            E, U, V = unit_axes(b.axis)
            rel = [X[i]-M(c[i]) for i in range(3)]
            h = sum(rel[i]*E[i] for i in range(3))
            t = mp.atan2(sum(rel[i]*V[i] for i in range(3)), sum(rel[i]*U[i] for i in range(3)))
            rho = mp.sqrt(max(mp.mpf(0), sum(q*q for q in rel)-h*h))
            ds = [d_sph, abs(h-inp.lo)]
            for phi in inp.phis:
                ds.append(abs(rho*mp.cos(t-phi)))
            worst = max(worst, min(ds)/size)
        if b.divide:
            rim = [v for v in vertices
                   if abs(sum((v[i]-CENTRE[i])*AXES[b.axis][0][i] for i in range(3))-float(inp.lo)) <= 1e-9
                   and abs(sum((v[i]-CENTRE[i])**2 for i in range(3))-RADIUS**2) <= 1e-9]
            assert len(rim) == 2, (b.name, 'a divided rim of two arcs', len(rim))
    assert worst <= 1e-12, ('a stored vertex off its construction', mp.nstr(worst, 3))
    return worst, read


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--workers', type=int, default=max(1, min(4, (os.cpu_count() or 2)-2)))
    parser.add_argument('--samples', type=int, default=100000, help='Monte-Carlo points per group')
    parser.add_argument('--only', help='one group, printed (no files)')
    args = parser.parse_args()
    listed = all_cases()
    validate(listed)
    todo = jobs(args.samples)
    if args.only:
        for r in base.run([j for j in todo if j[0] == args.only], evaluate_group, 1):
            print(r)
        return
    results = base.run(todo, evaluate_group, args.workers)
    worst, covered = reference_checks(results)
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
    print('reference checks (largest deviation, relative to the case size):',
          ', '.join(f'{k} {mp.nstr(v, 3)}' for k, v in sorted(worst.items())))
    print('cases per check:', ', '.join(f'{k} {v}' for k, v in sorted(covered.items())))
    print(f'bodies: {read} of {len(BODIES)} files read, every stored vertex within {mp.nstr(deviation, 3)} of the '
          'size of its construction\'s surfaces, every stored plane on or across the axis exactly')


if __name__ == '__main__':
    main()
