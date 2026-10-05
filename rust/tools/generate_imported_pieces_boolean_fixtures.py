#!/usr/bin/env python3
"""Fixtures for S9e.4b.3a of REVIEW_NOTES.md: imported plane pieces of a
sphere, a cylinder or a cone (bodies of one such face and plane faces OCCT
wrote to `.brep` files) given to Booleans, decided as their primitive common
the half-spaces of their planes.

The bodies are OCCT's own output, as S9e.4a's
(`generate_imported_boolean_fixtures.py`): `boolean-imported-pieces-
bodies.txt` holds one `write NAME imported/NAME.brep` block per body for
`occt_boolean_oracle.cpp` (a primitive's rows, a `boolean` row and a box's
rows: the Boolean's one solid written), which
`compare_imported_pieces_boolean.py --write-bodies` runs to write
`rust/fixtures/imported/NAME.brep` (format version 1, no triangulations).
Nothing here reads the kernel; the files are read only to check that each
is the body it claims to be and that it is this step's (below).

`boolean-imported-pieces-cases.txt` lists each case in the Boolean protocol,
an imported input as its one `brep imported/NAME.brep` row, a chain's
further `then` rows as S9e.1's; `boolean-imported-pieces-expected.tsv`
gives per case `expect KIND S9e.4b.3a CLASS` (`solid`, `empty`,
`degenerate`: a `Degenerate` of S9's rules; `unsupported`: `OutOfDomain`, a
later sub-step's; the class `sphere`: an imported piece and a
construction, `both`: two imported inputs, `chain`: an imported
input's result given to another Boolean), a degenerate or unsupported
case's `reason TEXT`, then `result N volume area cx cy cz` (totals over the
N solids, world coordinates) or `empty`; `boolean-imported-pieces-
frames.tsv` the stored axes of every solid the kernel builds from its rows.

The bodies are spheres' pieces only: OCCT writes a sphere's section by a
plane other than a meridian or a parallel of its stored frame with a
B-spline pcurve the converter cannot certify on the sphere
(`UncertifiedPcurveOffEdge`), and a cylinder's or a cone's oblique section
as an ellipse record the `.brep` reader does not read (`Ellipse`, the
import track's), so their pieces' planes are each the sphere's meridian
planes or parallels' planes, as the DRAW survey's `so1` to `so7` are. The
frames are rational rotations none of whose axes is a world axis or normal
to one (the reference's sphere families are the world's meridians and
parallels, which a plane holding the world's `z` through the centre, or
normal to it, would carry whole) nor one of the kernel's whole sphere's
own axes. Each body is `P common B` (`bitten` `P cut B`), `P` a whole
sphere on a frame and `B` a box or prism on the same frame's axes (its
planes the sphere's meridian planes and parallels' planes): `octant`, the
ball of radius 5 about `(5, 5, 4)` on `SKEW` common the box of side 8 at
its centre (the corner of its axes' positive combinations, three planes
through the centre, the stored pole where two of them meet the sphere);
`octant_tilt`, the same on `SKEW4` about `(21/4, 17/2, 7/4)`; `upper`, the ball of
radius 5 about `(3, 4, 1)` on `SKEW2` common the box of side 8 at the point
of its axis 2 above its centre (the wedge above a parallel's plane between
two meridian planes, `so5`'s); `lune`, the ball of radius 4 about `(5, 5,
5)` on `SKEW` between two meridian planes; `half`, the ball of radius 4
about `(5, 5, 5)` on one side of a meridian plane of `SKEW2` (its sphere's
frame's `x` reversed, its seam outside it: OCCT splits a sphere face its
seam crosses); `zone_wedge`, the ball of radius 5 about `(5, 5, 4)` on
`SKEW3` between two meridian planes and two parallels' planes at `-3/2`
and `5/2` along its axis (not symmetric: the reference's events where a
wall's generatrix leaves the sphere at both ends would coincide); declared S9e.4b.3c's: `octant_low`, the ball of `octant` below a
parallel's plane (S9e.4a's cap, a piece of the same sphere), and `bitten`,
the ball of `half` less the box of side 8 at its centre (not convex in its
planes).

Cases (each the three operations): `octant_box`, `box_octant`, the octant
and a box across its sphere face, either way; `tilt_rod`, the turned octant
and an upright rod through its sphere face; `upper_slab`, the upper wedge
and a `TILT` slab; `lune_ball`, `lune_cone`, the lune and a ball, an upright
cone across its sphere face; `half_rod`, the half ball and a rod across its
meridian disc; `zone_box`, the zone's wedge and a box; `pieces`, the
turned octant and the zone's wedge, both imported (two spheres);
`half_box`, the half ball and S9e.4a's imported box (its top across it);
`chain_octant`, the octant less an upright rod, then with a `TILT` slab.
Declared `degenerate`: `tilt_flush`, a box in the turned octant's frame on
its base plane (one plane in the construction, within the resolution of it
in the file); `half_touch`, a box whose face touches the half ball's
sphere at a point inside its face. `one_sphere`, the octant and
`octant_low` (faces of both on one sphere), declared `unsupported` until
S9e.4b.3c.1's kernel decided them, is solid; so is `bitten_box`, the bitten
ball and a box, declared `unsupported` until S9e.4b.3c.3a's kernel decided
the bitten ball as a bite (its ball less the hull of its planes turned
over).

The reference is the construction OCCT was given, through S9e.3a's chained
reference (`chained_curved_boolean_reference.py`): each imported piece its
first Boolean, `(P common B) op C` (swapped, the piece the tool), two
pieces `(P common B) op (Q common D)`, a chain `((P common B) op1 C) op2
D`, with S9e.4a's checks (the two families, each solid's closed form, the
pair identities on the last Boolean's arguments, Monte Carlo, solid counts
by rays at two resolutions, each piece one solid, the meetings' sines and
the events' spacing outside the declared cases) and each piece's closed
form where it has one (`corner_closed`: the corner's solid angle by Van
Oosterom and Strackee; `lune_closed`: the dihedral angle; the half). With
`--check`, where the bodies' files exist, each is read (`stored_records`):
its faces' kinds one sphere and planes, every stored vertex within 1e-12
of the size on the construction's surfaces; and each construction is no
S9e.4a construction (a plane of its box not normal to its sphere's axis):
the reason each body is this step's.
"""
import argparse
from fractions import Fraction as F
import os
import struct
import zlib

import mpmath as mp

from identity_reference import Case, native_case
from curve_surface_reference import stored_axes
import generate_imported_boolean_fixtures as base
from generate_curved_boolean_fixtures import square, disc

ROOT = base.ROOT
OPS = base.OPS
OPERATIONS, BOOLEANS = base.OPERATIONS, base.BOOLEANS
CLASSES = ('sphere', 'both', 'chain')
BODIES_FILE = 'boolean-imported-pieces-bodies.txt'
STEP = 'S9e.4b.3a'
# Frames of rational rotations, none of whose axes is a world axis or
# normal to one: the reference sweeps a sphere by meridians about the
# world's `z` and parallels normal to it, which a plane holding the world's
# `z` through the centre, or normal to it, would carry whole.
# Nor any of the kernel's whole sphere's own axes (S9d.1's rational
# rotations, `(2, 3, 6)`, `(1, 4, 8)`, `(2, 6, 9)`, `(1, 2, 2)`), whose split
# would run through the planes' line. Each normalizes to the same stored
# axes with a correctly rounded `hypot` and with macOS's (a partner on
# `SKEW3`'s frame is built by the kernel too).
FRAMES = dict(base.FRAMES, SKEW=(8.0, 4.0, 1.0, -1.0, 4.0, -8.0), SKEW2=(4.0, 4.0, 7.0, 1.0, -8.0, 4.0),
              SKEW3=(10.0, 11.0, 2.0, -10.0, 10.0, -5.0), SKEW4=(4.0, 1.0, 8.0, -7.0, -4.0, 4.0))
prism, sphere, cone, construction = base.prism, base.sphere, base.cone, base.construction
ref = base.ref


def at(name, origin):
    return tuple(float(c) for c in origin)+FRAMES[name]


# ------------------------------------------------------------------ bodies

class Body:
    """An imported piece: `primitive op box`, OCCT's rows for it, its class,
    its closed form `(volume, area)` where it has one, and whether it is
    S9e.4a's construction instead (`earlier`)."""

    def __init__(self, name, klass, primitive, box, op='common', closed=None, earlier=False):
        self.name, self.klass, self.primitive, self.box, self.op = name, klass, primitive, box, op
        self.closed, self.earlier = closed, earlier
        self.path = f'imported/{name}.brep'

    def specs(self):
        return [self.primitive, self.box]

    def native_rows(self):
        rows = lambda spec: native_case(construction(spec, 0)).split('\n')[1:-1]
        return rows(self.primitive)+[f'boolean {self.op}']+rows(self.box)


def box(frame, size, low=0.0):
    """`MakeBox`'s prism: `[0, DX] x [0, DY]` on the frame over `[low, low +
    DZ]`."""
    return prism([square(0.0, 0.0, size[0], size[1])], frame, low, low+size[2])


def axes(frame):
    """A frame's stored axes `x`, `y`, `n` as unit mpf vectors."""
    _, x, y, n = stored_axes(frame)
    q = lambda c: mp.mpf(F(c).numerator)/F(c).denominator
    unit = lambda v: [q(c)/mp.sqrt(sum(q(d)**2 for d in v)) for c in v]
    return unit(x), unit(y), unit(n)


def angle(a, b):
    """The angle between two unit vectors (`atan2` of the cross product's
    length and the dot product: no rounding near 0 or pi)."""
    c = [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
    return mp.atan2(mp.sqrt(sum(v*v for v in c)), sum(p*q for p, q in zip(a, b)))


def corner_closed(r, frame):
    """A ball of radius `r` about the frame's origin common the cone of its
    stored axes' positive combinations (three planes through the centre):
    its volume `Omega r^3 / 3` (`Omega` the corner's solid angle by Van
    Oosterom and Strackee), its area `Omega r^2` on the sphere and the three
    sectors between the axes."""
    R = mp.mpf(r)
    x, y, n = axes(frame)
    triple = abs(x[0]*(y[1]*n[2]-y[2]*n[1])-x[1]*(y[0]*n[2]-y[2]*n[0])+x[2]*(y[0]*n[1]-y[1]*n[0]))
    dot = lambda a, b: sum(p*q for p, q in zip(a, b))
    omega = 2*mp.atan2(triple, 1+dot(x, y)+dot(y, n)+dot(n, x))
    return omega*R**3/3, omega*R**2+(angle(x, y)+angle(y, n)+angle(n, x))*R**2/2


def lune_closed(r, frame):
    """A ball of radius `r` about the frame's origin between its stored
    planes `u = 0` and `v = 0` (`u, v >= 0`, both holding the stored `n`):
    the dihedral angle `alpha` between `x` and `y` projected normal to `n`,
    the volume `2 alpha r^3 / 3`, the area `2 alpha r^2` on the sphere and
    two half discs."""
    R = mp.mpf(r)
    x, y, n = axes(frame)
    perp = lambda v: [v[i]-sum(p*q for p, q in zip(v, n))*n[i] for i in range(3)]
    unit = lambda v: [c/mp.sqrt(sum(d*d for d in v)) for c in v]
    alpha = angle(unit(perp(x)), unit(perp(y)))
    return 2*alpha*R**3/3, 2*alpha*R**2+mp.pi*R**2


def half_closed(r):
    """A ball of radius `r` on one side of a plane through its centre."""
    R = mp.mpf(r)
    return 2*mp.pi*R**3/3, 3*mp.pi*R**2


OCTANT = at('SKEW', (5, 5, 4))
OCTANT_TILT = at('SKEW4', (5.25, 8.5, 1.75))
LUNE = at('SKEW', (5, 5, 5))
HALF = at('SKEW2', (5, 5, 5))
# The half ball's sphere on `SKEW2` with its `x` reversed: its seam (the
# frame's `x`) outside the half `u >= 0` (OCCT splits a sphere face its seam
# crosses: two faces on one sphere, S9e.4b.3c).
HALF_SPHERE = (5.0, 5.0, 5.0, 4.0, 4.0, 7.0, -1.0, 8.0, -4.0)
ZONE = at('SKEW3', (5, 5, 4))
# The upper wedge's box at the point of its axis 2 above its centre.
UPPER_BALL = at('SKEW2', (3, 4, 1))
UPPER_BOX = tuple(float(F(UPPER_BALL[i])+2*F(stored_axes(UPPER_BALL)[3][i])) for i in range(3))+FRAMES['SKEW2']


def make_bodies():
    out = [
        Body('octant', 'sphere', sphere(5.0, OCTANT), box(OCTANT, (8.0, 8.0, 8.0)),
             closed=lambda: corner_closed(5, OCTANT)),
        Body('octant_tilt', 'sphere', sphere(5.0, OCTANT_TILT), box(OCTANT_TILT, (8.0, 8.0, 8.0)),
             closed=lambda: corner_closed(5, OCTANT_TILT)),
        Body('upper', 'sphere', sphere(5.0, UPPER_BALL), box(UPPER_BOX, (8.0, 8.0, 8.0))),
        Body('lune', 'sphere', sphere(4.0, LUNE), prism([square(0.0, 0.0, 6.0, 6.0)], LUNE, -5.0, 5.0),
             closed=lambda: lune_closed(4, LUNE)),
        Body('half', 'sphere', sphere(4.0, HALF_SPHERE), prism([square(0.0, -6.0, 6.0, 6.0)], HALF, -6.0, 6.0),
             closed=lambda: half_closed(4)),
        Body('zone_wedge', 'sphere', sphere(5.0, ZONE), prism([square(0.0, 0.0, 8.0, 8.0)], ZONE, -1.5, 2.5)),
        Body('octant_low', 'sphere', sphere(5.0, OCTANT), prism([square(-6.0, -6.0, 6.0, 6.0)], OCTANT, -6.0, 2.0),
             earlier=True),
        Body('bitten', 'sphere', sphere(4.0, HALF), box(HALF, (8.0, 8.0, 8.0)), op='cut'),
    ]
    return {b.name: b for b in out}


BODIES = make_bodies()


def imported(name):
    return ('imported', name)


def flatten(items):
    """The reference's solids of a case's inputs (an imported piece its
    primitive and box, S9e.4a's imported body its construction) and each
    input's expression over them."""
    specs, exprs = [], []
    for it in items:
        if it[0] == 'imported':
            b = BODIES[it[1]]
            k = len(specs)
            specs += b.specs()
            exprs.append((b.op, k, k+1))
        elif it[0] == 'base':
            exprs.append(len(specs))
            specs.append(base.BODIES[it[1]].spec)
        else:
            exprs.append(len(specs))
            specs.append(it)
    return specs, exprs


def spec_frame(spec):
    return spec[{'sphere': 2, 'prism': 2, 'cone': 4}[spec[0]]]


class Imported:
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
            if s[0] in ('imported', 'base'):
                path = (BODIES if s[0] == 'imported' else base.BODIES)[s[1]].path
                c = Case(name, 1e-7, OPERATIONS[k], at('XY', (0, 0, 0)), 0.0, 0.0, [], brep=path)
            else:
                c = construction(s, OPERATIONS[k])
                c.name = name
            self.solid_cases.append(c)

    @property
    def last(self):
        return self.stages[-1][0] if self.stages else self.op1

    def expr(self):
        _, e = flatten(self.items)
        if not self.stages:
            return (self.op1, e[0], e[1])
        (op2, sw2) = self.stages[0]
        first = (self.op1, e[0], e[1])
        return (op2, e[2], first) if sw2 else (op2, first, e[2])

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


def group(name, klass, items, outcomes, first=None, reason=None, kind=None):
    """Cases of one group: the last operation varied over `outcomes` (`{op:
    solids}`, 0 for empty)."""
    out = []
    for op, n in outcomes.items():
        k = kind or ('degenerate' if reason else ('solid' if n else 'empty'))
        stages = [(op, False)] if first else []
        out.append(Imported(f'{name}_{op}', klass, items, first or op, stages, k, n, reason))
    return out


TOUCH = 'a box\'s face touching the imported half ball\'s sphere at a point inside its face'
FLUSH = ('a box in the turned octant\'s frame on the plane of its base (one plane in the construction, within '
         'the resolution of it in the file)')
ROD = prism([disc(5.5, 7.375, 0.75)], at('XY', (0, 0, -2)), 0.0, 14.0)
SLAB = prism([square(-4.0, -15.0, 14.0, 15.0)], at('TILT', (0, 0, 8.125)), 0.0, 1.5)
# The half ball's touching box: on the plane normal to `SKEW2`'s `x` through
# the sphere's point there.
TOUCH_FRAME = tuple(float(F(HALF[i])+4*F(stored_axes(HALF)[1][i])) for i in range(3))+(1.0, -8.0, 4.0, 4.0, 4.0, 7.0)


def cases():
    out = []
    three = {'fuse': 1, 'cut': 1, 'common': 1}
    octant_box = prism([square(4.5, 8.5, 7.5, 11.25)], at('XY', (0, 0, 1.75)), 0.0, 3.0)
    out += group('octant_box', 'sphere', [imported('octant'), octant_box], three)
    out += group('box_octant', 'sphere', [octant_box, imported('octant')], three)
    out += group('tilt_rod', 'sphere', [imported('octant_tilt'), prism([disc(5.5, 5.25, 0.75)],
                                                                       at('XY', (0, 0, -2.5)), 0.0, 11.5)], three)
    out += group('upper_slab', 'sphere', [imported('upper'), prism([square(-4.0, -15.0, 14.0, 15.0)],
                                                                    at('TILT', (0, 0, 4.625)), 0.0, 1.5)],
                 {'fuse': 1, 'cut': 2, 'common': 1})
    out += group('lune_ball', 'sphere', [imported('lune'), sphere(1.5, at('XY', (3.25, 8.875, 3.625)))], three)
    out += group('lune_cone', 'sphere', [imported('lune'), cone(1.5, 0.5, 5.0, at('XY', (3.75, 8.25, 1.5)))], three)
    out += group('half_rod', 'sphere', [imported('half'), prism([disc(5.25, -5.75, 1.0)],
                                                                (0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0), 0.0, 12.0)], three)
    out += group('zone_box', 'sphere', [imported('zone_wedge'), prism([square(0.25, 6.25, 3.0, 9.25)],
                                                                      at('XY', (0, 0, 4.75)), 0.0, 2.5)], three)
    out += group('pieces', 'both', [imported('octant'), imported('octant_tilt')], three)
    out += group('half_box', 'both', [imported('half'), ('base', 'box')], three)
    out += group('chain_octant', 'chain', [imported('octant'), ROD, SLAB], {'fuse': 1, 'cut': 2, 'common': 1},
                 first='cut')
    out += group('tilt_flush', 'sphere', [imported('octant_tilt'),
                                          prism([square(-3.0, -3.0, 3.0, 3.0)], OCTANT_TILT, -2.0, 0.0)],
                 {'fuse': 1, 'cut': 1, 'common': 0}, reason=FLUSH)
    out += group('half_touch', 'sphere', [imported('half'), prism([square(-2.0, -2.0, 2.0, 2.0)], TOUCH_FRAME,
                                                                  0.0, 2.0)],
                 {'fuse': 1, 'cut': 1, 'common': 0}, reason=TOUCH)
    out += group('one_sphere', 'both', [imported('octant'), imported('octant_low')], three)
    out += group('bitten_box', 'sphere', [imported('bitten'), prism([square(6.0, 5.5, 9.5, 8.5)],
                                                                    at('XY', (0, 0, 3.5)), 0.0, 3.0)],
                 three)
    return out


def all_cases():
    return cases()


def validate(listed):
    names = [c.name for c in listed]
    assert len(names) == len(set(names)), 'duplicate case names'
    for c in listed:
        assert len(c.items) == 2+len(c.stages), c.name
        assert c.last in OPS, c.name
        assert any(s[0] == 'imported' for s in c.items), c.name
        if c.klass == 'both':
            assert sum(s[0] in ('imported', 'base') for s in c.items) == 2, c.name
        if c.klass == 'chain':
            assert c.stages, c.name
        specs, _ = flatten(c.items)
        for s in specs:
            assert s[0] in ('sphere', 'prism', 'cone'), c.name


# ------------------------------------------------------------------ the reference

def owners(items, sub):
    """Each reference solid's case input."""
    out = {}
    for k, e in enumerate(sub):
        for j in ref.leaves(e):
            out[j] = k
    return out


def evaluate(job):
    """One group on the chained reference: rows, results and checks."""
    name, items, op1, swapped, chained, ops, degenerate, mc_n = job
    specs, sub = flatten(items)
    chain = ref.Chain([construction(s, 0) for s in specs])
    size = chain.size
    first = (op1, sub[0], sub[1]) if chained else None

    def expr_of(op):
        if not chained:
            return (op, sub[0], sub[1])
        return (op, sub[2], first) if swapped else (op, first, sub[2])
    exprs = {op: (expr_of(op), solids) for op, solids in ops.items()}
    rows, res = {}, {}
    for op, (expr, solids) in exprs.items():
        rows[op] = ref.rows(chain, expr, solids)
        res[op] = chain.result(expr, solids)
    checks = {}
    fam = mp.mpf(0)
    for op, (expr, _) in exprs.items():
        a, b = chain.measures(expr, 0), chain.measures(expr, 1)
        fam = max(fam, abs(a[0]-b[0])/size**3, abs(a[2]-b[2])/size**2,
                  max(abs(a[1][i]-b[1][i]) for i in range(3))/size**4)
    checks['families'] = fam
    cdev = mp.mpf(0)
    for k, I in enumerate(chain.inputs):
        V, mom, A = I.closed()
        v, m, a = chain.measures(k)
        cdev = max(cdev, abs(v-V)/size**3, abs(a-A)/size**2, max(abs(m[i]-mom[i]) for i in range(3))/size**4)
    checks['solid_closed_forms'] = cdev
    # Each imported piece one solid, its closed form where it has one.
    for it, e in zip(items, sub):
        if it[0] != 'imported':
            continue
        assert ref.count_solids(chain, e, 48) == 1, f'{name}: {it[1]} is not one solid'
        closed = BODIES[it[1]].closed
        if closed is not None:
            V, A = closed()
            v, _, a = chain.measures(e)
            d = max(abs(v-V)/size**3, abs(a-A)/size**2)
            checks['piece_closed_forms'] = max(checks.get('piece_closed_forms', mp.mpf(0)), d)
    # The pair identities on the last Boolean's arguments.
    if chained:
        assert ref.count_solids(chain, first, 48) == 1, f'{name}: the given result is not one solid'
        x_expr, c_expr = (sub[2], first) if swapped else (first, sub[2])
    else:
        x_expr, c_expr = sub[0], sub[1]
    vx, mx, ax = chain.measures(x_expr)
    vc, mc, ac = chain.measures(c_expr)
    m = {op: chain.measures(e) for op, (e, _) in exprs.items()}
    if set(m) == set(OPS):
        vf, mf, af = m['fuse']
        vn, mn, an = m['common']
        vt, mt, at_ = m['cut']
        dev = abs(vf+vn-vx-vc)/size**3
        for i in range(3):
            dev = max(dev, abs(mf[i]+mn[i]-mx[i]-mc[i])/size**4)
        dev = max(dev, abs(vt-(vx-vn))/size**3, max(abs(mt[i]-(mx[i]-mn[i])) for i in range(3))/size**4)
        checks['pair_identities'] = dev
        # Faces of two inputs of the case on one surface (not a piece's own
        # primitive and box) leave the area identity.
        owner = owners(items, sub)
        shared = any(owner[i] != owner[sw.others.surfs[k][0]] for i, _, sw in chain.sweeps(0) for k in sw.skip)
        if not shared:
            checks['area_identity'] = abs(af+an-ax-ac)/size**2
    for op, (expr, solids) in exprs.items():
        for n in (48, 71):
            got = ref.count_solids(chain, expr, n)
            assert got == solids, f'{name} {op}: {got} solids by rays ({n}), {solids} declared'
    mc = ref.monte_carlo(chain, [e for e, _ in exprs.values()], mc_n, zlib.crc32(name.encode()))
    z = 0.0
    for op, (expr, _) in exprs.items():
        vol, c, sv, sc = mc[expr]
        v, mom, _ = chain.measures(expr)
        z = max(z, abs(vol-float(v))/max(sv, 1e-12))
        if c is not None and v > 0:
            cent = [float(mom[i]/v) for i in range(3)]
            z = max(z, max(abs(c[i]-cent[i])/max(sc[i], 1e-12) for i in range(3)))
    checks['monte_carlo_sigma'] = mp.mpf(z)
    sine, gap, _ = chain.margins()
    margins = {'sine': sine, 'gap': gap, 'spacing': base.spacing(chain)}
    stats = {'quadrature': chain.quad_error()/size**4, 'missed': chain.missed(), 'refined': chain.refined()}
    return name, rows, res, checks, margins, stats, degenerate


# Groups whose margins (the meetings' sines, the events' spacing) are not
# checked though solid: `bitten_box`, declared `unsupported` until
# S9e.4b.3c.3a, its box's edge from the ball's centre along the frame's
# normal through the sphere's stored pole within rounding (the reference's
# events there 6.2e-18 of their range apart; the kernel takes the pole as a
# vertex of the edge).
MARGINS_UNCHECKED = {'bitten_box'}


def jobs(mc_n):
    """One job per group."""
    groups = {}
    for c in all_cases():
        e = groups.setdefault(c.group, [c.items, c.op1 if c.stages else None, bool(c.stages and c.stages[0][1]),
                                        bool(c.stages), {},
                                        c.kind in ('degenerate', 'unsupported') or c.group in MARGINS_UNCHECKED])
        e[4][c.last] = c.solids
    return [(name, items, op1, swapped, chained, ops, degenerate, mc_n)
            for name, (items, op1, swapped, chained, ops, degenerate) in groups.items()]


def cached(job):
    """`evaluate`, its result kept in `--cache` (a development aid)."""
    cache = os.environ.get('IMPORTED_PIECES_CACHE')
    if not cache:
        return evaluate(job)
    import hashlib
    import inspect
    import pickle
    from pathlib import Path
    key = hashlib.sha256((repr(job)+repr(flatten(job[1])[0])+inspect.getsource(evaluate)
                          + inspect.getsource(corner_closed)+Path(ref.__file__).read_text()).encode()).hexdigest()[:24]
    path = Path(cache)/f'{job[0]}-{key}.pickle'
    if path.exists():
        return pickle.loads(path.read_bytes())
    try:
        out = evaluate(job)
    except Exception as e:
        raise RuntimeError(f'{job[0]}: {e!r}') from e
    path.write_bytes(pickle.dumps(out))
    return out


LIMITS = dict(base.LIMITS, piece_closed_forms=1e-30)


def reference_checks(results):
    worst, covered = {}, {}

    def note(key, value):
        worst[key] = max(worst.get(key, mp.mpf(0)), value)
        covered[key] = covered.get(key, 0)+1

    for name, rows, _, checks, margins, stats, degenerate in results:
        for key, value in checks.items():
            note(key, value)
        note('quadrature_estimate', stats['quadrature'])
        assert stats['missed'] == 0, f'{name}: quadrature nodes off their structure'
        if not degenerate:
            note('least_sine', 1/max(margins['sine'], mp.mpf(10)**-30))
            note('event_spacing', 1/max(margins['spacing'], mp.mpf(10)**-40))
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
    out = [f'# case\trow ({STEP}; chained_curved_boolean_reference.py: expect KIND {STEP} CLASS, reason TEXT for a '
           'degenerate or unsupported case, then result N volume area cx cy cz or empty)']
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
        for k, (item, c) in enumerate(zip(case.items, case.solid_cases)):
            if item[0] in ('imported', 'base'):
                continue
            _, x, y, n = stored_axes(c.frame)
            for key, v in (('n', n), ('x', x), ('y', y)):
                frames.append(f'{case.name}\tsolid{k} {key}\t'+' '.join(struct.pack('>d', q).hex() for q in v))
    return {'boolean-imported-pieces-cases.txt': '\n'.join(blocks)+'\n',
            'boolean-imported-pieces-expected.tsv': '\n'.join(out)+'\n',
            'boolean-imported-pieces-frames.tsv': '\n'.join(frames)+'\n',
            BODIES_FILE: bodies_text()}


def native_input():
    """Every case's native rows (`compare_imported_pieces_boolean.py`)."""
    return '\n'.join(c.native() for c in all_cases())+'\n'


def case_size(case):
    """The case's size: its reference solids' reach from the origin."""
    specs, _ = flatten(case.items)
    chain = ref.Chain([construction(s, 0) for s in specs])
    return float(chain.size)


# ------------------------------------------------------------------ the files' check

def axis(spec):
    """A primitive's stored axis (its frame's normal), as floats."""
    return stored_axes(spec_frame(spec))[3]


def oblique(body):
    """Whether a body is no S9e.4a construction: a plane of its box not
    normal to a sphere's axis, or oblique to a cylinder's or a cone's
    (neither normal to nor along it)."""
    a = axis(body.primitive)
    _, x, y, n = stored_axes(spec_frame(body.box))
    dot = lambda u, v: sum(p*q for p, q in zip(u, v))
    for m in (x, y, n):
        c = abs(dot(m, a))
        if body.primitive[0] == 'sphere' and c < 1-1e-9:
            return True
        if body.primitive[0] != 'sphere' and 1e-9 < c < 1-1e-9:
            return True
    return False


def check_bodies():
    """Each existing body file: its faces' kinds one sphere, cylinder or
    cone and planes, each stored vertex on the construction's surfaces within
    1e-12 of the body's size; each construction no S9e.4a construction.
    Returns the largest deviation and the files read."""
    worst, read = mp.mpf(0), 0
    for b in BODIES.values():
        assert b.earlier or oblique(b), (b.name, 'an S9e.4a construction')
        path = ROOT/'fixtures'/b.path
        if not path.exists():
            continue
        read += 1
        kinds, vertices = base.stored_records(path.read_text())
        curved = [k for k in kinds if k != 'plane']
        want = 'cylinder' if b.primitive[0] == 'prism' else b.primitive[0]
        assert curved == [want], (b.name, kinds)
        surfs = []
        for spec in b.specs():
            surfs += base.construction_surfaces(spec)[1]
        size = max([1.0]+[abs(x) for v in vertices for x in v])
        for v in vertices:
            X = tuple(mp.mpf(c) for c in v)
            d = min(abs(s.f(X))/max(base.tc.norm(s.grad(X)), mp.mpf(10)**-30) for s in surfs)
            worst = max(worst, d/size)
    assert worst <= 1e-12, ('a stored vertex off its construction', mp.nstr(worst, 3))
    return worst, read


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--workers', type=int, default=max(1, min(4, (os.cpu_count() or 2)-2)))
    parser.add_argument('--samples', type=int, default=200000, help='Monte-Carlo points per group')
    parser.add_argument('--only', help='one group, printed (no files)')
    parser.add_argument('--cache', help='a directory keeping each group\'s result (development)')
    args = parser.parse_args()
    if args.cache:
        os.environ['IMPORTED_PIECES_CACHE'] = args.cache
    listed = all_cases()
    validate(listed)
    todo = jobs(args.samples)
    if args.only:
        todo = [j for j in todo if j[0] == args.only]
        for r in base.run(todo, cached, 1):
            print(r)
        return
    results = base.run(todo, cached, args.workers)
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
    print('reference checks (largest deviation, relative to the case size; least_sine and event_spacing as '
          'reciprocals):', ', '.join(f'{k} {mp.nstr(v, 3)}' for k, v in sorted(worst.items())))
    print('cases per check:', ', '.join(f'{k} {v}' for k, v in sorted(covered.items())))
    print(f'bodies: {read} of {len(BODIES)} files read, every stored vertex within {mp.nstr(deviation, 3)} of the '
          'size of its construction\'s surfaces')


if __name__ == '__main__':
    main()
