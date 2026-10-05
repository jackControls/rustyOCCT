#!/usr/bin/env python3
"""Fixtures for S9e.4b.3b of REVIEW_NOTES.md: the kernel's own plane pieces
(S8's split pieces: a prism's oblique piece, `Clipped`, and a cone's, a
zone's or a torus's piece, `Half`) given to Booleans with curved faces,
decided as their primitive common the half-space of their plane (S9e.4b.3a's
model).

A split piece is a solid's rows and one row `split ox oy oz nx ny nz xx xy
xz below|above` (`identity_reference.Case.split`): the kernel builds the
solid, splits it by the plane through the frame's origin normal to its
normal (`Solid::split_by_plane`) and keeps the one piece on that side (along
the normal: `above`); `occt_boolean_oracle.cpp` keeps the solid common
`BRepPrimAPI_MakeHalfSpace` of the plane on that side. Nothing here reads
the kernel.

`boolean-split-pieces-cases.txt` lists each case in the Boolean protocol, a
chain's further `then` rows as S9e.1's; `boolean-split-pieces-expected.tsv`
gives per case `expect KIND S9e.4b.3b CLASS` (`solid`, `empty`,
`degenerate`: a `Degenerate` of S9's rules; `unsupported`: `OutOfDomain`, a
later sub-step's; the class `prism`: a prism's oblique piece and a
construction, `revolved`: a cone's, a zone's or a torus's piece and a
construction, `both`: two split pieces, `chain`: a split piece's result given
to another Boolean), a degenerate or unsupported case's `reason TEXT`, then
`result N volume area cx cy cz` (totals over the N solids, world
coordinates) or `empty`; `boolean-split-pieces-frames.tsv` the stored axes
of every solid the kernel builds from its rows and of every split's plane.

The pieces: `cyl_low`, a cylinder of radius 3 below a plane leaning across
its wall (an ellipse, the cut through the axis at `z = 3`); `dee_low`, a
prism of lines and an arc (S9e.4a's `dee`) below a tilted plane crossing
both caps; `block_high`, a box above a leaning plane crossing both caps
(plane faces only, against a curved partner); `frustum_low`, a frustum below
a tilted plane across its wall (an ellipse); `zone_half`, a zone on `SKEW2`
on one side of the plane through its axis normal to its frame's `x`
(S8c.2's half); `zone_cut`, a zone on `SKEW` above a tilted plane (a circle
across its wall, S8d.2); `band`, a whole torus on `XY` above a plane
normal to its axis (S8d.1's band). Declared `degenerate`: `cone_axis`, a
frustum's half by a plane through its axis (through its virtual apex, S9d.3a's
rule) against a box; `cyl_flush`, a box on the cylinder piece's cut plane
(one plane in the construction, within the resolution of it in the kernel's
two models). `one_sphere`, the zone's half against the whole ball of its
sphere, declared `unsupported` until S9e.4b.3c.1's kernel decided it, is
solid (its cut empty).

Cases (each the three operations): `cyl_box`, `box_cyl`, the cylinder's
piece and a box across its wall and cut face, either way; `cyl_ball` and
`cyl_rod`, against a ball and a rod along the world's `y` across its cut
face; `dee_rod`, the dee's piece and an upright rod across its arc wall;
`block_rod`, the box's piece (plane faces only) and a rod along `y`;
`frustum_ball`, the frustum's piece and a ball across its section;
`zone_box` and `zone_cone`, the zone's half and a box, and a cone along the
zone's axis inside its band, across its cut face and wall; `cut_rod`, the
`SKEW` zone's piece and a rod along `y`; `band_ball`, the band and a ball
biting its top; `pair`, the cylinder's and the dee's pieces (two split
pieces); `chain_cyl`, the cylinder's piece less a rod along `y`, then with a
level slab across it.

The reference is the construction OCCT is given, through S9e.3a's chained
reference (`chained_curved_boolean_reference.py`): each split piece its
first Boolean `P common H`, `H` a box on the plane's frame on the kept side
reaching past the solid (its face at the frame's origin the plane, its
other faces apart from the solid), `(P common H) op C`, swapped, two pieces
`(P common H) op (Q common K)`, a chain `((P common H) op1 C) op2 D` (a
zone's piece its whole sphere common `H` common the slab between its ends'
parallels, at the heights the kernel stores, `r sin(lat)`: the chained
reference takes whole spheres), with
S9e.4a's checks (the two families, each solid's closed form, the pair
identities on the last Boolean's arguments, Monte Carlo, solid counts by
rays at two resolutions, each piece one solid, the meetings' sines and the
events' spacing outside the declared cases) and the cylinder's piece's
closed form (`cylinder_cut_closed`: its volume the disc's area times the
plane's height over the axis, its area the disc, the wall's mean height and
the ellipse's `pi r^2 |m| / |m_z|`, `m` its box's face's normal). Each
split is checked to be the one the kernel takes as a piece of its class: a
prism's plane oblique to its axis (S8a.2), a revolved solid's plane through
its axis within a quarter of the resolution (S8c.2) or neither (S8d.2), a
torus's normal to its axis (S8d.1).
"""
import argparse
from fractions import Fraction as F
import os
import struct
import zlib

import mpmath as mp

from brep_reference import sin_rn
from curve_surface_reference import stored_axes
import generate_imported_boolean_fixtures as base
import generate_imported_pieces_boolean_fixtures as pieces
from generate_curved_boolean_fixtures import square, disc

ROOT = base.ROOT
OPS = base.OPS
OPERATIONS, BOOLEANS = base.OPERATIONS, base.BOOLEANS
CLASSES = ('prism', 'revolved', 'both', 'chain')
STEP = 'S9e.4b.3b'
FRAMES = pieces.FRAMES
prism, sphere, cone, torus, construction = base.prism, base.sphere, base.cone, base.torus, base.construction
ref = base.ref
HP = base.HP


def at(name, origin):
    return tuple(float(c) for c in origin)+FRAMES[name]


# ------------------------------------------------------------------ pieces

class Piece:
    """A split piece: a solid split by a plane (a frame: its origin and
    normal), the side kept, its class (`clipped`, `half`, `conic`, `band`),
    its closed form `(volume, area)` where it has one."""

    def __init__(self, name, kind, solid, plane, side, reach, closed=None):
        assert side in ('below', 'above')
        self.name, self.kind, self.solid, self.plane, self.side = name, kind, solid, plane, side
        self.reach, self.closed = reach, closed

    def halfspace(self):
        """The box on the plane's frame on the kept side, reaching `reach`
        from the plane's origin each way (past the solid: its farthest point
        lies within four fifths of it, `test_halfspace_boxes`)."""
        r = self.reach
        low, high = (0.0, r) if self.side == 'above' else (-r, 0.0)
        return prism([square(-r, -r, r, r)], self.plane, low, high)

    def specs(self):
        """The reference's solids whose common the piece is: the solid and
        the half-space box, a zone's whole sphere and its slab between its
        ends' parallels (the chained reference takes whole spheres)."""
        spec = self.solid
        if spec[0] != 'sphere' or (spec[3], spec[4]) == (-HP, HP):
            return [spec, self.halfspace()]
        _, r, frame, low, high = spec
        # The ends' heights as the kernel stores them: `r sin(lat)`, the
        # sine rounded once and the product once.
        w = lambda lat, pole: pole*2*r if lat == pole*HP else r*sin_rn(lat)
        slab = prism([square(-self.reach, -self.reach, self.reach, self.reach)], frame, w(low, -1), w(high, 1))
        return [sphere(r, frame), self.halfspace(), slab]

    def expr(self, k):
        """The piece's expression over its solids from index `k`."""
        out = ('common', k, k+1)
        return ('common', out, k+2) if len(self.specs()) == 3 else out

    def case(self, name, op):
        c = construction(self.solid, op)
        c.name = name
        c.split = (self.plane, self.side)
        return c


def farthest(piece):
    """An upper bound of the distance from the piece's plane's origin to
    its solid: a prism's profile's bounding points (an arc's or a circle's
    its centre's square of side `2 r`) at both ends, a revolved solid's
    centre (a cone's ends' centres) plus its widest radius."""
    spec = piece.solid
    po = stored_axes(piece.plane)[0]
    frame = pieces.spec_frame(spec) if spec[0] != 'torus' else spec[3]
    o, x, y, n = stored_axes(frame)
    at = lambda u, v, w: [o[i]+u*x[i]+v*y[i]+w*n[i] for i in range(3)]
    dist = lambda p: mp.sqrt(sum((mp.mpf(p[i])-mp.mpf(po[i]))**2 for i in range(3)))
    if spec[0] == 'prism':
        _, boundaries, _, start, end = spec
        pts = []
        for b in boundaries:
            if b.circle is not None:
                cx, cy, r = b.circle
                pts += [(cx+sx*r, cy+sy*r) for sx in (-1, 1) for sy in (-1, 1)]
            else:
                pts += list(b.points)
                for seg in b.segments or []:
                    if seg is not None:
                        cx, cy, r = seg[:3]
                        pts += [(cx+sx*r, cy+sy*r) for sx in (-1, 1) for sy in (-1, 1)]
        return max(dist(at(u, v, w)) for u, v in pts for w in (start, end))
    if spec[0] == 'cone':
        _, r0, r1, h, _ = spec
        return max(dist(at(0, 0, 0))+r0, dist(at(0, 0, h))+r1)
    if spec[0] == 'torus':
        return dist(at(0, 0, 0))+spec[1]+spec[2]
    return dist(at(0, 0, 0))+spec[1]


def exact(v):
    return [mp.mpf(F(c).numerator)/F(c).denominator for c in v]


def cylinder_cut_closed(r, axis_point, plane):
    """A cylinder of radius `r` about the upright axis through `axis_point`
    (the low end's centre) below a plane crossing its wall only: the disc's
    area times the plane's height over the axis, its area the disc, the wall
    (its mean height that height) and the ellipse `pi r^2 |m| / |m_z|`, `m`
    the plane's normal as the half-space box's face holds it, its stored
    axes' `x * y`."""
    R = mp.mpf(r)
    o, x, y, _ = stored_axes(plane)
    x, y = exact(x), exact(y)
    m = [x[1]*y[2]-x[2]*y[1], x[2]*y[0]-x[0]*y[2], x[0]*y[1]-x[1]*y[0]]
    p = exact(o)
    a = exact(axis_point)
    # The plane's height over the axis: m . (X - p) = 0 at X = (a0, a1, z).
    z = p[2]-(m[0]*(a[0]-p[0])+m[1]*(a[1]-p[1]))/m[2]
    h = z-a[2]
    disc_area = mp.pi*R**2
    return disc_area*h, disc_area+2*mp.pi*R*h+disc_area*mp.sqrt(sum(c*c for c in m))/abs(m[2])


CYL = prism([disc(0.0, 0.0, 3.0)], at('XY', (5, 5, -1)), 0.0, 8.0)
CYL_PLANE = at('LEAN', (5, 5, 3))
DEE = prism([base.dee(0.0, 5.0)], at('XY', (0, 0, 0)), 0.0, 4.0)
BLOCK = prism([square(1.0, 1.0, 9.0, 7.0)], at('XY', (0, 0, 0)), 0.0, 5.0)
FRUSTUM = cone(3.0, 1.5, 6.0, at('XY', (5, 5, 0)))
ZONE = sphere(4.0, at('SKEW2', (5, 5, 5)), -0.5, 0.75)
# The plane through the zone's axis normal to its frame's `x` (`SKEW2`'s
# `x` its normal, its `n` the hint).
ZONE_PLANE = (5.0, 5.0, 5.0, 1.0, -8.0, 4.0, 4.0, 4.0, 7.0)
ZONE2 = sphere(4.5, at('SKEW', (5, 5, 4)), -0.75, 1.0)
TORUS = torus(3.0, 1.25, at('XY', (5, 5, 4)))
# Normal to the torus's axis (the world's `z`), through a point above its
# centre.
BAND_PLANE = at('XY', (5, 5, 4.5))


def make_pieces():
    out = [
        Piece('cyl_low', 'clipped', CYL, CYL_PLANE, 'below', 7.5,
              closed=lambda: cylinder_cut_closed(3, (5, 5, -1), CYL_PLANE)),
        Piece('dee_low', 'clipped', DEE, at('TILT', (7, 5, 2)), 'below', 12.625),
        Piece('block_high', 'clipped', BLOCK, at('LEAN', (5, 4, 2.5)), 'above', 7.5),
        Piece('frustum_low', 'conic', FRUSTUM, at('TILT', (5, 5, 3)), 'below', 8.0),
        Piece('zone_half', 'half', ZONE, ZONE_PLANE, 'above', 5.5),
        Piece('zone_cut', 'conic', ZONE2, at('TILT', (5, 5, 5)), 'above', 7.0),
        Piece('band', 'band', TORUS, BAND_PLANE, 'above', 6.0),
        Piece('cone_axis', 'half', FRUSTUM, at('SIDE', (5, 5, 3)), 'above', 8.0),
    ]
    return {p.name: p for p in out}


PIECES = make_pieces()


def split(name):
    return ('split', name)


def flatten(items):
    """The reference's solids of a case's inputs (a split piece its solid
    and half-space box) and each input's expression over them."""
    specs, exprs = [], []
    for it in items:
        if it[0] == 'split':
            k = len(specs)
            specs += PIECES[it[1]].specs()
            exprs.append(PIECES[it[1]].expr(k))
        else:
            exprs.append(len(specs))
            specs.append(it)
    return specs, exprs


class Split:
    """One case: its inputs (`split(NAME)` or a construction), the first
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
            if s[0] == 'split':
                c = PIECES[s[1]].case(name, OPERATIONS[k])
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
        out.append(Split(f'{name}_{op}', klass, items, first or op, stages, k, n, reason))
    return out


APEX = 'a frustum\'s half by a plane through its axis (through its virtual apex, S9d.3a\'s rule)'
FLUSH = ('a box on the cylinder piece\'s cut plane (one plane in the construction, within the resolution of it in '
         'the kernel\'s models)')
# A rod along the world's `y` (its frame's normal `y`, `x` the world's `x`:
# a profile point `(u, v)` lies at `x = u`, `z = -v`).
ALONG_Y = (0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0)


def rod_y(x, z, r, y0=0.0, y1=12.0):
    return prism([disc(x, -z, r)], ALONG_Y, y0, y1)


ROD = rod_y(6.5, 2.25, 0.875)
SLAB = prism([square(1.5, 1.0, 9.0, 9.5)], at('XY', (0, 0, 0.875)), 0.0, 1.25)


def cases():
    out = []
    three = {'fuse': 1, 'cut': 1, 'common': 1}
    cyl_box = prism([square(6.25, 3.5, 9.5, 6.75)], at('XY', (0, 0, -0.25)), 0.0, 3.0)
    out += group('cyl_box', 'prism', [split('cyl_low'), cyl_box], three)
    out += group('box_cyl', 'prism', [cyl_box, split('cyl_low')], three)
    out += group('cyl_ball', 'prism', [split('cyl_low'), sphere(2.0, at('SKEW', (3.25, 6.5, 3.75)))], three)
    out += group('cyl_rod', 'prism', [split('cyl_low'), ROD], three)
    out += group('dee_rod', 'prism', [split('dee_low'), prism([disc(14.0, 6.5, 1.25)], at('XY', (0, 0, -1)),
                                                               0.0, 6.0)], three)
    out += group('block_rod', 'prism', [split('block_high'), rod_y(6.25, 3.25, 1.0, -1.0, 9.0)], three)
    out += group('frustum_ball', 'revolved', [split('frustum_low'), sphere(1.75, at('SKEW', (6.375, 3.3125, 4.1875)))],
                 three)
    out += group('zone_box', 'revolved', [split('zone_half'), prism([square(5.75, 1.5, 9.25, 4.75)],
                                                                    at('XY', (0, 0, 3.25)), 0.0, 3.25)], three)
    out += group('zone_cone', 'revolved', [split('zone_half'), cone(1.25, 1.0, 3.0, at('SKEW2', (7.875, 3.4375, 3.0)))],
                 three)
    out += group('cut_rod', 'revolved', [split('zone_cut'), rod_y(4.5, 5.25, 1.0, -1.0, 11.0)], three)
    out += group('band_ball', 'revolved', [split('band'), sphere(0.75, at('SKEW', (7.875, 5.25, 5.125)))], three)
    out += group('pair', 'both', [split('cyl_low'), split('dee_low')], three)
    out += group('chain_cyl', 'chain', [split('cyl_low'), ROD, SLAB], {'fuse': 1, 'cut': 2, 'common': 1},
                 first='cut')
    out += group('cone_axis', 'revolved', [split('cone_axis'), prism([square(3.5, 4.25, 6.75, 7.5)],
                                                                     at('XY', (0, 0, 1.25)), 0.0, 2.5)],
                 three, reason=APEX)
    out += group('cyl_flush', 'prism', [split('cyl_low'), prism([square(-1.5, -1.25, 1.75, 1.5)], CYL_PLANE,
                                                                 0.0, 2.0)],
                 {'fuse': 1, 'cut': 1, 'common': 0}, reason=FLUSH)
    out += group('one_sphere', 'revolved', [split('zone_half'), sphere(4.0, at('SKEW2', (5, 5, 5)))],
                 {'fuse': 1, 'cut': 0, 'common': 1})
    return out


def all_cases():
    return cases()


def validate(listed):
    names = [c.name for c in listed]
    assert len(names) == len(set(names)), 'duplicate case names'
    for c in listed:
        assert len(c.items) == 2+len(c.stages), c.name
        assert c.last in OPS, c.name
        assert any(s[0] == 'split' for s in c.items), c.name
        if c.klass == 'both':
            assert sum(s[0] == 'split' for s in c.items) == 2, c.name
        if c.klass == 'chain':
            assert c.stages, c.name
        specs, _ = flatten(c.items)
        for s in specs:
            assert s[0] in ('sphere', 'prism', 'cone', 'torus'), c.name


# ------------------------------------------------------------------ the splits' classes

def dot(a, b):
    return sum(p*q for p, q in zip(a, b))


def split_class(piece):
    """The kernel's class of a split, from the solid's and the plane's
    stored axes exactly: a prism's plane oblique to its axis (`clipped`), a
    cone's or zone's through its axis within a quarter of the resolution
    over its ends (`half`) or not (`conic`), a torus's normal to its axis
    (`band`)."""
    spec = piece.solid
    frame = pieces.spec_frame(spec) if spec[0] != 'torus' else spec[3]
    o, x, y, n = (tuple(F(c) for c in v) for v in stored_axes(frame))
    po, _, _, m = (tuple(F(c) for c in v) for v in stored_axes(piece.plane))
    a, b, c = dot(m, x), dot(m, y), dot(m, n)
    d = dot(m, tuple(o[i]-po[i] for i in range(3)))
    tol = F(1, 10**7)
    if spec[0] == 'prism':
        return 'clipped' if c != 0 and (a, b) != (0, 0) else 'other'
    if spec[0] == 'torus':
        R, r = F(spec[1]), F(spec[2])
        normal = 16*(a*a+b*b)*(R+r)**2 <= tol*tol*c*c
        return 'band' if normal else 'other'
    m2 = a*a+b*b+c*c
    if spec[0] == 'cone':
        ends = (F(0), F(spec[3]))
    else:
        r = spec[1]
        # The zone's ends' heights as the kernel stores them: `r sin(lat)`
        # rounded once.
        ends = tuple(F(float(mp.mpf(r)*mp.sin(mp.mpf(lat)))) for lat in spec[3:5])
    through = all(16*(c*w+d)**2 <= tol*tol*m2 for w in ends)
    return 'half' if through else 'conic'


# ------------------------------------------------------------------ the reference

def owners(items, sub):
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
    # Each split piece one solid, its closed form where it has one.
    for it, e in zip(items, sub):
        if it[0] != 'split':
            continue
        assert ref.count_solids(chain, e, 48) == 1, f'{name}: {it[1]} is not one solid'
        closed = PIECES[it[1]].closed
        if closed is not None:
            V, A = closed()
            v, _, a = chain.measures(e)
            d = max(abs(v-V)/size**3, abs(a-A)/size**2)
            checks['piece_closed_forms'] = max(checks.get('piece_closed_forms', mp.mpf(0)), d)
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


def jobs(mc_n):
    """One job per group."""
    groups = {}
    for c in all_cases():
        e = groups.setdefault(c.group, [c.items, c.op1 if c.stages else None, bool(c.stages and c.stages[0][1]),
                                        bool(c.stages), {}, c.kind in ('degenerate', 'unsupported')])
        e[4][c.last] = c.solids
    out = [(name, items, op1, swapped, chained, ops, degenerate, mc_n)
           for name, (items, op1, swapped, chained, ops, degenerate) in groups.items()]
    # The costliest first (more solids, spheres' wide sweeps), for the
    # workers' balance: the results are keyed by group.
    cost = lambda job: (sum(sp[0] == 'sphere' for sp in flatten(job[1])[0]), len(flatten(job[1])[0]))
    return sorted(out, key=cost, reverse=True)


def cached(job):
    """`evaluate`, its result kept in `--cache` (a development aid)."""
    cache = os.environ.get('SPLIT_PIECES_CACHE')
    if not cache:
        return evaluate(job)
    import hashlib
    import inspect
    import pickle
    from pathlib import Path
    key = hashlib.sha256((repr(job)+repr(flatten(job[1])[0])+inspect.getsource(evaluate)
                          + inspect.getsource(cylinder_cut_closed)+Path(ref.__file__).read_text()).encode()).hexdigest()[:24]
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
    frames = ['# case\tframe (solid k, or split k its plane) axis (n, x, y)\tstored unit vector (stored_axes, as hex '
              'bits)']
    for case in all_cases():
        for k, (item, c) in enumerate(zip(case.items, case.solid_cases)):
            listed = [(f'solid{k}', c.frame)]
            if item[0] == 'split':
                listed.append((f'split{k}', PIECES[item[1]].plane))
            for which, frame in listed:
                _, x, y, n = stored_axes(frame)
                for key, v in (('n', n), ('x', x), ('y', y)):
                    frames.append(f'{case.name}\t{which} {key}\t'+' '.join(struct.pack('>d', q).hex() for q in v))
    return {'boolean-split-pieces-cases.txt': '\n'.join(blocks)+'\n',
            'boolean-split-pieces-expected.tsv': '\n'.join(out)+'\n',
            'boolean-split-pieces-frames.tsv': '\n'.join(frames)+'\n'}


def native_input():
    """Every case's native rows (`compare_split_pieces_boolean.py`)."""
    return '\n'.join(c.native() for c in all_cases())+'\n'


def case_size(case):
    """The case's size: its reference solids' reach from the origin."""
    specs, _ = flatten(case.items)
    chain = ref.Chain([construction(s, 0) for s in specs])
    return float(chain.size)


def check_splits():
    """Each piece's split is its class's."""
    for p in PIECES.values():
        assert split_class(p) == p.kind, (p.name, split_class(p), p.kind)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--workers', type=int, default=max(1, min(4, (os.cpu_count() or 2)-2)))
    parser.add_argument('--samples', type=int, default=50000, help='Monte-Carlo points per group')
    parser.add_argument('--only', help='one group, printed (no files)')
    parser.add_argument('--cache', help='a directory keeping each group\'s result (development)')
    args = parser.parse_args()
    if args.cache:
        os.environ['SPLIT_PIECES_CACHE'] = args.cache
    listed = all_cases()
    validate(listed)
    check_splits()
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
    kinds = {}
    for c in listed:
        kinds[c.kind] = kinds.get(c.kind, 0)+1
    print(len(listed), 'cases:', ', '.join(f'{sum(1 for c in listed if c.last == op)} {op}' for op in OPS),
          '-', ', '.join(f'{v} {k}' for k, v in sorted(kinds.items())),
          '-', ', '.join(f'{sum(1 for c in listed if c.klass == k)} {k}' for k in CLASSES))
    print('reference checks (largest deviation, relative to the case size; least_sine and event_spacing as '
          'reciprocals):', ', '.join(f'{k} {mp.nstr(v, 3)}' for k, v in sorted(worst.items())))
    print('cases per check:', ', '.join(f'{k} {v}' for k, v in sorted(covered.items())))


if __name__ == '__main__':
    main()
