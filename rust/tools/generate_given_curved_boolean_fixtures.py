#!/usr/bin/env python3
"""Fixtures for S9e.3a of REVIEW_NOTES.md: a Boolean's result given to
further Booleans where the solids are spheres, cones or frusta and whole tori
besides prisms with arcs, the given result's edges any its first arrangement
makes (a sphere's circles, a cone's and a torus's plane sections, meetings of
two curved faces the partner does not reach), deeper chains (a result of a
result given again), and given results against a sphere, a cone and a torus.

`boolean-given-curved-cases.txt` lists each case in the Boolean protocol
with its `then` rows (`identity_reference.encode_chained_case`, S9e.3: a
deeper chain's further `then` rows and solids); `boolean-given-curved-
expected.tsv` gives per case, from `chained_curved_boolean_reference.py`:

* `expect KIND S9e.3a CLASS`: the declared outcome (`solid`: one or more
  solids; `empty`; `degenerate`: `Degenerate` in the decisions) and the
  chain's class (`sphere`, `cone`, `torus`: the given result holds one;
  `procedural`: its meeting of two curved faces, which the partner does not
  reach; `deep`: a result of a result given; `partner`: a given result
  against a sphere, a cone or a torus), then for a degenerate case `reason
  TEXT`;
* `result N volume area cx cy cz` (totals over the N solids, world
  coordinates) or `empty`.

`boolean-given-curved-frames.tsv` records every solid's stored axes. Frames
are the curved generator's whose stored axes the kernel gives bit for bit
(`XY`, `SIDE`, `TILT`, `R125`).

The chains (`B` the box `[0, 10]^2 x [0, 4]` unless said): `dome_tilt`, `B`
fused with a sphere of radius 3 about the centre of its top face, then with a
slab of `TILT` across the dome and the top face's circle; `dome_drill`, a
cylinder of radius 1 through the dome's circle with the domed box (swapped:
the domed box the tool); DRAW's `bcut_simple/G9` body (a cylinder of radius
9 and height 3 fused with a frustum of radii 7 and 6 and height 4 on its
base) with G9's `pcylinder` moved clear of the frustum's top circle
(`g9_clear`, the axis at `x = 4.5`), with a box whose walls cross the
frustum's section by the cylinder's top (`g9_box`), and with a coaxial
cylinder of radius 2 (`g9_bore`, closed forms); `B` less a frustum from
radius 1 at height 1 to radius 3 at its top (a countersink) and a coaxial
cylinder of radius 1/2 (`countersink_drill`, closed forms); a box less a
torus of radii 3 and 1 whose tube its top face cuts above the equator (a
groove of surd radii), with a `TILT` slab across the groove (`groove_tilt`)
and a cylinder through its outer circle (`groove_drill`); a sphere fused
with a cylinder off its centre (their meeting a `Rise` curve), then with a
box slicing the sphere's bottom clear of it (`peg_clear`); deeper chains:
`B` of height 5 less a hole along `x` and a vertical hole, then with a
`TILT` slab across both (`holes_tilt`), and the domed box drilled, then with
an `R125` box whose wall crosses the dome clear of the drill (`dome_deep`);
`B` of height 5 less the hole along `x` against a sphere about its top face's
centre (`holed_ball`), a frustum standing in it through the hole
(`holed_cone`) and a torus about a vertical axis whose tube crosses the hole
(`holed_torus`); declared `degenerate`: DRAW's `G9` itself (the `pcylinder`
of radius 1 about `(5, 0)` touching the frustum's top circle of radius 6 at
`(6, 0, 4)`) and a box whose bottom face touches the dome's top point
(`dome_touch`).

Before writing, `reference_checks` compares the reference with independent
results (limits relative to the case's size): coaxial closed forms (the
cross sections' rings about the axis integrated exactly by Simpson's rule
between their kinks, their horizontal faces and walls where the chain's set
function changes across them) within 1e-30; the two families of curves of
every face (independent roots and events) within 1e-30; each solid's
closed form; for the given result `X` and the last solid `C` the pair
identities `V(X u C) + V(X n C) = V(X) + V(C)`, `V(X - C) = V(X) - V(X n
C)` (swapped `V(C - X) = V(C) - V(X n C)`), moments likewise, and `area(X u
C) + area(X n C) = area(X) + area(C)` where no face of `C` lies on another's
surface, within 1e-30; the given result one solid; Monte-Carlo estimates
(200,000 uniform points per chain) of every volume and centre within 5
standard errors; the declared solid counts given by rays at two
resolutions; no quadrature node off its interval's structure. Margins
outside the declared cases: every surface met at a sine of at least 0.05
on the scanned curves, every family's events at least 1e-6 of its range
apart (no near coincidence). Chains run in worker processes (`--workers`).
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
import os
from pathlib import Path
import struct
import zlib

import mpmath as mp

from identity_reference import Case, encode_chained_case
from curve_surface_reference import stored_axes
import chained_curved_boolean_reference as ref
from generate_curved_boolean_fixtures import FRAMES, at, disc, square
from torus_boolean_reference import TWO_PI

ROOT = Path(__file__).resolve().parents[1]
OPERATIONS = (91, 92, 94, 96)      # the solids' construction ids
BOOLEANS = (93, 95, 97)            # the Booleans' ids
CLASSES = ('sphere', 'cone', 'torus', 'procedural', 'deep', 'partner')
HP = 1.5707963267948966
OPS = ('fuse', 'cut', 'common')


def solid(spec, op):
    kind = spec[0]
    if kind == 'torus':
        _, R, r, frame = spec
        return Case('', 1e-7, op, frame, 0.0, 0.0, [], torus=(float(R), float(r), 0.0, TWO_PI, TWO_PI))
    if kind == 'sphere':
        _, r, frame = spec
        return Case('', 1e-7, op, frame, 0.0, 0.0, [], sphere=(float(r), -HP, HP))
    if kind == 'cone':
        _, r0, r1, h, frame = spec
        return Case('', 1e-7, op, frame, 0.0, 0.0, [], cone=(float(r0), float(r1), float(h)))
    _, boundaries, frame, start, end = spec
    return Case('', 1e-7, op, frame, float(start), float(end), boundaries)


def torus(R, r, frame):
    return ('torus', R, r, frame)


def sphere(r, frame):
    return ('sphere', r, frame)


def cone(r0, r1, h, frame):
    return ('cone', r0, r1, h, frame)


def prism(boundaries, frame, start, end):
    return ('prism', boundaries, frame, start, end)


def spec_frame(spec):
    return spec[{'torus': 3, 'sphere': 2, 'cone': 4, 'prism': 2}[spec[0]]]


class Given:
    """One case: the solids, the first Boolean, the fixed stages and the
    last stage's operation."""

    def __init__(self, name, klass, specs, op1, stages, kind, solids, reason=None):
        assert klass in CLASSES, klass
        self.name, self.klass, self.op1, self.stages = name, klass, op1, stages
        self.kind, self.solids, self.reason = kind, solids, reason
        self.chain_name = name.rsplit('_', 1)[0]
        self.specs = specs
        self.solid_cases = [solid(s, OPERATIONS[k]) for k, s in enumerate(specs)]
        for c in self.solid_cases:
            c.name = name
        self.frames = [spec_frame(s) for s in specs]

    @property
    def last(self):
        return self.stages[-1][0]

    def expr(self):
        (op2, sw2), more = self.stages[0], self.stages[1:]
        return ref.chain_expr(self.op1, op2, sw2, more)

    def given_expr(self):
        """The result given to the last Boolean."""
        expr = self.expr()
        return expr[1] if isinstance(expr[2], int) else expr[2]

    def encode(self):
        c = self.solid_cases
        (op2, sw2) = self.stages[0]
        more = [(op, c[3+k], BOOLEANS[2+k], sw, None) for k, (op, sw) in enumerate(self.stages[1:])]
        return encode_chained_case(c[0], self.op1, c[1], BOOLEANS[0], op2, c[2], BOOLEANS[1], sw2, None,
                                   more)


def group(name, klass, specs, op1, stages, outcomes, swapped=False, reason=None):
    """Cases of one chain: `stages` the fixed stages after the first
    Boolean (`(op, swapped)`), the last stage's operation varied over
    `outcomes` (`{op: solids}`, 0 for empty; a string kind for a declared
    degenerate case)."""
    out = []
    for op, n in outcomes.items():
        kind = 'degenerate' if reason else ('solid' if n else 'empty')
        out.append(Given(f'{name}_{op}', klass, specs, op1, list(stages)+[(op, swapped)], kind, n, reason))
    return out


# The box B and the domed box (a sphere of radius 3 about B's top face's
# centre), the dome's circle of radius 3 on the top face.
BOX = prism([square(0.0, 0.0, 10.0, 10.0)], at('XY', (0, 0, 0)), 0.0, 4.0)
BALL = sphere(3.0, at('XY', (5, 5, 4)))
# DRAW's bcut_simple/G9 body: pcylinder cyl 9 3, pcone kone 7 6 4, bfuse.
G9_CYL = prism([disc(0.0, 0.0, 9.0)], at('XY', (0, 0, 0)), 0.0, 3.0)
G9_CONE = cone(7.0, 6.0, 4.0, at('XY', (0, 0, 0)))
# B of height 5 less a hole of radius 3/2 along x about (y, z) = (5, 5/2).
BOX5 = prism([square(0.0, 0.0, 10.0, 10.0)], at('XY', (0, 0, 0)), 0.0, 5.0)
HOLE_X = prism([disc(5.0, 2.5, 1.5)], at('SIDE', (0, 0, 0)), -1.0, 11.0)
# A box less a torus of radii 3 and 1 about the origin, its top at z = 1/2
# (the groove's circles of radii 3 -+ sqrt(3)/2).
GROOVE = (prism([square(-5.0, -5.0, 5.0, 5.0)], at('XY', (0, 0, 0)), -2.0, 0.5),
          torus(3.0, 1.0, at('XY', (0, 0, 0))))
TOUCH = 'a third box\'s bottom face tangent to the given result\'s sphere at its top point'
G9_TANGENT = ('the third cylinder of radius 1 about (5, 0) touching the given frustum\'s top circle of radius 6 '
              'at (6, 0, 4) from inside (DRAW\'s bcut_simple/G9 and H3)')


def cases():
    out = []
    out += group('dome_tilt', 'sphere', [BOX, BALL, prism([square(-2.0, -15.0, 12.0, 15.0)],
                                                          at('TILT', (0, 0, 4)), 1.5, 4.0)],
                 'fuse', [], {'fuse': 1, 'cut': 2, 'common': 1})
    out += group('dome_drill', 'sphere', [BOX, BALL, prism([disc(7.5, 5.0, 1.0)], at('XY', (0, 0, 0)), -1.0, 8.0)],
                 'fuse', [], {'fuse': 1, 'cut': 2, 'common': 1}, swapped=True)
    out += group('g9_clear', 'cone', [G9_CYL, G9_CONE, prism([disc(0.0, 0.0, 1.0)], at('XY', (4.5, 0, -2)), 0.0, 9.0)],
                 'fuse', [], {'fuse': 1, 'cut': 1, 'common': 1})
    out += group('g9_box', 'cone', [G9_CYL, G9_CONE, prism([square(3.0, -2.0, 12.0, 2.0)], at('XY', (0, 0, 0)),
                                                           -1.0, 5.0)],
                 'fuse', [], {'fuse': 1, 'cut': 1, 'common': 1})
    out += group('g9_bore', 'cone', [G9_CYL, G9_CONE, prism([disc(0.0, 0.0, 2.0)], at('XY', (0, 0, 0)), -1.0, 5.0)],
                 'fuse', [], {'fuse': 1, 'cut': 1, 'common': 1})
    out += group('countersink_drill', 'cone', [BOX, cone(1.0, 3.0, 3.0, at('XY', (5, 5, 1))),
                                               prism([disc(5.0, 5.0, 0.5)], at('XY', (0, 0, 0)), -1.0, 6.0)],
                 'cut', [], {'fuse': 1, 'cut': 1, 'common': 1})
    out += group('groove_tilt', 'torus', [*GROOVE, prism([square(-6.0, -15.0, 6.0, 15.0)], at('TILT', (0, 0, 0)),
                                                         1.8, 2.4)],
                 'cut', [], {'fuse': 1, 'cut': 2, 'common': 1})
    out += group('groove_drill', 'torus', [*GROOVE, prism([disc(4.0, 0.0, 0.6)], at('XY', (0, 0, 0)), -3.0, 2.0)],
                 'cut', [], {'fuse': 1, 'cut': 1, 'common': 1})
    out += group('peg_clear', 'procedural', [sphere(3.0, at('XY', (0, 0, 0))),
                                             prism([disc(1.25, 0.0, 1.0)], at('XY', (0, 0, 0)), 0.0, 6.0),
                                             prism([square(-4.0, -4.0, 4.0, 4.0)], at('XY', (0, 0, 0)), -5.0, -1.5)],
                 'fuse', [], {'fuse': 1, 'cut': 1, 'common': 1})
    out += group('holes_tilt', 'deep', [BOX5, HOLE_X, prism([disc(3.0, 1.5, 1.0)], at('XY', (0, 0, 0)), -1.0, 6.0),
                                        prism([square(-2.0, -15.0, 12.0, 15.0)], at('TILT', (0, 0, 0)), 2.0, 4.5)],
                 'cut', [('cut', False)], {'fuse': 1, 'cut': 2, 'common': 1})
    out += group('dome_deep', 'deep', [BOX, BALL, prism([disc(7.5, 5.0, 1.0)], at('XY', (0, 0, 0)), -1.0, 8.0),
                                       prism([square(-20.0, -20.0, 6.0, 20.0)], at('R125', (0, 0, -1)), 0.0, 10.0)],
                 'fuse', [('cut', False)], {'fuse': 1, 'cut': 1, 'common': 1})
    out += group('holed_ball', 'partner', [BOX5, HOLE_X, sphere(3.0, at('XY', (5, 5, 5)))],
                 'cut', [], {'fuse': 1, 'cut': 1, 'common': 1})
    out += group('holed_cone', 'partner', [BOX5, HOLE_X, cone(2.0, 0.5, 4.0, at('XY', (5, 5, 3)))],
                 'cut', [], {'fuse': 1, 'cut': 1, 'common': 1})
    out += group('holed_torus', 'partner', [BOX5, HOLE_X, torus(3.0, 1.0, at('XY', (5, 5, 4.5)))],
                 'cut', [], {'fuse': 1, 'cut': 1, 'common': 1})
    out += group('g9', 'cone', [G9_CYL, G9_CONE, prism([disc(0.0, 0.0, 1.0)], at('XY', (5, 0, -2)), 0.0, 9.0)],
                 'fuse', [], {'fuse': 1, 'cut': 1, 'common': 1}, reason=G9_TANGENT)
    out += group('dome_touch', 'sphere', [BOX, BALL, prism([square(3.0, 3.0, 7.0, 7.0)], at('XY', (0, 0, 0)),
                                                           7.0, 9.0)],
                 'fuse', [], {'fuse': 2, 'cut': 1, 'common': 0}, reason=TOUCH)
    return out


def validate(listed):
    names = [c.name for c in listed]
    assert len(names) == len(set(names)), 'duplicate case names'
    for c in listed:
        for frame in c.frames:
            assert frame_name(frame) in FRAMES, c.name
        assert len(c.specs) == 2+len(c.stages), c.name
        kinds = {s[0] for s in c.specs}
        if c.klass in ('sphere', 'cone', 'torus'):
            assert c.klass in kinds, c.name
        if c.klass == 'deep':
            assert len(c.stages) > 1, c.name
        if c.klass == 'partner':
            assert c.specs[-1][0] in ('sphere', 'cone', 'torus'), c.name


def frame_name(frame):
    for k, v in FRAMES.items():
        if tuple(frame[3:]) == v:
            return k
    raise KeyError(frame)


# ------------------------------------------------------------------ coaxial closed forms

class Coaxial:
    """Solids about one vertical axis through `(cx, cy)`: discs of radius
    linear in `z` (a cylinder, a cone or frustum) over `[z0, z1]`, and a
    square of side `s` centred on the axis holding every disc."""

    def __init__(self, cx, cy):
        self.cx, self.cy = mp.mpf(cx), mp.mpf(cy)
        self.discs = []    # (index, r0, r1, z0, z1)
        self.square = None

    def disc(self, index, r0, r1, z0, z1):
        self.discs.append((index, mp.mpf(r0), mp.mpf(r1), mp.mpf(z0), mp.mpf(z1)))

    def box(self, index, side, z0, z1):
        self.square = (index, mp.mpf(side), mp.mpf(z0), mp.mpf(z1))

    def radius(self, d, z):
        _, r0, r1, z0, z1 = d
        return r0+(r1-r0)*(z-z0)/(z1-z0)

    def levels(self):
        zs = set()
        for _, _, _, z0, z1 in self.discs:
            zs |= {z0, z1}
        if self.square:
            zs |= {self.square[2], self.square[3]}
        # Where two discs' radii cross.
        for a in self.discs:
            for b in self.discs:
                if a is b:
                    continue
                lo, hi = max(a[3], b[3]), min(a[4], b[4])
                if hi <= lo:
                    continue
                fa, fb = self.radius(a, lo)-self.radius(b, lo), self.radius(a, hi)-self.radius(b, hi)
                if fa*fb < 0:
                    zs.add(lo+(hi-lo)*fa/(fa-fb))
        return sorted(zs)

    def rings(self, z, side):
        """The section's rings at `z` (just above it for `side` 1, below for
        -1): [(inner, outer, memberships)]; outer None for the square's rest."""
        n = 1+max([d[0] for d in self.discs]+([self.square[0]] if self.square else []))
        eps = mp.mpf(10)**-30*side
        zz = z+eps
        live = [d for d in self.discs if d[3] < zz < d[4]]
        radii = sorted({self.radius(d, z) for d in live} | {mp.mpf(0)})
        out = []
        bounds = list(zip(radii, radii[1:]))
        if self.square:
            bounds.append((radii[-1], None))
        for a, b in bounds:
            m = [False]*n
            rmid = (a+b)/2 if b is not None else None
            for d in live:
                if rmid is not None and rmid < self.radius(d, zz):
                    m[d[0]] = True
            if self.square:
                m[self.square[0]] = self.square[2] < zz < self.square[3]
            out.append((a, b, m))
        return out

    def ring_area(self, a, b):
        if b is None:
            return self.square[1]**2-mp.pi*a*a
        return mp.pi*(b*b-a*a)

    def measures(self, expr):
        """(volume, centre, area) of the expression's set."""
        zs = self.levels()
        V = Mz = A = mp.mpf(0)
        # Volume and the z moment by Simpson's rule between levels (the
        # section's area is quadratic in z there, z times it cubic).
        for z0, z1 in zip(zs, zs[1:]):
            zm = (z0+z1)/2
            # The interval's rings (their order and memberships at its
            # middle), their radii at each of Simpson's heights.
            s0, sm, s1 = (sec_over(self, expr, zm, z) for z in (z0, zm, z1))
            V += (z1-z0)*(s0+4*sm+s1)/6
            Mz += (z1-z0)*(z0*s0+4*zm*sm+z1*s1)/6
            # Walls: each disc's wall where the set changes across it.
            for d in self.discs:
                if not (d[3] <= z0 and z1 <= d[4]):
                    continue
                mid = self.radius(d, zm)
                ins, outs = members(self, zm, mid*(1-mp.mpf(10)**-20)), members(self, zm, mid*(1+mp.mpf(10)**-20))
                if ref.evaluate(expr, ins) != ref.evaluate(expr, outs):
                    r0, r1 = self.radius(d, z0), self.radius(d, z1)
                    k = (d[2]-d[1])/(d[4]-d[3])
                    A += mp.pi*(r0+r1)*(z1-z0)*mp.sqrt(1+k*k)
            if self.square and self.square[2] <= z0 and z1 <= self.square[3]:
                half = self.square[1]/2
                ins, outs = members(self, zm, None, half*(1-mp.mpf(10)**-20)), members(self, zm, None, half*2)
                if ref.evaluate(expr, ins) != ref.evaluate(expr, outs):
                    A += 4*self.square[1]*(z1-z0)
        # Horizontal faces at each level: rings where the set differs above
        # and below.
        for z in zs:
            below = {(a, b): ref.evaluate(expr, m) for a, b, m in self.rings(z, -1)}
            above = {(a, b): ref.evaluate(expr, m) for a, b, m in self.rings(z, 1)}
            radii = sorted({x for a, b in list(below)+list(above) for x in (a, b) if x is not None})
            pieces = list(zip(radii, radii[1:]))+([(radii[-1], None)] if self.square else [])
            for a, b in pieces:
                mid = (a+b)/2 if b is not None else None
                if ring_value(below, mid) != ring_value(above, mid):
                    A += self.ring_area(a, b)
        return V, (self.cx, self.cy, Mz/V), A


def members(co, z, r, sq=None):
    """The memberships at height `z` and radius `r` (or at `sq` along x from
    the axis, the square's side region)."""
    n = 1+max([d[0] for d in co.discs]+([co.square[0]] if co.square else []))
    m = [False]*n
    rho = r if r is not None else sq
    for d in co.discs:
        if d[3] < z < d[4] and rho < co.radius(d, z):
            m[d[0]] = True
    if co.square:
        half = co.square[1]/2
        m[co.square[0]] = co.square[2] < z < co.square[3] and rho < half
    return m


def sec_over(co, expr, zm, z):
    """The section's area at `z` with the rings of the interval holding
    `zm`: the discs live there in their order of radius at `zm`, each
    ring's membership at `zm`, its radii at `z`."""
    live = sorted((d for d in co.discs if d[3] < zm < d[4]), key=lambda d: co.radius(d, zm))
    tot = mp.mpf(0)
    inner = (mp.mpf(0), mp.mpf(0))
    bounds = [(inner, (co.radius(d, zm), co.radius(d, z))) for d in live[:1]]
    bounds += [((co.radius(a, zm), co.radius(a, z)), (co.radius(b, zm), co.radius(b, z)))
               for a, b in zip(live, live[1:])]
    if co.square:
        last = (co.radius(live[-1], zm), co.radius(live[-1], z)) if live else inner
        bounds.append((last, None))
    for a, b in bounds:
        rmid = (a[0]+b[0])/2 if b is not None else None
        m = members(co, zm, rmid, None if rmid is not None else co.square[1]/2*(1-mp.mpf(10)**-20))
        if ref.evaluate(expr, m):
            tot += co.ring_area(a[1], b[1] if b is not None else None)
    return tot


def ring_value(table, mid):
    for (a, b), v in table.items():
        if mid is None:
            if b is None:
                return v
        elif b is not None and a < mid < b:
            return v
        elif b is None and mid > a:
            return v
    return False


def closed_forms():
    """{(chain, op): (volume, area, centre)}: the coaxial chains."""
    out = {}
    g9 = Coaxial(0, 0)
    g9.disc(0, 9, 9, 0, 3)
    g9.disc(1, 7, 6, 0, 4)
    g9.disc(2, 2, 2, -1, 5)
    sink = Coaxial(5, 5)
    sink.box(0, 10, 0, 4)
    sink.disc(1, 1, 3, 1, 4)
    sink.disc(2, mp.mpf(1)/2, mp.mpf(1)/2, -1, 6)
    for name, co, op1 in (('g9_bore', g9, 'fuse'), ('countersink_drill', sink, 'cut')):
        for op in OPS:
            V, c, A = co.measures((op, (op1, 0, 1), 2))
            out[(name, op)] = (V, A, c)
    return out


# ------------------------------------------------------------------ the chains' work

def evaluate(job):
    """One chain: the reference's rows for its last operations and every
    check's deviations (run in a worker process)."""
    name, specs, op1, stages, ops, last_swapped, degenerate, mc_n = job
    cases_ = [solid(s, OPERATIONS[k]) for k, s in enumerate(specs)]
    chain = ref.Chain(cases_)
    size = chain.size
    exprs = {}
    for op, solids in ops.items():
        g = Given(name+'_'+op, 'deep', specs, op1, stages+[(op, last_swapped)], 'solid', 1)
        exprs[op] = (g.expr(), solids)
    given = g.given_expr()
    rows, res = {}, {}
    for op, (expr, solids) in exprs.items():
        rows[op] = ref.rows(chain, expr, solids)
        res[op] = chain.result(expr, solids)
    checks = {}
    # Two families.
    fam = mp.mpf(0)
    for op, (expr, _) in exprs.items():
        a, b = chain.measures(expr, 0), chain.measures(expr, 1)
        fam = max(fam, abs(a[0]-b[0])/size**3, abs(a[2]-b[2])/size**2,
                  max(abs(a[1][i]-b[1][i]) for i in range(3))/size**4)
    checks['families'] = fam
    # Each solid's closed form.
    cdev = mp.mpf(0)
    for k, I in enumerate(chain.inputs):
        V, mom, A = I.closed()
        v, m, a = chain.measures(k)
        cdev = max(cdev, abs(v-V)/size**3, abs(a-A)/size**2, max(abs(m[i]-mom[i]) for i in range(3))/size**4)
    checks['solid_closed_forms'] = cdev
    # Pair identities for the given result X and the last solid C.
    last = len(specs)-1
    vx, mx, ax = chain.measures(given)
    vc, mc, ac = chain.measures(last)
    m = {op: chain.measures(e) for op, (e, _) in exprs.items()}
    dev = mp.mpf(0)
    if set(m) == set(OPS):
        vf, mf, af = m['fuse']
        vn, mn, an = m['common']
        vt, mt, at_ = m['cut']
        dev = abs(vf+vn-vx-vc)/size**3
        for i in range(3):
            dev = max(dev, abs(mf[i]+mn[i]-mx[i]-mc[i])/size**4)
        if exprs['cut'][0][1] == last:
            dev = max(dev, abs(vt-(vc-vn))/size**3, max(abs(mt[i]-(mc[i]-mn[i])) for i in range(3))/size**4)
        else:
            dev = max(dev, abs(vt-(vx-vn))/size**3, max(abs(mt[i]-(mx[i]-mn[i])) for i in range(3))/size**4)
        checks['pair_identities'] = dev
        # A face of C on another solid's surface, or another's on C's.
        shared = any((i == last and sw.skip) or any(sw.others.surfs[k][0] == last for k in sw.skip)
                     for i, _, sw in chain.sweeps(0))
        if not shared:
            checks['area_identity'] = abs(af+an-ax-ac)/size**2
    # The given result one solid.
    assert ref.count_solids(chain, given, 48) == 1, f'{name}: the given result is not one solid'
    # Solid counts at two resolutions.
    for op, (expr, solids) in exprs.items():
        for n in (48, 71):
            got = ref.count_solids(chain, expr, n)
            assert got == solids, f'{name} {op}: {got} solids by rays ({n}), {solids} declared'
    # Monte Carlo.
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
    sine, gap, spacing = chain.margins()
    margins = {'sine': sine, 'gap': gap, 'spacing': spacing}
    stats = {'quadrature': chain.quad_error()/size**4, 'missed': chain.missed(), 'refined': chain.refined()}
    return name, rows, res, checks, margins, stats, degenerate


def cached(job):
    """`evaluate`, its result kept in `--cache` (a development aid: keyed by
    the job and the reference's and this generator's sources)."""
    cache = os.environ.get('GIVEN_CURVED_CACHE')
    if not cache:
        return evaluate(job)
    import hashlib
    import inspect
    import pickle
    key = hashlib.sha256((repr(job)+inspect.getsource(evaluate)
                          +(Path(ref.__file__).read_text())).encode()).hexdigest()[:24]
    path = Path(cache)/f'{job[0]}-{key}.pickle'
    if path.exists():
        return pickle.loads(path.read_bytes())
    try:
        out = evaluate(job)
    except Exception as e:
        raise RuntimeError(f'{job[0]}: {e!r}') from e
    path.write_bytes(pickle.dumps(out))
    return out


def run(jobs, fn, workers):
    if workers <= 1:
        return [fn(j) for j in jobs]
    with ProcessPoolExecutor(workers) as pool:
        return list(pool.map(fn, jobs))


def reference_checks(results):
    worst, covered = {}, {}

    def note(key, value):
        worst[key] = max(worst.get(key, mp.mpf(0)), value)
        covered[key] = covered.get(key, 0)+1

    by_chain = {r[0]: r for r in results}
    for (name, op), (V, S, Cn) in closed_forms().items():
        if name not in by_chain:
            continue
        res = by_chain[name][2]
        n, vol, area, centre = res[op]
        size = max(abs(V)**(mp.mpf(1)/3), 1)
        note('closed_forms', max(abs(vol-V)/abs(V), abs(area-S)/abs(S),
                                 max(abs(centre[i]-Cn[i]) for i in range(3))/size))
    for name, rows, _, checks, margins, stats, degenerate in results:
        for key, value in checks.items():
            note(key, value)
        note('quadrature_estimate', stats['quadrature'])
        assert stats['missed'] == 0, f'{name}: quadrature nodes off their structure'
        if not degenerate:
            note('least_sine', 1/max(margins['sine'], mp.mpf(10)**-30))
            note('event_spacing', 1/max(margins['spacing'], mp.mpf(10)**-40))
    return worst, covered


def generate(results):
    by_chain = {r[0]: r for r in results}
    blocks = []
    out = ['# case\trow (S9e.3a, chained_curved_boolean_reference.py: expect KIND S9e.3a CLASS, reason TEXT for '
           'a degenerate case, then result N volume area cx cy cz or empty)']
    for case in cases():
        blocks.append(case.encode())
        rows = by_chain[case.chain_name][1]
        row = rows[case.last]
        n = 0 if row[0] == 'empty' else int(row[0].split()[1])
        want = {'solid': n > 0, 'empty': n == 0, 'degenerate': True}[case.kind]
        assert want, f'{case.name}: declared {case.kind}, the reference gives {row[0]}'
        out.append(f'{case.name}\texpect {case.kind} S9e.3a {case.klass}')
        if case.kind == 'degenerate':
            out.append(f'{case.name}\treason {case.reason}')
        for r in row:
            out.append(f'{case.name}\t{r}')
    frames = ['# case\tframe (solid k) axis (n, x, y)\tstored unit vector (stored_axes, as hex bits)']
    for case in cases():
        for k, c in enumerate(case.solid_cases):
            _, x, y, n = stored_axes(c.frame)
            for key, v in (('n', n), ('x', x), ('y', y)):
                frames.append(f'{case.name}\tsolid{k} {key}\t'+' '.join(struct.pack('>d', q).hex() for q in v))
    return {'boolean-given-curved-cases.txt': '\n'.join(blocks)+'\n',
            'boolean-given-curved-expected.tsv': '\n'.join(out)+'\n',
            'boolean-given-curved-frames.tsv': '\n'.join(frames)+'\n'}


def jobs(mc_n):
    chains = {}
    for c in cases():
        e = chains.setdefault(c.chain_name, [c.specs, c.op1, c.stages[:-1], {}, c.stages[-1][1],
                                             c.kind == 'degenerate'])
        e[3][c.last] = c.solids
    return [(name, specs, op1, stages, ops, sw, degenerate, mc_n)
            for name, (specs, op1, stages, ops, sw, degenerate) in chains.items()]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--workers', type=int, default=max(1, min(4, (os.cpu_count() or 2)-2)))
    parser.add_argument('--samples', type=int, default=200000, help='Monte-Carlo points per chain')
    parser.add_argument('--only', help='one chain, printed (no files)')
    parser.add_argument('--cache', help='a directory keeping each chain\'s result (development)')
    args = parser.parse_args()
    if args.cache:
        os.environ['GIVEN_CURVED_CACHE'] = args.cache
    listed = cases()
    validate(listed)
    todo = jobs(args.samples)
    if args.only:
        todo = [j for j in todo if j[0] == args.only]
        for r in run(todo, cached, 1):
            print(r)
        return
    results = run(todo, cached, args.workers)
    worst, covered = reference_checks(results)
    limits = {'closed_forms': 1e-30, 'families': 1e-30, 'solid_closed_forms': 1e-30, 'pair_identities': 1e-30,
              'area_identity': 1e-30, 'monte_carlo_sigma': 5, 'quadrature_estimate': 1e-30,
              'least_sine': 20, 'event_spacing': 1e6}
    for key, value in worst.items():
        assert value <= limits[key], (key, mp.nstr(value, 3))
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
