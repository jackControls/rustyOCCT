#!/usr/bin/env python3
"""Fixtures for S9e.4a of REVIEW_NOTES.md: imported solids (bodies without a
construction, read from `.brep` files OCCT wrote) given to Booleans, decided
on the construction their stored surfaces give.

The bodies are OCCT's own output: `boolean-imported-bodies.txt` holds one
`write NAME imported/NAME.brep` block per body for `occt_boolean_oracle.cpp`
(`BRepPrimAPI_MakeBox`, `MakeCylinder`, `MakeSphere`, `MakeCone`, `MakeTorus`,
or `MakePrism` of a profile face), which `compare_imported_boolean.py
--write-bodies` runs to write `rust/fixtures/imported/NAME.brep` (format
version 1, no triangulations). Nothing here reads the kernel; the files are
read only to check that each is the body it claims to be (below).

`boolean-imported-cases.txt` lists each case in the Boolean protocol, an
imported input as its one `brep imported/NAME.brep` row
(`identity_reference.Case.brep`), a chain's further `then` rows as S9e.1's;
`boolean-imported-expected.tsv` gives per case:

* `expect KIND S9e.4a CLASS`: the declared outcome (`solid`: one or more
  solids; `empty`; `degenerate`: `Degenerate` in the decisions;
  `unsupported`: `OutOfDomain`, S9e.4b's) and the imported body's class
  (`prism`: a box or a prism of lines and arcs, `cylinder`, `sphere`,
  `cone`, `torus`; `both`: two imported inputs; `chain`: an imported
  input's result given to another Boolean), then for a degenerate or
  unsupported case `reason TEXT`;
* `result N volume area cx cy cz` (totals over the N solids, world
  coordinates) or `empty`.

`boolean-imported-frames.tsv` records the stored axes of every solid the
kernel builds from its rows (not of the imported ones, whose frames are the
converter's from the files).

The reference is the construction OCCT was given, exactly as the other
references take a construction (`stored_axes`, the profile's binary64
data): S9e.3a's chained reference (`chained_curved_boolean_reference.py`)
for every pair and the chain whose solids it takes (prisms of one convex
profile, whole spheres, cones and frusta, whole tori); S9d.1's sphere
reference (`sphere_boolean_reference.py`, through
`generate_sphere_boolean_fixtures.evaluate`) for the hemisphere against a
polyhedral box; the coaxial closed forms (`generate_given_curved_boolean_
fixtures.Coaxial`) alone for the plate with a hole (the chained reference
takes no hole), with a Monte-Carlo estimate of its own. The imported
bodies' stored data are OCCT's roundings of these constructions (vertices
to 15 significant digits), so the kernel's model of them is within the
resolution of the reference's, not equal to it.

The bodies (`B` a box): `box`, `MakeBox` of `[0, 10]^2 x [0, 4]`;
`box_tilt`, `MakeBox` of `6 x 5 x 4` at `(2, 1, 0)` in the `TILT` frame
(its corners rounded); `cyl`, `MakeCylinder` of radius 3 about `(5, 5)`
over `z` in `[-1, 5]`; `cyl_side`, of radius 3/2 along `x` through `(y, z)
= (5, 2)`; `dee`, `MakePrism` of a profile of lines and an arc tangent to
them (`(0, -5)` to `(10, -5)`, the arc about `(10, 0)` to `(10, 5)`, back
to `(0, 5)`), shifted to `(0, 5)`, height 6; `halves`, a circle of radius
5/2 about `(5, 5)` as two arcs (two cylinder faces, a periodic face split
at a seam), height 5; `plate`, `[-6, 6]^2 x [0, 3]` with a hole of radius
2; `ball`, `MakeSphere` of radius 3 about `(5, 5, 4)`; `dome`, the
hemisphere of radius 4 above the origin (a cap: one disc and a pole);
`frustum`, `MakeCone` of radii 3 and 1 and height 4 at `(5, 5, 0)`;
`spike`, of radius 2 and height 5 with its apex; `ring`, `MakeTorus` of
radii 5 and 3/2; `dee_turn`, the `dee` profile in a frame turned by 30
degrees about `z`.

Cases (each the three operations): `box_slab`, `B` and a `TILT` slab;
`box_drill`, `B` and a coaxial rod of radius 2 (closed forms);
`tilt_box`, the tilted box and an axis-aligned slab; `cyl_box`, the
cylinder and `B` (closed forms); `box_cyl`, `B` (built) and the imported
cylinder as the tool (closed forms); `rods`, the side cylinder and a
vertical rod of radius 2 (their meetings two cylinders' curves);
`dee_slab`, the tangent profile and a `TILT` slab; `halves_box`, the
seam-split cylinder and a box across its seam's edges; `plate_bore`, the
plate and a coaxial rod of radius 3 (closed forms); `ball_box`, the ball
and `B` (hemispheres in closed form); `dome_box`, the hemisphere and a box
in the `TILT` frame; `frustum_bore`, the frustum and a coaxial rod of
radius 2 (closed forms); `frustum_box`, the frustum and a box whose wall
cuts it off its axis (hyperbolas); `spike_ball`, the cone with its apex and
a sphere of radius 3/2 about the apex; `ring_box`, the torus and a box whose bottom cuts
its tube; `ring_pin`, the torus and a rod through its tube; `both`, the
imported box and the imported ball (closed forms); `slab_ball`, a `TILT`
slab and the imported ball as the tool; `chain_drill`, the imported box
less a rod, then with a `TILT` slab. Declared `degenerate`: `cyl_tangent`,
the cylinder and a box whose wall touches it along a generatrix;
`ball_touch`, the ball and a box whose bottom touches its top point;
`box_kiss`, `B` and a rod touching its wall. Declared `unsupported`:
`dee_turn`, the turned profile (its arc's ends rounded off its circle in
its cap's frame) and a box.

Before writing, the checks (limits relative to the case's size): the
coaxial and hemisphere closed forms within 1e-30; the chained reference's
two families, each solid's closed form and the pair identities (`V(A u B) +
V(A n B) = V(A) + V(B)`, `V(A - B) = V(A) - V(A n B)`, moments likewise,
`area(A u B) + area(A n B) = area(A) + area(B)` where no faces share a
surface) within 1e-30; the sphere reference's own checks (its slicing two
ways and against the closed forms, inclusion and exclusion, the face
classes); Monte Carlo within 5 standard errors; the declared solid counts
by rays at two resolutions (or by the sphere reference's adjacency); every
scanned meeting at a sine of at least 0.05 and events at least 1e-6 of
their range apart outside the declared cases (an event found twice, within
1e-25, counted once: a tangent profile's joint, a coaxial rod's circle).
With `--check`, where the bodies' files exist, each is read
(`stored_records`: `brep_io_reference`'s records and locations, every face
and vertex instance as placed): its faces' surfaces' kinds are the
construction's faces' and every stored vertex lies within 1e-12 of the
case's size on the construction's surfaces.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
import math
import os
from pathlib import Path
import struct
import types
import zlib

import mpmath as mp

from identity_reference import Boundary, Case, encode_case, native_case, BOOLEAN_OPERATIONS
from curve_surface_reference import stored_axes
import chained_curved_boolean_reference as ref
from generate_curved_boolean_fixtures import FRAMES as CURVED_FRAMES, square, disc
from generate_given_curved_boolean_fixtures import Coaxial
from torus_boolean_reference import TWO_PI, number
import torus_curved_boolean_reference as tc
import brep_io_reference as bio

ROOT = Path(__file__).resolve().parents[1]
OPERATIONS = (91, 92, 94)          # the solids' construction (import) ids
BOOLEANS = (93, 95)                # the Booleans' ids
CLASSES = ('prism', 'cylinder', 'sphere', 'cone', 'torus', 'both', 'chain')
HP = 1.5707963267948966
OPS = ('fuse', 'cut', 'common')
FRAMES = dict(CURVED_FRAMES, TURN30=(0.0, 0.0, 1.0, 0.8660254037844386, 0.5, 0.0))
BODIES_FILE = 'boolean-imported-bodies.txt'


def at(name, origin):
    return tuple(float(c) for c in origin)+FRAMES[name]


# ------------------------------------------------------------------ solids

def torus(R, r, frame):
    return ('torus', R, r, frame)


def sphere(r, frame, low=-HP, high=HP):
    return ('sphere', r, frame, low, high)


def cone(r0, r1, h, frame):
    return ('cone', r0, r1, h, frame)


def prism(boundaries, frame, start, end):
    return ('prism', boundaries, frame, start, end)


def construction(spec, op):
    """The identity case of a construction spec."""
    kind = spec[0]
    if kind == 'torus':
        _, R, r, frame = spec
        return Case('', 1e-7, op, frame, 0.0, 0.0, [], torus=(float(R), float(r), 0.0, TWO_PI, TWO_PI))
    if kind == 'sphere':
        _, r, frame, low, high = spec
        return Case('', 1e-7, op, frame, 0.0, 0.0, [], sphere=(float(r), low, high))
    if kind == 'cone':
        _, r0, r1, h, frame = spec
        return Case('', 1e-7, op, frame, 0.0, 0.0, [], cone=(float(r0), float(r1), float(h)))
    _, boundaries, frame, start, end = spec
    return Case('', 1e-7, op, frame, float(start), float(end), boundaries)


class Body:
    """An imported body: its construction (the reference's) and the rows
    OCCT builds it from (`primitive`: a `box` or `cylinder` row, else the
    construction's own `native_case` rows)."""

    def __init__(self, name, spec, klass, primitive=None):
        self.name, self.spec, self.klass, self.primitive = name, spec, klass, primitive
        self.path = f'imported/{name}.brep'

    def native_rows(self):
        if self.primitive is not None:
            return [self.primitive]
        return native_case(construction(self.spec, 0)).split('\n')[1:-1]


def numbers(*values):
    return ' '.join(repr(float(v)) for v in values)


def box_body(name, frame, size):
    """`BRepPrimAPI_MakeBox(gp_Ax2(origin, normal, x), DX, DY, DZ)`: the prism
    of `[0, DX] x [0, DY]` on the frame over `[0, DZ]`."""
    o, n, x = frame[:3], frame[3:6], frame[6:9]
    return Body(name, prism([square(0.0, 0.0, size[0], size[1])], frame, 0.0, size[2]), 'prism',
                'box '+numbers(*o, *n, *x, *size))


def cylinder_body(name, frame, r, h):
    """`BRepPrimAPI_MakeCylinder(gp_Ax2(origin, normal, x), R, H)`: the prism
    of the circle of radius `R` about the frame's origin over `[0, H]`."""
    o, n, x = frame[:3], frame[3:6], frame[6:9]
    return Body(name, prism([disc(0.0, 0.0, r)], frame, 0.0, h), 'cylinder',
                'cylinder '+numbers(*o, *n, *x, r, h))


def dee(x0, y0):
    """A profile of lines and an arc tangent to them: `(x0, y0 - 5)` to
    `(x0 + 10, y0 - 5)`, the arc about `(x0 + 10, y0)` to `(x0 + 10, y0 + 5)`,
    back to `(x0, y0 + 5)`."""
    return Boundary(points=[(x0, y0-5.0), (x0+10.0, y0-5.0), (x0+10.0, y0+5.0), (x0, y0+5.0)],
                    segments=[None, (x0+10.0, y0, 5.0, True), None, None])


def halves(cx, cy, r):
    """A circle as two arcs (two faces on one cylinder)."""
    return Boundary(points=[(cx+r, cy), (cx-r, cy)], segments=[(cx, cy, r, True), (cx, cy, r, True)])


BODIES = {b.name: b for b in [
    box_body('box', at('XY', (0, 0, 0)), (10.0, 10.0, 4.0)),
    box_body('box_tilt', at('TILT', (2, 1, 0)), (6.0, 5.0, 4.0)),
    cylinder_body('cyl', at('XY', (5, 5, -1)), 3.0, 6.0),
    cylinder_body('cyl_side', at('SIDE', (-1, 5, 2)), 1.5, 12.0),
    Body('dee', prism([dee(0.0, 5.0)], at('XY', (0, 0, 0)), 0.0, 6.0), 'prism'),
    Body('halves', prism([halves(5.0, 5.0, 2.5)], at('XY', (0, 0, 0)), 0.0, 5.0), 'cylinder'),
    Body('plate', prism([square(-6.0, -6.0, 6.0, 6.0), disc(0.0, 0.0, 2.0)], at('XY', (0, 0, 0)), 0.0, 3.0),
         'prism'),
    Body('ball', sphere(3.0, at('XY', (5, 5, 4))), 'sphere'),
    Body('dome', sphere(4.0, at('XY', (0, 0, 0)), 0.0, HP), 'sphere'),
    Body('frustum', cone(3.0, 1.0, 4.0, at('XY', (5, 5, 0))), 'cone'),
    Body('spike', cone(2.0, 0.0, 5.0, at('XY', (0, 0, 0))), 'cone'),
    Body('ring', torus(5.0, 1.5, at('XY', (0, 0, 0))), 'torus'),
    Body('dee_turn', prism([dee(0.0, 5.0)], at('TURN30', (0, 0, 0)), 0.0, 6.0), 'prism'),
]}


def imported(name):
    return ('imported', name)


def resolve(spec):
    """A spec's construction: an imported body's, or the spec itself."""
    return BODIES[spec[1]].spec if spec[0] == 'imported' else spec


B = prism([square(0.0, 0.0, 10.0, 10.0)], at('XY', (0, 0, 0)), 0.0, 4.0)
SLAB = prism([square(-2.0, -15.0, 12.0, 15.0)], at('TILT', (0, 0, 1.5)), 0.0, 1.5)
# A `TILT` slab across the ball's middle.
BALL_SLAB = prism([square(-2.0, -15.0, 12.0, 15.0)], at('TILT', (0, 0, 7.5)), 0.0, 1.5)


class Imported:
    """One case: its solids (`imported(NAME)` or a construction), the first
    Boolean, the fixed stages and the last stage's operation."""

    def __init__(self, name, klass, specs, op1, stages, kind, solids, reason=None, reference='chain'):
        assert klass in CLASSES, klass
        self.name, self.klass, self.op1, self.stages = name, klass, op1, stages
        self.kind, self.solids, self.reason, self.reference = kind, solids, reason, reference
        self.group = name.rsplit('_', 1)[0]
        self.specs = specs
        self.constructions = [construction(resolve(s), OPERATIONS[k]) for k, s in enumerate(specs)]
        self.solid_cases = []
        for k, s in enumerate(specs):
            c = construction(resolve(s), OPERATIONS[k])
            if s[0] == 'imported':
                c.brep = BODIES[s[1]].path
            c.name = name
            self.solid_cases.append(c)
        for c in self.constructions:
            c.name = name

    @property
    def last(self):
        return self.stages[-1][0] if self.stages else self.op1

    def expr(self):
        if not self.stages:
            return (self.op1, 0, 1)
        (op2, sw2), more = self.stages[0], self.stages[1:]
        return ref.chain_expr(self.op1, op2, sw2, more)

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


def group(name, klass, specs, outcomes, first=None, reason=None, kind=None, reference='chain'):
    """Cases of one group: the last operation varied over `outcomes` (`{op:
    solids}`, 0 for empty): a pair's Boolean, or with `first` a chain's
    second Boolean on the first's result."""
    out = []
    for op, n in outcomes.items():
        k = kind or ('degenerate' if reason else ('solid' if n else 'empty'))
        if first:
            out.append(Imported(f'{name}_{op}', klass, specs, first, [(op, False)], k, n, reason, reference))
        else:
            out.append(Imported(f'{name}_{op}', klass, specs, op, [], k, n, reason, reference))
    return out


TANGENT = 'a box\'s wall tangent to the imported cylinder along a generatrix'
TOUCH = 'a box\'s bottom face touching the imported sphere at its top point'
KISS = 'a rod tangent to the imported box\'s wall along a generatrix'
TURNED = ('S9e.4b: the imported prism\'s arc whose ends, rounded into its cap\'s frame turned by 30 degrees, '
          'lie off its circle')


def cases():
    out = []
    three = {'fuse': 1, 'cut': 1, 'common': 1}
    out += group('box_slab', 'prism', [imported('box'), SLAB], {'fuse': 1, 'cut': 2, 'common': 1})
    out += group('box_drill', 'prism', [imported('box'), prism([disc(5.0, 5.0, 2.0)], at('XY', (0, 0, -1)), 0.0, 6.0)],
                 three)
    out += group('tilt_box', 'prism', [imported('box_tilt'), prism([square(0.0, 0.0, 10.0, 10.0)],
                                                                    at('XY', (0, 0, 1)), 0.0, 1.5)],
                 {'fuse': 1, 'cut': 2, 'common': 1})
    out += group('cyl_box', 'cylinder', [imported('cyl'), B], {'fuse': 1, 'cut': 2, 'common': 1})
    out += group('box_cyl', 'cylinder', [B, imported('cyl')], three)
    out += group('rods', 'cylinder', [imported('cyl_side'), prism([disc(5.0, 5.0, 2.0)], at('XY', (0, 0, -1)),
                                                                  0.0, 6.0)],
                 {'fuse': 1, 'cut': 2, 'common': 1})
    out += group('dee_slab', 'prism', [imported('dee'), SLAB], {'fuse': 1, 'cut': 1, 'common': 1})
    out += group('halves_box', 'cylinder', [imported('halves'), prism([square(5.5, 2.0, 9.0, 8.0)],
                                                                      at('XY', (0, 0, 1)), 0.0, 3.0)],
                 three)
    out += group('plate_bore', 'prism', [imported('plate'), prism([disc(0.0, 0.0, 3.0)], at('XY', (0, 0, -1)),
                                                                  0.0, 5.0)],
                 three, reference='coaxial')
    out += group('ball_box', 'sphere', [imported('ball'), B], three)
    out += group('dome_box', 'sphere', [imported('dome'), prism([square(-2.0, -6.0, 7.0, 6.0)],
                                                                at('TILT', (0, 0, 1)), 0.0, 6.0)],
                 three, reference='sphere')
    out += group('frustum_bore', 'cone', [imported('frustum'), prism([disc(5.0, 5.0, 2.0)], at('XY', (0, 0, -1)),
                                                                     0.0, 6.0)],
                 three)
    out += group('frustum_box', 'cone', [imported('frustum'), prism([square(6.0, 0.0, 12.0, 10.0)],
                                                                    at('XY', (0, 0, 1)), 0.0, 2.0)],
                 three)
    out += group('spike_ball', 'cone', [imported('spike'), sphere(1.5, at('XY', (0, 0, 5)))],
                 three)
    out += group('ring_box', 'torus', [imported('ring'), prism([square(-8.0, -8.0, 8.0, 8.0)], at('XY', (0, 0, 0.5)),
                                                               0.0, 3.5)],
                 three)
    out += group('ring_pin', 'torus', [imported('ring'), prism([disc(5.0, 0.0, 0.75)], at('XY', (0, 0, -3)),
                                                               0.0, 6.0)],
                 {'fuse': 1, 'cut': 1, 'common': 1})
    out += group('both', 'both', [imported('box'), imported('ball')], three)
    out += group('slab_ball', 'sphere', [BALL_SLAB, imported('ball')], three)
    out += group('chain_drill', 'chain', [imported('box'), prism([disc(5.0, 5.0, 2.0)], at('XY', (0, 0, -1)),
                                                                 0.0, 6.0), SLAB],
                 {'fuse': 1, 'cut': 2, 'common': 1}, first='cut')
    out += group('cyl_tangent', 'cylinder', [imported('cyl'), prism([square(8.0, 0.0, 12.0, 10.0)],
                                                                    at('XY', (0, 0, 0)), 0.0, 4.0)],
                 {'fuse': 1, 'cut': 1, 'common': 0}, reason=TANGENT)
    out += group('ball_touch', 'sphere', [imported('ball'), prism([square(3.0, 3.0, 7.0, 7.0)],
                                                                  at('XY', (0, 0, 7)), 0.0, 2.0)],
                 {'fuse': 2, 'cut': 1, 'common': 0}, reason=TOUCH)
    out += group('box_kiss', 'prism', [imported('box'), prism([disc(12.0, 5.0, 2.0)], at('XY', (0, 0, -1)),
                                                              0.0, 6.0)],
                 {'fuse': 1, 'cut': 1, 'common': 0}, reason=KISS)
    out += group('dee_turn', 'prism', [imported('dee_turn'), B], three, kind='unsupported', reason=TURNED)
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
            frame = resolve(s)[{'torus': 3, 'sphere': 2, 'cone': 4, 'prism': 2}[resolve(s)[0]]]
            assert tuple(frame[3:]) in FRAMES.values(), c.name


# ------------------------------------------------------------------ closed forms

def coaxial_forms():
    """{(group, op): (volume, area, centre)} by the coaxial sections, the
    expressions over the forms' own indices."""
    out = {}

    def add(group, co, expr_of):
        for op in OPS:
            V, c, A = co.measures(expr_of(op))
            out[(group, op)] = (V, A, c)
    drill = Coaxial(5, 5)
    drill.box(0, 10, 0, 4)
    drill.disc(1, 2, 2, -1, 5)
    add('box_drill', drill, lambda op: (op, 0, 1))
    cyl = Coaxial(5, 5)
    cyl.disc(0, 3, 3, -1, 5)
    cyl.box(1, 10, 0, 4)
    add('cyl_box', cyl, lambda op: (op, 0, 1))
    add('box_cyl', cyl, lambda op: (op, 1, 0))
    plate = Coaxial(0, 0)
    plate.box(0, 12, 0, 3)
    plate.disc(1, 2, 2, 0, 3)
    plate.disc(2, 3, 3, -1, 4)
    add('plate_bore', plate, lambda op: (op, ('cut', 0, 1), 2))
    cone_ = Coaxial(5, 5)
    cone_.disc(0, 3, 1, 0, 4)
    cone_.disc(1, 2, 2, -1, 5)
    add('frustum_bore', cone_, lambda op: (op, 0, 1))
    return out


def hemisphere_forms():
    """{(group, op): (volume, area, centre)}: the ball of radius 3 about
    `(5, 5, 4)` against `[0, 10]^2 x [0, 4]`, whose top face holds the ball's
    centre (hemispheres; the box's volume 400 and area 360)."""
    pi = mp.pi
    half = 18*pi                            # (2/3) pi 3^3
    zc = mp.mpf(9)/8                        # a hemisphere's centroid, 3 r / 8
    out = {}
    five = mp.mpf(5)
    for group in ('ball_box', 'both'):
        out[(group, 'common')] = (half, 27*pi, (five, five, 4-zc))
        V = 400+half
        out[(group, 'fuse')] = (V, 360+9*pi, (five, five, (400*2+half*(4+zc))/V))
    # The ball less the box (`ball_box`: the upper hemisphere), the box less
    # the ball (`both`, the box first: the lower hemisphere's hollow).
    out[('ball_box', 'cut')] = (half, 27*pi, (five, five, 4+zc))
    V = 400-half
    out[('both', 'cut')] = (V, 360+9*pi, (five, five, (400*2-half*(4-zc))/V))
    return out


# ------------------------------------------------------------------ the groups' work

def evaluate_chain(job):
    """One group on the chained reference: rows, results and checks."""
    name, specs, op1, stages, ops, degenerate, mc_n = job
    cases_ = [construction(resolve(s), OPERATIONS[k]) for k, s in enumerate(specs)]
    chain = ref.Chain(cases_)
    size = chain.size
    exprs = {}
    for op, solids in ops.items():
        if stages:
            exprs[op] = (ref.chain_expr(op1, op, False), solids)
        else:
            exprs[op] = ((op, 0, 1), solids)
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
    # Pair identities: X the first argument of the last Boolean (the first
    # result in a chain), C its partner.
    last = len(specs)-1
    x_expr = (op1, 0, 1) if stages else 0
    c_expr = last if stages else 1
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
        shared = any(sw.skip for _, _, sw in chain.sweeps(0))
        if not shared:
            checks['area_identity'] = abs(af+an-ax-ac)/size**2
    if stages:
        assert ref.count_solids(chain, x_expr, 48) == 1, f'{name}: the given result is not one solid'
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
    margins = {'sine': sine, 'gap': gap, 'spacing': spacing(chain)}
    stats = {'quadrature': chain.quad_error()/size**4, 'missed': chain.missed(), 'refined': chain.refined()}
    return name, rows, res, checks, margins, stats, degenerate


def spacing(chain):
    """`Chain.margins`'s events' least spacing relative to their family's
    range, an event found twice (within 1e-25 of the range: one event from
    two sources, as a tangent profile's joint is both a vertex and its arc's
    extreme along a chord, or a coaxial rod's circle on a frustum both
    surfaces' event) counted once."""
    least = mp.mpf(1)
    for k in chain._sweeps:
        for _, _, sw in chain._sweeps[k]:
            span = sw.fam.a1-sw.fam.a0
            pts = sorted(set(sw.events))
            for a, b in zip(pts, pts[1:]):
                if (b-a)/span > mp.mpf(10)**-25:
                    least = min(least, (b-a)/span)
    return least


def evaluate_sphere(job):
    """The hemisphere's group on S9d.1's sphere reference, with its checks."""
    import generate_sphere_boolean_fixtures as sg
    name, specs, op1, stages, ops, degenerate, mc_n = job
    assert not stages
    obj, tool = (construction(resolve(s), OPERATIONS[k]) for k, s in enumerate(specs))
    first = types.SimpleNamespace(obj=obj, tool=tool, coplanar=False, pair_name=name)
    _, rows, res, checks, near, stats = sg.evaluate((name, first, tuple(ops), mc_n))
    assert not near, (name, near)
    out = {}
    for op, solids in ops.items():
        n = res[op][0]
        assert n == solids, f'{name} {op}: {n} solids, {solids} declared'
        out[op] = res[op]
    checks = {k: v for k, v in checks.items()}
    stats = {'quadrature': stats['volume_quadrature'], 'missed': 0, 'refined': 0}
    margins = {'sine': mp.mpf(1), 'gap': mp.mpf(1), 'spacing': mp.mpf(1)}
    return name, {op: rows[op] for op in ops}, out, checks, margins, stats, degenerate


def plate_contains(X):
    x, y, z = X
    return abs(x) < 6 and abs(y) < 6 and 0 < z < 3 and x*x+y*y > 4


def bore_contains(X):
    x, y, z = X
    return x*x+y*y < 9 and -1 < z < 4


def evaluate_coaxial(job):
    """The plate's group: the coaxial closed forms and a Monte-Carlo estimate
    of each operation (the chained reference takes no hole)."""
    import random
    name, specs, op1, stages, ops, degenerate, mc_n = job
    forms = coaxial_forms()
    rows, res = {}, {}
    for op, solids in ops.items():
        V, A, c = forms[(name, op)]
        res[op] = (solids, V, A, c)
        rows[op] = ['result {} {} {} {} {} {}'.format(solids, *(number(v) for v in (V, A, *c)))]
    rnd = random.Random(zlib.crc32(name.encode()))
    lo, hi = (-6.5, -6.5, -1.5), (6.5, 6.5, 4.5)
    vol = 13.0*13.0*6.0
    sets = {'fuse': lambda a, b: a or b, 'cut': lambda a, b: a and not b, 'common': lambda a, b: a and b}
    acc = {op: [0, [0.0]*3, [0.0]*3] for op in ops}
    for _ in range(mc_n):
        X = [lo[i]+(hi[i]-lo[i])*rnd.random() for i in range(3)]
        a, b = plate_contains(X), bore_contains(X)
        for op in ops:
            if sets[op](a, b):
                e = acc[op]
                e[0] += 1
                for i in range(3):
                    e[1][i] += X[i]
                    e[2][i] += X[i]*X[i]
    z = 0.0
    for op, (k, s, s2) in acc.items():
        p = k/mc_n
        sv = vol*math.sqrt(max(p*(1-p), 0.0)/mc_n)
        V, A, c = forms[(name, op)]
        z = max(z, abs(vol*p-float(V))/sv)
        cm = [s[i]/k for i in range(3)]
        sc = [math.sqrt(max(s2[i]/k-cm[i]*cm[i], 0.0)/k) for i in range(3)]
        z = max(z, max(abs(cm[i]-float(c[i]))/max(sc[i], 1e-12) for i in range(3)))
    checks = {'monte_carlo_sigma': mp.mpf(z)}
    # Inclusion and exclusion on the closed forms.
    (vf, af, cf), (vt, at_, ct), (vn, an, cn) = (forms[(name, op)] for op in OPS)
    stats = {'quadrature': mp.mpf(0), 'missed': 0, 'refined': 0}
    margins = {'sine': mp.mpf(1), 'gap': mp.mpf(1), 'spacing': mp.mpf(1)}
    vp = 12*12*3-mp.pi*4*3
    vb = mp.pi*9*5
    checks['pair_identities'] = max(abs(vf+vn-vp-vb), abs(vt-(vp-vn)))/mp.mpf(13)**3
    return name, rows, res, checks, margins, stats, degenerate


def evaluate(job):
    kind = job[-1]
    job = job[:-1]
    return {'chain': evaluate_chain, 'sphere': evaluate_sphere, 'coaxial': evaluate_coaxial}[kind](job)


def cached(job):
    """`evaluate`, its result kept in `--cache` (a development aid: keyed by
    the job and the references' and this generator's sources)."""
    cache = os.environ.get('IMPORTED_CACHE')
    if not cache:
        return evaluate(job)
    import hashlib
    import inspect
    import pickle
    key = hashlib.sha256((repr(job)+inspect.getsource(evaluate_chain)+inspect.getsource(evaluate_sphere)
                          +inspect.getsource(evaluate_coaxial)+Path(ref.__file__).read_text()).encode()).hexdigest()[:24]
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


def jobs(mc_n):
    groups = {}
    for c in all_cases():
        e = groups.setdefault(c.group, [c.specs, c.op1 if c.stages else None, bool(c.stages), {},
                                        c.kind in ('degenerate', 'unsupported'), c.reference])
        e[3][c.last] = c.solids
    out = []
    for name, (specs, op1, stages, ops, degenerate, reference) in groups.items():
        out.append((name, specs, op1, stages, ops, degenerate, mc_n, reference))
    return out


def reference_checks(results):
    worst, covered = {}, {}

    def note(key, value):
        worst[key] = max(worst.get(key, mp.mpf(0)), value)
        covered[key] = covered.get(key, 0)+1

    by_group = {r[0]: r for r in results}
    forms = dict(coaxial_forms())
    forms.update(hemisphere_forms())
    for (name, op), (V, S, Cn) in forms.items():
        if name not in by_group:
            continue
        n, vol, area, centre = by_group[name][2][op]
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


LIMITS = {'closed_forms': 1e-30, 'families': 1e-30, 'solid_closed_forms': 1e-30, 'pair_identities': 1e-30,
          'area_identity': 1e-30, 'monte_carlo_sigma': 5, 'quadrature_estimate': 1e-30, 'least_sine': 20,
          'event_spacing': 1e6,
          # S9d.1's sphere reference's own checks (`generate_sphere_boolean_fixtures.LIMITS`).
          'inputs_sliced': 1e-30, 'common_two_ways': 1e-30, 'inclusion_exclusion': 1e-30, 'face_classes': 1e-30,
          'shared_both_sides': 1e-30, 'input_areas': 1e-30, 'second_direction': 1e-30}


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
    out = ['# case\trow (S9e.4a; chained_curved_boolean_reference.py, sphere_boolean_reference.py or the coaxial '
           'closed forms: expect KIND S9e.4a CLASS, reason TEXT for a degenerate or unsupported case, then result '
           'N volume area cx cy cz or empty)']
    for case in all_cases():
        blocks.append(case.encode())
        row = by_group[case.group][1][case.last]
        n = 0 if row[0] == 'empty' else int(row[0].split()[1])
        want = {'solid': n > 0, 'empty': n == 0, 'degenerate': True, 'unsupported': True}[case.kind]
        assert want, f'{case.name}: declared {case.kind}, the reference gives {row[0]}'
        out.append(f'{case.name}\texpect {case.kind} S9e.4a {case.klass}')
        if case.kind in ('degenerate', 'unsupported'):
            out.append(f'{case.name}\treason {case.reason}')
        for r in row:
            out.append(f'{case.name}\t{r}')
    frames = ['# case\tframe (solid k) axis (n, x, y)\tstored unit vector (stored_axes, as hex bits; the solids '
              'the kernel builds from rows, not the imported ones)']
    for case in all_cases():
        for k, (spec, c) in enumerate(zip(case.specs, case.constructions)):
            if spec[0] == 'imported':
                continue
            _, x, y, n = stored_axes(c.frame)
            for key, v in (('n', n), ('x', x), ('y', y)):
                frames.append(f'{case.name}\tsolid{k} {key}\t'+' '.join(struct.pack('>d', q).hex() for q in v))
    return {'boolean-imported-cases.txt': '\n'.join(blocks)+'\n',
            'boolean-imported-expected.tsv': '\n'.join(out)+'\n',
            'boolean-imported-frames.tsv': '\n'.join(frames)+'\n',
            BODIES_FILE: bodies_text()}


def native_input():
    """Every case's native rows (`compare_imported_boolean.py`)."""
    return '\n'.join(c.native() for c in all_cases())+'\n'


# ------------------------------------------------------------------ the files' check

def stored_records(text):
    """The faces' surface kinds and the vertices' points of a `.brep` OCCT
    wrote, as placed: `brep_io_reference.read`'s records and locations (a
    prism's top shares its bottom's records under a translation, as
    `MakePrism` without a copy builds it), every face and vertex instance
    reached from the root once (by record and composed location), each
    vertex record's point (its tolerance's next three numbers, the records'
    order the reader's) under its instance's location; the kinds
    `brep_io_reference`'s names (`plane`, `cylinder`, `cone`, `sphere`,
    `torus`)."""
    locations, tables, shapes, root = bio.read(text)
    words = text.split()
    at = words.index('TShapes')
    points = iter(tuple(float(w) for w in words[j+2:j+5]) for j in range(at, len(words)) if words[j] == 'Ve')
    point_of = {i: next(points) for i, (kind, _, _) in enumerate(shapes) if kind == 'Ve'}
    loc = lambda i: bio.IDENTITY if i == 0 else locations[i-1]
    faces, vertices, seen = [], [], set()

    def key(index, t):
        return (index, tuple(round(x, 12) for x in t))

    def walk(index, t):
        kind, data, subs = shapes[index]
        if key(index, t) in seen:
            return
        seen.add(key(index, t))
        if kind == 'Fa':
            faces.append(bio.kind_of(tables['Surfaces'][data[0]-1]))
        if kind == 'Ve':
            x = point_of[index]
            vertices.append(tuple(t[4*i]*x[0]+t[4*i+1]*x[1]+t[4*i+2]*x[2]+t[4*i+3] for i in range(3)))
            return
        for _, sub, l in subs:
            walk(sub, bio.matmul(t, loc(l)))
    _, index, l = root
    walk(index, loc(l))
    return faces, vertices


def construction_surfaces(spec):
    """The construction's faces' kinds and its surfaces' functions
    (`(f, grad)` in the world, mpf): a prism's per boundary (each a prism of
    one boundary's own), a sphere's, cone's or torus's, and a sphere's or a
    cone's end planes."""
    kind = spec[0]
    case = construction(spec, 0)
    if kind == 'prism':
        kinds, surfs = ['plane', 'plane'], []
        for b in spec[1]:
            one = construction(prism([b], spec[2], spec[3], spec[4]), 0)
            I = ref.make_input(one)
            surfs += I.surfs
            pts = b.circle is None and b.points
            if b.circle is not None:
                kinds.append('cylinder')
            else:
                kinds += ['plane' if s is None else 'cylinder' for s in b.segments or [None]*len(pts)]
        return kinds, surfs
    fr = tc.Frame(case.frame)
    if kind == 'sphere':
        _, r, frame, low, high = spec
        whole = construction(sphere(r, frame), 0)
        surfs = list(ref.make_input(whole).surfs)
        kinds = ['sphere']
        for lat in (low, high):
            if abs(lat) != HP:
                h = mp.mpf(r)*mp.sin(mp.mpf(lat))
                surfs.append(tc.plane_surf('end', fr, (tc.Z, tc.Z, mp.mpf(1)), h))
                kinds.append('plane')
        return kinds, surfs
    if kind == 'cone':
        _, r0, r1, h, _ = spec
        surfs = list(ref.make_input(case).surfs)
        kinds = ['cone']+['plane']*sum(1 for r in (r0, r1) if r > 0)
        surfs.append(tc.plane_surf('bottom', fr, (tc.Z, tc.Z, mp.mpf(1)), 0))
        surfs.append(tc.plane_surf('top', fr, (tc.Z, tc.Z, mp.mpf(1)), h))
        return kinds, surfs
    return ['torus'], list(ref.make_input(case).surfs)


def check_bodies():
    """Each existing body file: its stored surfaces' kinds the
    construction's faces' and each stored vertex on the construction's
    surfaces within 1e-12 of the body's size. Returns the largest
    deviation and the files read."""
    worst, read = mp.mpf(0), 0
    for b in BODIES.values():
        path = ROOT/'fixtures'/b.path
        if not path.exists():
            continue
        read += 1
        kinds, vertices = stored_records(path.read_text())
        want, surfs = construction_surfaces(b.spec)
        # Every face of the construction, by its surface's kind.
        assert sorted(kinds) == sorted(want), (b.name, kinds, want)
        I = construction(b.spec, 0)
        size = max([1.0]+[abs(v) for v in I.frame[:3]]+[abs(x) for v in vertices for x in v])
        for v in vertices:
            X = tuple(mp.mpf(c) for c in v)
            d = min(abs(s.f(X))/max(tc.norm(s.grad(X)), mp.mpf(10)**-30) for s in surfs)
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
        os.environ['IMPORTED_CACHE'] = args.cache
    listed = all_cases()
    validate(listed)
    todo = jobs(args.samples)
    if args.only:
        todo = [j for j in todo if j[0] == args.only]
        for r in run(todo, cached, 1):
            print(r)
        return
    results = run(todo, cached, args.workers)
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
