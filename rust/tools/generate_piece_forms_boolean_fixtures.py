#!/usr/bin/env python3
"""Fixtures for S9e.4b.3c.3a of REVIEW_NOTES.md: imported bodies of one
sphere, cylinder or cone face and plane faces that are another Boolean of
their primitive and the convex hull of their planes than S9e.4b.3a's common
(bodies OCCT wrote to `.brep` files) given to Booleans.

The bodies are OCCT's own output, as S9e.4b.3a's
(`generate_imported_pieces_boolean_fixtures.py`): `boolean-piece-forms-
bodies.txt` holds one `write NAME imported/NAME.brep` block per body for
`occt_boolean_oracle.cpp` (one solid's rows, a `boolean` row and another
solid's rows: the Boolean's one solid written), which
`compare_piece_forms_boolean.py --write-bodies` runs to write
`rust/fixtures/imported/NAME.brep` (format version 1, no triangulations).
Nothing here reads the kernel; the files are read only to check that each is
the body it claims to be (below).

`boolean-piece-forms-cases.txt` lists each case in the Boolean protocol, an
imported input as its one `brep imported/NAME.brep` row, a chain's further
`then` rows as S9e.1's; `boolean-piece-forms-expected.tsv` gives per case
`expect KIND S9e.4b.3c.3a CLASS` (`solid`, `empty`, `degenerate`: a
`Degenerate` of S9's rules; the class `pieces`: an imported body and a
construction, `both`: two imported inputs, `chain`: an imported input's
result given to another Boolean), a degenerate case's `reason TEXT`, then
`result N volume area cx cy cz` (totals over the N solids, world
coordinates) or `empty`; `boolean-piece-forms-frames.tsv` the stored axes of
every solid the kernel builds from its rows.

The bodies, each `X op Y` of a box (a prism of a square) and a primitive
whose sections by the box's planes are circles and lines (the `.brep`
reader takes no ellipse, and a sphere's section off its stored frame's
parallels and meridians has a pcurve the converter does not certify,
S9e.4b.3a's evidence), written to `imported/form_NAME.brep`: reproducing
the DRAW survey's `bcut_complex/G4`, `boss`, a box on the world's axes
fused with a cylinder of radius 5/2 along the world's `y` from its face `y
= 0` to a cap at `y = -4` (the part, a box with a cylindrical boss whose
cap and whose end on the box's face lie on planes of the world's axes);
`scoop`, a box on `SKEW4` less a ball of radius 9/2 centred on its top face
(the ball's frame normal along the box's `y`, its `x` the box's normal: the
top face a meridian plane, in two faces on one plane, the side faces its
parallels' planes, its seam above the body); `slot`, a box on the world's
axes less a cylinder of radius 5/4 along `x` whose axis lies 1/2 above its
top face, ending inside the box at its cap (a groove through a box is
S9e.4a's prism); `dimple`, a box on `SKEW` less a ball of radius 5/2 about
a point 3/2 above its top face's middle (the ball's frame `SKEW`'s, its
section a parallel); `ball_boss`, a box on `SKEW2` fused with a ball of
radius 5/2 about a point 1 below its top face's middle; `sink`, a box on
`SKEW4` less a frustum along its normal through it (radii 1 and 3 over 6
from 1 below it: a conical hole); `bite`, a cylinder of radius 3 on the
world's axes less a box across its wall between its caps (two planes along
its axis, two normal to it, its caps whole: in a turned frame a plane along
a cylinder's axis is within rounding of it, S9's refusal); declared
`degenerate`: `notch`, reproducing `bcut_complex/I6`'s tool, a box on the
world's axes less a cylinder of radius 3 along `z` about the middle of one
face, its wall tangent to the two faces it meets (its frame's `x` reversed:
its seam off the body); `quarter`, reproducing `shading_132`, a frustum of
radii 2 and 1 over 2 on the world's axes less the quarter between two
planes through its axis (its frame's `x` between them: its seam off the
body). S9e.4b.3a's `bitten` (a ball less a box's corner), declared
S9e.4b.3c's there, is this step's too (its set's generator declares it
solid).

Cases (each the three operations): `boss_rod`, the boss and a cylinder of
radius 2 along `y` whose caps lie on the boss's cap plane and on the box's
face (`G4`'s tool: two cylinders crossing along two lines, faces of both on
two planes); `boss_slab`, the boss and a `TILT` slab across its cylinder
and the box; `scoop_rod`, the scoop and an upright rod through its floor;
`slab_scoop`, a `TILT` slab less the scoop (the body the tool);
`slot_rod`, the slot and an upright rod through its floor; `dimple_ball`,
the dimple and a ball across its sphere face; `ball_boss_slab`, the ball's
boss and a `TILT` slab across its sphere; `sink_rod`, the conical hole and
an upright rod across its wall; `bite_box`, the bitten cylinder and a
`TILT` box across the bite; `pieces`, the boss and the slot, both imported;
`chain_slot`, the slot less a rod, then with a `TILT` slab. Declared
`degenerate`: `notch_box`, the notch (its wall tangent to its own faces)
and a box; `quarter_ball`, the three-quarter frustum (its planes through
its virtual apex) and a ball.

The reference is the construction OCCT was given, through S9e.3a's chained
reference (`chained_curved_boolean_reference.py`): each imported body its
first Boolean, `(X op Y) op2 C` (swapped, the body the tool), two bodies
`(X op Y) op2 (Z op W)`, a chain `((X op Y) op2 C) op3 D`, with S9e.4b.3a's
checks (the two families, each solid's closed form, the pair identities on
the last Boolean's arguments, Monte Carlo, solid counts by rays at two
resolutions, each body one solid, the meetings' sines and the events'
spacing outside the declared groups) and each body's closed form where it
has one (`closed`). With `--check`, where the bodies' files exist, each is
read (`stored_records`): its faces' kinds one sphere, cylinder or cone and
planes, every stored vertex within 1e-12 of the size on the construction's
surfaces, and none of S9e.4b.3a's pieces (its primitive's material outside
it, or a plane of the box cutting no edge of the body's primitive's
common): the reason each body is this step's.
"""
import argparse
from fractions import Fraction as F
import os
import struct
import zlib

import mpmath as mp

from identity_reference import Case, native_case
from curve_surface_reference import stored_axes
import generate_imported_pieces_boolean_fixtures as pieces
from generate_curved_boolean_fixtures import square, disc

base = pieces.base
ROOT = base.ROOT
OPS = base.OPS
OPERATIONS, BOOLEANS = base.OPERATIONS, base.BOOLEANS
CLASSES = ('pieces', 'both', 'chain')
PREFIX = 'boolean-piece-forms'
BODIES_FILE = f'{PREFIX}-bodies.txt'
STEP = 'S9e.4b.3c.3a'
# `ALONGY`: the world's `y` as the normal, its `z` as the `x` (G4's files'
# cylinders): a frame's local `(u, v, w)` the world's `(z, x, y)`.
FRAMES = dict(pieces.FRAMES, ALONGY=(0.0, 1.0, 0.0, 0.0, 0.0, 1.0))
prism, sphere, cone, construction = base.prism, base.sphere, base.cone, base.construction
ref = base.ref


def at(name, origin):
    return tuple(float(c) for c in origin)+FRAMES[name]


def q(c):
    return mp.mpf(F(c).numerator)/F(c).denominator


def point(frame, u, v, w):
    """The frame's point `o + u x + v y + w n` (its stored axes), rounded
    once: a solid's origin placed in another's frame."""
    o, x, y, n = stored_axes(frame)
    return tuple(float(F(o[i])+F(u)*F(x[i])+F(v)*F(y[i])+F(w)*F(n[i])) for i in range(3))


def moved(frame, origin):
    """The frame's axes about another origin."""
    return tuple(float(c) for c in origin)+frame[3:]


def reversed_x(frame):
    """The frame turned half a turn about its normal (its `x` reversed)."""
    return frame[:6]+tuple(-c for c in frame[6:9])


# ------------------------------------------------------------------ bodies

class Body:
    """An imported body `first op second`, OCCT's rows for it, its closed
    form `(volume, area)` where it has one, and the form the kernel decides
    it as (`groove`, `boss`, `bite`; `None` where declared degenerate, its
    reason `reason`)."""

    def __init__(self, name, first, op, second, form, closed=None, reason=None):
        self.name, self.first, self.op, self.second = name, first, op, second
        self.form, self.closed, self.reason = form, closed, reason
        self.path = f'imported/form_{name}.brep'

    def specs(self):
        return [self.first, self.second]

    def primitive(self):
        """The body's curved primitive (a sphere, a cone or a cylinder's
        prism of a disc)."""
        for s in self.specs():
            if s[0] in ('sphere', 'cone') or (s[0] == 'prism' and s[1][0].circle is not None):
                return s
        raise AssertionError(self.name)

    def box(self):
        return next(s for s in self.specs() if s is not self.primitive())

    def native_rows(self):
        rows = lambda spec: native_case(construction(spec, 0)).split('\n')[1:-1]
        return rows(self.first)+[f'boolean {self.op}']+rows(self.second)


def cylinder(frame, cx, cy, r, h0, h1):
    return prism([disc(cx, cy, r)], frame, h0, h1)


def block(frame, x0, y0, x1, y1, h0, h1):
    return prism([square(x0, y0, x1, y1)], frame, h0, h1)


def segment_area(r, d):
    """The area of a disc of radius `r` beyond a chord at distance `d` from
    its centre."""
    R, D = mp.mpf(r), mp.mpf(d)
    return R**2*mp.acos(D/R)-D*mp.sqrt(R**2-D**2)


def cap_volume(r, d):
    """The volume of a ball of radius `r` beyond a plane at distance `d`
    from its centre: `pi t^2 (3 r - t) / 3`, `t = r - d`."""
    R, t = mp.mpf(r), mp.mpf(r)-d
    return mp.pi*t**2*(3*R-t)/3


def box_volume(frame, v):
    """A prism's volume `v` in its chart, in the world: times the stored
    axes' determinant (a unit within rounding)."""
    _, x, y, n = (tuple(F(c) for c in a) for a in stored_axes(frame))
    det = (x[0]*(y[1]*n[2]-y[2]*n[1])-x[1]*(y[0]*n[2]-y[2]*n[0])+x[2]*(y[0]*n[1]-y[1]*n[0]))
    return q(F(v)*det)


def plane_distance(frame, w, centre):
    """The signed distance of a point from the frame's plane at local height
    `w` (the prism's cap: its chart's `w`, through the stored axes' exact
    inverse, `cross(x, y) . (X - o) / det = w`)."""
    o, x, y, n = (tuple(F(v) for v in a) for a in stored_axes(frame))
    cross = lambda a, b: (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])
    row = cross(x, y)
    det = sum(n[i]*row[i] for i in range(3))
    level = sum(row[i]*(F(centre[i])-o[i]) for i in range(3))/det-F(w)
    length = mp.sqrt(q(sum(c*c for c in row)))/q(det)
    return q(level)/length


# The bodies' frames and places.
BOSS_FRAME = at('ALONGY', (0, 0, 0))
# The scoop's slab on `SKEW4` (its normal near the world's `z`, along the
# reference's rays: its floor no sliver to them), its ball's frame normal along the slab's `y`
# and its `x` the slab's normal (the seam above the top face, out of the
# body; the top face a meridian plane of the ball's frame through its
# centre).
SCOOP_FRAME = at('SKEW4', (1, 2, 0.5))
SCOOP_CENTRE = point(SCOOP_FRAME, 6, 3, 6)
SCOOP_BALL = SCOOP_CENTRE+(4.0, -8.0, -1.0, 4.0, 1.0, 8.0)
SLOT_FRAME = at('XY', (4, 1.25, 2.5))
SLOT_AXIS = at('SIDE', (3, 1.25, 2.5))
DIMPLE_FRAME = at('SKEW', (1, 2, 0.5))
DIMPLE_CENTRE = point(DIMPLE_FRAME, 4, 4, 5.5)
BALL_BOSS_FRAME = at('SKEW2', (0.5, 1, 1.5))
BALL_BOSS_CENTRE = point(BALL_BOSS_FRAME, 4, 4, 3)
SINK_FRAME = at('SKEW4', (2, 1, 0.5))
SINK_BASE = point(SINK_FRAME, 4, 4, -1)
BITE_FRAME = at('XY', (0, 0, 0))
NOTCH_FRAME = at('XY', (0, 0, 0))
# The notch's cylinder with its `x` reversed: its seam off the body.
NOTCH_AXIS = at('XY', (0, 0, 0))[:6]+(-1.0, 0.0, 0.0)
# shading_132's frustum on the world's axes (its planes through its axis
# exactly), its frame's `x` between the removed quarter's planes: its seam
# inside that quarter, off the body.
QUARTER_FRAME = at('XY', (1, 1, 1))
QUARTER_CONE = QUARTER_FRAME[:6]+(1.0, 1.0, 0.0)


def make_bodies():
    out = [
        # G4's part: the boss on the box's face `y = 0`, from `y = -4`.
        Body('boss', block(BOSS_FRAME, 0.0, 0.0, 10.0, 10.0, 0.0, 8.0), 'fuse',
             cylinder(BOSS_FRAME, 6.0, 4.0, 2.5, -4.0, 0.0), 'boss',
             closed=lambda: (800+mp.pi*25, 520+20*mp.pi)),
        # The top face in two faces on one plane, the ball's half below it.
        Body('scoop', block(SCOOP_FRAME, 0.0, 0.0, 12.0, 6.0, 0.0, 6.0), 'cut',
             sphere(4.5, SCOOP_BALL), 'groove'),
        # A groove ending inside the box at the cylinder's cap.
        Body('slot', block(SLOT_FRAME, 0.0, 0.0, 8.0, 6.0, 0.0, 4.0), 'cut',
             cylinder(SLOT_AXIS, 3.0, 4.5, 1.25, 0.0, 6.0), 'groove',
             closed=lambda: (192-5*segment_area(1.25, 0.5), None)),
        Body('dimple', block(DIMPLE_FRAME, 0.0, 0.0, 8.0, 8.0, 0.0, 4.0), 'cut',
             sphere(2.5, moved(DIMPLE_FRAME, DIMPLE_CENTRE)), 'groove',
             closed=lambda: (box_volume(DIMPLE_FRAME, 256)-cap_volume(2.5, plane_distance(DIMPLE_FRAME, 4, DIMPLE_CENTRE)),
                             None)),
        Body('ball_boss', block(BALL_BOSS_FRAME, 0.0, 0.0, 8.0, 8.0, 0.0, 4.0), 'fuse',
             sphere(2.5, moved(BALL_BOSS_FRAME, BALL_BOSS_CENTRE)), 'boss',
             closed=lambda: (box_volume(BALL_BOSS_FRAME, 256)
                             + cap_volume(2.5, -plane_distance(BALL_BOSS_FRAME, 4, BALL_BOSS_CENTRE)),
                             None)),
        Body('sink', block(SINK_FRAME, 0.0, 0.0, 8.0, 8.0, 0.0, 4.0), 'cut',
             cone(1.0, 3.0, 6.0, moved(SINK_FRAME, SINK_BASE)), 'groove'),
        # Its caps whole, the bite between them.
        Body('bite', cylinder(BITE_FRAME, 4.0, 3.0, 3.0, 0.0, 6.0), 'cut',
             block(BITE_FRAME, 5.0, 4.0, 12.0, 12.0, 2.0, 4.5), 'bite'),
        # I6's tool: the wall tangent to `y = 0` and `y = 6` at `x = 0`.
        Body('notch', block(NOTCH_FRAME, 0.0, 0.0, 8.0, 6.0, 0.0, 3.0), 'cut',
             cylinder(NOTCH_AXIS, 0.0, -3.0, 3.0, -1.0, 4.0), None,
             closed=lambda: (144-mp.pi*27/2, None),
             reason='an imported plane piece whose curved face is tangent to its plane faces'),
        # shading_132: the quarter between two planes through the axis.
        Body('quarter', cone(2.0, 1.0, 2.0, QUARTER_CONE), 'cut',
             block(QUARTER_FRAME, 0.0, 0.0, 5.0, 5.0, -1.0, 3.0), None,
             closed=lambda: (mp.pi*14/3*3/4, None),
             reason='a plane through a cone\'s apex'),
    ]
    return {b.name: b for b in out}


BODIES = make_bodies()


def imported(name):
    return ('imported', name)


def flatten(items):
    """The reference's solids of a case's inputs (an imported body its two
    solids) and each input's expression over them."""
    specs, exprs = [], []
    for it in items:
        if it[0] == 'imported':
            b = BODIES[it[1]]
            k = len(specs)
            specs += b.specs()
            exprs.append((b.op, k, k+1))
        else:
            exprs.append(len(specs))
            specs.append(it)
    return specs, exprs


class Imported(pieces.Imported):
    """One case, its imported inputs this step's bodies."""

    def __init__(self, name, klass, items, op1, stages, kind, solids, reason=None):
        assert klass in CLASSES, klass
        self.name, self.klass, self.op1, self.stages = name, klass, op1, stages
        self.kind, self.solids, self.reason = kind, solids, reason
        self.group = name.rsplit('_', 1)[0]
        self.items = items
        self.solid_cases = []
        for k, s in enumerate(items):
            if s[0] == 'imported':
                c = Case(name, 1e-7, OPERATIONS[k], at('XY', (0, 0, 0)), 0.0, 0.0, [], brep=BODIES[s[1]].path)
            else:
                c = construction(s, OPERATIONS[k])
                c.name = name
            self.solid_cases.append(c)

    def expr(self):
        _, e = flatten(self.items)
        if not self.stages:
            return (self.op1, e[0], e[1])
        (op2, sw2) = self.stages[0]
        first = (self.op1, e[0], e[1])
        return (op2, e[2], first) if sw2 else (op2, first, e[2])


def group(name, klass, items, outcomes, first=None, reason=None, kind=None):
    """Cases of one group: the last operation varied over `outcomes` (`{op:
    solids}`, 0 for empty)."""
    out = []
    for op, n in outcomes.items():
        k = kind or ('degenerate' if reason else ('solid' if n else 'empty'))
        stages = [(op, False)] if first else []
        out.append(Imported(f'{name}_{op}', klass, items, first or op, stages, k, n, reason))
    return out


# G4's tool: along `y` from the boss's cap plane to the box's face.
BOSS_ROD = cylinder(BOSS_FRAME, 4.5, 5.5, 2.0, -4.0, 0.0)
BOSS_SLAB = block(at('TILT', (0, 0, 0)), -1.0, -20.0, 12.0, 20.0, 3.5, 6.0)
SCOOP_ROD = cylinder(at('XY', point(SCOOP_FRAME, 5.5, 2.5, 0)), 0.0, 0.0, 0.75, -4.0, 10.0)
SCOOP_SLAB = block(at('TILT', point(SCOOP_FRAME, 6, 3, 3)), -12.0, -12.0, 12.0, 12.0, -0.75, 0.75)
SLOT_ROD = cylinder(SLOT_FRAME, 3.25, 2.75, 0.75, -1.0, 6.0)
SLOT_SLAB = block(at('TILT', (0, 0, 0)), -2.0, -10.0, 10.0, 12.0, 4.75, 6.0)
DIMPLE_BALL = sphere(1.75, at('XY', point(DIMPLE_FRAME, 5.25, 4.5, 4.25)))
BALL_BOSS_SLAB = block(at('TILT', point(BALL_BOSS_FRAME, 4, 4, 4)), -12.0, -12.0, 12.0, 12.0, -0.5, 1.0)
SINK_ROD = cylinder(at('XY', (0, 0, -6)), *point(SINK_FRAME, 4.0, 5.75, 2)[:2], 0.625, 0.0, 20.0)
BITE_BOX = block(at('TILT', (0, 0, 0)), 3.5, -2.0, 9.0, 12.0, 3.0, 5.5)
CHAIN_ROD = cylinder(SLOT_FRAME, 6.25, 3.25, 0.5, -1.0, 6.0)
NOTCH_BOX = block(NOTCH_FRAME, 2.0, -1.0, 5.0, 7.0, 1.0, 2.0)
QUARTER_BALL = sphere(1.5, at('XY', point(QUARTER_FRAME, -1.0, -0.5, 1.5)))


def cases():
    out = []
    three = {'fuse': 1, 'cut': 1, 'common': 1}
    out += group('boss_rod', 'pieces', [imported('boss'), BOSS_ROD], three)
    out += group('boss_slab', 'pieces', [imported('boss'), BOSS_SLAB], {'fuse': 1, 'cut': 2, 'common': 1})
    out += group('scoop_rod', 'pieces', [imported('scoop'), SCOOP_ROD], three)
    out += group('slab_scoop', 'pieces', [SCOOP_SLAB, imported('scoop')], three)
    out += group('slot_rod', 'pieces', [imported('slot'), SLOT_ROD], three)
    out += group('dimple_ball', 'pieces', [imported('dimple'), DIMPLE_BALL], three)
    out += group('ball_boss_slab', 'pieces', [imported('ball_boss'), BALL_BOSS_SLAB],
                 {'fuse': 1, 'cut': 2, 'common': 1})
    out += group('sink_rod', 'pieces', [imported('sink'), SINK_ROD], three)
    out += group('bite_box', 'pieces', [imported('bite'), BITE_BOX], three)
    out += group('pieces', 'both', [imported('boss'), imported('slot')], three)
    out += group('chain_slot', 'chain', [imported('slot'), CHAIN_ROD, SLOT_SLAB], three, first='cut')
    out += group('notch_box', 'pieces', [imported('notch'), NOTCH_BOX], three, reason=BODIES['notch'].reason)
    out += group('quarter_ball', 'pieces', [imported('quarter'), QUARTER_BALL], three,
                 reason=BODIES['quarter'].reason)
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
            assert sum(s[0] == 'imported' for s in c.items) == 2, c.name
        if c.klass == 'chain':
            assert c.stages, c.name
        # A declared body's refusal declares its cases.
        for s in c.items:
            if s[0] == 'imported' and BODIES[s[1]].reason:
                assert c.kind == 'degenerate' and c.reason == BODIES[s[1]].reason, c.name
        specs, _ = flatten(c.items)
        for s in specs:
            assert s[0] in ('sphere', 'prism', 'cone'), c.name


# ------------------------------------------------------------------ the reference

def count_solids(chain, expr, n):
    """S9e.4b.3c.1's count of solids by rays (`generate_one_sphere_boolean_
    fixtures.count_solids`: each ray's intervals joined to a neighbouring
    ray's within two grid spacings): a face nearly along the rays (the
    scoop's side faces on `SKEW4`) leaves rays grazing it whose short
    intervals the next ray's do not overlap, which the chained reference's
    strict overlap counts as further solids."""
    import generate_one_sphere_boolean_fixtures as one_sphere
    return one_sphere.count_solids(chain, expr, n)


def evaluate(job):
    """One group on the chained reference: rows, results and checks (S9e.4b.3a's,
    each body's closed form where it has one)."""
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
    # Each imported body one solid, its closed form where it has one.
    for it, e in zip(items, sub):
        if it[0] != 'imported':
            continue
        assert count_solids(chain, e, 48) == 1, f'{name}: {it[1]} is not one solid'
        closed = BODIES[it[1]].closed
        if closed is not None:
            V, A = closed()
            v, _, a = chain.measures(e)
            d = abs(v-V)/size**3
            if A is not None:
                d = max(d, abs(a-A)/size**2)
            checks['body_closed_forms'] = max(checks.get('body_closed_forms', mp.mpf(0)), d)
    # The pair identities on the last Boolean's arguments.
    if chained:
        assert count_solids(chain, first, 48) == 1, f'{name}: the given result is not one solid'
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
        # Faces of two inputs of the case on one surface (not a body's own
        # two solids') leave the area identity.
        owner = pieces.owners(items, sub)
        shared = any(owner[i] != owner[sw.others.surfs[k][0]] for i, _, sw in chain.sweeps(0) for k in sw.skip)
        if not shared:
            checks['area_identity'] = abs(af+an-ax-ac)/size**2
    for op, (expr, solids) in exprs.items():
        for n in (48, 71):
            got = count_solids(chain, expr, n)
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
                                        bool(c.stages), {}, c.kind == 'degenerate'])
        e[4][c.last] = c.solids
    return [(name, items, op1, swapped, chained, ops, degenerate, mc_n)
            for name, (items, op1, swapped, chained, ops, degenerate) in groups.items()]


def cached(job):
    """`evaluate`, its result kept in `--cache` (a development aid)."""
    cache = os.environ.get('PIECE_FORMS_CACHE')
    if not cache:
        return evaluate(job)
    import hashlib
    import inspect
    import pickle
    from pathlib import Path
    key = hashlib.sha256((repr(job)+repr(flatten(job[1])[0])+inspect.getsource(evaluate)
                          + Path(ref.__file__).read_text()).encode()).hexdigest()[:24]
    path = Path(cache)/f'{job[0]}-{key}.pickle'
    if path.exists():
        return pickle.loads(path.read_bytes())
    try:
        out = evaluate(job)
    except Exception as e:
        raise RuntimeError(f'{job[0]}: {e!r}') from e
    path.write_bytes(pickle.dumps(out))
    return out


LIMITS = dict(base.LIMITS, body_closed_forms=1e-30)
LIMITS.update({'degenerate_'+k: 1e-15 for k in list(LIMITS) if k not in ('monte_carlo_sigma', 'least_sine',
                                                                          'event_spacing')})


def reference_checks(results):
    worst, covered = {}, {}

    def note(key, value):
        worst[key] = max(worst.get(key, mp.mpf(0)), value)
        covered[key] = covered.get(key, 0)+1

    for name, rows, _, checks, margins, stats, degenerate in results:
        # A declared degenerate group's reference meets its body's own
        # tangency or planes through its apex: its checks kept apart, looser.
        for key, value in checks.items():
            note(('degenerate_'+key) if degenerate and key != 'monte_carlo_sigma' else key, value)
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
           'degenerate case, then result N volume area cx cy cz or empty)']
    for case in all_cases():
        blocks.append(case.encode())
        row = by_group[case.group][1][case.last]
        n = 0 if row[0] == 'empty' else int(row[0].split()[1])
        want = {'solid': n > 0, 'empty': n == 0, 'degenerate': True}[case.kind]
        assert want, f'{case.name}: declared {case.kind}, the reference gives {row[0]}'
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
            for key, v in (('n', n), ('x', x), ('y', y)):
                frames.append(f'{case.name}\tsolid{k} {key}\t'+' '.join(struct.pack('>d', q).hex() for q in v))
    return {f'{PREFIX}-cases.txt': '\n'.join(blocks)+'\n',
            f'{PREFIX}-expected.tsv': '\n'.join(out)+'\n',
            f'{PREFIX}-frames.tsv': '\n'.join(frames)+'\n',
            BODIES_FILE: bodies_text()}


def native_input():
    """Every case's native rows (`compare_piece_forms_boolean.py`)."""
    return '\n'.join(c.native() for c in all_cases())+'\n'


def case_size(case):
    """The case's size: its reference solids' reach from the origin."""
    specs, _ = flatten(case.items)
    chain = ref.Chain([construction(s, 0) for s in specs])
    return float(chain.size)


# ------------------------------------------------------------------ the files' check

def axis(spec):
    """A primitive's stored axis (its frame's normal)."""
    return stored_axes(pieces.spec_frame(spec))[3]


def this_steps(body):
    """Why a body is this step's and not S9e.4b.3a's piece: its primitive's
    material outside it (the box less the primitive), the primitive fused
    with the box, or the primitive less a box, each a Boolean other than a
    common."""
    p = body.primitive()
    return (body.op, body.first is p) in (('cut', False), ('fuse', True), ('fuse', False), ('cut', True))


def check_bodies():
    """Each existing body file: its faces' kinds one sphere, cylinder or
    cone and planes, each stored vertex on the construction's surfaces
    within 1e-12 of the body's size. Returns the largest deviation and the
    files read."""
    worst, read = mp.mpf(0), 0
    for b in BODIES.values():
        assert this_steps(b), (b.name, 'S9e.4b.3a\'s piece')
        path = ROOT/'fixtures'/b.path
        if not path.exists():
            continue
        read += 1
        kinds, vertices = base.stored_records(path.read_text())
        curved = [k for k in kinds if k != 'plane']
        p = b.primitive()
        want = 'cylinder' if p[0] == 'prism' else p[0]
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
    parser.add_argument('--samples', type=int, default=100000, help='Monte-Carlo points per group')
    parser.add_argument('--only', help='one group, printed (no files)')
    parser.add_argument('--cache', help='a directory keeping each group\'s result (development)')
    args = parser.parse_args()
    if args.cache:
        os.environ['PIECE_FORMS_CACHE'] = args.cache
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
