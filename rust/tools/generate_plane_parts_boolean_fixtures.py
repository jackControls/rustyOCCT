#!/usr/bin/env python3
"""Fixtures for S9e.4b.4b.2b.1 of REVIEW_NOTES.md: imported bodies of several
sphere, cylinder and cone faces whose other plane faces are a primitive's
flat, a prism of one cap (a boss, its other cap hidden in the part it stands
on) or a pocket (its planes' hull turned over, less the primitives standing in
it), the chain of the leaves, the primitives and those parts (bodies OCCT
wrote to `.brep` files), given to Booleans.

The bodies are OCCT's own output, as S9e.4b.4b.2a's
(`generate_prism_leaves_boolean_fixtures.py`): `boolean-plane-parts-bodies.txt`
holds one `write NAME imported/part_NAME.brep` block per body for
`occt_boolean_oracle.cpp` (the first solid's rows, a `boolean` row and the
second's: the Boolean's one solid written), which
`compare_plane_parts_boolean.py --write-bodies` runs to write
`rust/fixtures/imported/part_NAME.brep` (format version 1, no
triangulations). Nothing here reads the kernel; the files are read only to
check that each is the body it claims to be (below).

`boolean-plane-parts-cases.txt` lists each case in the Boolean protocol, an
imported input as its one `brep imported/part_NAME.brep` row, a chain's
further `then` rows as S9e.1's; `boolean-plane-parts-expected.tsv` gives per
case `expect KIND S9e.4b.4b.2b.1 CLASS` (`solid`, `empty`, `degenerate`: a
`Degenerate` of S9's rules or the body's own; `unsupported`: `OutOfDomain`, a
later step's; the class `pieces`: an imported body and a construction,
`both`: two imported inputs, `chain`: an imported input's result given to
another Boolean), a degenerate or unsupported case's `reason TEXT`, then
`result N volume area cx cy cz` (totals over the N solids, world
coordinates) or `empty`; `boolean-plane-parts-frames.tsv` the stored axes of
every solid the kernel builds from its rows.

Each body is one Boolean OCCT was given, every section a circle or a line
(the `.brep` reader takes no ellipse): `flat` (the `boolean` fuzz target's
cone on a ball's half), a ball of radius 3 on `SKEW` below its equator's
plane (a zone from its south pole to its equator: its disc a face) fused with
a frustum along its frame's normal about `(3/4, -1/2)` from radius 5/4 at
height -1 (inside the ball) to radius 3/4 at height 5/2, its axis off the
ball's: the disc an annulus between the ball's rim and the frustum's circle,
a flat of the ball, the frustum's base hidden; `stack` (the target's
pentagon under a stadium), a hexagon of rational corners over `[0, 1]` on
`SKEW2` fused with a stadium (half discs of radius 5/2 about `(0, 0)` and
`(6, 0)`) over `[1, 3]` on the same frame, its footprint holding the
hexagon's: the hexagon a prism of one cap, its top hidden in the stadium;
`cake`, a box `9 x 7 x 2` on `SKEW4` fused with a stadium boss (half discs
of radius 3/2 about `(5/2, 7/2)` and `(13/2, 7/2)`) from inside it, height
3/2, to its cap at 4: the boss a prism of one cap with arcs, its bottom
hidden in the box; `slot` (the target's plate with a hole less its moved
copy), the square `[0, 8]^2` with a hole of radius 2 about `(4, 4)` over
`[0, 3]` on the world's axes less the square `[3/4, 35/4] x [5/4, 37/4]`
with a hole of radius 3/2 about `(11/2, 4)` over `[1, 2]`, on the world's
axes turned half a turn about `z`: a slot across the plate (its floor, roof
and two walls a pocket) leaving a crescent of the second hole's disc outside
the first (a tooth: the second hole's cylinder standing in the pocket,
meeting the first's in two lines; each hole's seam inside the other's disc,
so OCCT keeps the crescent's walls and the plate's hole's wall one face
each, no seam edge in them).

Cases (each the three operations): `slot_rod`, the slot and an upright rod
through the tooth and the plate's hole; `chain_cake`, the cake less a `TILT`
rod through the boss, then with a `TILT` box; `cake_stack`, the cake and the
stack both imported (the box across the stadium); `flat_box`, the flat and a
`TILT` box across its disc and frustum; `stack_rod`, the stack and a `TILT`
rod through both levels; `cake_ball`, the cake and a ball across the boss's
rim; `rod_flat`, a `TILT` rod less the flat (the body the tool). Declared
`degenerate`: `flat_touch`, the flat and a ball resting on its frustum's top
disc at its centre (`Degenerate("a plane crossing a sphere within the
resolution of tangency (S9d.1)")`). The reference's cost grows with a case's
faces: a first set with the flat and the stack both imported, the stack's
chain with a slab and the slot with a `TILT` box ran 18 to 26 minutes of a
worker a group, so the pair is the cake and the stack, the chain the cake's
with a box and the slot's partner an upright rod; a `TILT` box less the
flat put two of the reference's events within 6e-18 of each other wherever
it was, so the body is the tool of a rod.

The reference is the construction OCCT was given, through S9e.3a's chained
reference (`chained_curved_boolean_reference.py`): each imported body its
Boolean, the ball's half its ball common a cylinder below its equator's
plane (`Half`), each plate with a hole its box less its hole's cylinder
(`Holed`: the reference models convex profiles), with S9e.4b.3c.3a's checks
(the two families, each solid's closed form, the pair identities on the last
Boolean's arguments, Monte Carlo, solid counts by rays at two resolutions,
each body one solid, the meetings' sines and the events' spacing outside the
declared groups) and each body's closed form (times the stored axes'
determinant). With `--check`, where the bodies' files exist, each is read
(`stored_records`): its curved faces' kinds its constructions', every stored
vertex within 1e-12 of the size on the construction's surfaces.
"""
import argparse
from fractions import Fraction as F
import os
import struct
import zlib

import mpmath as mp

from identity_reference import Boundary, Case, native_case
from curve_surface_reference import stored_axes
import generate_imported_pieces_boolean_fixtures as pieces
import generate_piece_forms_boolean_fixtures as forms
import generate_piece_trees_boolean_fixtures as trees
import generate_primitive_chains_boolean_fixtures as chains
import generate_prism_leaves_boolean_fixtures as leaves
from generate_curved_boolean_fixtures import square, disc

base = pieces.base
ROOT = base.ROOT
OPS = base.OPS
OPERATIONS, BOOLEANS = base.OPERATIONS, base.BOOLEANS
CLASSES = ('pieces', 'both', 'chain')
PREFIX = 'boolean-plane-parts'
BODIES_FILE = f'{PREFIX}-bodies.txt'
STEP = 'S9e.4b.4b.2b.1'
FRAMES = forms.FRAMES
HP = base.HP
prism, sphere, cone, construction = base.prism, base.sphere, base.cone, base.construction
ref = base.ref
q, point, moved = forms.q, forms.point, forms.moved
cylinder, block = forms.cylinder, forms.block
det, stadium, stadium_area, height = chains.det, leaves.stadium, leaves.stadium_area, leaves.height


def at(name, origin):
    return tuple(float(c) for c in origin)+FRAMES[name]


# ------------------------------------------------------------------ bodies

class Half:
    """A ball's half below its equator's plane for OCCT (`native`: a sphere's
    row from its south pole to its equator) and the reference's same set, its
    whole ball common a cylinder on its frame below its equator's plane, of
    twice its radius, four radii deep."""

    def __init__(self, r, frame):
        self.native = sphere(r, frame, -HP, 0.0)
        self.parts = [sphere(r, frame), cylinder(frame, 0.0, 0.0, 2.0*r, -4.0*r, 0.0)]
        self.op = 'common'


class Holed:
    """A prism of a square with a round hole for OCCT (`native`: one prism of
    both boundaries) and the reference's same set, the square's box less the
    hole's cylinder (past both caps)."""

    def __init__(self, frame, h0, h1, rect, hole):
        x0, y0, x1, y1 = rect
        cx, cy, r = hole
        self.native = prism([square(x0, y0, x1, y1), disc(cx, cy, r)], frame, h0, h1)
        self.parts = [block(frame, x0, y0, x1, y1, h0, h1), cylinder(frame, cx, cy, r, h0-1.0, h1+1.0)]
        self.op = 'cut'


class Body:
    """An imported body `first op second`, OCCT's rows for it, its closed
    form `(volume, area)` (the area `None` where not checked) and its
    declared refusal (`kind` `degenerate` or `unsupported`, its reason
    `reason`); `Half` and `Holed` solids the reference's sets of their own
    (`ref_specs`, `ref_expr`)."""

    def __init__(self, name, first, op, second, closed=None, kind=None, reason=None):
        self.name, self.first, self.op, self.second = name, first, op, second
        self.closed, self.kind, self.reason = closed, kind, reason
        self.path = f'imported/part_{name}.brep'
        specs, exprs = [], []
        for s in (first, second):
            if isinstance(s, (Half, Holed)):
                k = len(specs)
                specs += s.parts
                exprs.append((s.op, k, k+1))
            else:
                exprs.append(len(specs))
                specs.append(s)
        self.ref_specs, self.ref_expr = specs, (op, exprs[0], exprs[1])

    def natives(self):
        return [s.native if isinstance(s, (Half, Holed)) else s for s in (self.first, self.second)]

    def specs(self):
        return self.natives()

    def native_rows(self):
        rows = lambda spec: native_case(construction(spec, 0)).split('\n')[1:-1]
        first, second = self.natives()
        return rows(first)+[f'boolean {self.op}']+rows(second)

    def curved_kinds(self):
        return chains.Body.curved_kinds(self)


def polygon_area(points):
    """A polygon's area by the shoelace formula, exactly."""
    p = [(F(x), F(y)) for x, y in points]
    n = len(p)
    return abs(sum(p[i][0]*p[(i+1) % n][1]-p[(i+1) % n][0]*p[i][1] for i in range(n)))/2


def frustum_volume(r0, r1, h):
    return mp.pi*h/3*(r0**2+r0*r1+r1**2)


def lens_area(r1, r2, d):
    """The area common two discs of radii `r1` and `r2` whose centres lie `d`
    apart."""
    r1, r2, d = (x if isinstance(x, mp.mpf) else q(x) for x in (r1, r2, d))
    a = r1**2*mp.acos((d**2+r1**2-r2**2)/(2*d*r1))
    b = r2**2*mp.acos((d**2+r2**2-r1**2)/(2*d*r2))
    return a+b-mp.sqrt((-d+r1+r2)*(d+r1-r2)*(d-r1+r2)*(d+r1+r2))/2


# The bodies' frames and places.
STACK_FRAME = at('SKEW2', (3, 4.5, 5))
# The ball's centre above the stack's stadium, on the stack's frame.
FLAT_FRAME = at('SKEW', forms.point(STACK_FRAME, 6.5, 1.25, 4.5))
FLAT_CONE = cone(1.25, 0.75, 3.5, moved(FLAT_FRAME, point(FLAT_FRAME, 0.75, -0.5, -1.0)))
HEXAGON = [(1.0, -1.5), (4.0, -1.5), (5.5, 0.0), (4.0, 1.5), (1.0, 1.5), (-0.5, 0.0)]
# The cake across the stack's stadium (`cake_stack`).
CAKE_FRAME = at('SKEW4', point(STACK_FRAME, 4.0, -1.0, 2.0))
CAKE_BOSS = moved(CAKE_FRAME, point(CAKE_FRAME, 2.5, 3.5, 0.0))
SLOT_FRAME = at('XY', (0, 0, 0))
# The moved copy half a turn about `z` (local `(u, v)` the world's `(10 - x,
# 10 - y)`): its hole's seam toward the world's `-x`, inside the plate's
# hole, and the plate's hole's seam inside its disc.
SLOT_TOOL = forms.reversed_x(at('XY', (10, 10, 0)))


def flat_volume():
    """The ball's half (`18 pi`) and the frustum above its disc: the disc's
    plane is the frustum's chart's height `1` within rounding (its origin's
    height `-1` in the ball's chart, exactly), its radius there `5/4 - w/7`;
    the part above it, to `7/2`, in the frustum's chart times its axes'
    determinant."""
    w = -height(FLAT_FRAME, FLAT_CONE[4][:3])
    r = lambda h: q(F(5, 4))-q(h)/7
    return 18*mp.pi+frustum_volume(r(w), q(F(3, 4)), q(F(7, 2))-q(w))*det(FLAT_FRAME)


def boss_height():
    """The cake's boss's frame's origin's height in the box's chart (0
    within rounding): the boss's top, 4 in its own chart, above the box's
    top 2 by `2 + h`."""
    return q(height(CAKE_FRAME, CAKE_BOSS[:3]))


def slot_volume():
    """The plate `3 (64 - 4 pi)` less the slot: over its height 1 the
    rectangle `[3/4, 8] x [5/4, 8]` less both holes' discs (each inside it),
    their union the two discs less their lens (their centres 3/2 apart)."""
    union = mp.pi*(4+q(F(9, 4)))-lens_area(2, F(3, 2), F(3, 2))
    return 3*(64-4*mp.pi)-(q(F(29, 4))*q(F(27, 4))-union)


def make_bodies():
    out = [
        Body('flat', Half(3.0, FLAT_FRAME), 'fuse', FLAT_CONE, closed=lambda: (flat_volume(), None)),
        Body('stack', prism([Boundary(points=HEXAGON)], STACK_FRAME, 0.0, 1.0), 'fuse',
             prism([stadium(6, 2.5)], STACK_FRAME, 1.0, 3.0),
             closed=lambda: ((q(polygon_area(HEXAGON))+2*stadium_area(6, 2.5))*det(STACK_FRAME), None)),
        Body('cake', block(CAKE_FRAME, 0.0, 0.0, 9.0, 7.0, 0.0, 2.0), 'fuse',
             prism([stadium(4, 1.5)], CAKE_BOSS, 1.5, 4.0),
             closed=lambda: ((126+(2+boss_height())*stadium_area(4, 1.5))*det(CAKE_FRAME), None)),
        Body('slot', Holed(SLOT_FRAME, 0.0, 3.0, (0.0, 0.0, 8.0, 8.0), (4.0, 4.0, 2.0)), 'cut',
             Holed(SLOT_TOOL, 1.0, 2.0, (1.25, 0.75, 9.25, 8.75), (4.5, 6.0, 1.5)),
             closed=lambda: (slot_volume(), None)),
    ]
    return {b.name: b for b in out}


BODIES = make_bodies()


def imported(name):
    return ('imported', name)


def flatten(items):
    """The reference's solids of a case's inputs (an imported body its own
    solids) and each input's expression over them."""
    specs, exprs = [], []
    for it in items:
        if it[0] == 'imported':
            b = BODIES[it[1]]
            k = len(specs)
            specs += b.ref_specs
            exprs.append(trees.shifted(b.ref_expr, k))
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


def group(name, klass, items, outcomes, first=None, kind=None, reason=None):
    """Cases of one group: the last operation varied over `outcomes` (`{op:
    solids}`, 0 for empty); a declared body's refusal, or the group's own
    (`kind`, `reason`), declares its cases."""
    out = []
    declared = [BODIES[s[1]] for s in items if s[0] == 'imported' and BODIES[s[1]].kind]
    for op, n in outcomes.items():
        if declared:
            k, why = declared[0].kind, declared[0].reason
        elif kind:
            k, why = kind, reason
        else:
            k, why = ('solid' if n else 'empty'), None
        stages = [(op, False)] if first else []
        out.append(Imported(f'{name}_{op}', klass, items, first or op, stages, k, n, why))
    return out


FLAT_BOX = block(at('TILT', point(FLAT_FRAME, 1.5, -0.5, 0.25)), -1.25, -1.5, 1.25, 1.5, -1.0, 1.0)
STACK_ROD = cylinder(at('TILT', point(STACK_FRAME, 4.25, 0.75, 1.25)), 0.0, 0.0, 1.0, -5.0, 5.0)
CAKE_BALL = sphere(1.75, at('XY', point(CAKE_FRAME, 8.0, 3.0, 2.5)))
SLOT_ROD = cylinder(SLOT_FRAME, 6.25, 5.0, 0.75, -1.0, 4.0)
ROD_FLAT = cylinder(at('TILT', point(FLAT_FRAME, 1.25, 0.75, 0.625)), 0.0, 0.0, 1.0, -5.0, 5.0)
CHAIN_ROD = cylinder(at('TILT', point(CAKE_FRAME, 5.25, 3.25, 2.0)), 0.0, 0.0, 0.75, -5.0, 5.0)
CHAIN_BOX = block(at('TILT', point(CAKE_FRAME, 6.5, 4.25, 2.75)), -1.5, -1.5, 1.5, 1.5, -1.0, 1.0)
FLAT_TOUCH = sphere(1.0, at('XY', point(FLAT_CONE[4], 0.0, 0.0, 4.5)))


def cases():
    """The cases, the groups slowest for the reference first (its workers
    take them in this order)."""
    out = []
    three = {'fuse': 1, 'cut': 1, 'common': 1}
    out += group('slot_rod', 'pieces', [imported('slot'), SLOT_ROD], three)
    out += group('chain_cake', 'chain', [imported('cake'), CHAIN_ROD, CHAIN_BOX], three, first='cut')
    out += group('cake_stack', 'both', [imported('cake'), imported('stack')], three)
    out += group('flat_box', 'pieces', [imported('flat'), FLAT_BOX], three)
    out += group('stack_rod', 'pieces', [imported('stack'), STACK_ROD], three)
    out += group('cake_ball', 'pieces', [imported('cake'), CAKE_BALL], three)
    out += group('rod_flat', 'pieces', [ROD_FLAT, imported('flat')], three)
    out += group('flat_touch', 'pieces', [imported('flat'), FLAT_TOUCH], {'fuse': 1, 'cut': 1, 'common': 0},
                 kind='degenerate',
                 reason='a plane crossing a sphere within the resolution of tangency (S9d.1)')
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
        specs, _ = flatten(c.items)
        for s in specs:
            assert s[0] in ('sphere', 'prism', 'cone'), c.name


# ------------------------------------------------------------------ the reference

def evaluate(job):
    """One group on the chained reference: rows, results and checks
    (S9e.4b.3c.3a's, each body's closed form)."""
    name, items, op1, swapped, chained, ops, declared, mc_n = job
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
    # Each imported body one solid, its closed form.
    for it, e in zip(items, sub):
        if it[0] != 'imported':
            continue
        assert forms.count_solids(chain, e, 48) == 1, f'{name}: {it[1]} is not one solid'
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
        assert forms.count_solids(chain, first, 48) == 1, f'{name}: the given result is not one solid'
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
        # solids') leave the area identity.
        owner = pieces.owners(items, sub)
        shared = any(owner[i] != owner[sw.others.surfs[k][0]] for i, _, sw in chain.sweeps(0) for k in sw.skip)
        if not shared:
            checks['area_identity'] = abs(af+an-ax-ac)/size**2
    for op, (expr, _) in exprs.items():
        for n in (48, 71):
            got = forms.count_solids(chain, expr, n)
            assert got == ops[op], f'{name} {op}: {got} solids by rays ({n}), {ops[op]} declared'
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
    return name, rows, res, checks, margins, stats, declared


def jobs(mc_n):
    """One job per group."""
    groups = {}
    for c in all_cases():
        e = groups.setdefault(c.group, [c.items, c.op1 if c.stages else None, bool(c.stages and c.stages[0][1]),
                                        bool(c.stages), {}, c.kind in ('degenerate', 'unsupported')])
        e[4][c.last] = c.solids
    return [(name, items, op1, swapped, chained, ops, declared, mc_n)
            for name, (items, op1, swapped, chained, ops, declared) in groups.items()]


def cached(job):
    """`evaluate`, its result kept in `--cache` (a development aid)."""
    cache = os.environ.get('PLANE_PARTS_CACHE')
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


LIMITS = forms.LIMITS
# Groups whose events' spacing is not checked (their sines and gaps are).
SPACING_UNCHECKED = set()


def reference_checks(results):
    worst, covered = {}, {}

    def note(key, value):
        worst[key] = max(worst.get(key, mp.mpf(0)), value)
        covered[key] = covered.get(key, 0)+1

    for name, rows, _, checks, margins, stats, declared in results:
        # A declared group's reference meets a tangency: its checks kept
        # apart, looser.
        for key, value in checks.items():
            note(('degenerate_'+key) if declared and key != 'monte_carlo_sigma' else key, value)
        note('quadrature_estimate', stats['quadrature'])
        assert stats['missed'] == 0, f'{name}: quadrature nodes off their structure'
        if not declared:
            note('least_sine', 1/max(margins['sine'], mp.mpf(10)**-30))
            if name not in SPACING_UNCHECKED:
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
            if item[0] == 'imported':
                continue
            _, x, y, n = stored_axes(c.frame)
            for key, v in (('n', n), ('x', x), ('y', y)):
                frames.append(f'{case.name}\tsolid{k} {key}\t'+' '.join(struct.pack('>d', q_).hex() for q_ in v))
    return {f'{PREFIX}-cases.txt': '\n'.join(blocks)+'\n',
            f'{PREFIX}-expected.tsv': '\n'.join(out)+'\n',
            f'{PREFIX}-frames.tsv': '\n'.join(frames)+'\n',
            BODIES_FILE: bodies_text()}


def native_input():
    """Every case's native rows (`compare_plane_parts_boolean.py`)."""
    return '\n'.join(c.native() for c in all_cases())+'\n'


def case_size(case):
    """The case's size: its reference solids' reach from the origin."""
    specs, _ = flatten(case.items)
    chain = ref.Chain([construction(s, 0) for s in specs])
    return float(chain.size)


# ------------------------------------------------------------------ the files' check

def check_bodies():
    """Each existing body file: its curved faces' kinds its constructions'
    (compared as sets, at least two curved faces), each stored vertex on the
    constructions' surfaces within 1e-12 of the body's size. Returns the
    largest deviation and the files read."""
    worst, read = mp.mpf(0), 0
    for b in BODIES.values():
        path = ROOT/'fixtures'/b.path
        if not path.exists():
            continue
        read += 1
        kinds, vertices = base.stored_records(path.read_text())
        curved = sorted(set(k for k in kinds if k != 'plane'))
        assert curved == sorted(set(b.curved_kinds())), (b.name, kinds)
        assert len([k for k in kinds if k != 'plane']) >= 2, (b.name, kinds)
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
        os.environ['PLANE_PARTS_CACHE'] = args.cache
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
