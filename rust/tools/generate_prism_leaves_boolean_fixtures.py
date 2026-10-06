#!/usr/bin/env python3
"""Fixtures for S9e.4b.4b.2a of REVIEW_NOTES.md: imported bodies of several
sphere, cylinder and cone faces led by a prism leaf (a prism of lines and
arcs with both its caps, its walls' joints tangent where its corners are
rounded), the chain of the leaf and the other primitives (bodies OCCT wrote
to `.brep` files), given to Booleans.

The bodies are OCCT's own output, as S9e.4b.4b.1's
(`generate_primitive_chains_boolean_fixtures.py`): `boolean-prism-leaves-
bodies.txt` holds one `write NAME imported/leaf_NAME.brep` block per body
for `occt_boolean_oracle.cpp` (the prism's rows, a `boolean` row and the
primitive's rows: the Boolean's one solid written), which
`compare_prism_leaves_boolean.py --write-bodies` runs to write
`rust/fixtures/imported/leaf_NAME.brep` (format version 1, no
triangulations). Nothing here reads the kernel; the files are read only to
check that each is the body it claims to be (below).

`boolean-prism-leaves-cases.txt` lists each case in the Boolean protocol, an
imported input as its one `brep imported/leaf_NAME.brep` row, a chain's
further `then` rows as S9e.1's; `boolean-prism-leaves-expected.tsv` gives
per case `expect KIND S9e.4b.4b.2a CLASS` (`solid`, `empty`, `degenerate`: a
`Degenerate` of S9's rules or the body's own; `unsupported`: `OutOfDomain`,
a later step's; the class `pieces`: an imported body and a construction,
`both`: two imported inputs, `chain`: an imported input's result given to
another Boolean), a degenerate or unsupported case's `reason TEXT`, then
`result N volume area cx cy cz` (totals over the N solids, world
coordinates) or `empty`; `boolean-prism-leaves-frames.tsv` the stored axes of
every solid the kernel builds from its rows.

Each body is one Boolean OCCT was given of a prism of lines and arcs (a
square of side 9 its corners rounded by arcs of radius 2, as
`bfuse_complex/K1`'s part; or a stadium, two half discs of radius 5/2 about
`(0, 0)` and `(6, 0)` joined by lines) and a primitive, every section a
circle or a line (the `.brep` reader takes no ellipse, and a sphere's
section off its stored frame's parallels has a pcurve the converter does not
certify, S9e.4b.3a's evidence): `bore` (K1's part), the rounded square over
`[0, 6]` on the world's axes about `(1, 2, 1/2)` less a cylinder of radius
3/2 along `y` about `(9/2, *, 3)` through both its walls `y = 0` and `y = 9`
(on a turned frame its axis would lie within rounding of the walls `x = 0`
and `x = 9`, which S9's rules refuse); `boss`, the stadium over `[0, 2]` on
`SKEW2` fused with a cylinder of radius 5/4 along its axis about `(3, 0)`
over `[1, 5]` (from inside the plate to a cap above it); `dimple`, the
rounded square over `[0, 3]` on the world's axes less a ball of radius 5/2
about `(9/2, 9/2, 9/2)` (meeting the top cap at radius 2); `pocket`, the
stadium over `[0, 3]` on `SKEW` (across the boss's plate) less a frustum
along its axis about `(3, 0)` from radius 1/2 at height 1 to radius 2 at
height 4 (a conical pocket with a flat floor); `dome`, the rounded square over `[0, 3]`
on `SKEW2` fused with a ball of radius 2 about `(9/2, 9/2, 5/2)` (meeting
the top cap alone); declared `degenerate`: `post`, the rounded square over
`[0, 3]` on the world's axes fused with a cylinder of radius 2 about the
corner `(2, 2)` over `[1, 5]` (a post on the fillet's circle: the fillet's
face and the post's above it on one surface along their arc,
`Degenerate("an imported body of several primitives whose faces are
tangent along an edge")`, given to no case); declared `unsupported`:
`notch`, the rounded square over `[0, 2]` on the world's axes less a rod of
radius 3/2 along `y` about `(9/2, *, 1)` from `y = -1` to `y = 4` (a blind
notch through both caps at the rim: no cap's loop the prism's,
S9e.4b.4b.2b's).

Cases (each the three operations): `bore_rod`, the bore and a `TILT` rod
across it; `boss_box`, the boss and a `TILT` box across the boss and the
plate's rim; `dimple_ball`, the dimple and a ball across its rim;
`pocket_rod`, the pocket and a `TILT` rod through its floor; `dome_box`, a
`TILT` box less the dome (the body the tool); `boss_pocket`, the boss and
the pocket both imported (their plates crossing); `chain_dimple`, the dimple
less a `TILT` rod through its rounded corner, then with a `TILT` slab.
Declared `degenerate`: `bore_touch`, the bore and a rod of radius 1 along `x`
about `(*, 9/2, 11/2)`, tangent to the bore from above at `(9/2, 9/2, 9/2)`
(`Degenerate("a tangency between the inputs (S9c)")`; K1's own
configuration, a rod of the bore's radius whose axis crosses it, the kernel
does not refuse alike in every Boolean, so it is left to the DRAW trial). The post and the
notch are given to no case (the reference's cost grows with a case's faces:
a first set with the bore and the boss both imported, the bore's chain with a
box and a slab and the notch against a box ran over twenty minutes a group),
so the kernel's tests refuse them on import alone.

The reference is the construction OCCT was given, through S9e.3a's chained
reference (`chained_curved_boolean_reference.py`): each imported body its
Boolean, `(X op Y) op2 C` (swapped, the body the tool), two bodies `(X op Y)
op2 (Z op W)`, a chain `((X op Y) op2 C) op3 D`, with S9e.4b.3c.3a's checks
(the two families, each solid's closed form, the pair identities on the
last Boolean's arguments, Monte Carlo, solid counts by rays at two
resolutions, each body one solid, the meetings' sines and the events'
spacing outside the declared groups) and each body's closed form (its
prism's section in its chart, less or plus the primitive's part in it,
times the stored axes' determinant). With `--check`, where the bodies'
files exist, each is read (`stored_records`): its curved faces' kinds its
prism's arcs' and its primitive's, every stored vertex within 1e-12 of the
size on the construction's surfaces.
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
import generate_primitive_chains_boolean_fixtures as chains

base = pieces.base
ROOT = base.ROOT
OPS = base.OPS
OPERATIONS, BOOLEANS = base.OPERATIONS, base.BOOLEANS
CLASSES = ('pieces', 'both', 'chain')
PREFIX = 'boolean-prism-leaves'
BODIES_FILE = f'{PREFIX}-bodies.txt'
STEP = 'S9e.4b.4b.2a'
FRAMES = forms.FRAMES
HP = base.HP
prism, sphere, cone, construction = base.prism, base.sphere, base.cone, base.construction
ref = base.ref
q, point, moved = forms.q, forms.point, forms.moved
cylinder, block, cap_volume = forms.cylinder, forms.block, forms.cap_volume
det, rounded_square = chains.det, chains.rounded_square


def at(name, origin):
    return tuple(float(c) for c in origin)+FRAMES[name]


# ------------------------------------------------------------------ bodies

class Body(chains.Body):
    """An imported body `prism op primitive`, OCCT's rows for it, its closed
    form `(volume, area)` (the area `None` where not checked), and its
    declared refusal (`kind` `degenerate` or `unsupported`, its reason
    `reason`)."""

    def __init__(self, name, first, op, second, closed=None, kind=None, reason=None):
        super().__init__(name, first, op, second, closed, kind, reason)
        self.path = f'imported/leaf_{name}.brep'


def stadium(length, r):
    """Two half discs of radius `r` about `(0, 0)` and `(length, 0)` joined
    by lines, counter-clockwise."""
    points = [(0.0, -r), (length, -r), (length, r), (0.0, r)]
    segments = [None, (float(length), 0.0, float(r), True), None, (0.0, 0.0, float(r), True)]
    return Boundary(points=[(float(x), float(y)) for x, y in points], segments=segments)


def rounded_area(side, r):
    """The rounded square's area: the square less its corners' `(4 - pi)
    r^2`."""
    return q(side)**2-(4-mp.pi)*q(r)**2


def stadium_area(length, r):
    return q(length)*2*q(r)+mp.pi*q(r)**2


def frustum_volume(r0, r1, h):
    return mp.pi*q(h)/3*(q(r0)**2+q(r0)*q(r1)+q(r1)**2)


def height(frame, p):
    """A point's height `w` in a frame's chart on its stored axes, exactly
    (`p = o + u x + v y + w n`, Cramer's rule): a primitive's rounded origin
    placed in its plate's chart."""
    o, x, y, n = (tuple(F(c) for c in a) for a in stored_axes(frame))
    d = tuple(F(p[i])-o[i] for i in range(3))

    def det3(a, b, c):
        return a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0])
    return det3(x, y, d)/det3(x, y, n)


def axes(frame):
    """A frame's stored origin and axes as exact rationals."""
    return [tuple(F(c) for c in a) for a in stored_axes(frame)]


def bore_volume():
    """The bore's part between the plate's walls `v = 0` and `v = 9`: its
    frame's chart map `o' + u x' + v y' + w n'` (stored axes, exact) carries
    the disc of radius 3/2 about `(3, 9/2)` along `n'`; the walls bound `w`
    where the plate's chart's `v` (a linear function of the bore's chart) is
    0 and 9, a range of `w` alike over the disc (the walls parallel), so the
    part is the disc's area times that range times the bore's axes'
    determinant."""
    o, x, y, n = axes(BORE_FRAME)
    o2, x2, y2, n2 = axes(BORE_HOLE)

    def det3(a, b, c):
        return a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0])
    # The plate's `v` of a world vector `d`: Cramer's rule on its axes.
    dv = det3(x, n2, n)/det3(x, y, n)
    return mp.pi*q(F(9, 4))*q(9/abs(dv))*q(det3(x2, y2, n2))


def pocket_depth():
    """The pocket's frustum base's height in its plate's chart (1 within
    rounding)."""
    return height(POCKET_FRAME, POCKET_CONE[4][:3])


def dome_centre():
    """The dome's ball centre's height in its plate's chart (5/2 within
    rounding)."""
    return height(DOME_FRAME, DOME_BALL[2][:3])


def dome_gap():
    """The world distance from the dome's ball centre to the plate's top
    cap, the plane `w = 3` of its chart: `(3 - c) det / |x * y|` on the
    stored axes (the chart's `w` a unit along the normal within
    rounding)."""
    _, x, y, n = axes(DOME_FRAME)
    xy = (x[1]*y[2]-x[2]*y[1], x[2]*y[0]-x[0]*y[2], x[0]*y[1]-x[1]*y[0])
    d = xy[0]*n[0]+xy[1]*n[1]+xy[2]*n[2]
    return q((3-dome_centre())*d)/mp.sqrt(q(sum(c*c for c in xy)))


# The bodies' frames and places.
BORE_FRAME = at('XY', (1, 2, 0.5))
BORE_HOLE = moved(at('ALONGY', (0, 0, 0)), BORE_FRAME[:3])
BOSS_FRAME = at('SKEW2', (3, 4.5, 5))
DIMPLE_FRAME = at('XY', (0, 0, 0))
POCKET_FRAME = at('SKEW', point(BOSS_FRAME, 1.5, -0.5, -1.25))
DOME_FRAME = at('SKEW2', (1, -1, 0.5))
POCKET_CONE = cone(0.5, 2.0, 3.0, at('SKEW', point(POCKET_FRAME, 3.0, 0.0, 1.0)))
DOME_BALL = sphere(2.0, at('SKEW2', point(DOME_FRAME, 4.5, 4.5, 2.5)))
POST_FRAME = at('XY', (0, 0, 0))
NOTCH_ROD = at('ALONGY', (0, 0, 0))


def make_bodies():
    pi = mp.pi
    out = [
        # The rounded square over 6 less the bore across its 9.
        Body('bore', prism([rounded_square(9, 2)], BORE_FRAME, 0.0, 6.0), 'cut',
             cylinder(BORE_HOLE, 3.0, 4.5, 1.5, -1.0, 10.0),
             closed=lambda: (6*rounded_area(9, 2)*det(BORE_FRAME)-bore_volume(), None)),
        # The stadium over 2 and the boss above it, `[2, 5]`.
        Body('boss', prism([stadium(6, 2.5)], BOSS_FRAME, 0.0, 2.0), 'fuse',
             cylinder(BOSS_FRAME, 3.0, 0.0, 1.25, 1.0, 5.0),
             closed=lambda: ((2*stadium_area(6, 2.5)+3*pi*q(F(25, 16)))*det(BOSS_FRAME), None)),
        # The rounded square over 3 less the ball's cap below it, 3/2 from
        # its centre.
        Body('dimple', prism([rounded_square(9, 2)], DIMPLE_FRAME, 0.0, 3.0), 'cut',
             sphere(2.5, at('XY', (4.5, 4.5, 4.5))),
             closed=lambda: (3*rounded_area(9, 2)-cap_volume(q(F(5, 2)), q(F(3, 2))), None)),
        # The stadium over 3 less the frustum's part below its top, from
        # radius 1/2 at its base's height (1 within rounding) up to 3, its
        # slope 1/2.
        Body('pocket', prism([stadium(6, 2.5)], POCKET_FRAME, 0.0, 3.0), 'cut', POCKET_CONE,
             closed=lambda: ((3*stadium_area(6, 2.5)-frustum_volume(F(1, 2), F(1, 2)+(3-pocket_depth())/2,
                                                                     3-pocket_depth()))*det(POCKET_FRAME), None)),
        # The rounded square over 3 and the ball's cap above it, its centre
        # (5/2 high within rounding) below the cap.
        Body('dome', prism([rounded_square(9, 2)], DOME_FRAME, 0.0, 3.0), 'fuse', DOME_BALL,
             closed=lambda: (3*rounded_area(9, 2)*det(DOME_FRAME)+cap_volume(2, dome_gap()), None)),
        # A post on a rounded corner, on its fillet's circle: the fillet's
        # face and the post's above it on one surface along their arc.
        Body('post', prism([rounded_square(9, 2)], POST_FRAME, 0.0, 3.0), 'fuse',
             cylinder(POST_FRAME, 2.0, 2.0, 2.0, 1.0, 5.0),
             kind='degenerate', reason='an imported body of several primitives whose faces are tangent along an edge'),
        Body('notch', prism([rounded_square(9, 2)], POST_FRAME, 0.0, 2.0), 'cut',
             cylinder(NOTCH_ROD, 1.0, 4.5, 1.5, -1.0, 4.0),
             kind='unsupported',
             reason='an imported body of several primitives with plane faces other than their ends or a prism\'s '
                    '(S9e.4b.4b.2b)'),
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
            specs += b.ref_specs
            exprs.append(chains.shifted(b.ref_expr, k))
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


BORE_ROD = cylinder(at('TILT', point(BORE_FRAME, 6.25, 3.5, 3.25)), 0.0, 0.0, 1.0, -6.0, 6.0)
BORE_TOUCH = cylinder(moved(at('SIDE', (0, 0, 0)), BORE_FRAME[:3]), 4.5, 5.5, 1.0, -1.0, 10.0)
BOSS_BOX = block(at('TILT', point(BOSS_FRAME, 4.25, 1.5, 2.25)), -2.0, -1.75, 2.0, 2.0, -1.5, 1.5)
DIMPLE_BALL = sphere(1.75, at('XY', (6.0, 4.75, 3.25)))
POCKET_ROD = cylinder(at('TILT', point(POCKET_FRAME, 3.25, 0.5, 1.0)), 0.0, 0.0, 0.75, -3.0, 3.0)
DOME_BOX = block(at('TILT', point(DOME_FRAME, 5.5, 4.0, 3.0)), -3.0, -3.0, 3.0, 3.0, -1.25, 1.0)
CHAIN_ROD = cylinder(at('TILT', (7.25, 7.5, 1.5)), 0.0, 0.0, 1.0, -4.0, 4.0)
CHAIN_SLAB = block(at('TILT', (4.25, 4.75, 2.125)), -8.0, -8.0, 8.0, 8.0, -0.75, 0.5)


def cases():
    """The cases, the groups slowest for the reference first (its workers
    take them in this order)."""
    out = []
    three = {'fuse': 1, 'cut': 1, 'common': 1}
    out += group('boss_pocket', 'both', [imported('boss'), imported('pocket')], three)
    out += group('chain_dimple', 'chain', [imported('dimple'), CHAIN_ROD, CHAIN_SLAB],
                 {'fuse': 1, 'cut': 2, 'common': 1}, first='cut')
    out += group('bore_rod', 'pieces', [imported('bore'), BORE_ROD], three)
    out += group('boss_box', 'pieces', [imported('boss'), BOSS_BOX], three)
    out += group('dimple_ball', 'pieces', [imported('dimple'), DIMPLE_BALL], three)
    out += group('pocket_rod', 'pieces', [imported('pocket'), POCKET_ROD], three)
    out += group('dome_box', 'pieces', [DOME_BOX, imported('dome')], three)
    out += group('bore_touch', 'pieces', [imported('bore'), BORE_TOUCH], three, kind='degenerate',
                 reason='a tangency between the inputs (S9c)')
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
        for s in c.items:
            if s[0] == 'imported' and BODIES[s[1]].kind:
                assert c.kind == BODIES[s[1]].kind and c.reason == BODIES[s[1]].reason, c.name
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
        # two solids') leave the area identity.
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
    cache = os.environ.get('PRISM_LEAVES_CACHE')
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
# Groups whose events' spacing is not checked (their sines and gaps are):
# the chain's slab crosses the dimple's rim circle, where its top cap and
# its ball meet the slab's face at one point, found as several events of the
# slab's chords within 1e-21 of each other.
SPACING_UNCHECKED = {'chain_dimple'}


def reference_checks(results):
    worst, covered = {}, {}

    def note(key, value):
        worst[key] = max(worst.get(key, mp.mpf(0)), value)
        covered[key] = covered.get(key, 0)+1

    for name, rows, _, checks, margins, stats, declared in results:
        # A declared group's reference meets a tangency or the later step's
        # planes: its checks kept apart, looser.
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
    """Every case's native rows (`compare_prism_leaves_boolean.py`)."""
    return '\n'.join(c.native() for c in all_cases())+'\n'


def case_size(case):
    """The case's size: its reference solids' reach from the origin."""
    specs, _ = flatten(case.items)
    chain = ref.Chain([construction(s, 0) for s in specs])
    return float(chain.size)


# ------------------------------------------------------------------ the files' check

def check_bodies():
    """Each existing body file: its curved faces' kinds its prism's arcs'
    and its primitive's (compared as sets), each stored vertex on the
    construction's surfaces within 1e-12 of the body's size. Returns the
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
        os.environ['PRISM_LEAVES_CACHE'] = args.cache
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
