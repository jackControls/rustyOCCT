#!/usr/bin/env python3
"""Fixtures for S9e.4b.4b.1 of REVIEW_NOTES.md: imported bodies of several
sphere, cylinder and cone faces whose plane faces are all ends of their
primitives, a Boolean chain of those primitives (bodies OCCT wrote to
`.brep` files), given to Booleans.

The bodies are OCCT's own output, as S9e.4b.3c.3b's
(`generate_piece_trees_boolean_fixtures.py`): `boolean-primitive-chains-
bodies.txt` holds one `write NAME imported/chain_NAME.brep` block per body
for `occt_boolean_oracle.cpp` (one solid's rows, a `boolean` row and
another solid's rows: the Boolean's one solid written), which
`compare_primitive_chains_boolean.py --write-bodies` runs to write
`rust/fixtures/imported/chain_NAME.brep` (format version 1, no
triangulations). Nothing here reads the kernel; the files are read only to
check that each is the body it claims to be (below).

`boolean-primitive-chains-cases.txt` lists each case in the Boolean
protocol, an imported input as its one `brep imported/chain_NAME.brep` row,
a chain's further `then` rows as S9e.1's; `boolean-primitive-chains-
expected.tsv` gives per case `expect KIND S9e.4b.4b.1 CLASS` (`solid`,
`empty`, `degenerate`: a `Degenerate` of S9's rules or the body's own;
`unsupported`: `OutOfDomain`, a later step's; the class `pieces`: an
imported body and a construction, `both`: two imported inputs, `chain`: an
imported input's result given to another Boolean), a degenerate or
unsupported case's `reason TEXT`, then `result N volume area cx cy cz`
(totals over the N solids, world coordinates) or `empty`;
`boolean-primitive-chains-frames.tsv` the stored axes of every solid the
kernel builds from its rows.

Each body is one Boolean OCCT was given of two primitives (a cylinder, the
prism of a disc; a cone; a ball), coaxial, so every section is a circle (a
coaxial meeting, or a plane normal to the axis: the `.brep` reader takes no
ellipse, and a sphere's section off its stored frame's parallels has a
pcurve the converter does not certify, S9e.4b.3a's evidence): `shaft`
(`bfuse_complex/E5`'s part), a cylinder of radius 3 over `[0, 4]` on
`SKEW` fused with a coaxial one of radius 3/2 over `[2, 9]` (from inside it
to a cap: its step a ring); `cup`, a cylinder of radius 3 over `[0, 5]` on
`SKEW2` less a coaxial bore of radius 2 over `[1, 7]` (a blind hole, its
floor a cap of the bore); `dome` (`bcut_complex/G9`'s `cts21128c`), a
frustum of radii 7/2 and 5/2 over 4 on the world's axes common a ball of
radius 5 about `(0, 0, -2)` (meeting at radius 3, height 2; its pole a
vertex loop); `pin` (G9's `cts21128d`), a rod of radius 7/5 over `[0, 6]`
common that ball (its rim about 14/5 high, its disc on the dome's base
plane); `bead`, a ball of radius 5 on `SKEW4` less a coaxial rod of radius 3
through it (its rims at heights `+-4`, no plane face); `knob`, a ball of
radius 2 on the world's axes about `(1, 1, 1)` fused with a coaxial rod of
radius 6/5 over `[0, 5]` (from its centre to a cap, its rim at height
8/5); declared `degenerate`: `capsule`, a rod of radius 2 over `[0, 4]`
fused with a ball of its radius about its top's centre (tangent along the
rim: `Degenerate("an imported body of several primitives whose faces are
tangent along an edge")`); declared `unsupported`: `rounded`
(`bfuse_complex/K1`'s part), a prism of a square of side 8 its corners
rounded by arcs of radius 2, less a cylinder of radius 3/2 along the world's
`y` through it (plane faces other than its primitives' ends, S9e.4b.4b.2's).

Cases (each the three operations): `shaft_box`, the shaft and a `TILT` box
across its step; `cup_ball`, the cup and a ball across its rim between its
wall and its bore; `dome_pin`, the dome and the pin, both imported (G9's
configuration: one sphere, the pin's disc inside the dome's on one plane);
`dome_slab`, the dome and a `TILT` slab; `bead_slab`, the bead and a `TILT`
slab; `knob_box`, a `TILT` box less the knob (the body the tool);
`chain_dome`, the dome less a `TILT` box biting its side, then with a
`TILT` slab. Declared `unsupported`: `rounded_rod`. The capsule is given to no
case: its own tangency along its rim keeps the reference's sweeps over twenty
minutes a group, so the kernel's tests refuse it on import alone.

The reference is the construction OCCT was given, through S9e.3a's chained
reference (`chained_curved_boolean_reference.py`): each imported body its
Boolean, `(X op Y) op2 C` (swapped, the body the tool), two bodies `(X op Y)
op2 (Z op W)`, a chain `((X op Y) op2 C) op3 D`, with S9e.4b.3c.3a's checks
(the two families, each solid's closed form, the pair identities on the
last Boolean's arguments, Monte Carlo, solid counts by rays at two
resolutions, each body one solid, the meetings' sines and the events'
spacing outside the declared groups) and each body's closed form where it
has one (`closed`: the shaft and the cup, prisms on one frame, their
volumes in their chart times the stored axes' determinant; the dome, the
pin and the knob on the world's axes). With `--check`, where the bodies'
files exist, each is read (`stored_records`): its curved faces' kinds its
primitives', every stored vertex within 1e-12 of the size on the
construction's surfaces.
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

base = pieces.base
ROOT = base.ROOT
OPS = base.OPS
OPERATIONS, BOOLEANS = base.OPERATIONS, base.BOOLEANS
CLASSES = ('pieces', 'both', 'chain')
PREFIX = 'boolean-primitive-chains'
BODIES_FILE = f'{PREFIX}-bodies.txt'
STEP = 'S9e.4b.4b.1'
FRAMES = forms.FRAMES
HP = base.HP
prism, sphere, cone, construction = base.prism, base.sphere, base.cone, base.construction
ref = base.ref
at, q, point, moved = forms.at, forms.q, forms.point, forms.moved
cylinder, block, cap_volume = forms.cylinder, forms.block, forms.cap_volume


# ------------------------------------------------------------------ bodies

class Body:
    """An imported body `first op second` of two primitives, OCCT's rows for
    it, its closed form `(volume, area)` where it has one (the area `None`
    where not checked), and its declared refusal (`kind` `degenerate` or
    `unsupported`, its reason `reason`)."""

    def __init__(self, name, first, op, second, closed=None, kind=None, reason=None):
        self.name, self.first, self.op, self.second = name, first, op, second
        self.closed, self.kind, self.reason = closed, kind, reason
        self.path = f'imported/chain_{name}.brep'
        self.ref_specs, self.ref_expr = [first, second], (op, 0, 1)

    def specs(self):
        return [self.first, self.second]

    def native_rows(self):
        rows = lambda spec: native_case(construction(spec, 0)).split('\n')[1:-1]
        return rows(self.first)+[f'boolean {self.op}']+rows(self.second)

    def curved_kinds(self):
        """The kinds of its primitives' curved faces' surfaces, sorted (a
        prism's the cylinders of its arcs and circles)."""
        out = []
        for s in self.specs():
            if s[0] == 'prism':
                for b in s[1]:
                    if b.circle is not None:
                        out.append('cylinder')
                    else:
                        out += ['cylinder' for seg in (b.segments or []) if seg is not None]
            else:
                out.append(s[0])
        return sorted(out)


def det(frame):
    """The frame's stored axes' determinant (a unit within rounding)."""
    _, x, y, n = (tuple(F(c) for c in a) for a in stored_axes(frame))
    return q(x[0]*(y[1]*n[2]-y[2]*n[1])-x[1]*(y[0]*n[2]-y[2]*n[0])+x[2]*(y[0]*n[1]-y[1]*n[0]))


def rounded_square(side, r):
    """A square `[0, side]^2` its corners rounded by arcs of radius `r`,
    counter-clockwise."""
    s = side
    points = [(r, 0), (s-r, 0), (s, r), (s, s-r), (s-r, s), (r, s), (0, s-r), (0, r)]
    arcs = [None, (s-r, r, r, True), None, (s-r, s-r, r, True), None, (r, s-r, r, True), None, (r, r, r, True)]
    return Boundary(points=[(float(x), float(y)) for x, y in points],
                    segments=[None if a is None else tuple(float(c) for c in a[:3])+(a[3],) for a in arcs])


def rod_in_ball(r, R, above):
    """A rod of radius `r` from a plane `above` over a ball's centre (of
    radius `R`) common the ball: the rod up to its rim, at `d = sqrt(R^2 -
    r^2)` above the centre, and the ball's cap beyond it."""
    d = mp.sqrt(R**2-r**2)
    return (mp.pi*r**2*(d-above)+cap_volume(R, d), None)


def knob_volume(r, R, top):
    """A ball of radius `R` fused with a rod of radius `r` from its centre
    to a cap `top` above it: the ball, and the rod above its rim at `d =
    sqrt(R^2 - r^2)` less the ball's cap beyond it."""
    d = mp.sqrt(R**2-r**2)
    return (mp.pi*mp.mpf(4)/3*R**3+mp.pi*r**2*(top-d)-cap_volume(R, d), None)


# The bodies' frames and places.
SHAFT_FRAME = at('SKEW', (1, 2, 0.5))
CUP_FRAME = at('SKEW2', (0.5, 1, 1.5))
DOME_FRAME = at('XY', (0, 0, 0))
DOME_BALL = at('XY', (0, 0, -2))
BEAD_FRAME = at('SKEW4', (2, 1, 0.5))
KNOB_FRAME = at('XY', (1, 1, 1))
CAPSULE_FRAME = at('XY', (0, 0, 0))
ROUNDED_FRAME = at('XY', (0, 0, 0))
ROUNDED_HOLE = at('ALONGY', (0, 0, 0))


def make_bodies():
    pi = mp.pi
    out = [
        Body('shaft', cylinder(SHAFT_FRAME, 0.0, 0.0, 3.0, 0.0, 4.0), 'fuse',
             cylinder(SHAFT_FRAME, 0.0, 0.0, 1.5, 2.0, 9.0),
             closed=lambda: (q(36+F(9, 4)*5)*pi*det(SHAFT_FRAME), None)),
        Body('cup', cylinder(CUP_FRAME, 0.0, 0.0, 3.0, 0.0, 5.0), 'cut',
             cylinder(CUP_FRAME, 0.0, 0.0, 2.0, 1.0, 7.0),
             closed=lambda: ((45-16)*pi*det(CUP_FRAME), None)),
        # The frustum over `[0, 2]` and the ball's cap above height 2 (its
        # centre's 4).
        Body('dome', cone(3.5, 2.5, 4.0, DOME_FRAME), 'common', sphere(5.0, DOME_BALL),
             closed=lambda: (pi*q(F(2, 3)*(F(49, 4)+F(21, 2)+9))+cap_volume(5, 4), None)),
        # The rod up to its rim, `d - 2` above its base (`d` the rim's
        # height over the ball's centre, about 24/5 for the binary64 radius
        # 1.4), and the ball's cap above it.
        Body('pin', cylinder(DOME_FRAME, 0.0, 0.0, 1.4, 0.0, 6.0), 'common', sphere(5.0, DOME_BALL),
             closed=lambda: rod_in_ball(q(1.4), 5, 2)),
        Body('bead', sphere(5.0, BEAD_FRAME), 'cut', cylinder(BEAD_FRAME, 0.0, 0.0, 3.0, -6.0, 6.0)),
        # The ball, and the rod above its rim (about 8/5 above the centre for
        # the binary64 radius 1.2) less the ball's cap above it.
        Body('knob', sphere(2.0, KNOB_FRAME), 'fuse', cylinder(KNOB_FRAME, 0.0, 0.0, 1.2, 0.0, 5.0),
             closed=lambda: knob_volume(q(1.2), 2, 5)),
        Body('capsule', cylinder(CAPSULE_FRAME, 0.0, 0.0, 2.0, 0.0, 4.0), 'fuse', sphere(2.0, at('XY', (0, 0, 4))),
             kind='degenerate', reason='an imported body of several primitives whose faces are tangent along an edge'),
        Body('rounded', prism([rounded_square(8, 2)], ROUNDED_FRAME, 0.0, 8.0), 'cut',
             cylinder(ROUNDED_HOLE, 4.0, 4.0, 1.5, -1.0, 9.0),
             kind='unsupported',
             reason='an imported body of several primitives with plane faces other than their ends (S9e.4b.4b.2)'),
    ]
    return {b.name: b for b in out}


BODIES = make_bodies()


def imported(name):
    return ('imported', name)


def shifted(expr, k):
    """An expression over a body's own solids, its leaves moved by `k`."""
    if isinstance(expr, int):
        return expr+k
    return (expr[0], shifted(expr[1], k), shifted(expr[2], k))


def flatten(items):
    """The reference's solids of a case's inputs (an imported body its two
    solids) and each input's expression over them."""
    specs, exprs = [], []
    for it in items:
        if it[0] == 'imported':
            b = BODIES[it[1]]
            k = len(specs)
            specs += b.ref_specs
            exprs.append(shifted(b.ref_expr, k))
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


def group(name, klass, items, outcomes, first=None):
    """Cases of one group: the last operation varied over `outcomes` (`{op:
    solids}`, 0 for empty); a declared body's refusal declares its cases."""
    out = []
    declared = [BODIES[s[1]] for s in items if s[0] == 'imported' and BODIES[s[1]].kind]
    for op, n in outcomes.items():
        if declared:
            kind, reason = declared[0].kind, declared[0].reason
        else:
            kind, reason = ('solid' if n else 'empty'), None
        stages = [(op, False)] if first else []
        out.append(Imported(f'{name}_{op}', klass, items, first or op, stages, kind, n, reason))
    return out


SHAFT_BOX = block(at('TILT', point(SHAFT_FRAME, 1.0, 0.5, 4.25)), -2.5, -2.0, 2.0, 2.5, -1.25, 1.75)
CUP_BALL = sphere(1.75, at('XY', point(CUP_FRAME, 2.5, 0.25, 5.0)))
DOME_SLAB = block(at('TILT', (0.5, 0.25, 1.5)), -7.0, -7.0, 7.0, 7.0, -0.75, 0.625)
BEAD_SLAB = block(at('TILT', point(BEAD_FRAME, 0.5, 1.0, 1.0)), -7.0, -7.0, 7.0, 7.0, -1.0, 1.5)
KNOB_BOX = block(at('TILT', point(KNOB_FRAME, 0.5, -0.25, 1.5)), -3.0, -3.0, 3.0, 3.0, -0.5, 1.25)
CHAIN_BOX = block(at('TILT', (2.5, 0.5, 0.75)), -1.5, -1.5, 1.5, 1.5, -1.0, 1.0)
CHAIN_SLAB = block(at('TILT', (0.25, -0.5, 2.25)), -6.0, -6.0, 6.0, 6.0, -0.5, 0.625)
ROUNDED_ROD = cylinder(ROUNDED_FRAME, 6.5, 6.5, 0.75, -1.0, 9.0)


def cases():
    """The cases, the groups slowest for the reference first (its workers
    take them in this order)."""
    out = []
    three = {'fuse': 1, 'cut': 1, 'common': 1}
    out += group('cup_ball', 'pieces', [imported('cup'), CUP_BALL], three)
    out += group('shaft_box', 'pieces', [imported('shaft'), SHAFT_BOX], three)
    out += group('dome_pin', 'both', [imported('dome'), imported('pin')], three)
    out += group('chain_dome', 'chain', [imported('dome'), CHAIN_BOX, CHAIN_SLAB], {'fuse': 1, 'cut': 2, 'common': 1},
                 first='cut')
    out += group('dome_slab', 'pieces', [imported('dome'), DOME_SLAB], {'fuse': 1, 'cut': 2, 'common': 1})
    out += group('bead_slab', 'pieces', [imported('bead'), BEAD_SLAB], {'fuse': 1, 'cut': 2, 'common': 1})
    out += group('knob_box', 'pieces', [KNOB_BOX, imported('knob')], three)
    out += group('rounded_rod', 'pieces', [imported('rounded'), ROUNDED_ROD], three)
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
    (S9e.4b.3c.3a's, each body's closed form where it has one)."""
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
    # Each imported body one solid, its closed form where it has one.
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
    for op, (expr, solids) in exprs.items():
        for n in (48, 71):
            got = forms.count_solids(chain, expr, n)
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
    cache = os.environ.get('PRIMITIVE_CHAINS_CACHE')
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
# the bead's bore, a cylinder on `SKEW4`'s stored axes (an affine model, its
# axes orthonormal within rounding), meets its ball in its rims, where the
# reference finds each rim as both surfaces' event, a rounding apart.
SPACING_UNCHECKED = {'bead_slab'}


def reference_checks(results):
    worst, covered = {}, {}

    def note(key, value):
        worst[key] = max(worst.get(key, mp.mpf(0)), value)
        covered[key] = covered.get(key, 0)+1

    for name, rows, _, checks, margins, stats, declared in results:
        # A declared group's reference meets its body's own tangency or the
        # later step's planes: its checks kept apart, looser.
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
    """Every case's native rows (`compare_primitive_chains_boolean.py`)."""
    return '\n'.join(c.native() for c in all_cases())+'\n'


def case_size(case):
    """The case's size: its reference solids' reach from the origin."""
    specs, _ = flatten(case.items)
    chain = ref.Chain([construction(s, 0) for s in specs])
    return float(chain.size)


# ------------------------------------------------------------------ the files' check

def check_bodies():
    """Each existing body file: its curved faces' kinds its primitives'
    (each primitive's one face but a sphere's or a cylinder's that OCCT
    keeps whole: compared as sets), each stored vertex on the construction's
    surfaces within 1e-12 of the body's size. Returns the largest deviation
    and the files read."""
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
        os.environ['PRIMITIVE_CHAINS_CACHE'] = args.cache
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
