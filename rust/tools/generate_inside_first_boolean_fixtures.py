#!/usr/bin/env python3
"""Fixtures for S9e.4b.4c.2b.2 of REVIEW_NOTES.md: imported bodies of
several primitives that are a Boolean chain only in the inside-first order
(the primitives whose material lies inside their quadric first, widest
first, then those outside, widest first, the prisms of one cap after them),
given to Booleans.

The bodies are OCCT's own output, as S9e.4b.4c.2b.1's
(`generate_nested_leaves_boolean_fixtures.py`):
`boolean-inside-first-bodies.txt` holds one `write NAME
imported/inside_NAME.brep` block per body for `occt_boolean_oracle.cpp` (a
solid's rows, a `boolean` row and the other's: the Boolean's one solid
written), which `compare_inside_first_boolean.py --write-bodies` runs to
write `rust/fixtures/imported/inside_NAME.brep` (format version 1, no
triangulations). Nothing here reads the kernel; the files are read only to
check that each is the body it claims to be (below).

`boolean-inside-first-cases.txt` lists each case in the Boolean protocol,
an imported input as its one `brep imported/inside_NAME.brep` row, a chain's
further `then` rows as S9e.1's; `boolean-inside-first-expected.tsv` gives
per case `expect KIND S9e.4b.4c.2b.2 CLASS` (`solid`, `empty`, `degenerate`:
a `Degenerate` of S9's rules; `unsupported`: `OutOfDomain`, a later step's;
the class `pieces`: an imported body and a construction, `both`: two
imported inputs, `chain`: an imported input's result given to another
Boolean), a degenerate or unsupported case's `reason TEXT`, then `result N
volume area cx cy cz` (totals over the N solids, world coordinates) or
`empty`; `boolean-inside-first-frames.tsv` the stored axes of every solid
the kernel builds from its rows.

The bodies, every section a circle or a line, each about the world's `z`
through its origin, its turned part OCCT's `BRepPrimAPI_MakeRevol` of its
profile (a body OCCT gives a Boolean of a rod and a ball stores a meeting
the reader's validator refuses, `PcurveOffEdge`): `dimple` (the dataset's
`OCC14536-sh1`'s form), a rod of radius 3 over `[0, 6]` less a ball of
radius 5 about its axis's point at height `-4`, which meets the rod's wall
exactly at its bottom rim (`3^2 + 4^2 = 5^2`: no bottom face, the dimple the
rod's whole end); `cup`, a rod of radius 3 over `[0, 6]` less a ball of
radius 5 about its axis's point at height `-9/2` (an annulus at the bottom
about a dimple, of radius `sqrt 19 / 2`); in both the ball, the widest
primitive, has its material outside its sphere, which S9e.4b.4b.1's chain
(widest first) refuses; `flange` (the dataset's `CTO904_cts21347a`'s form),
a disc of radius 4 over `[0, 3/2]` with a recess between radii 2 and 3 from
its bottom to `1/2` and a spigot of radius 2 down to `-1/2`, fused with a
lug standing on its top over `[3/2, 5/2]`, the prism of the outline from
`(1, -3/4)` to `(2, -3/2)` along the ray `(4, -3)`, the arc of radius `5/2`
about the axis to `(2, 3/2)`, back along the ray `(4, 3)` to `(1, 3/4)` and
across (not turned, and no leaf: a prism of one cap, as
`cts21347a`'s lug), whose chains build the spigot, a prism of one cap,
before the recess's cylinder is cut.

Cases (each the three operations): `dimple_rod`, the dimple and an upright
rod on its frame across its wall and its dimple; `cup_box`, the cup and a
box on its frame across its wall, its bottom annulus and its dimple;
`flange_rod`, the flange and an upright rod across the spigot's wall, the
recess's wall and the lug; `flange_box`, the flange and a box across the
lug and the disc's rim; `rod_dimple`, an upright rod less the dimple (the
body the tool) across its wall; `flange_cyl`, the flange and S9e.4a's
imported cylinder `cyl` (`generate_imported_boolean_fixtures.py`: radius 3
about the world's `(5, 5)` over `[-1, 5]`) across the disc, the recess and
the lug; `chain_cup`, the cup less an upright rod across its wall, then
with a box on its frame across the rod's hole and the dimple. Declared
`degenerate`: `flange_touch`, the flange and a ball of radius 1 resting on
the lug's top at its centre, exactly (`Degenerate("a tangency between the
inputs (S9c)")`).

The reference is the construction OCCT was given, through S9e.3a's chained
reference (`chained_curved_boolean_reference.py`): the dimple and the cup
their rod less their ball (the dimple's rod from 1 below its rim, inside the
ball there, so that no plane meets the rim), the flange its disc less the recess's cylinder
(radius 3 over `[-1, 1/2]`) fused with the spigot's (radius 2 over `[-1/2,
1]`) and the lug (the turned profile as a Boolean of its cylinders),
`(X op Y) op2 C` (swapped, the body the tool), two bodies `(...) op2 Z`, a
chain `((...) op2 C) op3 D`, with S9e.4b.3c.3a's checks and each body's
closed form (`148 pi / 3`, `54 pi - 29 pi / 24`, each times its stored
axes' determinant; `47 pi / 2 + 2`). With `--check`, where the bodies'
files exist, each is read (`stored_records`): its curved faces' kinds its
primitives', every stored vertex within 1e-12 of the size on the
constructions' surfaces.
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
import generate_crossed_leaves_boolean_fixtures as crossed
import generate_turned_profiles_boolean_fixtures as profiles

base = pieces.base
ROOT = base.ROOT
OPS = base.OPS
OPERATIONS, BOOLEANS = base.OPERATIONS, base.BOOLEANS
CLASSES = ('pieces', 'both', 'chain')
PREFIX = 'boolean-inside-first'
BODIES_FILE = f'{PREFIX}-bodies.txt'
STEP = 'S9e.4b.4c.2b.2'
FRAMES = crossed.FRAMES
prism, sphere, construction = base.prism, base.sphere, base.construction
ref = base.ref
q, point = forms.q, forms.point
cylinder, block = forms.cylinder, forms.block
det = chains.det
num, TURN = profiles.num, profiles.TURN


def at(name, origin):
    return tuple(float(c) for c in origin)+FRAMES[name]


# ------------------------------------------------------------------ bodies

class Dimpled:
    """A rod of radius 3 over `[0, 6]` about the world's `z` through
    `origin` less a ball of radius 5 about its axis's point at height `c`
    (its top above the rod's bottom, its section there inside the rod):
    OCCT's rows (the profile revolved), the reference's solids and
    expression, its closed form."""

    def __init__(self, name, origin, c):
        self.name, self.origin, self.c = name, tuple(F(x) for x in origin), F(c)
        self.path = f'imported/inside_{name}.brep'
        self.frame = at('XY', self.origin)
        # The dimple's rod reaches 1 below its rim inside the ball (the ball
        # through the rod's bottom rim and its bottom disc: a triple meeting
        # the reference does not take), the same set less the ball.
        low = -1.0 if self.c == -4 else 0.0
        self.ref_specs = [cylinder(self.frame, 0.0, 0.0, 3.0, low, 6.0),
                          sphere(5.0, at('XY', point(self.frame, 0, 0, self.c)))]
        self.ref_expr = ('cut', 0, 1)

    def specs(self):
        return self.ref_specs

    def cap(self):
        """The ball's height above the rod's bottom `h = c + 5`."""
        return self.c+5

    def native_rows(self):
        o, c = self.origin, self.c
        place = lambda r, w: ' '.join(num(x) for x in (o[0]+F(r), o[1], o[2]+F(w)))
        # The profile from the ball's top on the axis along its arc
        # (clockwise in the profile's plane) to its section at the rod's
        # bottom (the rim, or the annulus's inner circle at `sqrt(25 - c^2)`
        # rounded once), out to the rim, up the wall, in to the axis.
        if c == -4:
            bottom = [place(3, 0)]
        else:
            with mp.workdps(50):
                x = float(q(o[0])+mp.sqrt(25-q(c)**2))
            bottom = [f'{x!r} {num(o[1])} {num(o[2])}', place(3, 0)]
        points = [place(0, c+5)]+bottom+[place(3, 6), place(0, 6)]
        segments = [f'A {place(0, c)} 5 0']+['L']*(len(points)-1)
        wire = [f'wire S {len(points)}']+[f'{p} {g}' for p, g in zip(points, segments)]
        return [f'plane {num(o[0])} {num(o[1])} {num(o[2])} 0 -1 0 1 0 0', ' '.join(wire),
                f'revolve {num(o[0])} {num(o[1])} {num(o[2])} 0 0 1 {TURN}']

    def closed(self):
        """The rod `54 pi` less the ball's cap above the rod's bottom, of
        height `h`: `pi h^2 (15 - h) / 3`."""
        h = q(self.cap())
        return 54*mp.pi-mp.pi*h**2*(15-h)/3, None

    def curved_kinds(self):
        return ['cylinder', 'sphere']


class Flange:
    """The turned flange (its profile revolved by OCCT about the world's `z`
    through `origin`) fused with its lug: OCCT's rows, the reference's
    solids (the disc, the recess's and the spigot's cylinders, the lug) and
    expression, its closed form."""

    # The profile `(r, z)`: the spigot's bottom, its wall to the recess's
    # top, the recess to its outer wall and down to the disc's bottom, the
    # disc's rim and top.
    PROFILE = [(0, F(-1, 2)), (2, F(-1, 2)), (2, F(1, 2)), (3, F(1, 2)), (3, 0), (4, 0), (4, F(3, 2)),
               (0, F(3, 2))]

    def __init__(self, name, origin):
        self.name, self.origin = name, tuple(F(c) for c in origin)
        self.path = f'imported/inside_{name}.brep'
        frame = at('XY', self.origin)
        self.frame = frame
        self.lug = prism([LUG], frame, 1.5, 2.5)
        self.ref_specs = [cylinder(frame, 0.0, 0.0, 4.0, 0.0, 1.5), cylinder(frame, 0.0, 0.0, 3.0, -1.0, 0.5),
                          cylinder(frame, 0.0, 0.0, 2.0, -0.5, 1.0), self.lug]
        self.ref_expr = ('fuse', ('fuse', ('cut', 0, 1), 2), 3)

    def specs(self):
        return self.ref_specs

    def native_rows(self):
        o = self.origin
        wire = [f'wire S {len(self.PROFILE)}']
        for r, w in self.PROFILE:
            wire.append(' '.join(num(c) for c in (o[0]+r, o[1], o[2]+w)))
            wire.append('L')
        lug = native_case(construction(self.lug, 0)).split('\n')[1:-1]
        return [f'plane {num(o[0])} {num(o[1])} {num(o[2])} 0 -1 0 1 0 0', ' '.join(wire),
                f'revolve {num(o[0])} {num(o[1])} {num(o[2])} 0 0 1 {TURN}', 'boolean fuse', *lug]

    def closed(self):
        """The disc `24 pi` less the recess `9 pi / 2`, the spigot's `4 pi`
        beyond them, the lug its sector of radius `5/2` over the angle `2
        atan(3/4)` less the triangle to its chord, `25 atan(3/4) / 4 - 3/4`,
        over a unit height."""
        return 47*mp.pi/2+q(F(25, 4))*mp.atan(q(F(3, 4)))-q(F(3, 4)), None

    def curved_kinds(self):
        return ['cylinder', 'cylinder', 'cylinder', 'cylinder']


# The lug's outline on the flange's frame, counter-clockwise.
LUG = Boundary(points=[(1.0, -0.75), (2.0, -1.5), (2.0, 1.5), (1.0, 0.75)],
               segments=[None, (0.0, 0.0, 2.5, True), None, None])
DIMPLE_ORIGIN = (1, 2, 0.5)
CUP_ORIGIN = (-1, 0.5, 1)
FLANGE_ORIGIN = (3, 3, 0.5)


def make_bodies():
    out = [
        Dimpled('dimple', DIMPLE_ORIGIN, -4),
        Dimpled('cup', CUP_ORIGIN, F(-9, 2)),
        Flange('flange', FLANGE_ORIGIN),
    ]
    return {b.name: b for b in out}


BODIES = make_bodies()
DIMPLE_FRAME = BODIES['dimple'].frame
CUP_FRAME = BODIES['cup'].frame
FLANGE_FRAME = BODIES['flange'].frame


def imported(name):
    return ('imported', name)


def earlier(name):
    """S9e.4a's body `name` (its file read again)."""
    return ('base', name)


def path_of(item):
    return (BODIES[item[1]] if item[0] == 'imported' else base.BODIES[item[1]]).path


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
        elif it[0] == 'base':
            exprs.append(len(specs))
            specs.append(base.BODIES[it[1]].spec)
        else:
            exprs.append(len(specs))
            specs.append(it)
    return specs, exprs


class Imported(crossed.Imported):
    """One case, its imported inputs this step's bodies or S9e.4a's
    cylinder."""

    def __init__(self, name, klass, items, op1, stages, kind, solids, reason=None):
        assert klass in CLASSES, klass
        self.name, self.klass, self.op1, self.stages = name, klass, op1, stages
        self.kind, self.solids, self.reason = kind, solids, reason
        self.group = name.rsplit('_', 1)[0]
        self.items = items
        self.solid_cases = []
        for k, s in enumerate(items):
            if s[0] in ('imported', 'base'):
                c = Case(name, 1e-7, OPERATIONS[k], at('XY', (0, 0, 0)), 0.0, 0.0, [], brep=path_of(s))
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
    solids}`, 0 for empty); the group's own refusal (`kind`, `reason`)
    declares its cases."""
    out = []
    for op, n in outcomes.items():
        if kind:
            k, why = kind, reason
        else:
            k, why = ('solid' if n else 'empty'), None
        stages = [(op, False)] if first else []
        out.append(Imported(f'{name}_{op}', klass, items, first or op, stages, k, n, why))
    return out


DIMPLE_ROD = cylinder(DIMPLE_FRAME, 2.5, 0.25, 1.0, -1.0, 3.0)
CUP_BOX = block(CUP_FRAME, 1.25, -0.75, 4.25, 1.25, -1.0, 2.0)
FLANGE_ROD = cylinder(FLANGE_FRAME, 2.25, 0.125, 0.875, -1.0, 3.0)
FLANGE_BOX = block(FLANGE_FRAME, 2.25, -1.75, 5.0, 0.25, 0.75, 2.0)
ROD_DIMPLE = cylinder(DIMPLE_FRAME, 2.5, 0.5, 1.0, 1.5, 4.0)
CHAIN_ROD = cylinder(CUP_FRAME, -2.5, 0.5, 1.0, -1.0, 7.0)
CHAIN_BOX = block(CUP_FRAME, -3.75, -1.25, 0.25, 1.75, -1.5, 1.5)
FLANGE_TOUCH = sphere(1.0, at('XY', point(FLANGE_FRAME, 1.75, 0.0, 3.5)))
TOUCH = 'a tangency between the inputs (S9c)'


def cases():
    """The cases, the groups slowest for the reference first (its workers
    take them in this order)."""
    out = []
    three = {'fuse': 1, 'cut': 1, 'common': 1}
    out += group('flange_cyl', 'both', [imported('flange'), earlier('cyl')], three)
    out += group('chain_cup', 'chain', [imported('cup'), CHAIN_ROD, CHAIN_BOX], three, first='cut')
    out += group('dimple_rod', 'pieces', [imported('dimple'), DIMPLE_ROD], three)
    out += group('cup_box', 'pieces', [imported('cup'), CUP_BOX], three)
    out += group('flange_rod', 'pieces', [imported('flange'), FLANGE_ROD], three)
    out += group('flange_box', 'pieces', [imported('flange'), FLANGE_BOX], three)
    out += group('rod_dimple', 'pieces', [ROD_DIMPLE, imported('dimple')], three)
    out += group('flange_touch', 'pieces', [imported('flange'), FLANGE_TOUCH], {'fuse': 1, 'cut': 1, 'common': 0},
                 kind='degenerate', reason=TOUCH)
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
        V, A = BODIES[it[1]].closed()
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
    cache = os.environ.get('INSIDE_FIRST_CACHE')
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
            if item[0] in ('imported', 'base'):
                continue
            _, x, y, n = stored_axes(c.frame)
            for key, v in (('n', n), ('x', x), ('y', y)):
                frames.append(f'{case.name}\tsolid{k} {key}\t'+' '.join(struct.pack('>d', q_).hex() for q_ in v))
    return {f'{PREFIX}-cases.txt': '\n'.join(blocks)+'\n',
            f'{PREFIX}-expected.tsv': '\n'.join(out)+'\n',
            f'{PREFIX}-frames.tsv': '\n'.join(frames)+'\n',
            BODIES_FILE: bodies_text()}


def native_input():
    """Every case's native rows (`compare_inside_first_boolean.py`)."""
    return '\n'.join(c.native() for c in all_cases())+'\n'


def case_size(case):
    """The case's size: its reference solids' reach from the origin."""
    specs, _ = flatten(case.items)
    chain = ref.Chain([construction(s, 0) for s in specs])
    return float(chain.size)


# ------------------------------------------------------------------ the files' check

def check_bodies():
    """Each existing body file: its curved faces' kinds its primitives'
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
    parser.add_argument('--bodies', action='store_true', help='write the bodies\' blocks alone (no reference)')
    args = parser.parse_args()
    if args.cache:
        os.environ['INSIDE_FIRST_CACHE'] = args.cache
    if args.bodies:
        (ROOT/'fixtures'/BODIES_FILE).write_text(bodies_text())
        return
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
