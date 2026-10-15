#!/usr/bin/env python3
"""Fixtures for REVIEW_NOTES.md's "Moved turned bodies refined": imported
turned bodies about an axis off the world's (`(1, 2, 2) / 3`), OCCT storing
every surface's frame rounded, given to Booleans with an imported box whose
planes OCCT rounded too: a coaxial cylinder's and cone's meeting on such a
frame, which the kernel's rulings made a quartic graph over the frame's
rounded data, met by the box's rounded plane passed S9e.3b's eliminant's
budget (`ComputationLimit`) before the kernel took it as a plane section.

The bodies are OCCT's own output, as S9e.4b.4d.2a's and 2b's: `boolean-
turned-tilted-bodies.txt` holds one `write NAME imported/tilted_NAME.brep`
block per body for `occt_boolean_oracle.cpp` (for a turned body a `plane`
row through the axis, the profile's `wire` of lines, a `revolve` row turning
it a whole turn about the axis through `ORIGIN` along `(1, 2, 2)`; for the
box a `plane` row, its bottom's `wire` and a `prism` row; each with a
`tolerance` row setting every tolerance to 1e-3), which
`compare_turned_tilted_boolean.py --write-bodies` runs to write
`rust/fixtures/imported/tilted_NAME.brep` (format version 1, no
triangulations), their coordinates off the dyadic grid rounded once in the
rows' text. Nothing here reads the kernel; the files are read only to check
that each is the body it claims to be (below).

`boolean-turned-tilted-cases.txt` lists each case in the Boolean protocol, an
imported input as its one `brep imported/tilted_NAME.brep` row;
`boolean-turned-tilted-expected.tsv` gives per case `expect KIND STEP CLASS`
(`solid`; the class `pieces`), then `result N volume area cx cy cz` (totals
over the N solids, world coordinates); `boolean-turned-tilted-frames.tsv`
the stored axes of every solid the kernel builds from its rows.

Each turned body's profile is a path of dyadic points `(r, w)` from the axis
out and back to it, its point at radius `r` and height `w` the world's
`ORIGIN + r x + w (1, 2, 2)` (`x` = `(2, 1, -2) / 3`, the height's unit 3:
every primitive's frame at a dyadic point on the axis), as S9e.4b.4d.2b's an
OCCT profile and a model profile (near-flat segments level at their wider
ends). The bodies: `waist` (corners only): a disc of radius 6, a cylinder to
height 1, a cone narrowing by 2 to radius 4 at 2, a cylinder to 3 and its
disc; `pin` (S9e.4b.4d.2b's kind): a bottom whose apex stands `D` = 2^-17
above its rim of radius 6 (3 D along the axis), a cylinder to 1, a step in
to radius `4 - D / 2` rising `D` (taken as the annulus at 1), a draft cone
narrowing by 1/2 to radius 3 at 3 and its disc; `skewbox`, a box on the
`SKEW` frame (`x` = `(-1, 4, -8) / 9`, normal `(8, 4, 1) / 9`) about the
world's `(11/2, 15/4, -2)`, the bodies' point at radius 6 and height 1.

Cases (each the three operations): `waist_box` and `pin_box`, each turned
body and the box across its first corner (the waist's cylinder and cone, the
pin's step); `waist_ball`, the waist and a ball across its upper cylinder.
No partner's face is parallel or normal to the axis (S9's rule refuses a
plane within rounding of a cylinder's direction once the frames are
rounded).

The reference is the model's construction through S9e.3a's chained
reference (`chained_curved_boolean_reference.py`): each body the max-min
form of its model profile's lines (S9e.4b.4d.2a's `Lattice`), each primitive
on the tilted frame at its dyadic point on the axis (the reference takes the
kernel's stored axes of a frame, so the primitives' axes lie a rounding of
their heights apart from one line: each body's closed form within 1e-15),
the box its `SKEW` spec, then `(X) op C`, with S9e.4b.3c.3a's checks (the
two families, each solid's closed form, the pair identities, Monte Carlo,
each body one solid, solid counts by rays at two resolutions, the meetings'
sines and the events' spacing) and each body's closed form by Pappus over
its model profile. With `--check`, where the bodies' files exist, each is
read (`stored_records`): its faces' kinds its OCCT profile's (the box's six
planes), every stored vertex within 1e-12 of the size on the OCCT profile's
surfaces (the box's at its corners).
"""
import argparse
from fractions import Fraction as F
import os
import struct
import zlib

import mpmath as mp

from identity_reference import Case
from curve_surface_reference import stored_axes
import generate_imported_pieces_boolean_fixtures as pieces
import generate_piece_forms_boolean_fixtures as forms
import generate_turned_profiles_boolean_fixtures as profiles

base = pieces.base
ROOT = base.ROOT
OPS = base.OPS
OPERATIONS, BOOLEANS = base.OPERATIONS, base.BOOLEANS
CLASSES = ('pieces',)
PREFIX = 'boolean-turned-tilted'
BODIES_FILE = f'{PREFIX}-bodies.txt'
STEP = 'S9e.4b.4d.2b-moved'
prism, sphere, cone, construction = base.prism, base.sphere, base.cone, base.construction
ref = base.ref
q = forms.q
cylinder, block = forms.cylinder, forms.block
Lattice, num, at, TURN = profiles.Lattice, profiles.num, profiles.at, profiles.TURN
# A near-flat segment's rise.
D = F(1, 2**17)
# The bodies' resolution (the oracle's `tolerance` row).
RESOLUTION = '0.001'
# The axis: through `ORIGIN` along `AXIS` (length 3, the height's unit),
# the radius along `X`.
ORIGIN = (F(1, 2), F(-1, 4), F(0))
AXIS = (F(1), F(2), F(2))
X = (F(2, 3), F(1, 3), F(-2, 3))
# The axis frame: its normal along the axis, its `x` along `X`.
GEN = (1.0, 2.0, 2.0, 2.0, 1.0, -2.0)


def at_axis(origin):
    return tuple(float(c) for c in origin)+GEN


def fnum(x):
    """A coordinate's shortest round-trip text, rounded once (the profile's
    points off the dyadic grid along `X`)."""
    return repr(float(F(x)))


# ------------------------------------------------------------------ bodies

def near(a, b):
    """Whether a segment is near flat: within `4 D` of level, not level."""
    (ra, wa), (rb, wb) = a, b
    return wa != wb and abs(wb-wa) <= 4*D and ra != rb


class Body:
    """An imported turned body about the tilted axis: OCCT's profile, the
    model's (its lattice), the relation between them checked."""

    def __init__(self, name, occt, model):
        self.name, self.origin = name, ORIGIN
        self.path = f'imported/tilted_{name}.brep'
        self.occt = [(F(r), F(w)) for r, w in occt]
        self.model = [(F(r), F(w)) for r, w in model]
        self.relate()
        self.lattice = Lattice(self.model)
        self.ref_specs = self.specs()
        self.ref_expr = self.lattice.expr()

    def specs(self):
        """Each part's reference solid on the tilted axis frame: a cylinder
        (the prism of its circle) or a frustum from its lower height, the
        heights 3 times the profile's."""
        out = []
        L = self.lattice
        for p in L.parts:
            lo, hi = p['range']
            if p['k'] == 0:
                out.append(prism([forms.disc(0.0, 0.0, float(p['r0']))], self.frame_at(None), float(3*lo),
                                 float(3*hi)))
            else:
                r_lo, r_hi = L.at(p, lo), L.at(p, hi)
                for x in (r_lo, r_hi, 3*(hi-lo)):
                    assert float(x) == x, ('a frustum off the dyadic grid', x)
                out.append(cone(float(r_lo), float(r_hi), float(3*(hi-lo)), self.frame_at(lo)))
        return out

    def relate(self):
        """The model's profile is OCCT's with each near-flat segment level at
        its wider end's height, every other segment on its OCCT line; every
        point within `4 D` of OCCT's."""
        o, m = self.occt, self.model
        assert len(o) == len(m), self.name
        for i, ((a, b), (x, y)) in enumerate(zip(zip(o, o[1:]), zip(m, m[1:]))):
            if near(a, b):
                wide = a if a[0] > b[0] else b
                assert x[1] == y[1] == wide[1], (self.name, i)
            else:
                for p in (x, y):
                    assert (b[0]-a[0])*(p[1]-a[1]) == (b[1]-a[1])*(p[0]-a[0]), (self.name, i)
        for p, s in zip(o, m):
            assert abs(p[0]-s[0]) <= 4*D and abs(p[1]-s[1]) <= 4*D, (self.name, p, s)

    def flats(self):
        """OCCT's near-flat segments."""
        return [(a, b) for a, b in zip(self.occt, self.occt[1:]) if near(a, b)]

    def place(self, r, w):
        o = self.origin
        return tuple(o[i]+F(r)*X[i]+F(w)*AXIS[i] for i in range(3))

    def frame_at(self, w):
        return at_axis(self.place(0, 0 if w is None else w))

    def native_rows(self):
        o = self.origin
        wire = [f'wire S {len(self.occt)}']
        for r, w in self.occt:
            wire.append(' '.join(fnum(c) for c in self.place(r, w)))
            wire.append('L')
        a = ' '.join(num(c) for c in AXIS)
        return [f'plane {num(o[0])} {num(o[1])} {num(o[2])} 2 -2 1 2 1 -2', ' '.join(wire),
                f'revolve {num(o[0])} {num(o[1])} {num(o[2])} {a} {TURN}', f'tolerance {RESOLUTION}']

    def closed(self):
        """Pappus over the model profile, its heights 3 times its `w`."""
        vol, area = mp.mpf(0), mp.mpf(0)
        for (ra, wa), (rb, wb) in zip(self.model, self.model[1:]):
            if wb > wa:
                vol += q((wb-wa)*(ra*ra+ra*rb+rb*rb))
                area += q(ra+rb)*mp.sqrt(q((rb-ra)**2+9*(wb-wa)**2))
            else:
                area += q(abs(rb*rb-ra*ra))
        return mp.pi*vol, mp.pi*area

    def expr(self, k=0):
        """Its expression over the specs from index `k`."""
        return self.lattice.expr(k)

    def kinds(self):
        """Its stored faces' surface kinds: OCCT's profile's segments'."""
        out = []
        for (ra, wa), (rb, wb) in zip(self.occt, self.occt[1:]):
            out.append('plane' if wa == wb else ('cylinder' if ra == rb else 'cone'))
        return sorted(out)


class BoxBody:
    """An imported box OCCT wrote (`MakePrism` of its bottom face along its
    frame's normal) in the `SKEW` frame (`x` = `(-1, 4, -8) / 9`, `y` = `(-4,
    7, 4) / 9`, normal `(8, 4, 1) / 9`): its reference spec exact, OCCT's
    frames rounded."""

    AXES = ((F(-1, 9), F(4, 9), F(-8, 9)), (F(-4, 9), F(7, 9), F(4, 9)), (F(8, 9), F(4, 9), F(1, 9)))

    def __init__(self, name, origin, extent):
        self.name, self.origin = name, tuple(F(c) for c in origin)
        self.extent = tuple(F(c) for c in extent)
        self.path = f'imported/tilted_{name}.brep'
        x0, y0, x1, y1, h0, h1 = self.extent
        self.spec = block(at('SKEW', self.origin), float(x0), float(y0), float(x1), float(y1), float(h0),
                          float(h1))
        self.ref_specs = [self.spec]
        self.ref_expr = 0

    def point(self, u, v, h):
        x, y, n = self.AXES
        return tuple(self.origin[i]+u*x[i]+v*y[i]+h*n[i] for i in range(3))

    def corners(self):
        x0, y0, x1, y1, h0, h1 = self.extent
        return [self.point(u, v, h) for h in (h0, h1) for u, v in ((x0, y0), (x1, y0), (x1, y1), (x0, y1))]

    def expr(self, k=0):
        return k

    def native_rows(self):
        x0, y0, x1, y1, h0, h1 = self.extent
        n = self.AXES[2]
        return ['plane '+' '.join(fnum(c) for c in self.point(0, 0, h0))+' 8 4 1 -1 4 -8',
                'wire P 4 '+' '.join(' '.join(fnum(c) for c in v) for v in self.corners()[:4]),
                'prism '+' '.join(fnum((h1-h0)*c) for c in n), f'tolerance {RESOLUTION}']

    def closed(self):
        x0, y0, x1, y1, h0, h1 = self.extent
        a, b, c = x1-x0, y1-y0, h1-h0
        return q(a*b*c), q(2*(a*b+b*c+c*a))

    def kinds(self):
        return ['plane']*6


def make_bodies():
    waist = [(0, 0), (6, 0), (6, 1), (4, 2), (4, 3), (0, 3)]
    out = [
        Body('waist', waist, waist),
        Body('pin', [(0, D), (6, 0), (6, 1), (4-D/2, 1+D), (3, 3), (0, 3)],
             [(0, 0), (6, 0), (6, 1), (4, 1), (3, 3), (0, 3)]),
    ]
    out.append(BoxBody('skewbox', (F(11, 2), F(15, 4), F(-2)),
                       (F(-5, 4), F(-3, 2), F(7, 4), F(5, 4), F(-7, 4), F(3, 2))))
    return {b.name: b for b in out}


BODIES = make_bodies()


def imported(name):
    return ('imported', name)


def flatten(items):
    """The reference's solids of a case's inputs (an imported body its
    lattice's) and each input's expression over them."""
    specs, exprs = [], []
    for it in items:
        if it[0] == 'imported':
            b = BODIES[it[1]]
            k = len(specs)
            specs += b.ref_specs
            exprs.append(b.expr(k))
        else:
            exprs.append(len(specs))
            specs.append(it)
    return specs, exprs


class Imported(pieces.Imported):
    """One case, its imported input this step's body."""

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
        return (self.op1, e[0], e[1])


def group(name, items, outcomes, kind=None, reason=None):
    """Cases of one group: the operation varied over `outcomes` (`{op:
    solids}`, 0 for empty); the group's own refusal declares its cases."""
    out = []
    for op, n in outcomes.items():
        k = kind or ('solid' if n else 'empty')
        out.append(Imported(f'{name}_{op}', 'pieces', items, op, [], k, n, reason if kind else None))
    return out


# A ball across the waist's cylinder of radius 4 (about its point at radius
# 3.42 and height 2.53, the height's unit 3).
WAIST_BALL = sphere(1.25, at('XY', (5.25, 6.0, 2.75)))


def cases():
    """The cases, the groups slowest for the reference first (its workers
    take them in this order)."""
    three = {'fuse': 1, 'cut': 1, 'common': 1}
    out = []
    out += group('waist_box', [imported('waist'), imported('skewbox')], three)
    out += group('pin_box', [imported('pin'), imported('skewbox')], three)
    out += group('waist_ball', [imported('waist'), WAIST_BALL], three)
    return out


def all_cases():
    return cases()


def validate(listed):
    names = [c.name for c in listed]
    assert len(names) == len(set(names)), 'duplicate case names'
    for c in listed:
        assert len(c.items) == 2 and not c.stages, c.name
        assert c.last in OPS, c.name
        assert c.items[0][0] == 'imported', c.name
        specs, _ = flatten(c.items)
        for s in specs:
            assert s[0] in ('sphere', 'prism', 'cone'), c.name


# ------------------------------------------------------------------ the reference

def evaluate(job):
    """One group on the chained reference: rows, results and checks
    (S9e.4b.3c.3a's, each body's closed form)."""
    name, items, ops, declared, mc_n = job
    specs, sub = flatten(items)
    chain = ref.Chain([construction(s, 0) for s in specs])
    size = chain.size
    exprs = {op: ((op, sub[0], sub[1]), solids) for op, solids in ops.items()}
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
    # The imported body one solid, its closed form.
    it, e = items[0], sub[0]
    assert forms.count_solids(chain, e, 48) == 1, f'{name}: {it[1]} is not one solid'
    V, A = BODIES[it[1]].closed()
    v, _, a = chain.measures(e)
    checks['body_closed_forms'] = max(abs(v-V)/size**3, abs(a-A)/size**2)
    # The pair identities.
    vx, mx, ax = chain.measures(sub[0])
    vc, mc, ac = chain.measures(sub[1])
    m = {op: chain.measures(e) for op, (e, _) in exprs.items()}
    vf, mf, af = m['fuse']
    vn, mn, an = m['common']
    vt, mt, _ = m['cut']
    dev = abs(vf+vn-vx-vc)/size**3
    for i in range(3):
        dev = max(dev, abs(mf[i]+mn[i]-mx[i]-mc[i])/size**4)
    dev = max(dev, abs(vt-(vx-vn))/size**3, max(abs(mt[i]-(mx[i]-mn[i])) for i in range(3))/size**4)
    checks['pair_identities'] = dev
    checks['area_identity'] = abs(af+an-ax-ac)/size**2
    # Solid counts by rays.
    for op, (expr, _) in exprs.items():
        for n in (48, 71):
            got = forms.count_solids(chain, expr, n)
            assert got == ops[op], f'{name} {op}: {got} solids by rays ({n}), {ops[op]} declared'
    mcr = ref.monte_carlo(chain, [e for e, _ in exprs.values()], mc_n, zlib.crc32(name.encode()))
    z = 0.0
    for op, (expr, _) in exprs.items():
        vol, c, sv, sc = mcr[expr]
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
        e = groups.setdefault(c.group, [c.items, {}, c.kind == 'degenerate'])
        e[1][c.last] = c.solids
    return [(name, items, ops, declared, mc_n) for name, (items, ops, declared) in groups.items()]


def cached(job):
    """`evaluate`, its result kept in `--cache` (a development aid)."""
    cache = os.environ.get('TURNED_TILTED_CACHE')
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


# Each body's closed form (Pappus, exactly on one axis line) against its
# reference solids, which stand on the kernel's stored axes of the tilted frame
# (its rounded unit normal through each primitive's dyadic origin on the
# exact axis): apart from one line by a rounding of their heights.
LIMITS = dict(profiles.LIMITS, body_closed_forms=1e-15)


def reference_checks(results):
    worst, covered = {}, {}

    def note(key, value):
        worst[key] = max(worst.get(key, mp.mpf(0)), value)
        covered[key] = covered.get(key, 0)+1

    for name, rows, _, checks, margins, stats, declared in results:
        # (No declared group here.)
        for key, value in checks.items():
            note(('degenerate_'+key) if declared and key != 'monte_carlo_sigma' else key, value)
        note('quadrature_estimate', stats['quadrature'])
        assert stats['missed'] == 0, f'{name}: quadrature nodes off their structure'
        if not declared:
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
        want = {'solid': n > 0, 'empty': n == 0, 'degenerate': n == case.solids}[case.kind]
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
                frames.append(f'{case.name}\tsolid{k} {key}\t'+' '.join(struct.pack('>d', q_).hex() for q_ in v))
    return {f'{PREFIX}-cases.txt': '\n'.join(blocks)+'\n',
            f'{PREFIX}-expected.tsv': '\n'.join(out)+'\n',
            f'{PREFIX}-frames.tsv': '\n'.join(frames)+'\n',
            BODIES_FILE: bodies_text()}


def native_input():
    """Every case's native rows (`compare_turned_tilted_boolean.py`)."""
    return '\n'.join(c.native() for c in all_cases())+'\n'


def case_size(case):
    """The case's size: its reference solids' reach from the origin."""
    specs, _ = flatten(case.items)
    chain = ref.Chain([construction(s, 0) for s in specs])
    return float(chain.size)


# ------------------------------------------------------------------ the files' check

def profile_distance(body, X):
    """A world point's distance from the body's OCCT profile's surfaces (each
    segment's cone, cylinder or plane, unbounded), in its axis frame (world
    heights 3 times the profile's)."""
    o = [q(c) for c in body.origin]
    d = [X[i]-o[i] for i in range(3)]
    n = [q(c)/3 for c in AXIS]
    w = sum(d[i]*n[i] for i in range(3))
    r = mp.sqrt(sum((d[i]-w*n[i])**2 for i in range(3)))
    best = None
    for (ra, wa), (rb, wb) in zip(body.occt, body.occt[1:]):
        a, b = q(rb-ra), -3*q(wb-wa)
        dist = abs(a*(w-3*q(wa))+b*(r-q(ra)))/mp.sqrt(a*a+b*b)
        best = dist if best is None else min(best, dist)
    return best


def check_bodies():
    """Each existing body file: its faces' kinds its OCCT profile's, each
    stored vertex on that profile's surfaces within 1e-12 of the body's size.
    Returns the largest deviation and the files read."""
    worst, read = mp.mpf(0), 0
    for b in BODIES.values():
        path = ROOT/'fixtures'/b.path
        if not path.exists():
            continue
        read += 1
        kinds, vertices = base.stored_records(path.read_text())
        assert sorted(kinds) == b.kinds(), (b.name, kinds)
        size = max([1.0]+[abs(x) for v in vertices for x in v])
        for v in vertices:
            X = tuple(mp.mpf(c) for c in v)
            if isinstance(b, BoxBody):
                d = min(mp.sqrt(sum((X[i]-q(c[i]))**2 for i in range(3))) for c in b.corners())
            else:
                d = profile_distance(b, X)
            worst = max(worst, d/size)
    assert worst <= 1e-12, ('a stored vertex off its profile', mp.nstr(worst, 3))
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
        os.environ['TURNED_TILTED_CACHE'] = args.cache
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
          '-', ', '.join(f'{v} {k}' for k, v in sorted(kinds.items())))
    print('reference checks (largest deviation, relative to the case size; least_sine and event_spacing as '
          'reciprocals):', ', '.join(f'{k} {mp.nstr(v, 3)}' for k, v in sorted(worst.items())))
    print('cases per check:', ', '.join(f'{k} {v}' for k, v in sorted(covered.items())))
    print(f'bodies: {read} of {len(BODIES)} files read, every stored vertex within {mp.nstr(deviation, 3)} of the '
          'size of its profile\'s surfaces')


if __name__ == '__main__':
    main()
