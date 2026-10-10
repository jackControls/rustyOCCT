#!/usr/bin/env python3
"""Fixtures for S9e.4b.4c.1 of REVIEW_NOTES.md: imported polyhedra (bodies of
plane faces and line edges OCCT wrote to `.brep` files, decided on their
stored vertices, S9e.4b.2's) against curved faces, their results given to
further Booleans, and imported polyhedra with a cavity.

The bodies are OCCT's own output: S9e.4b.2's (`generate_imported_polyhedra_
boolean_fixtures.py`, its files under `rust/fixtures/imported/` read again,
not written again), one of S9e.4a's (`generate_imported_boolean_fixtures.py`,
the ball) and this step's own, `boolean-polyhedra-curved-bodies.txt` holding
one `write NAME imported/poly_NAME.brep` block per body for
`occt_boolean_oracle.cpp` (a Boolean of two `box` rows: its result's one
solid written), which `compare_polyhedra_curved_boolean.py --write-bodies`
runs to write `rust/fixtures/imported/poly_NAME.brep` (format version 1, no
triangulations). Nothing here reads the kernel; the files are read only to
check that each is the body it claims to be.

`boolean-polyhedra-curved-cases.txt` lists each case in the Boolean protocol,
an imported input as its one `brep imported/NAME.brep` row, a chain's
further `then` rows as S9e.1's; `boolean-polyhedra-curved-expected.tsv` gives
per case `expect KIND S9e.4b.4c.1 CLASS` (`solid`, `empty`, `degenerate`: a
`Degenerate` of S9's rules; `unsupported`: `OutOfDomain`; the class
`curved`: an imported polyhedron and a construction with curved faces,
`both`: an imported polyhedron and an imported body with curved faces,
`chain`: an imported polyhedron's result given to another Boolean,
`cavity`: an imported polyhedron with a cavity), a degenerate or unsupported
case's `reason TEXT`, then `result N volume area cx cy cz` (totals over the
N solids, world coordinates) or `empty`; `boolean-polyhedra-curved-
frames.tsv` the stored axes of every solid the kernel builds from its rows.

The bodies: S9e.4b.2's `pyramid` (`MakeWedge` in the `TILT` frame, an apex
of four faces), `truncated` (a pyramid's frustum in the `TURN30` frame),
`wedge` (slanted on three sides, in the `R125` frame), `tetra` (a
tetrahedron of points rounded in the `TURN30` frame), `octa` (an octahedron
of exact points), `notched` (a box less a box in a skew frame, not convex)
and `hollow` (the box `[0, 10]^3` less the box `[3, 7]^3`: a cavity); S9e.4a's
`ball` (a ball of radius 3 about `(5, 5, 4)`); and this step's `cavity`, a
box `10 x 8 x 6` in the `TURN30` frame at `(1, 2, 0)` less a box `3 x 2.5 x
2` in the `TILT` frame at `(1.75, 6.375, 3)` inside it (a cavity whose faces' stored corners are
roundings in two turned frames, so each face is two triangles).

Cases (each the three operations): `pyramid_rod`, the pyramid and an upright
rod through its tilted base and its lower face; `wedge_ball`, the wedge and
a ball across its slanted top and its back face; `octa_cone`, the
octahedron and an upright frustum through two of its faces; `rod_truncated`,
a `TILT` rod less the truncated pyramid (the body the tool); `notched_ball`,
the notched box and a ball in its notch; `tetra_ball`, the tetrahedron and
S9e.4a's imported ball (both imported); `chain_pyramid`, the pyramid less a
box across one of its corners (S9b's polyhedral result), then with a ball;
`chain_wedge`, the wedge less an upright rod (a curved result), then with a
`TILT` box; `hollow_ball`, the hollow box and a ball through its wall into
its cavity (the fuse keeps the cavity, the cut opens it); `cavity_slab`, the
turned cavity and a `TILT` slab through it (the fuse's cavity in two, the
cut two solids). Declared `degenerate`: `octa_touch`, the octahedron and a
ball tangent to one of its faces within rounding (`Degenerate("a plane
crossing a sphere within the resolution of tangency (S9d.1)")`). Declared
`unsupported`: `hollow_inner_fuse`, the hollow box and a ball inside its
cavity (two solids, one of them with a cavity: `OutOfDomain("a cavity among
several solids (S9c)")`, by design).

The reference is the construction OCCT was given
(`polyhedra_curved_boolean_reference.py`: S9e.3a's chained reference with
convex hulls of exact points: a wedge's corners on its frame's stored axes,
a polyhedron's binary64 points; a Boolean of boxes its two prisms), with
S9e.4b.3c.3a's checks (the two families, each solid's closed form, each
body's closed form, the pair identities on the last Boolean's arguments, the
area identity where no two inputs share a surface, Monte Carlo, solid counts
by rays at two resolutions, each body one solid, the meetings' sines and the
events' spacing outside the declared groups). With `--check`, where the
bodies' files exist, each is read (`stored_records`): its faces' kinds
planes, every stored vertex within 1e-12 of the size on the construction's
planes.
"""
import argparse
from fractions import Fraction as F
import os
import struct
import zlib

import mpmath as mp

from identity_reference import Case
from curve_surface_reference import stored_axes
import generate_imported_boolean_fixtures as base
import generate_imported_pieces_boolean_fixtures as pieces
import generate_imported_polyhedra_boolean_fixtures as polyhedra
import generate_piece_forms_boolean_fixtures as forms
import generate_primitive_chains_boolean_fixtures as chains
import polyhedra_curved_boolean_reference as ref
import polyhedral_reference as pr

ROOT = base.ROOT
OPS = base.OPS
OPERATIONS, BOOLEANS = base.OPERATIONS, base.BOOLEANS
CLASSES = ('curved', 'both', 'chain', 'cavity')
PREFIX = 'boolean-polyhedra-curved'
BODIES_FILE = f'{PREFIX}-bodies.txt'
STEP = 'S9e.4b.4c.1'
FRAMES = polyhedra.FRAMES
HP = base.HP
prism, sphere, cone = base.prism, base.sphere, base.cone
q, point = forms.q, forms.point
cylinder, block = forms.cylinder, forms.block
det = chains.det


def at(name, origin):
    return tuple(float(c) for c in origin)+FRAMES[name]


# ------------------------------------------------------------------ bodies

def hull(points):
    """A convex hull's spec: its exact points."""
    return ('hull', tuple(tuple(F(c) for c in p) for p in points))


def construction(spec, op):
    """The reference's case of a spec: a hull's points, or S9e.4a's
    construction."""
    if spec[0] == 'hull':
        return ref.HullCase(spec[1])
    return base.construction(spec, op)


def cell_points(body):
    """S9e.4b.2's convex body's points (its one cell's)."""
    cells = body.body().cells
    assert len(cells) == 1, body.name
    out = []
    for c in cells[0]:
        for p in c:
            if p not in out:
                out.append(p)
    return out


class Body:
    """An imported body: its file, the reference's specs and expression over
    them, its closed-form volume (or None) and, for this step's own, OCCT's
    rows."""

    def __init__(self, name, path, specs, expr, closed=None, rows=None):
        self.name, self.path, self.specs, self.expr = name, path, specs, expr
        self.closed, self.rows = closed, rows

    def native_rows(self):
        return self.rows


def earlier(name):
    """An S9e.4b.2 body: a convex one its hull, a Boolean of boxes its two
    prisms."""
    b = polyhedra.BODIES[name]
    if any(r.startswith('boolean ') for r in b.rows):
        a, op, c = BOXES[name]
        return Body(name, b.path, [polyhedra.box_spec(*a), polyhedra.box_spec(*c)], (op, 0, 1),
                    closed=CLOSED.get(name))
    closed = b.closed

    def volume():
        return pr.M(F(closed()))
    return Body(name, b.path, [hull(cell_points(b))], 0, closed=volume if closed else None)


# S9e.4b.2's Booleans of boxes (`generate_imported_polyhedra_boolean_fixtures`).
BOXES = {
    'notched': ((at('XY', (0, 0, 0)), (10.0, 10.0, 4.0)), 'cut', (at('SKEW', (8.2, 7.9, 2.4)), (6.0, 6.0, 6.0))),
    'hollow': ((at('XY', (0, 0, 0)), (10.0, 10.0, 10.0)), 'cut', (at('XY', (3, 3, 3)), (4.0, 4.0, 4.0))),
}
CLOSED = {'hollow': lambda: mp.mpf(936)}
# This step's cavity: a turned box less a box in another turned frame inside
# it.
CAVITY_OUTER = (at('TURN30', (1, 2, 0)), (10.0, 8.0, 6.0))
CAVITY_INNER = (at('TILT', (1.75, 6.375, 3.0)), (3.0, 2.5, 2.0))


def cavity_volume():
    """The outer box less the inner one, each times its stored axes'
    determinant."""
    (fa, sa), (fb, sb) = CAVITY_OUTER, CAVITY_INNER
    return q(F(sa[0])*F(sa[1])*F(sa[2]))*det(fa)-q(F(sb[0])*F(sb[1])*F(sb[2]))*det(fb)


def make_bodies():
    out = [earlier(n) for n in ('pyramid', 'truncated', 'wedge', 'tetra', 'octa', 'notched', 'hollow')]
    ball = base.BODIES['ball']
    out.append(Body('ball', ball.path, [ball.spec], 0))
    out.append(Body('cavity', 'imported/poly_cavity.brep',
                    [polyhedra.box_spec(*CAVITY_OUTER), polyhedra.box_spec(*CAVITY_INNER)], ('cut', 0, 1),
                    closed=cavity_volume,
                    rows=[polyhedra.box_row(*CAVITY_OUTER), 'boolean cut', polyhedra.box_row(*CAVITY_INNER)]))
    return {b.name: b for b in out}


BODIES = make_bodies()
OWN = [b for b in BODIES.values() if b.rows is not None]


def imported(name):
    return ('imported', name)


def shifted(expr, k):
    if isinstance(expr, int):
        return expr+k
    return (expr[0], shifted(expr[1], k), shifted(expr[2], k))


def flatten(items):
    """The reference's solids of a case's inputs (an imported body its own
    solids) and each input's expression over them."""
    specs, exprs = [], []
    for it in items:
        if it[0] == 'imported':
            b = BODIES[it[1]]
            k = len(specs)
            specs += b.specs
            exprs.append(shifted(b.expr, k))
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
                c = base.construction(s, OPERATIONS[k])
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
    solids}`, 0 for empty); a group's declared refusal (`kind`, `reason`)
    declares its cases."""
    out = []
    for op, n in outcomes.items():
        k = kind or ('solid' if n else 'empty')
        stages = [(op, False)] if first else []
        out.append(Imported(f'{name}_{op}', klass, items, first or op, stages, k, n, reason))
    return out


PYRAMID_ROD = cylinder(at('XY', (0, 0, 0)), 4.75, 4.375, 1.25, -3.0, 9.0)
WEDGE_BALL = sphere(2.0, at('XY', (4.25, 5.5, 6.625)))
OCTA_CONE = cone(2.25, 0.75, 10.0, at('XY', (4.375, 5.75, -1.0)))
ROD_TRUNCATED = cylinder(at('TILT', (3.5, 3.75, 1.5)), 0.0, 0.0, 1.0, -6.0, 6.0)
NOTCHED_BALL = sphere(2.0, at('XY', (7.5, 8.25, 2.75)))
CHAIN_BOX = block(at('XY', (0, 0, 0)), 6.5, -1.0, 11.0, 4.25, -2.0, 8.0)
CHAIN_BALL = sphere(2.5, at('XY', (4.875, 3.75, 2.5)))
CHAIN_ROD = cylinder(at('XY', (0, 0, 0)), 5.25, 4.5, 1.25, -1.0, 9.0)
CHAIN_TILT = block(at('TILT', (6.125, 5.5, 3.0)), -2.0, -1.5, 2.0, 1.5, -1.5, 1.5)
HOLLOW_BALL = sphere(2.0625, at('XY', (8.375, 5.125, 4.875)))
# Through the cavity, parallel to its faces normal to `TILT`'s `n` (the
# cavity's `w` from 6.225 to 8.225, the slab's from 6.75 to 7.75).
CAVITY_SLAB = block(at('TILT', (0, 0, 8.4375)), -3.5, 2.5, 10.5, 17.0, 0.0, 1.0)
# A ball resting on the tetrahedron's base (its corners' heights `-1/2`
# exactly: `TURN30`'s normal is the world's `z`) at `(4, 9/2, -1/2)`.
TETRA_TOUCH = sphere(1.5, at('XY', (4.0, 4.5, 1.0)))
HOLLOW_INNER = sphere(1.25, at('XY', (5.125, 4.875, 5.25)))
TANGENT = 'a tangency between the inputs (S9c)'
SEVERAL = 'a cavity among several solids (S9c)'


def cases():
    """The cases, the groups slowest for the reference first (its workers
    take them in this order)."""
    out = []
    three = {'fuse': 1, 'cut': 1, 'common': 1}
    out += group('hollow_ball', 'cavity', [imported('hollow'), HOLLOW_BALL], three)
    out += group('cavity_slab', 'cavity', [imported('cavity'), CAVITY_SLAB], {'fuse': 1, 'cut': 2, 'common': 1})
    out += group('chain_pyramid', 'chain', [imported('pyramid'), CHAIN_BOX, CHAIN_BALL], three, first='cut')
    out += group('chain_wedge', 'chain', [imported('wedge'), CHAIN_ROD, CHAIN_TILT], {'fuse': 1, 'cut': 1, 'common': 2},
                 first='cut')
    out += group('notched_ball', 'curved', [imported('notched'), NOTCHED_BALL], three)
    out += group('tetra_ball', 'both', [imported('tetra'), imported('ball')], three)
    out += group('pyramid_rod', 'curved', [imported('pyramid'), PYRAMID_ROD], three)
    out += group('wedge_ball', 'curved', [imported('wedge'), WEDGE_BALL], three)
    out += group('octa_cone', 'curved', [imported('octa'), OCTA_CONE], three)
    out += group('rod_truncated', 'curved', [ROD_TRUNCATED, imported('truncated')], {'fuse': 1, 'cut': 2, 'common': 1})
    out += group('tetra_touch', 'curved', [imported('tetra'), TETRA_TOUCH], three, kind='degenerate', reason=TANGENT)
    out += group('hollow_inner', 'cavity', [imported('hollow'), HOLLOW_INNER], {'fuse': 2}, kind='unsupported',
                 reason=SEVERAL)
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
        if c.klass == 'cavity':
            assert any(s[0] == 'imported' and s[1] in ('hollow', 'cavity') for s in c.items), c.name
        specs, _ = flatten(c.items)
        for s in specs:
            assert s[0] in ('sphere', 'prism', 'cone', 'hull'), c.name


# ------------------------------------------------------------------ the reference

def evaluate(job):
    """One group on the chained reference with hulls: rows, results and
    checks (S9e.4b.3c.3a's, each body's closed form)."""
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
            v, _, _ = chain.measures(e)
            d = abs(v-closed())/size**3
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
    cache = os.environ.get('POLYHEDRA_CURVED_CACHE')
    if not cache:
        return evaluate(job)
    import hashlib
    import inspect
    import pickle
    from pathlib import Path
    key = hashlib.sha256((repr(job)+repr(flatten(job[1])[0])+inspect.getsource(evaluate)
                          + Path(ref.__file__).read_text()
                          + Path(ref.ch.__file__).read_text()).encode()).hexdigest()[:24]
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
# the notched box's own corners where its notch's wall meets the box's edges
# are events of two of the box's planes on that wall at one height, found
# twice a rounding apart.
SPACING_UNCHECKED = {'notched_ball'}


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
    """The writer's blocks: each of this step's own bodies' rows under
    `write NAME PATH`."""
    out = []
    for b in OWN:
        out += [f'write {b.name} {b.path}', *b.native_rows(), 'end']
    return '\n'.join(out)+'\n'


def generate(results):
    by_group = {r[0]: r for r in results}
    blocks = []
    out = [f'# case\trow ({STEP}; polyhedra_curved_boolean_reference.py: expect KIND {STEP} CLASS, reason TEXT for '
           'a degenerate or unsupported case, then result N volume area cx cy cz or empty)']
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
    """Every case's native rows (`compare_polyhedra_curved_boolean.py`)."""
    return '\n'.join(c.native() for c in all_cases())+'\n'


def case_size(case):
    """The case's size: its reference solids' reach from the origin."""
    specs, _ = flatten(case.items)
    chain = ref.Chain([construction(s, 0) for s in specs])
    return float(chain.size)


# ------------------------------------------------------------------ the files' check

def planes_of(spec):
    """A polyhedral body's planes as `(n, d)` (exact): a hull's faces', a
    prism's caps' and walls'."""
    if spec[0] == 'hull':
        I = ref.Hull(ref.HullCase(spec[1]))
        return [(tuple(F(x) for x in n), F(d)) for n, d in I.planes]
    out = []
    for s in base.construction_surfaces(spec)[1]:
        g = s.grad((mp.mpf(0),)*3)
        out.append((g, -s.f((mp.mpf(0),)*3)))
    return out


def check_bodies():
    """Each existing body file: its faces' kinds planes (but the ball's), each
    stored vertex on its construction's planes within 1e-12 of the body's
    size. Returns the largest deviation and the files read."""
    worst, read = mp.mpf(0), 0
    for b in BODIES.values():
        path = ROOT/'fixtures'/b.path
        if not path.exists():
            continue
        read += 1
        kinds, vertices = base.stored_records(path.read_text())
        if b.name == 'ball':
            assert kinds == ['sphere'], (b.name, kinds)
            continue
        assert set(kinds) == {'plane'}, (b.name, kinds)
        planes = []
        for spec in b.specs:
            planes += planes_of(spec)
        size = max([1.0]+[abs(x) for v in vertices for x in v])
        for v in vertices:
            X = tuple(mp.mpf(c) for c in v)
            d = min(abs(sum(base.tc.M(n[i])*X[i] for i in range(3))-base.tc.M(dd))
                    / mp.sqrt(sum(base.tc.M(n[i])**2 for i in range(3))) for n, dd in planes)
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
        os.environ['POLYHEDRA_CURVED_CACHE'] = args.cache
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
          'size of its construction\'s planes')


if __name__ == '__main__':
    main()
