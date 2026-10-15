#!/usr/bin/env python3
"""Fixtures for S9e.4b.4d.2b of REVIEW_NOTES.md: imported turned bodies of
cylinders, cones and planes about one axis whose steps or ends are cones
within the resolution of a plane (`OCC497c`'s and `OCC497e`'s kind: a step
or an end turned a few millionths from flat, an apex dimpling or bulging a
disc), each such cone taken as the plane through its wider rim, given to
Booleans.

The bodies are OCCT's own output, as S9e.4b.4d.2a's
(`generate_turned_profiles_boolean_fixtures.py`): `boolean-turned-flats-
bodies.txt` holds one `write NAME imported/flat_NAME.brep` block per body for
`occt_boolean_oracle.cpp` (a `plane` row, the profile's `wire` of lines, a
`revolve` row turning it a whole turn about the axis, and a `tolerance` row
setting every tolerance to 1e-3 as a translated part stores its own), which
`compare_turned_flats_boolean.py --write-bodies` runs to write
`rust/fixtures/imported/flat_NAME.brep` (format version 1, no
triangulations). Nothing here reads the kernel; the files are read only to
check that each is the body it claims to be (below). OCCT's revolve stores a
step rising less than 1e-7 as a plane, and a wall widening by up to 2^-10
over 6 (1.6e-4 rad) as a cylinder (a cone within the resolution of a
cylinder does not arise from it, so no body here has one).

`boolean-turned-flats-cases.txt` lists each case in the Boolean protocol, an
imported input as its one `brep imported/flat_NAME.brep` row;
`boolean-turned-flats-expected.tsv` gives per case `expect KIND S9e.4b.4d.2b
CLASS` (`solid`, `empty`, `degenerate`: a `Degenerate` of S9's rules; the
class `pieces`: an imported body and a construction), a degenerate case's
`reason TEXT`, then `result N volume area cx cy cz` (totals over the N
solids, world coordinates) or `empty`; `boolean-turned-flats-frames.tsv` the
stored axes of every solid the kernel builds from its rows.

Each body has two profiles of dyadic points `(r, w)` (`r` the radius, `w`
the height along the axis), from the axis out and back to it: the one OCCT
is given (`occt`), whose near-flat segments rise or fall by `D` = 2^-17
(7.6e-6) over 2 to 6 (a cone 1.3e-6 to 3.8e-6 rad from flat, 0.76 % of the
resolution), and the model's (`model`), each such segment taken level at its
wider end's height (a plane), its neighbours' lines kept and met there: the
lattice of S9e.4b.4d.2a over the model's profile. The bodies, each about the
world's `z` through the origin: `pin` (`OCC497c`'s kind): a bottom whose
apex stands `D` above its rim of radius 6 (a dimple, taken as the disc at
height 0, the apex a vertex within the resolution of its centre), a cylinder
to 4, a step in to radius `4 - D/8` rising `D` (taken as the annulus at 4),
a draft cone narrowing by 1/8 to radius 3 at 12 and its disc; `collar`
(`OCC497e`'s kind): a bottom whose apex lies `D` below its rim of radius 3
(a bulge, taken as the disc at `D`), a cylinder to 6, a step out to radius 6
falling `D` (taken as the annulus at `6 - D` facing down), a cylinder to 10,
a chamfer narrowing by 1 to radius 4 at 12 and its disc. `rod` (no Boolean):
a cylinder of radius 2 from a disc at 0 to height 5 and a top cone rising
`D` to its apex: one curved face beside its flat, no lattice, its import
refused as S9e.4b.4d.1's cone within the resolution of a plane.

Cases (each the three operations): `pin_box`, the pin and a box across its
step and its draft cone; `pin_ball`, the pin and a ball across its bottom
about the dimple's apex; `collar_rod`, the collar and a rod along the
world's `x` across its two cylinders and its step out; and declared
`degenerate`, `pin_rest`, the pin and a box over its step 2^-12 above it
(within the resolution of the step's plane: `bug497_2`'s configuration,
`OCC497d`'s floor 1.9e-5 above `OCC497c`'s step; `Degenerate("two faces
within the resolution of one plane")`).

The reference is the model's construction through S9e.3a's chained
reference (`chained_curved_boolean_reference.py`): each body the max-min
form of its model profile's lines (S9e.4b.4d.2a's `Lattice`), then `(X) op
C`, with S9e.4b.3c.3a's checks (the two families, each solid's closed form,
the pair identities, Monte Carlo, each body one solid, and outside the
declared group solid counts by rays at two resolutions, the meetings' sines
and the events' spacing) and each body's closed form by Pappus over its
model profile. With `--check`, where the bodies' files exist, each is read
(`stored_records`): its faces' kinds its OCCT profile's, every stored vertex
within 1e-12 of the size on the OCCT profile's surfaces.
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
PREFIX = 'boolean-turned-flats'
BODIES_FILE = f'{PREFIX}-bodies.txt'
STEP = 'S9e.4b.4d.2b'
prism, sphere, cone, construction = base.prism, base.sphere, base.cone, base.construction
ref = base.ref
q = forms.q
cylinder, block = forms.cylinder, forms.block
Lattice, num, at, TURN = profiles.Lattice, profiles.num, profiles.at, profiles.TURN
# A near-flat segment's rise.
D = F(1, 2**17)
# The bodies' resolution (the oracle's `tolerance` row).
RESOLUTION = '0.001'


# ------------------------------------------------------------------ bodies

def near(a, b):
    """Whether a segment is near flat: within `4 D` of level, not level."""
    (ra, wa), (rb, wb) = a, b
    return wa != wb and abs(wb-wa) <= 4*D and ra != rb


class Body:
    """An imported turned body about the world's `z` through `origin`: OCCT's
    profile, the model's (its lattice, `None` for a body that is none), the
    relation between them checked."""

    def __init__(self, name, occt, model, origin=(0, 0, 0)):
        self.name, self.axis, self.origin = name, 'z', tuple(F(c) for c in origin)
        self.path = f'imported/flat_{name}.brep'
        self.occt = [(F(r), F(w)) for r, w in occt]
        self.model = [(F(r), F(w)) for r, w in model] if model else None
        if self.model:
            self.relate()
            self.lattice = Lattice(self.model)
            self.ref_specs = self.lattice.specs(self.frame_at)
            self.ref_expr = self.lattice.expr()

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
        return (o[0]+F(r), o[1], o[2]+F(w))

    def frame_at(self, w):
        return at('XY', self.place(0, 0 if w is None else w))

    def native_rows(self):
        o = self.origin
        wire = [f'wire S {len(self.occt)}']
        for r, w in self.occt:
            wire.append(' '.join(num(c) for c in self.place(r, w)))
            wire.append('L')
        return [f'plane {num(o[0])} {num(o[1])} {num(o[2])} 0 -1 0 1 0 0', ' '.join(wire),
                f'revolve {num(o[0])} {num(o[1])} {num(o[2])} 0 0 1 {TURN}', f'tolerance {RESOLUTION}']

    def closed(self):
        return self.lattice.closed()

    def kinds(self):
        """Its stored faces' surface kinds: OCCT's profile's segments'."""
        out = []
        for (ra, wa), (rb, wb) in zip(self.occt, self.occt[1:]):
            out.append('plane' if wa == wb else ('cylinder' if ra == rb else 'cone'))
        return sorted(out)


def make_bodies():
    out = [
        Body('pin', [(0, D), (6, 0), (6, 4), (4-D/8, 4+D), (3, 12), (0, 12)],
             [(0, 0), (6, 0), (6, 4), (4, 4), (3, 12), (0, 12)]),
        Body('collar', [(0, 0), (3, D), (3, 6), (6, 6-D), (6, 10), (4, 12), (0, 12)],
             [(0, D), (3, D), (3, 6-D), (6, 6-D), (6, 10), (4, 12), (0, 12)]),
        Body('rod', [(0, 0), (2, 0), (2, 5), (0, 5+D)], None),
    ]
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
            exprs.append(b.lattice.expr(k))
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


# The pin's box across its step (radii 4 to 6 at height 4) and its draft
# cone, from the cylinder below.
PIN_BOX = block(at('TILT', (2.5, -0.75, 2.5)), 0.0, 0.0, 3.0, 2.0, 0.0, 3.5)
# A ball across the pin's bottom about the dimple's apex (the origin).
PIN_BALL = sphere(2.0, at('XY', (0.5, 0.25, -1.0)))
# The collar's rod along the world's x about (y, z) = (3/4, 21/4), across
# its cylinder of radius 3 (below 6 - D), its step out and its cylinder of
# radius 6.
COLLAR_ROD = cylinder(at('SIDE', (-8, 0.75, 0)), 0.0, 5.25, 1.25, 0.0, 16.0)
# A box over the pin's step, its floor 2^-12 above it.
PIN_REST = block(at('XY', (0, 0, 0)), 4.5, -1.0, 7.5, 1.0, 4+2.0**-12, 7.0)
REST = 'two faces within the resolution of one plane'


def cases():
    """The cases, the groups slowest for the reference first (its workers
    take them in this order)."""
    three = {'fuse': 1, 'cut': 1, 'common': 1}
    out = []
    out += group('collar_rod', [imported('collar'), COLLAR_ROD], three)
    out += group('pin_box', [imported('pin'), PIN_BOX], three)
    out += group('pin_ball', [imported('pin'), PIN_BALL], three)
    out += group('pin_rest', [imported('pin'), PIN_REST], {'fuse': 2, 'cut': 1, 'common': 0},
                 kind='degenerate', reason=REST)
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
    # Solid counts by rays (not the declared group's: its box's gap of
    # 2^-12 lies below the rays' spacing).
    for op, (expr, _) in exprs.items():
        for n in () if declared else (48, 71):
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
    cache = os.environ.get('TURNED_FLATS_CACHE')
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


LIMITS = profiles.LIMITS


def reference_checks(results):
    worst, covered = {}, {}

    def note(key, value):
        worst[key] = max(worst.get(key, mp.mpf(0)), value)
        covered[key] = covered.get(key, 0)+1

    for name, rows, _, checks, margins, stats, declared in results:
        # The declared group's box rests within the resolution of the step:
        # its checks kept apart, looser.
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
    """Every case's native rows (`compare_turned_flats_boolean.py`)."""
    return '\n'.join(c.native() for c in all_cases())+'\n'


def case_size(case):
    """The case's size: its reference solids' reach from the origin."""
    specs, _ = flatten(case.items)
    chain = ref.Chain([construction(s, 0) for s in specs])
    return float(chain.size)


# ------------------------------------------------------------------ the files' check

def profile_distance(body, X):
    """A world point's distance from the body's OCCT profile's surfaces (each
    segment's cone, cylinder or plane, unbounded), in its axis frame."""
    o = [q(c) for c in body.origin]
    d = [X[i]-o[i] for i in range(3)]
    r, w = mp.sqrt(d[0]**2+d[1]**2), d[2]
    best = None
    for (ra, wa), (rb, wb) in zip(body.occt, body.occt[1:]):
        a, b = q(rb-ra), -q(wb-wa)
        dist = abs(a*(w-q(wa))+b*(r-q(ra)))/mp.sqrt(a*a+b*b)
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
            worst = max(worst, profile_distance(b, X)/size)
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
        os.environ['TURNED_FLATS_CACHE'] = args.cache
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
