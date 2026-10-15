#!/usr/bin/env python3
"""Fixtures for S9e.4b.4d.2c.2 of REVIEW_NOTES.md: imported turned bodies
given to Booleans whose partners cross their smooth edges (the tangent
circles where a fillet or a round meets its plane, cylinder or cone),
among them the known failure of "S9e.4b.4d.2c.1a's consistency fixed" (a
rib whose side planes cross a round's tangent circle, on a body of
`OCC485a`'s form with a fillet and a round).

The bodies are OCCT's own output, as S9e.4b.4d.2c.1a's and 2c.1b's: this
step's own, `ribbed`, is the one `write ribbed imported/crossing_ribbed.brep`
block of `boolean-turned-crossings-bodies.txt` for `occt_boolean_oracle.cpp`
(a `plane` row, the profile's `wire` of lines and arcs, a `revolve` row
turning it a whole turn about the world's `z`, and a `tolerance` row setting
every tolerance to 3e-3, `OCC485a`'s), which
`compare_turned_crossings_boolean.py --write-bodies` runs to write
`rust/fixtures/imported/crossing_ribbed.brep`; the others are the earlier
steps' files, read as they are: S9e.4b.4d.2c.1a's `fillet_step` and
`fillet_collar` (`generate_turned_fillets_boolean_fixtures.py`) and
2c.1b's `cone_fillet_round` and `cone_fillet_fillet`
(`generate_turned_cone_fillets_boolean_fixtures.py`). Nothing here reads
the kernel; the files are read only to check that each is the body it
claims to be. The cases' files are 2c.1a's, under the prefix
`boolean-turned-crossings`, the class `pieces`.

`ribbed` (`OCC485a`'s form on exact data at its resolution, its profile
`(r, w)` about the world's `z` from the axis out and back): a bottom disc of
radius 24, the cylinder to 32, a chamfer in to 21 at 35, the step's annulus
in to 10, a fillet of radius 1/2 about `(10, 71/2)` to the shaft of radius
19/2, the shaft to 41, a round of radius 3 about `(13/2, 41)` to the top's
disc at 44; OCCT's profile is the model's (every number dyadic).

Cases (each the three operations), every partner crossing a smooth edge
transversally:
`ribbed_rib`, the known failure: the rib over the triangle `(8, 69/2)`,
`(20, 69/2)`, `(8, 42)` in `(x, z)` across `|y| <= 1` (its side planes
across the fillet's two tangent circles and the round's lower one, its
hypotenuse across the round's torus, its top at 42 over the shaft's end at
41); `ribbed_post`, a rod parallel to the axis about `(77/8, 0)` of radius
5/8 from 33 to 47 (across the step's plane, the fillet, the shaft and the
round, the fillet's two circles and the round's lower one); `ribbed_box`,
the box `[8, 11] x [-1, 1] x [34, 43]` (across the same circles, the round's
torus by its top); `ribbed_skew`, a cube of side 2 on the skew frame `SKEW`
about the round's lower circle's point `(19/2, 0, 41)` (none of its faces
parallel or normal to the axis, so that it evaluates turned); `step_ring`, `fillet_step` and a rod along the world's
`x` about `(y, z) = (0, 225/32)` of radius 1/4 (across the round's lower
circle, the shaft's end at 7); `collar_box`, `fillet_collar` and the box
`[5/2, 13/2] x [-1/2, 1/2] x [15/4, 15/2]` (across the fillet's two circles and
the upper round's lower one); `round_rib`, `cone_fillet_round` and the rib
over `(3, 3)`, `(8, 3)`, `(3, 15/2)` across `|y| <= 1/2` (its side planes
across the round's circles on the cylinder and on the cone); `fillet_rod`,
`cone_fillet_fillet` and a rod along `x` about `(y, z) = (0, 81/32)` of
radius 1/4 (across the fillet's circle on its cone, at `43/17`). Declared
`degenerate`: `step_touch`, `fillet_step` and the box `[7/2, 9/2] x [-1/2,
1/2] x [5/2, 5]`, its side `x = 9/2` touching the fillet's circle on the
step (of radius 9/2 at 3) at `(9/2, 0, 3)` (the plane through the circle's
tangent there: `Degenerate("a tangency between the inputs (S9c)")`, contact
only; the same box on `ribbed`, its side at the fillet's circle of radius
10, took the reference past 20 minutes).

The reference is 2c.1a's and 2c.1b's: S9e.3a's chained reference with each
body one input turned from its model profile (`TurnedInput`, its faces the
profile's own, a cone's band among them), S9e.4b.3c.3a's checks and each
body's closed form by Pappus over its model profile. With `--check`, where
the bodies' files exist, each is read (`stored_records`): its faces' kinds
its OCCT profile's, every stored vertex within 1e-12 of the size on the OCCT
profile's lines and circles.
"""
import argparse
from fractions import Fraction as F
import os
import struct
import zlib

import mpmath as mp

from identity_reference import Boundary, Case
from curve_surface_reference import stored_axes
import generate_turned_fillets_boolean_fixtures as fillets
import generate_turned_cone_fillets_boolean_fixtures as cones

base = fillets.base
forms = fillets.forms
ROOT = base.ROOT
OPS = base.OPS
OPERATIONS = base.OPERATIONS
CLASSES = ('pieces',)
PREFIX = 'boolean-turned-crossings'
BODIES_FILE = f'{PREFIX}-bodies.txt'
STEP = 'S9e.4b.4d.2c.2'
construction = base.construction
ref = base.ref
q = forms.q
cylinder, block, prism = forms.cylinder, forms.block, base.prism
at = forms.at
L, arc = fillets.L, fillets.arc
flatten = fillets.flatten
validate = fillets.validate


# ------------------------------------------------------------------ bodies

class Model(fillets.Model):
    """S9e.4b.4d.2c.1a's model profile (each torus between a plane and a
    cylinder; cones elsewhere at sharp corners), with each corner's kite
    (its centre, tangent points and corner) for 2c.1b's binary64
    membership."""

    def kite(self, c):
        (cr, cw), (pr, pw) = c['centre'], c['corner']
        # The tangent points: the centre's feet on the corner's two lines
        # (a plane's at the corner's height, a cylinder's at its radius).
        return [(cr, cw), (pr, cw), (pr, pw), (cr, pw)]


class Body(fillets.Body):
    """This step's own body about the world's `z`: OCCT's profile the
    model's."""

    def __init__(self, name, profile, tolerance=None):
        self.name, self.axis, self.origin = name, 'z', (F(0), F(0), F(0))
        self.path = f'imported/crossing_{name}.brep'
        self.occt, self.tolerance = profile, tolerance
        self.model = Model(profile)
        self.lattice = self.model.lattice
        self.ref_specs = [('turned', name)]
        self.ref_expr = 0


def make_ribbed():
    return [(0, 0, L), (24, 0, L), (24, 32, L), (21, 35, L), (10, 35, arc(10, F(71, 2), F(1, 2), False)),
            (F(19, 2), F(71, 2), L), (F(19, 2), 41, arc(F(13, 2), 41, 3, True)), (F(13, 2), 44, L), (0, 44, L)]


OWN = {'ribbed': Body('ribbed', make_ribbed(), tolerance='0.003')}
# Every body of the cases: this step's and the earlier steps' files.
BODIES = dict(OWN, step=fillets.BODIES['step'], collar=fillets.BODIES['collar'],
              round=cones.BODIES['round'], fillet=cones.BODIES['fillet'])


def imported(name):
    return ('imported', name)


class Imported(fillets.Imported):
    """One case, its imported input one of `BODIES`."""

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


def group(name, items, outcomes, kind=None, reason=None):
    """Cases of one group: the operation varied over `outcomes` (`{op:
    solids}`, 0 for empty); the group's own refusal declares its cases."""
    out = []
    for op, n in outcomes.items():
        k = kind or ('solid' if n else 'empty')
        out.append(Imported(f'{name}_{op}', 'pieces', items, op, [], k, n, reason if kind else None))
    return out


def rib(points, y0, y1):
    """The prism over a triangle in `(x, z)` (`points` as `(x, z)`, counter-
    clockwise in `(z, x)`) across `y0 <= y <= y1`: on the frame `ALONGY`
    (its normal `y`, its `x` the world's `z`, its `y` the world's `x`)."""
    return prism([Boundary(points=[(float(z), float(x)) for x, z in points])], at('ALONGY', (0, 0, 0)),
                 float(y0), float(y1))


# The known failure: its side planes across the fillet's circles (radius 10
# at 35, 19/2 at 71/2) and the round's lower one (19/2 at 41), its
# hypotenuse from (20, 69/2) to (8, 42) across the round's torus.
RIBBED_RIB = rib([(8, F(69, 2)), (8, 42), (20, F(69, 2))], -1, 1)
# A rod parallel to the axis about (77/8, 0) of radius 5/8 (radii 9 to
# 41/4 from the axis): across the step's plane beyond 10, the fillet, the
# shaft and the round below its radius 19/2.
RIBBED_POST = cylinder(at('XY', (0, 0, 33)), 9.625, 0.0, 0.625, 0.0, 14.0)
RIBBED_BOX = block(at('XY', (0, 0, 0)), 8.0, -1.0, 11.0, 1.0, 34.0, 43.0)
# A cube of side 2 on the skew frame `SKEW` (no face parallel or normal to
# the axis) about the ribbed body's round's lower circle's point (19/2, 0,
# 41): across the shaft, that circle and the round.
RIBBED_SKEW = block(at('SKEW', (9.5, 0, 41)), -1.0, -1.0, 1.0, 1.0, -1.0, 1.0)
# The rod along x across the step body's round's lower circle (radius 4 at
# 7), its axis 1/32 above it (at the circle's height every point of the
# rod's generatrix through it is both the shaft's and the round's event, a
# rounding apart for the reference).
STEP_RING = cylinder(at('SIDE', (-10, 0, 0)), 0.0, 7.03125, 0.25, 0.0, 20.0)
# The collar's box: its sides y = +-1/2 across the fillet's circles (radius
# 3 at 9/2, 7/2 at 5) and the upper round's lower one (6 at 7).
COLLAR_BOX = block(at('XY', (0, 0, 0)), 2.5, -0.5, 6.5, 0.5, 3.75, 7.5)
# The cone round's rib: its sides across the round's circles on the
# cylinder (6 at 4) and on the cone (98/17 at 84/17).
ROUND_RIB = rib([(3, 3), (3, F(15, 2)), (8, 3)], F(-1, 2), F(1, 2))
# The rod along x about (0, 81/32) across the cone fillet's circle on its
# cone (70/17 at 43/17).
FILLET_ROD = cylinder(at('SIDE', (-10, 0, 0)), 0.0, 2.53125, 0.25, 0.0, 20.0)
# The box whose side x = 9/2 touches the step body's fillet's circle on its
# step (radius 9/2 at 3) at (9/2, 0, 3).
STEP_TOUCH = block(at('XY', (0, 0, 0)), 3.5, -0.5, 4.5, 0.5, 2.5, 5.0)
TOUCH = 'a tangency between the inputs (S9c)'


def cases():
    """The cases, the groups slowest for the reference first (its workers
    take them in this order)."""
    three = {'fuse': 1, 'cut': 1, 'common': 1}
    out = []
    out += group('ribbed_post', [imported('ribbed'), RIBBED_POST], three)
    out += group('ribbed_rib', [imported('ribbed'), RIBBED_RIB], three)
    out += group('ribbed_box', [imported('ribbed'), RIBBED_BOX], three)
    out += group('ribbed_skew', [imported('ribbed'), RIBBED_SKEW], three)
    out += group('collar_box', [imported('collar'), COLLAR_BOX], three)
    out += group('step_ring', [imported('step'), STEP_RING], three)
    out += group('fillet_rod', [imported('fillet'), FILLET_ROD], three)
    out += group('round_rib', [imported('round'), ROUND_RIB], three)
    out += group('step_touch', [imported('step'), STEP_TOUCH], {'fuse': 1, 'cut': 1, 'common': 1},
                 kind='degenerate', reason=TOUCH)
    return out


def all_cases():
    return cases()


# ------------------------------------------------------------------ the reference's turned input

# 2c.1b's chain (its own bodies' and 2c.1a's inputs): this step's own body
# 2c.1b's `TurnedInput` (a cone's band among its faces), its binary64 form
# 2c.1b's `FloatTurned` (by `make_float`, which takes 2c.1b's inputs).
_make_input = ref.make_input


def make_input(case):
    body = getattr(case, 'turned', None)
    return cones.TurnedInput(body) if isinstance(body, Body) else _make_input(case)


ref.make_input = make_input
TurnedCase = fillets.TurnedCase


def reference_case(spec):
    return TurnedCase(BODIES[spec[1]]) if spec[0] == 'turned' else construction(spec, 0)


# ------------------------------------------------------------------ the reference

def evaluate(job):
    """One group on the chained reference: rows, results and checks
    (S9e.4b.3c.3a's, each body's closed form)."""
    name, items, ops, declared, mc_n = job
    specs, sub = flatten(items)
    chain = ref.Chain([reference_case(s) for s in specs])
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
    # Solid counts by rays (not the declared group's: its box touches the
    # fillet's circle at a point).
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
    sine, _, _ = chain.margins()
    margins = {'sine': sine, 'spacing': base.spacing(chain)}
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
    cache = os.environ.get('TURNED_CROSSINGS_CACHE')
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


LIMITS = fillets.LIMITS
reference_checks = fillets.reference_checks


# ------------------------------------------------------------------ the files

def bodies_text():
    """The writer's blocks: this step's own bodies' rows under `write NAME
    PATH` (the earlier steps' files are theirs)."""
    out = []
    for b in OWN.values():
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
                frames.append(f'{case.name}\tsolid{k} {key}\t'+' '.join(struct.pack('>d', q_).hex() for q_ in v))
    return {f'{PREFIX}-cases.txt': '\n'.join(blocks)+'\n',
            f'{PREFIX}-expected.tsv': '\n'.join(out)+'\n',
            f'{PREFIX}-frames.tsv': '\n'.join(frames)+'\n',
            BODIES_FILE: bodies_text()}


def native_input():
    """Every case's native rows (`compare_turned_crossings_boolean.py`)."""
    return '\n'.join(c.native() for c in all_cases())+'\n'


def case_size(case):
    """The case's size: its reference solids' reach from the origin."""
    specs, _ = flatten(case.items)
    chain = ref.Chain([reference_case(s) for s in specs])
    return float(chain.size)


# ------------------------------------------------------------------ the files' check

def check_bodies():
    """Each existing body file: its faces' kinds its OCCT profile's, each
    stored vertex on that profile's lines and circles within 1e-12 of the
    body's size. Returns the largest deviation and the files read."""
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
            worst = max(worst, fillets.profile_distance(b, X)/size)
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
        os.environ['TURNED_CROSSINGS_CACHE'] = args.cache
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
          'size of its profile\'s lines and circles')


if __name__ == '__main__':
    main()
