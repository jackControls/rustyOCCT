#!/usr/bin/env python3
"""Fixtures for S9e.3b of REVIEW_NOTES.md: a Boolean's result given to
another Boolean whose solid's faces meet the given result's edges on
meetings of two curved faces (a cylinder's and a sphere's `Rise`, two
cylinders' or a cone's and a cylinder's `Meet`, a torus's and a cylinder's
`Toric`) or on a plane's general section of a cone or of a torus: three
surfaces, two or three of them curved.

`boolean-given-met-cases.txt` lists each case in the Boolean protocol with
its `then` row (`identity_reference.encode_chained_case`);
`boolean-given-met-expected.tsv` gives per case, from
`chained_curved_boolean_reference.py` (S9e.3a's reference, unchanged: it
decides the chain by its solids' memberships):

* `expect KIND S9e.3b CLASS PARTNER`: the declared outcome (`solid`: one or
  more solids; `empty`; `degenerate`: `Degenerate` in the decisions), the
  given edges' class (`rise`, `meet`, `toric`, `cone`: a cone's section,
  `spiric`: a torus's section) and the partner's surface met there
  (`plane`, `cylinder`, `sphere`, `cone`), then for a degenerate case
  `reason TEXT`;
* `meet N SINE`: the triple points found on the given result's edges inside
  the partner's faces (`given_met_reference.py`) and their least sine;
* `result N volume area cx cy cz` (totals over the N solids, world
  coordinates) or `empty`.

`boolean-given-met-frames.tsv` records every solid's stored axes. Frames
are the curved generator's whose stored axes the kernel gives bit for bit
(`XY`, `SIDE`, `TILT`).

The chains (the given result first, then the partner): a sphere of radius 5
about the origin fused with a peg of radius 1 about `(0, 2)` (their meeting
a `Rise` loop from `(0, 3, 4)` up to `(0, 1, sqrt 24)`), with a `TILTX` slab
(normal `(0, -4, 3) / 5`) across the loop (`peg_tilt`), a box whose wall `y
= 2` holds the peg's axis (parallel to its rulings, `peg_wall`, swapped), a
pipe along `x` (`peg_pipe`), a ball (`peg_ball`) and a frustum along `x`
(`peg_cone`) across it; a rod of radius 2 along `x` fused with a rod of
radius 1 along `z` off its axis by 1/2 (two `Meet` rings over the thinner
rod), with a `TILT` slab across the upper ring (`cross_tilt`), a box whose
wall `x = 0` holds the thin rod's axis (`cross_wall`) and a ball across the
ring (`cross_ball`); a frustum of radii 3 and 1 fused with a pipe along `x`
through its wall (a cone's and a cylinder's `Meet` loops), with a wall
across a loop (`frustum_pipe_wall`); a torus of radii 3 and 1 fused with a
rod of radius 1/2 through its tube (two `Toric` loops), with a `TILTX` slab
across the upper loop (`torus_rod_tilt`), and with the rod about `(3.1, 0)`
and a wall holding its axis (`torus_rod_wall`); the frustum less a `TILT`
box above `w = 2.2` (an elliptic section of the cone), with a rod
(`cone_cut_rod`) and a ball (`cone_cut_ball`) across the section; the torus
less a box beyond the wall `x = 5/2` (a spiric section), with a pipe along
`x` (`torus_cut_pipe`, its common two solids) and a ball (`torus_cut_ball`)
across the section; declared `degenerate`: a box whose bottom face `z = 4`
touches the peg's loop at its lowest point `(0, 3, 4)` (`peg_touch`), and a
ball of radius 13/16 about `(0, 15/4, 69/16)` through that point, whose
normal there, `(0, -12, -5) / 13`, is dependent on the sphere's and the
peg's though no two of the three surfaces touch (`peg_kiss`, fuse and cut:
the common's pieces touch at the point, which rays count unsteadily). The
slabs end just past their chains' solids (Monte Carlo samples the solids'
common box).

Before writing, the chained reference's checks run as S9e.3a's
(`generate_given_curved_boolean_fixtures.evaluate`: the two families of
curves of every face, each solid's closed form, the pair identities for the
given result and the partner, the area identity, the given result one solid,
the declared solid counts by rays at two resolutions, Monte Carlo), and
`given_met_reference.py`'s triple points: every chain outside the declared
cases holds at least one on the given result's edges of its class, every one
found at least twice (from both solids' faces or both families), every
sine at least 0.05; each declared case's tangency at its declared point
(the three surfaces through it within 1e-30, their normals' determinant
within 1e-30, on the given result's edge and the partner's face). The
declared tangencies' sections touch on the faces through the point, so for
them alone the chained reference takes roots of different surfaces within
1e-30 of a curve's range as one (`MERGE`: else the noise-classified piece
between them flips the structure, an event found again and again) and
integrates to 1e-22 of the size to the fourth (`QUAD`: the touching
pieces' endpoint singularities converge slowly), its checks there at 1e-18
(their rows are rounded to binary64; the kernel refuses them). Chains run in
worker processes (`--workers`).
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
import os
from pathlib import Path
import struct

import mpmath as mp

from curve_surface_reference import stored_axes
import chained_curved_boolean_reference as ref
from generate_curved_boolean_fixtures import FRAMES, at, disc, square
import generate_given_curved_boolean_fixtures as gc
from generate_given_curved_boolean_fixtures import cone, prism, sphere, torus, frame_name
import given_met_reference as gm
from identity_reference import encode_chained_case

ROOT = Path(__file__).resolve().parents[1]
CLASSES = ('rise', 'meet', 'toric', 'cone', 'spiric')
PARTNERS = ('plane', 'cylinder', 'sphere', 'cone')
OPS = gc.OPS


class Met:
    """One case: the solids, the first Boolean and the last one's operation
    (`swapped`: the partner the object)."""

    def __init__(self, name, klass, partner, specs, op1, op2, swapped, kind, solids, reason=None, touch=None):
        assert klass in CLASSES and partner in PARTNERS, (klass, partner)
        self.name, self.klass, self.partner = name, klass, partner
        self.op1, self.op2, self.swapped = op1, op2, swapped
        self.kind, self.solids, self.reason, self.touch = kind, solids, reason, touch
        self.chain_name = name.rsplit('_', 1)[0]
        self.specs = specs
        self.solid_cases = [gc.solid(s, gc.OPERATIONS[k]) for k, s in enumerate(specs)]
        for c in self.solid_cases:
            c.name = name
        self.frames = [gc.spec_frame(s) for s in specs]
        self.stages = [(op2, swapped)]

    @property
    def last(self):
        return self.op2

    def expr(self):
        return ref.chain_expr(self.op1, self.op2, self.swapped)

    def encode(self):
        c = self.solid_cases
        return encode_chained_case(c[0], self.op1, c[1], gc.BOOLEANS[0], self.op2, c[2], gc.BOOLEANS[1],
                                   self.swapped, None, [])


def group(name, klass, partner, specs, op1, outcomes, swapped=False, reason=None, touch=None):
    """Cases of one chain, the last Boolean's operation over `outcomes`
    (`{op: solids}`, 0 for empty)."""
    out = []
    for op, n in outcomes.items():
        kind = 'degenerate' if reason else ('solid' if n else 'empty')
        out.append(Met(f'{name}_{op}', klass, partner, specs, op1, op, swapped, kind, n, reason, touch))
    return out


# The sphere and the peg: their meeting a Rise loop from (0, 3, 4) (where
# the peg's circle is farthest from the sphere's axis) up to (0, 1, sqrt 24).
BALL5 = sphere(5.0, at('XY', (0, 0, 0)))
PEG = prism([disc(0.0, 2.0, 1.0)], at('XY', (0, 0, 0)), 0.0, 8.0)
# Two rods: radius 2 along x, radius 1 along z about (0, 1/2): two rings
# over the thinner rod's angle.
ROD_X = prism([disc(0.0, 0.0, 2.0)], at('SIDE', (0, 0, 0)), -6.0, 6.0)
ROD_Z = prism([disc(0.0, 0.5, 1.0)], at('XY', (0, 0, 0)), -6.0, 6.0)
# A frustum of radii 3 and 1 and height 4.
FRUSTUM = cone(3.0, 1.0, 4.0, at('XY', (0, 0, 0)))
# A torus of radii 3 and 1 about the z axis.
RING = torus(3.0, 1.0, at('XY', (0, 0, 0)))
TOUCH = ('a box\'s bottom face z = 4 tangent to the given sphere\'s and peg\'s meeting at its lowest point '
         '(0, 3, 4)')
KISS = ('a ball of radius 13/16 about (0, 15/4, 69/16) through the given meeting\'s lowest point (0, 3, 4), its '
        'normal there (0, -12, -5) / 13 dependent on the sphere\'s and the peg\'s (the ball tangent to the meeting, '
        'no two of the surfaces tangent)')


def big_box(x0, y0, z0, x1, y1, z1):
    return prism([square(x0, y0, x1, y1)], at('XY', (0, 0, 0)), z0, z1)


def tilt_slab(w0, w1, frame='TILT', across=(-8.0, -8.0, 8.0, 8.0)):
    """A `TILT` box over `w` in `[w0, w1]` along the frame's normal `(0, 3,
    4) / 5` (`TILTX`'s `(0, -4, 3) / 5`), over the rectangle `across` of its
    profile (each chain's just past its solids: Monte Carlo samples the
    solids' common box)."""
    return prism([square(*across)], at(frame, (0, 0, 0)), w0, w1)


# The frustum less a TILT box across it above w = 2.2 (an elliptic section).
FRUSTUM_CUT = tilt_slab(2.2, 5.5, 'TILT', (-4.0, -5.0, 4.0, 5.0))


def cases():
    out = []
    out += group('peg_tilt', 'rise', 'plane', [BALL5, PEG, tilt_slab(1.0, 6.0, 'TILTX', (-6.0, -6.0, 6.0, 9.0))],
                 'fuse', {'fuse': 1, 'cut': 1, 'common': 1})
    out += group('peg_wall', 'rise', 'plane', [BALL5, PEG, big_box(-8.0, 2.0, -8.0, 8.0, 8.0, 10.0)], 'fuse',
                 {'fuse': 1, 'cut': 1, 'common': 1}, swapped=True)
    out += group('peg_pipe', 'rise', 'cylinder',
                 [BALL5, PEG, prism([disc(2.9, 4.1, 0.5)], at('SIDE', (0, 0, 0)), -7.0, 7.0)], 'fuse',
                 {'fuse': 1, 'cut': 1, 'common': 1})
    out += group('peg_ball', 'rise', 'sphere', [BALL5, PEG, sphere(1.0, at('XY', (0.25, 2.75, 4.25)))], 'fuse',
                 {'fuse': 1, 'cut': 1, 'common': 1})
    out += group('peg_cone', 'rise', 'cone', [BALL5, PEG, cone(1.2, 0.4, 3.0, at('SIDE', (-1.5, 3.2, 4.2)))],
                 'fuse', {'fuse': 1, 'cut': 1, 'common': 1})
    out += group('cross_tilt', 'meet', 'plane', [ROD_X, ROD_Z, tilt_slab(1.5, 6.5, 'TILT', (-7.0, -6.0, 7.0, 6.0))],
                 'fuse', {'fuse': 1, 'cut': 1, 'common': 1})
    out += group('cross_wall', 'meet', 'plane', [ROD_X, ROD_Z, big_box(0.0, -8.0, -8.0, 8.0, 8.0, 8.0)], 'fuse',
                 {'fuse': 1, 'cut': 1, 'common': 1})
    out += group('cross_ball', 'meet', 'sphere', [ROD_X, ROD_Z, sphere(0.6, at('XY', (0.95, 0.6, 1.9)))], 'fuse',
                 {'fuse': 1, 'cut': 1, 'common': 1})
    out += group('frustum_pipe_wall', 'meet', 'plane',
                 [FRUSTUM, prism([disc(0.5, 2.0, 0.6)], at('SIDE', (0, 0, 0)), -5.0, 5.0),
                  big_box(1.85, -8.0, -8.0, 8.0, 8.0, 8.0)], 'fuse', {'fuse': 1, 'cut': 1, 'common': 1})
    out += group('torus_rod_tilt', 'toric', 'plane',
                 [RING, prism([disc(3.0, 0.0, 0.5)], at('XY', (0, 0, 0)), -3.0, 3.0),
                  tilt_slab(0.3, 5.5, 'TILTX', (-5.0, -6.0, 5.0, 6.0))], 'fuse', {'fuse': 1, 'cut': 1, 'common': 1})
    out += group('torus_rod_wall', 'toric', 'plane',
                 [RING, prism([disc(3.1, 0.0, 0.5)], at('XY', (0, 0, 0)), -3.0, 3.0),
                  big_box(3.1, -8.0, -4.0, 8.0, 8.0, 4.0)], 'fuse', {'fuse': 1, 'cut': 1, 'common': 1})
    out += group('cone_cut_rod', 'cone', 'cylinder',
                 [FRUSTUM, FRUSTUM_CUT, prism([disc(0.1, 2.5, 0.4)], at('XY', (0, 0, 0)), -1.0, 6.0)],
                 'cut', {'fuse': 1, 'cut': 1, 'common': 1})
    out += group('cone_cut_ball', 'cone', 'sphere',
                 [FRUSTUM, FRUSTUM_CUT, sphere(0.6, at('XY', (0.2, -1.2, 3.6)))], 'cut',
                 {'fuse': 1, 'cut': 1, 'common': 1})
    out += group('torus_cut_pipe', 'spiric', 'cylinder',
                 [RING, big_box(2.5, -8.0, -4.0, 8.0, 8.0, 4.0),
                  prism([disc(1.6, 1.0, 0.4)], at('SIDE', (0, 0, 0)), -6.0, 6.0)], 'cut',
                 {'fuse': 1, 'cut': 1, 'common': 2})
    out += group('torus_cut_ball', 'spiric', 'sphere',
                 [RING, big_box(2.5, -8.0, -4.0, 8.0, 8.0, 4.0), sphere(0.5, at('XY', (2.6, -1.6, -0.9)))],
                 'cut', {'fuse': 1, 'cut': 1, 'common': 1})
    out += group('peg_touch', 'rise', 'plane', [BALL5, PEG, big_box(-8.0, -8.0, 4.0, 8.0, 8.0, 12.0)], 'fuse',
                 {'fuse': 1, 'cut': 1, 'common': 1}, reason=TOUCH, touch=(0, 3, 4))
    # Its common is left out: the ball's and the given result's pieces touch
    # at the point, which rays count unsteadily.
    out += group('peg_kiss', 'rise', 'sphere', [BALL5, PEG, sphere(0.8125, at('XY', (0, 3.75, 4.3125)))], 'fuse',
                 {'fuse': 1, 'cut': 1}, reason=KISS, touch=(0, 3, 4))
    return out


def validate(listed):
    names = [c.name for c in listed]
    assert len(names) == len(set(names)), 'duplicate case names'
    for c in listed:
        for frame in c.frames:
            assert frame_name(frame) in FRAMES, c.name
        assert len(c.specs) == 3, c.name
        assert (c.kind == 'degenerate') == (c.touch is not None), c.name


# ------------------------------------------------------------------ the chains' work

def evaluate(job):
    """One chain: S9e.3a's rows and checks (`generate_given_curved_boolean_
    fixtures.evaluate`) and the triple points."""
    name, specs, op1, ops, swapped, degenerate, touch, mc_n = job
    # A declared tangency's sections touch on the faces through it: roots of
    # their surfaces within rounding along a curve are one there, and the
    # quadrature takes 1e-22 (the touching pieces' endpoint singularities
    # converge slowly; the kernel refuses these cases, their rows are
    # rounded to binary64).
    ref.MERGE = degenerate
    ref.QUAD = mp.mpf(10)**(-22 if degenerate else -32)
    out = gc.evaluate((name, specs, op1, [], ops, swapped, degenerate, mc_n))
    cases_ = [gc.solid(s, gc.OPERATIONS[k]) for k, s in enumerate(specs)]
    chain = ref.Chain(cases_)
    tracer = gm.Tracer(chain, (op1, 0, 1), 2)
    points = tracer.points()
    declared = tracer.declared(touch) if touch is not None else None
    return out+(points, declared)


def cached(job):
    """`evaluate`, its result kept in `--cache` (a development aid)."""
    cache = os.environ.get('GIVEN_MET_CACHE')
    if not cache:
        return evaluate(job)
    import hashlib
    import inspect
    import pickle
    key = hashlib.sha256((repr(job)+inspect.getsource(evaluate)+inspect.getsource(gc.evaluate)
                          + Path(ref.__file__).read_text()+Path(gm.__file__).read_text()).encode()).hexdigest()[:24]
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


def klass_points(points, klass):
    """The triple points of the chain's class: on a meeting of two curved
    faces, or on a section met by a curved face (`given_met_reference.in_scope`)."""
    return [p for p in points if gm.in_scope(*p[2])]


def reference_checks(results, listed):
    worst, covered = {}, {}

    def note(key, value):
        worst[key] = max(worst.get(key, mp.mpf(0)), value)
        covered[key] = covered.get(key, 0)+1

    chains = {c.chain_name: c for c in listed}
    for name, rows, _, checks, margins, stats, degenerate, points, declared in results:
        tag = '' if declared is None else '_declared'
        for key, value in checks.items():
            note(key+tag, value)
        note('quadrature_estimate'+tag, stats['quadrature'])
        assert stats['missed'] == 0, f'{name}: quadrature nodes off their structure'
        if declared is not None:
            values, d, edge, face = declared
            assert edge and face, f'{name}: the declared tangency off the given edge or the partner\'s face'
            note('declared_values', values)
            note('declared_determinant', d)
            continue
        note('least_sine', 1/max(margins['sine'], mp.mpf(10)**-30))
        note('event_spacing', 1/max(margins['spacing'], mp.mpf(10)**-40))
        mine = klass_points(points, chains[name].klass)
        assert mine, f'{name}: no triple point of its class on the given edges'
        for X, sine, kinds, found in points:
            assert found >= 2, f'{name}: a triple point found once at {[mp.nstr(x, 8) for x in X]}'
            note('triple_sine', 1/sine)
        partner = chains[name].partner
        assert any(k[2] == partner for _, _, k, _ in mine), f'{name}: no triple point on the partner\'s {partner}'
    return worst, covered


def generate(results, listed):
    by_chain = {r[0]: r for r in results}
    blocks = []
    out = ['# case\trow (S9e.3b, chained_curved_boolean_reference.py and given_met_reference.py: expect KIND S9e.3b '
           'CLASS PARTNER, reason TEXT for a degenerate case, meet N SINE, then result N volume area cx cy cz or '
           'empty)']
    for case in listed:
        blocks.append(case.encode())
        r = by_chain[case.chain_name]
        rows, points = r[1], r[7]
        row = rows[case.last]
        n = 0 if row[0] == 'empty' else int(row[0].split()[1])
        want = {'solid': n > 0, 'empty': n == 0, 'degenerate': True}[case.kind]
        assert want, f'{case.name}: declared {case.kind}, the reference gives {row[0]}'
        out.append(f'{case.name}\texpect {case.kind} S9e.3b {case.klass} {case.partner}')
        if case.kind == 'degenerate':
            out.append(f'{case.name}\treason {case.reason}')
        else:
            mine = klass_points(points, case.klass)
            out.append(f'{case.name}\tmeet {len(mine)} {mp.nstr(min(p[1] for p in mine), 3)}')
        for x in row:
            out.append(f'{case.name}\t{x}')
    frames = ['# case\tframe (solid k) axis (n, x, y)\tstored unit vector (stored_axes, as hex bits)']
    for case in listed:
        for k, c in enumerate(case.solid_cases):
            _, x, y, n = stored_axes(c.frame)
            for key, v in (('n', n), ('x', x), ('y', y)):
                frames.append(f'{case.name}\tsolid{k} {key}\t'+' '.join(struct.pack('>d', q).hex() for q in v))
    return {'boolean-given-met-cases.txt': '\n'.join(blocks)+'\n',
            'boolean-given-met-expected.tsv': '\n'.join(out)+'\n',
            'boolean-given-met-frames.tsv': '\n'.join(frames)+'\n'}


def jobs(listed, mc_n):
    chains = {}
    for c in listed:
        e = chains.setdefault(c.chain_name, [c.specs, c.op1, {}, c.swapped, c.kind == 'degenerate', c.touch])
        e[2][c.last] = c.solids
    return [(name, specs, op1, ops, sw, degenerate, touch, mc_n)
            for name, (specs, op1, ops, sw, degenerate, touch) in chains.items()]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--workers', type=int, default=max(1, min(4, (os.cpu_count() or 2)-2)))
    parser.add_argument('--samples', type=int, default=200000, help='Monte-Carlo points per chain')
    parser.add_argument('--only', help='one chain, printed (no files)')
    parser.add_argument('--cache', help='a directory keeping each chain\'s result (development)')
    args = parser.parse_args()
    if args.cache:
        os.environ['GIVEN_MET_CACHE'] = args.cache
    listed = cases()
    validate(listed)
    todo = jobs(listed, args.samples)
    if args.only:
        todo = [j for j in todo if j[0] == args.only]
        for r in run(todo, cached, 1):
            print(r)
        return
    results = run(todo, cached, args.workers)
    worst, covered = reference_checks(results, listed)
    limits = {'families': 1e-30, 'solid_closed_forms': 1e-30, 'pair_identities': 1e-30, 'area_identity': 1e-30,
              'monte_carlo_sigma': 5, 'quadrature_estimate': 1e-30, 'least_sine': 20, 'event_spacing': 1e6,
              'triple_sine': 20, 'declared_values': 1e-30, 'declared_determinant': 1e-30}
    # The declared tangencies' checks at their quadrature's coarser tolerance.
    for key in ('families', 'solid_closed_forms', 'pair_identities', 'area_identity', 'quadrature_estimate'):
        limits[key+'_declared'] = 1e-18
    limits['monte_carlo_sigma_declared'] = 5
    for key, value in worst.items():
        assert value <= limits[key], (key, mp.nstr(value, 3))
    files = generate(results, listed)
    for name, contents in files.items():
        target = ROOT/'fixtures'/name
        if args.check:
            if target.read_text() != contents:
                raise SystemExit(f'{target} is stale')
        else:
            target.write_text(contents)
    kinds = {}
    for c in listed:
        kinds[c.kind] = kinds.get(c.kind, 0)+1
    print(len(listed), 'cases:', ', '.join(f'{sum(1 for c in listed if c.last == op)} {op}' for op in OPS),
          '-', ', '.join(f'{v} {k}' for k, v in sorted(kinds.items())),
          '-', ', '.join(f'{sum(1 for c in listed if c.klass == k)} {k}' for k in CLASSES),
          '-', ', '.join(f'{sum(1 for c in listed if c.partner == k)} {k}' for k in PARTNERS))
    print('reference checks (largest deviation, relative to the case size; least_sine, event_spacing and '
          'triple_sine as reciprocals):', ', '.join(f'{k} {mp.nstr(v, 3)}' for k, v in sorted(worst.items())))
    print('cases per check:', ', '.join(f'{k} {v}' for k, v in sorted(covered.items())))


if __name__ == '__main__':
    main()
