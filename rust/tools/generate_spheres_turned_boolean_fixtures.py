#!/usr/bin/env python3
"""Fixtures for S9d.2c of REVIEW_NOTES.md: Booleans of a sphere against a
prism with circles where a frame is turned, before any of its kernel code:
a cap's circles (its rim and the split's great circle, of unequal axes on a
turned frame's stored axis) against a cylinder, and a sphere meeting a
cylinder in a turned frame in a loop.

`boolean-spheres-turned-cases.txt` lists each case in the Boolean protocol
(`identity_reference.encode_boolean_case`); `boolean-spheres-turned-
expected.tsv` gives per case `expect KIND S9d.2c` (`solid`, `empty` or
`degenerate`, then `reason TEXT`) and `result N volume area cx cy cz` or
`empty` from `spheres_boolean_reference.py` (a cap's end planes as the
kernel's `Ball` reads them, `AxisSphere`: through `o + h n` normal to the
stored axis); `boolean-spheres-turned-frames.tsv` the stored axes' bits.
Frames are the curved generator's (`FRAMES`). Caps are hemispheres only
(latitude 0, a height of exactly 0 whatever the platform's `sin`), centred
at the world's origin.

Pairs. Turned caps against cylinders in an exact frame: a hemisphere in
`TILT` against a coaxial pipe through its disc and dome, the rim and the
split crossing the pipe's wall between its rings (`dome_tilt_pipe`), the
lower hemisphere in `LEAN` against S9d.2's bite (`dome_lean_bite`, all three
operations: the equator crossing the loop), a hemisphere in `TILTX` against
S9d.2's rod (`dome_tiltx_rings`: both rings cut by the equator). Loops in a
turned frame: S9d.2's bite in `LEAN` (`bite_lean`, all three operations)
and its fuse and common in `TILT` (`bite_tilt`; the cut is S9d.2's
fixture), a long loop round most of a thick cylinder in `TILTX` about an
offset sphere (`graze_tiltx`), a thin cylinder turned about the axis
(`bite_r125`, `R125`). Both at once: a hemisphere in `TILT` against the
bite in `LEAN` (`dome_bite_turned`). And `degenerate`: a hemisphere in
`TILT` whose rim is tangent (within rounding) to a cylinder crossing the
sphere (`rim_tangent`, its axis through the rim's tangent point's outward
normal of the projected rim, the centre rounded once).

Before writing, `reference_checks` compares the reference with independent
results (limits relative to the case's size):

* closed forms where available: a whole sphere against a cylinder in a
  turned frame by S9d.2's lens of two discs integrated along the axis
  (`generate_spheres_boolean_fixtures.prism_forms`); a hemisphere against a
  coaxial pipe covering it (`dome_tilt_pipe`) by its sections along the
  pipe's axis, each a disc cut by the hemisphere's plane (a circular
  segment's area and moment in closed form, one quadrature between the
  kinks), its faces by central symmetry (half the ball's sphere face and
  wall inside) and the disc inside the pipe as a circle and a concentric
  ellipse (closed form); within 1e-15 (the frames' ideal axes against their
  stored ones);
* halves: a cap and its complement in the same frame (their end planes one
  plane) against the prism sum to the whole sphere sliced along the
  prism's own axis, volumes, moments, the sphere face's area inside, each
  prism face's area inside; within 1e-30;
* S9d.2's checks (`generate_spheres_boolean_fixtures.evaluate`): inclusion
  and exclusion, both inputs against their closed forms, the area
  identity, every face's classes summing to its area, a second slicing
  direction for whole spheres, Monte-Carlo estimates within 5 standard
  errors, and a scan for near coincidences and cap circles near tangency
  (none but in the declared degenerate pair).

Pairs run in worker processes (`--workers`).
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
import math
import os
from pathlib import Path
import struct

import mpmath as mp

from curve_surface_reference import stored_axes
import generate_spheres_boolean_fixtures as g2
import spheres_boolean_reference as ref
from generate_curved_boolean_fixtures import FRAMES
from generate_sphere_boolean_fixtures import frame_name, ideal_axes

ROOT = Path(__file__).resolve().parents[1]
PREFIX = 'boolean-spheres-turned'
STEP = 'S9d.2c'
HP = math.pi/2
EXACT = {'XY', 'SIDE', 'DOWN', 'TURN'}

at, circle, sphere, prism, group = g2.at, g2.circle, g2.sphere, g2.prism, g2.group

RIM_TANGENT = 'a cap\'s rim tangent to a cylinder'


def tangent_rod(R, a):
    """A cylinder of radius `a` along `z` tangent (ideally) to the rim of a
    hemisphere of radius `R` in `TILT` at its point of angle `(3/5, 4/5)`:
    the rim projects onto `(R cos t, 4/5 R sin t)`, whose outward normal
    there is `(cos t / R, sin t / (4/5 R))`; the centre the tangent point
    moved `a` along it, rounded once."""
    R, a = mp.mpf(R), mp.mpf(a)
    c, s = mp.mpf(3)/5, mp.mpf(4)/5
    p = (R*c, R*s*4/5)
    nrm = (c/R, s/(R*4/5))
    k = a/mp.sqrt(nrm[0]**2+nrm[1]**2)
    return float(p[0]+k*nrm[0]), float(p[1]+k*nrm[1]), float(a)


def cases():
    """S9d.2c's pairs (see the module's docstring)."""
    O = (0, 0, 0)
    out = []
    out += group('dome_tilt_pipe', sphere(2, at('TILT', O), 0.0, HP),
                 prism([circle(0, 0, 1.75)], at('XY', O), -1.5, 3), {'cut': 'solid', 'common': 'solid'})
    bite = prism([circle(2, 0.5, 1)], at('XY', O), -3, 3)
    out += group('dome_lean_bite', sphere(2, at('LEAN', O), -HP, 0.0), bite,
                 {'fuse': 'solid', 'cut': 'solid', 'common': 'solid'})
    out += group('dome_tiltx_rings', sphere(2, at('TILTX', O), 0.0, HP),
                 prism([circle(0.5, 0.25, 0.75)], at('XY', O), -3, 3), {'cut': 'solid', 'common': 'solid'})
    out += group('bite_lean', sphere(2, at('XY', O)), prism([circle(2, 0.5, 1)], at('LEAN', O), -3, 3),
                 {'fuse': 'solid', 'cut': 'solid', 'common': 'solid'})
    out += group('bite_tilt', sphere(2, at('XY', O)), prism([circle(2, 0.5, 1)], at('TILT', O), -3, 3),
                 {'fuse': 'solid', 'common': 'solid'})
    out += group('graze_tiltx', sphere(2, at('XY', (0.25, -0.25, 0.5))),
                 prism([circle(1.5, 0, 1.75)], at('TILTX', O), -3.5, 3.5), {'cut': 'solid', 'common': 'solid'})
    out += group('bite_r125', sphere(2, at('XY', O)), prism([circle(1.75, -0.5, 0.75)], at('R125', O), -2.5, 2.5),
                 {'cut': 'solid'})
    out += group('dome_bite_turned', sphere(2, at('TILT', O), 0.0, HP),
                 prism([circle(2, 0.5, 1)], at('LEAN', O), -3, 3), {'cut': 'solid', 'common': 'solid'})
    out += group('rim_tangent', sphere(2, at('TILT', O), 0.0, HP),
                 prism([circle(*tangent_rod(2, 0.5))], at('XY', O), -3, 3),
                 {'common': ('degenerate', RIM_TANGENT)})
    return out


def is_cap(spec):
    return spec[0] == 'sphere' and (spec[3] != -HP or spec[4] != HP)


# ------------------------------------------------------------------ closed forms

Z = mp.mpf(0)


def segment(rho, y0):
    """A disc of radius `rho` about the origin above the chord `y = y0`: its
    area and first moment in `y`."""
    if y0 <= -rho:
        return mp.pi*rho*rho, Z
    if y0 >= rho:
        return Z, Z
    h = mp.sqrt(rho*rho-y0*y0)
    return rho*rho*mp.acos(y0/rho)-y0*h, 2*h**3/3


def disc_ellipse(R, a, b):
    """The area common to a circle of radius `R` and a concentric ellipse of
    semi-axes `a < R < b` (the first along the circle's `x`)."""
    tan2 = (1/a**2-1/R**2)/(1/R**2-1/b**2)
    phi = mp.atan(mp.sqrt(tan2))
    return 4*(a*b/2*mp.atan(a/b*mp.tan(phi))+R*R/2*(mp.pi/2-phi))


def dome_pipe_forms(first):
    """A hemisphere (radius `R` about the origin, `n . X >= 0`, `n` the
    frame's ideal axis with no `x` part) against a pipe along `z` of radius
    `a` about the origin whose ends miss it: sections along `z` are discs of
    radius `min(a, sqrt(R^2 - z^2))` above the chord `n_y y + n_z z = 0`;
    faces by central symmetry (the ball less the pipe is symmetric through
    the centre, the hemisphere's plane holds it), the flat face inside the
    pipe a circle and an ellipse of semi-axes `a` and `a / n_z`."""
    s, p = first.specs
    R = mp.mpf(s[1])
    _, _, n = ideal_axes(frame_name(s[2]))
    assert n[0] == 0 and s[3] == 0.0 and s[4] == HP
    (cx, cy, a), = [b.circle for b in p[1]]
    assert cx == 0 and cy == 0
    a = mp.mpf(a)
    z0, z1 = mp.mpf(p[3]), mp.mpf(p[4])
    assert z0 < -R*n[1] and z1 > R, 'the pipe covers the hemisphere'
    h0 = mp.sqrt(R*R-a*a)
    rho = lambda z: min(a, mp.sqrt(max(R*R-z*z, Z)))
    y0 = lambda z: -n[2]*z/n[1]
    kinks = sorted({-R, R, -h0, h0}|{k for k in (a*n[1]/n[2], -a*n[1]/n[2], R*n[1], -R*n[1]) if -R < k < R})
    V = mp.quad(lambda z: segment(rho(z), y0(z))[0], kinks)
    My = mp.quad(lambda z: segment(rho(z), y0(z))[1], kinks)
    Mz = mp.quad(lambda z: z*segment(rho(z), y0(z))[0], kinks)
    sphere_in = 2*mp.pi*R*(R-h0)
    disc_in = disc_ellipse(R, a, a/n[2])
    wall_in = 2*mp.pi*a*h0
    VH = 2*mp.pi*R**3/3
    H = (VH, tuple(VH*3*R/8*n[i] for i in range(3)), 3*mp.pi*R*R)
    VP = mp.pi*a*a*(z1-z0)
    P = (VP, (Z, Z, VP*(z0+z1)/2), 2*mp.pi*a*(z1-z0)+2*mp.pi*a*a)
    return g2.assemble(H, P, (V, (Z, My, Mz)), sphere_in+disc_in, wall_in)


def closed_forms(first):
    """{op: (volume, area, centre) or None}, or None where no closed form is
    taken."""
    obj, tool = first.specs
    if first.pair_name == 'dome_tilt_pipe':
        return dome_pipe_forms(first)
    if is_cap(obj) or is_cap(tool):
        return None
    return g2.prism_forms(first)


# ------------------------------------------------------------------ halves

def complement(spec):
    """The other half of a hemisphere, in the same frame."""
    _, r, frame, low, high = spec
    assert (low, high) in ((0.0, HP), (-HP, 0.0))
    return ('sphere', r, frame, -HP, 0.0) if low == 0.0 else ('sphere', r, frame, 0.0, HP)


def halves(first):
    """A cap and its complement against the prism against the whole sphere
    sliced along the prism's axis: the largest deviation of the common's
    volume and moments, the sphere face inside, each prism face's class
    inside (relative to the case's size)."""
    obj, tool = first.specs
    if not (is_cap(obj) or is_cap(tool)):
        return None
    swap = is_cap(tool)
    s, p = (tool, obj) if swap else (obj, tool)
    whole = ('sphere', s[1], s[2], -HP, HP)
    pairs = []
    for spec in (s, complement(s), whole):
        a, b = g2.make(spec, 91), g2.make(p, 92)
        pairs.append(ref.Pair(a, b))
    size = pairs[2].size
    res = [x.sliced() for x in pairs]
    dev = Z
    common = [r['common'] for r in res]
    dev = max(dev, abs(common[0][0]+common[1][0]-common[2][0])/size**3)
    for i in range(3):
        dev = max(dev, abs(common[0][1][i]+common[1][1][i]-common[2][1][i])/size**4)
    sph = [r[('sphere', 'A')]['in'] for r in res]
    dev = max(dev, abs(sph[0]+sph[1]-sph[2])/size**2)
    faces = [{tag: cls for w, tag, cls, _ in x.faces() if w == 'B'} for x in pairs]
    for tag in faces[2]:
        dev = max(dev, abs(faces[0][tag]['in']+faces[1][tag]['in']-faces[2][tag]['in'])/size**2)
    return dev


# ------------------------------------------------------------------ the pairs' work

def evaluate(job):
    """S9d.2's checks on a pair, and the halves."""
    out = g2.evaluate(job)
    name, first = job[0], job[1]
    h = halves(first)
    if h is not None:
        out[3]['halves'] = h
    return out


def run(jobs, fn, workers):
    if workers <= 1:
        return [fn(j) for j in jobs]
    with ProcessPoolExecutor(workers) as pool:
        return list(pool.map(fn, jobs))


def closed_job(first):
    return first.pair_name, closed_forms(first)


LIMITS = {'closed_forms_turned_frames': 1e-15, 'inputs_sliced': 1e-30, 'inclusion_exclusion': 1e-30,
          'face_classes': 1e-30, 'shared_both_sides': 1e-30, 'area_identity': 1e-30, 'input_areas': 1e-30,
          'second_direction': 1e-30, 'halves': 1e-30, 'monte_carlo_sigma': 5, 'volume_quadrature_estimate': 1e-30}


def reference_checks(results, forms):
    """Largest deviations of each check, relative to the case's size."""
    worst, covered = {}, {}

    def note(key, value):
        worst[key] = max(worst.get(key, Z), value)
        covered[key] = covered.get(key, 0)+1
    by_pair = {r[0]: r for r in results}
    for name, form in forms.items():
        if form is None:
            continue
        _, _, res, _, _, _ = by_pair[name]
        for op, want in form.items():
            n, vol, area, centre = res[op]
            if want is None:
                assert n == 0, (name, op, 'the closed form is empty')
                continue
            V, A, C = want
            size = max(abs(V)**(mp.mpf(1)/3), 1)
            dev = max(abs(vol-V)/abs(V), abs(area-A)/abs(A), max(abs(centre[i]-C[i]) for i in range(3))/size)
            note('closed_forms_turned_frames', dev)
    for name, rows, res, checks, near, stats in results:
        for key, value in checks.items():
            note(key, value)
        note('volume_quadrature_estimate', stats['volume_quadrature'])
    return worst, covered


def generate(results):
    by_pair = {r[0]: r for r in results}
    blocks = []
    out = [f'# case\trow ({STEP}, spheres_boolean_reference.py: expect KIND {STEP}, reason TEXT for a degenerate '
           'case, then result N volume area cx cy cz or empty)']
    for case in cases():
        blocks.append(case.encode())
        _, rows, _, _, _, _ = by_pair[case.pair_name]
        row = rows[case.operation]
        n = 0 if row[0] == 'empty' else int(row[0].split()[1])
        want = {'solid': n > 0, 'empty': n == 0, 'degenerate': True}[case.kind]
        assert want, f'{case.name}: declared {case.kind}, the reference gives {row[0]}'
        out.append(f'{case.name}\texpect {case.kind} {STEP}')
        if case.kind == 'degenerate':
            out.append(f'{case.name}\treason {case.reason}')
        for r in row:
            out.append(f'{case.name}\t{r}')
    frames = ['# case\tframe (obj, tool) axis (n, x, y)\tstored unit vector (stored_axes, as hex bits)']
    for case in cases():
        for label, c in (('obj', case.obj), ('tool', case.tool)):
            _, x, y, n = stored_axes(c.frame)
            for key, v in (('n', n), ('x', x), ('y', y)):
                frames.append(f'{case.name}\t{label} {key}\t'+' '.join(struct.pack('>d', q).hex() for q in v))
    return {f'{PREFIX}-cases.txt': '\n'.join(blocks)+'\n',
            f'{PREFIX}-expected.tsv': '\n'.join(out)+'\n',
            f'{PREFIX}-frames.tsv': '\n'.join(frames)+'\n'}


def jobs(mc_n):
    pairs = {}
    for c in cases():
        pairs.setdefault(c.pair_name, [c, []])[1].append(c.operation)
    return [(name, first, ops, mc_n) for name, (first, ops) in pairs.items()]


def validate(listed):
    names = [c.name for c in listed]
    assert len(names) == len(set(names)), 'duplicate case names'
    for c in listed:
        specs = c.specs
        assert any(s[0] == 'sphere' for s in specs) and any(s[0] == 'prism' for s in specs), \
            f'{c.name}: a sphere against a prism'
        turned = [frame_name(s[2]) not in EXACT for s in specs if s[0] == 'prism' or is_cap(s)]
        assert any(turned), f'{c.name}: a turned cap or a turned cylinder'
        for s in specs:
            if is_cap(s):
                assert s[2][:3] == (0.0, 0.0, 0.0) and 0.0 in (s[3], s[4]), f'{c.name}: a hemisphere at the origin'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--workers', type=int, default=max(1, min(6, (os.cpu_count() or 2)-2)))
    parser.add_argument('--samples', type=int, default=200000, help='Monte-Carlo points per pair')
    parser.add_argument('--only', help='evaluate only this pair (no files written)')
    args = parser.parse_args()
    listed = cases()
    validate(listed)
    todo = jobs(args.samples)
    if args.only:
        todo = [j for j in todo if j[0] == args.only]
    results = run(todo, evaluate, args.workers)
    degenerate = {c.pair_name for c in listed if c.kind == 'degenerate'}
    near = [f'{name}: near coincidences {n}' for name, _, _, _, n, _ in results if n and name not in degenerate]
    firsts = {}
    for c in listed:
        firsts.setdefault(c.pair_name, c)
    done = {r[0] for r in results}
    forms = dict(run([f for k, f in firsts.items() if k in done], closed_job, args.workers))
    worst, covered = reference_checks(results, forms)
    print('reference checks (largest deviation, relative to the case size):',
          ', '.join(f'{k} {mp.nstr(v, 3)}' for k, v in sorted(worst.items())))
    print('cases per check:', ', '.join(f'{k} {v}' for k, v in sorted(covered.items())))
    for name, rows, _, _, n, stats in results:
        print(name, 'breaks', stats['breaks'], 'clearance',
              None if stats['clearance'] is None else mp.nstr(stats['clearance'], 3),
              '; '.join(f'{op} {r[0]}' for op, r in rows.items()))
        if n:
            print('  near coincidences:', '; '.join(n[:4]))
    if near:
        raise SystemExit('\n'.join(near))
    for key, value in worst.items():
        assert value <= LIMITS[key], (key, mp.nstr(value, 3))
    if args.only:
        return
    files = generate(results)
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
    print(len(listed), 'cases:', ', '.join(f'{sum(1 for c in listed if c.operation == op)} {op}'
                                           for op in ('fuse', 'cut', 'common')),
          '-', ', '.join(f'{v} {k}' for k, v in sorted(kinds.items())))


if __name__ == '__main__':
    main()
