#!/usr/bin/env python3
"""Fixtures for S9d.3c of REVIEW_NOTES.md: Booleans of a cone or frustum
(`Solid::cone_with`) against a cylinder or another cone meeting it in loops
(the discriminant of both inputs' rulings against the other's quadric
changing sign over their angles), in exact and turned frames, of a turned
cone against a sphere meeting it in a loop, and of a turned cap (a
hemisphere, its circles of unequal axes on its stored axis) against a cone,
before any of its kernel code.

`boolean-cones-loops-cases.txt` lists each case in the Boolean protocol
(`identity_reference.encode_boolean_case`, as S9d.3b's);
`boolean-cones-loops-expected.tsv` gives per case `expect KIND S9d.3c`
(`solid`, `empty` or `degenerate`, then `reason TEXT`) and `result N volume
area cx cy cz` or `empty` from `cones_boolean_reference.py` (a cap's end
planes as the kernel's `Ball` reads them, S9d.2c's `AxisSphere`);
`boolean-cones-loops-frames.tsv` the stored axes' bits. Frames are the
curved generator's (`FRAMES`). Caps are hemispheres centred at the world's
origin (a height of exactly 0 whatever the platform's `sin`).

Pairs, with `F` the frustum of radii 2 and 1 and height 2 and `A` the cone
of radius 3/2 and height 3 on `XY` at the origin unless said. A cone and a
cylinder in loops: a rod across the axis grazing `F`'s wall
(`rod_graze`, all three operations; the rod less `F`, `rod_bitten`), a rod
across `A`'s tip region beside its apex (`tip_graze`), a rod in `LEAN`
grazing `F`'s wall (`rod_lean`), a vertical rod against `F` in `TILT`
(`tilt_rod`). Two cones in loops: a thin frustum across the axis grazing
`F`'s wall (`cones_graze`, all three), `F` in `TILT` against a frustum in
`LEAN` (`cones_turned`), and a cone in `LEAN` beside `F` whose direction
cone crosses `F`'s, its curve running through infinity (`cones_asymptotic`).
A turned cone against a sphere in a loop: `F` in `TILT` (`ball_tilt`, all
three; S9d.3b's pair), `A` in `LEAN` (`ball_lean`), `F` turned about its
axis (`ball_r125`). A turned hemisphere against a cone: a coaxial cone
through the hemisphere in `TILT` (`dome_tilt_cone`: its rim and disc cut
by the wall, the dome in a circle), the lower hemisphere in `LEAN` against
a frustum off the axis (`dome_lean_frustum`: the rim across the loop), a
hemisphere in `TILT` against a frustum in `LEAN` (`dome_cone_turned`:
both). And `degenerate`: a cylinder exactly along a cone's ruling
(`rod_ruling_exact`: the cone's slope `fl(0.6) / fl(0.8)` of `LEAN`'s
stored axis, so the cylinder's `A` is zero), a hemisphere in `TILT` whose
rim is tangent to a cone within rounding (`rim_tangent_cone`).

Before writing, `reference_checks` compares the reference with independent
results (limits relative to the case's size):

* closed forms by one quadrature where available: a rod across a cone's
  axis by S9d.3b's strip (`generate_cones_boolean_fixtures.strip_form`, the
  disc of the cone's section against the rod's strip along the axis); a
  whole sphere against a turned cone by S9d.3b's lens along the cone's
  ideal axis (`axial_form`); a hemisphere against a coaxial cone by its
  sections along the axis, discs of radius the smaller of the cone's and
  the sphere's cut by the hemisphere's plane (circular segments), its
  sphere face and the cone's wall by their circles' arcs above the chord,
  the flat disc inside the cone by its chords in its own plane
  (`dome_cone_forms`); within 1e-30 in exact frames and 1e-15 in turned
  ones (their stored axes not exactly orthonormal);
* splits: an input cut in two along its axis (a prism at a height between
  its caps in its own frame, a cone in an exact frame at a dyadic height,
  both parts exactly the input's model) against the other, the parts'
  commons summing to the whole's (volume, moments), the other's faces
  inside summing, the cone's wall inside summing; within 1e-30;
* halves: a hemisphere and its complement against the cone summing to the
  whole sphere's (S9d.2c's check, the cones' reference); within 1e-30;
* S9d.3b's checks (`generate_cones_boolean_fixtures.evaluate`): inclusion
  and exclusion, both inputs against their closed forms, the area
  identity, every face's classes summing to its area, the cone's wall two
  ways where the chart is its own, a second slicing direction (pairs
  without a cap or a declared degeneracy), Monte-Carlo estimates within 5
  standard errors, and a scan for near coincidences and edges or vertices
  within 1e-3 of tangency with or incidence on the other's surfaces (none
  but in the declared degenerate pairs).

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
import cones_boolean_reference as ref
import generate_cones_boolean_fixtures as g3
from generate_sphere_boolean_fixtures import frame_name, ideal_axes

ROOT = Path(__file__).resolve().parents[1]
PREFIX = 'boolean-cones-loops'
STEP = 'S9d.3c'
HP = math.pi/2
EXACT = g3.EXACT
Z = mp.mpf(0)

at, cone, sphere, prism, circle, group = g3.at, g3.cone, g3.sphere, g3.prism, g3.circle, g3.group
ALL = {'fuse': 'solid', 'cut': 'solid', 'common': 'solid'}
BOTH = {'cut': 'solid', 'common': 'solid'}

RULING = 'a cylinder along a cone\'s ruling exactly'
RIM_TANGENT = 'a cap\'s rim tangent to a cone'


def tangent_cone(R, k, rho, z0, h):
    """A cone along `z` (slope `k`, height `h` from `z0`) tangent (ideally)
    to the rim of a hemisphere of radius `R` in `TILT` at its point of angle
    `(3/5, 4/5)`, its radius there `rho`: the cone's normal `(h_xy / rho,
    -k)` at the point normal to the rim's tangent `T`, so the horizontal
    offset `h_xy` from the axis makes the angle `acos(k T_z / |T_xy|)` with
    `T_xy`, on the side facing out of the rim. Its base centre and radius
    rounded once."""
    R, k, rho, z0 = mp.mpf(R), mp.mpf(k), mp.mpf(rho), mp.mpf(z0)
    c, s = mp.mpf(3)/5, mp.mpf(4)/5
    x, y = (1, 0, 0), (0, mp.mpf(4)/5, -mp.mpf(3)/5)
    p = tuple(R*(c*x[i]+s*y[i]) for i in range(3))
    T = tuple(R*(-s*x[i]+c*y[i]) for i in range(3))
    txy = mp.sqrt(T[0]**2+T[1]**2)
    phi = -mp.acos(k*T[2]/txy)
    u = (T[0]/txy, T[1]/txy)
    hx = (u[0]*mp.cos(phi)-u[1]*mp.sin(phi), u[0]*mp.sin(phi)+u[1]*mp.cos(phi))
    axis = (p[0]-rho*hx[0], p[1]-rho*hx[1])
    r0 = rho-k*(p[2]-z0)
    return (float(axis[0]), float(axis[1]), float(z0)), float(r0), float(r0+k*h), float(h)


def cases():
    """S9d.3c's pairs (see the module's docstring)."""
    O = (0, 0, 0)
    F_ = cone(2, 1, 2, at('XY', O))
    A_ = cone(1.5, 0, 3, at('XY', O))
    out = []
    rod = prism([circle(1.5, 1, 0.375)], at('SIDE', O), -3, 3)
    out += group('rod_graze', F_, rod, ALL, ('strip',))
    out += group('rod_bitten', rod, F_, {'cut': 'solid'}, ('strip',))
    out += group('tip_graze', A_, prism([circle(0.5, 2.25, 0.375)], at('SIDE', O), -3, 3), BOTH, ('strip',))
    out += group('rod_lean', F_, prism([circle(1.625, 0, 0.375)], at('LEAN', (0, 0, 1)), -3, 3), BOTH, None)
    out += group('tilt_rod', cone(2, 1, 2, at('TILT', O)), prism([circle(1.25, 1.25, 0.5)], at('XY', O), -1, 3),
                 BOTH, None)
    out += group('cones_graze', F_, cone(0.5, 0.25, 4.5, at('SIDE', (-2.25, 1.5, 1))), ALL, None)
    out += group('cones_turned', cone(2, 1, 2, at('TILT', O)), cone(1, 0.25, 3, at('LEAN', (1, 1.5, 0))), BOTH,
                 None)
    out += group('cones_asymptotic', F_, cone(1.5, 0, 3, at('LEAN', (0, 2.5, 0))), BOTH, None)
    out += group('ball_tilt', cone(2, 1, 2, at('TILT', O)), sphere(1, at('XY', (1.5, -0.25, 1.25))), ALL,
                 ('axial',))
    out += group('ball_lean', cone(1.5, 0, 3, at('LEAN', O)), sphere(0.5, at('XY', (0.9, 0.75, 1.2))), BOTH,
                 ('axial',))
    out += group('ball_r125', cone(2, 1, 2, at('R125', O)), sphere(1, at('XY', (1.25, 0.75, 1))), {'cut': 'solid'},
                 ('axial',))
    out += group('dome_tilt_cone', sphere(2, at('TILT', O), 0.0, HP), cone(1.25, 0, 4.5, at('XY', (0, 0, -2.25))),
                 BOTH, ('dome',))
    out += group('dome_lean_frustum', sphere(2, at('LEAN', O), -HP, 0.0), cone(1.25, 0.75, 3, at('XY', (1.5, 0.5, -1))),
                 BOTH, None)
    out += group('dome_cone_turned', sphere(2, at('TILT', O), 0.0, HP),
                 cone(1.25, 0.5, 3, at('LEAN', (0.5, 0.25, -1.5))), BOTH, None)
    out += group('rod_ruling_exact', cone(2.4, 0, 3.2, at('XY', O)),
                 prism([circle(0, 0, 0.5)], at('LEAN', (-1.25, 0.5, 1.5)), -3, 3),
                 {'common': ('degenerate', RULING)}, None)
    o, r0, r1, h = tangent_cone(2, -0.25, 0.75, -2.5, 4)
    out += group('rim_tangent_cone', sphere(2, at('TILT', O), 0.0, HP), cone(r0, r1, h, at('XY', o)),
                 {'common': ('degenerate', RIM_TANGENT)}, None)
    return out


# Splits: the input cut in two along its axis ('obj' or 'tool') at a height.
SPLITS = {'rod_graze': ('tool', 0.5), 'tip_graze': ('tool', -0.25), 'rod_lean': ('tool', 0.5),
          'tilt_rod': ('tool', 1.25), 'cones_graze': ('obj', 1), 'cones_asymptotic': ('obj', 1)}


def is_cap(spec):
    return spec[0] == 'sphere' and (spec[3] != -HP or spec[4] != HP)


# ------------------------------------------------------------------ closed forms

def arc(r, y0):
    """The angle of a circle of radius `r` about the origin above the chord
    `y = y0`."""
    return g3.segment(r, y0)[2]


def dome_cone_forms(first):
    """A hemisphere (radius `R` about the origin, `n . X >= 0`, `n` its
    frame's ideal axis with no `x` part) against a cone along `z` about the
    origin's vertical (radius `rho(z)` for `z` in `[z0, z1]`): sections along
    `z` are discs of radius `min(rho(z), s(z))`, `s(z) = sqrt(R^2 - z^2)`,
    above the chord `n_y y + n_z z = 0`; the sphere face inside where `s <=
    rho` (Archimedes: `R` per unit angle and height), the wall inside where
    `rho <= s` (`sqrt(1 + k^2) rho` per unit angle and height), the cone's
    end discs inside the ball above the chord, and the hemisphere's disc
    inside the cone by its chords along `x` in its own plane."""
    s_, c_ = first.specs
    R = mp.mpf(s_[1])
    _, _, n = ideal_axes(frame_name(s_[2]))
    assert n[0] == 0 and s_[3] == 0.0 and s_[4] == HP
    _, r0, r1, h, frame = c_
    assert frame_name(frame) == 'XY' and frame[0] == 0 and frame[1] == 0
    r0, r1, h, z0 = mp.mpf(r0), mp.mpf(r1), mp.mpf(h), mp.mpf(frame[2])
    z1 = z0+h
    k = (r1-r0)/h
    rho = lambda z: max(r0+k*(z-z0), Z)
    s = lambda z: mp.sqrt(max(R*R-z*z, Z))
    y0 = lambda z: -n[2]*z/n[1]
    lo, hi = max(z0, -R), min(z1, R)
    kinks = {lo, hi}
    kinks.update(g3.roots_on(lambda z: rho(z)-s(z), lo, hi))
    for sg in (1, -1):
        kinks.update(g3.roots_on(lambda z, sg=sg: y0(z)-sg*rho(z), lo, hi))
        kinks.update(g3.roots_on(lambda z, sg=sg: y0(z)-sg*s(z), lo, hi))
    pts = sorted(kinks)
    m = lambda z: min(rho(z), s(z))
    sec = lambda z: g3.segment(m(z), y0(z))
    V = mp.quad(lambda z: sec(z)[0], pts)
    My = mp.quad(lambda z: sec(z)[1], pts)
    Mz = mp.quad(lambda z: z*sec(z)[0], pts)
    sphere_in = mp.quad(lambda z: R*arc(s(z), y0(z)) if s(z) <= rho(z) else Z, pts)
    sk = mp.sqrt(1+k*k)
    wall_in = mp.quad(lambda z: sk*rho(z)*arc(rho(z), y0(z)) if rho(z) <= s(z) else Z, pts)
    discs = Z
    for z, r in ((z0, r0), (z1, r1)):
        if r > 0 and -R < z < R:
            discs += g3.segment(min(r, s(z)), y0(z))[0]
    # The hemisphere's disc: X = x e_x + t e, e = (0, n_z, -n_y), inside the
    # ball (x^2 + t^2 <= R^2) and the cone (x^2 + (n_z t)^2 <= rho(-n_y t)^2,
    # the height within the cone's).
    zt = lambda t: -n[1]*t
    yt = lambda t: n[2]*t

    def half(t):
        if not z0 <= zt(t) <= z1:
            return Z
        return mp.sqrt(max(min(R*R-t*t, rho(zt(t))**2-yt(t)**2), Z))
    tk = {-R, R}
    for z in (z0, z1):
        t = -z/n[1]
        if -R < t < R:
            tk.add(t)
    tk.update(g3.roots_on(lambda t: (R*R-t*t)-(rho(zt(t))**2-yt(t)**2), -R, R))
    tk.update(g3.roots_on(lambda t: rho(zt(t))**2-yt(t)**2, -R, R))
    disc_in = mp.quad(lambda t: 2*half(t), sorted(tk))
    VH = 2*mp.pi*R**3/3
    H = (VH, tuple(VH*3*R/8*n[i] for i in range(3)), 3*mp.pi*R*R)
    VK, MwK, AK = g3.cone_measures(r0, r1, h)
    K = (VK, (Z, Z, VK*z0+MwK), AK)
    return g3.assemble(H, K, (V, (Z, My, Mz)), sphere_in+disc_in, wall_in+discs)


def closed_forms(first):
    """{op: (volume, area, centre) or None}, or None where no closed form is
    taken."""
    if first.form is None:
        return None
    if first.form[0] == 'dome':
        return dome_cone_forms(first)
    return g3.closed_forms(first)


# ------------------------------------------------------------------ splits and halves

def split_specs(spec, at_):
    """An input cut in two along its axis at the height `at_`: a prism's
    heights in its own frame, a cone (in an exact frame) at the dyadic
    height, its upper part's base the stored origin moved along the stored
    exact axis."""
    if spec[0] == 'prism':
        _, bounds, frame, w0, w1 = spec
        assert w0 < at_ < w1
        return prism(bounds, frame, w0, at_), prism(bounds, frame, at_, w1)
    _, r0, r1, h, frame = spec
    assert frame_name(frame) in EXACT and 0 < at_ < h
    rm = r0+(r1-r0)*at_/h
    assert float(rm) == rm and rm > 0
    n = frame[3:6]
    base = tuple(frame[i]+at_*n[i] for i in range(3))
    assert all(float(c) == c for c in base)
    return cone(r0, rm, at_, frame), cone(rm, r1, h-at_, tuple(base)+frame[3:])


def face_in(pair, role):
    """{tag: area inside the other} of one input's faces."""
    return {tag: cls['in'] for w, tag, cls, _ in pair.faces() if w == role}


def splits(first):
    """The whole input against the other beside its two parts: the largest
    deviation of the commons' volume and moments, the other's faces inside
    and the cone's wall inside (relative to the case's size)."""
    which, at_ = SPLITS[first.pair_name]
    obj, tool = first.specs
    whole = obj if which == 'obj' else tool
    other = tool if which == 'obj' else obj
    role, orole = ('A', 'B') if which == 'obj' else ('B', 'A')
    parts = split_specs(whole, at_)

    def pair_of(x):
        a, b = (x, other) if which == 'obj' else (other, x)
        return ref.Pair(g3.make(a, 91), g3.make(b, 92))
    pairs = [pair_of(whole)]+[pair_of(p) for p in parts]
    size = pairs[0].size
    res = [p.sliced()['common'] for p in pairs]
    dev = abs(res[1][0]+res[2][0]-res[0][0])/size**3
    for i in range(3):
        dev = max(dev, abs(res[1][1][i]+res[2][1][i]-res[0][1][i])/size**4)
    other_in = [face_in(p, orole) for p in pairs]
    for tag, v in other_in[0].items():
        dev = max(dev, abs(other_in[1][tag]+other_in[2][tag]-v)/size**2)
    if whole[0] == 'cone':
        walls = [face_in(p, role)[('wall',)] for p in pairs]
        dev = max(dev, abs(walls[1]+walls[2]-walls[0])/size**2)
    return dev


def complement(spec):
    """The other half of a hemisphere, in the same frame."""
    _, r, frame, low, high = spec
    assert (low, high) in ((0.0, HP), (-HP, 0.0))
    return ('sphere', r, frame, -HP, 0.0) if low == 0.0 else ('sphere', r, frame, 0.0, HP)


def halves(first):
    """A hemisphere and its complement against the cone against the whole
    sphere: the largest deviation of the commons' volume and moments, the
    sphere face inside and each of the cone's faces inside (relative to the
    case's size)."""
    obj, tool = first.specs
    swap = is_cap(tool)
    s, k = (tool, obj) if swap else (obj, tool)
    whole = ('sphere', s[1], s[2], -HP, HP)
    pairs = []
    for spec in (s, complement(s), whole):
        pairs.append(ref.Pair(g3.make(spec, 91), g3.make(k, 92)))
    size = pairs[2].size
    common = [p.sliced()['common'] for p in pairs]
    dev = abs(common[0][0]+common[1][0]-common[2][0])/size**3
    for i in range(3):
        dev = max(dev, abs(common[0][1][i]+common[1][1][i]-common[2][1][i])/size**4)
    sph = [face_in(p, 'A')[('sphere',)] for p in pairs]
    dev = max(dev, abs(sph[0]+sph[1]-sph[2])/size**2)
    faces = [face_in(p, 'B') for p in pairs]
    for tag in faces[2]:
        dev = max(dev, abs(faces[0][tag]+faces[1][tag]-faces[2][tag])/size**2)
    return dev


# ------------------------------------------------------------------ the pairs' work

def evaluate(job):
    """S9d.3b's checks on a pair, the splits and the halves."""
    name, first, ops, mc_n = job
    try:
        out = g3.evaluate(job)
        if name in SPLITS:
            out[3]['splits'] = splits(first)
        obj, tool = first.specs
        if is_cap(obj) or is_cap(tool):
            out[3]['halves'] = halves(first)
        return out
    except Exception:
        import traceback
        return ('error', name, traceback.format_exc())


def run(jobs, fn, workers):
    if workers <= 1:
        return [fn(j) for j in jobs]
    with ProcessPoolExecutor(workers) as pool:
        return list(pool.map(fn, jobs))


def closed_job(first):
    return first.pair_name, closed_forms(first)


PER_PAIR = {}
LIMITS = {'closed_forms_exact_frames': 1e-30, 'closed_forms_turned_frames': 1e-15,
          'inputs_sliced': 1e-30, 'inclusion_exclusion': 1e-30, 'face_classes': 1e-30,
          'shared_both_sides': 1e-30, 'area_identity': 1e-30, 'input_areas': 1e-30, 'wall_two_ways': 1e-30,
          'second_direction': 1e-30, 'splits': 1e-30, 'halves': 1e-30, 'monte_carlo_sigma': 5,
          'volume_quadrature_estimate': 1e-30}


def exact_pair(case):
    """Whether the closed forms hold exactly: no cone, prism or cap in a
    turned frame (a whole sphere's model does not depend on its frame's
    axes)."""
    return all(frame_name(s[4] if s[0] == 'cone' else s[2]) in EXACT
               for s in case.specs if s[0] != 'sphere' or is_cap(s))


def reference_checks(results, forms):
    """Largest deviations of each check, relative to the case's size."""
    worst, covered = {}, {}

    def note(key, value):
        worst[key] = max(worst.get(key, Z), value)
        covered[key] = covered.get(key, 0)+1
    by_pair = {r[0]: r for r in results}
    firsts = {}
    for c in cases():
        firsts.setdefault(c.pair_name, c)
    for name, form in forms.items():
        if form is None:
            continue
        _, _, res, _, _, _ = by_pair[name]
        kind = 'exact' if exact_pair(firsts[name]) else 'turned'
        for op, want in form.items():
            n, vol, area, centre = res[op]
            if want is None:
                assert n == 0, (name, op, 'the closed form is empty')
                continue
            assert n > 0, (name, op, 'the reference is empty')
            V, A, C = want
            size = max(abs(V)**(mp.mpf(1)/3), 1)
            dev = max(abs(vol-V)/abs(V), abs(area-A)/abs(A), max(abs(centre[i]-C[i]) for i in range(3))/size)
            note(f'closed_forms_{kind}_frames', dev)
            PER_PAIR[name] = max(PER_PAIR.get(name, Z), dev)
    for name, rows, res, checks, near, stats in results:
        for key, value in checks.items():
            note(key, value)
        note('volume_quadrature_estimate', stats['volume_quadrature'])
    return worst, covered


def generate(results):
    by_pair = {r[0]: r for r in results}
    blocks = []
    out = [f'# case\trow ({STEP}, cones_boolean_reference.py: expect KIND {STEP}, reason TEXT for a degenerate '
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
        assert c.obj.cone is not None or c.tool.cone is not None, f'{c.name}: a cone'
        for s in c.specs:
            if is_cap(s):
                assert s[2][:3] == (0.0, 0.0, 0.0) and 0.0 in (s[3], s[4]), f'{c.name}: a hemisphere at the origin'
        for f in c.frames:
            frame_name(f)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--workers', type=int, default=max(1, min(6, (os.cpu_count() or 2)-2)))
    parser.add_argument('--samples', type=int, default=200000, help='Monte-Carlo points per pair')
    parser.add_argument('--only', help='evaluate only these pairs, comma separated (no files written)')
    args = parser.parse_args()
    listed = cases()
    validate(listed)
    todo = jobs(args.samples)
    if args.only:
        wanted = set(args.only.split(','))
        todo = [j for j in todo if j[0] in wanted]
    results = run(todo, evaluate, args.workers)
    errors = [r for r in results if r[0] == 'error']
    for _, name, tb in errors:
        print('FAILED', name, tb[-1500:])
    if errors:
        raise SystemExit(f'{len(errors)} pairs failed')
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
    for name, rows, _, checks, n, stats in results:
        print(name, 'closed', mp.nstr(PER_PAIR[name], 3) if name in PER_PAIR else None,
              'splits', mp.nstr(checks['splits'], 3) if 'splits' in checks else None,
              'halves', mp.nstr(checks['halves'], 3) if 'halves' in checks else None,
              'chart', stats['chart'], 'breaks', stats['breaks'], 'sweep', stats['sweep_breaks'],
              'clearance', None if stats['clearance'] is None else mp.nstr(stats['clearance'], 3),
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
          '-', ', '.join(f'{v} {k}' for k, v in sorted(kinds.items())),
          '-', f'{sum(1 for c in listed if exact_pair(c))} in exact frames (whole spheres in any)')


if __name__ == '__main__':
    main()
