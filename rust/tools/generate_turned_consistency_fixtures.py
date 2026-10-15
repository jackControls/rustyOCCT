#!/usr/bin/env python3
"""Fixtures for REVIEW_NOTES.md's "S9e.4b.4d.2c.1a's consistency fixed": an
imported turned body of `OCC485a`'s form whose stored geometry lies off its
model within the resolution, given to Booleans with a rib through its step,
fillet, shaft and round and with boxes across each of its cylinders, its
fillet and its round, every result checked against the others and against
the model's exact measures.

The body is OCCT's own output, as S9e.4b.4d.2a's to 2c.1a's:
`boolean-turned-consistency-bodies.txt` holds one `write NAME
imported/consistency_NAME.brep` block for `occt_boolean_oracle.cpp` (a
`plane` row, the profile's `wire` of lines and arcs as
`generate_turned_fillets_boolean_fixtures.py` writes them, a `revolve` row
turning it a whole turn about the axis through the origin along `(EPS, 0,
1)` and a `tolerance` row setting every tolerance to 3e-3, `OCC485a`'s),
which `--write-bodies` runs to write
`rust/fixtures/imported/consistency_NAME.brep` (format version 1, no
triangulations). Nothing here reads the kernel.

The body `ribbed`, about an axis tilted by `EPS` = 2^-27 rad from the
world's `z` in the `(x, z)` plane (`OCC485a`'s is tilted by 1.08e-8: the
kernel's model takes the axis along the world's `z` through the bottom's
apex, every stored surface off it by up to 3.4e-7), its profile `(r, w)`
(the point at `(r + EPS w, 0, w - EPS r)`, exactly): a bottom cone whose
rim of radius 24 stands `2^-15` above its apex (a near-flat cone, taken as
the plane through its rim), a cylinder to 32, a chamfer in to radius 21 at
35, a step cone in to radius 10 rising `2^-17` (taken as the plane through
its wider rim at 35), a fillet of radius 1/2 about `(10, 35.5 + 2^-17)` to
the shaft of radius 19/2 (snapped tangent to the step's plane: its centre
`(10, 35.5)`, moved by `2^-17`), the shaft to 41, a round about `(13/2, 41 +
D)`, `D` = 2^-18, through the shaft's end (radius `sqrt(9 + D^2)` rounded
once, crossing the shaft at an angle `D / 3`: snapped tangent to the shaft
and the top, its centre `(13/2, 41 + D)` and radius 3), the top's disc at
`41 + D + sqrt(9 + D^2)` (rounded once; its plane at `44 + D` on the
kernel's grid). Before the fix the kernel's results kept the stored
surfaces and edges beside the model's new edges, so the cut of the rib and
the body and their common did not sum to the rib (`OCC485a`'s by 4.7e-5),
the two orders of the fuse differed, a result symmetric in `y` had its
centre off `y = 0`, and boxes across the cylinders, the fillet and the
round failed the validator.

`boolean-turned-consistency-expected.tsv`: a `body NAME PATH volume V` row
(the model's volume by Pappus: its lines' frusta and annuli and each arc's
band by Green's theorem in the meridian, exactly but for pi), then per
partner `case NAME KIND ... volume V common C`: the partner (`rib X0 W0 X1
W1 Y0 Y1`, the prism over the triangle `(X0, W0)`, `(X1, W0)`, `(X0, W1)` in
`(x, z)` from `y = Y0` to `Y1`; `box X Y Z DX DY DZ`, `Solid::box_at`'s), its
volume and its common with the model (the integral over the heights of each
section's area inside the model's radius, by mpmath at 40 digits with every
kink of the integrand a breakpoint: the profile's nodes and each height
where the model's radius meets a corner's or an edge line's distance from
the axis). Every partner is symmetric in `y` but the last box.

`boolean-turned-consistency-crossings.tsv` (S9e.4b.4d.2c.2): the same
`body` row, then partners across both the fillet and the round, each
failing "an open Boolean of arcs in any position" before S9e.4b.4d.2c.2
(the arrangement's meetings found once per plane and per quadric, handing
the round the fillet's sections): `rib_high`, the known failure (the rib
over `(8, 34.5)`, `(20, 34.5)`, `(8, 42)`, its side planes across the
round's tangent circle); `box_both`, a box across the fillet's circles and
the round's lower one; and `post CX CY R Z0 Z1`, a rod parallel to the axis
about `(77/8, 0)` of radius 5/8 from 33 to 47, its common each height's
lens of two discs (the kinks where the model's radius is the rod's nearest
or farthest distance from the axis).

`boolean-turned-consistency-corners.tsv` (S9e.4b.4d.2c.1c): a `body` row
for `imported/corner_step.brep` (`generate_turned_corners_boolean_fixtures.py`'s
`step`: a cylinder of radius 4 to 3, a round of radius 1 about `(3, 3)`
crossing the shaft of radius 7/2 at the corner `3 + sqrt(3)/2`, the shaft
to 6), its volume by Pappus (the round's band between its tangent circle
and its corner in closed form), then partners across its corner ring: the
box `step_box` of that generator, a box off `y = 0` across the ring and the
round, and a rod parallel to the axis about `(7/2, 0)` of radius 1/4 from
9/4 to 5, each common by the same integral (the corner a kink of the model's
radius); every partner's plane clear of the tube's top and bottom, of the
planes `x` or `y` at a cylinder's radius and of the tube's equators (S9c's
tangencies).

The integrator is checked against the closed form first (a box holding the
model is the model, within 1e-25; so is a rod about the axis holding it,
and a lens is its chords' widths integrated). `--check` regenerates the texts and
fails where they differ, or where the body's file exists and its surfaces'
kinds are not the profile's (a near-flat cone, a cylinder, the chamfer's
cone, a near-flat cone, the fillet's torus, a cylinder, the round's torus, a
plane); `--write-bodies --occt-root PREFIX --sdk-manifest FILE` writes the body with
the pinned SDK.
"""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import sys

import mpmath as mp

ROOT = Path(__file__).resolve().parents[2]
FIXTURES = ROOT/'rust/fixtures'
PREFIX = 'boolean-turned-consistency'
BODIES_FILE = f'{PREFIX}-bodies.txt'
EXPECTED_FILE = f'{PREFIX}-expected.tsv'
CROSSINGS_FILE = f'{PREFIX}-crossings.tsv'
CORNERS_FILE = f'{PREFIX}-corners.tsv'

EPS = F(1, 2**27)
BOTTOM = F(1, 2**15)
STEP = F(1, 2**17)
D = F(1, 2**18)
TOLERANCE = '0.003'
mp.mp.dps = 40


def num(x):
    """A coordinate's shortest round-trip text."""
    x = F(x)
    return str(x.numerator) if x.denominator == 1 else repr(float(x))


def m(x):
    """A Fraction (or an mpf) as an mpf at the working precision."""
    if isinstance(x, mp.mpf):
        return x
    x = F(x)
    return mp.mpf(x.numerator)/x.denominator


def rounded(x):
    """`x` (an mpf) rounded once to binary64, as a Fraction."""
    return F(float(x))


# The round's radius and the top's height as OCCT is given them (rounded
# once).
ROUND = rounded(mp.sqrt(9 + m(D)**2))
TOP = rounded(m(41) + m(D) + mp.sqrt(9 + m(D)**2))


def occt_profile():
    """The profile OCCT revolves: `(r, w, segment)` to the next point, a
    segment `L` or `(A, cr, cw, radius, ccw)`."""
    return [
        (0, 0, 'L'),
        (24, BOTTOM, 'L'),
        (24, 32, 'L'),
        (21, 35, 'L'),
        (10, 35 + STEP, ('A', 10, F(71, 2) + STEP, F(1, 2), 0)),
        (F(19, 2), F(71, 2) + STEP, 'L'),
        (F(19, 2), 41, ('A', F(13, 2), 41 + D, ROUND, 1)),
        (F(13, 2), TOP, 'L'),
        (0, TOP, 'L'),
    ]


def world(r, w):
    """A meridian point on the tilted axis's frame, exactly."""
    return (F(r) + EPS*F(w), F(0), F(w) - EPS*F(r))


def bodies_text():
    rows = []
    for r, w, seg in occt_profile():
        x, y, z = world(r, w)
        text = f'{num(x)} {num(y)} {num(z)}'
        if seg == 'L':
            rows.append(f'{text} L')
        else:
            _, cr, cw, radius, ccw = seg
            cx, cy, cz = world(cr, cw)
            rows.append(f'{text} A {num(cx)} {num(cy)} {num(cz)} {num(radius)} {ccw}')
    wire = ' '.join(rows)
    return '\n'.join([
        'write ribbed imported/consistency_ribbed.brep',
        'plane 0 0 0 0 -1 0 1 0 0',
        f'wire S {len(rows)} {wire}',
        f'revolve 0 0 0 {num(EPS)} 0 1 6.283185307179586',
        f'tolerance {TOLERANCE}',
        'end',
    ]) + '\n'


# ------------------------------------------------------------------ the model

# The model's profile: segments `(w0, w1, kind, data)` from the bottom up:
# a line `r = r0 + k (w - w0)`, an arc `r = cr + s sqrt(t^2 - (w - cw)^2)`
# (`s` = -1 for the fillet, its material outside the tube; +1 for the round).
H = F(35)
TOP_MODEL = 44 + D
MODEL = [
    (BOTTOM, F(32), 'line', (F(24), F(0))),
    (F(32), H, 'line', (F(24), F(-1))),
    (H, H + F(1, 2), 'arc', (F(10), H + F(1, 2), F(1, 2), -1)),
    (H + F(1, 2), 41 + D, 'line', (F(19, 2), F(0))),
    (41 + D, TOP_MODEL, 'arc', (F(13, 2), 41 + D, F(3), 1)),
]


def radius(z, model=None):
    """The model's radius at height `z` (an mpf; zero outside it)."""
    for w0, w1, kind, data in model or MODEL:
        if m(w0) <= z <= m(w1):
            if kind == 'line':
                r0, k = data
                return m(r0) + m(k)*(z - m(w0))
            cr, cw, t, s = data
            d = m(t)**2 - (z - m(cw))**2
            return m(cr) + s*mp.sqrt(max(d, m(0)))
    return m(0)


def pappus():
    """The model's volume: each line's frustum, each arc's band (`pi int
    r^2 dw`, exactly: for `r = c + s sqrt(t^2 - u^2)` over a quarter,
    `c^2 t + s c pi t^2 / 2 + 2 t^3 / 3`)."""
    total = m(0)
    for w0, w1, kind, data in MODEL:
        if kind == 'line':
            r0, k = (m(x) for x in data)
            h = m(w1 - w0)
            r1 = r0 + k*h
            total += mp.pi*h*(r0*r0 + r0*r1 + r1*r1)/3
        else:
            c, _, t, s = data
            c, t = m(c), m(t)
            total += mp.pi*(c*c*t + s*c*mp.pi*t*t/2 + 2*t**3/3)
    return total


def segment_area(x0, x1, y0, y1, r):
    """The area of `[x0, x1] x [y0, y1]` inside the disc of radius `r`."""
    if r <= 0 or x1 <= x0 or y1 <= y0:
        return m(0)

    def width(y):
        if y*y >= r*r:
            return m(0)
        s = mp.sqrt(r*r - y*y)
        return max(m(0), min(x1, s) - max(x0, -s))
    cuts = {y0, y1}
    for x in (x0, x1):
        if r*r > x*x:
            for y in (mp.sqrt(r*r - x*x), -mp.sqrt(r*r - x*x)):
                cuts.add(y)
    cuts.update({r, -r, m(0)})
    points = sorted(c for c in cuts if y0 <= c <= y1)
    return mp.quad(width, points)


def roots(f, a, b, samples=64):
    """Every sign change of `f` on `[a, b]`, by bisection."""
    out = []
    xs = [a + (b - a)*i/samples for i in range(samples + 1)]
    for p, q in zip(xs, xs[1:]):
        fp, fq = f(p), f(q)
        if fp == 0:
            out.append(p)
        if fp*fq < 0:
            for _ in range(200):
                m = (p + q)/2
                fm = f(m)
                if fp*fm <= 0:
                    q = m
                else:
                    p, fp = m, fm
            out.append((p + q)/2)
    return out


def common(section, z0, z1, distances, model=None):
    """The volume of the partner (each height's section `x0, x1, y0, y1`)
    inside the model over `[z0, z1]`; `distances(z)` the distances from the
    axis where the section's area has a kink."""
    lo, hi = m(z0), m(z1)
    breaks = {lo, hi}
    for w0, w1, _, _ in model or MODEL:
        for w in (w0, w1):
            if lo < m(w) < hi:
                breaks.add(m(w))
    nodes = sorted(breaks)
    for a, b in zip(nodes, nodes[1:]):
        for g in distances:
            breaks.update(roots(lambda z: radius(z, model) - g(z), a, b))
    points = sorted(breaks)
    return mp.quad(lambda z: segment_area(*section(z), radius(z, model)), points)


def box_common(x, y, z, dx, dy, dz, model=None):
    x0, x1 = m(x), m(x) + m(dx)
    y0, y1 = m(y), m(y) + m(dy)
    distances = [lambda z, d=abs(c): d for c in (x0, x1, y0, y1)]
    distances += [lambda z, d=mp.sqrt(a*a + b*b): d for a in (x0, x1) for b in (y0, y1)]
    return common(lambda z: (x0, x1, y0, y1), z, z + dz, distances, model)


def rib_common(x0, w0, x1, w1, y0, y1):
    """The prism over the triangle `(x0, w0)`, `(x1, w0)`, `(x0, w1)` (its
    hypotenuse from `(x1, w0)` to `(x0, w1)`)."""
    x0, w0, x1, w1, y0, y1 = (m(v) for v in (x0, w0, x1, w1, y0, y1))

    def far(z):
        return x1 + (x0 - x1)*(z - w0)/(w1 - w0)
    distances = [lambda z, d=abs(c): d for c in (x0, y0, y1)]
    distances += [lambda z, d=mp.sqrt(x0*x0 + b*b): d for b in (y0, y1)]
    distances += [far] + [lambda z, b=b: mp.sqrt(far(z)**2 + b*b) for b in (y0, y1)]
    return common(lambda z: (x0, far(z), y0, y1), w0, w1, distances)


def lens(d, rho, r):
    """The area of the disc of radius `rho` about a point `d` from the axis
    inside the disc of radius `r` about the axis."""
    if r <= 0:
        return m(0)
    if d >= r+rho:
        return m(0)
    if d <= abs(r-rho):
        return mp.pi*min(r, rho)**2
    a = r*r*mp.acos((d*d+r*r-rho*rho)/(2*d*r))
    b = rho*rho*mp.acos((d*d+rho*rho-r*r)/(2*d*rho))
    return a+b-mp.sqrt((-d+r+rho)*(d+r-rho)*(d-r+rho)*(d+r+rho))/2


def post_common(cx, cy, rho, z0, z1, model=None):
    """A rod parallel to the axis about `(cx, cy)` of radius `rho` from `z0`
    to `z1` inside the model: each height's lens, its kinks where the
    model's radius is the rod's nearest or farthest distance from the axis."""
    d = mp.sqrt(m(cx)**2+m(cy)**2)
    rho = m(rho)
    lo, hi = m(z0), m(z1)
    breaks = {lo, hi}
    for w0, w1, _, _ in model or MODEL:
        for w in (w0, w1):
            if lo < m(w) < hi:
                breaks.add(m(w))
    nodes = sorted(breaks)
    for a, b in zip(nodes, nodes[1:]):
        for g in (d-rho, d+rho):
            breaks.update(roots(lambda z: radius(z, model)-g, a, b))
    return mp.quad(lambda z: lens(d, rho, radius(z, model)), sorted(breaks))


RIB = (8, F(69, 2), 20, F(83, 2), -1, 1)
BOXES = [
    ('base_x', (20, -1, 10, 10, 2, 5)),
    ('base_y', (-1, 20, 10, 2, 10, 5)),
    ('chamfer', (21, -1, 30, 4, 2, 4)),
    ('fillet', (9, -1, F(139, 4), 2, 2, F(11, 8))),
    ('shaft', (8, -1, 37, 4, 2, 2)),
    ('round', (8, -1, F(81, 2), 2, 2, F(5, 2))),
    ('round_skew', (7, 6, 40, 3, 3, F(9, 2))),
]


# S9e.4b.4d.2c.2's crossings (`boolean-turned-consistency-crossings.tsv`):
# partners across both the fillet and the round, the known failure of
# "S9e.4b.4d.2c.1a's consistency fixed" first (the rib whose side planes
# cross the round's tangent circle, its top at 42 over the shaft's end at
# 41), each failing "an open Boolean of arcs in any position" before
# S9e.4b.4d.2c.2: the box across the fillet's circles and the round's lower
# one (its top across the round), and a rod parallel to the axis about
# `(77/8, 0)` of radius 5/8 across the step's plane, the fillet, the shaft
# and the round (`post CX CY R Z0 Z1`).
CROSSING_RIB = ('rib_high', (8, F(69, 2), 20, 42, -1, 1))
CROSSING_BOX = ('box_both', (8, -1, 34, 3, 2, 9))
CROSSING_POST = ('post', (F(77, 8), 0, F(5, 8), 33, 47))


# S9e.4b.4d.2c.1c's corner body (`boolean-turned-consistency-corners.tsv`):
# `corner_step`'s model, its round's band from its tangent circle at 3 to
# its corner at `3 + sqrt(3)/2` (an mpf: `radius` and `common` take one).
CORNER = m(3) + mp.sqrt(3)/2
CORNER_MODEL = [
    (F(0), F(3), 'line', (F(4), F(0))),
    (F(3), CORNER, 'arc', (F(3), F(3), F(1), 1)),
    (CORNER, F(6), 'line', (F(7, 2), F(0))),
]
CORNER_BOXES = [
    ('step_box', (3, F(-1, 2), F(13, 4), F(3, 2), 1, F(5, 4))),
    ('box_off', (F(9, 4), F(9, 4), F(13, 4), 1, 1, F(5, 4))),
]
CORNER_POST = ('post', (F(7, 2), 0, F(1, 4), F(9, 4), 5))


def pappus_of(model):
    """A model's volume: each line's frustum, each arc's band between its
    heights in closed form (`pi int (c + s sqrt(t^2 - u^2))^2 du`, `u` the
    height off the centre)."""
    total = m(0)
    for w0, w1, kind, data in model:
        if kind == 'line':
            r0, k = (m(x) for x in data)
            h = m(w1) - m(w0)
            r1 = r0 + k*h
            total += mp.pi*h*(r0*r0 + r0*r1 + r1*r1)/3
        else:
            c, cw, t, s = data
            c, t = m(c), m(t)

            def prim(u):
                root = mp.sqrt(max(t*t - u*u, m(0)))
                return c*c*u + s*c*(u*root + t*t*mp.asin(u/t)) + t*t*u - u**3/3
            total += mp.pi*(prim(m(w1) - m(cw)) - prim(m(w0) - m(cw)))
    return total


def self_check():
    """The integrator against the closed form: a box holding the model is
    the model, by Pappus."""
    whole = box_common(-25, -25, -1, 50, 50, 46)
    assert abs(whole - pappus()) < mp.mpf(10)**-25, (whole, pappus())
    # A rod about the axis holding the model is the model; a lens against
    # its chords' widths integrated across the rod.
    whole = post_common(0, 0, 30, -1, 46)
    assert abs(whole - pappus()) < mp.mpf(10)**-25, (whole, pappus())
    # The corner body: a box holding it is it; its closed form against a
    # direct quadrature of its radius.
    whole = box_common(-5, -5, -1, 10, 10, 8, CORNER_MODEL)
    assert abs(whole - pappus_of(CORNER_MODEL)) < mp.mpf(10)**-25, (whole, pappus_of(CORNER_MODEL))
    direct = mp.pi*mp.quad(lambda z: radius(z, CORNER_MODEL)**2, [0, 3, CORNER, 6])
    assert abs(direct - pappus_of(CORNER_MODEL)) < mp.mpf(10)**-25
    d, rho = m(F(77, 8)), m(F(5, 8))
    for r in (m(9), m(F(19, 2)), m(10)):
        def width(x):
            s = rho*rho-(x-d)**2
            if s <= 0:
                return m(0)
            h = mp.sqrt(s)
            inside = r*r-x*x
            if inside <= 0:
                return m(0)
            return min(h, mp.sqrt(inside))*2
        cuts = sorted({d-rho, d+rho, r, *(c for c in [(r*r-rho*rho+d*d)/(2*d)] if d-rho < c < d+rho)})
        direct = mp.quad(width, cuts)
        assert abs(direct - lens(d, rho, r)) < mp.mpf(10)**-25, (r, direct, lens(d, rho, r))


# The stored surfaces' kinds in the file's order (OCCT's codes: 1 a plane,
# 2 a cylinder, 3 a cone, 5 a torus): the bottom's near-flat cone, the
# cylinder, the chamfer, the step's near-flat cone, the fillet, the shaft,
# the round, the top.
KINDS = ['3', '2', '3', '3', '5', '2', '5', '1']


def stored_kinds(path):
    """The surfaces' kinds of a `.brep` text, in order."""
    lines = path.read_text().splitlines()
    start = next(i for i, l in enumerate(lines) if l.startswith('Surfaces '))
    count = int(lines[start].split()[1])
    kinds, i = [], start + 1
    while len(kinds) < count:
        kind = lines[i].split()[0]
        kinds.append(kind)
        # A cone's half-angle and a cylinder's or a torus's radii follow on
        # its row; a cone's half-angle on the next.
        i += 2 if kind == '3' else 1
    return kinds


def expected_text():
    self_check()
    lines = [
        '# generated by rust/tools/generate_turned_consistency_fixtures.py',
        f'body\tribbed imported/consistency_ribbed.brep volume {mp.nstr(pappus(), 20)}',
    ]
    x0, w0, x1, w1, y0, y1 = RIB
    rib = F(1, 2)*(x1 - x0)*(w1 - w0)*(y1 - y0)
    spec = ' '.join(num(v) for v in RIB)
    lines.append(f'case\trib rib {spec} volume {mp.nstr(m(rib), 20)} '
                 f'common {mp.nstr(rib_common(*RIB), 20)}')
    for name, b in BOXES:
        volume = F(b[3])*F(b[4])*F(b[5])
        spec = ' '.join(num(v) for v in b)
        lines.append(f'case\t{name} box {spec} volume {mp.nstr(m(volume), 20)} '
                     f'common {mp.nstr(box_common(*b), 20)}')
    return '\n'.join(lines) + '\n'


def crossings_text():
    """S9e.4b.4d.2c.2's partners of the body, in `expected_text`'s rows."""
    self_check()
    lines = [
        '# generated by rust/tools/generate_turned_consistency_fixtures.py (S9e.4b.4d.2c.2)',
        f'body\tribbed imported/consistency_ribbed.brep volume {mp.nstr(pappus(), 20)}',
    ]
    name, r = CROSSING_RIB
    x0, w0, x1, w1, y0, y1 = r
    rib = F(1, 2)*(x1 - x0)*(w1 - w0)*(y1 - y0)
    lines.append(f'case\t{name} rib {" ".join(num(v) for v in r)} volume {mp.nstr(m(rib), 20)} '
                 f'common {mp.nstr(rib_common(*r), 20)}')
    name, b = CROSSING_BOX
    volume = F(b[3])*F(b[4])*F(b[5])
    lines.append(f'case\t{name} box {" ".join(num(v) for v in b)} volume {mp.nstr(m(volume), 20)} '
                 f'common {mp.nstr(box_common(*b), 20)}')
    name, p = CROSSING_POST
    cx, cy, rho, z0, z1 = p
    volume = mp.pi*m(rho)**2*m(z1-z0)
    lines.append(f'case\t{name} post {" ".join(num(v) for v in p)} volume {mp.nstr(volume, 20)} '
                 f'common {mp.nstr(post_common(*p), 20)}')
    return '\n'.join(lines) + '\n'


def corners_text():
    """S9e.4b.4d.2c.1c's corner body and its partners, in `expected_text`'s
    rows."""
    self_check()
    lines = [
        '# generated by rust/tools/generate_turned_consistency_fixtures.py (S9e.4b.4d.2c.1c)',
        f'body\tstep imported/corner_step.brep volume {mp.nstr(pappus_of(CORNER_MODEL), 20)}',
    ]
    for name, b in CORNER_BOXES:
        volume = F(b[3])*F(b[4])*F(b[5])
        lines.append(f'case\t{name} box {" ".join(num(v) for v in b)} volume {mp.nstr(m(volume), 20)} '
                     f'common {mp.nstr(box_common(*b, CORNER_MODEL), 20)}')
    name, p = CORNER_POST
    cx, cy, rho, z0, z1 = p
    volume = mp.pi*m(rho)**2*m(z1-z0)
    lines.append(f'case\t{name} post {" ".join(num(v) for v in p)} volume {mp.nstr(volume, 20)} '
                 f'common {mp.nstr(post_common(*p, CORNER_MODEL), 20)}')
    return '\n'.join(lines) + '\n'


def write_bodies(args):
    sys.path.insert(0, str(ROOT/'rust/tools'))
    import compare_boolean as base
    from compare_degree_elevation import verify_sdk
    output = (ROOT/'target/boolean-turned-consistency-oracle').resolve()
    output.mkdir(parents=True, exist_ok=True)
    prefix = args.occt_root.resolve()
    verify_sdk(prefix, args.sdk_manifest)
    executable, env, _, _ = base.build(prefix, output, base.SOURCE_FILE, 'boolean-oracle',
                                       ('TKBO', 'TKPrim', 'TKShHealing'),
                                       ['TKBO', 'TKPrim', 'TKShHealing']+base.TOOLKITS)
    env = dict(env, OCCT_BOOLEAN_FIXTURES=str(FIXTURES))
    (FIXTURES/'imported').mkdir(exist_ok=True)
    record = base.run(executable, (FIXTURES/BODIES_FILE).read_text(), env)
    written = [w for w in record['stdout'].splitlines() if w.endswith(' written')]
    if record['exit_code'] != 0 or len(written) != 1:
        raise SystemExit('the body was not written: '+json.dumps(record)[:2000])
    print(f'{len(written)} body written')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--write-bodies', action='store_true')
    parser.add_argument('--occt-root', type=Path)
    parser.add_argument('--sdk-manifest', type=Path)
    args = parser.parse_args()
    if args.write_bodies:
        if not (args.occt_root and args.sdk_manifest):
            parser.error('--write-bodies needs --occt-root and --sdk-manifest')
        write_bodies(args)
        return
    texts = {BODIES_FILE: bodies_text(), EXPECTED_FILE: expected_text(), CROSSINGS_FILE: crossings_text(),
             CORNERS_FILE: corners_text()}
    stale = []
    body = FIXTURES/'imported/consistency_ribbed.brep'
    if body.exists() and stored_kinds(body) != KINDS:
        stale.append(f'{body.name} ({stored_kinds(body)})')
    for name, text in texts.items():
        path = FIXTURES/name
        if args.check:
            if not path.exists() or path.read_text() != text:
                stale.append(name)
        else:
            path.write_text(text)
    if stale:
        raise SystemExit(f'stale: {", ".join(stale)}')
    print('checked' if args.check else 'written', ', '.join(texts))


if __name__ == '__main__':
    main()
