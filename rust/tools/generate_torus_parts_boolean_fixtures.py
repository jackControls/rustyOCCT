#!/usr/bin/env python3
"""Fixtures for S9d.4c of REVIEW_NOTES.md: Booleans of a sphere's cap or zone
against a whole torus, and of a torus v-segment or wedge against a prism
with arcs, a sphere (whole or a cap), a cone or frustum and a whole torus,
before any of its kernel code.

`boolean-torus-parts-cases.txt` lists each case in the Boolean protocol
(`identity_reference.encode_boolean_case`: a torus's block is its `frame` and
`torus MAJOR MINOR LOW HIGH ANGLE` rows, a sphere's `frame` and `sphere R LOW
HIGH`, a cone's `frame` and `cone BOTTOM TOP HEIGHT`);
`boolean-torus-parts-expected.tsv` gives per case `expect KIND S9d.4c`
(`solid`, `empty` or `degenerate`, then `reason TEXT`) and `result N volume
area cx cy cz` or `empty` from `torus_parts_boolean_reference.py`, as
S9d.4b.2's; `boolean-torus-parts-frames.tsv` the stored axes' bits. Frames
are the curved generator's.

`T` is the torus of radii 5/2 and 1 on `XY` at the origin (its hole of
radius 3/2, its outer equator of radius 7/2); `O` its outer half (latitudes
`-pi/2` to `pi/2`, the solid from the axis to the tube's outer side, its
end discs of radius 5/2 tangent to the torus along its top and bottom
circles), `I` its inner half (`pi/2` to `3 pi/2`, inside out), `B` the band
from latitude 1/2 to 9/4 (its rims of surd radius `R +- sqrt(r^2 - z^2)`),
`Q`, `H` and `W` the wedges of a quarter, a half and three quarters of a
turn (their ends' half-planes on the rounded directions `(cos, sin)` of the
turn). Caps and zones against `T`: a zone of radius 11/4 about the centre
between latitudes -1/2 and 1/2 meeting the tube in two circles
(`zone_coax`, all three), a cap of that sphere above latitude 1/4 whose
disc cuts the tube (`cap_coax`, and `coax_cap` the cap first), a cap of
radius 3/4 above latitude 1/2 on the tube's top, its rim of surd radius
across the tube (`cap_top`, all three), a zone about `x` (`SIDE`) across
the outer side, its disc cutting the tube and its rim crossing the disc's
spiric loop (`zone_side`), a hemisphere in `LEAN` on the tube's top (its
rim on a basis of unequal lengths, `dome_lean`) and a zone against the
torus in `TILT` (`zone_tilt`). Parts: `O` and `I` against coaxial pipes
(`oh_pipe`, all three; `ih_pipe` cutting its top disc), `I` against a
coaxial cone and sphere (`ih_cone`, `ih_ball`), `O` against a sphere
about a point of its top rim (`oh_ball`), a coaxial lower hemisphere
(`oh_dome`) and a small torus about `x` ringing the tube's outer side
(`oh_torus`, all three); `B` against a sphere across its lower rim
(`band_ball`, all three), a frustum across its upper rim (`band_cone`) and
the ringing torus across its lower rim (`band_torus`), and in `TILT`
against a coaxial pipe (`band_tilt_pipe`); `Q` against a pipe across its
end disc (`qw_pipe`, all three) and a small torus of a parallel axis across
it (`qw_torus`), `H` against a sphere across its start disc and rim (`hw_ball`),
`W` against a frustum across its end disc (`tw_cone`), and `Q` in `LEAN`
against a sphere across its start disc's rim (`qw_lean_ball`). `degenerate` with
reasons: a cap's rim within the resolution of tangency to the torus's outer
equator (`rim_touch`) and `Q`'s start rim touching a sphere at its outer
equator's point (`qw_rim_touch`: the sphere's centre off the rim's plane,
the surfaces crossing there).

Before writing, `reference_checks` compares the reference with independent
results (limits relative to the case's size):

* closed forms of the coaxial pairs by one quadrature along the torus's
  axis in its ideal frame: each slice's radial intervals of both inputs
  (the torus's annulus, a segment's disc or annulus, a wedge's annulus over
  its turn; a pipe's, a cap's or zone's, a cone's discs, a torus's annulus),
  the common by their intersections (volume, moments, the walls' areas
  inside the other by their own elements: the tori's by their latitudes, a
  sphere's by Archimedes, a cone's by `rho sqrt(1 + k^2)`, the discs by
  their radial intervals, a wedge's meridian discs by their chords), every
  other operation from the inputs' closed forms; within 1e-30 in exact
  frames, 1e-15 in turned ones;
* two-way checks: a half and the other half (the rest of the tube's circle
  between the same end planes, weighted by their orientation: the inner
  half inside out) against the tool as the whole torus against it
  (S9d.4b.2's reference), and a cap or zone and the rest of its sphere
  against the torus as the whole sphere, the common's volume and moments
  within 1e-30;
* every operation's volume, moments and area two ways (each face's two
  families of curves), both inputs' measures from their faces against
  their closed forms, inclusion and exclusion, the area identity, every
  face's two families' totals equal, within 1e-30; every quadrature node's
  structure its interval's; the declared degenerate pairs within 1e-18;
* Monte-Carlo estimates (200,000 uniform points per pair) of every volume
  and centre within 5 standard errors;
* the solids by the sweep of the torus's (the part's) normal slices, the
  fuse's by the rule (one where the inputs overlap);
* a scan for near coincidences (as S9d.4b.2's) must find none but in the
  declared degenerate pairs; every other pair's margins (S9d.4b.2's with
  every edge, the rims too) are reported, and every plane of one input
  must lie at least 1e-3 of the case's size from tangency to the other's
  spheres and tori (the kernel refuses a plane tangent to a sphere or a
  torus wherever it touches, S9d.1's and S9d.4a's rules).

Pairs run in worker processes (`--workers`).
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
import math
import os
from pathlib import Path
import random
import struct
from types import SimpleNamespace
import zlib

import mpmath as mp

from brep_reference import cos_rn, sin_rn
from identity_reference import Boundary, Case, encode_boolean_case
from curve_surface_reference import stored_axes
import torus_curved_boolean_reference as tc
import torus_parts_boolean_reference as ref
from generate_curved_boolean_fixtures import FRAMES
from generate_sphere_boolean_fixtures import frame_name, ideal_axes
from generate_cones_boolean_fixtures import roots_on

ROOT = Path(__file__).resolve().parents[1]
PREFIX = 'boolean-torus-parts'
STEP = 'S9d.4c'
BOOLEAN_OPERATION = 94
EXACT = {'XY', 'SIDE', 'DOWN', 'TURN'}
TWO_PI = tc.TWO_PI
HP = math.pi/2
Z = mp.mpf(0)


def at(name, origin):
    return tuple(float(c) for c in origin)+FRAMES[name]


def circle(cx, cy, r):
    return Boundary(circle=(float(cx), float(cy), float(r)))


def make(spec, op):
    kind = spec[0]
    if kind == 'torus':
        _, R, r, frame, low, high, angle = spec
        return Case('', 1e-7, op, frame, 0.0, 0.0, [], torus=(float(R), float(r), low, high, angle))
    if kind == 'sphere':
        _, r, frame, low, high = spec
        return Case('', 1e-7, op, frame, 0.0, 0.0, [], sphere=(float(r), low, high))
    if kind == 'cone':
        _, r0, r1, h, frame = spec
        return Case('', 1e-7, op, frame, 0.0, 0.0, [], cone=(float(r0), float(r1), float(h)))
    _, boundaries, frame, start, end = spec
    return Case('', 1e-7, op, frame, float(start), float(end), boundaries)


def torus(R, r, frame, low=0.0, high=TWO_PI, angle=TWO_PI):
    return ('torus', R, r, frame, low, high, angle)


def sphere(r, frame, low=-HP, high=HP):
    return ('sphere', r, frame, low, high)


def cone(r0, r1, h, frame):
    return ('cone', r0, r1, h, frame)


def prism(boundaries, frame, start, end):
    return ('prism', boundaries, frame, start, end)


def spec_frame(spec):
    return {'torus': 3, 'sphere': 2, 'cone': 4, 'prism': 2}[spec[0]]


def whole(spec):
    if spec[0] == 'torus':
        return spec[5]-spec[4] == TWO_PI and spec[6] == TWO_PI
    if spec[0] == 'sphere':
        return spec[3] == -HP and spec[4] == HP
    return True


class Boolean:
    def __init__(self, name, operation, obj, tool, form, kind='solid', reason=None):
        self.name, self.operation, self.kind, self.reason = name, operation, kind, reason
        self.form = form
        self.pair_name = name.rsplit('_', 1)[0]
        self.obj = make(obj, 91)
        self.obj.name = name
        self.tool = make(tool, 92)
        self.tool.name = name
        self.specs = (obj, tool)
        self.frames = (obj[spec_frame(obj)], tool[spec_frame(tool)])

    def encode(self):
        return encode_boolean_case(self.obj, self.operation, self.tool, BOOLEAN_OPERATION)


def group(name, obj, tool, ops, form):
    out = []
    for op, kind in ops.items():
        reason = None
        if isinstance(kind, tuple):
            kind, reason = kind
        out.append(Boolean(f'{name}_{op}', op, obj, tool, form, kind, reason))
    return out


RIM_TOUCH = "a sphere's circle within the resolution of tangency to a torus"
RIM_BALL = "a torus part's rim tangent to the other input's surface"
ALL = {'fuse': 'solid', 'cut': 'solid', 'common': 'solid'}
CC = {'cut': 'solid', 'common': 'solid'}
# The pairs whose checks include a complement: a half's other half (a
# band's rest of the tube bounds its region between the arc and the axis
# with a boundary crossing itself, whose parity set is no signed
# complement), a cap's rest of its sphere.
COMPLEMENTS = {'oh_pipe', 'oh_ball', 'ih_ball', 'cap_top', 'zone_side'}


def rim_touch_centre(a, lat):
    """A cap of radius `a` below latitude `lat` whose rim lies on the plane
    `w = 0` (the centre's height `a sin(lat)` rounded once, the kernel's rim
    `o + h n` exactly on it) and touches the circle of radius 7/2 there
    within rounding: its centre's `x` is `7/2 - sqrt(a^2 - h^2)` rounded."""
    h = a*sin_rn(lat)
    assert h == a*math.sin(lat)
    rho = mp.sqrt(mp.mpf(a)**2-mp.mpf(h)**2)
    return float(mp.mpf(3.5)-rho), h


def cases():
    """S9d.4c's pairs (see the module's docstring)."""
    O0 = (0, 0, 0)
    XY = at('XY', O0)
    T = torus(2.5, 1, XY)
    PI = math.pi
    O = torus(2.5, 1, XY, -HP, HP)
    I = torus(2.5, 1, XY, HP, 3*HP)
    B = torus(2.5, 1, XY, 0.5, 2.25)
    Q = torus(2.5, 1, XY, 0.0, TWO_PI, HP)
    H = torus(2.5, 1, XY, 0.0, TWO_PI, PI)
    W = torus(2.5, 1, XY, 0.0, TWO_PI, 3*HP)
    ring = torus(1.25, 0.5, at('SIDE', (0, 2.5, 0)))
    out = []
    # Caps and zones against the whole torus.
    out += group('zone_coax', T, sphere(2.75, XY, -0.5, 0.5), ALL, 'coaxial')
    cap = sphere(2.75, XY, 0.25, HP)
    out += group('cap_coax', T, cap, CC, 'coaxial')
    out += group('coax_cap', cap, T, {'cut': 'solid'}, 'coaxial')
    out += group('cap_top', T, sphere(0.75, at('XY', (2.5, 0, 0.5)), 0.5, HP), ALL, None)
    out += group('zone_side', T, sphere(1.25, at('SIDE', (3.3, 0, 0)), -0.5, 0.5), CC, None)
    out += group('dome_lean', T, sphere(0.75, at('LEAN', (2.5, 0, 0.75)), 0.0, HP), CC, None)
    out += group('zone_tilt', torus(2.5, 1, at('TILT', O0)), sphere(1, at('XY', (2.5, -0.5, 1)), -0.5, 0.5), CC,
                 None)
    x0, h = rim_touch_centre(0.75, 0.5)
    out += group('rim_touch', T, sphere(0.75, at('XY', (x0, 0, h)), -0.5, HP), {'fuse': ('degenerate', RIM_TOUCH)},
                 None)
    # The halves.
    out += group('oh_pipe', O, prism([circle(0, 0, 3)], XY, -2, 2), ALL, 'coaxial')
    out += group('ih_pipe', I, prism([circle(0, 0, 2)], XY, -0.5, 2), CC, 'coaxial')
    out += group('ih_cone', I, cone(4, 0, 4, at('XY', (0, 0, -2))), CC, 'coaxial')
    out += group('ih_ball', I, sphere(2, XY), CC, 'coaxial')
    out += group('oh_ball', O, sphere(0.6, at('XY', (2.5, 0, 1))), CC, None)
    out += group('oh_dome', O, sphere(3, at('XY', (0, 0, 0.25)), -HP, 0.0), CC, 'coaxial')
    out += group('oh_torus', O, torus(0.9, 0.35, at('SIDE', (0, 3.3, 0))), ALL, None)
    # The band.
    out += group('band_ball', B, sphere(0.5, at('XY', (3.3, 0, 0.45))), ALL, None)
    out += group('band_cone', B, cone(0.5, 0.25, 2, at('XY', (1.9, 0, -0.5))), CC, None)
    out += group('band_torus', B, ring, CC, None)
    TILT = at('TILT', O0)
    out += group('band_tilt_pipe', torus(2.5, 1, TILT, 0.5, 2.25), prism([circle(0, 0, 3)], TILT, -2, 2), CC,
                 'coaxial')
    # The wedges.
    out += group('qw_pipe', Q, prism([circle(0, 2.5, 0.5)], XY, -2, 2), ALL, None)
    out += group('qw_torus', Q, torus(0.75, 0.3, at('XY', (0.25, 2.45, 0.5))), CC, None)
    out += group('hw_ball', H, sphere(0.6, at('XY', (2.5, 0, 0.75))), CC, None)
    out += group('tw_cone', W, cone(0.75, 0.25, 3, at('XY', (0.3, -2.4, -1.5))), CC, None)
    out += group('qw_rim_touch', Q, sphere(0.625, at('XY', (4, 0.375, 0))), {'common': ('degenerate', RIM_BALL)},
                 None)
    out += group('qw_lean_ball', torus(2.5, 1, at('LEAN', O0), 0.0, TWO_PI, HP),
                 sphere(0.75, at('XY', (-0.4, 3.2, 0.25))), CC, None)
    return out


def exact_pair(case):
    """Whether the closed forms hold exactly: no input in a turned frame (a
    whole sphere's model does not depend on its frame's axes)."""
    return all(frame_name(f) in EXACT for f, s in zip(case.frames, case.specs)
               if not (s[0] == 'sphere' and whole(s)))


def check_latitudes(listed):
    """Every latitude's sine and cosine as the kernel's on this host
    (`scaled_sin`, `scaled_cos`: the platform's; the model rounds once)."""
    for c in listed:
        for s in c.specs:
            if s[0] == 'torus':
                lats = s[4:6]
            elif s[0] == 'sphere':
                lats = s[3:5]
            else:
                continue
            for lat in lats:
                assert math.sin(lat) == sin_rn(lat) and math.cos(lat) == cos_rn(lat), (c.name, lat)


# ------------------------------------------------------------------ closed forms: coaxial pairs

class Ideal:
    """An input's ideal frame: origin and orthonormal axes (mpf)."""

    def __init__(self, frame):
        self.o = tuple(mp.mpf(v) for v in frame[:3])
        self.x, self.y, self.n = ideal_axes(frame_name(frame))

    def world(self, u, v, w):
        return tuple(self.o[i]+u*self.x[i]+v*self.y[i]+w*self.n[i] for i in range(3))

    def local(self, X):
        d = [X[i]-self.o[i] for i in range(3)]
        return tuple(sum(d[i]*a[i] for i in range(3)) for a in (self.x, self.y, self.n))


def dot(a, b):
    return sum(x*y for x, y in zip(a, b))


class RadialTorus:
    """A torus (whole or a part) about its axis: its radial intervals at a
    level, its turn, its wall's latitudes and its end discs."""

    def __init__(self, spec):
        _, R, r, frame, low, high, angle = spec
        self.fr = Ideal(frame)
        self.R, self.r = mp.mpf(R), mp.mpf(r)
        self.turn = 2*mp.pi
        self.arc = (-mp.pi, mp.pi)
        self.ends = []
        self.segment = self.wedge = False
        if angle != TWO_PI:
            self.wedge = True
            self.turn = mp.atan2(mp.mpf(sin_rn(angle)), mp.mpf(cos_rn(angle))) % (2*mp.pi)
        elif high-low != TWO_PI:
            self.segment = True
            vs = []
            for lat in (low, high):
                z = mp.mpf(r*sin_rn(lat))
                side = 1 if cos_rn(lat) > 0 else -1
                q = mp.sqrt(max(self.r**2-z*z, Z))
                v = mp.atan2(z, side*q)
                v += mp.nint((mp.mpf(lat)-v)/(2*mp.pi))*2*mp.pi
                vs.append(v)
                self.ends.append((z, self.R+side*q))
            self.arc = tuple(vs)
        self.levels = [-self.r, self.r]+[z for z, _ in self.ends]

    def on_arc(self, v):
        lo, hi = self.arc
        k = mp.ceil((lo-v)/(2*mp.pi))
        return v+k*2*mp.pi <= hi

    def intervals(self, s):
        if abs(s) >= self.r:
            return []
        q = mp.sqrt(self.r**2-s*s)
        if not self.segment:
            return [(self.R-q, self.R+q)]
        a = mp.asin(s/self.r)
        sides = [sd for sd, ang in ((1, a), (-1, mp.pi-a)) if self.on_arc(ang)]
        if len(sides) == 1:
            return [(Z, self.R+sides[0]*q)]
        if len(sides) == 2:
            return [(self.R-q, self.R+q)]
        return []

    def radii(self):
        q = lambda s: mp.sqrt(max(self.r**2-s*s, Z))
        return [lambda s: self.R+q(s), lambda s: self.R-q(s)]

    def measures(self):
        """Volume, chart moments and area, closed (Pappus on the meridian
        section: the arc's Green integrals `t^2 dw`, `t^2 w dw`, the wall's
        `r rho dphi`, the end discs)."""
        R, r = self.R, self.r
        lo, hi = self.arc
        if not self.segment:
            A = self.turn*(2*mp.pi*R*r)
            V = self.turn*mp.pi*r*r*R
            I2 = mp.pi*r*r*(R*R+r*r/4)
            mom = (I2*mp.sin(self.turn), I2*(1-mp.cos(self.turn)), Z)
            if self.wedge:
                A += 2*mp.pi*r*r
            return V, mom, A
        # t = R + r cos phi, w = r sin phi: V = pi int t^2 dw (signed), Mw =
        # pi int t^2 w dw, over the arc (the end segments add nothing).
        f = lambda p: (R+r*mp.cos(p))**2*r*mp.cos(p)
        g = lambda p: (R+r*mp.cos(p))**2*r*mp.sin(p)*r*mp.cos(p)
        V = mp.pi*mp.quad(f, [lo, hi])
        Mw = mp.pi*mp.quad(g, [lo, hi])
        sgn = 1 if V > 0 else -1
        A = 2*mp.pi*r*(R*(hi-lo)+r*(mp.sin(hi)-mp.sin(lo)))
        A += sum(mp.pi*rho*rho for _, rho in self.ends)
        return sgn*V, (Z, Z, sgn*Mw), A

    def wall_inside(self, inside, extra):
        """The wall's area inside `inside(rho, s)` (by latitudes)."""
        R, r = self.R, self.r
        lo, hi = self.arc
        pts = {lo, hi}
        for g in extra:
            pts |= set(roots_on(lambda p, g=g: g(R+r*mp.cos(p), r*mp.sin(p)), lo, hi, n=1500))
        for p in [-mp.pi/2, mp.pi/2, 3*mp.pi/2, -3*mp.pi/2]:
            if lo < p < hi:
                pts.add(p)
        f = lambda p: r*(R+r*mp.cos(p))*(1 if inside(R+r*mp.cos(p), r*mp.sin(p)) else 0)
        return self.turn*mp.quad(f, sorted(pts))


class RadialOther:
    """The coaxial other input: its radial intervals at a level (its axis on
    the torus's), its caps, its wall element."""

    def __init__(self, spec, K):
        self.kind = spec[0]
        frame = spec[spec_frame(spec)]
        fr = Ideal(frame)
        c = K.fr.local(fr.o)
        assert abs(c[0])+abs(c[1]) < mp.mpf(10)**-12, 'coaxial'
        sg = dot(fr.n, K.fr.n)
        assert abs(abs(sg)-1) < mp.mpf(10)**-12, 'coaxial'
        self.sg = 1 if sg > 0 else -1
        self.c = c[2]
        self.caps = []
        if self.kind == 'sphere':
            _, rad, _, low, high = spec
            self.rad = mp.mpf(rad)
            zlo = -self.rad if low == -HP else mp.mpf(rad*sin_rn(low))
            zhi = self.rad if high == HP else mp.mpf(rad*sin_rn(high))
            ends = sorted((self.c+self.sg*zlo, self.c+self.sg*zhi))
            self.range = tuple(ends)
            for z, cut in ((zlo, low != -HP), (zhi, high != HP)):
                if cut:
                    self.caps.append((self.c+self.sg*z, mp.sqrt(self.rad**2-z*z)))
        elif self.kind == 'torus':
            _, R2, r2, _, _, _, _ = spec
            self.R2, self.r2 = mp.mpf(R2), mp.mpf(r2)
            self.range = (self.c-self.r2, self.c+self.r2)
        elif self.kind == 'cone':
            _, r0, r1, h, _ = spec
            self.r0, self.r1, self.h = mp.mpf(r0), mp.mpf(r1), mp.mpf(h)
            self.k = (self.r1-self.r0)/self.h
            ends = (self.c, self.c+self.sg*self.h)
            self.range = (min(ends), max(ends))
            self.caps = [(ends[0], self.r0), (ends[1], self.r1)]
        else:
            _, boundaries, _, lo, hi = spec
            self.a = mp.mpf(boundaries[0].circle[2])
            ends = (self.c+self.sg*mp.mpf(lo), self.c+self.sg*mp.mpf(hi))
            self.range = (min(ends), max(ends))
            self.caps = [(e, self.a) for e in ends]

    def rho(self, s):
        if self.kind == 'sphere':
            return mp.sqrt(max(self.rad**2-(s-self.c)**2, Z))
        if self.kind == 'cone':
            return self.r0+self.k*self.sg*(s-self.c)
        return self.a

    def intervals(self, s):
        lo, hi = self.range
        if not lo < s < hi:
            return []
        if self.kind == 'torus':
            q = mp.sqrt(self.r2**2-(s-self.c)**2)
            return [(self.R2-q, self.R2+q)]
        return [(Z, self.rho(s))]

    def radii(self):
        if self.kind == 'torus':
            q = lambda s: mp.sqrt(max(self.r2**2-(s-self.c)**2, Z))
            return [lambda s: self.R2+q(s), lambda s: self.R2-q(s)]
        return [self.rho]

    def measures(self):
        """Volume, chart moment along the axis and area, closed."""
        if self.kind == 'sphere':
            lo, hi = sorted((self.range[0]-self.c, self.range[1]-self.c))
            R = self.rad
            V = mp.pi*((R*R*hi-hi**3/3)-(R*R*lo-lo**3/3))
            Mz = mp.pi*((R*R*hi**2/2-hi**4/4)-(R*R*lo**2/2-lo**4/4))
            A = 2*mp.pi*R*(hi-lo)+sum(mp.pi*rr*rr for _, rr in self.caps)
            return V, self.c*V+Mz, A
        if self.kind == 'torus':
            V = 2*mp.pi**2*self.R2*self.r2**2
            return V, self.c*V, 4*mp.pi**2*self.R2*self.r2
        if self.kind == 'cone':
            r0, r1, h = self.r0, self.r1, self.h
            V = mp.pi*h*(r0*r0+r0*r1+r1*r1)/3
            zc = h*(r0*r0+2*r0*r1+3*r1*r1)/(4*(r0*r0+r0*r1+r1*r1))
            A = mp.pi*(r0+r1)*mp.sqrt(h*h+(r1-r0)**2)+mp.pi*(r0*r0+r1*r1)
            return V, (self.c+self.sg*zc)*V, A
        lo, hi = self.range
        a = self.a
        V = mp.pi*a*a*(hi-lo)
        return V, (lo+hi)/2*V, 2*mp.pi*a*(hi-lo)+2*mp.pi*a*a

    def wall_inside(self, K, inside, extra):
        """Its lateral wall's area inside the torus's section (per the torus's
        turn)."""
        lo, hi = self.range
        if self.kind == 'torus':
            R2, r2, c = self.R2, self.r2, self.c
            pts = {-mp.pi, mp.pi, -mp.pi/2, mp.pi/2}
            for g in extra:
                pts |= set(roots_on(lambda p, g=g: g(R2+r2*mp.cos(p), c+r2*mp.sin(p)), -mp.pi, mp.pi, n=1500))
            f = lambda p: r2*(R2+r2*mp.cos(p))*(1 if inside(R2+r2*mp.cos(p), c+r2*mp.sin(p)) else 0)
            return K.turn*mp.quad(f, sorted(pts))
        pts = {lo, hi}
        for g in extra:
            pts |= set(roots_on(lambda s, g=g: g(self.rho(s), s), lo, hi, n=1500))
        pts |= {p for p in K.levels if lo < p < hi}
        if self.kind == 'sphere':
            el = lambda s: self.rad
        elif self.kind == 'cone':
            el = lambda s: self.rho(s)*mp.sqrt(1+self.k**2)
        else:
            el = lambda s: self.a
        f = lambda s: el(s)*(1 if inside(self.rho(s), s) else 0)
        return K.turn*mp.quad(f, sorted(pts))


def iv_meet(A, B):
    out = []
    for a0, a1 in A:
        for b0, b1 in B:
            lo, hi = max(a0, b0), min(a1, b1)
            if hi > lo:
                out.append((lo, hi))
    return out


def within(ivs, x):
    return any(a < x < b for a, b in ivs)


def coaxial_form(case):
    """The torus (whole or a part) against a coaxial input: every operation's
    volume, area and centre (world, ideal frames)."""
    tspec = next(s for s in case.specs if s[0] == 'torus' and not whole(s)) if any(
        s[0] == 'torus' and not whole(s) for s in case.specs) else next(s for s in case.specs if s[0] == 'torus')
    pspec = case.specs[1] if case.specs[0] is tspec else case.specs[0]
    k_first = case.specs[0] is tspec
    K = RadialTorus(tspec)
    P = RadialOther(pspec, K)
    # Kinks along the axis: the inputs' levels and the levels where a
    # boundary radius of one meets one of the other.
    lo = max(-K.r, P.range[0])
    hi = min(K.r, P.range[1])
    pts = {lo, hi}
    for z in K.levels+[c[0] for c in P.caps]:
        if lo < z < hi:
            pts.add(z)
    for f1 in K.radii():
        for f2 in P.radii():
            pts |= set(roots_on(lambda s, f1=f1, f2=f2: f1(s)-f2(s), lo, hi, n=1500))
    for f2 in P.radii():
        for z, rho in K.ends:
            pts |= set(roots_on(lambda s, f2=f2, rho=rho: f2(s)-rho, lo, hi, n=400))
    pts = sorted(pts)

    def g(s, k):
        ivs = iv_meet(K.intervals(s), P.intervals(s))
        if k == 0:
            return sum((b*b-a*a)/2 for a, b in ivs)
        if k == 1:
            return s*sum((b*b-a*a)/2 for a, b in ivs)
        return sum((b**3-a**3)/3 for a, b in ivs)
    if hi > lo:
        A0, A1, A2 = (mp.quad(lambda s: g(s, k), pts) for k in range(3))
    else:
        A0 = A1 = A2 = Z
    V = K.turn*A0
    mom = (A2*mp.sin(K.turn), A2*(1-mp.cos(K.turn)), K.turn*A1)
    # The walls inside the other, the discs by their radial intervals.
    kin = K.wall_inside(lambda rho, s: within(P.intervals(s), rho),
                        [lambda rho, s, f=f: rho-f(s) for f in P.radii()]
                        + [lambda rho, s, z=c[0]: s-z for c in P.caps])
    for z, rho in K.ends:
        kin += 2*mp.pi*sum((b*b-a*a)/2 for a, b in iv_meet([(Z, rho)], P.intervals(z)))
    if K.wedge:
        # Its meridian discs: the chord [R - q, R + q] inside the other at
        # each level, twice.
        kpts = sorted({p for p in pts if -K.r <= p <= K.r} | {-K.r, K.r})
        kin += 2*mp.quad(lambda s: sum(b-a for a, b in iv_meet(K.intervals(s), P.intervals(s))), kpts)
    pin = P.wall_inside(K, lambda rho, s: within(K.intervals(s), rho),
                        [lambda rho, s, f=f: rho-f(s) for f in K.radii()]
                        + [lambda rho, s, z=z: s-z for z in K.levels])
    for z, rad in P.caps:
        pin += K.turn*sum((b*b-a*a)/2 for a, b in iv_meet([(Z, rad)], K.intervals(z)))
    VK, MK, AK = K.measures()
    VP, MPw, AP = P.measures()
    MP = (Z, Z, MPw)
    to_world = lambda m: tuple(sum(m[j]*(K.fr.x, K.fr.y, K.fr.n)[j][i] for j in range(3)) for i in range(3))
    common = (V, mom, kin+pin)
    fuse = (VK+VP-V, tuple(a+b-c for a, b, c in zip(MK, MP, mom)), (AK-kin)+(AP-pin))
    if k_first:
        cut = (VK-V, tuple(a-c for a, c in zip(MK, mom)), (AK-kin)+pin)
    else:
        cut = (VP-V, tuple(a-c for a, c in zip(MP, mom)), (AP-pin)+kin)
    out = {}
    for op, (v, m, a) in (('common', common), ('fuse', fuse), ('cut', cut)):
        if v < mp.mpf(10)**-30:
            out[op] = None
            continue
        mw = to_world(m)
        out[op] = (v, a, tuple(K.fr.o[i]+mw[i]/v for i in range(3)))
    return out


def pair_forms(case):
    if case.form == 'coaxial':
        return coaxial_form(case)
    return None


# ------------------------------------------------------------------ two-way checks

def pair_of(A, B):
    """A reference pair of two built inputs."""
    p = ref.Pair.__new__(ref.Pair)
    p.A, p.B = A, B
    tori = [S for S in (A, B) if S.kind == 'torus']
    parts = [S for S in tori if ref.is_part(S)]
    p.K = parts[0] if parts else tori[0]
    p.P = B if p.K is A else A
    p.size = mp.mpf(tc.box_size([A.box(), B.box()]))
    p.c0 = p.K.fr.om
    p._sweeps = None
    p._counts = {}
    return p


def common_vm(pair):
    V, m, _ = pair.measures('common')
    return V, m


def complement_check(first):
    """A segment and its complement against the tool as the whole torus
    against it (weighted by their orientation), or a cap or zone and the
    rest of its sphere against the torus as the whole sphere: the common's
    volume and moments, the largest deviation relative to the size."""
    obj, tool = first.obj, first.tool
    A, B = ref.make_input(obj), ref.make_input(tool)
    if ref.is_part(A) or ref.is_part(B):
        k_first = ref.is_part(A)
        K = A if k_first else B
        O = B if k_first else A
        Kc = ref.TorusPart(K.case, ref.complement_part(K.part))
        parts = [(K.part.sigma, pair_of(K, O) if k_first else pair_of(O, K)),
                 (Kc.part.sigma, pair_of(Kc, O) if k_first else pair_of(O, Kc))]
        wcase = SimpleNamespace(**vars(K.case))
        R, r, _, _, _ = K.case.torus
        wcase.torus = (R, r, 0.0, TWO_PI, TWO_PI)
        whole_pair = tc.Pair(wcase, tool) if k_first else tc.Pair(obj, wcase)
    else:
        k_first = ref.is_cap(B)
        C = B if k_first else A
        T = A if k_first else B
        radius, low, high = C.case.sphere
        pieces = [(low, high)]+[p for p in ((-HP, low), (high, HP)) if p[1] > p[0]]
        parts = []
        for lo_, hi_ in pieces:
            c = SimpleNamespace(**vars(C.case))
            c.sphere = (radius, lo_, hi_)
            S = ref.Cap(c) if (lo_, hi_) != (-HP, HP) else tc.Sphere(c)
            parts.append((1, pair_of(T, S) if k_first else pair_of(S, T)))
        wcase = SimpleNamespace(**vars(C.case))
        wcase.sphere = (radius, -HP, HP)
        whole_pair = tc.Pair(T.case, wcase) if k_first else tc.Pair(wcase, T.case)
    Vw, mw = common_vm(whole_pair)
    Vs, ms = Z, [Z, Z, Z]
    for sgn, p in parts:
        V, m = common_vm(p)
        Vs += sgn*V
        ms = [a+sgn*b for a, b in zip(ms, m)]
    size = whole_pair.size
    return max(abs(Vs-Vw)/size**3, max(abs(a-b) for a, b in zip(ms, mw))/size**4)


# ------------------------------------------------------------------ plane tangencies

def plane_margin(pair):
    """The least distance (relative to the size) from tangency of every plane
    of one input to the other's spheres (the centre's distance against the
    radius) and tori (a linear function's critical values on the torus in
    closed form: on `o + (R + r cos v)(cos u x + sin u y) + r sin v n` it is
    `C + (R + r cos v)(alpha cos u + beta sin u) + r gamma sin v`, stationary
    at `C +- R rho +- r sqrt(rho^2 + gamma^2)`, `rho = |(alpha, beta)|`)."""
    out = mp.inf
    for S, O in ((pair.A, pair.B), (pair.B, pair.A)):
        for pl in (s for s in S.surfs if s.degree == 1):
            if O.kind not in ('sphere', 'torus'):
                continue
            g = pl.grad(O.fr.om)
            gn = tc.norm(g)
            if O.kind == 'sphere':
                out = min(out, abs(abs(pl.f(O.cm))/gn-O.Rm)/pair.size)
                continue
            fr = O.fr
            C = pl.f(fr.om)
            al, be, ga = dot(g, fr.xm), dot(g, fr.ym), dot(g, fr.nm)
            rho = mp.sqrt(al*al+be*be)
            for s1 in (1, -1):
                for s2 in (1, -1):
                    v = C+s1*O.Rm*rho+s2*O.rm*mp.sqrt(rho*rho+ga*ga)
                    out = min(out, abs(v)/gn/pair.size)
    return out


# ------------------------------------------------------------------ Monte Carlo

def monte_carlo(pair, n, seed):
    """Volume and centre estimates of the three operations and their
    standard errors, from `n` uniform points in the inputs' box."""
    fa, fb = ref.float_model(pair.A), ref.float_model(pair.B)
    ba, bb = pair.A.box(), pair.B.box()
    box0 = [min(ba[0][i], bb[0][i]) for i in range(3)]
    box1 = [max(ba[1][i], bb[1][i]) for i in range(3)]
    vbox = (box1[0]-box0[0])*(box1[1]-box0[1])*(box1[2]-box0[2])
    rng = random.Random(seed)
    sets = {'fuse': lambda a, b: a or b, 'cut': lambda a, b: a and not b, 'common': lambda a, b: a and b}
    acc = {op: [0, [0.0]*3, [0.0]*3] for op in ref.OPS}
    for _ in range(n):
        X = [box0[i]+(box1[i]-box0[i])*rng.random() for i in range(3)]
        a, b = fa.contains(X), fb.contains(X)
        for op in ref.OPS:
            if sets[op](a, b):
                e = acc[op]
                e[0] += 1
                for i in range(3):
                    e[1][i] += X[i]
                    e[2][i] += X[i]*X[i]
    out = {}
    for op, (k, s1, s2) in acc.items():
        p = k/n
        vol = vbox*p
        sv = vbox*math.sqrt(max(p*(1-p), 1e-300)/n)
        if k > 1:
            c = [s/k for s in s1]
            sc = [math.sqrt(max(s2[i]/k-c[i]*c[i], 0)/k) for i in range(3)]
        else:
            c, sc = None, None
        out[op] = (vol, sv, c, sc)
    return out


# ------------------------------------------------------------------ the pairs' work

def dev3(a, b, size):
    """Largest relative deviation of (volume, world moments, area)."""
    (va, ma, aa), (vb, mb, ab) = a, b
    return max(abs(va-vb)/size**3, max(abs(x-y) for x, y in zip(ma, mb))/size**4, abs(aa-ab)/size**2)


def evaluate(job):
    """One pair: the reference's rows for its operations and every check's
    deviations (run in a worker process)."""
    name, first, ops, mc_n = job
    obj, tool = first.obj, first.tool
    pair = ref.Pair(obj, tool)
    size = pair.size
    res = {op: pair.result(op) for op in ref.OPS}
    rows = {op: ref.rows(obj, op, tool, pair)[0] for op in ops}
    checks = {}
    closed = {'A': pair.A.closed(), 'B': pair.B.closed()}
    checks['inputs'] = max(dev3(pair.input_measures(role, k), closed[role], size) for role in 'AB' for k in (0, 1))
    checks['second_family'] = max(dev3(pair.measures(op, 0), pair.measures(op, 1), size) for op in ref.OPS)
    (VA, MA, AA), (VB, MB, AB) = closed['A'], closed['B']
    dev = Z
    for k in (0, 1):
        Vc, Mc, Ac = pair.measures('common', k)
        Vf, Mf, Af = pair.measures('fuse', k)
        Vx, Mx, _ = pair.measures('cut', k)
        dev = max(dev, dev3((Vf, Mf, Z), (VA+VB-Vc, tuple(a+b-c for a, b, c in zip(MA, MB, Mc)), Z), size),
                  dev3((Vx, Mx, Z), (VA-Vc, tuple(a-c for a, c in zip(MA, Mc)), Z), size))
        checks['area_identity'] = max(checks.get('area_identity', Z), abs(Af+Ac-AA-AB)/size**2)
    checks['inclusion_exclusion'] = dev
    face = Z
    sw = pair.sweeps()
    for (k, role, fname), s in sw.items():
        if k == 0:
            other = sw[(1, role, fname)]
            face = max(face, abs(sum(s.result[c][0] for c in s.result)-sum(other.result[c][0] for c in other.result)))
    checks['face_classes'] = face/size**2
    mcr = monte_carlo(pair, mc_n, zlib.crc32(name.encode()))
    z = 0.0
    for op in ref.OPS:
        vol, sv, c, sc = mcr[op]
        V, mm = pair.volume(op)
        v = float(V)
        z = max(z, abs(vol-v)/sv)
        if v > 1e-9 and c is not None:
            cent = [float(x/V) for x in mm]
            z = max(z, max(abs(c[i]-cent[i])/max(sc[i], 1e-12) for i in range(3)))
    checks['monte_carlo_sigma'] = mp.mpf(z)
    if name in COMPLEMENTS:
        checks['complement'] = complement_check(first)
    near = ref.near_coincidences(pair)
    counts = {k: pair.swept_count(k) for k in ('common', 'A-B', 'B-A', 'fuse')}
    overlap = pair.volume('common')[0] > mp.mpf(10)**-20*size**3
    rule = 1 if overlap else 2
    assert counts['fuse'] == rule, (name, 'the fuse swept', counts['fuse'], 'by the rule', rule)
    margins = ref.margins(pair)
    margins['plane_tangency'] = plane_margin(pair)
    stats = {'missed': pair.missed(), 'quadrature': pair.quad_error()/size**4, 'margins': margins,
             'solids': counts, 'events': {f'{k}{r}:{f}': len(s.events) for (k, r, f), s in sw.items()}}
    return name, rows, res, checks, near, stats


def run(jobs, fn, workers):
    if workers <= 1:
        return [fn(j) for j in jobs]
    with ProcessPoolExecutor(workers) as pool:
        return list(pool.map(fn, jobs))


def closed_job(first):
    return first.pair_name, pair_forms(first)


LIMITS = {'closed_forms_exact_frames': 1e-30, 'closed_forms_turned_frames': 1e-15, 'degenerate_pairs': 1e-18,
          'inputs': 1e-30, 'second_family': 1e-30, 'inclusion_exclusion': 1e-30, 'area_identity': 1e-30,
          'face_classes': 1e-30, 'complement': 1e-30, 'monte_carlo_sigma': 5, 'volume_quadrature_estimate': 1e-28}
PLANE_MARGIN = 1e-3


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
        kind = 'exact_frames' if exact_pair(firsts[name]) else 'turned_frames'
        for op, want in form.items():
            n, vol, area, centre = res[op]
            if want is None:
                assert n == 0, (name, op, 'the closed form is empty')
                continue
            V, A, C = want
            size = max(abs(V)**(mp.mpf(1)/3), 1)
            dev = max(abs(vol-V)/abs(V), abs(area-A)/abs(A), max(abs(centre[i]-C[i]) for i in range(3))/size)
            note(f'closed_forms_{kind}', dev)
    degenerate = {c.pair_name for c in cases() if c.kind == 'degenerate'}
    for name, rows, res, checks, near, stats in results:
        for key, value in checks.items():
            if name in degenerate and key != 'monte_carlo_sigma':
                note('degenerate_pairs', value)
            else:
                note(key, value)
        note('volume_quadrature_estimate', stats['quadrature'])
        if name not in degenerate:
            assert stats['missed'] == 0, (name, 'quadrature nodes off their interval\'s structure', stats['missed'])
            assert stats['margins']['plane_tangency'] > PLANE_MARGIN, (name, 'a plane near tangency',
                                                                      stats['margins']['plane_tangency'])
    return worst, covered


def generate(results):
    by_pair = {r[0]: r for r in results}
    blocks = []
    out = [f'# case\trow ({STEP}, torus_parts_boolean_reference.py: expect KIND {STEP}, reason TEXT for a '
           'degenerate case, then result N volume area cx cy cz or empty)']
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


def jobs(mc_n, only=None):
    pairs = {}
    for c in cases():
        pairs.setdefault(c.pair_name, [c, []])[1].append(c.operation)
    return [(name, first, ops, mc_n) for name, (first, ops) in pairs.items() if only is None or name in only]


def validate(listed):
    names = [c.name for c in listed]
    assert len(names) == len(set(names)), 'duplicate case names'
    for c in listed:
        assert c.obj.torus is not None or c.tool.torus is not None, f'{c.name}: a torus'
        assert not all(whole(s) for s in c.specs), f'{c.name}: a part or a cap'
        for f in c.frames:
            frame_name(f)
    check_latitudes(listed)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--workers', type=int, default=max(1, min(4, (os.cpu_count() or 2)-2)))
    parser.add_argument('--samples', type=int, default=200000, help='Monte-Carlo points per pair')
    parser.add_argument('--only', nargs='*', help='evaluate these pairs only (no files written)')
    args = parser.parse_args()
    listed = cases()
    validate(listed)
    if args.only:
        results = run(jobs(args.samples, set(args.only)), evaluate, args.workers)
        firsts = {c.pair_name: c for c in listed}
        for name, rows, res, checks, near, stats in results:
            print(name, rows, {k: mp.nstr(v, 3) for k, v in checks.items()}, near,
                  {k: (mp.nstr(v, 3) if not isinstance(v, (dict, tuple, int)) else v) for k, v in stats.items()})
            form = pair_forms(firsts[name])
            if form:
                for op, want in form.items():
                    n, vol, area, centre = res[op]
                    if want is None:
                        print('  form', op, 'empty; reference', n)
                    else:
                        print('  form', op, mp.nstr(vol-want[0], 3), mp.nstr(area-want[1], 3),
                              [mp.nstr(centre[i]-want[2][i], 3) for i in range(3)])
        return
    results = run(jobs(args.samples), evaluate, args.workers)
    degenerate = {c.pair_name for c in listed if c.kind == 'degenerate'}
    near = [f'{name}: near coincidences {n}' for name, _, _, _, n, _ in results if n and name not in degenerate]
    if near:
        raise SystemExit('\n'.join(near))
    firsts = {}
    for c in listed:
        firsts.setdefault(c.pair_name, c)
    forms = dict(run(list(firsts.values()), closed_job, args.workers))
    worst, covered = reference_checks(results, forms)
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
    kinds = {}
    for c in listed:
        kinds[c.kind] = kinds.get(c.kind, 0)+1
    print(len(listed), 'cases:', ', '.join(f'{sum(1 for c in listed if c.operation == op)} {op}'
                                           for op in ('fuse', 'cut', 'common')),
          '-', ', '.join(f'{v} {k}' for k, v in sorted(kinds.items())),
          '-', f'{sum(1 for c in listed if exact_pair(c))} in exact frames')
    print('reference checks (largest deviation, relative to the case size):',
          ', '.join(f'{k} {mp.nstr(v, 3)}' for k, v in sorted(worst.items())))
    print('cases per check:', ', '.join(f'{k} {v}' for k, v in sorted(covered.items())))
    print('closed forms:', sum(1 for f in forms.values() if f is not None), 'of', len(forms), 'pairs')
    by_pair = {r[0]: r for r in results}
    least = {}
    for name, _, _, _, n, stats in results:
        if n:
            print('declared near coincidences:', name, '; '.join(n[:4]))
        if name not in degenerate:
            for k, v in stats['margins'].items():
                if v < least.get(k, (mp.inf, None))[0]:
                    least[k] = (v, name)
    print('least margins (non-degenerate pairs):', ', '.join(f'{k} {mp.nstr(v, 3)} ({n})'
                                                             for k, (v, n) in sorted(least.items())))
    for name, _, _, _, _, stats in results:
        print('solids', name, stats['solids'], 'events', sum(stats['events'].values()))
    for c in listed:
        print(c.name, by_pair[c.pair_name][1][c.operation][0])


if __name__ == '__main__':
    main()
