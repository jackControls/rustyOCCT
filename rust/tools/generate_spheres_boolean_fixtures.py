#!/usr/bin/env python3
"""Fixtures for S9d.2 of REVIEW_NOTES.md: Booleans of a sphere (whole or a
cap, `Solid::sphere_with`) against a prism whose profile holds arcs and
circles (cylindrical walls), and of two spheres, before any of its kernel
code.

`boolean-spheres-cases.txt` lists each case in the Boolean protocol
(`identity_reference.encode_boolean_case`: a sphere's block is its `frame`
and `sphere R LOW HIGH` rows, either input or both); `boolean-spheres-
expected.tsv` gives per case `expect KIND S9d.2` (`solid`, `empty` or
`degenerate`, then `reason TEXT`) and `result N volume area cx cy cz` or
`empty` from `spheres_boolean_reference.py`, as the curved fixtures;
`boolean-spheres-frames.tsv` the stored axes' bits. Frames are the curved
generator's (`FRAMES`, their stored axes the kernel's bit for bit). A cap
only in an exact frame centred at the world's origin, its axis the prism's
(S9d.1's rule; the reference slices along it).

Pairs. Two spheres: crossing (`spheres_cross`, all three operations), a
small sphere inside a big one (`spheres_nested`, the cut a cavity), apart
(`spheres_apart`: the fuse two solids, the common empty), crossing in
turned frames (`spheres_turned`). A sphere and a cylinder with the centre on
its axis, meeting in parallels at rational heights: a pipe through a
sphere (`pipe_ring`: the cut the spherical ring, `pi h^3 / 6`;
`pipe_ends`, the pipe less the sphere, two solids), in `SIDE`
(`pipe_side`), a sphere against a box with a coaxial round hole whose walls
cut it (`hole_box`: the sphere less the box its core and four bulges, five
solids; `holed_box`, the box less the sphere), a hemisphere against a pipe
through its disc (`dome_pipe`). A sphere and a cylinder off its centre (the
quartic): a rod through the sphere, two rings (`rings`, `rod_ends` two
solids; `rings_tilt` in `TILT`), a partial bite, one loop (`bite`, all
three operations; `bite_tilt` in `TILT`), a cylinder ending inside the
sphere, its cap circle crossing it (`cap_cross`). A stadium through a
sphere (`stadium`: the sphere less it two solids; `stadium_ends`). And
`degenerate` with reasons: two spheres tangent (`spheres_tangent`), a
cylinder tangent to the sphere outside (`rod_tangent`) and inside
(`rod_inside_tangent`).

Before writing, `reference_checks` compares the reference with independent
results (limits relative to the case's size):

* closed forms where available, by textbook formulas or one quadrature
  (mpmath's tanh-sinh between the kinks) of exact 2D forms: two spheres'
  lens as the sum of two caps (`pi h^2 (3R - h)/3`, `2 pi R h`); a sphere
  against a coaxial cylinder by the axial integrals of `pi min(a^2, R^2 -
  z^2)` (a ring `pi h^3 / 6` whatever the radius), the sphere's face inside
  by Archimedes, the wall inside `2 pi a` times its height inside; a box
  with a coaxial hole as S9d.1's box less the coaxial core; an off-axis
  cylinder by the lens of two discs (areas, first moments by circular
  segments, the angles of each circle inside the other) integrated along
  the axis; within 1e-30 in exact frames and 1e-15 in turned ones (their
  stored axes not exactly orthonormal);
* `fuse = A + B - common` and `cut = A - common` (each operation sliced
  apart), both inputs' slicings against their closed forms, the area
  identity `area(fuse) + area(common) + 2 opp = area(A) + area(B)`, every
  face's classes summing to its area (a cylindrical wall's by mpmath's
  quadrature of its element), both inputs' areas, a second slicing
  direction `(2, -3, 5)` (the prism then cut obliquely: its circles
  ellipses in the slices; pairs without a cap), within 1e-30;
* Monte-Carlo estimates (200,000 uniform points per pair) of every volume
  and centre within 5 standard errors;
* a scan for near coincidences (a class of a face or a result's volume
  positive but below 1e-9 of the case's size, slicing breakpoints closer
  than that), and every prism cap circle at least 1e-3 (of the squared
  radius) from tangency with the sphere (its critical values of `|X - c|^2
  - R^2`), a cap's rim likewise from each profile circle: none but in the
  declared degenerate pairs.

Pairs run in worker processes (`--workers`).
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
import math
import os
from pathlib import Path
import random
import struct
import zlib

import mpmath as mp

from identity_reference import Boundary, Case, encode_boolean_case
from curve_surface_reference import stored_axes
import spheres_boolean_reference as ref
from generate_curved_boolean_fixtures import FRAMES, stadium
from generate_sphere_boolean_fixtures import box_sphere, frame_name, ideal_axes, sphere_measures

ROOT = Path(__file__).resolve().parents[1]
PREFIX = 'boolean-spheres'
STEP = 'S9d.2'
BOOLEAN_OPERATION = 93
HP = math.pi/2
EXACT = {'XY', 'SIDE', 'DOWN', 'TURN'}
SECOND = (2, -3, 5)


def at(name, origin):
    return tuple(float(c) for c in origin)+FRAMES[name]


def square(x0, y0, x1, y1):
    return Boundary(points=[(x0, y0), (x1, y0), (x1, y1), (x0, y1)])


def circle(cx, cy, r):
    return Boundary(circle=(float(cx), float(cy), float(r)))


def make(spec, op):
    if spec[0] == 'sphere':
        _, r, frame, low, high = spec
        return Case('', 1e-7, op, frame, 0.0, 0.0, [], sphere=(float(r), low, high))
    _, boundaries, frame, start, end = spec
    return Case('', 1e-7, op, frame, float(start), float(end), boundaries)


def sphere(r, frame, low=-HP, high=HP):
    return ('sphere', r, frame, low, high)


def prism(boundaries, frame, start, end):
    return ('prism', boundaries, frame, start, end)


class Boolean:
    def __init__(self, name, operation, obj, tool, kind='solid', reason=None):
        self.name, self.operation, self.kind, self.reason = name, operation, kind, reason
        self.pair_name = name.rsplit('_', 1)[0]
        self.obj = make(obj, 91)
        self.obj.name = name
        self.tool = make(tool, 92)
        self.tool.name = name
        self.frames = (obj[2], tool[2])
        self.specs = (obj, tool)

    def encode(self):
        return encode_boolean_case(self.obj, self.operation, self.tool, BOOLEAN_OPERATION)


def group(name, obj, tool, ops):
    out = []
    for op, kind in ops.items():
        reason = None
        if isinstance(kind, tuple):
            kind, reason = kind
        out.append(Boolean(f'{name}_{op}', op, obj, tool, kind, reason))
    return out


TANGENT_SPHERES = 'two spheres tangent'
TANGENT_CYLINDER = 'a cylinder tangent to the sphere'
HOLED = [square(-2.25, -2.25, 2.25, 2.25), circle(0, 0, 1.5)]


def cases():
    """S9d.2's pairs (see the module's docstring)."""
    O = (0, 0, 0)
    out = []
    out += group('spheres_cross', sphere(2, at('XY', O)), sphere(1.5, at('XY', (2.25, 0.5, 0.25))),
                 {'fuse': 'solid', 'cut': 'solid', 'common': 'solid'})
    out += group('spheres_nested', sphere(2, at('XY', O)), sphere(0.75, at('XY', (0.5, -0.25, 0.5))),
                 {'cut': 'solid'})
    out += group('spheres_apart', sphere(1, at('XY', O)), sphere(1.25, at('XY', (1.5, 2, 0.5))),
                 {'fuse': 'solid', 'common': 'empty'})
    out += group('spheres_turned', sphere(2, at('TILTX', (1, 2, -1))), sphere(1.25, at('LEAN', (2.5, 1.25, 0))),
                 {'common': 'solid'})
    pipe = prism([circle(0, 0, 1.5)], at('XY', O), -3, 3)
    out += group('pipe_ring', sphere(2.5, at('XY', O)), pipe, {'cut': 'solid', 'common': 'solid'})
    out += group('pipe_ends', pipe, sphere(2.5, at('XY', O)), {'cut': 'solid'})
    out += group('pipe_side', sphere(2.5, at('XY', (0.5, 0.25, -0.25))),
                 prism([circle(0.25, -0.25, 2)], at('SIDE', O), -3, 3.5), {'common': 'solid'})
    holed = prism(HOLED, at('XY', O), -3, 3.5)
    out += group('hole_box', sphere(2.5, at('XY', O)), holed, {'cut': 'solid', 'common': 'solid'})
    out += group('holed_box', holed, sphere(2.5, at('XY', O)), {'cut': 'solid'})
    out += group('dome_pipe', sphere(2.5, at('XY', O), 0.0, HP), prism([circle(0, 0, 1.5)], at('XY', O), -1, 4),
                 {'cut': 'solid', 'common': 'solid'})
    rod = prism([circle(0.5, 0.25, 0.75)], at('XY', O), -3, 3)
    out += group('rings', sphere(2, at('XY', O)), rod, {'cut': 'solid', 'common': 'solid'})
    out += group('rod_ends', rod, sphere(2, at('XY', O)), {'cut': 'solid'})
    out += group('rings_tilt', sphere(2, at('XY', O)), prism([circle(0.5, 0.25, 0.75)], at('TILT', O), -3, 3),
                 {'common': 'solid'})
    out += group('bite', sphere(2, at('XY', O)), prism([circle(2, 0.5, 1)], at('XY', O), -3, 3),
                 {'fuse': 'solid', 'cut': 'solid', 'common': 'solid'})
    out += group('bite_tilt', sphere(2, at('XY', O)), prism([circle(2, 0.5, 1)], at('TILT', O), -3, 3),
                 {'cut': 'solid'})
    out += group('cap_cross', sphere(2, at('XY', O)), prism([circle(1.5, 0, 1)], at('XY', O), -3, 1),
                 {'cut': 'solid', 'common': 'solid'})
    stad = prism([stadium(2.0, 1.0)], at('XY', O), -3, 3)
    out += group('stadium', sphere(2, at('XY', (0.25, 0.25, 0))), stad, {'cut': 'solid', 'common': 'solid'})
    out += group('stadium_ends', stad, sphere(2, at('XY', (0.25, 0.25, 0))), {'cut': 'solid'})
    out += group('spheres_tangent', sphere(2, at('XY', O)), sphere(1, at('XY', (3, 0, 0))),
                 {'fuse': ('degenerate', TANGENT_SPHERES), 'common': ('degenerate', TANGENT_SPHERES)})
    out += group('rod_tangent', sphere(2, at('XY', O)), prism([circle(3, 0, 1)], at('XY', O), -3, 3),
                 {'common': ('degenerate', TANGENT_CYLINDER)})
    out += group('rod_inside_tangent', sphere(2, at('XY', O)), prism([circle(1.25, 0, 0.75)], at('XY', O), -3, 3),
                 {'cut': ('degenerate', TANGENT_CYLINDER)})
    return out


def exact_pair(case):
    """Whether the closed forms hold exactly: no prism in a turned frame (a
    whole sphere's model does not depend on its frame's axes)."""
    return all(frame_name(f) in EXACT for f, s in zip(case.frames, case.specs) if s[0] == 'prism')


# ------------------------------------------------------------------ closed forms

Z = mp.mpf(0)


def cap_volume(R, h):
    return mp.pi*h*h*(3*R-h)/3


def cap_centroid(R, h):
    """The distance of a cap's centroid from its sphere's centre."""
    return 3*(2*R-h)**2/(4*(3*R-h))


def lens(r1, r2, d):
    """Two discs (radii `r1`, `r2`, centres `d` apart): the common's area,
    its first moment along the line of centres about the first centre, the
    angle of the first circle inside the second disc and of the second
    inside the first."""
    tau = 2*mp.pi
    if r1 <= 0:
        return Z, Z, Z, Z
    if d >= r1+r2:
        return Z, Z, Z, Z
    if d <= r2-r1:
        return mp.pi*r1*r1, Z, tau, Z
    if d <= r1-r2:
        return mp.pi*r2*r2, mp.pi*r2*r2*d, Z, tau
    x1 = (d*d+r1*r1-r2*r2)/(2*d)
    x2 = d-x1
    seg = lambda r, x: r*r*mp.acos(x/r)-x*mp.sqrt(r*r-x*x)
    mom = lambda r, x: 2*(r*r-x*x)**mp.mpf(1.5)/3
    A1, A2 = seg(r1, x1), seg(r2, x2)
    return A1+A2, mom(r1, x1)+A2*d-mom(r2, x2), 2*mp.acos(x1/r1), 2*mp.acos(x2/r2)


def assemble(A, B, common, A_in, B_in):
    """Every operation by inclusion and exclusion: `A`, `B` (volume,
    moments, area), the common (volume, moments), each input's area inside
    the other."""
    VA, MA, AA = A
    VB, MB, AB = B
    Vc, Mc = common
    A_out, B_out = AA-A_in, AB-B_in
    ops = {'common': (Vc, Mc, A_in+B_in),
           'fuse': (VA+VB-Vc, tuple(a+b-c for a, b, c in zip(MA, MB, Mc)), A_out+B_out),
           'cut': (VA-Vc, tuple(a-c for a, c in zip(MA, Mc)), A_out+B_in)}
    out = {}
    for op, (V, m, area) in ops.items():
        out[op] = None if V <= mp.mpf(10)**-30 else (V, area, tuple(x/V for x in m))
    return out


def spheres_forms(sa, sb):
    """Two whole spheres in world coordinates."""
    ca = tuple(mp.mpf(v) for v in sa[2][:3])
    cb = tuple(mp.mpf(v) for v in sb[2][:3])
    Ra, Rb = mp.mpf(sa[1]), mp.mpf(sb[1])
    d = mp.sqrt(sum((x-y)**2 for x, y in zip(ca, cb)))
    u = tuple((y-x)/d for x, y in zip(ca, cb))
    VA, VB = 4*mp.pi*Ra**3/3, 4*mp.pi*Rb**3/3
    A = (VA, tuple(VA*x for x in ca), 4*mp.pi*Ra*Ra)
    B = (VB, tuple(VB*x for x in cb), 4*mp.pi*Rb*Rb)
    if d >= Ra+Rb:
        return assemble(A, B, (Z, (Z, Z, Z)), Z, Z)
    if d <= Ra-Rb:
        return assemble(A, B, (VB, B[1]), Z, B[2])
    if d <= Rb-Ra:
        return assemble(A, B, (VA, A[1]), A[2], Z)
    x = (d*d+Ra*Ra-Rb*Rb)/(2*d)
    ha, hb = Ra-x, Rb-(d-x)
    Va, Vb = cap_volume(Ra, ha), cap_volume(Rb, hb)
    ga, gb = cap_centroid(Ra, ha), cap_centroid(Rb, hb)
    M = tuple(Va*(ca[i]+ga*u[i])+Vb*(cb[i]-gb*u[i]) for i in range(3))
    return assemble(A, B, (Va+Vb, M), 2*mp.pi*Ra*ha, 2*mp.pi*Rb*hb)


def coax(R, a, z0, z1, zone=None):
    """A ball (radius `R`, `zone` its axial range, from its centre) against
    the coaxial cylinder of radius `a` between `z0` and `z1`: the common's
    volume and axial moment, the sphere's face inside the cylinder, the
    cylinder's wall inside the ball, its caps inside the ball and the zone's
    end discs inside the cylinder."""
    R, a, z0, z1 = (mp.mpf(v) for v in (R, a, z0, z1))
    zlo, zhi = (-R, R) if zone is None else (mp.mpf(zone[0]), mp.mpf(zone[1]))
    h0 = mp.sqrt(R*R-a*a) if a < R else Z
    lo, hi = max(z0, zlo, -R), min(z1, zhi, R)
    V = Mz = Z
    if hi > lo:
        cuts = sorted({lo, hi}|{z for z in (-h0, h0) if lo < z < hi})
        for p, q in zip(cuts, cuts[1:]):
            if abs((p+q)/2) < h0:
                V += mp.pi*a*a*(q-p)
                Mz += mp.pi*a*a*(q*q-p*p)/2
            else:
                V += mp.pi*(R*R*(q-p)-(q**3-p**3)/3)
                Mz += mp.pi*(R*R*(q*q-p*p)/2-(q**4-p**4)/4)

    def length(a_, b_, inner):
        """The length of [a_, b_] where |z| < h0 (inner) or > h0."""
        if b_ <= a_:
            return Z
        mid = max(Z, min(b_, h0)-max(a_, -h0))
        return mid if inner else (b_-a_)-mid
    sphere_in = 2*mp.pi*R*length(max(z0, zlo, -R), min(z1, zhi, R), False)
    wall_in = 2*mp.pi*a*length(max(z0, zlo), min(z1, zhi), True)
    disc = lambda z: mp.pi*min(a*a, R*R-z*z) if abs(z) < R else Z
    caps_in = sum((disc(z) for z in (z0, z1) if zlo < z < zhi), Z)
    ends = [z for z, whole in ((zlo, zone is None or zlo <= -R), (zhi, zone is None or zhi >= R)) if not whole]
    ends_in = sum((disc(z) for z in ends if z0 < z < z1), Z)
    return {'volume': V, 'moment': Mz, 'sphere_in': sphere_in, 'wall_in': wall_in, 'caps_in': caps_in,
            'ends_in': ends_in}


def local_frame(spec):
    """A prism's ideal orthonormal axes and origin."""
    x, y, n = ideal_axes(frame_name(spec[2]))
    return tuple(mp.mpf(v) for v in spec[2][:3]), (x, y, n)


def to_local(o, axes, X):
    return tuple(sum((X[j]-o[j])*axes[i][j] for j in range(3)) for i in range(3))


def to_world(o, axes, V, m):
    return tuple(o[j]*V+sum(m[i]*axes[i][j] for i in range(3)) for j in range(3))


def sphere_input(s, zone_local=None):
    """A ball's volume, moments (world) and area."""
    R = mp.mpf(s[1])
    c = tuple(mp.mpf(v) for v in s[2][:3])
    zone = zone_of(s)
    V, Mz, A = sphere_measures(R, zone)
    if zone is None:
        return V, tuple(V*x for x in c), A
    _, _, n = ideal_axes(frame_name(s[2]))
    return V, tuple(V*c[i]+Mz*n[i] for i in range(3)), A


def zone_of(s):
    _, r, _, low, high = s
    lo = None if low == -HP else r*math.sin(low)
    hi = None if high == HP else r*math.sin(high)
    if lo is None and hi is None:
        return None
    return (-r if lo is None else lo, r if hi is None else hi)


def prism_forms(first):
    """A sphere (or cap) against a prism of one circle (coaxial or not) or a
    square with a coaxial round hole, in the prism's ideal frame."""
    obj, tool = first.specs
    s, p = (obj, tool) if obj[0] == 'sphere' else (tool, obj)
    o, axes = local_frame(p)
    R = mp.mpf(s[1])
    c = tuple(mp.mpf(v) for v in s[2][:3])
    cl = to_local(o, axes, c)
    w0, w1 = mp.mpf(p[3]), mp.mpf(p[4])
    bounds = p[1]
    zone = zone_of(s)
    if zone is not None:
        sn = ideal_axes(frame_name(s[2]))[2]
        assert abs(ref.dot(sn, axes[2])-1) < mp.mpf(10)**-30, 'a cap along the prism axis'
    S = sphere_input(s)
    if len(bounds) == 1 and bounds[0].circle is not None:
        cx, cy, a = (mp.mpf(v) for v in bounds[0].circle)
        e = mp.sqrt((cx-cl[0])**2+(cy-cl[1])**2)
        VP = mp.pi*a*a*(w1-w0)
        P = (VP, to_world(o, axes, VP, (cx*VP, cy*VP, (w0+w1)/2*VP)), 2*mp.pi*a*(w1-w0)+2*mp.pi*a*a)
        if e == 0:
            f = coax(R, a, w0-cl[2], w1-cl[2], zone)
            Vc = f['volume']
            Mc = to_world(o, axes, Vc, (cl[0]*Vc, cl[1]*Vc, cl[2]*Vc+f['moment']))
            s_in, p_in = f['sphere_in']+f['ends_in'], f['wall_in']+f['caps_in']
        else:
            assert zone is None, 'an off-axis cylinder against a whole sphere'
            u = ((cx-cl[0])/e, (cy-cl[1])/e)
            rho = lambda w: mp.sqrt(max(R*R-(w-cl[2])**2, Z))
            lo, hi = max(w0, cl[2]-R), min(w1, cl[2]+R)
            kinks = {lo, hi}
            for k in (e+a, abs(e-a)):
                if k <= R:
                    dz = mp.sqrt(R*R-k*k)
                    kinks.update(z for z in (cl[2]-dz, cl[2]+dz) if lo < z < hi)
            pts = sorted(kinks)
            L = lambda w: lens(rho(w), a, e)
            Vc = mp.quad(lambda w: L(w)[0], pts)
            mu = mp.quad(lambda w: L(w)[1], pts)
            Mw = mp.quad(lambda w: w*L(w)[0], pts)
            Mc = to_world(o, axes, Vc, (cl[0]*Vc+u[0]*mu, cl[1]*Vc+u[1]*mu, Mw))
            s_in = R*mp.quad(lambda w: L(w)[2], pts)
            p_in = a*mp.quad(lambda w: L(w)[3], pts)
            p_in += sum((lens(rho(w), a, e)[0] for w in (w0, w1) if abs(w-cl[2]) < R), Z)
    else:
        # A square with a coaxial round hole.
        (sq, hole) = bounds
        pts = sq.points
        (u0, v0), (u1, v1) = pts[0], pts[2]
        hx, hy, a = (mp.mpf(v) for v in hole.circle)
        assert hx == cl[0] and hy == cl[1] and zone is None
        box = [(u0, u1), (v0, v1), (w0, w1)]
        f = box_sphere(box, cl, R)
        core = coax(R, a, w0-cl[2], w1-cl[2])
        Vc = f['volume']-core['volume']
        m = tuple(f['moments'][i]-cl[i]*core['volume']-(core['moment'] if i == 2 else 0) for i in range(3))
        Mc = to_world(o, axes, Vc, m)
        s_in = f['sphere_in']-core['sphere_in']
        p_in = f['box_in']-core['caps_in']+core['wall_in']
        side_u, side_v = mp.mpf(u1-u0), mp.mpf(v1-v0)
        Ap = side_u*side_v-mp.pi*a*a
        h = w1-w0
        VP = Ap*h
        mu = mp.mpf(u0+u1)/2*side_u*side_v-hx*mp.pi*a*a
        mv = mp.mpf(v0+v1)/2*side_u*side_v-hy*mp.pi*a*a
        P = (VP, to_world(o, axes, VP, (mu*h, mv*h, (w0+w1)/2*VP)), 2*Ap+2*(side_u+side_v)*h+2*mp.pi*a*h)
    if obj[0] == 'sphere':
        return assemble(S, P, (Vc, Mc), s_in, p_in)
    return assemble(P, S, (Vc, Mc), p_in, s_in)


def closed_forms(first):
    """{op: (volume, area, centre) or None}, or None where no closed form is
    taken (the stadium)."""
    obj, tool = first.specs
    if obj[0] == 'sphere' and tool[0] == 'sphere':
        return spheres_forms(obj, tool)
    p = obj if obj[0] == 'prism' else tool
    if any(b.segments is not None for b in p[1]):
        return None
    return prism_forms(first)


# ------------------------------------------------------------------ Monte Carlo

class FloatBall:
    def __init__(self, S):
        self.c = [float(v) for v in S.c]
        self.r2 = float(S.r2)
        self.planes = [([float(v) for v in a], float(b)) for a, b, _ in S.planes]
        r = float(S.r)
        self.box = ([v-r for v in self.c], [v+r for v in self.c])

    def contains(self, X):
        if sum((X[i]-self.c[i])**2 for i in range(3)) >= self.r2:
            return False
        return all(a[0]*X[0]+a[1]*X[1]+a[2]*X[2] > b for a, b in self.planes)


class FloatArcPrism:
    def __init__(self, P):
        self.o = [float(v) for v in P.o]
        self.inv = [[float(c) for c in row] for row in P.Kinv]
        self.lo, self.hi = float(P.lo), float(P.hi)
        self.els = []
        us, vs = [], []
        for e in P.elements:
            if e[0] == 'seg':
                self.els.append(('seg', float(e[1][0]), float(e[1][1]), float(e[2][0]-e[1][0]), float(e[2][1]-e[1][1])))
                us += [float(e[1][0])]
                vs += [float(e[1][1])]
            else:
                _, C, r, t0, t1, full, _, _ = e
                self.els.append(('arc', float(C[0]), float(C[1]), float(r), float(t0), float(t1), full))
                us += [float(C[0]-r), float(C[0]+r)]
                vs += [float(C[1]-r), float(C[1]+r)]
        pts = [[float(c) for c in P.world((F_(u), F_(v)), w)] for u in (min(us), max(us)) for v in (min(vs), max(vs))
               for w in (P.lo, P.hi)]
        self.box = ([min(p[i] for p in pts) for i in range(3)], [max(p[i] for p in pts) for i in range(3)])

    def contains(self, X):
        d = (X[0]-self.o[0], X[1]-self.o[1], X[2]-self.o[2])
        w = sum(self.inv[2][i]*d[i] for i in range(3))
        if not self.lo < w < self.hi:
            return False
        u = sum(self.inv[0][i]*d[i] for i in range(3))
        v = sum(self.inv[1][i]*d[i] for i in range(3))
        ru, rv = 0.8150192, 0.5794340
        count = 0
        for el in self.els:
            if el[0] == 'seg':
                _, px, py, ex, ey = el
                den = ru*ey-rv*ex
                if den == 0:
                    continue
                wx, wy = px-u, py-v
                t = (wx*ey-wy*ex)/den
                s = (wx*rv-wy*ru)/den
                if t > 0 and 0 <= s < 1:
                    count += 1
            else:
                _, cx, cy, r, t0, t1, full = el
                Wx, Wy = u-cx, v-cy
                b = ru*Wx+rv*Wy
                disc_ = b*b-(Wx*Wx+Wy*Wy-r*r)
                if disc_ <= 0:
                    continue
                sq = math.sqrt(disc_)
                for t in (-b-sq, -b+sq):
                    if t > 0:
                        if full:
                            count += 1
                            continue
                        phi = math.atan2(v+t*rv-cy, u+t*ru-cx)
                        span = abs(t1-t0)
                        rel = (phi-t0) % (2*math.pi) if t1 > t0 else (t0-phi) % (2*math.pi)
                        if rel <= span:
                            count += 1
        return count % 2 == 1


def F_(x):
    from fractions import Fraction
    return Fraction(x)


def monte_carlo(pair, n, seed):
    """Volume and centre estimates of the three operations and their
    standard errors, from `n` uniform points in the inputs' box."""
    fa = FloatBall(pair.A.S) if isinstance(pair.A, ref.Ball) else FloatArcPrism(pair.A)
    fb = FloatBall(pair.B.S) if isinstance(pair.B, ref.Ball) else FloatArcPrism(pair.B)
    box0 = [min(fa.box[0][i], fb.box[0][i]) for i in range(3)]
    box1 = [max(fa.box[1][i], fb.box[1][i]) for i in range(3)]
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


# ------------------------------------------------------------------ clearances

def cap_circle_clearance(pair):
    """The smallest critical value of `|X - c|^2 - R^2` over every prism
    cap circle's arc (relative to `R^2`), and of `|q - C|^2 - r^2` over a
    cap's rim against each profile circle (relative to `r^2`)."""
    if pair.two_balls:
        return None
    ball, P = pair.D, pair.P
    best = None

    def note(v):
        nonlocal best
        best = v if best is None else min(best, v)
    for e in P.elements:
        if e[0] != 'arc':
            continue
        _, C, r, t0, t1, full, _, _ = e
        for w in (P.lo, P.hi):
            g = ref.cap_circle_trig(P, C, r, w, ball)
            dg = g.deriv()
            if dg.is_zero():
                note(abs(g.value(mp.mpf(0)))/ball.rm**2)
                continue
            for th in dg.angles():
                if full or ref.span_param(t0, t1, th) is not None:
                    note(abs(g.value(th))/ball.rm**2)
        # A cap's rims (slices) against this circle, sampled and refined.
        for (k, b) in ball.sbounds:
            s = M(b/k)
            rim = ball.section(s, ignore_zone=True).base[0]
            uv = lambda t: P.uv_of(s, rim.point(t))
            h = lambda t: ((uv(t)[0]-M(C[0]))**2+(uv(t)[1]-M(C[1]))**2-M(r)**2)/M(r)**2
            n = 720
            vals = [h(2*mp.pi*j/n) for j in range(n)]
            for j in range(n):
                a_, b_, c_ = vals[j-1], vals[j], vals[(j+1) % n]
                if (b_-a_)*(c_-b_) <= 0:
                    lo_, hi_ = 2*mp.pi*(j-1)/n, 2*mp.pi*(j+1)/n
                    for _ in range(60):
                        m1, m2 = lo_+(hi_-lo_)/3, hi_-(hi_-lo_)/3
                        if abs(h(m1)) < abs(h(m2)):
                            hi_ = m2
                        else:
                            lo_ = m1
                    note(abs(h((lo_+hi_)/2)))
    return best


def M(x):
    return ref.M(x)


# ------------------------------------------------------------------ the pairs' work

def evaluate(job):
    """One pair: the reference's rows for its operations and every check's
    deviations (run in a worker process)."""
    name, first, ops, mc_n = job
    obj, tool = first.obj, first.tool
    pair = ref.Pair(obj, tool)
    res = {op: pair.result(op) for op in ref.OPS}
    rows = {op: ref.rows(obj, op, tool, pair)[0] for op in ops}
    size = pair.size
    r = pair.sliced()
    closed = {}
    for x in (pair.D, pair.P):
        closed[x] = x.closed()
    VD, MD = closed[pair.D][:2]
    VP, MP = closed[pair.P][:2]
    checks = {}
    dev = mp.mpf(0)
    for (v, m), (V, Mm) in ((r['D'], (VD, MD)), (r['P'], (VP, MP))):
        dev = max(dev, abs(v-V)/size**3, max(abs(x-y) for x, y in zip(m, Mm))/size**4)
    checks['inputs_sliced'] = dev
    vc, mc = r['common']
    dev = mp.mpf(0)
    for key, (V, Mm) in (('fuse', (VD+VP-vc, [a+b-c for a, b, c in zip(MD, MP, mc)])),
                         ('D-P', (VD-vc, [a-c for a, c in zip(MD, mc)])),
                         ('P-D', (VP-vc, [a-c for a, c in zip(MP, mc)]))):
        v, m = r[key]
        dev = max(dev, abs(v-V)/size**3, max(abs(x-y) for x, y in zip(m, Mm))/size**4)
    checks['inclusion_exclusion'] = dev
    faces = pair.faces()
    checks['face_classes'] = max(abs(sum(cls.values())-area) for _, _, cls, area in faces)/size**2
    ca, cb = pair.classes('A'), pair.classes('B')
    checks['shared_both_sides'] = max(abs(ca['same']-cb['same']), abs(ca['opp']-cb['opp']))/size**2
    area_inputs = {}
    for x in (pair.A, pair.B):
        area_inputs[x.role] = input_area(x)
    checks['area_identity'] = abs(pair.area('fuse')+pair.area('common')+2*ca['opp']
                                  - area_inputs['A']-area_inputs['B'])/size**2
    checks['input_areas'] = max(abs(sum(sum(cls.values()) for w, _, cls, _ in faces if w == k)-A)
                                for k, A in area_inputs.items())/size**2
    if not any(isinstance(x, ref.Ball) and x.planes for x in (pair.A, pair.B)):
        second = ref.Pair(obj, tool, d=SECOND).sliced()
        dev = mp.mpf(0)
        for key in ('D', 'P', 'common', 'fuse', 'D-P', 'P-D'):
            dev = max(dev, abs(second[key][0]-r[key][0])/size**3,
                      max(abs(x-y) for x, y in zip(second[key][1], r[key][1]))/size**4)
        for key in r:
            if isinstance(key, tuple):
                for k in ('in', 'out'):
                    dev = max(dev, abs(second[key][k]-r[key][k])/size**2)
        checks['second_direction'] = dev
    mc = monte_carlo(pair, mc_n, zlib.crc32(name.encode()))
    z = 0.0
    for op in ref.OPS:
        vol, sv, c, sc = mc[op]
        V, m = pair.volume(op)
        v = float(V)
        z = max(z, abs(vol-v)/sv)
        if v > 1e-9 and c is not None:
            cent = [float(x/V) for x in m]
            z = max(z, max(abs(c[i]-cent[i])/max(sc[i], 1e-12) for i in range(3)))
    checks['monte_carlo_sigma'] = mp.mpf(z)
    near = near_coincidences(pair)
    clearance = cap_circle_clearance(pair)
    if clearance is not None and clearance < mp.mpf(10)**-3:
        near.append(f'a cap circle within {mp.nstr(clearance, 3)} of tangency')
    stats = {'volume_quadrature': pair.quad_error/size**4, 'breaks': len(pair.breaks), 'clearance': clearance}
    return name, rows, res, checks, near, stats


def input_area(x):
    """An input's area in closed form: a sphere's by `Sphere.closed`, a
    prism's caps by its profile's area and its walls by their own elements
    (a cylindrical one by mpmath's quadrature)."""
    if isinstance(x, ref.Ball):
        return x.closed()[2]
    A = x.profile_moments()[0]
    xy = ref.cross(x.x, x.y)
    total = 2*A*mp.sqrt(M(ref.dot(xy, xy)))
    h = M(x.hi-x.lo)
    xm, ym, nm = ref.Mv(x.x), ref.Mv(x.y), ref.Mv(x.n)
    for e in x.elements:
        if e[0] == 'seg':
            p, q = e[1], e[2]
            ew = ref.add(ref.scale(x.x, q[0]-p[0]), ref.scale(x.y, q[1]-p[1]))
            c = ref.cross(ew, x.n)
            total += mp.sqrt(M(ref.dot(c, c)))*h
        else:
            _, C, r, t0, t1, _, _, _ = e

            def el(th):
                t = tuple(-mp.sin(th)*xm[i]+mp.cos(th)*ym[i] for i in range(3))
                c = ref.cross(t, nm)
                return mp.sqrt(ref.dot(c, c))
            total += M(r)*h*mp.quad(el, [min(t0, t1), max(t0, t1)])
    return total


def near_coincidences(pair):
    """A face class, or a result's volume, positive but thinner than 1e-9 of
    the case's size; slicing breakpoints closer than that (and farther than
    1e-30)."""
    size = pair.size
    out = []
    for w, tag, cls, _ in pair.faces():
        for c, v in cls.items():
            if 0 < v < mp.mpf(10)**-9*size**2:
                out.append(f'{w} {tag} {c} {mp.nstr(v, 3)}')
    for op in ref.OPS:
        v = pair.volume(op)[0]
        if mp.mpf(10)**-25*size**3 < v < mp.mpf(10)**-9*size**3:
            out.append(f'{op} volume {mp.nstr(v, 3)}')
    b = pair.raw_breaks
    for p, q in zip(b, b[1:]):
        if mp.mpf(10)**-30*size < q-p < mp.mpf(10)**-9*size:
            out.append(f'breakpoints {mp.nstr(p, 12)} and {mp.nstr(q, 12)}')
    return out


def run(jobs, fn, workers):
    if workers <= 1:
        return [fn(j) for j in jobs]
    with ProcessPoolExecutor(workers) as pool:
        return list(pool.map(fn, jobs))


def closed_job(first):
    return first.pair_name, closed_forms(first)


LIMITS = {'closed_forms_exact_frames': 1e-30, 'closed_forms_turned_frames': 1e-15,
          'inputs_sliced': 1e-30, 'inclusion_exclusion': 1e-30, 'face_classes': 1e-30,
          'shared_both_sides': 1e-30, 'area_identity': 1e-30, 'input_areas': 1e-30,
          'second_direction': 1e-30, 'monte_carlo_sigma': 5, 'volume_quadrature_estimate': 1e-30}


def reference_checks(results, forms):
    """Largest deviations of each check, relative to the case's size."""
    worst, covered = {}, {}

    def note(key, value):
        worst[key] = max(worst.get(key, mp.mpf(0)), value)
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
            V, A, C = want
            size = max(abs(V)**(mp.mpf(1)/3), 1)
            dev = max(abs(vol-V)/abs(V), abs(area-A)/abs(A), max(abs(centre[i]-C[i]) for i in range(3))/size)
            note(f'closed_forms_{kind}_frames', dev)
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
        assert c.obj.sphere is not None or c.tool.sphere is not None, f'{c.name}: a sphere'
        for s in (c.obj, c.tool):
            if s.sphere is not None and (s.sphere[1] != -HP or s.sphere[2] != HP):
                assert frame_name(s.frame) in EXACT and s.frame[:3] == (0.0, 0.0, 0.0), \
                    f'{c.name}: a cap or zone in an exact frame at the origin'
        for f in c.frames:
            frame_name(f)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--workers', type=int, default=max(1, (os.cpu_count() or 2)-2))
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
          '-', ', '.join(f'{v} {k}' for k, v in sorted(kinds.items())),
          '-', f'{sum(1 for c in listed if exact_pair(c))} in exact frames (whole spheres in any)')


if __name__ == '__main__':
    main()
