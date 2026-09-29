#!/usr/bin/env python3
"""Fixtures for S9d.4a of REVIEW_NOTES.md: Booleans of a whole torus
(`Solid::torus_with`, the full tube and turn) against a polyhedral prism in
any relative position, before any of its kernel code.

`boolean-torus-cases.txt` lists each case in the Boolean protocol
(`identity_reference.encode_boolean_case`: a torus's block is its `frame`
and `torus MAJOR MINOR 0 2pi 2pi` rows); `boolean-torus-expected.tsv` gives
per case `expect KIND S9d.4a` (`solid`, `empty` or `degenerate`, then
`reason TEXT`) and `result N volume area cx cy cz` or `empty` from
`torus_boolean_reference.py`, as the cone fixtures;
`boolean-torus-frames.tsv` the stored axes' bits. Frames are the curved
generator's (`FRAMES`, their stored axes the kernel's bit for bit).

Pairs, with `T` the torus of radii 5/2 and 3/2 (its hole of radius 1, its
tube from 1 to 4 about the axis and from -3/2 to 3/2 along it) on `XY` at
the origin unless said: a bar through the hole across the whole torus, its
walls cutting the tube in loops about the tube (`bar`), and a thin strip
along the equator cut in three (`strip`, the strip first); a slab normal to
the axis cutting loops about the axis (`band`, and `slab` the slab first:
the hole's disc and the outside two solids), a half-space above the
equator (`up`); a cap cut from the tube's outside by a wall (`cap_out`), a
box through the hole whose four walls cut such caps (`frame`), a box in the
hole whose vertical edges cut caps from the tube's inside (`corner_caps`,
contractible loops over two walls) and a wedge whose one edge does
(`wedge_cap`); a half-space through the axis (`half`, loops about the
tube); a box inside the tube (`in_tube`), the torus inside a box
(`in_box`), a box in the hole (`hole`); turned frames: `T` in `TILT` above
a plane inclined as a Villarceau plane but 3/4 from the centre (`tilt_up`,
one contractible loop), `T` in `TILTX` across a slab between its saddle
levels (`tiltx_band`, two loops about the tube in each plane), a bar in
`LEAN` across `T` (`lean_bar`); and `degenerate` with reasons: a face
tangent to the torus along its top circle (`top_tangent`), a wall tangent
to its inner equator (`inner_tangent`, a figure eight), a face on a
Villarceau plane (`villarceau`, tangent at two points) and a box's vertex
on the torus (`vertex_on`).

Before writing, `reference_checks` compares the reference with independent
results (limits relative to the case's size):

* closed forms of every pair but the wedge, in the torus's ideal
  (orthonormal) frame, by one quadrature (mpmath's tanh-sinh between the
  kinks) of exact 2D forms along the axis, each slice's annulus the outer
  disc less the inner: an aligned box by S9d.1's rectangle-in-disc
  antiderivatives (the walls inside by their chords `|[v0, v1] n [-Q, Q]|`,
  `Q^2 = rho^2 - u^2`, the wall's element `r rho / q` times the circle's
  angle inside the rectangle), a half-space by circular segments of both
  circles (normal to the axis in closed form: `4 pi R q` per slice, the
  wall `4 pi R r` times the angle of latitude), a slab as the torus less
  two half-spaces' complements; so the common's volume, moments and both
  inputs' surfaces inside the other, and every operation by inclusion and
  exclusion; within 1e-30 in exact frames and 1e-15 in turned ones (their
  stored axes not exactly orthonormal), the declared degenerate pairs within
  1e-25 (their kinks within rounding of each other: a vertex 2e-17 off the
  torus leaves features of its square root's width);
* `fuse = A + B - common` and `cut = A - common` (each operation sliced
  apart, `A` and `B` in closed form), the common a second way (each piece
  clipped by both discs), both inputs' slicings against their closed forms
  (the torus's `2 pi^2 R r^2` and `4 pi^2 R r`), the area identity
  `area(fuse) + area(common) = area(A) + area(B)`, every face's classes
  summing to its area, a second slicing direction (the meridian
  half-planes about the axis, every section the tube's disc, the volume by
  the cylindrical element: every operation's volume and moments), the wall
  two ways (in its own `(theta, phi)` by the meridians, and by the normal
  slices' circles with the element `r rho / q`), within 1e-30;
* Monte-Carlo estimates (200,000 uniform points per pair) of every volume
  and centre within 5 standard errors;
* a scan for near coincidences (a class of a face or a result's volume
  positive but below 1e-9 of the case's size, slicing breakpoints of either
  direction closer than that) must find none but in the declared
  degenerate pairs; and every non-degenerate pair's margins (a vertex's
  distance from the torus, a face's plane's distance from the tangent
  planes of its normal whose points of contact lie on the face, the sine
  of the angle at which an edge crosses the torus) are reported.

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
import torus_boolean_reference as ref
from generate_curved_boolean_fixtures import FRAMES
from generate_sphere_boolean_fixtures import FloatPrism, arc_in_rect, frame_name, ideal_axes, rect_disc
from generate_cone_boolean_fixtures import prism_measures, segment

ROOT = Path(__file__).resolve().parents[1]
PREFIX = 'boolean-torus'
STEP = 'S9d.4a'
BOOLEAN_OPERATION = 93
EXACT = {'XY', 'SIDE', 'DOWN', 'TURN'}
TWO_PI = ref.TWO_PI


def at(name, origin):
    return tuple(float(c) for c in origin)+FRAMES[name]


def square(x0, y0, x1, y1):
    return Boundary(points=[(float(x0), float(y0)), (float(x1), float(y0)), (float(x1), float(y1)),
                            (float(x0), float(y1))])


def polygon(*pts):
    return Boundary(points=[(float(u), float(v)) for u, v in pts])


def make(spec, op):
    if spec[0] == 'torus':
        _, R, r, frame = spec
        return Case('', 1e-7, op, frame, 0.0, 0.0, [], torus=(float(R), float(r), 0.0, TWO_PI, TWO_PI))
    _, boundaries, frame, start, end = spec
    return Case('', 1e-7, op, frame, float(start), float(end), boundaries)


def torus(R, r, frame):
    return ('torus', R, r, frame)


def prism(boundaries, frame, start, end):
    return ('prism', boundaries, frame, start, end)


def spec_frame(spec):
    return spec[3] if spec[0] == 'torus' else spec[2]


class Boolean:
    def __init__(self, name, operation, obj, tool, form, kind='solid', reason=None):
        self.name, self.operation, self.kind, self.reason = name, operation, kind, reason
        self.form = form
        self.pair_name = name.rsplit('_', 1)[0]
        self.obj = make(obj, 91)
        self.obj.name = name
        self.tool = make(tool, 92)
        self.tool.name = name
        self.frames = (spec_frame(obj), spec_frame(tool))

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


TOP_TANGENT = 'a face of the prism tangent to the torus along a circle'
POINT_TANGENT = 'a face of the prism tangent to the torus at a point'
VILLARCEAU = 'a face of the prism tangent to the torus at two points (a Villarceau plane)'
VERTEX_ON = 'a vertex of the prism on the torus'


def cases():
    """S9d.4a's pairs (see the module's docstring)."""
    O = (0, 0, 0)
    T_ = torus(2.5, 1.5, at('XY', O))
    BIG = square(-5, -5, 5, 5)
    both = {'cut': 'solid', 'common': 'solid'}
    out = []
    out += group('bar', T_, prism([square(-0.5, -5, 0.5, 5)], at('XY', O), -2, 2),
                 {'fuse': 'solid', 'cut': 'solid', 'common': 'solid'}, ('box',))
    out += group('strip', prism([square(-5, -0.5, 5, 0.5)], at('XY', O), -0.5, 0.5), T_, {'cut': 'solid'}, ('box',))
    out += group('band', T_, prism([BIG], at('XY', O), -0.5, 0.5), both, ('box',))
    out += group('slab', prism([BIG], at('XY', O), -0.5, 0.5), T_, {'cut': 'solid'}, ('box',))
    out += group('up', T_, prism([BIG], at('XY', O), 0.5, 3), both, ('box',))
    out += group('cap_out', T_, prism([square(3.25, -5, 6, 5)], at('XY', O), -2, 2), both, ('box',))
    out += group('frame', T_, prism([square(-3, -3, 3, 3)], at('XY', O), -2, 2), both, ('box',))
    out += group('corner_caps', T_, prism([square(-0.875, -0.875, 0.875, 0.875)], at('XY', O), -2, 2), both,
                 ('box',))
    out += group('wedge_cap', T_, prism([polygon((1.5, 0), (-0.5, 0.75), (-0.5, -0.75))], at('XY', O), -2, 2),
                 {'common': 'solid'}, None)
    out += group('half', T_, prism([square(0, -5, 5, 5)], at('XY', O), -2, 2),
                 {'fuse': 'solid', 'cut': 'solid', 'common': 'solid'}, ('box',))
    out += group('in_tube', T_, prism([square(2.25, -0.25, 2.75, 0.25)], at('XY', O), -0.25, 0.25), both,
                 ('box',))
    out += group('in_box', prism([BIG], at('XY', O), -2, 2), T_, both, ('box',))
    out += group('hole', T_, prism([square(-0.625, -0.625, 0.625, 0.625)], at('XY', O), -2, 2),
                 {'fuse': 'solid', 'common': 'empty'}, ('box',))
    out += group('tilt_up', torus(2.5, 1.5, at('TILT', O)), prism([square(-6, -6, 6, 6)], at('XY', O), 0.75, 5),
                 both, ('planes', [('cap', 0)]))
    out += group('tiltx_band', torus(2.5, 1.5, at('TILTX', (0.25, 0.5, 0))),
                 prism([square(-6, -6, 6, 6)], at('XY', O), -0.25, 0.25), both, ('planes', [('cap', 0), ('cap', 1)]))
    out += group('lean_bar', prism([square(-0.5, -6, 0.5, 6)], at('LEAN', O), -6, 6), T_, both,
                 ('planes', [('wall', 1), ('wall', 3)]))
    out += group('top_tangent', T_, prism([BIG], at('XY', O), 1.5, 3), {'fuse': ('degenerate', TOP_TANGENT)},
                 ('box',))
    out += group('inner_tangent', T_, prism([square(1, -5, 5, 5)], at('XY', O), -2, 2),
                 {'common': ('degenerate', POINT_TANGENT)}, ('box',))
    out += group('villarceau', T_, prism([polygon((-8, -10), (8, -10), (8, 6), (-8, -6))], at('SIDE', O), -6, 6),
                 {'common': ('degenerate', VILLARCEAU)}, ('planes', [('wall', 2)]))
    out += group('vertex_on', T_, prism([square(0, 0, 3.4, 2)], at('XY', O), 0, 1.2),
                 {'cut': ('degenerate', VERTEX_ON)}, ('box',))
    return out


def exact_pair(case):
    return all(frame_name(f) in EXACT for f in case.frames)


# ------------------------------------------------------------------ closed forms

def torus_measures(R, r):
    """Volume and area of a torus."""
    R, r = mp.mpf(R), mp.mpf(r)
    return 2*mp.pi**2*R*r*r, 4*mp.pi**2*R*r


def halfspace_form(R, r, a, b):
    """The torus (ideal frame) and the half-space `a . p >= b`: the common's
    volume and moments, the torus's surface inside, the plane's face
    inside."""
    R, r = mp.mpf(R), mp.mpf(r)
    au, av, aw = (mp.mpf(c) for c in a)
    b = mp.mpf(b)
    nuv = mp.sqrt(au*au+av*av)
    q = lambda s: mp.sqrt(max(r*r-s*s, mp.mpf(0)))
    if nuv < mp.mpf(10)**-30:
        c = b/aw
        lo, hi = (max(c, -r), r) if aw > 0 else (-r, min(c, r))
        zero = {'volume': mp.mpf(0), 'moments': (mp.mpf(0),)*3, 'torus_in': mp.mpf(0), 'face_in': mp.mpf(0)}
        if hi <= lo:
            return zero
        Q1 = lambda s: (s*q(s)+r*r*mp.asin(s/r))/2
        Q2 = lambda s: -q(s)**3/3
        V = 4*mp.pi*R*(Q1(hi)-Q1(lo))
        Mw = 4*mp.pi*R*(Q2(hi)-Q2(lo))
        wall = 4*mp.pi*R*r*(mp.asin(hi/r)-mp.asin(lo/r))
        face = 4*mp.pi*R*q(c) if abs(c) < r else mp.mpf(0)
        return {'volume': V, 'moments': (mp.mpf(0), mp.mpf(0), Mw), 'torus_in': wall, 'face_in': face}
    mu, mv = au/nuv, av/nuv
    d = lambda s: (b-aw*s)/nuv
    kinks = [-r, r]
    # |d| = R +- q: 4 R^2 N D^2 = (D^2 + N (R^2 - r^2 + s^2))^2, D = b - aw s.
    N = nuv*nuv
    D2 = [aw*aw, -2*b*aw, b*b]
    E = [D2[0]+N, D2[1], D2[2]+N*(R*R-r*r)]
    poly = [E[0]*E[0], 2*E[0]*E[1], E[1]*E[1]+2*E[0]*E[2], 2*E[1]*E[2], E[2]*E[2]]
    poly = [poly[0], poly[1], poly[2]-4*R*R*N*D2[0], poly[3]-4*R*R*N*D2[1], poly[4]-4*R*R*N*D2[2]]
    for s in mp.polyroots(poly, maxsteps=400, extraprec=400):
        s = mp.mpc(s)
        if abs(s.imag) < mp.mpf(10)**-15 and -r < s.real < r:
            kinks.append(s.real)
    kinks = sorted(kinks)

    def parts(s):
        qs = q(s)
        out = []
        for rho in (R+qs, R-qs):
            out.append(segment(rho, d(s)))
        return out, qs

    def f(k):
        def g(s):
            (so, si), qs = parts(s)
            if k == 'area':
                return so[0]-si[0]
            if k == 'mom':
                return so[1]-si[1]
            return so[3]-si[3]
        return g

    def wall_phi(phi):
        # The wall by its latitude, s = r sin(phi): the element r rho /
        # q ds is r rho dphi.
        (so, si), qs = parts(r*mp.sin(phi))
        return r*((R+qs)*so[2]+(R-qs)*si[2])
    V = mp.quad(f('area'), kinks)
    Mn = mp.quad(f('mom'), kinks)
    Mw = mp.quad(lambda s: s*f('area')(s), kinks)
    wall = mp.quad(wall_phi, [mp.asin(max(min(s/r, mp.mpf(1)), mp.mpf(-1))) for s in kinks])
    face = mp.sqrt(au*au+av*av+aw*aw)/nuv*mp.quad(f('chord'), kinks)
    return {'volume': V, 'moments': (mu*Mn, mv*Mn, Mw), 'torus_in': wall, 'face_in': face}


def box_form(R, r, box):
    """The torus (ideal frame) and an aligned box `[u0, u1] x [v0, v1] x [w0,
    w1]`: the common's volume and moments, the torus's surface inside and
    the box's faces inside."""
    R, r = mp.mpf(R), mp.mpf(r)
    (u0, u1), (v0, v1), (w0, w1) = [tuple(mp.mpf(x) for x in rng) for rng in box]
    q = lambda s: mp.sqrt(max(r*r-s*s, mp.mpf(0)))
    wa, wb = max(w0, -r), min(w1, r)
    out = {'volume': mp.mpf(0), 'moments': (mp.mpf(0),)*3, 'torus_in': mp.mpf(0), 'box_in': mp.mpf(0)}
    if wb <= wa:
        return out
    kinks = set([wa, wb])
    vals = [abs(u0), abs(u1), abs(v0), abs(v1)]+[mp.sqrt(x*x+y*y) for x in (u0, u1) for y in (v0, v1)]
    for val in vals:
        dq = abs(val-R)
        if dq <= r:
            for s in (-mp.sqrt(r*r-dq*dq), mp.sqrt(r*r-dq*dq)):
                if wa < s < wb:
                    kinks.add(s)
    pts = sorted(kinks)

    def rd(s):
        qs = q(s)
        o = rect_disc(u0, u1, v0, v1, R+qs)
        i = rect_disc(u0, u1, v0, v1, R-qs)
        return tuple(x-y for x, y in zip(o, i))
    V = mp.quad(lambda s: rd(s)[0], pts)
    Mu = mp.quad(lambda s: rd(s)[1], pts)
    Mv_ = mp.quad(lambda s: rd(s)[2], pts)
    Mw = mp.quad(lambda s: s*rd(s)[0], pts)
    out['volume'], out['moments'] = V, (Mu, Mv_, Mw)

    def wall(phi):
        # By the latitude, s = r sin(phi): the element r rho / q ds is r rho
        # dphi.
        qs = q(r*mp.sin(phi))
        return r*((R+qs)*arc_in_rect(u0, u1, v0, v1, R+qs)+(R-qs)*arc_in_rect(u0, u1, v0, v1, R-qs))
    out['torus_in'] = mp.quad(wall, [mp.asin(max(min(s/r, mp.mpf(1)), mp.mpf(-1))) for s in pts])

    def chord(c, lo, hi, s):
        total = mp.mpf(0)
        qs = q(s)
        for rho, sign in ((R+qs, 1), (R-qs, -1)):
            Q2 = rho*rho-c*c
            if Q2 <= 0:
                continue
            Q = mp.sqrt(Q2)
            total += sign*max(mp.mpf(0), min(hi, Q)-max(lo, -Q))
        return total
    box_in = mp.mpf(0)
    for u in (u0, u1):
        box_in += mp.quad(lambda s: chord(u, v0, v1, s), pts)
    for v in (v0, v1):
        box_in += mp.quad(lambda s: chord(v, u0, u1, s), pts)
    for w in (w0, w1):
        if -r < w < r:
            box_in += rd(w)[0]
    out['box_in'] = box_in
    return out


def pair_forms(case):
    """{op: (volume, area, centre)} in closed form for a pair (see the
    module's docstring), in the torus's ideal frame; None for a pair
    without one."""
    if case.form is None:
        return None
    K_role = 'A' if case.obj.torus is not None else 'B'
    kc, pc = (case.obj, case.tool) if K_role == 'A' else (case.tool, case.obj)
    X, Y, Nn = ideal_axes(frame_name(kc.frame))
    ko = tuple(mp.mpf(v) for v in kc.frame[:3])
    R, r = kc.torus[0], kc.torus[1]
    local = lambda P: tuple(ref.dot(ref.sub(P, ko), ax) for ax in (X, Y, Nn))
    world_m = lambda V, m: tuple(ko[j]*V+sum(m[i]*(X, Y, Nn)[i][j] for i in range(3)) for j in range(3))
    VK, AK = torus_measures(R, r)
    MK = world_m(VK, (0, 0, 0))
    VP, MP, AP = prism_measures(pc)
    px, py, pn = ideal_axes(frame_name(pc.frame))
    po = tuple(mp.mpf(v) for v in pc.frame[:3])
    pworld = lambda u, v, w: tuple(po[i]+u*px[i]+v*py[i]+w*pn[i] for i in range(3))
    kind = case.form[0]
    if kind == 'box':
        pts = pc.boundaries[0].points
        corners = [local(pworld(mp.mpf(u), mp.mpf(v), mp.mpf(w))) for u, v in pts for w in (pc.start, pc.end)]
        box = []
        for i in range(3):
            vals = sorted(c[i] for c in corners)
            assert vals[3]-vals[0] < mp.mpf(10)**-14 and vals[7]-vals[4] < mp.mpf(10)**-14, 'an aligned box'
            box.append((vals[0], vals[7]))
        f = box_form(R, r, box)
        Vc, Mc, K_in, P_in = f['volume'], f['moments'], f['torus_in'], f['box_in']
    else:
        # The torus less the parts beyond each listed face (parallel faces,
        # their outer half-spaces apart).
        _, faces = case.form
        Vc, Mc, K_in, P_in = VK, (mp.mpf(0),)*3, AK, mp.mpf(0)
        pts = [(mp.mpf(u), mp.mpf(v)) for u, v in pc.boundaries[0].points]
        for tag, idx in faces:
            if tag == 'wall':
                p, q = pts[idx], pts[(idx+1) % len(pts)]
                e = (q[0]-p[0], q[1]-p[1])
                # The outer normal of a counter-clockwise profile's edge.
                a = tuple(e[1]*px[i]-e[0]*py[i] for i in range(3))
                P0 = pworld(p[0], p[1], mp.mpf(0))
            else:
                sg = -1 if idx == 0 else 1
                a = tuple(sg*c for c in pn)
                P0 = pworld(mp.mpf(0), mp.mpf(0), mp.mpf(pc.start if idx == 0 else pc.end))
            ac = tuple(ref.dot(a, ax) for ax in (X, Y, Nn))
            f = halfspace_form(R, r, ac, ref.dot(a, P0)-ref.dot(a, ko))
            Vc -= f['volume']
            Mc = tuple(m-x for m, x in zip(Mc, f['moments']))
            K_in -= f['torus_in']
            P_in += f['face_in']
    McW = world_m(Vc, Mc)
    parts = {'K': (VK, MK, AK, K_in), 'P': (VP, MP, AP, P_in)}
    first, second = ('K', 'P') if K_role == 'A' else ('P', 'K')
    VA, MA, AA, A_in = parts[first]
    VB, MB, AB, B_in = parts[second]
    ops = {'common': (Vc, McW, A_in+B_in),
           'fuse': (VA+VB-Vc, tuple(p+q-c for p, q, c in zip(MA, MB, McW)), AA-A_in+AB-B_in),
           'cut': (VA-Vc, tuple(p-c for p, c in zip(MA, McW)), AA-A_in+B_in)}
    out = {}
    for op, (V, m, A) in ops.items():
        if V <= mp.mpf(10)**-30:
            out[op] = None
            continue
        out[op] = (V, A, tuple(v/V for v in m))
    return out


# ------------------------------------------------------------------ Monte Carlo

class FloatTorus:
    def __init__(self, T):
        self.o = [float(v) for v in T.o]
        self.inv = [[float(c) for c in row] for row in T.inv]
        self.R, self.r = float(T.R), float(T.r)
        x, y, n = T.x, T.y, T.n
        e = self.R+self.r
        pts = []
        for su in (-1, 1):
            for sv in (-1, 1):
                for sw in (-1, 1):
                    pts.append([self.o[i]+su*e*float(x[i])+sv*e*float(y[i])+sw*self.r*float(n[i]) for i in range(3)])
        self.box = ([min(p[i] for p in pts) for i in range(3)], [max(p[i] for p in pts) for i in range(3)])

    def contains(self, X):
        d = [X[i]-self.o[i] for i in range(3)]
        u, v, w = (sum(row[i]*d[i] for i in range(3)) for row in self.inv)
        s = u*u+v*v+w*w+self.R*self.R-self.r*self.r
        return s*s < 4*self.R*self.R*(u*u+v*v)


def monte_carlo(pair, n, seed):
    """Volume and centre estimates of the three operations and their
    standard errors, from `n` uniform points in the inputs' box."""
    fk, fp = FloatTorus(pair.T), FloatPrism(pair.P)
    box0 = [min(fk.box[0][i], fp.box[0][i]) for i in range(3)]
    box1 = [max(fk.box[1][i], fp.box[1][i]) for i in range(3)]
    vbox = (box1[0]-box0[0])*(box1[1]-box0[1])*(box1[2]-box0[2])
    rng = random.Random(seed)
    sets = {'fuse': lambda a, b: a or b, 'cut': lambda a, b: a and not b, 'common': lambda a, b: a and b}
    acc = {op: [0, [0.0]*3, [0.0]*3] for op in ref.OPS}
    for _ in range(n):
        X = [box0[i]+(box1[i]-box0[i])*rng.random() for i in range(3)]
        s, p = fk.contains(X), fp.contains(X)
        a, b = (s, p) if pair.torus_first else (p, s)
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


# ------------------------------------------------------------------ margins

def margins(pair):
    """The pair's distances from tangency and incidence, relative to its
    size: the prism's vertices from the torus, each face's plane from the
    tangent planes of its normal whose points of contact lie on the face,
    and the least sine of the angle at which an edge crosses the torus."""
    T, P, size = pair.T, pair.P, pair.size
    R, r = T.Rm, T.rm
    pc = {k: T.chart(X) for k, X in P.P.items()}
    pm = {k: ref.Mv(p) for k, p in pc.items()}
    vertex = min(abs(mp.sqrt((mp.sqrt(p[0]**2+p[1]**2)-R)**2+p[2]**2)-r) for p in pm.values())/size
    face = mp.inf
    for piece in P.pieces:
        for f in piece.faces:
            loop = [pm[k] for k in f['loop']]
            a = ref.cross(ref.sub(loop[1], loop[0]), ref.sub(loop[2], loop[0]))
            na = mp.sqrt(ref.dot(a, a))
            N = tuple(c/na for c in a)
            h = ref.dot(N, loop[0])
            sb = mp.sqrt(N[0]**2+N[1]**2)
            if sb < mp.mpf(10)**-30:
                face = min(face, min(abs(h-c) for c in (r, -r))/size)
                continue
            g = mp.atan2(N[1], N[0])
            cb = N[2]
            cands = []
            for th in (g, g+mp.pi):
                ct = mp.cos(th-g)
                for sg in (1, -1):
                    # The normal (cos phi cos th, cos phi sin th, sin phi) = sg N.
                    cphi, sphi = sg*sb*ct, sg*cb
                    X = ((R+r*cphi)*mp.cos(th), (R+r*cphi)*mp.sin(th), r*sphi)
                    cands.append(X)
            for X in cands:
                dist = ref.dot(N, X)-h
                Xp = tuple(X[i]-dist*N[i] for i in range(3))
                # On the face: inside every edge of the loop (convex).
                n = len(loop)
                inside = all(ref.dot(ref.cross(ref.sub(loop[(k+1) % n], loop[k]), ref.sub(Xp, loop[k])), N) >= 0
                             for k in range(n)) or \
                    all(ref.dot(ref.cross(ref.sub(loop[(k+1) % n], loop[k]), ref.sub(Xp, loop[k])), N) <= 0
                        for k in range(n))
                if inside:
                    face = min(face, abs(dist)/size)
    edge = mp.inf
    edges = set(e for piece in P.pieces for e in piece.edges)
    for i, j in edges:
        vi, vj = pc[i], pc[j]
        e = ref.sub(vj, vi)
        em = ref.Mv(e)
        for t in ref.real_roots(T.line_quartic(vi, e)):
            if not 0 <= t <= 1:
                continue
            X = tuple(pm[i][k]+t*em[k] for k in range(3))
            s2 = X[0]**2+X[1]**2+X[2]**2+R*R-r*r
            grad = (4*s2*X[0]-8*R*R*X[0], 4*s2*X[1]-8*R*R*X[1], 4*s2*X[2])
            ng = mp.sqrt(ref.dot(grad, grad))
            ne = mp.sqrt(ref.dot(em, em))
            edge = min(edge, abs(ref.dot(grad, em))/(ng*ne))
    return {'vertex': vertex, 'face': face, 'edge_crossing': edge}


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
    m = pair.swept()
    VK, MK, AK = pair.T.closed()
    VP, MP, AP = pair.P.closed()
    checks = {}
    dev = mp.mpf(0)
    for res_ in (r, m):
        for (v, mm), (V, Mm) in ((res_['K'], (VK, MK)), (res_['P'], (VP, MP))):
            dev = max(dev, abs(v-V)/size**3, max(abs(x-y) for x, y in zip(mm, Mm))/size**4)
    checks['inputs_sliced'] = dev
    vc, mc = r['common2']
    dev = abs(r['common'][0]-vc)/size**3
    dev = max(dev, max(abs(x-y) for x, y in zip(r['common'][1], mc))/size**4)
    checks['common_two_ways'] = dev
    dev = mp.mpf(0)
    for key, (V, Mm) in (('fuse', (VK+VP-vc, [a+b-c for a, b, c in zip(MK, MP, mc)])),
                         ('K-P', (VK-vc, [a-c for a, c in zip(MK, mc)])),
                         ('P-K', (VP-vc, [a-c for a, c in zip(MP, mc)]))):
        v, mm = r[key]
        dev = max(dev, abs(v-V)/size**3, max(abs(x-y) for x, y in zip(mm, Mm))/size**4)
    checks['inclusion_exclusion'] = dev
    faces = pair.faces()
    checks['face_classes'] = max(abs(sum(cls.values())-area) for _, _, cls, area in faces)/size**2
    checks['area_identity'] = abs(pair.area('fuse')+pair.area('common')-AK-AP)/size**2
    checks['input_areas'] = max(abs(sum(sum(cls.values()) for w, _, cls, _ in faces if w == k)-A)
                                for k, A in (('K', AK), ('P', AP)))/size**2
    dev = mp.mpf(0)
    for key in ('K', 'P', 'common', 'fuse', 'K-P', 'P-K'):
        dev = max(dev, abs(m[key][0]-r[key][0])/size**3,
                  max(abs(x-y) for x, y in zip(m[key][1], r[key][1]))/size**4)
    checks['second_direction'] = dev
    checks['wall_two_ways'] = max(abs(r['wall'][k]-m['wall'][k]) for k in ('in', 'out'))/size**2
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
    near = near_coincidences(pair)
    stats = {'volume_quadrature': max(pair.normal.quad_error, pair.meridian.quad_error)/size**4,
             'breaks': (len(pair.normal.breaks), len(pair.meridian.breaks)), 'margins': margins(pair),
             'solids': pair.meridian.solids()}
    return name, rows, res, checks, near, stats


def near_coincidences(pair):
    """A face class, or a result's volume, positive but thinner than 1e-9 of
    the case's size; slicing breakpoints of either direction closer than
    that (and farther than 1e-30)."""
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
    for label, b, unit in (('levels', pair.normal.raw_breaks, size), ('angles', pair.meridian.raw_breaks, 1)):
        b = sorted(b)
        for p, q in zip(b, b[1:]):
            if mp.mpf(10)**-30*unit < q-p < mp.mpf(10)**-9*unit:
                out.append(f'{label} {mp.nstr(p, 12)} and {mp.nstr(q, 12)}')
    return out


def run(jobs, fn, workers):
    if workers <= 1:
        return [fn(j) for j in jobs]
    with ProcessPoolExecutor(workers) as pool:
        return list(pool.map(fn, jobs))


def closed_job(first):
    return first.pair_name, pair_forms(first)


LIMITS = {'closed_forms_exact_frames': 1e-30, 'closed_forms_turned_frames': 1e-15, 'closed_forms_degenerate': 1e-25,
          'inputs_sliced': 1e-30, 'common_two_ways': 1e-30, 'inclusion_exclusion': 1e-30,
          'face_classes': 1e-30, 'area_identity': 1e-30, 'input_areas': 1e-30,
          'second_direction': 1e-30, 'wall_two_ways': 1e-30, 'monte_carlo_sigma': 5,
          'volume_quadrature_estimate': 1e-30}


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
        if any(c.kind == 'degenerate' for c in cases() if c.pair_name == name):
            # Kinks of the forms within rounding of each other (a vertex
            # 2e-17 off the torus: features of its square root's width).
            kind = 'degenerate'
        for op, want in form.items():
            n, vol, area, centre = res[op]
            if want is None:
                assert n == 0, (name, op, 'the closed form is empty')
                continue
            V, A, C = want
            size = max(abs(V)**(mp.mpf(1)/3), 1)
            dev = max(abs(vol-V)/abs(V), abs(area-A)/abs(A), max(abs(centre[i]-C[i]) for i in range(3))/size)
            note(f'closed_forms_{kind}_frames' if kind != 'degenerate' else 'closed_forms_degenerate', dev)
    for name, rows, res, checks, near, stats in results:
        for key, value in checks.items():
            note(key, value)
        note('volume_quadrature_estimate', stats['volume_quadrature'])
    return worst, covered


def generate(results):
    by_pair = {r[0]: r for r in results}
    blocks = []
    out = [f'# case\trow ({STEP}, torus_boolean_reference.py: expect KIND {STEP}, reason TEXT for a degenerate '
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


def jobs(mc_n, only=None):
    pairs = {}
    for c in cases():
        pairs.setdefault(c.pair_name, [c, []])[1].append(c.operation)
    return [(name, first, ops, mc_n) for name, (first, ops) in pairs.items() if only is None or name in only]


def validate(listed):
    names = [c.name for c in listed]
    assert len(names) == len(set(names)), 'duplicate case names'
    for c in listed:
        assert (c.obj.torus is None) != (c.tool.torus is None), f'{c.name}: a torus and a prism'
        for f in c.frames:
            frame_name(f)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--workers', type=int, default=max(1, (os.cpu_count() or 2)-2))
    parser.add_argument('--samples', type=int, default=200000, help='Monte-Carlo points per pair')
    parser.add_argument('--only', nargs='*', help='evaluate these pairs only (no files written)')
    args = parser.parse_args()
    listed = cases()
    validate(listed)
    if args.only:
        results = run(jobs(args.samples, set(args.only)), evaluate, args.workers)
        for name, rows, res, checks, near, stats in results:
            print(name, rows, {k: mp.nstr(v, 3) for k, v in checks.items()}, near, stats)
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
    for c in listed:
        print(c.name, by_pair[c.pair_name][1][c.operation][0])


if __name__ == '__main__':
    main()
