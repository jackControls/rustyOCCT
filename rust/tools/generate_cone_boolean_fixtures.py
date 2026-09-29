#!/usr/bin/env python3
"""Fixtures for S9d.3a of REVIEW_NOTES.md: Booleans of a cone or frustum
(`Solid::cone_with`) against a polyhedral prism in any relative position,
before any of its kernel code.

`boolean-cone-cases.txt` lists each case in the Boolean protocol
(`identity_reference.encode_boolean_case`: a cone's block is its `frame`
and `cone BOTTOM TOP HEIGHT` rows); `boolean-cone-expected.tsv` gives per
case `expect KIND S9d.3a` (`solid`, `empty` or `degenerate`, then `reason
TEXT`) and `result N volume area cx cy cz` or `empty` from
`cone_boolean_reference.py`, as the sphere fixtures;
`boolean-cone-frames.tsv` the stored axes' bits. Frames are the curved
generator's (`FRAMES`, their stored axes the kernel's bit for bit).

Pairs, with `F` the frustum of radii 2 and 1 and height 2 and `A` the cone
of radius 1.5 and height 3 (its apex on top), both on `XY` at the origin
unless said: `F` cut by a prism's face oblique to its axis in an ellipse
(`ellipse`: a prism along `x` whose top face is `z = 1 + y/4`), by a face
parallel to a ruling in a parabola clipped by both end planes
(`parabola`: `z = 5/2 - 2y`), by a box's wall parallel to the axis in a
hyperbola (`hyperbola`, the box first); a box through `A`'s apex region,
its bottom face a circle and a wall a hyperbola (`apex_box`), a box's
corner below the apex (`apex_corner`, two hyperbolic walls); a cone inside
a box (`cone_in_box`, a cavity) and a box inside a frustum
(`box_in_frustum`, a cavity); `F` and a half-space normal to its axis, a
frustum (`normal_cut`); a slab parallel to the axis cutting `F` in two
(`slab`); a bar along `x` through `F` (`bar`, the bar less the cone two
solids); a box on `F`'s top disc (`on_top`: coplanar, opposite; the common
empty) and a box standing on its base plane (`base_wall`: coplanar, the
same orientation, a hyperbolic wall); turned frames: `F` in `TILT` below a
box's face (`tilt_ellipse`, the section crossing the base), `A` in `LEAN`
with its apex in a box (`lean_apex`, an oblique cone of an ellipse), a
frustum in `TILTX` inside a box (`inside_tilt`), a slab across `F` in
`TILTX` (`tilt_slab`, two solids); and `degenerate` with reasons: a
box's wall through `A`'s apex (`apex_plane`, lines), a face tangent to `F`
along a ruling (`ruling_tangent`), a box's vertex on `F`'s wall
(`vertex_on`).

Before writing, `reference_checks` compares the reference with independent
results (limits relative to the case's size):

* closed forms of every pair, in the cone's ideal (orthonormal) frame, by
  one quadrature (mpmath's tanh-sinh between the kinks) of exact 2D forms
  along the axis: an aligned box by S9d.1's rectangle-in-disc
  antiderivatives and angles of a circle inside a rectangle (the box's
  walls inside by their chords, `|[v0, v1] n [-q, q]|`, `q^2 = r(w)^2 -
  u^2`: hyperbolas), a half-space by circular segments (`r^2 acos(d/r) - d
  sqrt(r^2 - d^2)`, its moment `2/3 (r^2 - d^2)^(3/2)`, the wall's angle `2
  acos(d/r)`, the face's chord `2 sqrt(r^2 - d^2)` times `|a| / |a_uv|`)
  or, normal to the axis, a frustum's closed form, a slab as one
  half-space less another, a solid inside the other; the wall's element `sqrt(1 + k^2) r dtheta dw`; so the common's
  volume, moments and both inputs' surfaces inside the other, and every
  operation by inclusion and exclusion (coplanar faces by hand); within
  1e-30 in exact frames and 1e-15 in turned ones (their stored axes not
  exactly orthonormal);
* `fuse = A + B - common` and `cut = A - common` (each operation sliced
  apart, `A` and `B` in closed form), the common a second way (each piece
  clipped by the disc), both inputs' slicings against their closed forms,
  the area identity `area(fuse) + area(common) + 2 opp = area(A) +
  area(B)`, every face's classes summing to its area, both sides' shared
  areas equal, a second slicing direction `(1, -2, m)` in the cone's chart
  (`m` the least integer with `2 sqrt(5) |k| < m`: every section an
  ellipse cut by the end planes' lines; pairs without coplanar faces),
  within 1e-30;
* Monte-Carlo estimates (200,000 uniform points per pair) of every volume
  and centre within 5 standard errors;
* a scan for near coincidences (a class of a face or a result's volume
  positive but below 1e-9 of the case's size, slicing breakpoints closer
  than that) must find none but in the declared degenerate pairs.

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
import cone_boolean_reference as ref
from generate_curved_boolean_fixtures import FRAMES
from generate_sphere_boolean_fixtures import FloatPrism, arc_in_rect, frame_name, ideal_axes, rect_disc

ROOT = Path(__file__).resolve().parents[1]
PREFIX = 'boolean-cone'
STEP = 'S9d.3a'
BOOLEAN_OPERATION = 93
EXACT = {'XY', 'SIDE', 'DOWN', 'TURN'}


def at(name, origin):
    return tuple(float(c) for c in origin)+FRAMES[name]


def square(x0, y0, x1, y1):
    return Boundary(points=[(x0, y0), (x1, y0), (x1, y1), (x0, y1)])


def polygon(*pts):
    return Boundary(points=[(float(u), float(v)) for u, v in pts])


def make(spec, op):
    if spec[0] == 'cone':
        _, r0, r1, h, frame = spec
        return Case('', 1e-7, op, frame, 0.0, 0.0, [], cone=(float(r0), float(r1), float(h)))
    _, boundaries, frame, start, end = spec
    return Case('', 1e-7, op, frame, float(start), float(end), boundaries)


def cone(r0, r1, h, frame):
    return ('cone', r0, r1, h, frame)


def prism(boundaries, frame, start, end):
    return ('prism', boundaries, frame, start, end)


class Boolean:
    def __init__(self, name, operation, obj, tool, form, kind='solid', reason=None, coplanar=False):
        self.name, self.operation, self.kind, self.reason = name, operation, kind, reason
        self.form, self.coplanar = form, coplanar
        self.pair_name = name.rsplit('_', 1)[0]
        self.obj = make(obj, 91)
        self.obj.name = name
        self.tool = make(tool, 92)
        self.tool.name = name
        self.frames = (obj[4] if obj[0] == 'cone' else obj[2], tool[4] if tool[0] == 'cone' else tool[2])

    def encode(self):
        return encode_boolean_case(self.obj, self.operation, self.tool, BOOLEAN_OPERATION)


def group(name, obj, tool, ops, form, **kw):
    out = []
    for op, kind in ops.items():
        reason = None
        if isinstance(kind, tuple):
            kind, reason = kind
        out.append(Boolean(f'{name}_{op}', op, obj, tool, form, kind, reason, **kw))
    return out


APEX_PLANE = "a face of the prism through the cone's apex"
RULING = 'a face of the prism tangent to the cone along a ruling'
VERTEX_ON = "a vertex of the prism on the cone's wall"


def cases():
    """S9d.3a's pairs (see the module's docstring)."""
    O = (0, 0, 0)
    F_ = cone(2, 1, 2, at('XY', O))
    A_ = cone(1.5, 0, 3, at('XY', O))
    out = []
    out += group('ellipse', F_, prism([polygon((-8, -1), (4, -1), (4, 2))], at('SIDE', O), -4, 4),
                 {'fuse': 'solid', 'cut': 'solid', 'common': 'solid'}, ('half', ('wall', 2)))
    out += group('parabola', F_, prism([polygon((2.75, -3), (4, -3), (4, 3.5), (-0.5, 3.5))], at('SIDE', O), -4, 4),
                 {'cut': 'solid', 'common': 'solid'}, ('half', ('wall', 3)))
    out += group('hyperbola', prism([square(1.25, -4, 4, 4)], at('XY', O), -1, 3), F_,
                 {'cut': 'solid', 'common': 'solid'}, ('half', ('wall', 3)))
    out += group('apex_box', A_, prism([square(-0.25, -1, 1, 1)], at('XY', O), 1.5, 4),
                 {'cut': 'solid', 'common': 'solid'}, ('box',))
    out += group('apex_corner', A_, prism([square(0.125, -0.25, 2, 2)], at('XY', O), 2, 4),
                 {'common': 'solid'}, ('box',))
    out += group('cone_in_box', prism([square(-2, -2, 3, 3)], at('XY', O), -2, 2),
                 cone(1, 0, 2, at('XY', (0.5, 0.25, -1))), {'cut': 'solid'}, ('inside',))
    out += group('box_in_frustum', cone(3, 2, 3, at('XY', O)), prism([square(-1, -1, 1, 1)], at('XY', O), 0.5, 2),
                 {'cut': 'solid'}, ('contains',))
    out += group('normal_cut', F_, prism([square(-3, -3, 3, 3)], at('XY', O), 1, 5),
                 {'fuse': 'solid', 'cut': 'solid', 'common': 'solid'}, ('half', ('cap', 0)))
    out += group('slab', F_, prism([square(-0.25, -3, 0.25, 3)], at('XY', O), -1, 3), {'cut': 'solid'}, ('box',))
    out += group('bar', prism([square(-0.25, 0.75, 0.25, 1.25)], at('SIDE', O), -4, 4), F_, {'cut': 'solid'},
                 ('box',))
    out += group('on_top', F_, prism([square(-3, -3, 3, 3)], at('XY', O), 2, 3),
                 {'fuse': 'solid', 'common': 'empty'}, ('box',), coplanar=True)
    out += group('base_wall', prism([square(0.5, -3, 3, 3)], at('XY', O), 0, 1), F_, {'common': 'solid'}, ('box',),
                 coplanar=True)
    out += group('tilt_ellipse', cone(2, 1, 2, at('TILT', O)), prism([square(-5, -5, 5, 5)], at('XY', O), -4, 0.5),
                 {'cut': 'solid', 'common': 'solid'}, ('half', ('cap', 1)))
    out += group('lean_apex', cone(1.5, 0, 3, at('LEAN', (1, 2, -1))),
                 prism([square(0, -1, 5, 5)], at('XY', O), 0.5, 4), {'cut': 'solid', 'common': 'solid'},
                 ('half', ('cap', 0)))
    out += group('inside_tilt', prism([square(-3, -3, 3, 3)], at('XY', O), -3, 3),
                 cone(1.5, 0.5, 2, at('TILTX', (0.25, 0.5, -0.5))), {'cut': 'solid'}, ('inside',))
    out += group('tilt_slab', cone(2, 1, 2, at('TILTX', O)), prism([square(-5, -5, 5, 5)], at('XY', O), 0.25, 0.75),
                 {'cut': 'solid'}, ('slab',))
    out += group('apex_plane', A_, prism([square(0, -3, 3, 3)], at('XY', O), -1, 4),
                 {'common': ('degenerate', APEX_PLANE)}, ('box',))
    out += group('ruling_tangent', F_, prism([polygon((3.5, -3), (5, -3), (5, 5), (-0.5, 5))], at('SIDE', O), -4, 3),
                 {'fuse': ('degenerate', RULING), 'common': ('degenerate', RULING)}, ('half', ('wall', 3)))
    out += group('vertex_on', prism([square(0.75, 1, 4, 4)], at('XY', O), -1, 1.5), F_,
                 {'cut': ('degenerate', VERTEX_ON)}, ('box',))
    return out


def exact_pair(case):
    return all(frame_name(f) in EXACT for f in case.frames)


# ------------------------------------------------------------------ closed forms

def cone_measures(r0, r1, h):
    """Volume, axial moment (from the base) and wall's and total area of a
    cone or frustum."""
    r0, r1, h = mp.mpf(r0), mp.mpf(r1), mp.mpf(h)
    V = mp.pi*h*(r0*r0+r0*r1+r1*r1)/3
    Mw = mp.pi*h*h*(r0*r0+2*r0*r1+3*r1*r1)/12
    wall = mp.pi*(r0+r1)*mp.sqrt(h*h+(r1-r0)**2)
    return V, Mw, wall, wall+mp.pi*(r0*r0+r1*r1)


def segment(r, d):
    """The disc of radius `r` beyond the line at signed distance `d` from its
    centre: area, first moment along the line's normal, the circle's angle
    beyond it, the chord's length."""
    if r <= 0 or d >= r:
        return mp.mpf(0), mp.mpf(0), mp.mpf(0), mp.mpf(0)
    if d <= -r:
        return mp.pi*r*r, mp.mpf(0), 2*mp.pi, mp.mpf(0)
    q = mp.sqrt(r*r-d*d)
    return r*r*mp.acos(d/r)-d*q, 2*q**3/3, 2*mp.acos(d/r), 2*q


def halfspace_form(r0, r1, h, a, b):
    """The cone (ideal frame) and the half-space `a . p >= b`: the common's
    volume and moments, the cone's surface inside, the plane's face inside."""
    r0, r1, h = mp.mpf(r0), mp.mpf(r1), mp.mpf(h)
    k = (r1-r0)/h
    r = lambda w: r0+k*w
    au, av, aw = a
    b = mp.mpf(b)
    nuv = mp.sqrt(au*au+av*av)
    sk = mp.sqrt(1+k*k)
    if nuv < mp.mpf(10)**-30:
        w0 = b/aw
        lo, hi = (max(w0, 0), h) if aw > 0 else (mp.mpf(0), min(w0, h))
        if hi <= lo:
            return {'volume': mp.mpf(0), 'moments': (mp.mpf(0),)*3, 'cone_in': mp.mpf(0), 'face_in': mp.mpf(0)}
        V = mp.pi*((r0*r0*hi+r0*k*hi**2+k*k*hi**3/3)-(r0*r0*lo+r0*k*lo**2+k*k*lo**3/3))
        Mw = mp.pi*((r0*r0*hi**2/2+2*r0*k*hi**3/3+k*k*hi**4/4)-(r0*r0*lo**2/2+2*r0*k*lo**3/3+k*k*lo**4/4))
        cone_in = sk*mp.pi*(r(lo)+r(hi))*(hi-lo)
        for w in (mp.mpf(0), h):
            if lo <= w <= hi and (w != w0):
                cone_in += mp.pi*r(w)**2
        face_in = mp.pi*r(w0)**2 if 0 < w0 < h else mp.mpf(0)
        return {'volume': V, 'moments': (mp.mpf(0), mp.mpf(0), Mw), 'cone_in': cone_in, 'face_in': face_in}
    mu, mv = au/nuv, av/nuv
    d = lambda w: (b-aw*w)/nuv
    kinks = [mp.mpf(0), h]
    for sg in (1, -1):
        den = aw+sg*nuv*k
        if den != 0:
            w = (b-sg*nuv*r0)/den
            if 0 < w < h:
                kinks.append(w)
    kinks.sort()
    seg = lambda w: segment(r(w), d(w))
    V = mp.quad(lambda w: seg(w)[0], kinks)
    Mn = mp.quad(lambda w: seg(w)[1], kinks)
    Mw = mp.quad(lambda w: w*seg(w)[0], kinks)
    wall = sk*mp.quad(lambda w: r(w)*seg(w)[2], kinks)
    face = mp.sqrt(au*au+av*av+aw*aw)/nuv*mp.quad(lambda w: seg(w)[3], kinks)
    ends = segment(r0, d(0))[0]+segment(r1, d(h))[0]
    return {'volume': V, 'moments': (mu*Mn, mv*Mn, Mw), 'cone_in': wall+ends, 'face_in': face}


def box_form(r0, r1, h, box):
    """The cone (ideal frame) and an aligned box `[u0, u1] x [v0, v1] x [w0,
    w1]`: the common's volume and moments, the cone's surface inside, the
    box's faces inside, and the areas of the box's faces on the cone's end
    planes inside its discs ({(w, outward sign): area})."""
    r0, r1, h = mp.mpf(r0), mp.mpf(r1), mp.mpf(h)
    (u0, u1), (v0, v1), (w0, w1) = [tuple(mp.mpf(x) for x in rng) for rng in box]
    k = (r1-r0)/h
    r = lambda w: max(r0+k*w, mp.mpf(0))
    sk = mp.sqrt(1+k*k)
    wa, wb = max(w0, mp.mpf(0)), min(w1, h)
    out = {'volume': mp.mpf(0), 'moments': (mp.mpf(0),)*3, 'cone_in': mp.mpf(0), 'box_in': mp.mpf(0),
           'coplanar': {}}
    for w, sign in ((w0, -1), (w1, 1)):
        if w in (0, h):
            out['coplanar'][(w, sign)] = rect_disc(u0, u1, v0, v1, r(w))[0]
    if wb <= wa:
        return out
    kinks = set([wa, wb])
    vals = [abs(u0), abs(u1), abs(v0), abs(v1)]+[mp.sqrt(x*x+y*y) for x in (u0, u1) for y in (v0, v1)]
    for val in vals:
        w = (val-r0)/k
        if wa < w < wb:
            kinks.add(w)
    pts = sorted(kinks)
    rd = lambda w: rect_disc(u0, u1, v0, v1, r(w))
    V = mp.quad(lambda w: rd(w)[0], pts)
    Mu = mp.quad(lambda w: rd(w)[1], pts)
    Mv = mp.quad(lambda w: rd(w)[2], pts)
    Mw = mp.quad(lambda w: w*rd(w)[0], pts)
    out['volume'], out['moments'] = V, (Mu, Mv, Mw)
    cone_in = sk*mp.quad(lambda w: r(w)*arc_in_rect(u0, u1, v0, v1, r(w)), pts)
    for w, rr in ((mp.mpf(0), r0), (h, r1)):
        if w0 < w < w1:
            cone_in += rect_disc(u0, u1, v0, v1, rr)[0]
    out['cone_in'] = cone_in

    def chord(c, lo, hi, w):
        q2 = r(w)**2-c*c
        if q2 <= 0:
            return mp.mpf(0)
        q = mp.sqrt(q2)
        return max(mp.mpf(0), min(hi, q)-max(lo, -q))
    box_in = mp.mpf(0)
    for u in (u0, u1):
        box_in += mp.quad(lambda w: chord(u, v0, v1, w), pts)
    for v in (v0, v1):
        box_in += mp.quad(lambda w: chord(v, u0, u1, w), pts)
    for w, sign in ((w0, -1), (w1, 1)):
        if 0 < w < h:
            box_in += rect_disc(u0, u1, v0, v1, r(w))[0]
    out['box_in'] = box_in
    return out


def prism_measures(pc):
    """Volume, moments and area (world, ideal frame) of a prism of one
    polygon."""
    x, y, n = ideal_axes(frame_name(pc.frame))
    o = tuple(mp.mpf(v) for v in pc.frame[:3])
    pts = [(mp.mpf(u), mp.mpf(v)) for u, v in pc.boundaries[0].points]
    N = len(pts)
    cr = [pts[k-1][0]*pts[k][1]-pts[k][0]*pts[k-1][1] for k in range(N)]
    A = sum(cr)/2
    cu = sum((pts[k-1][0]+pts[k][0])*cr[k] for k in range(N))/(6*A)
    cv = sum((pts[k-1][1]+pts[k][1])*cr[k] for k in range(N))/(6*A)
    per = sum(mp.sqrt((pts[k][0]-pts[k-1][0])**2+(pts[k][1]-pts[k-1][1])**2) for k in range(N))
    L = mp.mpf(pc.end)-mp.mpf(pc.start)
    V = A*L
    cw = (mp.mpf(pc.start)+mp.mpf(pc.end))/2
    c = tuple(o[i]+cu*x[i]+cv*y[i]+cw*n[i] for i in range(3))
    return V, tuple(V*ci for ci in c), 2*A+per*L


def pair_forms(case):
    """{op: (volume, area, centre)} in closed form for a pair (see the
    module's docstring), in the cone's ideal frame."""
    K_role = 'A' if case.obj.cone is not None else 'B'
    kc, pc = (case.obj, case.tool) if K_role == 'A' else (case.tool, case.obj)
    X, Y, Nn = ideal_axes(frame_name(kc.frame))
    ko = tuple(mp.mpf(v) for v in kc.frame[:3])
    r0, r1, h = kc.cone
    local = lambda P: tuple(ref.dot(ref.sub(P, ko), ax) for ax in (X, Y, Nn))
    world_m = lambda V, m: tuple(ko[j]*V+sum(m[i]*(X, Y, Nn)[i][j] for i in range(3)) for j in range(3))
    VK, MwK, _, AK = cone_measures(r0, r1, h)
    MK = world_m(VK, (0, 0, MwK))
    VP, MP, AP = prism_measures(pc)
    px, py, pn = ideal_axes(frame_name(pc.frame))
    po = tuple(mp.mpf(v) for v in pc.frame[:3])
    pworld = lambda u, v, w: tuple(po[i]+u*px[i]+v*py[i]+w*pn[i] for i in range(3))
    on_same = on_opp = mp.mpf(0)
    kind = case.form[0]
    if kind == 'inside':
        Vc, Mc, K_in, P_in = VK, (0, 0, MwK), AK, mp.mpf(0)
    elif kind == 'contains':
        Vc, Mc, K_in, P_in = None, None, mp.mpf(0), AP
    elif kind == 'box':
        pts = pc.boundaries[0].points
        corners = [local(pworld(mp.mpf(u), mp.mpf(v), mp.mpf(w))) for u, v in pts for w in (pc.start, pc.end)]
        box = []
        for i in range(3):
            vals = sorted(c[i] for c in corners)
            assert vals[3]-vals[0] < mp.mpf(10)**-14 and vals[7]-vals[4] < mp.mpf(10)**-14, 'an aligned box'
            box.append((vals[0], vals[7]))
        f = box_form(r0, r1, h, box)
        Vc, Mc, K_in, P_in = f['volume'], f['moments'], f['cone_in'], f['box_in']
        for (w, sign), area in f['coplanar'].items():
            # The cone's end's outward normal: -w at 0, +w at h.
            if (sign > 0) == (w != 0):
                on_same += area
            else:
                on_opp += area
    elif kind == 'slab':
        # Both caps cross the cone, the walls outside: the common is the
        # half-space above the start cap less the one beyond the end cap.
        halves = []
        for level in (pc.start, pc.end):
            P0 = pworld(mp.mpf(0), mp.mpf(0), mp.mpf(level))
            ac = tuple(ref.dot(pn, ax) for ax in (X, Y, Nn))
            halves.append(halfspace_form(r0, r1, h, ac, ref.dot(pn, P0)-ref.dot(pn, ko)))
        lo, hi = halves
        Vc = lo['volume']-hi['volume']
        Mc = tuple(p-q for p, q in zip(lo['moments'], hi['moments']))
        K_in = lo['cone_in']-hi['cone_in']
        P_in = lo['face_in']+hi['face_in']
    else:
        _, (tag, idx) = case.form
        if tag == 'wall':
            pts = [(mp.mpf(u), mp.mpf(v)) for u, v in pc.boundaries[0].points]
            p, q = pts[idx], pts[(idx+1) % len(pts)]
            e = (q[0]-p[0], q[1]-p[1])
            a = tuple(-(e[1]*px[i]-e[0]*py[i]) for i in range(3))
            P0 = pworld(p[0], p[1], mp.mpf(0))
        else:
            sg = 1 if idx == 0 else -1
            a = tuple(sg*c for c in pn)
            P0 = pworld(mp.mpf(0), mp.mpf(0), mp.mpf(pc.start if idx == 0 else pc.end))
        b = ref.dot(a, P0)
        ac = tuple(ref.dot(a, ax) for ax in (X, Y, Nn))
        f = halfspace_form(r0, r1, h, ac, b-ref.dot(a, ko))
        Vc, Mc, K_in, P_in = f['volume'], f['moments'], f['cone_in'], f['face_in']
    if Vc is None:
        Vc, McW = VP, MP
    else:
        McW = world_m(Vc, Mc)
    parts = {'K': (VK, MK, AK, K_in), 'P': (VP, MP, AP, P_in)}
    first, second = ('K', 'P') if K_role == 'A' else ('P', 'K')
    VA, MA, AA, A_in = parts[first]
    VB, MB, AB, B_in = parts[second]
    A_out = AA-A_in-on_same-on_opp
    B_out = AB-B_in-on_same-on_opp
    ops = {'common': (Vc, McW, A_in+B_in+on_same),
           'fuse': (VA+VB-Vc, tuple(p+q-c for p, q, c in zip(MA, MB, McW)), A_out+B_out+on_same),
           'cut': (VA-Vc, tuple(p-c for p, c in zip(MA, McW)), A_out+on_opp+B_in)}
    out = {}
    for op, (V, m, A) in ops.items():
        if V <= mp.mpf(10)**-30:
            out[op] = None
            continue
        out[op] = (V, A, tuple(v/V for v in m))
    return out


# ------------------------------------------------------------------ Monte Carlo

class FloatCone:
    def __init__(self, K):
        self.o = [float(v) for v in K.o]
        self.inv = [[float(c) for c in row] for row in K.inv]
        self.r0, self.k, self.h = float(K.r0), float(K.k), float(K.h)
        x, y, n = K.x, K.y, K.n
        rmax = float(max(K.r0, K.r1))
        pts = []
        for w in (0.0, self.h):
            for su in (-1, 1):
                for sv in (-1, 1):
                    pts.append([self.o[i]+su*rmax*float(x[i])+sv*rmax*float(y[i])+w*float(n[i]) for i in range(3)])
        self.box = ([min(p[i] for p in pts) for i in range(3)], [max(p[i] for p in pts) for i in range(3)])

    def contains(self, X):
        d = [X[i]-self.o[i] for i in range(3)]
        u, v, w = (sum(row[i]*d[i] for i in range(3)) for row in self.inv)
        if not 0 < w < self.h:
            return False
        r = self.r0+self.k*w
        return u*u+v*v < r*r


def monte_carlo(pair, n, seed):
    """Volume and centre estimates of the three operations and their
    standard errors, from `n` uniform points in the inputs' box."""
    fk, fp = FloatCone(pair.K), FloatPrism(pair.P)
    box0 = [min(fk.box[0][i], fp.box[0][i]) for i in range(3)]
    box1 = [max(fk.box[1][i], fp.box[1][i]) for i in range(3)]
    vbox = (box1[0]-box0[0])*(box1[1]-box0[1])*(box1[2]-box0[2])
    rng = random.Random(seed)
    sets = {'fuse': lambda a, b: a or b, 'cut': lambda a, b: a and not b, 'common': lambda a, b: a and b}
    acc = {op: [0, [0.0]*3, [0.0]*3] for op in ref.OPS}
    for _ in range(n):
        X = [box0[i]+(box1[i]-box0[i])*rng.random() for i in range(3)]
        s, p = fk.contains(X), fp.contains(X)
        a, b = (s, p) if pair.cone_first else (p, s)
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

def second_direction(K):
    """`(1, -2, m)`: every slice's section of the cone an ellipse."""
    m = int(math.floor(2*math.sqrt(5)*abs(float(K.k))))+1
    while not 5*K.k*K.k < ref.F(m*m, 4):
        m += 1
    return (1, -2, m)


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
    VK, MK, AK = pair.K.closed()
    VP, MP, AP = pair.P.closed()
    checks = {}
    dev = mp.mpf(0)
    for (v, m), (V, Mm) in ((r['K'], (VK, MK)), (r['P'], (VP, MP))):
        dev = max(dev, abs(v-V)/size**3, max(abs(x-y) for x, y in zip(m, Mm))/size**4)
    checks['inputs_sliced'] = dev
    vc, mc = r['common2']
    dev = abs(r['common'][0]-vc)/size**3
    dev = max(dev, max(abs(x-y) for x, y in zip(r['common'][1], mc))/size**4)
    dev = max(dev, abs(sum(r['pieces'])-vc)/size**3)
    checks['common_two_ways'] = dev
    dev = mp.mpf(0)
    for key, (V, Mm) in (('fuse', (VK+VP-vc, [a+b-c for a, b, c in zip(MK, MP, mc)])),
                         ('K-P', (VK-vc, [a-c for a, c in zip(MK, mc)])),
                         ('P-K', (VP-vc, [a-c for a, c in zip(MP, mc)]))):
        v, m = r[key]
        dev = max(dev, abs(v-V)/size**3, max(abs(x-y) for x, y in zip(m, Mm))/size**4)
    checks['inclusion_exclusion'] = dev
    faces = pair.faces()
    checks['face_classes'] = max(abs(sum(cls.values())-area) for _, _, cls, area in faces)/size**2
    ck, cp = pair.classes('K'), pair.classes('P')
    checks['shared_both_sides'] = max(abs(ck['same']-cp['same']), abs(ck['opp']-cp['opp']))/size**2
    checks['area_identity'] = abs(pair.area('fuse')+pair.area('common')+2*ck['opp']-AK-AP)/size**2
    checks['input_areas'] = max(abs(sum(sum(cls.values()) for w, _, cls, _ in faces if w == k)-A)
                                for k, A in (('K', AK), ('P', AP)))/size**2
    if not first.coplanar:
        second = ref.Slicing(pair.K, pair.P, second_direction(pair.K), size).results()
        dev = mp.mpf(0)
        for key in ('K', 'P', 'common', 'fuse', 'K-P', 'P-K'):
            dev = max(dev, abs(second[key][0]-r[key][0])/size**3,
                      max(abs(x-y) for x, y in zip(second[key][1], r[key][1]))/size**4)
        checks['second_direction'] = dev
    mcr = monte_carlo(pair, mc_n, zlib.crc32(name.encode()))
    z = 0.0
    for op in ref.OPS:
        vol, sv, c, sc = mcr[op]
        V, m = pair.volume(op)
        v = float(V)
        z = max(z, abs(vol-v)/sv)
        if v > 1e-9 and c is not None:
            cent = [float(x/V) for x in m]
            z = max(z, max(abs(c[i]-cent[i])/max(sc[i], 1e-12) for i in range(3)))
    checks['monte_carlo_sigma'] = mp.mpf(z)
    near = near_coincidences(pair)
    stats = {'volume_quadrature': pair.slicing.quad_error/size**4, 'breaks': len(pair.slicing.breaks)}
    return name, rows, res, checks, near, stats


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
    b = pair.slicing.raw_breaks
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
    return first.pair_name, pair_forms(first)


LIMITS = {'closed_forms_exact_frames': 1e-30, 'closed_forms_turned_frames': 1e-15,
          'inputs_sliced': 1e-30, 'common_two_ways': 1e-30, 'inclusion_exclusion': 1e-30,
          'face_classes': 1e-30, 'shared_both_sides': 1e-30, 'area_identity': 1e-30, 'input_areas': 1e-30,
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
    out = [f'# case\trow ({STEP}, cone_boolean_reference.py: expect KIND {STEP}, reason TEXT for a degenerate '
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
        assert (c.obj.cone is None) != (c.tool.cone is None), f'{c.name}: a cone and a prism'
        for f in c.frames:
            frame_name(f)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--workers', type=int, default=max(1, (os.cpu_count() or 2)-2))
    parser.add_argument('--samples', type=int, default=200000, help='Monte-Carlo points per pair')
    args = parser.parse_args()
    listed = cases()
    validate(listed)
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
    by_pair = {r[0]: r for r in results}
    for name, _, _, _, n, stats in results:
        if n:
            print('declared near coincidences:', name, '; '.join(n[:4]))
    for c in listed:
        print(c.name, by_pair[c.pair_name][1][c.operation][0])


if __name__ == '__main__':
    main()
