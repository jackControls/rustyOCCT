#!/usr/bin/env python3
"""Independent reference for S9d.4b.1 of REVIEW_NOTES.md: Booleans of a torus
v-segment or wedge (`Solid::torus_with` other than a whole torus: S3's
segments between two latitudes of the tube, a full turn, and wedges of the
whole tube over a partial turn) against a polyhedral prism (a profile of
lines, every face a plane) in any relative position.

**The inputs' exact models.** In rationals from the stored binary64 data, in
the torus's chart `(u, v, w)` of its stored frame (`stored_axes`, affine where
the axes are not exactly orthonormal), as S9d.4a's (`torus_boolean_
reference.Torus`, whose chart, wall element and weights this module uses).
S3 builds (`topology.rs`, `Topology::torus`, as `BRepPrimAPI_MakeTorus(gp_Ax2,
R, r, low, high, angle)`):

* **A v-segment** (latitudes `low < high <= low + 2 pi`, a full turn) is the
  region between the tube's arc from `low` to `high` and the axis, revolved:
  its end faces are planar discs normal to the axis at the heights `z =
  r sin(latitude)` (`math::scaled_sin`, one binary64 product), each from the
  axis to the arc's end, not a surface swept by the tube's point. In the
  meridian half-plane `(t, w)` (`t >= 0` the distance from the axis) the
  section `M` is bounded by the arc, the two horizontal segments from its
  ends to the axis and the axis between them; the model takes the planes
  `w = z` at the stored heights and the torus exactly, the arc's ends where
  the planes meet the tube on the side of the latitude's cosine (`rho = R
  +- sqrt(r^2 - z^2)`), its angles `v_lo < v_hi` those points' own. `M` is
  the set with an odd number of the arc's points at its height to its right
  (one between the end heights, where `M` reaches the axis: `[0, c]`; none
  or two elsewhere: `[R - q, R + q]`, the tube's cap). Its boundary runs
  counter-clockwise exactly when S3's wall is forward (the arc bulging away
  from the axis); only the arc carries Green's integrals of the volume and
  moments (`t^2 dw`: the axis has `t = 0`, the end segments `dw = 0`).
  The outer half's (`+-pi/2`) end planes `w = +-r` are tangent to the torus
  along its rings (the wall meets its end discs tangentially), the inner
  half's likewise (and OCCT builds that half inside out).
* **A wedge** (the whole tube from `v = 0`, `0 < angle < 2 pi`) is the torus
  `(|p|^2 + R^2 - r^2)^2 <= 4 R^2 (u^2 + v^2)` in the sector from the half-
  plane of `x` (`theta = 0`) to the half-plane of the chart direction `e =
  (cos angle, sin angle)` rounded (`cos_rn`, `sin_rn`: the kernel's end disc
  lies on the frame's `x cos + y sin`, rounded), `theta_U = atan2(e)`: the
  half-turn's end lies 1.2e-16 off the plane of `x`, the quarter's 6.1e-17
  off the plane of `y`. Its end faces are the tube's discs in those
  half-planes.

**Volumes and moments, two ways.** (1) Normal slices `w = s`: `K`'s section
is the union of `M`'s intervals revolved (a disc `[0, c]` or an annulus), cut
to the sector for a wedge (two convex sectors when `theta_U > pi`), each
prism piece's section a convex polygon; the common by clipping each polygon
to the sector and by the discs (S9d.1's `convex_disc`, signed), `K` and `P`
in closed form per slice, the other operations by inclusion and exclusion.
(2) Meridian half-planes `theta`: the section is `M` (a wedge's: the tube's
disc for `theta` in the turn, nothing outside), the prism's a convex polygon,
the boundaries cut and classified (the arc by the polygon's lines, the
polygon's edges by the tube's circle and the end heights, each piece at its
midpoint: `M` by its parity rule), the volume by the cylindrical element as
Green's integrals of `t`, `t^2` and `t w`, every operation apart. Breakpoints:
S9d.4a's in both directions and the end heights, rings, end half-planes and
their circles meeting the prism's edges and faces.

**Areas, every face in its own parameters.** The wall: by the meridians in
its own `(theta, phi)` (the arcs of the tube's circle inside the section
polygon) and by the normal slices (each circle's angle inside the pieces
and the sector, the element `r rho / q`). A segment's end discs: exactly at
their levels (each piece's section there, or a piece's face on the plane:
`same` or `opp` by the pieces' side against `K`'s), and by the meridians (the
radial segment `[0, rho_end]` at the end's height inside the polygon, `t dt
dtheta`). A wedge's end discs: exactly in their half-planes (each piece's
section, clipped at the axis) and by the normal slices (the end ray's chord
`[R - q, R + q]` inside the pieces). The prism's faces: normal to the axis by
the regions of `K` just above and below the level (`in` where both, `same`
or `opp` where one: a face on an end disc's plane), any other by the lines
of constant `w` across it (the chord's parameters inside the discs and the
sector), a face on a wedge's end half-plane by its chords on the end's ray.

**Solids.** By sweeping the meridian half-planes: in each section the
operation's region is cut into bands of `w` between its critical levels
(the end heights, `+-r`, the polygon's vertices, the polygon's edges meeting
the tube's circle), each band's intervals of `t` joined across a level where
they overlap; sections within one interval of `theta` between breakpoints
are sampled densely and joined where their components' sample points lie in
each other (a one-to-one correspondence required, else the samples refined:
no event is missed silently), across a breakpoint where they do also in the
section at the breakpoint itself (regions meeting along a curve or at a
point stay apart: the regularized Boolean), around the whole turn, and
through the axis where two components touch it along overlapping heights.
The fuse is one solid where the inputs overlap or share a face (`opp`), two
otherwise (both inputs connected); the sweep's fuse must agree.

`rows(obj, operation, tool)` gives `result N volume area cx cy cz` (totals
over the N solids, world coordinates) or `empty`, as S9d.4a's. Nothing here
uses the kernel, a surface/surface intersection or an arrangement; shared
with other references are `stored`, `stored_axes`, S9d.1's prism, 2D
regions and quadrature, S9d.3a's classification and S9d.4a's chart, wall
element, polynomial roots and meridian sections.
"""
from bisect import bisect_right
from fractions import Fraction as F
import itertools
import math
from types import SimpleNamespace

import mpmath as mp

from brep_reference import cos_rn, sin_rn
from sphere_boolean_reference import (
    M, Mv, Disc2, Polys2, Prism, add, clip_poly, cross, cross2, dot, green_seg, hull2, integrate, line_of,
    roots2, scale, sub, tau, zero3, vadd)
from cone_boolean_reference import ccw2, classify, convex_disc, seg_circle_params, circle_line_angles
import torus_boolean_reference as whole
from torus_boolean_reference import (
    Torus, circle_tangency_quartic, green_mer_arc, green_mer_seg, padd, pmul, pscale, real_roots, number, OPS,
    TWO_PI)

mp.mp.dps = 40


def in_range(v, lo, hi):
    """Whether the angle `v` (mod 2 pi) lies in `[lo, hi]`."""
    T = tau()
    k = mp.ceil((lo-v)/T)
    return v+k*T <= hi


def iv_intersect(A, B):
    out = []
    for a0, a1 in A:
        for b0, b1 in B:
            lo, hi = max(a0, b0), min(a1, b1)
            if hi > lo:
                out.append((lo, hi))
    return sorted(out)


def iv_diff(A, B):
    out = []
    for a0, a1 in A:
        pieces = [(a0, a1)]
        for b0, b1 in B:
            nxt = []
            for p0, p1 in pieces:
                if b1 <= p0 or b0 >= p1:
                    nxt.append((p0, p1))
                    continue
                if b0 > p0:
                    nxt.append((p0, b0))
                if b1 < p1:
                    nxt.append((b1, p1))
            pieces = nxt
        out += pieces
    return sorted(out)


def iv_union(A, B):
    out = []
    for a0, a1 in sorted(list(A)+list(B)):
        if out and a0 <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], a1))
        else:
            out.append((a0, a1))
    return out


def quad_roots_mp(a, b, c):
    """Real roots (mpf) of `a x^2 + b x + c`."""
    if abs(a) < mp.mpf(10)**-35:
        return [] if b == 0 else [-c/b]
    d = b*b-4*a*c
    if d < 0:
        return []
    s = mp.sqrt(d)
    return [(-b-s)/(2*a), (-b+s)/(2*a)]


# ------------------------------------------------------------------ the part

class Part:
    """A torus v-segment or wedge on its exact model (see the module's note),
    in its chart; `T` the whole torus of its radii on the same frame."""

    def __init__(self, case):
        R, r, low, high, angle = case.torus
        tube, turn = high-low == TWO_PI, angle == TWO_PI
        assert not (tube and turn), 'S9d.4b.1: a v-segment or a wedge, not a whole torus'
        assert tube or turn, 'S3: no segment of a partial turn'
        self.case = case
        self.kind = 'wedge' if tube else 'segment'
        self.T = T = Torus(SimpleNamespace(frame=case.frame, torus=(R, r, 0.0, TWO_PI, TWO_PI)))
        self.R, self.r, self.Rm, self.rm = T.R, T.r, T.Rm, T.rm
        self.o, self.x, self.y, self.n = T.o, T.x, T.y, T.n
        if self.kind == 'segment':
            ends = []
            for lat in (low, high):
                z = F(r*sin_rn(lat))
                side = 1 if cos_rn(lat) > 0 else -1
                zm = M(z)
                q = mp.sqrt(max(self.rm**2-zm**2, mp.mpf(0)))
                v = mp.atan2(zm, side*q)
                v += mp.nint((mp.mpf(lat)-v)/tau())*tau()
                ends.append({'lat': lat, 'z': z, 'zm': zm, 'side': side, 'q': q, 'rho': self.Rm+side*q, 'v': v})
            self.ends = ends
            self.vlo, self.vhi = ends[0]['v'], ends[1]['v']
            assert self.vlo < self.vhi <= self.vlo+tau()
            R_, r_ = self.Rm, self.rm
            area = r_*(R_*(mp.sin(self.vhi)-mp.sin(self.vlo))
                       + r_*((self.vhi-self.vlo)/2+(mp.sin(2*self.vhi)-mp.sin(2*self.vlo))/4))
            self.sigma = 1 if area > 0 else -1
            self.zband = sorted((ends[0]['zm'], ends[1]['zm']))
            assert self.zband[1] > self.zband[0]
            for k, e in enumerate(ends):
                # The material lies toward the other end's height.
                e['kside'] = 1 if ends[1-k]['zm'] > e['zm'] else -1
            self.sectors = [[]]
        else:
            assert low == 0.0, 'S3: a wedge starts its tube at v = 0'
            cu, su = F(cos_rn(angle)), F(sin_rn(angle))
            self.eu = (cu, su)
            Um = mp.atan2(M(su), M(cu))
            if Um <= 0:
                Um += tau()
            self.Um = Um
            self.sigma = 1
            lines0 = (mp.mpf(0), mp.mpf(1), mp.mpf(0))
            linesU = (M(su), -M(cu), mp.mpf(0))
            if Um <= mp.pi:
                self.sectors = [[lines0, linesU]]
            else:
                self.sectors = [[lines0], [(mp.mpf(0), mp.mpf(-1), mp.mpf(0)), linesU]]
            # The end half-planes: (c, s) of the direction, K's side of the
            # plane `-s u + c v` (positive past theta = 0, negative before
            # theta_U).
            self.end_planes = [((F(1), F(0)), 1), ((cu, su), -1)]

    # -- the chart (the whole torus's)

    def chart(self, X):
        return self.T.chart(X)

    def to_world_moments(self, V, m):
        return self.T.to_world_moments(V, m)

    def face_factor(self, a):
        return self.T.face_factor(a)

    # -- the model

    def in_turn(self, theta):
        return 0 < theta % tau() < self.Um

    def in_sector(self, u, v):
        cu, su = self.eu
        a, b = v > 0, u*M(su)-v*M(cu) > 0
        return (a and b) if self.Um <= mp.pi else (a or b)

    def sides_at(self, w):
        """The arc's points at height `w`: 1 on the tube's outer side (`R +
        q`), -1 on its inner side."""
        if abs(w) >= self.rm:
            return ()
        a = mp.asin(w/self.rm)
        out = []
        if in_range(a, self.vlo, self.vhi):
            out.append(1)
        if in_range(mp.pi-a, self.vlo, self.vhi):
            out.append(-1)
        return tuple(out)

    def wall_sides(self, w):
        if abs(w) >= self.rm:
            return ()
        return self.sides_at(w) if self.kind == 'segment' else (1, -1)

    def rho_intervals(self, w, sides=None, at=None):
        """`K`'s section at height `w` as intervals of the distance from the
        axis (a wedge's: in its sector); `sides` given, the radii at `at`."""
        at = w if at is None else at
        if abs(at) > self.rm or (sides is None and abs(at) == self.rm):
            return []
        q = mp.sqrt(self.rm**2-at*at)
        if self.kind == 'wedge':
            return [(self.Rm-q, self.Rm+q)] if abs(w) < self.rm else []
        s = self.sides_at(w) if sides is None else sides
        inside = self.zband[0] < w < self.zband[1]
        if sides is None:
            assert (len(s) == 1) == inside, ('the arc crossed an odd number of times outside its end heights', w, s)
        if len(s) == 1:
            return [(mp.mpf(0), self.Rm+s[0]*q)]
        if len(s) == 2:
            return [(self.Rm-q, self.Rm+q)]
        return []

    def contains(self, p):
        """Strictly inside (mpf chart point)."""
        t = mp.sqrt(p[0]**2+p[1]**2)
        if self.kind == 'wedge':
            return self.T.contains(p) and self.in_sector(p[0], p[1])
        return any(lo < t < hi for lo, hi in self.rho_intervals(p[2]))

    def in_M(self, X):
        """`(t, w)` strictly inside the meridian section (a wedge's: the tube's
        disc)."""
        t, w = X
        if self.kind == 'wedge':
            return (t-self.Rm)**2+w*w < self.rm**2
        return any(lo < t < hi for lo, hi in self.rho_intervals(w))

    def turn_pieces(self, t0, t1):
        """The parts of the angle interval `[t0, t1]` in the turn."""
        if self.kind == 'segment':
            return [(t0, t1)]
        out = []
        T = tau()
        for k in range(-1, 3):
            lo, hi = max(t0, k*T), min(t1, k*T+self.Um)
            if hi > lo:
                out.append((lo, hi))
        return out

    def end_units(self):
        """A wedge's end directions (mpf unit vectors)."""
        out = []
        for (c, s), _ in self.end_planes:
            nn = mp.sqrt(M(c*c+s*s))
            out.append((M(c)/nn, M(s)/nn))
        return out

    # -- closed forms

    def wall_area(self):
        T = self.T
        R, r = self.Rm, self.rm
        if self.kind == 'segment':
            if T.orthonormal:
                return T.N0*2*mp.pi*r*(R*(self.vhi-self.vlo)+r*(mp.sin(self.vhi)-mp.sin(self.vlo)))
            m = 48
            h = 2*mp.pi/m
            return sum(T.meridian_weight(i*h, self.vlo, self.vhi) for i in range(m))*h
        if T.orthonormal:
            return T.N0*self.Um*2*mp.pi*R*r
        return T._gl(lambda th: T.meridian_weight(th, mp.mpf(0), 2*mp.pi), mp.mpf(0), self.Um)

    def end_areas(self):
        if self.kind == 'segment':
            f = self.face_factor((F(0), F(0), F(1)))
            return [mp.pi*e['rho']**2*f for e in self.ends]
        return [mp.pi*self.rm**2*self.face_factor((-s, c, F(0))) for (c, s), _ in self.end_planes]

    def closed(self):
        """Volume, first moments and area (world) in closed form: Pappus's
        theorems on the meridian section."""
        R, r = self.Rm, self.rm
        if self.kind == 'segment':
            g = green_mer_arc(R, r, self.vlo, self.vhi)
            V = self.sigma*2*mp.pi*g[0]
            mom = (mp.mpf(0), mp.mpf(0), self.sigma*2*mp.pi*g[2])
        else:
            V = self.Um*mp.pi*r*r*R
            I2 = mp.pi*r*r*(R*R+r*r/4)
            mom = (I2*mp.sin(self.Um), I2*(1-mp.cos(self.Um)), mp.mpf(0))
        Vw, momw = self.to_world_moments(V, mom)
        return Vw, momw, self.wall_area()+sum(self.end_areas())


# ------------------------------------------------------------------ normal slices

def chord_in_lines(P1, P2, lines):
    """The parameter interval of the segment from `P1` to `P2` where every
    `al x + be y >= ga` holds, or None."""
    lo, hi = mp.mpf(0), mp.mpf(1)
    for al, be, ga in lines:
        f0 = al*P1[0]+be*P1[1]-ga
        f1 = al*P2[0]+be*P2[1]-ga
        if f0 >= 0 and f1 >= 0:
            continue
        if f0 < 0 and f1 < 0:
            return None
        t = f0/(f0-f1)
        if f0 < 0:
            lo = max(lo, t)
        else:
            hi = min(hi, t)
    return (lo, hi) if hi > lo else None


def chord_in_disc(P1, P2, rho2):
    """The parameter interval of the segment inside the disc of squared radius
    `rho2` about the origin, or None."""
    ts = seg_circle_params(P1, P2, (mp.mpf(0), mp.mpf(0)), rho2)
    if len(ts) != 2:
        return None
    lo, hi = max(ts[0], mp.mpf(0)), min(ts[1], mp.mpf(1))
    return (lo, hi) if hi > lo else None


def chord_inside(K, P1, P2, s, skip=None):
    """The parameter length of the chord inside `K`'s section at height `s`
    (its intervals, in its sector); `skip` a parameter interval left out."""
    total = mp.mpf(0)
    ivs = K.rho_intervals(s)
    for sector in K.sectors:
        sec = chord_in_lines(P1, P2, sector)
        if sec is None:
            continue
        pieces = [sec]
        if skip is not None:
            pieces = [(a, b) for a, b in iv_diff([sec], [skip])]
        for lo, hi in ivs:
            for a, b in pieces:
                for rad, sign in ((hi, 1), (lo, -1)):
                    if rad <= 0:
                        continue
                    d = chord_in_disc(P1, P2, rad*rad)
                    if d is not None:
                        total += sign*max(mp.mpf(0), min(d[1], b)-max(d[0], a))
    return total


def poly_green(poly):
    out = zero3()
    for i in range(len(poly)):
        out = vadd(out, green_seg(poly[i], poly[(i+1) % len(poly)]))
    return out


def seg_in_convex(A, B, poly):
    """The parameter length of the segment from `A` to `B` inside a convex
    counter-clockwise polygon."""
    lines = [line_of(poly[i], poly[(i+1) % len(poly)]) for i in range(len(poly))]
    iv = chord_in_lines(A, B, lines)
    return mp.mpf(0) if iv is None else iv[1]-iv[0]


class Normal:
    """Both solids sliced by the chart planes `w = s`."""

    def __init__(self, K, P, size):
        self.K, self.T, self.P, self.size = K, K.T, P, size
        self.pc = {key: K.chart(X) for key, X in P.P.items()}
        self.pm = {key: Mv(p) for key, p in self.pc.items()}
        self.level = {key: p[2] for key, p in self.pc.items()}
        self.quad_error = mp.mpf(0)
        self.breaks = None
        self._res = None

    sections = whole.Normal.sections

    def integrand(self, s):
        K, T = self.K, self.T
        polys = [p for p in self.sections(s) if p]
        ivs = K.rho_intervals(s)
        if K.kind == 'segment':
            aK = sum(mp.pi*(b*b-a*a) for a, b in ivs)
            muK = mvK = mp.mpf(0)
        elif ivs:
            (a, b), = ivs
            aK = K.Um/2*(b*b-a*a)
            c3 = (b**3-a**3)/3
            muK, mvK = c3*mp.sin(K.Um), c3*(1-mp.cos(K.Um))
        else:
            aK = muK = mvK = mp.mpf(0)
        gP = zero3()
        for poly in polys:
            gP = vadd(gP, poly_green(poly))
        common = zero3()
        for sector in K.sectors:
            for poly in polys:
                cp = poly
                for ln in sector:
                    cp = clip_poly(cp, ln)
                    if not cp:
                        break
                if not cp:
                    continue
                for lo, hi in ivs:
                    common = vadd(common, convex_disc(cp, Disc2((mp.mpf(0), mp.mpf(0)), hi*hi)))
                    if lo > 0:
                        common = vadd(common, convex_disc(cp, Disc2((mp.mpf(0), mp.mpf(0)), lo*lo)), -1)
        out = [aK, s*aK, muK, mvK, gP[0], s*gP[0], gP[1], gP[2], common[0], s*common[0], common[1], common[2]]
        w_in = w_all = mp.mpf(0)
        if abs(s) < T.rm:
            q = mp.sqrt(T.rm**2-s*s)
            P2 = Polys2(polys)
            for side in K.wall_sides(s):
                rho = T.Rm+side*q
                phi = mp.atan2(s, side*q)
                arcs = classify(Disc2((mp.mpf(0), mp.mpf(0)), rho*rho), P2)[4]
                for t0, t1, inside in arcs:
                    for a, b in K.turn_pieces(t0, t1):
                        wgt = T.rm*rho/q*T.parallel_weight(phi, a, b)
                        w_all += wgt
                        if inside:
                            w_in += wgt
        out += [w_in, w_all]
        if K.kind == 'wedge':
            for e in K.end_units():
                if abs(s) < T.rm:
                    q = mp.sqrt(T.rm**2-s*s)
                    A = ((T.Rm-q)*e[0], (T.Rm-q)*e[1])
                    B = ((T.Rm+q)*e[0], (T.Rm+q)*e[1])
                    inside = sum((seg_in_convex(A, B, poly) for poly in polys), mp.mpf(0))*2*q
                    out += [inside, 2*q]
                else:
                    out += [mp.mpf(0), mp.mpf(0)]
        return out

    def breakpoints(self):
        K, T, P = self.K, self.T, self.P
        levels = set(self.level.values())
        levels |= {-T.r, T.r}
        polys = []
        edges = set(e for piece in P.pieces for e in piece.edges)
        for i, j in edges:
            vi, vj = self.pc[i], self.pc[j]
            e = sub(vj, vi)
            if e[2] == 0:
                continue
            quart = T.line_quartic(vi, e)
            lin = [-vi[2]/e[2], 1/e[2]]
            poly, powk = [F(0)], [F(1)]
            for c in quart:
                poly = padd(poly, pscale(powk, c))
                powk = pmul(powk, lin)
            polys.append(poly)
        faces = []
        for piece in P.pieces:
            for f in piece.faces:
                p0, p1, p2 = (self.pc[k] for k in f['loop'][:3])
                a = cross(sub(p1, p0), sub(p2, p0))
                faces.append((a, dot(a, p0)))
                if a[0] == 0 and a[1] == 0:
                    continue
                polys.append(circle_tangency_quartic(T, a, dot(a, p0)))
        pts = [M(v) for v in levels]
        if K.kind == 'segment':
            pts += [e['zm'] for e in K.ends]
        else:
            pts += wedge_levels(K, self.pc, edges, faces)
        for p in polys:
            pts += real_roots(p)
        lo = min([M(v) for v in self.level.values()]+[-T.rm])
        hi = max([M(v) for v in self.level.values()]+[T.rm])
        pts = sorted(p for p in pts if lo <= p <= hi)
        out = []
        gap = mp.mpf(10)**-30*self.size
        for p in [lo]+pts+[hi]:
            if not out or p-out[-1] > gap:
                out.append(p)
        self.raw_breaks = pts
        return out

    def measure(self):
        if self.breaks is None:
            self.breaks = self.breakpoints()
        tol = mp.mpf(10)**-33*self.size**4
        total = None
        for a, b in zip(self.breaks, self.breaks[1:]):
            est, diff = integrate(self.integrand, a, b, tol)
            if est is None:
                continue
            self.quad_error = max(self.quad_error, diff)
            total = est if total is None else [x+y for x, y in zip(total, est)]
        return total

    def results(self):
        """{name: (world volume, world moments)} for `K`, `P`, `common` (by
        clipping), `fuse`, `K-P`, `P-K` (inclusion and exclusion); the wall
        inside the prism (`wall`) and a wedge's end discs' chords (`ends`,
        chart areas inside and in all)."""
        if self._res is not None:
            return self._res
        v = self.measure()
        G = {name: v[4*k:4*k+4] for k, name in enumerate(('K', 'P', 'common'))}
        combos = {'K': [('K', 1)], 'P': [('P', 1)], 'common': [('common', 1)],
                  'fuse': [('K', 1), ('P', 1), ('common', -1)], 'K-P': [('K', 1), ('common', -1)],
                  'P-K': [('P', 1), ('common', -1)]}
        out = {}
        for name, parts in combos.items():
            A, sA, ma, mb = [sum(k*G[p][i] for p, k in parts) for i in range(4)]
            out[name] = self.K.to_world_moments(A, (ma, mb, sA))
        out['wall'] = {'in': v[12], 'out': v[13]-v[12], 'total': v[13]}
        out['ends'] = [(v[14+2*k], v[15+2*k]) for k in range((len(v)-14)//2)]
        self._res = out
        return out


def wedge_levels(K, pc, edges, faces):
    """Heights where a wedge's slice changes against the prism: the prism's
    edges crossing an end's plane, the faces' planes meeting an end's circle
    or the axis."""
    T = K.T
    out = []
    units = K.end_units()
    for ((c, s), _), e in zip(K.end_planes, units):
        sig = {k: -s*p[0]+c*p[1] for k, p in pc.items()}
        for i, j in edges:
            if (sig[i] < 0 < sig[j]) or (sig[j] < 0 < sig[i]):
                t = sig[i]/(sig[i]-sig[j])
                out.append(M(pc[i][2]+t*(pc[j][2]-pc[i][2])))
        for a, b in faces:
            am, bm = Mv(a), M(b)
            alpha = am[0]*e[0]+am[1]*e[1]
            if abs(alpha) > mp.mpf(10)**-30*(abs(am[0])+abs(am[1])+abs(am[2])):
                k, m = -am[2]/alpha, bm/alpha-T.Rm
                out += quad_roots_mp(k*k+1, 2*k*m, m*m-T.rm**2)
            elif am[2] != 0:
                out.append(bm/am[2])
    for a, b in faces:
        if a[2] != 0:
            out.append(M(b/a[2]))
    return out


# ------------------------------------------------------------------ meridian half-planes

class Meridian(whole.Meridian):
    """Both solids sliced by the half-planes `theta = const` about the axis:
    S9d.4a's sections of the prism against `M` (see the module's note)."""

    def __init__(self, K, P, size):
        super().__init__(K.T, P, size)
        self.K = K

    def plines(self, polys):
        out = []
        for poly in polys:
            n = len(poly)
            out.append([line_of(poly[k][0], poly[(k+1) % n][0]) for k in range(n)])
        return out

    def integrand(self, theta):
        K = self.K
        if K.kind == 'wedge':
            if K.in_turn(theta):
                return super().integrand(theta)
            polys = [p for p in self.sections(theta) if p]
            c, s = mp.cos(theta), mp.sin(theta)
            g = [mp.mpf(0)]*3
            for poly in polys:
                n = len(poly)
                for k in range(n):
                    g = [x+y for x, y in zip(g, green_mer_seg(poly[k][0], poly[(k+1) % n][0]))]
            return [mp.mpf(0)]*12+[g[0], c*g[1], s*g[1], g[2]]+[mp.mpf(0), mp.mpf(0)]
        return self.segment_integrand(theta)

    def segment_integrand(self, theta):
        K, T = self.K, self.T
        R, r = T.Rm, T.rm
        C = (R, mp.mpf(0))
        polys = [p for p in self.sections(theta) if p]
        pls = self.plines(polys)
        in_poly = lambda X: any(all(al*X[0]+be*X[1] > ga for al, be, ga in ls) for ls in pls)
        acc = {k: [mp.mpf(0)]*3 for k in ('din', 'dout', 'pin', 'pout')}
        w_in = w_all = mp.mpf(0)
        Tt = tau()
        angs = [K.vlo, K.vhi]
        for ls in pls:
            for ln in ls:
                for t in circle_line_angles(C, r, ln):
                    t = t+(mp.floor((K.vlo-t)/Tt)+1)*Tt
                    while t < K.vhi:
                        if t > K.vlo:
                            angs.append(t)
                        t += Tt
        angs.sort()
        for a, b in zip(angs, angs[1:]):
            if b-a <= 0:
                continue
            m = (a+b)/2
            inside = in_poly((R+r*mp.cos(m), r*mp.sin(m)))
            g = [K.sigma*x for x in green_mer_arc(R, r, a, b)]
            key = 'din' if inside else 'dout'
            acc[key] = [x+y for x, y in zip(acc[key], g)]
            wgt = T.meridian_weight(theta, a, b)
            w_all += wgt
            if inside:
                w_in += wgt
        zs = [e['zm'] for e in K.ends]
        for poly in polys:
            n = len(poly)
            for k in range(n):
                p, q = poly[k][0], poly[(k+1) % n][0]
                ts = [mp.mpf(0), mp.mpf(1)]
                ts += [t for t in seg_circle_params(p, q, C, r*r) if 0 < t < 1]
                for z in zs:
                    if (p[1] < z < q[1]) or (q[1] < z < p[1]):
                        ts.append((z-p[1])/(q[1]-p[1]))
                ts.sort()
                for t0, t1 in zip(ts, ts[1:]):
                    if t1-t0 <= 0:
                        continue
                    a = (p[0]+t0*(q[0]-p[0]), p[1]+t0*(q[1]-p[1]))
                    b = (p[0]+t1*(q[0]-p[0]), p[1]+t1*(q[1]-p[1]))
                    key = 'pin' if K.in_M(((a[0]+b[0])/2, (a[1]+b[1])/2)) else 'pout'
                    acc[key] = [x+y for x, y in zip(acc[key], green_mer_seg(a, b))]
        c, s = mp.cos(theta), mp.sin(theta)
        out = []
        for key in ('din', 'dout', 'pin', 'pout'):
            m1, m2, mtw = acc[key]
            out += [m1, c*m2, s*m2, mtw]
        out += [w_in, w_all]
        for e in K.ends:
            out += self.end_radial(e, polys, pls)
        return out

    def end_radial(self, e, polys, pls):
        """The end disc's radius `[0, rho_end]` at its height in this section:
        `t dt` inside the polygon, on a coplanar edge (`same`, `opp`), in all."""
        z, rho = e['zm'], e['rho']
        eps = mp.mpf(10)**-30*self.size
        ts = [mp.mpf(0), rho]
        flats = []
        for poly in polys:
            n = len(poly)
            cw = sum(p[0][1] for p in poly)/n
            for k in range(n):
                p, q = poly[k][0], poly[(k+1) % n][0]
                if abs(p[1]-z) < eps and abs(q[1]-z) < eps:
                    flats.append((min(p[0], q[0]), max(p[0], q[0]), 1 if cw > z else -1))
                    ts += [p[0], q[0]]
                elif (p[1] < z < q[1]) or (q[1] < z < p[1]):
                    ts.append(p[0]+(z-p[1])/(q[1]-p[1])*(q[0]-p[0]))
        ts = sorted(t for t in ts if 0 <= t <= rho)
        cls = {'in': mp.mpf(0), 'same': mp.mpf(0), 'opp': mp.mpf(0)}
        for t0, t1 in zip(ts, ts[1:]):
            if t1-t0 <= 0:
                continue
            tm = (t0+t1)/2
            wgt = (t1*t1-t0*t0)/2
            flat = [side for lo, hi, side in flats if lo < tm < hi]
            if flat:
                cls['same' if flat[0] == e['kside'] else 'opp'] += wgt
            elif any(all(al*tm+be*z > ga for al, be, ga in ls) for ls in pls):
                cls['in'] += wgt
        return [cls['in'], cls['same'], cls['opp'], rho*rho/2]

    def breakpoints(self):
        base = super().breakpoints()
        K, T = self.K, self.T
        pi = mp.pi
        angs = list(self.raw_breaks)
        if K.kind == 'segment':
            for e in K.ends:
                z, rho = e['z'], e['rho']
                for f in self.faces:
                    au, av, aw = f['am']
                    nuv = mp.sqrt(au*au+av*av)
                    if nuv == 0:
                        continue
                    val = (f['bm']-aw*M(z))/(rho*nuv)
                    if abs(val) <= 1:
                        phi, w = mp.atan2(av, au), mp.acos(val)
                        angs += [phi-w, phi+w]
                edges = set(ed for piece in self.P.pieces for ed in piece.edges)
                for i, j in edges:
                    vi, vj = self.pc[i], self.pc[j]
                    if (vi[2] < z < vj[2]) or (vj[2] < z < vi[2]):
                        t = (z-vi[2])/(vj[2]-vi[2])
                        X = add(vi, scale(sub(vj, vi), t))
                        if X[0] != 0 or X[1] != 0:
                            angs.append(mp.atan2(M(X[1]), M(X[0])))
        else:
            angs += [mp.mpf(0), K.Um]
        T2 = 2*pi
        norm = sorted(((a+pi) % T2)-pi for a in angs)
        out = []
        gap = mp.mpf(10)**-30
        for p in [-pi]+norm+[pi]:
            if not out or p-out[-1] > gap:
                out.append(p)
        if pi-out[-1] <= gap:
            out[-1] = pi
        else:
            out.append(pi)
        self.raw_breaks = norm
        return out

    def results(self):
        """As `Normal.results`, every operation apart: the volumes and
        moments by the cylindrical element (`M`'s arc signed by its
        orientation), the wall by its own parameters, a segment's end discs
        by their radii (`ends`: chart areas in, same, opp, all)."""
        if self._res is not None:
            return self._res
        v = self.measure()
        G = {name: v[4*k:4*k+4] for k, name in enumerate(('din', 'dout', 'pin', 'pout'))}
        combos = {'K': [('din', 1), ('dout', 1)], 'P': [('pin', 1), ('pout', 1)],
                  'common': [('din', 1), ('pin', 1)], 'fuse': [('dout', 1), ('pout', 1)],
                  'K-P': [('dout', 1), ('pin', -1)], 'P-K': [('pout', 1), ('din', -1)]}
        out = {}
        for name, parts in combos.items():
            V, mu, mv, mw = [sum(k*G[p][i] for p, k in parts) for i in range(4)]
            out[name] = self.K.to_world_moments(V, (mu, mv, mw))
        out['wall'] = {'in': v[16], 'out': v[17]-v[16], 'total': v[17]}
        out['ends'] = [tuple(v[18+4*k:22+4*k]) for k in range((len(v)-18)//4)]
        self._res = out
        return out


# ------------------------------------------------------------------ end faces, exactly

def segment_end_classes(K, prism):
    """[(end, {in, out, same, opp}, area)] for a segment's end discs: each
    piece's section at the end's height, or its face on the plane."""
    out = []
    factor = K.face_factor((F(0), F(0), F(1)))
    for k, e in enumerate(K.ends):
        z = e['z']
        D = Disc2((mp.mpf(0), mp.mpf(0)), e['rho']**2)
        area = mp.pi*e['rho']**2*factor
        groups = {'in': [], 'same': [], 'opp': []}
        for piece in prism.pieces:
            ch = {key: K.chart(prism.P[key]) for key in piece.keys}
            vals = {key: p[2]-z for key, p in ch.items()}
            pos = any(v > 0 for v in vals.values())
            neg = any(v < 0 for v in vals.values())
            zeros = [key for key, v in vals.items() if v == 0]
            if pos and neg:
                pts = [(M(ch[key][0]), M(ch[key][1])) for key in zeros]
                for i, j in piece.edges:
                    if (vals[i] < 0 < vals[j]) or (vals[j] < 0 < vals[i]):
                        t = vals[i]/(vals[i]-vals[j])
                        X = add(ch[i], scale(sub(ch[j], ch[i]), t))
                        pts.append((M(X[0]), M(X[1])))
                poly = hull2(pts)
                if poly:
                    groups['in'].append(poly)
            elif len(zeros) >= 3:
                poly = hull2([(M(ch[key][0]), M(ch[key][1])) for key in zeros])
                side = 1 if pos else -1
                groups['same' if side == e['kside'] else 'opp'].append(poly)
        cls = classes_in_disc(D, groups, area, factor)
        out.append((('end', k), cls, area))
    return out


def classes_in_disc(D, groups, area, factor):
    cls = {'in': mp.mpf(0), 'out': mp.mpf(0), 'same': mp.mpf(0), 'opp': mp.mpf(0)}
    allp = groups['in']+groups['same']+groups['opp']
    for key, polys in groups.items():
        for p in polys:
            cls[key] += convex_disc(p, D)[0]*factor
    if allp:
        din, dout, pin, pout, _ = classify(D, Polys2(allp))
        cls['out'] = (dout[0]-pin[0])*factor
    else:
        cls['out'] = area
    return cls


def wedge_end_classes(K, prism):
    """[(end, {in, out, same, opp}, area)] for a wedge's end discs: each
    piece's section by the end's half-plane in its own coordinates `(tau, w)`
    (`tau` along the end's unit direction), clipped at the axis, or its face
    on the plane."""
    out = []
    R, r = K.Rm, K.rm
    for k, (((c, s), kside), e) in enumerate(zip(K.end_planes, K.end_units())):
        factor = K.face_factor((-s, c, F(0)))
        D = Disc2((R, mp.mpf(0)), r*r)
        area = mp.pi*r*r*factor
        groups = {'in': [], 'same': [], 'opp': []}
        to2 = lambda X: (M(X[0])*e[0]+M(X[1])*e[1], M(X[2]))
        axis = (mp.mpf(1), mp.mpf(0), mp.mpf(0))
        for piece in prism.pieces:
            ch = {key: K.chart(prism.P[key]) for key in piece.keys}
            vals = {key: -s*p[0]+c*p[1] for key, p in ch.items()}
            pos = any(v > 0 for v in vals.values())
            neg = any(v < 0 for v in vals.values())
            zeros = [key for key, v in vals.items() if v == 0]
            if pos and neg:
                pts = [to2(ch[key]) for key in zeros]
                for i, j in piece.edges:
                    if (vals[i] < 0 < vals[j]) or (vals[j] < 0 < vals[i]):
                        t = vals[i]/(vals[i]-vals[j])
                        pts.append(to2(add(ch[i], scale(sub(ch[j], ch[i]), t))))
                poly = hull2(pts)
                poly = clip_poly(poly, axis) if poly else []
                if poly:
                    groups['in'].append(poly)
            elif len(zeros) >= 3:
                poly = hull2([to2(ch[key]) for key in zeros])
                poly = clip_poly(poly, axis) if poly else []
                if poly:
                    side = 1 if pos else -1
                    groups['same' if side == kside else 'opp'].append(poly)
        out.append((('end', k), classes_in_disc(D, groups, area, factor), area))
    return out


# ------------------------------------------------------------------ the prism's faces

def slanted_face(K, polys, a, b, size, coplanar=None):
    """(inside, coplanar, all) chart areas of convex chart polygons in the
    plane `a . p = b` (not normal to the axis), sliced by the lines of
    constant `w`; `coplanar` a wedge end's unit direction when the plane
    holds its half-plane (the chord on the end's ray inside the annulus is
    coplanar with the end disc)."""
    T = K.T
    au, av, aw = a
    nuv2 = au*au+av*av
    spacing = mp.sqrt(M(dot(a, a)))/mp.sqrt(M(nuv2))
    levels = {-T.r, T.r}
    polys_ = [circle_tangency_quartic(T, a, b)]
    edges, pts3 = [], []
    for poly in polys:
        for i in range(len(poly)):
            p, q = poly[i], poly[(i+1) % len(poly)]
            levels.add(p[2])
            edges.append((p, q))
            e = sub(q, p)
            if e[2] != 0:
                quart = T.line_quartic(p, e)
                lin = [-p[2]/e[2], 1/e[2]]
                acc, powk = [F(0)], [F(1)]
                for c in quart:
                    acc = padd(acc, pscale(powk, c))
                    powk = pmul(powk, lin)
                polys_.append(acc)
    lo = min(p[2] for poly in polys for p in poly)
    hi = max(p[2] for poly in polys for p in poly)
    pts = [M(v) for v in levels if lo <= v <= hi]
    for c in polys_:
        pts += [t for t in real_roots(c) if M(lo) <= t <= M(hi)]
    if K.kind == 'segment':
        pts += [e['zm'] for e in K.ends if M(lo) <= e['zm'] <= M(hi)]
    else:
        keyed = {('v', i): p for i, p in enumerate(p for poly in polys for p in poly)}
        idx = {}
        for (_, i), p in keyed.items():
            idx.setdefault(p, ('v', i))
        eds = [(idx[p], idx[q]) for p, q in edges]
        pc = {k: p for k, p in keyed.items()}
        pts += [t for t in wedge_levels(K, pc, eds, [(a, b)]) if M(lo) <= t <= M(hi)]
    pts.sort()
    brk = []
    for p in pts:
        if not brk or p-brk[-1] > mp.mpf(10)**-30*size:
            brk.append(p)
    pm = [[Mv(p) for p in poly] for poly in polys]
    Rm, rm = T.Rm, T.rm

    def f(s):
        Lin = Lco = Lall = mp.mpf(0)
        for poly in pm:
            ends = []
            for i in range(len(poly)):
                p, q = poly[i], poly[(i+1) % len(poly)]
                if (p[2] < s < q[2]) or (q[2] < s < p[2]):
                    t = (s-p[2])/(q[2]-p[2])
                    ends.append((p[0]+t*(q[0]-p[0]), p[1]+t*(q[1]-p[1])))
            if len(ends) != 2:
                continue
            P1, P2 = ends
            L = mp.sqrt((P2[0]-P1[0])**2+(P2[1]-P1[1])**2)
            Lall += L
            skip = None
            if coplanar is not None:
                # The chord's part on the end's ray (a positive component
                # along its direction), inside the annulus.
                d1 = P1[0]*coplanar[0]+P1[1]*coplanar[1]
                d2 = P2[0]*coplanar[0]+P2[1]*coplanar[1]
                ray = chord_in_lines(P1, P2, [(coplanar[0], coplanar[1], mp.mpf(0))]) if (d1 > 0 or d2 > 0) \
                    else None
                if ray is not None and abs(s) < rm:
                    skip = ray
                    q = mp.sqrt(rm*rm-s*s)
                    for rad, sign in ((Rm+q, 1), (Rm-q, -1)):
                        d = chord_in_disc(P1, P2, rad*rad)
                        if d is not None:
                            Lco += sign*L*max(mp.mpf(0), min(d[1], ray[1])-max(d[0], ray[0]))
                elif ray is not None:
                    skip = ray
            Lin += L*chord_inside(K, P1, P2, s, skip)
        return [Lin, Lco, Lall]
    tol = mp.mpf(10)**-33*size**3
    tin = tco = tall = mp.mpf(0)
    for x0, x1 in zip(brk, brk[1:]):
        est, _ = integrate(f, x0, x1, tol)
        if est is not None:
            tin += est[0]
            tco += est[1]
            tall += est[2]
    return spacing*tin, spacing*tco, spacing*tall


def horizontal_regions(K, w0):
    """`K`'s regions at the height `w0` as intervals of the distance from the
    axis: (both sides, above only, below only)."""
    if abs(w0) > K.rm:
        return [], [], []
    d = mp.mpf(10)**-35
    if K.kind == 'segment':
        above = K.rho_intervals(w0+d, K.sides_at(w0+d), w0)
        below = K.rho_intervals(w0-d, K.sides_at(w0-d), w0)
    else:
        above = below = K.rho_intervals(w0) if abs(w0) < K.rm else []
    return iv_intersect(above, below), iv_diff(above, below), iv_diff(below, above)


def area_in_annuli(K, p2d, ivs):
    """The chart area of convex polygons inside the intervals revolved (in
    the sector)."""
    total = mp.mpf(0)
    for sector in K.sectors:
        for poly in p2d:
            cp = poly
            for ln in sector:
                cp = clip_poly(cp, ln)
                if not cp:
                    break
            if not cp:
                continue
            for lo, hi in ivs:
                total += convex_disc(cp, Disc2((mp.mpf(0), mp.mpf(0)), hi*hi))[0]
                if lo > 0:
                    total -= convex_disc(cp, Disc2((mp.mpf(0), mp.mpf(0)), lo*lo))[0]
    return total


def face_classes(K, polys, outward, size):
    """{in, out, same, opp} (world areas) of a planar face of the prism
    (convex chart polygons, its chart outward normal) against the part."""
    p0, p1, p2 = polys[0][0], polys[0][1], polys[0][2]
    a = cross(sub(p1, p0), sub(p2, p0))
    if dot(a, outward) < 0:
        a = scale(a, -1)
    b = dot(a, p0)
    factor = K.face_factor(a)
    cls = {'in': mp.mpf(0), 'out': mp.mpf(0), 'same': mp.mpf(0), 'opp': mp.mpf(0)}
    if a[0] == 0 and a[1] == 0:
        w0 = M(b/a[2])
        p2d = [ccw2([(M(p[0]), M(p[1])) for p in poly]) for poly in polys]
        total = sum(abs(sum(cross2(p[i-1], p[i]) for i in range(len(p)))/2) for p in p2d)
        both, above, below = horizontal_regions(K, w0)
        cls['in'] = area_in_annuli(K, p2d, both)*factor
        # The prism lies below its face when the outward normal points up;
        # `K` only below (above) the plane: its end disc faces up (down).
        prism_below = a[2] > 0
        k_below = area_in_annuli(K, p2d, below)*factor
        k_above = area_in_annuli(K, p2d, above)*factor
        cls['same'] = k_below if prism_below else k_above
        cls['opp'] = k_above if prism_below else k_below
        cls['out'] = total*factor-cls['in']-cls['same']-cls['opp']
        return cls
    coplanar = None
    if K.kind == 'wedge' and a[2] == 0 and b == 0:
        for ((c, s), kside), e in zip(K.end_planes, K.end_units()):
            if a[0]*c+a[1]*s == 0:
                coplanar = (e, kside)
    inside, co, total = slanted_face(K, polys, a, b, size, coplanar[0] if coplanar else None)
    cls['in'] = inside*factor
    if coplanar is not None:
        (e, kside) = coplanar
        # The end disc's outward normal: -kside (-s, c); the prism's `a`.
        (c, s), _ = next(pl for pl in K.end_planes if pl[1] == kside)
        same = (a[0]*(-s)+a[1]*c)*(-kside) > 0
        cls['same' if same else 'opp'] = co*factor
    cls['out'] = (total-inside-co)*factor
    return cls


def prism_face_classes(K, prism, size):
    """[(tag, {in, out, same, opp}, exact area)] for the prism's faces."""
    out = []
    for tag, polys, piece, area in prism.faces():
        cp = [[K.chart(X) for X in poly] for poly in polys]
        p0, p1, p2 = cp[0][0], cp[0][1], cp[0][2]
        a = cross(sub(p1, p0), sub(p2, p0))
        if dot(a, K.chart(piece.centroid)) > dot(a, p0):
            a = scale(a, -1)
        out.append((tag, face_classes(K, cp, a, size), area))
    return out


# ------------------------------------------------------------------ solids by sweeping

class Section:
    """One operation's region in the meridian half-plane at `theta`, in bands
    of `w`: each band's intervals of `t`, joined into components."""

    def __init__(self, sweep, theta, op):
        self.theta = theta
        K, mer = sweep.K, sweep.mer
        polys = [p for p in mer.sections(theta) if p]
        self.polys = [[p[0] for p in poly] for poly in polys]
        tiny = sweep.tiny
        R, r = K.Rm, K.rm
        levels = [-r, r]
        if K.kind == 'segment':
            levels += [e['zm'] for e in K.ends]
        for poly in self.polys:
            n = len(poly)
            for k in range(n):
                p, q = poly[k], poly[(k+1) % n]
                levels.append(p[1])
                for t in seg_circle_params(p, q, (R, mp.mpf(0)), r*r):
                    if 0 < t < 1:
                        levels.append(p[1]+t*(q[1]-p[1]))
        levels.sort()
        levs = []
        for x in levels:
            if not levs or x-levs[-1] > mp.mpf(10)**-28*sweep.size:
                levs.append(x)
        self.levs = levs
        f = lambda w: [iv for iv in sweep.op_intervals(op, theta, self.polys, w) if iv[1] > iv[0]]
        self.f = f
        self.bands = [f((a+b)/2) for a, b in zip(levs, levs[1:])]
        parent = {(k, i): (k, i) for k, band in enumerate(self.bands) for i in range(len(band))}

        def find(x):
            while parent[x] != x:
                x = parent[x]
            return x
        for k in range(len(self.bands)-1):
            L = levs[k+1]
            eps = min(mp.mpf(10)**-30*sweep.size, (levs[k+1]-levs[k])/4, (levs[k+2]-levs[k+1])/4)
            below, above = f(L-eps), f(L+eps)
            assert len(below) == len(self.bands[k]) and len(above) == len(self.bands[k+1]), \
                ('a band changed between its levels', op, theta, L)
            for i, x in enumerate(below):
                for j, y in enumerate(above):
                    if min(x[1], y[1])-max(x[0], y[0]) > tiny:
                        parent[find((k, i))] = find((k+1, j))
        self.root = {key: find(key) for key in parent}
        self.comps = sorted(set(self.root.values()))
        self.samples = {c: [] for c in self.comps}
        self.axis = {c: [] for c in self.comps}
        for (k, i), c in self.root.items():
            lo, hi = self.bands[k][i]
            wm = (levs[k]+levs[k+1])/2
            self.samples[c].append(((lo+hi)/2, wm))
            if lo <= tiny:
                self.axis[c].append((levs[k], levs[k+1]))

    def locate(self, X):
        """The component holding `(t, w)`: the intervals at its height, in
        the order of its band's (constant between the levels)."""
        t, w = X
        k = bisect_right(self.levs, w)-1
        if k < 0 or k >= len(self.bands) or not (self.levs[k] < w < self.levs[k+1]):
            return None
        ivs = self.f(w)
        assert len(ivs) == len(self.bands[k]), ('a band changed between its levels', self.theta, w)
        for i, (lo, hi) in enumerate(ivs):
            if lo < t < hi:
                return self.root[(k, i)]
        return None


class Sweep:
    """Solid counts by sweeping the meridian half-planes (see the module's
    note)."""

    def __init__(self, pair):
        self.pair, self.K, self.mer = pair, pair.K, pair.meridian
        self.size = pair.size
        self.tiny = mp.mpf(10)**-13*self.size

    def c_intervals(self, polys, w):
        out = []
        for poly in polys:
            ts = []
            n = len(poly)
            for k in range(n):
                p, q = poly[k], poly[(k+1) % n]
                if p[1] == w:
                    ts.append(p[0])
                if (p[1] < w < q[1]) or (q[1] < w < p[1]):
                    ts.append(p[0]+(w-p[1])/(q[1]-p[1])*(q[0]-p[0]))
            if len(ts) >= 2 and max(ts) > min(ts):
                out.append((min(ts), max(ts)))
        return iv_union(out, [])

    def op_intervals(self, op, theta, polys, w):
        K = self.K
        ki = K.rho_intervals(w) if (K.kind == 'segment' or K.in_turn(theta)) else []
        ci = self.c_intervals(polys, w)
        if op == 'common':
            return iv_intersect(ki, ci)
        if op == 'K-P':
            return iv_diff(ki, ci)
        if op == 'P-K':
            return iv_diff(ci, ki)
        return iv_union(ki, ci)

    def links(self, A, B, at=None):
        pairs = set()
        for S1, S2, flip in ((A, B, False), (B, A, True)):
            for c, pts in S1.samples.items():
                for X in pts:
                    d = S2.locate(X)
                    if d is None or (at is not None and not at(X)):
                        continue
                    pairs.add((d, c) if flip else (c, d))
        return pairs

    def member(self, theta, op):
        """Membership of `(t, w)` in the operation's region at `theta` (a
        breakpoint: no bands built there)."""
        polys = [[p[0] for p in poly] for poly in self.mer.sections(theta) if poly]
        return lambda X: any(lo < X[0] < hi for lo, hi in self.op_intervals(op, theta, polys, X[1]))

    def bijective(self, A, B, pairs):
        left = {}
        right = {}
        for a, b in pairs:
            left.setdefault(a, set()).add(b)
            right.setdefault(b, set()).add(a)
        return set(left) == set(A.comps) and set(right) == set(B.comps) \
            and all(len(v) == 1 for v in left.values()) and all(len(v) == 1 for v in right.values())

    def count(self, op):
        mer = self.mer
        if mer.breaks is None:
            mer.breaks = mer.breakpoints()
        br = mer.breaks
        parent = {}

        def node(x):
            parent.setdefault(x, x)
            return x

        def find(x):
            while parent[x] != x:
                x = parent[x]
            return x

        def join(x, y):
            parent[find(x)] = find(y)
        sections = []   # per interval: [Section]
        for a, b in zip(br, br[1:]):
            d = min(mp.mpf(10)**-8, (b-a)/8)
            m = 6
            ths = [a+d]+[a+(b-a)*(1-mp.cos(mp.pi*j/m))/2 for j in range(1, m)]+[b-d]
            secs = [Section(self, th, op) for th in ths]
            secs = self.refine(secs, op, 0)
            sections.append(secs)
        key = lambda S, c: (id(S), c)
        for secs in sections:
            for S in secs:
                for c in S.comps:
                    node(key(S, c))
            for S1, S2 in zip(secs, secs[1:]):
                for c1, c2 in self.links(S1, S2):
                    join(key(S1, c1), key(S2, c2))
        n = len(sections)
        for i in range(n):
            S1, S2 = sections[i][-1], sections[(i+1) % n][0]
            for c1, c2 in self.links(S1, S2, self.member(br[i+1], op)):
                join(key(S1, c1), key(S2, c2))
        # Through the axis: components touching it along overlapping heights.
        touch = []
        for secs in sections:
            for S in secs:
                for c, rng in S.axis.items():
                    for lo, hi in rng:
                        touch.append((lo, hi, key(S, c)))
        for (a0, a1, x), (b0, b1, y) in itertools.combinations(touch, 2):
            if min(a1, b1)-max(a0, b0) > self.tiny and find(x) != find(y):
                join(x, y)
        return len({find(x) for x in parent})

    def refine(self, secs, op, depth):
        out = [secs[0]]
        for S1, S2 in zip(secs, secs[1:]):
            if not self.bijective(S1, S2, self.links(S1, S2)):
                assert depth < 7, ('components change between breakpoints', op, S1.theta, S2.theta)
                mid = Section(self, (S1.theta+S2.theta)/2, op)
                out += self.refine([S1, mid, S2], op, depth+1)[1:]
            else:
                out.append(S2)
        return out


# ------------------------------------------------------------------ the pair

def case_size(K, P):
    vals = [F(1)]
    vals += [abs(x) for p in P.P.values() for x in p]
    vals += [abs(x)+K.R+K.r for x in K.o]
    return M(max(vals))


class Pair:
    """A torus segment or wedge and a prism, in either order."""

    def __init__(self, obj, tool):
        self.part_first = obj.torus is not None
        k_case, p_case = (obj, tool) if self.part_first else (tool, obj)
        assert k_case.torus is not None and p_case.torus is None and p_case.sphere is None \
            and p_case.cone is None, 'S9d.4b.1: a torus segment or wedge and a prism'
        self.K, self.P = Part(k_case), Prism(p_case)
        self.size = case_size(self.K, self.P)
        self.normal = Normal(self.K, self.P, self.size)
        self.meridian = Meridian(self.K, self.P, self.size)
        self._faces = None
        self._solids = {}

    def sliced(self):
        return self.normal.results()

    def swept(self):
        return self.meridian.results()

    def end_classes(self):
        if self.K.kind == 'segment':
            return segment_end_classes(self.K, self.P)
        return wedge_end_classes(self.K, self.P)

    def faces(self):
        """[(input 'K' or 'P', tag, classes, exact area)]."""
        if self._faces is None:
            out = [('P', tag, cls, area) for tag, cls, area in prism_face_classes(self.K, self.P, self.size)]
            out += [('K', tag, cls, area) for tag, cls, area in self.end_classes()]
            wall = self.swept()['wall']
            out.append(('K', ('wall',), {'in': wall['in'], 'out': wall['out'], 'same': mp.mpf(0),
                                         'opp': mp.mpf(0)}, self.K.wall_area()))
            self._faces = out
        return self._faces

    def classes(self, which):
        tot = {'in': mp.mpf(0), 'out': mp.mpf(0), 'same': mp.mpf(0), 'opp': mp.mpf(0)}
        for w, _, cls, _ in self.faces():
            if w == which:
                for k in tot:
                    tot[k] += cls[k]
        return tot

    def roles(self):
        return ('K', 'P') if self.part_first else ('P', 'K')

    def area(self, op):
        a, b = self.roles()
        ca, cb = self.classes(a), self.classes(b)
        keep = {'fuse': (('out', 'same'), ('out',)), 'cut': (('out', 'opp'), ('in',)),
                'common': (('in', 'same'), ('in',))}[op]
        return sum(ca[k] for k in keep[0])+sum(cb[k] for k in keep[1])

    def volume(self, op):
        r = self.sliced()
        key = {'fuse': 'fuse', 'common': 'common', 'cut': 'K-P' if self.part_first else 'P-K'}[op]
        return r[key]

    def swept_count(self, key):
        """Components of `common`, `K-P`, `P-K` or `fuse` by the sweep."""
        if key not in self._solids:
            self._solids[key] = Sweep(self).count(key)
        return self._solids[key]

    def solids(self, op):
        if op == 'fuse':
            eps = mp.mpf(10)**-25
            opp = self.classes('K')['opp']
            return 1 if self.sliced()['common'][0] > eps*self.size**3 or opp > eps*self.size**2 else 2
        if op == 'common':
            return self.swept_count('common')
        return self.swept_count('K-P' if self.part_first else 'P-K')

    def result(self, op):
        vol, mom = self.volume(op)
        if vol <= mp.mpf(10)**-25*self.size**3:
            return 0, mp.mpf(0), mp.mpf(0), None
        n = self.solids(op)
        assert n > 0, 'a result of positive volume without solids'
        return n, vol, self.area(op), tuple(m/vol for m in mom)


def rows(obj, operation, tool, pair=None):
    """`result N volume area cx cy cz` or `empty`, and the pair (reused
    across the three operations)."""
    pair = pair or Pair(obj, tool)
    n, vol, surface, centre = pair.result(operation)
    if n == 0:
        return ['empty'], pair
    return [' '.join(['result', str(n), number(vol), number(surface)]+[number(c) for c in centre])], pair
