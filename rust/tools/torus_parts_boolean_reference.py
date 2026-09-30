#!/usr/bin/env python3
"""Independent reference for S9d.4c of REVIEW_NOTES.md: Booleans of a
sphere's cap or zone (`Solid::sphere_with` between two latitudes) against a
whole torus, and of a torus v-segment or wedge (`Solid::torus_with` other
than a whole torus, S9d.4b.1's parts) against a prism with arcs (planar and
cylindrical walls), a sphere (whole, a cap or a zone), a cone or frustum
and a whole torus, in any relative position. A cap's circles have surd radii
(and unequal axes on a turned frame's stored axis); a part's rims are
rings of its end planes' or half-planes' sections, of surd radius or on a
rounded direction.

**The inputs' exact models**, in rationals from the stored binary64 data
(`stored_axes`), as S9d.4b.2's reference (`torus_curved_boolean_reference`,
whose frames, surfaces, curve families, face sweeps, classes and
quadrature this module uses) with two more:

* **A torus part** is S9d.4b.1's (`torus_segment_boolean_reference.Part`):
  a v-segment the region between the tube's arc from `low` to `high` and
  the axis, revolved, its end planes `w = z` at the stored heights `r
  sin(latitude)` (correctly rounded, `sin_rn`, as the kernel's
  `scaled_sin` on every host the fixtures run), a point inside when an odd
  count of the arc's points at its height lie beyond it; a wedge the whole
  tube in the sector from the half-plane of `x` to that of the chart
  direction `(cos angle, sin angle)` rounded. Its surfaces are the torus's
  quartic and its end planes (a wedge's the planes of its half-planes);
  its faces the wall over its range (a segment's tube circles over the arc
  `[v_lo, v_hi]` and its parallels over that range, a wedge's tube circles
  over the turn and its parallels' arcs over it, oriented by the segment's
  `sigma`: the inner half is inside out) and its end discs (a segment's
  from the axis to the arc's end at its height, a wedge's the tube's discs
  in its half-planes), each by chords in two directions.
* **A cap or zone** is S9d.2c's `AxisSphere`: `|X - o|^2 <= R^2` with its
  end planes through `o + h n` normal to the stored axis `n` (`h = R
  sin(latitude)`), its faces the wall between the ends' latitudes about the
  unit axis `n / |n|` (meridians and parallels) and its end discs of radius
  `sqrt(R^2 - (h |n|)^2)`, by chords.

Every face of both inputs is swept by its two families of curves against the
other input's surfaces and membership, each piece classified and measured by
the divergence theorem (S9d.4b.2's `Sweep`): the results are the classes'
sums, every operation two ways. Nothing here uses the kernel, a surface/
surface intersection or an arrangement.

**Solids** by a sweep of the torus's (the part's) normal slices, as
S9d.4b.2's: in each slice the operation's intervals along rays from the
axis (the part's section its own intervals: a segment's disc `[0, c]` or
annulus, a wedge's annulus on the rays of its turn), joined where they
overlap on adjacent rays, at the axis and on the same ray of the adjacent
slice, or through a chain of intervals on that ray at levels between the
two slices (bisected: a thin region sheared between them, as the tip of a
half's ring outside a coaxial dome, is one solid, where overlaps alone, as
S9d.4b.2's sweep joins, would count its tip apart).

**Margins** (S9d.4b.2's, with both inputs' edges): the least sine between
the inputs' surfaces where they meet, the least distance from a face to a
surface of the other input it never crosses, the least sine at which any
edge (a prism's, a cone's rims, a cap's rims, a part's rims) crosses the
other input's surfaces and its least distance from them where stationary,
and the vertices' distances.

`rows(obj, operation, tool)` gives `result N volume area cx cy cz` (totals
over the N solids, world coordinates) or `empty`.
"""
import copy
from fractions import Fraction as F
import math

import mpmath as mp

import torus_curved_boolean_reference as tc
from torus_curved_boolean_reference import (
    Z, Face, Family, Frame, Surf, box_size, circle_curve, critical_points, curve_roots, disc_chords,
    distance_estimate, fpoly_add, fpoly_mul, frame_box, iv_ops, lin, line_curve, norm, plane_surf,
    sine)
from sphere_boolean_reference import M, Mv, cross, dot, scale
from spheres_boolean_reference import AxisSphere
from torus_segment_boolean_reference import Part
from torus_boolean_reference import number, OPS, TWO_PI

mp.mp.dps = 40

HALF_PI = math.pi/2


def unit(v):
    n = norm(v)
    return (v[0]/n, v[1]/n, v[2]/n)


def chart_disc(fr, centre, U, V, rad, direction, outward, name):
    """Chords of the disc `centre + a U + b V`, `a^2 + b^2 < rad^2` (a chart
    point and chart vectors orthonormal in the chart), along `direction` in
    `(U, V)`, the outer parameter across it; `outward` the face's outward
    normal as a chart covector (the family oriented by it)."""
    d = Mv(direction)
    dn = mp.sqrt(d[0]*d[0]+d[1]*d[1])
    d = (d[0]/dn, d[1]/dn)
    nperp = (-d[1], d[0])
    comb = lambda a, b: tuple(a*U[i]+b*V[i] for i in range(3))
    D = fr.vec(*comb(d[0], d[1]))
    dP = fr.vec(*comb(nperp[0], nperp[1]))
    sign = 1 if dot(cross(dP, D), fr.covector(outward)) > 0 else -1
    zero = (Z, Z, Z)

    def chord(t):
        q = mp.sqrt(max(rad*rad-t*t, Z))
        p = comb(t*nperp[0], t*nperp[1])
        return line_curve(fr.world(centre[0]+p[0], centre[1]+p[1], centre[2]+p[2]), D, dP, zero, -q, q)
    return Family(name, -rad, rad, False, chord, sign)


def world_disc(centre, U, V, rad, direction, outward, name):
    """Chords of the world disc `centre + a U + b V` (`U`, `V` orthonormal),
    oriented by the world vector `outward`."""
    d = Mv(direction)
    dn = mp.sqrt(d[0]*d[0]+d[1]*d[1])
    d = (d[0]/dn, d[1]/dn)
    nperp = (-d[1], d[0])
    comb = lambda a, b: tuple(a*U[i]+b*V[i] for i in range(3))
    D, dP = comb(d[0], d[1]), comb(nperp[0], nperp[1])
    sign = 1 if dot(cross(dP, D), outward) > 0 else -1
    zero = (Z, Z, Z)

    def chord(t):
        q = mp.sqrt(max(rad*rad-t*t, Z))
        p = comb(t*nperp[0], t*nperp[1])
        return line_curve(tuple(centre[i]+p[i] for i in range(3)), D, dP, zero, -q, q)
    return Family(name, -rad, rad, False, chord, sign)


# ------------------------------------------------------------------ a torus part

def complement_part(P):
    """The v-segment of the rest of the tube's circle: the same end planes,
    its arc from `P`'s high end to its low end a turn on (so `P` and it,
    weighted by their `sigma`, make the whole torus)."""
    assert P.kind == 'segment'
    Q = copy.copy(P)
    lo, hi = P.ends
    e0, e1 = dict(hi), dict(lo)
    e1['v'] = lo['v']+2*mp.pi
    Q.ends = [e0, e1]
    Q.vlo, Q.vhi = e0['v'], e1['v']
    R_, r_ = P.Rm, P.rm
    area = r_*(R_*(mp.sin(Q.vhi)-mp.sin(Q.vlo))
               + r_*((Q.vhi-Q.vlo)/2+(mp.sin(2*Q.vhi)-mp.sin(2*Q.vlo))/4))
    Q.sigma = 1 if area > 0 else -1
    Q.zband = sorted((e0['zm'], e1['zm']))
    for k, e in enumerate(Q.ends):
        e['kside'] = 1 if Q.ends[1-k]['zm'] > e['zm'] else -1
    return Q


class TorusPart:
    """A torus v-segment or wedge on S9d.4b.1's model (see the module's
    note)."""
    kind = 'torus'

    def __init__(self, case, part=None):
        R, r, low, high, angle = case.torus
        self.case = case
        self.fr = fr = Frame(case.frame)
        self.R, self.r = F(R), F(r)
        assert self.R > self.r > 0
        self.Rm, self.rm = M(self.R), M(self.r)
        self.part = P = part or Part(case)
        self.segment = P.kind == 'segment'
        R2, k = self.Rm**2, self.Rm**2-self.rm**2

        def f(X):
            p = fr.chart(X)
            s = p[0]*p[0]+p[1]*p[1]+p[2]*p[2]
            return (s+k)**2-4*R2*(p[0]*p[0]+p[1]*p[1])

        def grad(X):
            p = fr.chart(X)
            s = p[0]*p[0]+p[1]*p[1]+p[2]*p[2]
            g = 4*(s+k)
            return fr.covector((g*p[0]-8*R2*p[0], g*p[1]-8*R2*p[1], g*p[2]))
        self.surfs = [Surf('torus', 4, f, grad)]
        if self.segment:
            for e, name in zip(P.ends, ('low', 'high')):
                self.surfs.append(plane_surf(name, fr, (F(0), F(0), F(1)), e['z']))
        else:
            cu, su = P.eu
            self.surfs.append(plane_surf('start', fr, (F(0), F(1), F(0)), F(0)))
            self.surfs.append(plane_surf('end', fr, (-su, cu, F(0)), F(0)))
        self.faces = self.make_faces()

    def contains(self, X):
        return self.part.contains(self.fr.chart(X))

    def e(self, t):
        fr = self.fr
        c, s = mp.cos(t), mp.sin(t)
        return lin((c, fr.xm), (s, fr.ym)), lin((-s, fr.xm), (c, fr.ym))

    def make_faces(self):
        fr, R, r, P = self.fr, self.Rm, self.rm, self.part
        zero = (Z, Z, Z)
        if self.segment:
            vlo, vhi, sg = P.vlo, P.vhi, P.sigma

            def tube(th):
                e, de = self.e(th)
                return circle_curve(lin((1, fr.om), (R, e)), scale(e, r), scale(fr.nm, r),
                                    scale(de, R), scale(de, r), zero, (vlo, vhi))

            def parallel(ph):
                c, s = mp.cos(ph), mp.sin(ph)
                rho = R+r*c
                return circle_curve(lin((1, fr.om), (r*s, fr.nm)), scale(fr.xm, rho), scale(fr.ym, rho),
                                    scale(fr.nm, r*c), scale(fr.xm, -r*s), scale(fr.ym, -r*s))
            faces = [Face('wall', [Family('tube circles', Z, 2*mp.pi, True, tube, sg),
                                   Family('parallels', vlo, vhi, False, parallel, -sg)])]
            for e, name in zip(P.ends, ('low', 'high')):
                out = -e['kside']
                faces.append(Face(name, [disc_chords(fr, e['zm'], e['rho'], (F(1), F(0)), out, 'chords x'),
                                         disc_chords(fr, e['zm'], e['rho'], (F(2), F(7)), out, 'chords skew')]))
            return faces
        Um = P.Um

        def tube(th):
            e, de = self.e(th)
            return circle_curve(lin((1, fr.om), (R, e)), scale(e, r), scale(fr.nm, r),
                                scale(de, R), scale(de, r), zero)

        def parallel(ph):
            c, s = mp.cos(ph), mp.sin(ph)
            rho = R+r*c
            return circle_curve(lin((1, fr.om), (r*s, fr.nm)), scale(fr.xm, rho), scale(fr.ym, rho),
                                scale(fr.nm, r*c), scale(fr.xm, -r*s), scale(fr.ym, -r*s), (Z, Um))
        faces = [Face('wall', [Family('tube circles', Z, Um, False, tube, 1),
                               Family('parallels', Z, 2*mp.pi, True, parallel, -1)])]
        one, zero_ = mp.mpf(1), Z
        (c0, s0), (c1, s1) = P.end_units()
        n = (zero_, zero_, one)
        for (c, s), outward, name in (((c0, s0), (zero_, -one, zero_), 'start'), ((c1, s1), (-s1, c1, zero_), 'end')):
            U = (c, s, zero_)
            centre = (R*c, R*s, zero_)
            faces.append(Face(name, [chart_disc(fr, centre, U, n, r, (F(1), F(0)), outward, 'chords radial'),
                                     chart_disc(fr, centre, U, n, r, (F(2), F(7)), outward, 'chords skew')]))
        return faces

    def closed(self):
        return self.part.closed()

    def box(self):
        R, r = float(self.R), float(self.r)
        return frame_box(self.fr, [(u, v, w) for u in (-R-r, R+r) for v in (-R-r, R+r) for w in (-r, r)])

    def rims(self):
        """Its rims as curves: a segment's end parallels, a wedge's end tube
        circles."""
        wall = self.faces[0]
        if self.segment:
            fam = wall.families[1]
            return [fam.curve(fam.a0), fam.curve(fam.a1)]
        fam = wall.families[0]
        return [fam.curve(fam.a0), fam.curve(fam.a1)]

    def slice_intervals(self, s):
        """The part's section at the chart level `s` (binary64) as radial
        intervals, and whether it is cut to the turn (a wedge's)."""
        P = self.part
        if self.segment:
            return [(float(a), float(b)) for a, b in P.rho_intervals(mp.mpf(s))], None
        R, r = float(self.R), float(self.r)
        if abs(s) >= r:
            return [], float(P.Um)
        q = math.sqrt(r*r-s*s)
        return [(R-q, R+q)], float(P.Um)


# ------------------------------------------------------------------ a cap or zone

class Cap:
    """A sphere's cap or zone (S9d.2c's `AxisSphere` model)."""
    kind = 'sphere'

    def __init__(self, case):
        self.case = case
        S = self.S = AxisSphere(case)
        self.fr = fr = Frame(case.frame)
        self.c, self.R = S.c, S.r
        self.cm, self.Rm = Mv(self.c), M(self.R)
        cm, R2 = self.cm, self.Rm**2
        self.mhat = unit(fr.nm)
        x = fr.xm
        xm = dot(x, self.mhat)
        self.e1 = unit(tuple(x[i]-xm*self.mhat[i] for i in range(3)))
        self.e2 = cross(self.mhat, self.e1)
        self.zlo, self.zhi = S.z_range()
        self.blo = -mp.pi/2 if S.heights[0] is None else mp.asin(self.zlo/self.Rm)
        self.bhi = mp.pi/2 if S.heights[1] is None else mp.asin(self.zhi/self.Rm)
        self.surfs = [Surf('sphere', 2, lambda X: (X[0]-cm[0])**2+(X[1]-cm[1])**2+(X[2]-cm[2])**2-R2,
                           lambda X: (2*(X[0]-cm[0]), 2*(X[1]-cm[1]), 2*(X[2]-cm[2])))]
        self.planes = []
        for a, b, label in S.planes:
            am, bm = Mv(a), M(b)
            self.planes.append((am, bm))
            self.surfs.append(Surf(label, 1, lambda X, am=am, bm=bm: bm-dot(am, X),
                                   lambda X, am=am: (-am[0], -am[1], -am[2])))
        self.faces = self.make_faces()

    def contains(self, X):
        if self.surfs[0].f(X) >= 0:
            return False
        return all(dot(a, X) > b for a, b in self.planes)

    def make_faces(self):
        c, R, m, e1, e2 = self.cm, self.Rm, self.mhat, self.e1, self.e2
        blo, bhi = self.blo, self.bhi
        zero = (Z, Z, Z)

        def meridian(lam):
            cl, sl = mp.cos(lam), mp.sin(lam)
            return circle_curve(c, lin((R*cl, e1), (R*sl, e2)), scale(m, R), zero, lin((-R*sl, e1), (R*cl, e2)),
                                zero, (blo, bhi))

        def parallel(beta):
            cb, sb = mp.cos(beta), mp.sin(beta)
            return circle_curve(lin((1, c), (R*sb, m)), scale(e1, R*cb), scale(e2, R*cb), scale(m, R*cb),
                                scale(e1, -R*sb), scale(e2, -R*sb))
        faces = [Face('wall', [Family('meridians', Z, 2*mp.pi, True, meridian, 1),
                               Family('parallels', blo, bhi, False, parallel, -1)])]
        for z, h, name, out in ((self.zlo, self.S.heights[0], 'low', -1), (self.zhi, self.S.heights[1], 'high', 1)):
            if h is None:
                continue
            centre = lin((1, c), (z, m))
            rad = mp.sqrt(R*R-z*z)
            outward = scale(m, out)
            faces.append(Face(name, [world_disc(centre, e1, e2, rad, (F(1), F(0)), outward, 'chords x'),
                                     world_disc(centre, e1, e2, rad, (F(2), F(7)), outward, 'chords skew')]))
        return faces

    def closed(self):
        return self.S.closed()

    def box(self):
        c, R = [float(v) for v in self.c], float(self.R)
        return [v-R for v in c], [v+R for v in c]

    def rims(self):
        fam = self.faces[0].families[1]
        out = []
        if self.S.heights[0] is not None:
            out.append(fam.curve(self.blo))
        if self.S.heights[1] is not None:
            out.append(fam.curve(self.bhi))
        return out


def make_input(case):
    if case.torus is not None:
        R, r, low, high, angle = case.torus
        if high-low == TWO_PI and angle == TWO_PI:
            return tc.Torus(case)
        return TorusPart(case)
    if case.sphere is not None:
        radius, low, high = case.sphere
        if low == -HALF_PI and high == HALF_PI:
            return tc.Sphere(case)
        return Cap(case)
    if case.cone is not None:
        return tc.Cone(case)
    return tc.Prism(case)


def is_part(S):
    return isinstance(S, TorusPart)


def is_cap(S):
    return isinstance(S, Cap)


# ------------------------------------------------------------------ binary64 models

class FloatPart(tc.FloatModel):
    """A torus part in binary64: the torus's quartic and its end planes along
    a ray, its membership by S9d.4b.1's model."""

    def __init__(self, S):
        fr = S.fr
        self.S = S
        self.o = [float(v) for v in fr.om]
        self.inv = [[float(v) for v in row] for row in fr.invm]
        self.kind = 'part'
        self.R, self.r = float(S.R), float(S.r)
        P = S.part
        if S.segment:
            self.vlo, self.vhi = float(P.vlo), float(P.vhi)
            self.zband = [float(z) for z in P.zband]
            self.zs = [float(e['z']) for e in P.ends]
        else:
            self.eu = [float(c) for c in P.eu]
            self.Um = float(P.Um)

    def polys(self, P, D):
        p0, p1 = self.chart(P), self.chart_vec(D)
        lin_ = lambda i: [p0[i], p1[i]]
        sq = lambda i: fpoly_mul(lin_(i), lin_(i))
        s = fpoly_add(fpoly_add(sq(0), sq(1)), sq(2))
        s[0] += self.R*self.R-self.r*self.r
        h = fpoly_add(sq(0), sq(1))
        out = [fpoly_add(fpoly_mul(s, s), [-4*self.R*self.R*c for c in h])]
        if self.S.segment:
            out += [[p0[2]-z, p1[2]] for z in self.zs]
        else:
            cu, su = self.eu
            out += [lin_(1), [-su*p0[0]+cu*p0[1], -su*p1[0]+cu*p1[1]]]
        return out

    def contains(self, X):
        p = self.chart(X)
        R, r = self.R, self.r
        if not self.S.segment:
            s = p[0]*p[0]+p[1]*p[1]+p[2]*p[2]+R*R-r*r
            if not s*s < 4*R*R*(p[0]*p[0]+p[1]*p[1]):
                return False
            cu, su = self.eu
            a, b = p[1] > 0, p[0]*su-p[1]*cu > 0
            return (a and b) if self.Um <= math.pi else (a or b)
        w = p[2]
        if abs(w) >= r:
            return False
        t = math.hypot(p[0], p[1])
        a = math.asin(w/r)
        T = 2*math.pi

        def in_range(v):
            k = math.ceil((self.vlo-v)/T)
            return v+k*T <= self.vhi
        sides = [s for s, ang in ((1, a), (-1, math.pi-a)) if in_range(ang)]
        q = math.sqrt(r*r-w*w)
        if len(sides) == 1:
            return t < R+sides[0]*q
        if len(sides) == 2:
            return R-q < t < R+q
        return False


class FloatCap(tc.FloatModel):
    """A cap or zone in binary64."""

    def __init__(self, S):
        self.S = S
        self.kind = 'cap'
        self.c, self.R2 = [float(v) for v in S.cm], float(S.R)**2
        self.planes = [([float(v) for v in a], float(b)) for a, b in S.planes]

    def polys(self, P, D):
        q = [P[i]-self.c[i] for i in range(3)]
        out = [[sum(x*x for x in q)-self.R2, 2*sum(q[i]*D[i] for i in range(3)), sum(x*x for x in D)]]
        for a, b in self.planes:
            out.append([sum(a[i]*P[i] for i in range(3))-b, sum(a[i]*D[i] for i in range(3))])
        return out

    def contains(self, X):
        if sum((X[i]-self.c[i])**2 for i in range(3)) >= self.R2:
            return False
        return all(sum(a[i]*X[i] for i in range(3)) > b for a, b in self.planes)


def float_model(S):
    if is_part(S):
        return FloatPart(S)
    if is_cap(S):
        return FloatCap(S)
    return tc.FloatModel(S)


# ------------------------------------------------------------------ the pair

class Pair(tc.Pair):
    """A torus (whole or a part) and a prism with arcs, a sphere (whole, a
    cap or a zone), a cone or a whole torus, in either order: `A` the
    object, `B` the tool; `K` the part, else the torus."""

    def __init__(self, obj, tool):
        self.A, self.B = make_input(obj), make_input(tool)
        tori = [S for S in (self.A, self.B) if S.kind == 'torus']
        assert tori, 'S9d.4c: a torus or a torus part'
        parts = [S for S in tori if is_part(S)]
        assert parts or any(is_cap(S) for S in (self.A, self.B)), 'S9d.4c: a torus part or a cap'
        self.K = parts[0] if parts else tori[0]
        self.P = self.B if self.K is self.A else self.A
        self.size = mp.mpf(box_size([self.A.box(), self.B.box()]))
        self.c0 = self.K.fr.om
        self._sweeps = None
        self._counts = {}

    def swept_count(self, key):
        if key not in self._counts:
            self._counts.update(sweep_counts(self))
        return self._counts[key]


def sweep_counts(pair, levels=200, rays=512):
    """Solids of every operation by a sweep of `K`'s normal slices (S9d.4b.2's
    with the part's own sections). Intervals on one ray at adjacent levels
    join where they overlap, or where a chain of intervals at levels between
    them (bisected, to 2^-16 of the spacing; an interval below within twice
    the spacing) overlaps from one to the other:
    a thin region sheared between the levels (a wedge's tip, where two
    boundaries meet at different slopes) is still one solid."""
    K = pair.K
    fK = float_model(K)
    fA, fB = float_model(pair.A), float_model(pair.B)
    o = fK.o
    x, y, n = ([float(v) for v in K.fr.xm], [float(v) for v in K.fr.ym], [float(v) for v in K.fr.nm])
    corners = []
    for lo, hi in (pair.A.box(), pair.B.box()):
        corners += [fK.chart([(lo, hi)[i >> k & 1][k] for k in range(3)]) for i in range(8)]
    wmin, wmax = min(c[2] for c in corners)-0.01, max(c[2] for c in corners)+0.01
    tmax = max(math.hypot(c[0], c[1]) for c in corners)*1.5+1
    R, r = float(K.R), float(K.r)
    parent = {}

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb
    keys = ('common', 'fuse', 'A-B', 'B-A')
    slices = {}

    def k_slice(s):
        if s not in slices:
            if is_part(K):
                slices[s] = K.slice_intervals(s)
            else:
                slices[s] = ([(R-math.sqrt(r*r-s*s), R+math.sqrt(r*r-s*s))] if abs(s) < r else []), None
        return slices[s]

    def ops_at(s, i):
        kiv, turn = k_slice(s)
        th = 2*math.pi*(i+0.25)/rays
        c, sn = math.cos(th), math.sin(th)
        P = [o[k]+s*n[k] for k in range(3)]
        D = [c*x[k]+sn*y[k] for k in range(3)]
        iv = {}
        for role, fm in (('A', fA), ('B', fB)):
            if fm.S is K:
                iv[role] = kiv if turn is None or th < turn else []
            else:
                iv[role] = fm.intervals(P, D, tmax)
        return iv_ops(iv['A'], iv['B'])

    overlap = lambda u, v: min(u[1], v[1])-max(u[0], v[0]) > 1e-9

    def chained(lo_iv, s_lo, hi_iv, s_hi, i, key, depth=0):
        """Whether an interval at `s_lo` reaches one at `s_hi` on ray `i`
        through intervals at levels between them."""
        if overlap(lo_iv, hi_iv):
            return True
        if depth >= 16:
            return False
        s_mid = 0.5*s_lo+0.5*s_hi
        for mid in ops_at(s_mid, i)[key]:
            # One step at a time: a middle interval touching one end, then
            # the rest of the way from it.
            if overlap(lo_iv, mid) and chained(mid, s_mid, hi_iv, s_hi, i, key, depth+1):
                return True
            if overlap(mid, hi_iv) and chained(lo_iv, s_lo, mid, s_mid, i, key, depth+1):
                return True
        return False
    prev, s_prev = None, None
    for lv in range(levels):
        s = wmin+(wmax-wmin)*(lv+0.5)/levels
        row = []
        for i in range(rays):
            ops = ops_at(s, i)
            for key in keys:
                for j, _ in enumerate(ops[key]):
                    parent[(key, lv, i, j)] = (key, lv, i, j)
            row.append(ops)
        for key in keys:
            axis = None
            for i in range(rays):
                cur, nxt = row[i][key], row[(i+1) % rays][key]
                for j, (a, b) in enumerate(cur):
                    if a <= 1e-9:
                        if axis is not None:
                            union(axis, (key, lv, i, j))
                        axis = (key, lv, i, j)
                    for jj, (a2, b2) in enumerate(nxt):
                        if min(b, b2)-max(a, a2) > 1e-9:
                            union((key, lv, i, j), (key, lv, (i+1) % rays, jj))
                    if prev is not None:
                        below = prev[i][key]
                        hits = [jj for jj, v in enumerate(below) if overlap((a, b), v)]
                        if not hits:
                            # No overlap: a chain between the levels, for an
                            # interval below within twice the spacing.
                            gap = lambda v: max(v[0]-b, a-v[1])
                            hits = [jj for jj, v in enumerate(below)
                                    if gap(v) < 2*(s-s_prev) and chained(v, s_prev, (a, b), s, i, key)]
                        for jj in hits:
                            union((key, lv, i, j), (key, lv-1, i, jj))
        prev, s_prev = row, s
    counts = {}
    for key in keys:
        counts[key] = len({find(a) for a in parent if a[0] == key})
    return counts


# ------------------------------------------------------------------ margins

def edges_of(S):
    if is_part(S) or is_cap(S):
        return S.rims()
    return tc.edges_of(S)


def margins(pair):
    """S9d.4b.2's margins with every input's edges against the other's
    surfaces (a cap's and a part's rims too), relative to the case's
    size."""
    size = pair.size
    out = {'face_crossing': mp.inf, 'face_gap': mp.inf, 'edge_crossing': mp.inf, 'edge_gap': mp.inf,
           'vertex': mp.inf}
    for (k, role, name), sw in pair.sweeps().items():
        if k != 0:
            continue
        O = pair.B if role == 'A' else pair.A
        crossed = set()
        for a, cur, pcs in sw.samples:
            n = len(pcs)
            for i in range(n):
                p, q = pcs[i], pcs[(i+1) % n]
                closed = cur.kind == 'circle' and cur.rng is None
                if i+1 == n and not closed:
                    continue
                if p[2] == q[2] or p[4] in (None, 'end'):
                    continue
                j = p[4]
                crossed.add(j)
                X, Xa, Xb = cur.frame_at(p[1])
                out['face_crossing'] = min(out['face_crossing'], sine(cross(Xa, Xb), O.surfs[j].grad(X)))
        for j, surf in enumerate(O.surfs):
            if j in crossed:
                continue
            for a, cur, pcs in sw.samples[::4]:
                if cur.degenerate():
                    continue
                for b in critical_points(cur, surf):
                    X = cur.point(b)
                    if tc.on_face(O, surf, X, size):
                        out['face_gap'] = min(out['face_gap'], distance_estimate(surf, X)/size)
    for S in (pair.A, pair.B):
        O = pair.B if S is pair.A else pair.A
        for cur in edges_of(S):
            if cur.degenerate():
                continue
            for surf in O.surfs:
                for b in curve_roots(cur, surf):
                    X, _, Xb = cur.frame_at(b)
                    if not tc.on_face(O, surf, X, size):
                        continue
                    g = surf.grad(X)
                    out['edge_crossing'] = min(out['edge_crossing'], abs(dot(Xb, g))/(norm(Xb)*norm(g)))
                for b in critical_points(cur, surf):
                    X = cur.point(b)
                    if tc.on_face(O, surf, X, size):
                        out['edge_gap'] = min(out['edge_gap'], distance_estimate(surf, X)/size)
        for X in tc.vertices_of(S):
            for surf in O.surfs:
                out['vertex'] = min(out['vertex'], distance_estimate(surf, X)/size)
    return out


def near_coincidences(pair):
    return tc.near_coincidences(pair)


def rows(obj, operation, tool, pair=None):
    """`result N volume area cx cy cz` or `empty`, and the pair (reused
    across the three operations)."""
    pair = pair or Pair(obj, tool)
    n, vol, surface, centre = pair.result(operation)
    if n == 0:
        return ['empty'], pair
    return [' '.join(['result', str(n), number(vol), number(surface)]+[number(c) for c in centre])], pair
