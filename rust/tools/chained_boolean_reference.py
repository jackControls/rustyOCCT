#!/usr/bin/env python3
"""Independent reference for S9e.1 of REVIEW_NOTES.md: a Boolean's result
given to another Boolean, `(A op1 B) op2 C` or, swapped, `C op2 (A op1 B)`,
for three prisms of line, arc and circle profiles (planar and cylindrical
walls) in any relative position where every pair of faces meets in lines,
circles or ellipses (S9c.1's pairs).

It is S9c.1's reference (`curved_boolean_reference.py`: each prism its
construction's exact model on its frame's stored axes, sliced in planes
holding the prisms' axes, every face swept by the lines of its own
parameterisation) generalized from two prisms and a set function of two to
three prisms and the chained set function of three, `f(a, b, c) =
op2(op1(a, b), c)` (swapped: `op2(c, op1(a, b))`). Nothing here uses the
kernel, the first operation's result or a 3D arrangement: the chain is one
point set, decided at every point by the three prisms' memberships.

* **Slicing (volume, first moments, solids).** The planes `d . X = s` hold
  every prism's axis, so the three axes take at most two directions: `d =
  n_1 x n_2` (two directions, `g = n_1`, `h = n_2`) or `d = n x e` for a
  world axis `e` (one direction). Each prism's section is S9c.1's union of
  parallelograms. Each slice's first result is the convex pieces of
  `op1(A, B)` (S9c.1's clipping: common the pieces' intersections, cut the
  object's parallelograms less the tool's by successive half-plane
  complements, fuse the cut and the tool's), and the chain the convex
  pieces of `op2` between those pieces and the third's parallelograms (a
  convex piece's half-planes its edges); their areas and first moments by
  Green's theorem. The breakpoints are S9c.1's over the three sections'
  lines: two parallel lines of different sections coinciding, three lines
  not all of one section concurrent, each section's own structure, all
  roots of exact polynomials (the surds squared away, spurious roots
  filtered by the meeting lying on the closed boundaries of the sections
  whose lines meet). Gauss-Legendre between breakpoints as S9c.1's. Solids
  by S9c.1's combinatorics with the chained function.
* **Surface area.** Every face of each prism is swept by the lines of its
  own parameterisation against both other prisms at once (S9c.1's
  `FaceSweep` against each, their crossings merged along the line, the
  breakpoints each sweep's and those where a crossing of one other meets a
  crossing of the other); each piece of a line is classified against each
  other prism at its midpoint (inside, outside, or on a face of it with the
  same or the opposite orientation). A piece bounds the chain where the
  chained function differs across the face: in front the face's own prism
  is left, behind it is held, the others as classified (on a coincident
  face of the same orientation left in front and held behind, of the
  opposite orientation the reverse). A piece on faces of several prisms is
  counted once, by the first of them in the order object, tool, third.

`Chain(obj, tool, third, op1, swapped)` gives `result(op2)`: the solid
count, volume, area and centre, or no solid; `rows(...)` gives `result N
volume area cx cy cz` or `empty`, as the other references.
"""
import itertools

import mpmath as mp

import curved_boolean_reference as cref
from curved_boolean_reference import (
    F, M, Mv, P, Desc, FaceSweep, Prism, Section, clip, clip_f, area_f, corner, cross, cross2, det3,
    difference, dot, integrate, merge_breaks, overlap_area, padd, poly_moments, pscale, psub, pzero,
    real_roots, same_line, scale, sq_poly, sub)

mp.mp.dps = 40

OPS = cref.OPS
# Corners closer than this are one (the cases' coordinates are below 100).
EPS = mp.mpf(10)**-32
SET = cref.SET


def chained(op1, op2, swapped):
    """The chain's set function of the three memberships."""
    if swapped:
        return lambda a, b, c: SET[op2](c, SET[op1](a, b))
    return lambda a, b, c: SET[op2](SET[op1](a, b), c)


# ------------------------------------------------------------------ convex pieces

def halfplanes_of(poly):
    """A counter-clockwise convex polygon's half-planes `a x + b y <= c`."""
    out = []
    n = len(poly)
    for i in range(n):
        (px, py), (qx, qy) = poly[i], poly[(i+1) % n]
        a, b = qy-py, -(qx-px)
        out.append((a, b, a*px+b*py))
    return out


def ccw(poly):
    """Counter-clockwise, without repeated corners (a clip keeps a corner
    on the clipping line and adds the crossing there too)."""
    out = []
    for p in poly:
        if not out or abs(p[0]-out[-1][0])+abs(p[1]-out[-1][1]) > EPS:
            out.append(p)
    while len(out) > 1 and abs(out[0][0]-out[-1][0])+abs(out[0][1]-out[-1][1]) <= EPS:
        out.pop()
    if len(out) < 3:
        return []
    if poly_moments(out)[0] < 0:
        return list(reversed(out))
    return out


def piece_common(P_, Q_):
    """Disjoint convex pieces of the intersection of two sets given as
    disjoint convex pieces with their half-planes."""
    out = []
    for p, _ in P_:
        for _, hq in Q_:
            q = p
            for h in hq:
                q = clip(q, *h)
                if not q:
                    break
            q = ccw(q) if q else q
            if q:
                out.append((q, halfplanes_of(q)))
    return out


def piece_cut(P_, Q_):
    out = []
    for p, _ in P_:
        rest = [p]
        for _, hq in Q_:
            rest = [piece for r in rest for piece in difference(r, hq)]
            if not rest:
                break
        for r in rest:
            r = ccw(r)
            if r:
                out.append((r, halfplanes_of(r)))
    return out


def piece_op(op, P_, Q_):
    if op == 'common':
        return piece_common(P_, Q_)
    if op == 'cut':
        return piece_cut(P_, Q_)
    return piece_cut(P_, Q_)+list(Q_)


# ------------------------------------------------------------------ slicing

class Slicing3:
    """The three prisms sliced by planes holding their axes."""

    def __init__(self, prisms, axis=None):
        self.prisms = prisms
        dirs = []
        for p in prisms:
            if not any(cross(p.n, q) == (0, 0, 0) for q in dirs):
                dirs.append(p.n)
        assert len(dirs) <= 2, 'three axis directions: no slicing holds them all'
        if len(dirs) == 1:
            self.parallel = True
            n = dirs[0]
            axes = [(F(1), F(0), F(0)), (F(0), F(1), F(0)), (F(0), F(0), F(1))]
            e = axis if axis is not None else min(axes, key=lambda e: abs(dot(n, e)))
            d = cross(n, e)
            g, h = n, cross(n, d)
        else:
            self.parallel = False
            d, g, h = cross(dirs[0], dirs[1]), dirs[0], dirs[1]
        self.d, self.g, self.h = d, g, h
        self.dp = scale(d, 1/dot(d, d))
        self.P0 = prisms[0].o
        self.J = abs(det3(self.dp, g, h))
        self.secs = [Section(p, self.P0, self.dp, g, h, tag) for p, tag in zip(prisms, 'ABC')]
        self.ranges = [self.s_range(p) for p in prisms]
        self.size = max(max(abs(x) for x in bnd) for p in prisms for bnd in p.bounds())+1
        corners = []
        for p in prisms:
            b0, b1 = p.bounds()
            for i in range(8):
                X = tuple(b1[k] if (i >> k) & 1 else b0[k] for k in range(3))
                corners.append(self.plane_coords(X))
        self.box = (min(c[1] for c in corners), max(c[1] for c in corners),
                    min(c[2] for c in corners), max(c[2] for c in corners))

    def plane_coords(self, X):
        R = sub(X, self.P0)
        det = det3(self.dp, self.g, self.h)
        return det3(R, self.g, self.h)/det, det3(self.dp, R, self.h)/det, det3(self.dp, self.g, R)/det

    def s_range(self, prism):
        p = prism.profile
        dx, dy = dot(self.d, prism.x), dot(self.d, prism.y)
        base = dot(self.d, sub(prism.o, self.P0))
        vals = [base+v[0]*dx+v[1]*dy for v in p.vertices]
        for el in p.elements:
            if el.kind != 'seg':
                c = base+el.c[0]*dx+el.c[1]*dy
                vals += [c-el.r*(abs(dx)+abs(dy)), c+el.r*(abs(dx)+abs(dy))]
        return min(vals), max(vals)

    # ---- the integrand

    def pieces(self, s, op1, swapped):
        """Each second operation's convex pieces of the slice."""
        pa, pb, pc = (sec.parallelograms(s) for sec in self.secs)
        first = piece_op(op1, pa, pb)
        out = {}
        for op2 in OPS:
            pieces = piece_op(op2, pc, first) if swapped else piece_op(op2, first, pc)
            out[op2] = [q for q, _ in pieces]
        return out

    def integrand_for(self, op1, swapped):
        def integrand(s):
            out = []
            pieces = self.pieces(s, op1, swapped)
            for op in OPS:
                A = X = Y = 0
                for q in pieces[op]:
                    a, x, y = poly_moments(q)
                    A, X, Y = A+a, X+x, Y+y
                out += [A, s*A, X, Y]
            return out
        return integrand

    # ---- breakpoints

    def breakpoints(self):
        lo = min(r[0] for r in self.ranges)
        hi = max(r[1] for r in self.ranges)
        cand = [M(r[0]) for r in self.ranges]+[M(r[1]) for r in self.ranges]
        for sec in self.secs:
            for p in sec.structure:
                cand += real_roots(p)
        tol = mp.mpf(10)**-20*self.size
        lines = [l for sec in self.secs for l in sec.lines]
        found = []
        # Two parallel lines of different sections coinciding.
        for i in range(len(lines)):
            for j in range(i+1, len(lines)):
                (ab1, d1), (ab2, d2) = lines[i], lines[j]
                if d1.key[0] == d2.key[0] or cross2(ab1, ab2) != 0:
                    continue
                kappa = ab1[0]/ab2[0] if ab2[0] else ab1[1]/ab2[1]
                L = psub(d1.L, pscale(d2.L, kappa))
                poly = sq_poly(L, [(d1.k, d1.Q), (-kappa*d2.k, d2.Q)])
                if pzero(poly):
                    continue
                for s in real_roots(poly):
                    if M(lo)-tol <= s <= M(hi)+tol and self.relevant(s, (d1.key, d2.key), tol):
                        found.append(s)
        # Three pairwise non-parallel lines, not all of one section, concurrent.
        n = len(lines)
        for i in range(n):
            for j in range(i+1, n):
                for k in range(j+1, n):
                    trip = (lines[i], lines[j], lines[k])
                    if len({l[1].key[0] for l in trip}) < 2:
                        continue
                    abs_ = [l[0] for l in trip]
                    if any(cross2(abs_[p], abs_[q]) == 0 for p, q in ((0, 1), (0, 2), (1, 2))):
                        continue
                    (a1, b1), (a2, b2), (a3, b3) = abs_
                    cof = (a2*b3-a3*b2, -(a1*b3-a3*b1), a1*b2-a2*b1)
                    L = [F(0)]
                    terms = []
                    for cf, (_, d) in zip(cof, trip):
                        L = padd(L, pscale(d.L, cf))
                        if d.k:
                            terms.append((cf*d.k, d.Q))
                    # Three circles' lines (sq_poly squares two surds
                    # away): no fixture has them.
                    assert len(terms) <= 2, 'three circle lines concurrent'
                    poly = sq_poly(L, terms)
                    if pzero(poly):
                        continue
                    for s in real_roots(poly):
                        if M(lo)-tol <= s <= M(hi)+tol and \
                                self.relevant(s, tuple(l[1].key for l in trip), tol):
                            found.append(s)
        self.events = len(found)
        return merge_breaks(cand+found, M(lo), M(hi), mp.mpf(10)**-30*self.size)

    def relevant(self, s, keys, tol):
        """Whether the lines `keys` (some branch of each) meet at `s` on the
        closed boundaries of the sections they belong to."""
        inst = {}
        for sec in self.secs:
            for key, line, bad in sec.instances(s, tol):
                if bad <= tol:
                    inst.setdefault(key[:3], []).append((line, sec))
        options = []
        for key in keys:
            if key not in inst:
                return False
            options.append(inst[key])
        owners = {k[0] for k in keys}
        involved = [sec for sec in self.secs if sec.tag in owners]
        for combo in itertools.product(*options):
            ls = [c[0] for c in combo]
            if len(ls) == 2:
                (a1, b1, c1), (a2, b2, c2) = ls
                kappa = a1/a2 if abs(a2) > abs(b2) else b1/b2
                if abs(c1-kappa*c2) > tol*(1+abs(kappa)):
                    continue
                if self.portions_overlap(s, combo, involved, tol):
                    return True
                continue
            pair = None
            for p, q in ((0, 1), (0, 2), (1, 2)):
                if abs(ls[p][0]*ls[q][1]-ls[q][0]*ls[p][1]) > 0:
                    pair = (p, q)
                    break
            X = corner(ls[pair[0]], ls[pair[1]])
            r = [l for idx, l in enumerate(ls) if idx not in pair][0]
            if abs(r[0]*X[0]+r[1]*X[1]-r[2]) > tol*mp.sqrt(r[0]**2+r[1]**2):
                continue
            if all(sec.closure_bad(s, *X) <= tol for sec in involved):
                return True
        return False

    def portions_overlap(self, s, combo, involved, tol):
        (a, b, c), _ = combo[0]
        norm = a*a+b*b
        foot = (a*c/norm, b*c/norm)
        direction = (-b, a)
        ts = []
        for sec in self.secs:
            for key, (a2, b2, c2), bad in sec.instances(s, tol):
                det = a*b2-a2*b
                if abs(det) <= tol*mp.sqrt(norm*(a2*a2+b2*b2)) or bad > tol:
                    continue
                X = corner((a, b, c), (a2, b2, c2))
                ts.append(((X[0]-foot[0])*direction[0]+(X[1]-foot[1])*direction[1])/norm)
        ts.sort()
        ts = ts+[(p+q)/2 for p, q in zip(ts, ts[1:])]
        for t in ts:
            X = (foot[0]+t*direction[0], foot[1]+t*direction[1])
            if all(sec.closure_bad(s, *X) <= tol for sec in involved):
                return True
        return False

    # ---- measures

    def measure(self, op1, swapped):
        if not hasattr(self, 'breaks'):
            self.breaks = self.breakpoints()
        breaks = self.breaks
        tol = mp.mpf(10)**-33*self.size**4
        total = [mp.mpf(0)]*12
        err = mp.mpf(0)
        integrand = self.integrand_for(op1, swapped)
        for a, b in zip(breaks, breaks[1:]):
            est, diff = integrate(integrand, a, b, tol)
            if est is None:
                continue
            total = [x+y for x, y in zip(total, est)]
            err = max(err, diff)
        self.quad_error = err
        out = {}
        J = M(self.J)
        P0, dp, g, h = Mv(self.P0), Mv(self.dp), Mv(self.g), Mv(self.h)
        for k, op in enumerate(OPS):
            A, SA, X, Y = total[4*k:4*k+4]
            vol = J*A
            mom = tuple(J*(P0[i]*A+dp[i]*SA+g[i]*X+h[i]*Y) for i in range(3))
            out[op] = (vol, mom)
        return out

    # ---- solids

    def solids(self, fn):
        """The chain's number of solids for the set function `fn` (floats:
        combinatorics only)."""
        breaks = self.breaks
        box = [float(x) for x in self.box]
        margin = max(box[1]-box[0], box[3]-box[2])+1.0
        bx = [(box[0]-margin, box[2]-margin), (box[1]+margin, box[2]-margin),
              (box[1]+margin, box[3]+margin), (box[0]-margin, box[3]+margin)]
        thr = 1e-11*margin*margin
        if not hasattr(self, '_intervals'):
            canon = self.canonical_lines()
            intervals = []
            for a, b in zip(breaks, breaks[1:]):
                if b-a <= mp.mpf(10)**-9*self.size:
                    intervals.append(None)
                    continue
                intervals.append(self.faces_in(a, b, bx, thr, canon))
            self._intervals = intervals
        parent = {}

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        def union(x, y):
            parent[find(x)] = find(y)

        live = []
        for k, iv in enumerate(self._intervals):
            if iv is None:
                continue
            faces = [f for f in iv if fn(*f['in'])]
            for f in faces:
                parent[(k, f['id'])] = (k, f['id'])
            for x in range(len(faces)):
                for y in range(x+1, len(faces)):
                    sx, sy = faces[x]['signs'], faces[y]['signs']
                    if sum(1 for key in sx if sx[key] != sy.get(key)) == 1 and set(sx) == set(sy):
                        union((k, faces[x]['id']), (k, faces[y]['id']))
            live.append((k, faces))
        for (k1, f1), (k2, f2) in zip(live, live[1:]):
            for f in f1:
                if f['end'] is None:
                    continue
                for g in f2:
                    if g['start'] is None:
                        continue
                    if overlap_area(f['end'], g['start']) > thr:
                        union((k1, f['id']), (k2, g['id']))
        return len({find(x) for x in parent})

    def canonical_lines(self):
        items = []
        for sec in self.secs:
            for ab, d in sec.lines:
                for br in d.branches():
                    key = d.key+(br,) if d.key[1] == 'chord' else d.key+(0,)
                    items.append((key, ab, d, br))
        canon = {}
        for i, (key, ab, d, br) in enumerate(items):
            canon.setdefault(key, key)
            for key2, ab2, d2, br2 in items[:i]:
                if canon.get(key2) != key2:
                    continue
                if same_line(ab, d, br, ab2, d2, br2):
                    canon[key] = key2
                    break
        self.line_info = {key: (ab, d, br) for key, ab, d, br in items}
        return canon

    def line_at(self, key, s):
        ab, d, br = self.line_info[key]
        return (float(ab[0]), float(ab[1]), float(d.value(s, br, P(1))))

    def faces_in(self, a, b, bx, thr, canon):
        s = (a+b)/2
        lines = {}
        for sec in self.secs:
            for (lo_, ta), (hi_, tb) in sec.chords(s):
                for t, tag in ((lo_, ta), (hi_, tb)):
                    lines[canon[(sec.tag, 'chord', tag[0], tag[1])]] = None
            lines[canon[(sec.tag, 'height', 'lo', 0)]] = None
            lines[canon[(sec.tag, 'height', 'hi', 0)]] = None
        keys = sorted(lines)
        vals = {key: self.line_at(key, s) for key in keys}
        faces = [(bx, {})]
        for key in keys:
            la, lb, lc = vals[key]
            nxt = []
            for poly, signs in faces:
                below = clip_f(poly, la, lb, lc)
                above = clip_f(poly, -la, -lb, -lc)
                for part, sg in ((below, -1), (above, 1)):
                    if part and abs(area_f(part)) > thr:
                        ns = dict(signs)
                        ns[key] = sg
                        nxt.append((part, ns))
            faces = nxt
        out = []
        for idx, (poly, signs) in enumerate(faces):
            cx = sum(p[0] for p in poly)/len(poly)
            cy = sum(p[1] for p in poly)/len(poly)
            memb = tuple(sec.contains(s, mp.mpf(cx), mp.mpf(cy)) for sec in self.secs)
            ends = []
            for t in (a, b):
                q = bx
                for key in keys:
                    la, lb, lc = self.line_at(key, t)
                    sg = signs[key]
                    q = clip_f(q, la, lb, lc) if sg < 0 else clip_f(q, -la, -lb, -lc)
                    if not q:
                        break
                ends.append(q if q and abs(area_f(q)) > thr else None)
            out.append({'id': idx, 'signs': signs, 'in': memb, 'start': ends[0], 'end': ends[1]})
        return out


# ------------------------------------------------------------------ face areas

class MultiSweep:
    """One face of a prism against the two other prisms at once: the areas
    of its pieces by their classes against each (`in`, `out`, `same`,
    `opp`)."""

    def __init__(self, face, others, size):
        self.face = face
        self.sweeps = [FaceSweep(face, o, size) for o in others]
        self.size = size
        self.tol = mp.mpf(10)**-20*size

    def breakpoints(self):
        f = self.face
        lo, hi = M(f.domain[0]), M(f.domain[1])
        cand = []
        for sw in self.sweeps:
            cand += sw.breakpoints()
        found = []
        s1, s2 = self.sweeps
        for g1 in ('Oh', 'Oe', 'Ov'):
            for g2 in ('Oh', 'Oe', 'Ov'):
                for d1 in s1.groups[g1]:
                    for d2 in s2.groups[g2]:
                        poly = sq_poly(psub(d1.L, d2.L), [(d1.k, d1.Q), (-d2.k, d2.Q)])
                        if pzero(poly):
                            continue
                        for x in s1.params(poly):
                            for xx in s1.in_domain(x):
                                if self.relevant(xx, d1.key, d2.key):
                                    found.append(xx)
        self.events = len(found)
        return merge_breaks(cand+found, lo, hi, mp.mpf(10)**-30*max(1, hi-lo))

    def relevant(self, x, k1, k2):
        """A crossing of the first other and one of the second at one point
        of the line, within the face."""
        tol = self.tol
        s1, s2 = self.sweeps

        def match(key, tag):
            return tag[:len(key)] == key

        c1 = [t for t, tag, bad in s1.instances(x, tol) if match(k1, tag) and bad <= tol]
        c2 = [t for t, tag, bad in s2.instances(x, tol) if match(k2, tag) and bad <= tol]
        for t1 in c1:
            for t2 in c2:
                if abs(t1-t2) <= tol and self.face.own_bad(x, t1) <= tol:
                    return True
        return False

    def structure_at(self, x):
        """[(classes, tag0, tag1)] of the line's pieces."""
        f = self.face
        inst = []
        for k, sw in enumerate(self.sweeps):
            inst += [(t, (k,)+tag) for t, tag, bad in sw.instances(x) if bad == 0]
        X0 = f.base(x)
        out = []
        for (a, ta), (b, tb) in f.own_intervals(x):
            inner = sorted((t, tag) for t, tag in inst if a < t < b)
            pts = [(a, ta)]+inner+[(b, tb)]
            for (t0, g0), (t1, g1) in zip(pts, pts[1:]):
                tm = (t0+t1)/2
                X = tuple(X0[i]+tm*f.Em[i] for i in range(3))
                cls = tuple(sw.classify(X) for sw in self.sweeps)
                if out and out[-1][0] == cls and out[-1][2] == g0:
                    out[-1] = (cls, out[-1][1], g1)
                else:
                    out.append((cls, g0, g1))
        return out

    def values_at(self, x):
        f = self.face
        vals = {}
        for k, sw in enumerate(self.sweeps):
            for t, tag, bad in sw.instances(x):
                if bad == 0:
                    vals[(k,)+tag] = t
        for (a, ta), (b, tb) in f.own_intervals(x):
            vals[ta], vals[tb] = a, b
        return vals

    def areas(self):
        """{(class against the first other, against the second): area}."""
        f = self.face
        breaks = self.breakpoints()
        self.breaks = breaks
        total = {}
        self.quad_error = mp.mpf(0)
        tol = mp.mpf(10)**-33*self.size**2
        for a, b in zip(breaks, breaks[1:]):
            mid = (a+b)/2
            struct = self.structure_at(mid)
            if not struct:
                continue
            classes = sorted({cls for cls, _, _ in struct})

            def fn(x, struct=struct, classes=classes):
                vals = self.values_at(x)
                out = [mp.mpf(0)]*len(classes)
                m = f.metric(x)
                for cls, t0, t1 in struct:
                    if t0 not in vals or t1 not in vals:
                        raise AssertionError(f'a missed breakpoint on {cref.describe(f)} at {mp.nstr(x, 15)} '
                                             f'({mp.nstr(a, 15)}, {mp.nstr(b, 15)}): {t0} {t1}')
                    out[classes.index(cls)] += (vals[t1]-vals[t0])*m
                return out

            est, diff = integrate(fn, a, b, tol)
            self.quad_error = max(self.quad_error, diff)
            for c, v in zip(classes, est):
                total[c] = total.get(c, mp.mpf(0))+v
        return total


SIDES = {'in': (True, True), 'out': (False, False), 'same': (False, True), 'opp': (True, False)}


def bounds_chain(fn, own, others, classes):
    """Whether a piece of a face of prism `own` with the given classes
    against the prisms `others` bounds the chain: the chained function
    differs in front of the face (its prism left) and behind it (held)."""
    front, back = [None]*3, [None]*3
    front[own], back[own] = False, True
    for o, c in zip(others, classes):
        front[o], back[o] = SIDES[c]
    return fn(*front) != fn(*back)


# ------------------------------------------------------------------ the chain

class Chain:
    def __init__(self, obj, tool, third, op1, swapped, axis=None):
        self.prisms = [Prism(obj), Prism(tool), Prism(third)]
        self.op1, self.swapped = op1, swapped
        self.slicing = Slicing3(self.prisms, axis)
        self.size = self.slicing.size
        self._volumes = None
        self._areas = None
        self._solids = {}

    def fn(self, op2):
        return chained(self.op1, op2, self.swapped)

    def volumes(self):
        if self._volumes is None:
            self._volumes = self.slicing.measure(self.op1, self.swapped)
        return self._volumes

    def face_areas(self):
        """[(prism index, face, {classes: area}, sweep)]."""
        if self._areas is None:
            out = []
            for i, p in enumerate(self.prisms):
                others = [j for j in range(3) if j != i]
                for f in p.faces:
                    sweep = MultiSweep(f, [self.prisms[j] for j in others], self.size)
                    out.append((i, f, sweep.areas(), sweep))
            self._areas = out
        return self._areas

    def area(self, op2):
        fn = self.fn(op2)
        total = mp.mpf(0)
        for i, f, cls, _ in self.face_areas():
            others = [j for j in range(3) if j != i]
            for classes, a in cls.items():
                # On a face of an earlier prism: counted there.
                if any(c in ('same', 'opp') and o < i for o, c in zip(others, classes)):
                    continue
                if bounds_chain(fn, i, others, classes):
                    total += a
        return total

    def solids(self, op2):
        if op2 not in self._solids:
            self.volumes()
            self._solids[op2] = self.slicing.solids(self.fn(op2))
        return self._solids[op2]

    def result(self, op2):
        vol, mom = self.volumes()[op2]
        if vol <= mp.mpf(10)**-25*self.size**3:
            return 0, mp.mpf(0), mp.mpf(0), None
        return self.solids(op2), vol, self.area(op2), tuple(m/vol for m in mom)


def rows(chain, op2):
    n, vol, surface, centre = chain.result(op2)
    if n == 0:
        return ['empty']
    return [' '.join(['result', str(n), cref.number(vol), cref.number(surface)]+[cref.number(c) for c in centre])]
