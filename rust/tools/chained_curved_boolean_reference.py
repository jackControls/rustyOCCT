#!/usr/bin/env python3
"""Independent reference for S9e.3a of REVIEW_NOTES.md: a Boolean's result
given to further Booleans where the solids are prisms with arcs (planar and
cylindrical walls), spheres, cones or frusta and whole tori, in any relative
position: `(A op1 B) op2 C`, swapped `C op2 (A op1 B)`, and deeper chains
(`((A op1 B) op2 C) op3 D`, any tree of Booleans over the solids).

It is S9d.4b.2's reference (`torus_curved_boolean_reference.py`: each input
its construction's exact model from the stored binary64 data, every face
covered by two families of circles and lines, each curve cut at the roots of
the other inputs' surfaces, each piece classified and measured by the
divergence theorem) generalized from two inputs and an operation to several
inputs and a set function of their memberships. Nothing here uses the
kernel, a first result, a surface/surface intersection or an arrangement:
the chain is one point set, decided at every point by its solids'
memberships.

* **Pieces and their classes.** Every face of every input is swept by its
  two families of curves (S9d.4b.2's, whose inputs, families, roots along a
  curve and quadrature this module uses). Along one curve the roots of
  every surface of every other input cut it, and each piece is classified
  at a point a golden fraction along it: each other input's membership in
  front of the face and behind it. A surface of another input vanishing on
  the whole face (a plane or a quadric equal to the face's own, where two
  inputs' faces lie on one surface) gives no roots, and that input's
  memberships are taken at the point pushed off the face along its outward
  normal by `PUSH` of the case's size both ways (the others' at the point
  itself, the same in front and behind), so faces of several inputs on one
  surface need no special case. A piece bounds the chain where the chain's set function differs
  across the face (in front the face's own input left, behind it held, the
  others as pushed): with the face's orientation where the chain holds the
  back, reversed where it holds the front. A piece where an earlier input
  (in the inputs' order) has a boundary too (its memberships differ across
  the face) is counted by that input's face, not again.
* **Measures.** Each class's area, volume and first moments as S9d.4b.2's
  (the divergence theorem with the outward normal, the volume's and
  moments' integrands in closed form along a curve, the area element by
  Gauss-Legendre), summed over the outer parameter between the family's
  events by adaptive Gauss-Legendre to 1e-32 of the case's size to the
  fourth; an event is where the curve's structure (its pieces' classes,
  adjacent equal ones merged, with the surface bounding each) changes,
  found on a scan of `SCAN` curves and located by bisection; every
  quadrature node's structure is checked against its interval's: a node of
  another structure (a feature narrower than the scan's spacing) is an
  event the scan missed, bisected against the interval's middle, and the
  interval is measured again in parts (`refined`). One sweep of every face serves every set function of the
  same solids (the first result, the chain, each input alone).
* **Solids** by rays: a grid of rays along the world's `z` through the
  chain's box, each ray's intervals inside every input from the roots of
  its surfaces in binary64 (S9d.4b.2's `FloatModel`), combined by the set
  function, intervals joined where they overlap by more than 1e-9 on
  neighbouring rays; the components counted at two resolutions (the
  fixtures declare the count, which both must give).
* **Margins.** On the scanned curves, the least sine between a face and a
  surface of another input where the curve crosses it (a tangency of
  surfaces makes it small), the least gap between two roots of different
  surfaces along a curve relative to its range (a curve through a meeting
  of two surfaces of others: an edge of one on a face), and events' least
  spacing.

`Chain(cases)` holds the solids; `result(expr)` gives `(solids, volume,
area, centre)` for an expression tree (`('cut', ('fuse', 0, 1), 2)`, leaves
the solids' indices); `rows(...)` gives `result N volume area cx cy cz` or
`empty`, as the other references.
"""
import math

import mpmath as mp

import torus_curved_boolean_reference as tc
from torus_curved_boolean_reference import (
    GOLDEN, INNER, Z, FloatModel, box_size, curve_roots, fint, fmul, lin, make_input, multiples, norm,
    normal_jet, quarter, signature)
from sphere_boolean_reference import cross, dot, nodes, scale, sub
from torus_boolean_reference import number

mp.mp.dps = 40

SCAN = 240
# The pushes off a face, relative to the case's size: far above the
# working precision's rounding of a surface's value (1e-40 of the size to
# its degree), and the measure of a sliver they could misclassify far below
# the checks' 1e-30.
PUSH = mp.mpf(10)**-34
# A surface vanishing on a face: below this at every sample, relative to the
# size to its degree.
VANISH = mp.mpf(10)**-28
# Roots of different surfaces this near along a curve, relative to its range,
# are one (`merged`), where `MERGE` is set: S9e.3b's declared tangencies only,
# whose margins are not checked (at an ordinary crossing of two roots it
# makes two events this near, the spacing margin's near coincidence).
TINY = mp.mpf(10)**-30
MERGE = False
# The quadrature's tolerance relative to the size to the fourth (S9e.3b's
# declared tangencies take a coarser one: `QUAD`).
QUAD = mp.mpf(10)**-32

SET = {
    'fuse': lambda a, b: a or b,
    'cut': lambda a, b: a and not b,
    'common': lambda a, b: a and b,
}


def evaluate(expr, mem):
    """The set function of an expression tree at the memberships `mem`."""
    if isinstance(expr, int):
        return mem[expr]
    op, x, y = expr
    return SET[op](evaluate(x, mem), evaluate(y, mem))


def leaves(expr):
    if isinstance(expr, int):
        return {expr}
    return leaves(expr[1]) | leaves(expr[2])


def chain_expr(op1, op2, swapped, more=()):
    """`(A op1 B) op2 C` (swapped `C op2 (A op1 B)`), then each `(op,
    swapped)` of `more` with the next solid."""
    expr = (op2, 2, (op1, 0, 1)) if swapped else (op2, (op1, 0, 1), 2)
    for k, (op, sw) in enumerate(more):
        expr = (op, 3+k, expr) if sw else (op, expr, 3+k)
    return expr


# ------------------------------------------------------------------ the others

class Others:
    """The inputs other than a face's own: their surfaces, and their
    memberships at a point pushed off the face both ways."""

    def __init__(self, inputs, own, size):
        self.inputs = [(j, I) for j, I in enumerate(inputs) if j != own]
        self.surfs = [(j, s) for j, I in self.inputs for s in I.surfs]
        self.push = PUSH*size

    def classify(self, X, n, pushed):
        """(front, back): each other input's membership at `X + push n` and
        `X - push n` (`n` the unit outward normal) for an input in `pushed`
        (a surface of it on the face's), at `X` for the others."""
        f = tuple(X[i]+self.push*n[i] for i in range(3))
        b = tuple(X[i]-self.push*n[i] for i in range(3))
        front, back = [], []
        for j, I in self.inputs:
            if j in pushed:
                front.append(bool(I.contains(f)))
                back.append(bool(I.contains(b)))
            else:
                m = bool(I.contains(X))
                front.append(m)
                back.append(m)
        return tuple(front), tuple(back)


def unit_normal(cur, b, sign):
    _, Xa, Xb = cur.frame_at(b)
    N = cross(Xa, Xb)
    n = norm(N)
    assert n > 0, 'a degenerate normal at a classified point'
    k = 1 if sign > 0 else -1
    return tuple(k*x/n for x in N)


def curve_pieces(cur, others, skip, sign):
    """The curve cut at every root of the other inputs' surfaces (but those
    vanishing on the face, `skip`), each piece `(b0, b1, class, bound0,
    bound1)` classified at a golden fraction along it."""
    if cur.degenerate():
        return []
    roots = []
    closed = cur.kind == 'circle' and cur.rng is None
    span = 2*mp.pi if closed else cur.rng[1]-cur.rng[0]
    for k, (j, s) in enumerate(others.surfs):
        if k in skip:
            continue
        rts = curve_roots(cur, s)
        # A double root (a tangency) split by rounding: both dropped (S9d.4b.2's).
        keep = [True]*len(rts)
        n = len(rts)
        for i in range(n):
            m = (i+1) % n
            if (m == 0 and not closed) or n < 2 or not keep[i] or not keep[m]:
                continue
            gap = rts[m]-rts[i]+(span if m == 0 else 0)
            if gap < mp.mpf(10)**-15*span:
                keep[i] = keep[m] = False
        roots += [(t, k) for t, ok in zip(rts, keep) if ok]
    roots.sort(key=lambda x: x[0])
    if MERGE:
        roots = merged(roots, span)

    pushed = {others.surfs[k][0] for k in skip}

    def cls(t0, t1):
        b = t0+(t1-t0)*GOLDEN
        return others.classify(cur.point(b), unit_normal(cur, b, sign), pushed)
    out = []
    if closed:
        T = 2*mp.pi
        if not roots:
            return [(Z, T, cls(Z, T), None, None)]
        for i, (t, j) in enumerate(roots):
            t1, j1 = roots[(i+1) % len(roots)]
            if i+1 == len(roots):
                t1 += T
            if t1 > t:
                out.append((t, t1, cls(t, t1), j, j1))
        return out
    lo, hi = cur.rng
    # A root at an end within rounding (an edge of the face on another
    # input's surface: faces of two inputs sharing an edge) cuts nothing:
    # within 1e-30 of the range, or of the push's distance along the curve.
    near = span*mp.mpf(10)**-30+others.push/norm(cur.C)
    pts = [(lo, 'end')]+[r for r in roots if lo+near < r[0] < hi-near]+[(hi, 'end')]
    for (t0, j0), (t1, j1) in zip(pts, pts[1:]):
        if t1 > t0:
            out.append((t0, t1, cls(t0, t1), j0, j1))
    return out


def merged(roots, span):
    """Sorted roots `(t, surface)` with consecutive roots of different
    surfaces within `TINY` of the curve's range taken as one, bounded by
    both (S9e.3b's declared tangencies, `MERGE`: two surfaces whose sections
    are tangent at a point the curve passes within rounding of; the piece
    between them would be classified by rounding noise, its measure far
    below the checks, and their order flip the structure, an event found
    again and again)."""
    out = []
    for t, k in roots:
        if out and t-out[-1][0] <= TINY*span and out[-1][1] != k:
            prev = out[-1][1]
            both = tuple(sorted(set(prev if isinstance(prev, tuple) else (prev,)) | {k}, key=repr))
            out[-1] = (out[-1][0], both)
        else:
            out.append((t, k))
    return out


def integrate_pieces(cur, pieces, sign, c0):
    """Per class: area, volume and first moments (about `c0`) of the pieces
    (S9d.4b.2's `integrate_pieces`, its classes any key)."""
    out = {}
    jet = normal_jet(cur, sign)
    X0 = sub(cur.X0, c0)
    C, S = cur.C, cur.S
    if cur.kind == 'circle':
        N0, Nc, Ns, N2c, N2s = jet
        Ds = [([X0[i], C[i]], [Z, S[i]]) for i in range(3)]
        Nsr = [([N0[i], Nc[i], N2c[i]], [Z, Ns[i], N2s[i]]) for i in range(3)]
        vol = None
        mom = []
        for i in range(3):
            p = fmul(Ds[i], Nsr[i])
            vol = p if vol is None else ([x+y for x, y in zip(vol[0], p[0])], [x+y for x, y in zip(vol[1], p[1])])
            mom.append(fmul(Ds[i], p))
    else:
        N0, N1 = jet
    for b0, b1, key, _, _ in pieces:
        acc = out.setdefault(key, [Z]*5)
        if cur.kind == 'circle':
            e0, e1 = multiples(b0, 4), multiples(b1, 4)
            cs = [b1-b0]+[(e1[k][0]-e0[k][0], e1[k][1]-e0[k][1]) for k in range(1, 5)]
            acc[1] += fint(vol, cs)/3
            for i in range(3):
                acc[2+i] += fint(mom[i], cs)/2
            n = max(1, int(mp.ceil((b1-b0)/quarter())))
        else:
            def pint(coef):
                return sum(c*(b1**(k+1)-b0**(k+1))/(k+1) for k, c in enumerate(coef))
            v = [Z, Z, Z]
            for i in range(3):
                x0, x1, n0, n1 = X0[i], C[i], N0[i], N1[i]
                v[0] += x0*n0
                v[1] += x0*n1+x1*n0
                v[2] += x1*n1
                sq = (x0*x0, 2*x0*x1, x1*x1)
                acc[2+i] += pint([sq[0]*n0, sq[0]*n1+sq[1]*n0, sq[1]*n1+sq[2]*n0, sq[2]*n1])/2
            acc[1] += pint(v)/3
            n = 1
        step = (b1-b0)/n
        a0 = Z
        for i in range(n):
            lo = b0+i*step
            for xi, wi in nodes(INNER):
                b = lo+step*(xi+1)/2
                if cur.kind == 'line':
                    N = (N0[0]+b*N1[0], N0[1]+b*N1[1], N0[2]+b*N1[2])
                else:
                    c, s = mp.cos(b), mp.sin(b)
                    c2, s2 = c*c-s*s, 2*c*s
                    N = tuple(N0[j]+c*Nc[j]+s*Ns[j]+c2*N2c[j]+s2*N2s[j] for j in range(3))
                a0 += wi*mp.sqrt(N[0]*N[0]+N[1]*N[1]+N[2]*N[2])
        acc[0] += a0*step/2
    return out


# ------------------------------------------------------------------ a family swept

class Unconverged(Exception):
    pass


def capped(f, a, b, tol, depth=0):
    """S9d.1's adaptive Gauss-Legendre (`integrate`), its halving capped at
    40 levels (`Unconverged`)."""
    prev = None
    half = (b-a)/2
    for degree in (3, 4, 5, 6):
        est = None
        for xi, wi in nodes(degree):
            th = (xi+1)*mp.pi/2
            s = a+half*(1-mp.cos(th))
            jac = wi*half*mp.sin(th)*mp.pi/2
            v = f(s)
            est = [jac*c for c in v] if est is None else [e+jac*c for e, c in zip(est, v)]
        if prev is not None:
            diff = max(abs(x-y) for x, y in zip(est, prev))
            if diff <= tol:
                return est, diff
        prev = est
    if depth >= 40:
        raise Unconverged()
    m = (a+b)/2
    e1, d1 = capped(f, a, m, tol*3/5, depth+1)
    e2, d2 = capped(f, m, b, tol*3/5, depth+1)
    return [x+y for x, y in zip(e1, e2)], d1+d2


class Sweep:
    """A face's family against the other inputs: its events, its own kinks,
    each class's measures integrated between them (S9d.4b.2's `Sweep`)."""

    def __init__(self, fam, others, c0, size):
        self.fam, self.others, self.c0, self.size = fam, others, c0, size
        self.missed = 0
        self.refined = 0
        self.quad_error = Z
        self.events = None
        self.samples = []
        self.skip = self.vanishing()

    def vanishing(self):
        """The other inputs' surfaces vanishing on the face (sampled on three
        curves at three points each)."""
        fam = self.fam
        out = set()
        curves = []
        for f in (mp.mpf(1)/7, mp.mpf(1)/2, mp.mpf(5)/7):
            cur = fam.curve(fam.a0+(fam.a1-fam.a0)*f)
            if cur.degenerate():
                continue
            if cur.kind == 'circle' and cur.rng is None:
                bs = [mp.mpf(1), mp.mpf(3), mp.mpf(5)]
            else:
                lo, hi = cur.rng
                bs = [lo+(hi-lo)*g for g in (mp.mpf(1)/5, mp.mpf(1)/2, mp.mpf(4)/5)]
            curves.append((cur, bs))
        for k, (_, s) in enumerate(self.others.surfs):
            tol = VANISH*self.size**s.degree
            if curves and all(abs(s.f(cur.point(b))) <= tol for cur, bs in curves for b in bs):
                out.add(k)
        return out

    def structure(self, a):
        cur = self.fam.curve(a)
        pcs = curve_pieces(cur, self.others, self.skip, self.fam.sign)
        return signature(cur, pcs), cur, pcs

    def find_events(self, keep_samples=False):
        fam = self.fam
        a0, a1 = fam.a0, fam.a1
        span = a1-a0
        eps = span*mp.mpf(10)**-35
        if fam.periodic:
            grid = [a0+span*(j+mp.mpf(0.3819660112501051))/SCAN for j in range(SCAN)]
        else:
            delta = span*mp.mpf(10)**-25
            grid = [a0+delta]+[a0+span*j/SCAN for j in range(1, SCAN)]+[a1-delta]
        for k in fam.own:
            grid += [k-span*mp.mpf(10)**-25, k+span*mp.mpf(10)**-25]
        grid.sort()
        sigs = []
        for a in grid:
            s, cur, pcs = self.structure(a)
            sigs.append(s)
            if keep_samples:
                self.samples.append((a, cur, pcs))
        found = []
        bisect = lambda lo, slo, hi, shi: self.bisect(lo, slo, hi, shi, found)
        n = len(grid)
        pairs = list(zip(range(n-1), range(1, n)))
        if fam.periodic:
            pairs.append((n-1, 0))
        for i, j in pairs:
            if sigs[i] != sigs[j]:
                lo, hi = grid[i], grid[j]
                if j == 0:
                    hi += span
                bisect(lo, sigs[i], hi, sigs[j])
        out = []
        for e in found:
            if fam.periodic:
                e = a0+(e-a0) % span
            out.append(e)
        self.events = sorted(out)
        return self.events

    def bisect(self, lo, slo, hi, shi, found):
        """The events between two outer parameters of different structures
        (S9d.4b.2's bisection), appended to `found`."""
        span = self.fam.a1-self.fam.a0
        if hi-lo <= span*mp.mpf(10)**-35:
            found.append((lo+hi)/2)
            return
        m = (lo+hi)/2
        sm = self.structure(m)[0]
        if sm != slo and sm != shi and hi-lo <= span*mp.mpf(10)**-20:
            # Two surfaces tangent along an edge (S9d.4b.2's note).
            found.append(m)
            return
        if sm != slo:
            self.bisect(lo, slo, m, sm, found)
        if sm != shi:
            self.bisect(m, sm, hi, shi, found)

    def breakpoints(self):
        fam = self.fam
        pts = sorted(set([fam.a0, fam.a1]+list(self.events)+list(fam.own)))
        return [p for p in pts if fam.a0 <= p <= fam.a1]

    def measure(self):
        """{class: [area, volume, moments about c0]} of the face."""
        if self.events is None:
            self.find_events()
        tol = QUAD*self.size**4
        tot = {}
        pts = self.breakpoints()
        k = 0
        while k < len(pts)-1:
            lo, hi = pts[k], pts[k+1]
            k += 1
            if hi-lo <= (self.fam.a1-self.fam.a0)*mp.mpf(10)**-34:
                continue
            mid = (lo+hi)/2
            ref, _, pcs = self.structure(mid)
            keys = sorted({p[2] for p in pcs}, key=repr)
            if not keys:
                continue
            off = []

            def f(a, ref=ref, keys=keys, off=off):
                s, cur, pcs = self.structure(a)
                if s != ref:
                    off.append((a, s))
                m = integrate_pieces(cur, pcs, self.fam.sign, self.c0)
                out = []
                for key in keys:
                    out += m.get(key, [Z]*5)
                return out
            try:
                est, err = capped(f, lo, hi, tol)
            except Unconverged:
                est = None
            if off or est is None:
                # An event the scan missed (a feature narrower than its
                # spacing: two events close together whose outer
                # structures agree): bisected between a node of another
                # structure and the interval's middle, and the interval
                # measured again in parts.
                found = []
                for a, s in off[:4]:
                    if a < mid:
                        self.bisect(a, s, mid, ref, found)
                    else:
                        self.bisect(mid, ref, a, s, found)
                new = sorted(set(x for x in found if lo < x < hi))
                assert new, 'quadrature did not converge'
                self.events = sorted(set(self.events)|set(new))
                self.refined += len(new)
                pts = pts[:k]+new+pts[k:]
                k -= 1
                continue
            self.quad_error += err
            for i, key in enumerate(keys):
                acc = tot.setdefault(key, [Z]*5)
                for m in range(5):
                    acc[m] += est[5*i+m]
        return tot

    def margins(self):
        """From the scanned curves: the least sine between the face and a
        surface it crosses, and the least gap between roots of different
        surfaces relative to the curve's range."""
        sine, gap = mp.mpf(1), mp.mpf(1)
        for a, cur, pcs in self.samples:
            if not pcs:
                continue
            closed = cur.kind == 'circle' and cur.rng is None
            span = 2*mp.pi if closed else cur.rng[1]-cur.rng[0]
            bounds = []
            for b0, _, _, j0, _ in pcs:
                if isinstance(j0, int):
                    bounds.append((b0, j0))
            for b, k in bounds:
                _, s = self.others.surfs[k]
                X = cur.point(b)
                n = unit_normal(cur, b, 1)
                g = s.grad(X)
                gn = norm(g)
                if gn == 0:
                    sine = Z
                    continue
                sine = min(sine, norm(cross(n, g))/gn)
            bounds.sort(key=lambda x: x[0])
            for (b0, k0), (b1, k1) in zip(bounds, bounds[1:]):
                if k0 != k1:
                    gap = min(gap, (b1-b0)/span)
        return sine, gap


# ------------------------------------------------------------------ the chain

class Chain:
    """Solids (prisms of one convex profile of segments and arcs, whole
    spheres, cones or frusta, whole tori) and the sweeps of all their
    faces."""

    def __init__(self, cases):
        self.inputs = [make_input(c) for c in cases]
        boxes = [I.box() for I in self.inputs]
        self.size = mp.mpf(box_size(boxes))
        lo = [min(b[0][i] for b in boxes) for i in range(3)]
        hi = [max(b[1][i] for b in boxes) for i in range(3)]
        self.c0 = tuple(mp.mpf((lo[i]+hi[i])/2) for i in range(3))
        self._sweeps = {}

    def sweeps(self, k):
        """[(input, face name, Sweep)] of family `k`, measured."""
        if k not in self._sweeps:
            out = []
            for i, I in enumerate(self.inputs):
                others = Others(self.inputs, i, self.size)
                for face in I.faces:
                    sw = Sweep(face.families[k], others, self.c0, self.size)
                    sw.find_events(keep_samples=(k == 0))
                    sw.result = sw.measure()
                    out.append((i, face.name, sw))
            self._sweeps[k] = out
        return self._sweeps[k]

    def contribution(self, i, key, expr, used):
        """+1, -1 or 0: how a piece of input `i`'s face of class `key` bounds
        the expression's set (its orientation, reversed, or not at all)."""
        if i not in used:
            return 0
        front, back = key
        others = [j for j in range(len(self.inputs)) if j != i]
        for idx, j in enumerate(others):
            if j < i and j in used and front[idx] != back[idx]:
                return 0
        fm = [False]*len(self.inputs)
        bm = [False]*len(self.inputs)
        fm[i], bm[i] = False, True
        for idx, j in enumerate(others):
            fm[j], bm[j] = front[idx], back[idx]
        f, b = evaluate(expr, fm), evaluate(expr, bm)
        if b and not f:
            return 1
        if f and not b:
            return -1
        return 0

    def measures(self, expr, k=0):
        """(volume, world moments, area) of an expression's set by family `k`."""
        used = leaves(expr)
        tot = [Z]*5
        for i, _, sw in self.sweeps(k):
            for key, vec in sw.result.items():
                c = self.contribution(i, key, expr, used)
                if c == 0:
                    continue
                tot[0] += vec[0]
                for m in range(1, 5):
                    tot[m] += c*vec[m]
        V = tot[1]
        return V, tuple(self.c0[i]*V+tot[2+i] for i in range(3)), tot[0]

    def face_measures(self, expr, k=0):
        """{(input, face): area bounding the expression's set} by family `k`."""
        used = leaves(expr)
        out = {}
        for i, name, sw in self.sweeps(k):
            a = Z
            for key, vec in sw.result.items():
                if self.contribution(i, key, expr, used) != 0:
                    a += vec[0]
            out[(i, name)] = a
        return out

    def result(self, expr, solids):
        """(solids, volume, area, centre) with the declared solid count."""
        V, m, A = self.measures(expr)
        if solids == 0:
            return 0, Z, Z, (Z, Z, Z)
        return solids, V, A, tuple(x/V for x in m)

    def missed(self):
        return sum(sw.missed for k in self._sweeps for _, _, sw in self._sweeps[k])

    def refined(self):
        return sum(sw.refined for k in self._sweeps for _, _, sw in self._sweeps[k])

    def quad_error(self):
        return max([Z]+[sw.quad_error for k in self._sweeps for _, _, sw in self._sweeps[k]])

    def margins(self):
        """The least sine and root gap over the first family's scans, and the
        events' least spacing relative to their family's range."""
        sine, gap, spacing = mp.mpf(1), mp.mpf(1), mp.mpf(1)
        for _, _, sw in self.sweeps(0):
            s, g = sw.margins()
            sine, gap = min(sine, s), min(gap, g)
        for k in self._sweeps:
            for _, _, sw in self._sweeps[k]:
                span = sw.fam.a1-sw.fam.a0
                pts = sorted(set(sw.events))
                for a, b in zip(pts, pts[1:]):
                    spacing = min(spacing, (b-a)/span)
        return sine, gap, spacing


# ------------------------------------------------------------------ solids by rays

class Float(FloatModel):
    """S9d.4b.2's binary64 input, a prism's lower cap at its own height
    (`FloatModel` puts it at height zero, where every prism of S9d.4b.2's
    fixtures starts)."""

    def polys(self, P, D):
        out = super().polys(P, D)
        if self.kind == 'prism':
            p0, p1 = self.chart(P), self.chart_vec(D)
            out[0] = [p0[2]-self.lo, p1[2]]
        return out


def count_solids(chain, expr, n):
    """Components of the expression's set by a grid of `n` by `n` rays along
    the world's `z` through the chain's box (binary64), intervals joined where
    they overlap by more than 1e-9 on neighbouring rays."""
    used = sorted(leaves(expr))
    fms = {j: Float(chain.inputs[j]) for j in used}
    boxes = [chain.inputs[j].box() for j in used]
    lo = [min(b[0][i] for b in boxes) for i in range(3)]
    hi = [max(b[1][i] for b in boxes) for i in range(3)]
    pad = 1e-3*max(hi[i]-lo[i] for i in range(3))
    lo = [x-pad for x in lo]
    hi = [x+pad for x in hi]
    tmax = hi[2]-lo[2]
    # Grid offsets by irrational fractions: no ray along a face.
    xs = [lo[0]+(hi[0]-lo[0])*(i+0.3819660112501051)/n for i in range(n)]
    ys = [lo[1]+(hi[1]-lo[1])*(j+0.6180339887498949)/n for j in range(n)]
    rays = {}
    mem = [False]*len(chain.inputs)
    for i, x in enumerate(xs):
        for j, y in enumerate(ys):
            P, D = [x, y, lo[2]], [0.0, 0.0, 1.0]
            ivs = {k: fms[k].intervals(P, D, tmax) for k in used}
            pts = sorted(set([0.0, tmax]+[t for k in used for iv in ivs[k] for t in iv]))
            out = []
            for a, b in zip(pts, pts[1:]):
                m = (a+b)/2
                for k in used:
                    mem[k] = any(p < m < q for p, q in ivs[k])
                if evaluate(expr, mem):
                    if out and out[-1][1] == a:
                        out[-1] = (out[-1][0], b)
                    else:
                        out.append((a, b))
            rays[(i, j)] = out
    parent = {}

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a
    for (i, j), ivs in rays.items():
        for k in range(len(ivs)):
            parent[(i, j, k)] = (i, j, k)
    for (i, j), ivs in rays.items():
        for di, dj in ((1, 0), (0, 1)):
            nb = rays.get((i+di, j+dj))
            if not nb:
                continue
            for k, (a, b) in enumerate(ivs):
                for kk, (a2, b2) in enumerate(nb):
                    if min(b, b2)-max(a, a2) > 1e-9:
                        ra, rb = find((i, j, k)), find((i+di, j+dj, kk))
                        if ra != rb:
                            parent[ra] = rb
    return len({find(a) for a in parent})


def rows(chain, expr, solids):
    """The fixture rows of an expression with its declared solid count."""
    n, V, A, c = chain.result(expr, solids)
    if n == 0:
        return ['empty']
    return ['result {} {} {} {} {} {}'.format(n, number(V), number(A), *(number(x) for x in c))]


def monte_carlo(chain, exprs, samples, seed):
    """Monte-Carlo volumes and centres of several expressions over the
    chain's box (binary64 memberships): {expr: (volume, centre, standard
    error of the volume, of each centre coordinate)}."""
    import random
    rnd = random.Random(seed)
    fms = [Float(I) for I in chain.inputs]
    boxes = [I.box() for I in chain.inputs]
    lo = [min(b[0][i] for b in boxes) for i in range(3)]
    hi = [max(b[1][i] for b in boxes) for i in range(3)]
    vol = 1.0
    for i in range(3):
        vol *= hi[i]-lo[i]
    acc = {e: [0, [0.0]*3, [0.0]*3] for e in exprs}
    for _ in range(samples):
        X = [lo[i]+(hi[i]-lo[i])*rnd.random() for i in range(3)]
        mem = [fm.contains(X) for fm in fms]
        for e in exprs:
            if evaluate(e, mem):
                a = acc[e]
                a[0] += 1
                for i in range(3):
                    a[1][i] += X[i]
                    a[2][i] += X[i]*X[i]
    out = {}
    for e, (k, s, s2) in acc.items():
        p = k/samples
        V = vol*p
        se_v = vol*math.sqrt(max(p*(1-p), 0.0)/samples)
        if k == 0:
            out[e] = (V, None, se_v, None)
            continue
        c = [s[i]/k for i in range(3)]
        se_c = [math.sqrt(max(s2[i]/k-c[i]*c[i], 0.0)/k) for i in range(3)]
        out[e] = (V, c, se_v, se_c)
    return out


__all__ = ['Chain', 'chain_expr', 'count_solids', 'evaluate', 'leaves', 'monte_carlo', 'rows', 'tc']
