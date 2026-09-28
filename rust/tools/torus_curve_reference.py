#!/usr/bin/env python3
"""Independent reference for S7b.3b of REVIEW_NOTES.md: the intersection of a
torus with a cylinder or a cone off its axis.

The curve is parameterised on the torus's meridians, as in S7b.3a: the
point `C(phi) + r (cos t e + sin t a)`, `C = o + R e(phi)`, in the frame whose
unit `x` is the rational direction of the other surface's axis component
normal to the torus's axis (the direction to the other's origin when the axes
are parallel). The curve is the zero set of `G(phi, t) = f(p(phi, t))` on the
flat parameter torus, `f` the other surface's implicit function (a quadric,
so `G` is a trigonometric polynomial of degree two in `t`).

Method, independent of the kernel's certified subdivision:

* every meridian's roots are those of the polynomial `z^2 G` in
  `z = e^{it}` (80-digit `polyroots`, the unit-modulus ones), its Fourier
  coefficients from five exact samples;
* the critical meridians (a fold, `G = G_t = 0`, or a tangency) are the real
  roots of the resultant of `z^2 G` and `z^2 G_t`, a trigonometric
  polynomial of degree at most 16 in `phi` found from 64 samples, kept where
  a double root lies on the circle; folds are then polished by
  two-dimensional Newton (`findroot`);
* components by continuity of the roots between meridians placed densely
  and on both sides of every critical one, joined at the folds, with each
  component's winding numbers on the torus;
* tangencies (singular points of the curve) by exact algebra (sympy): a
  torus is the pipe of radius `r` about its spine circle, a cylinder the pipe
  of radius `r_c` about its axis, and they touch where a critical distance
  between the spine and the axis is `r + r_c` or `|r - r_c|` (a Groebner
  basis of three conics in the spine's plane); the local type (an isolated
  point or two crossing branches) by the sign of the Hessian's determinant
  of `G` there. A cone is never exactly tangent for a binary64 half-angle.

Rows, in canonical order:

* `empty`;
* `fold phi t x y z` for each fold, sorted by `phi` in `(-pi, pi]`, then `t`;
* `node phi t x y z crossing|isolated` for each tangency;
* `component folds w_phi w_t` for each smooth closed component (its winding
  numbers normalised to `w_phi > 0`, or `w_t >= 0` when `w_phi = 0`), and
  `cluster folds nodes` for each component through crossings, sorted;
* `ring phi t x y z` for each component without folds, by its point at
  `phi = 0` (sorted by `t`).
"""
from fractions import Fraction as F

import mpmath as mp

import analytic_intersection_reference as ana
from analytic_intersection_reference import cross, dot, mpf, scale, sub, zero
from procedural_intersection_reference import cone_form, cylinder_form, frame, unit
from identity_reference import frame_axes

mp.mp.dps = 80


class TorusCurve:
    def __init__(self, torus, other):
        ot, A = torus.axes()
        oo, ao = other.axes()
        toward = sub(ao, scale(A, dot(ao, A)/dot(A, A)))
        self.other = other
        if zero(toward):
            w = sub(oo, ot)
            toward = sub(w, scale(A, dot(w, A)/dot(A, A)))
        self.o, self.a, self.x, self.y = frame(ot, A, toward)
        self.R, self.r = mpf(torus.radius), mpf(torus.minor)
        # G's degree in cos t, sin t: two, for a torus too (on a meridian
        # circle |p - o2|^2 is affine in cos t and sin t).
        self.deg = 2
        if other.kind == 'cylinder':
            self.f = cylinder_form(oo, ao, other.radius)['f']
        elif other.kind == 'cone':
            self.f = cone_form(oo, ao, other.radius, other.angle)['f']
        else:
            self.f = torus_form(oo, ao, other.radius, other.minor)

    def point(self, phi, t):
        c, s = mp.cos(phi), mp.sin(phi)
        k = self.R+self.r*mp.cos(t)
        return [o+k*(c*x+s*y)+self.r*mp.sin(t)*a for o, x, y, a in zip(self.o, self.x, self.y, self.a)]

    def G(self, phi, t):
        return self.f(self.point(phi, t))

    def fourier(self, phi):
        """([c_0..c_d], [s_0..s_d]) of G(phi, t) = sum c_j cos jt + s_j sin jt,
        from 2d + 1 samples."""
        n = 2*self.deg+1
        ts = [2*mp.pi*k/n for k in range(n)]
        gs = [self.G(phi, t) for t in ts]
        c = [sum(gs)/n]+[2*sum(g*mp.cos(j*t) for g, t in zip(gs, ts))/n for j in range(1, self.deg+1)]
        s = [mp.mpf(0)]+[2*sum(g*mp.sin(j*t) for g, t in zip(gs, ts))/n for j in range(1, self.deg+1)]
        return c, s

    def roots(self, phi, loose=False):
        """The real roots t in (-pi, pi] of G(phi, .), sorted (with `loose`,
        also those within 1e-12 of the circle: a double root's neighbours)."""
        coeffs = zpoly(*self.fourier(phi))
        while coeffs and abs(coeffs[0]) < mp.mpf(10)**-60:
            coeffs = coeffs[1:]
        out = []
        for z in mp.polyroots(coeffs, maxsteps=200, extraprec=200):
            if abs(abs(z)-1) < (mp.mpf(10)**-12 if loose else mp.mpf(10)**-25):
                out.append(mp.atan2(mp.im(z), mp.re(z)))
        return sorted(out)

    def Gt(self, phi, t):
        return mp.diff(lambda s: self.G(phi, s), t)

    def fold(self, phi, t):
        """A fold near (phi, t), by Newton on (G, G_t)."""
        f = lambda p, s: [self.G(p, s), self.Gt(p, s)]
        p, s = mp.findroot(f, (phi, t))
        return mp.atan2(mp.sin(p), mp.cos(p)), mp.atan2(mp.sin(s), mp.cos(s))


def zpoly(c, s):
    """z^d G in z = e^{it}, highest power first."""
    d = len(c)-1
    return [(c[j]-1j*s[j])/2 for j in range(d, 0, -1)]+[c[0]]+[(c[j]+1j*s[j])/2 for j in range(1, d+1)]


def torus_form(o, a, R, r):
    """f(p) = (|w|^2 + R^2 - r^2)^2 - 4 R^2 (|w|^2 - (w . a)^2), w = p - o."""
    om, am = [mpf(v) for v in o], unit(a)
    Rm, rm = mpf(R), mpf(r)

    def f(p):
        w = [x-y for x, y in zip(p, om)]
        ww = sum(x*x for x in w)
        h = sum(x*y for x, y in zip(w, am))
        return (ww+Rm*Rm-rm*rm)**2-4*Rm*Rm*(ww-h*h)
    return f


def wrap(x):
    return mp.atan2(mp.sin(x), mp.cos(x))


def dist(a, b):
    return abs(wrap(a-b))


def critical_phis(cv):
    """Every meridian angle where G(phi, .) has a double root on the circle
    (a fold or a tangency): the real roots of the resultant of z^d G and
    z^d G_t in z = e^{it}, a trigonometric polynomial of degree at most
    4 d^2 in phi (16 for a quadric, 64 for a torus) found from exact samples,
    kept where a double root is on the circle."""
    d = cv.deg

    def res(phi):
        c, s = cv.fourier(phi)
        a = zpoly(c, s)
        b = zpoly([0]+[j*s[j] for j in range(1, d+1)], [0]+[-j*c[j] for j in range(1, d+1)])
        n = 2*d
        rows_ = [[0]*i+a+[0]*(n-1-i) for i in range(n)]+[[0]*i+b+[0]*(n-1-i) for i in range(n)]
        return mp.det(mp.matrix(rows_))
    top = 4*d*d
    N = 4*top
    samples = [res(2*mp.pi*k/N) for k in range(N)]
    # R(phi) = sum_{j=-top}^{top} r_j e^{i j phi}.
    coeffs = []
    for j in range(top, -top-1, -1):
        coeffs.append(sum(v*mp.expj(-j*2*mp.pi*k/N) for k, v in enumerate(samples))/N)
    scale = max(abs(c) for c in coeffs)
    while abs(coeffs[0]) < mp.mpf(10)**-50*scale:
        coeffs = coeffs[1:]
    while abs(coeffs[-1]) < mp.mpf(10)**-50*scale:
        coeffs = coeffs[:-1]
    out = []
    for w in mp.polyroots(coeffs, maxsteps=400, extraprec=400):
        if abs(abs(w)-1) > mp.mpf(10)**-20:
            continue
        phi = mp.atan2(mp.im(w), mp.re(w))
        ts = cv.roots(phi, loose=True)
        close = any(dist(ts[i], ts[(i+1) % len(ts)]) < mp.mpf(10)**-12 for i in range(len(ts))) if len(ts) > 1 else False
        if close and not any(dist(phi, q) < mp.mpf(10)**-30 for q in out):
            out.append(phi)
    return sorted(out)


def scan(cv, n=720):
    """Roots over meridians placed about the critical angles (at most one
    critical angle between neighbours, at least n per turn) and the folds
    between them: a graph whose vertices are (meridian, root) and whose edges
    join a root to its continuation on the next meridian, or the two roots
    of a pair that vanishes at a fold. Returns (phis, roots, folds, edges)."""
    from itertools import combinations
    crit = critical_phis(cv)
    marks = []
    for c in crit:
        marks += [c-mp.mpf(10)**-8, c+mp.mpf(10)**-8]
    grid = [-mp.pi+2*mp.pi*(k+mp.mpf(1)/3)/n for k in range(n)]
    phis = sorted(set(grid+marks), key=float)
    phis = [p for i, p in enumerate(phis) if i == 0 or p-phis[i-1] > mp.mpf(10)**-12]
    roots = [cv.roots(p) for p in phis]
    n = len(phis)
    folds, edges = [], []
    for k in range(n):
        j = (k+1) % n
        a, b = roots[k], roots[j]
        few, many, kf, km = (a, b, k, j) if len(a) <= len(b) else (b, a, j, k)
        m, f = len(many), len(few)
        best = None
        # Match `few` to a subset of `many` in cyclic order; the rest of
        # `many` pairs up into cyclically adjacent pairs that vanish at folds.
        for chosen in combinations(range(m), f):
            rest = [i for i in range(m) if i not in chosen]
            pairs = []
            if rest:
                # A vanishing pair has no root between its two: adjacent in
                # the cyclic order of `many`.
                for start_ in (0, 1):
                    rot = rest[start_:]+rest[:start_]
                    cand = [(rot[2*i], rot[2*i+1]) for i in range(len(rot)//2)]
                    if all(y == (x+1) % m for x, y in cand):
                        pairs = cand
                        break
                else:
                    continue
            for shift in range(max(f, 1)):
                cost = sum(dist(few[(i+shift) % f], many[chosen[i]]) for i in range(f))
                cost += sum(dist(many[x], many[y]) for x, y in pairs)
                if best is None or cost < best[0]:
                    best = (cost, chosen, shift, pairs)
        _, chosen, shift, pairs = best
        for i in range(f):
            edges.append(((kf, (i+shift) % f), (km, chosen[i])))
        for x, y in pairs:
            mid = (many[x]+many[y])/2
            if abs(many[x]-many[y]) > mp.pi:
                mid += mp.pi
            folds.append(cv.fold((phis[k]+phis[j]+(2*mp.pi if j == 0 else 0))/2, mid))
            edges.append(((km, x), (km, y)))
    return phis, roots, folds, edges


def components(phis, roots, edges):
    """Cycles of the scan's graph (every vertex has degree two), each with its
    winding numbers (w_phi, w_t), normalised so that w_phi > 0, or w_t >= 0
    when w_phi = 0; and its vertices."""
    adj = {}
    for x, y in edges:
        adj.setdefault(x, []).append(y)
        adj.setdefault(y, []).append(x)
    assert all(len(v) == 2 for v in adj.values()), 'a vertex of degree other than two'
    seen, out = set(), []
    for start in adj:
        if start in seen:
            continue
        cycle, prev, cur = [start], None, start
        wphi, wt = mp.mpf(0), mp.mpf(0)
        while True:
            seen.add(cur)
            nxt = adj[cur][0] if adj[cur][0] != prev else adj[cur][1]
            if adj[cur][0] == adj[cur][1]:
                nxt = adj[cur][0]
            if nxt[0] != cur[0]:
                wphi += wrap(phis[nxt[0]]-phis[cur[0]])
            wt += wrap(roots[nxt[0]][nxt[1]]-roots[cur[0]][cur[1]])
            prev, cur = cur, nxt
            if cur == start:
                break
            cycle.append(cur)
        w = (int(mp.nint(wphi/(2*mp.pi))), int(mp.nint(wt/(2*mp.pi))))
        if w[0] < 0 or (w[0] == 0 and w[1] < 0):
            w = (-w[0], -w[1])
        out.append((w, cycle))
    return out


def rows(torus, other):
    """The canonical rows of a torus and a cylinder or cone off its axis
    (see the module's documentation)."""
    if torus.kind != 'torus':
        torus, other = other, torus
    cv = TorusCurve(torus, other)
    nodes = []
    for _, p in tangencies(torus, other):
        phi, t = angles(cv, p)
        nodes.append((phi, t, p, node_type(cv, phi, t)))
    phis, roots, folds, edges = scan(cv)
    cycles = components(phis, roots, edges) if edges else []
    if not cycles and not nodes:
        return ['empty']
    out = []
    key = lambda x: (float(x[0]), float(x[1]))
    for phi, t in sorted(folds, key=key):
        out.append(('fold', phi, t, cv.point(phi, t)))
    for phi, t, p, kind in sorted(nodes, key=key):
        out.append(('node', phi, t, p, kind))
    # Each cycle's folds (its fold edges) and nodes (its vertices beside a
    # crossing); cycles through a common node form one component.
    near = lambda a, b: dist(a[0], b[0]) < mp.mpf(10)**-6 and dist(a[1], b[1]) < mp.mpf(10)**-3
    groups = []
    for w, cycle in cycles:
        verts = set(cycle)
        nf = sum(1 for x, y in edges if x in verts and x[0] == y[0])
        at = {i for i, (phi, t, _, kind) in enumerate(nodes)
              if kind == 'crossing' and any(near((phis[k], roots[k][j]), (phi, t)) for k, j in cycle)}
        for g in groups:
            if g['nodes'] & at:
                g['nodes'] |= at
                g['folds'] += nf
                g['cycles'].append(w)
                break
        else:
            groups.append({'nodes': at, 'folds': nf, 'cycles': [w]})
    rings = []
    for g in groups:
        if g['nodes']:
            out.append(('cluster', g['folds'], len(g['nodes'])))
        else:
            (w,) = g['cycles']
            out.append(('component', g['folds'], w[0], w[1]))
    for w, cycle in cycles:
        verts = set(cycle)
        if any(x in verts and x[0] == y[0] for x, y in edges):
            continue
        # A ring: its point at phi = 0, the exact meridian's root nearest
        # its vertex on the scanned meridian nearest phi = 0.
        k0 = min(range(len(phis)), key=lambda k: abs(phis[k]))
        (j0,) = [j for k, j in cycle if k == k0]
        t0 = min(cv.roots(mp.mpf(0)), key=lambda t: dist(t, roots[k0][j0]))
        rings.append(('ring', mp.mpf(0), t0, cv.point(mp.mpf(0), t0)))
    body = [r for r in out if r[0] in ('fold', 'node')]
    comps = sorted((r for r in out if r[0] in ('component', 'cluster')), key=lambda r: (r[0], r[1:]))
    return body+comps+sorted(rings, key=lambda r: float(r[2]))


def number(x):
    x = mpf(x)
    return '0.0' if abs(x) < mp.mpf(10)**-60 else mp.nstr(x, 30, min_fixed=-5, max_fixed=5)


def text(row):
    if isinstance(row, str):
        return row
    words = [row[0]]
    for part in row[1:]:
        if isinstance(part, str):
            words.append(part)
        elif isinstance(part, int):
            words.append(str(part))
        elif isinstance(part, (list, tuple)):
            words += [number(v) for v in part]
        else:
            words.append(number(part))
    return ' '.join(words)


def tangencies(torus, other):
    """Exact tangencies of a torus and a cylinder: critical points of the
    distance between the torus's spine circle and the cylinder's axis at the
    distance r + r_c or |r - r_c|, by a Groebner basis and sympy's exact
    roots. [(spine point, contact point)]; a cone gives none; another torus:
    `torus_tangencies`."""
    import sympy as sy
    if other.kind == 'torus':
        return torus_tangencies(torus, other)
    if other.kind != 'cylinder':
        return []
    ot, A = torus.axes()
    o2, a2 = other.axes()
    R, r, rc = F(torus.radius), F(torus.minor), F(other.radius)
    # A rational basis of the spine's plane.
    U = next(u for u in (cross(A, e) for e in ((1, 0, 0), (0, 1, 0), (0, 0, 1))) if not zero(u))
    V = cross(A, U)
    u, v = sy.symbols('u v', real=True)
    Q = lambda x: sy.Rational(x.numerator, x.denominator)
    s = [Q(ot[i])+u*Q(U[i])+v*Q(V[i]) for i in range(3)]
    w = [s[i]-Q(o2[i]) for i in range(3)]
    aa = Q(dot(a2, a2))
    h = sum(w[i]*Q(a2[i]) for i in range(3))/aa
    m = [w[i]-h*Q(a2[i]) for i in range(3)]
    rel = [s[i]-Q(ot[i]) for i in range(3)]
    tangent = [Q(A[1])*rel[2]-Q(A[2])*rel[1], Q(A[2])*rel[0]-Q(A[0])*rel[2], Q(A[0])*rel[1]-Q(A[1])*rel[0]]
    E0 = sy.expand(sum(x*x for x in rel)-Q(R)**2)
    E1 = sy.expand(sum(m[i]*tangent[i] for i in range(3)))
    d2 = sy.expand(sum(x*x for x in m))
    out = []
    for k in {(r+rc)**2, (r-rc)**2}:
        basis = sy.groebner([E0, E1, d2-Q(k)], u, v, order='lex')
        if list(basis) == [1]:
            continue
        for sol in sy.solve(list(basis), [u, v], dict=True):
            su, sv = sol[u], sol[v]
            if not (su.is_real and sv.is_real):
                continue
            sm = [mp.mpf(sy.N(x.subs({u: su, v: sv}), 90)) for x in s]
            lm = [mp.mpf(sy.N((s[i]-m[i]).subs({u: su, v: sv}), 90)) for i in range(3)]
            gap = [b-a for a, b in zip(sm, lm)]
            n = mp.sqrt(sum(x*x for x in gap))
            assert n > mp.mpf(10)**-40, 'the axis meets the spine'
            for sign in (1, -1):
                p = [a+sign*mpf(r)*g/n for a, g in zip(sm, gap)]
                w2 = [x-mpf(y) for x, y in zip(p, o2)]
                hh = sum(x*mpf(y) for x, y in zip(w2, a2))/mpf(dot(a2, a2))
                rho = mp.sqrt(sum((x-hh*mpf(y))**2 for x, y in zip(w2, a2)))
                if abs(rho-mpf(rc)) < mp.mpf(10)**-40:
                    out.append((sm, p))
    return out


def angles(cv, p):
    """(phi, t) of a point of the torus."""
    rel = [x-o for x, o in zip(p, cv.o)]
    phi = mp.atan2(sum(x*y for x, y in zip(rel, cv.y)), sum(x*y for x, y in zip(rel, cv.x)))
    h = sum(x*y for x, y in zip(rel, cv.a))
    rho = mp.sqrt(sum(x*x for x in rel)-h*h)
    return phi, mp.atan2(h, rho-cv.R)


def node_type(cv, phi, t):
    """'crossing' or 'isolated' by the sign of det Hess G at a tangency."""
    g = lambda a, b: cv.G(a, b)
    gpp = mp.diff(g, (phi, t), (2, 0))
    gtt = mp.diff(g, (phi, t), (0, 2))
    gpt = mp.diff(g, (phi, t), (1, 1))
    det = gpp*gtt-gpt*gpt
    assert abs(det) > mp.mpf(10)**-30, 'a degenerate contact'
    return 'crossing' if det < 0 else 'isolated'


def curve_samples(torus, other, every=6):
    """Float points along each component of the reference curve (every
    `every`-th scanned vertex), and the isolated tangency points, for
    comparisons with sampled native lines."""
    if torus.kind != 'torus':
        torus, other = other, torus
    cv = TorusCurve(torus, other)
    phis, roots, folds, edges = scan(cv)
    cycles = components(phis, roots, edges) if edges else []
    comps = [[[float(x) for x in cv.point(phis[k], roots[k][j])] for k, j in cycle[::every]+cycle[-1:]]
             for _, cycle in cycles]
    isolated = []
    for _, p in tangencies(torus, other):
        phi, t = angles(cv, p)
        if node_type(cv, phi, t) == 'isolated':
            isolated.append([float(x) for x in p])
    return comps, isolated


def torus_tangencies(t1, t2):
    """Exact tangencies of two tori: pairs of points (s1, s2) of the spine
    circles whose segment is normal to both circles, at the distance
    r1 + r2 or |r1 - r2|: five polynomial equations in the four coordinates
    of the points in rational bases of the spines' planes, by a Groebner basis
    and sympy's exact roots. [(spine point of t1, contact point)]."""
    import sympy as sy
    o1, A1 = t1.axes()
    o2, A2 = t2.axes()
    R1, r1, R2, r2 = F(t1.radius), F(t1.minor), F(t2.radius), F(t2.minor)
    Q = lambda x: sy.Rational(x.numerator, x.denominator)

    def plane(A):
        U = next(u for u in (cross(A, e) for e in ((1, 0, 0), (0, 1, 0), (0, 0, 1))) if not zero(u))
        return U, cross(A, U)
    (U1, V1), (U2, V2) = plane(A1), plane(A2)
    u1, v1, u2, v2 = sy.symbols('u1 v1 u2 v2', real=True)
    s1 = [Q(o1[i])+u1*Q(U1[i])+v1*Q(V1[i]) for i in range(3)]
    s2 = [Q(o2[i])+u2*Q(U2[i])+v2*Q(V2[i]) for i in range(3)]
    rel1 = [s1[i]-Q(o1[i]) for i in range(3)]
    rel2 = [s2[i]-Q(o2[i]) for i in range(3)]
    crossq = lambda a, b: [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
    T1 = crossq([Q(x) for x in A1], rel1)
    T2 = crossq([Q(x) for x in A2], rel2)
    g = [s1[i]-s2[i] for i in range(3)]
    base = [sy.expand(sum(x*x for x in rel1)-Q(R1)**2), sy.expand(sum(x*x for x in rel2)-Q(R2)**2),
            sy.expand(sum(g[i]*T1[i] for i in range(3))), sy.expand(sum(g[i]*T2[i] for i in range(3)))]
    out = []
    for k in {(r1+r2)**2, (r1-r2)**2}:
        system = base+[sy.expand(sum(x*x for x in g)-Q(k))]
        basis = sy.groebner(system, u1, v1, u2, v2, order='lex')
        if list(basis) == [1]:
            continue
        for sol in sy.solve(list(basis), [u1, v1, u2, v2], dict=True):
            if not all(sol[x].is_real for x in (u1, v1, u2, v2)):
                continue
            sm = [mp.mpf(sy.N(x.subs(sol), 90)) for x in s1]
            tm = [mp.mpf(sy.N(x.subs(sol), 90)) for x in s2]
            gap = [b-a for a, b in zip(sm, tm)]
            n = mp.sqrt(sum(x*x for x in gap))
            assert n > mp.mpf(10)**-40, 'the spines meet'
            # The contact is r1 from s1 on the segment, beyond it towards s2
            # unless inside a thicker tube (k = (r1 - r2)^2 with r2 > r1).
            sign = -1 if (k != (r1+r2)**2 and r2 > r1) else 1
            p = [a+sign*mpf(r1)*x/n for a, x in zip(sm, gap)]
            out.append((sm, p))
    return out
