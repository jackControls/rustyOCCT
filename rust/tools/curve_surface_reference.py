#!/usr/bin/env python3
"""Independent reference for S7c.1 of REVIEW_NOTES.md: the topology's lines
and circles against planes, cylinders, cones, spheres and tori.

A line is the exact set through its two rational points, parameterised
`p0 + s (p1 - p0)`; a circle the exact point set of its stored frame (the
plane through the origin normal to the stored normal) and radius,
parameterised by the angle from the stored frame's x axis (`frame_axes`, the
kernel's `Frame3::new`). A surface is the exact point set of S7a/S7b.

* A line against a plane, cylinder, sphere or torus: the surface's function
  along the line is a polynomial in `s` with rational coefficients (degree
  one, two or four), fitted exactly from rational samples; its real roots
  and their multiplicities by sympy (`real_roots`): a double root is a
  tangency, the zero polynomial containment.
* A circle against them: in a rational basis `(U, V)` of the circle's plane,
  its points solve `|u U + v V|^2 = r^2` and `f(o + u U + v V) = 0`, both
  rational; sympy's Groebner basis gives the solutions, a tangency where the
  two curves' gradients are parallel there (an exact zero determinant), and
  containment where `f` vanishes on the whole circle.
* A cone's `cos` and `sin` are transcendental: its function along a line is
  a quadratic and along a circle a trigonometric polynomial of degree two,
  both with 80-digit mpmath coefficients and `polyroots` (no exact
  tangency).

S7c.2 adds ellipses and hyperbolas of stored frames (`conic_rows`: their
rational parameterisations, exact, 80 digits against a cone) and rational
B-splines against tori and cones (`spline_rows`: exact span polynomials).

Rows, sorted by the curve parameter:

* `empty`, or `contained`;
* `point s x y z crossing|tangent` (a line's `s`, a circle's or an
  ellipse's angle in `(-pi, pi]`, a hyperbola's or a spline's parameter);
* `overlap s0 s1` (a spline's spans on the surface).
"""
from fractions import Fraction as F

import mpmath as mp

import analytic_intersection_reference as ana
from analytic_intersection_reference import cross, dot, mpf, scale, sub, zero
from identity_reference import frame_axes

mp.mp.dps = 80


def implicit(surface):
    """f(p) for a surface: exact on rationals for planes, cylinders, spheres
    and tori; an mpmath function for cones."""
    o, a = surface.axes()
    k = surface.kind
    if k == 'plane':
        return lambda p: dot(sub(p, o), a)
    if k == 'sphere':
        r2 = F(surface.radius)**2
        return lambda p: dot(sub(p, o), sub(p, o))-r2
    if k == 'cylinder':
        r2 = F(surface.radius)**2
        aa = dot(a, a)
        return lambda p: dot(sub(p, o), sub(p, o))-dot(sub(p, o), a)**2/aa-r2
    if k == 'torus':
        R2, r2 = F(surface.radius)**2, F(surface.minor)**2
        aa = dot(a, a)

        def f(p):
            w = sub(p, o)
            ww = dot(w, w)
            return (ww+R2-r2)**2-4*R2*(ww-dot(w, a)**2/aa)
        return f
    return ana.quadric(surface)


def exact(surface):
    return surface.kind != 'cone'


def line_rows(p0, p1, surface):
    import sympy as sy
    d = sub(p1, p0)
    f = implicit(surface)
    at = lambda s: tuple(x+s*y for x, y in zip(p0, d))
    point = lambda s: tuple(mpf(x)+mpf(s)*mpf(y) for x, y in zip(p0, d))
    if exact(surface):
        # Degree at most four: interpolate at five rational parameters.
        xs = [F(k) for k in range(-2, 3)]
        s = sy.symbols('s')
        poly = sy.interpolate([(sy.Rational(x.numerator, x.denominator),
                                sy.Rational(f(at(x)).numerator, f(at(x)).denominator)) for x in xs], s)
        poly = sy.Poly(sy.expand(poly), s)
        if poly.is_zero:
            return ['contained']
        rows = []
        roots = sy.real_roots(poly)
        for r in sorted(set(roots), key=lambda x: float(x)):
            mult = roots.count(r)
            sm = mp.mpf(sy.N(r, 90))
            rows.append(('point', sm, point(sm), 'tangent' if mult > 1 else 'crossing'))
        return rows or ['empty']
    # A cone: the quadratic A s^2 + 2 B s + C from three samples.
    g = lambda s: f(point(s))
    c0, cp, cm = g(0), g(1), g(-1)
    A, B = (cp+cm)/2-c0, (cp-cm)/4
    rows = []
    if abs(A) < mp.mpf(10)**-60:
        if abs(B) < mp.mpf(10)**-60:
            return ['empty'] if abs(c0) > mp.mpf(10)**-60 else ['contained']
        s0 = -c0/(2*B)
        return [('point', s0, point(s0), 'crossing')]
    disc = B*B-A*c0
    if disc < -mp.mpf(10)**-60:
        return ['empty']
    assert abs(disc) > mp.mpf(10)**-60, 'a cone is never exactly tangent to a line here'
    for s0 in sorted([(-B-mp.sqrt(disc))/A, (-B+mp.sqrt(disc))/A]):
        rows.append(('point', s0, point(s0), 'crossing'))
    return rows


def circle_rows(frame, radius, surface):
    import sympy as sy
    o, x, y, n = frame_axes(frame)
    om = [mpf(v) for v in o]
    xm, ym = [mpf(v) for v in x], [mpf(v) for v in y]
    rm = mpf(radius)
    angle = lambda p: mp.atan2(sum((a-b)*c for a, b, c in zip(p, om, ym)),
                               sum((a-b)*c for a, b, c in zip(p, om, xm)))
    point = lambda t: tuple(a+rm*(mp.cos(t)*b+mp.sin(t)*c) for a, b, c in zip(om, xm, ym))
    if exact(surface):
        Q = lambda v: sy.Rational(v.numerator, v.denominator)
        on = tuple(F(v) for v in o)
        nn = tuple(F(v) for v in n)
        U = next(u for u in (cross(nn, e) for e in ((1, 0, 0), (0, 1, 0), (0, 0, 1))) if not zero(u))
        V = cross(nn, U)
        u, v = sy.symbols('u v', real=True)
        p = [Q(on[i])+u*Q(U[i])+v*Q(V[i]) for i in range(3)]
        rel = [p[i]-Q(on[i]) for i in range(3)]
        E0 = sy.expand(sum(c*c for c in rel)-Q(F(radius))**2)
        f = implicit(surface)
        # f with sympy arithmetic: rebuild from the surface's exact data.
        E1 = sy.expand(sympy_implicit(surface, p))
        rem = sy.reduced(E1, [E0], v, u)[1] if E1 != 0 else 0
        # Contained: f vanishes on the whole circle (a remainder of zero
        # modulo E0 as a polynomial in v).
        if E1 == 0 or sy.expand(rem) == 0:
            return ['contained']
        basis = sy.groebner([E0, E1], u, v, order='lex')
        if list(basis) == [1]:
            return ['empty']
        rows = []
        for sol in sy.solve(list(basis), [u, v], dict=True):
            su, sv = sol[u], sol[v]
            if not (su.is_real and sv.is_real):
                continue
            jac = sy.Matrix([[sy.diff(E0, u), sy.diff(E0, v)], [sy.diff(E1, u), sy.diff(E1, v)]]).det()
            tangent = sy.simplify(jac.subs({u: su, v: sv})) == 0
            pt = tuple(mp.mpf(sy.N(c.subs({u: su, v: sv}), 90)) for c in p)
            rows.append(('point', angle(pt), pt, 'tangent' if tangent else 'crossing'))
        rows.sort(key=lambda r: float(r[1]))
        return rows or ['empty']
    # A cone: a trigonometric polynomial of degree two in the angle.
    f = implicit(surface)
    g = lambda t: f(point(t))
    ts = [2*mp.pi*k/5 for k in range(5)]
    gs = [g(t) for t in ts]
    c = [sum(gs)/5]+[2*sum(q*mp.cos(j*t) for q, t in zip(gs, ts))/5 for j in (1, 2)]
    s = [0]+[2*sum(q*mp.sin(j*t) for q, t in zip(gs, ts))/5 for j in (1, 2)]
    coeffs = [(c[2]-1j*s[2])/2, (c[1]-1j*s[1])/2, c[0], (c[1]+1j*s[1])/2, (c[2]+1j*s[2])/2]
    while coeffs and abs(coeffs[0]) < mp.mpf(10)**-60:
        coeffs = coeffs[1:]
    if not coeffs or all(abs(q) < mp.mpf(10)**-60 for q in coeffs):
        return ['contained']
    rows = []
    for z in mp.polyroots(coeffs, maxsteps=200, extraprec=200):
        if abs(abs(z)-1) < mp.mpf(10)**-25:
            t = mp.atan2(mp.im(z), mp.re(z))
            rows.append(('point', t, point(t), 'crossing'))
    rows.sort(key=lambda r: float(r[1]))
    return rows or ['empty']


def conic_rows(frame, major, minor, hyperbola, surface):
    """An ellipse `o + A cos t x + B sin t y` or a hyperbola's branch
    `o + A cosh t x + B sinh t y` (the stored axes `x`, `y`) against a
    surface: exactly by the conic's rational parameterisation substituted
    into the surface's function (planes, cylinders, spheres, tori; sympy's
    `real_roots` with multiplicities), in 80 digits for a cone (the
    ellipse's trigonometric polynomial of degree two; the hyperbola's quartic
    in `z = e^t`)."""
    import sympy as sy
    o, x, y, _ = frame_axes(frame)
    A, B = F(major), F(minor)
    om, xm, ym = [mpf(v) for v in o], [mpf(v) for v in x], [mpf(v) for v in y]
    Am, Bm = mpf(A), mpf(B)
    if hyperbola:
        point = lambda t: tuple(a+Am*mp.cosh(t)*b+Bm*mp.sinh(t)*c for a, b, c in zip(om, xm, ym))
        param = lambda xi, eta: mp.asinh(eta/Bm)
    else:
        point = lambda t: tuple(a+Am*mp.cos(t)*b+Bm*mp.sin(t)*c for a, b, c in zip(om, xm, ym))
        param = lambda xi, eta: mp.atan2(eta/Bm, xi/Am)
    if exact(surface):
        # The rational parameterisation `s = tan(t / 2)` (`tanh(t / 2)` on the
        # hyperbola's branch, `|s| < 1`): the surface's function times
        # `(1 +- s^2)^k` is a polynomial in `s` with rational coefficients,
        # its real roots and multiplicities exact; the ellipse's point
        # `s = oo` (`t = pi`) is checked apart.
        Q = lambda v: sy.Rational(v.numerator, v.denominator)
        sv = sy.symbols('s', real=True)
        den = 1-sv**2 if hyperbola else 1+sv**2
        xi = Q(A)*(1+sv**2 if hyperbola else 1-sv**2)
        eta = Q(B)*2*sv
        p = [Q(F(o[i]))*den+xi*Q(F(x[i]))+eta*Q(F(y[i])) for i in range(3)]
        f = sympy_implicit(surface, [c/den for c in p])
        # The function's degree in the point: its product with `den^k` is a
        # polynomial of degree at most `2 k`.
        k = {'plane': 1, 'sphere': 2, 'cylinder': 2, 'torus': 4}[surface.kind]
        num = sy.expand(sy.cancel(f*den**k))
        if num == 0:
            return ['contained']
        poly = sy.Poly(num, sv)
        rows = []
        roots = sy.real_roots(poly)
        for r in sorted(set(roots), key=lambda z: float(z)):
            if hyperbola and not (-1 < r < 1):
                continue
            rm = mp.mpf(sy.N(r, 90))
            t = 2*mp.atanh(rm) if hyperbola else 2*mp.atan(rm)
            rows.append(('point', t, point(t), 'tangent' if roots.count(r) > 1 else 'crossing'))
        if not hyperbola:
            # s = oo: the degree drop below 2 k is the order of the root at
            # infinity.
            drop = 2*k-poly.degree()
            if drop > 0:
                rows.append(('point', mp.pi, point(mp.pi), 'tangent' if drop > 1 else 'crossing'))
        rows.sort(key=lambda r: float(r[1]))
        return rows or ['empty']
    f = implicit(surface)
    tiny = mp.mpf(10)**-60
    if not hyperbola:
        g = lambda t: f(point(t))
        ts = [2*mp.pi*k/5 for k in range(5)]
        gs = [g(t) for t in ts]
        c = [sum(gs)/5]+[2*sum(q*mp.cos(j*t) for q, t in zip(gs, ts))/5 for j in (1, 2)]
        s = [0]+[2*sum(q*mp.sin(j*t) for q, t in zip(gs, ts))/5 for j in (1, 2)]
        coeffs = [(c[2]-1j*s[2])/2, (c[1]-1j*s[1])/2, c[0], (c[1]+1j*s[1])/2, (c[2]+1j*s[2])/2]
        roots = [mp.atan2(mp.im(z), mp.re(z)) for z in mp.polyroots(coeffs, maxsteps=200, extraprec=200)
                 if abs(abs(z)-1) < mp.mpf(10)**-25]
    else:
        # 4 z^2 G, G at t = ln z: a quartic in z, interpolated at five points.
        zs = [mp.mpf(k) for k in (1, 2, 3, 4, 5)]
        vals = [4*z*z*f(point(mp.log(z))) for z in zs]
        M = mp.matrix([[z**j for j in range(5)] for z in zs])
        cs = mp.lu_solve(M, mp.matrix(vals))
        coeffs = [cs[j] for j in range(4, -1, -1)]
        while coeffs and abs(coeffs[0]) < tiny:
            coeffs = coeffs[1:]
        roots = [mp.log(mp.re(z)) for z in mp.polyroots(coeffs, maxsteps=200, extraprec=200)
                 if abs(mp.im(z)) < mp.mpf(10)**-25 and mp.re(z) > 0]
    rows = [('point', t, point(t), 'crossing') for t in roots]
    rows.sort(key=lambda r: float(r[1]))
    return rows or ['empty']


def spline_spans(spline):
    """Each nonzero span `(a, b, H)` of a clamped rational B-spline `(degree,
    poles (x, y, z, w), knots, mults)`, `H` its homogeneous coordinates
    `(w x, w y, w z, w)` as exact sympy polynomials in the curve parameter,
    interpolated from exact de Boor evaluations."""
    import sympy as sy
    degree, poles, knots, mults = spline
    flat = [F(k) for k, m in zip(knots, mults) for _ in range(m)]
    hom = [tuple(F(c)*F(p[3]) for c in p[:3])+(F(p[3]),) for p in poles]

    def de_boor(k, t):
        d = [hom[j+k-degree] for j in range(degree+1)]
        for r in range(1, degree+1):
            for j in range(degree, r-1, -1):
                i = j+k-degree
                a = (t-flat[i])/(flat[i+degree+1-r]-flat[i])
                d[j] = tuple((1-a)*x+a*y for x, y in zip(d[j-1], d[j]))
        return d[degree]
    s = sy.symbols('s')
    out = []
    for k in range(degree, len(flat)-degree-1):
        a, b = flat[k], flat[k+1]
        if a == b:
            continue
        ts = [a+(b-a)*F(j, degree) if degree else a for j in range(degree+1)]
        pts = [de_boor(k, t) for t in ts]
        H = [sy.expand(sy.interpolate([(sy.Rational(t.numerator, t.denominator),
                                        sy.Rational(p[c].numerator, p[c].denominator)) for t, p in zip(ts, pts)], s))
             if degree else sy.Rational(pts[0][c].numerator, pts[0][c].denominator) for c in range(4)]
        out.append((a, b, H))
    return s, out


def spline_rows(spline, surface):
    """A rational B-spline against a torus (exact: the homogeneous quartic on
    each span, sympy's `real_roots` with multiplicities) or a cone (`W^2 F =
    Q0 + tau Q1 + tau^2 Q2`, exact `Q`s, 80-digit `tau` and `polyroots`; a
    span on the cone exactly when all three vanish, their common roots exact
    tangencies). Rows: `overlap s0 s1` for maximal runs of spans on the
    surface, `point s x y z crossing|tangent` elsewhere, by parameter."""
    import sympy as sy
    s, spans = spline_spans(spline)
    Q = lambda v: sy.Rational(v.numerator, v.denominator)
    o, a = surface.axes()
    o, a = [Q(v) for v in o], [Q(v) for v in a]
    aa = sum(c*c for c in a)
    rows, overlaps = [], []
    for lo, hi, H in spans:
        W = H[3]
        D = [sy.expand(H[i]-o[i]*W) for i in range(3)]
        DD = sy.expand(sum(c*c for c in D))
        Da = sy.expand(sum(D[i]*a[i] for i in range(3)))
        if surface.kind == 'torus':
            R2, r2 = Q(F(surface.radius))**2, Q(F(surface.minor))**2
            polys = [sy.expand((DD+(R2-r2)*W**2)**2-4*R2*W**2*(DD-Da**2/aa))]
        else:
            r = Q(F(surface.radius))
            polys = [sy.expand(DD-Da**2/aa-r**2*W**2), sy.expand(-2*r*Da*W), sy.expand(-Da**2)]
        if all(p == 0 for p in polys):
            if overlaps and overlaps[-1][1] == lo:
                overlaps[-1][1] = hi
            else:
                overlaps.append([lo, hi])
            continue
        inside = lambda x: Q(lo) <= x <= Q(hi)
        point = lambda sm: tuple(mp.mpf(sy.N(H[i].subs(s, sm), 90))/mp.mpf(sy.N(W.subs(s, sm), 90)) for i in range(3))
        if surface.kind == 'torus':
            poly = sy.Poly(polys[0], s)
            roots = sy.real_roots(poly)
            for r in sorted(set(roots), key=lambda x: float(x)):
                if inside(r):
                    rows.append(('point', r, point(r), 'tangent' if roots.count(r) > 1 else 'crossing'))
            continue
        common = sy.gcd(sy.gcd(polys[0], polys[1]), polys[2])
        exact_roots = [r for r in set(sy.real_roots(sy.Poly(common, s))) if inside(r)] if sy.Poly(common, s).degree() > 0 else []
        for r in exact_roots:
            rows.append(('point', r, point(r), 'tangent'))
        rest = [sy.quo(p, common, s) if common != 1 else p for p in polys]
        tau = mp.tan(mpf(surface.angle))/mp.sqrt(mpf(F(aa.p, aa.q)))
        coeffs = [mp.mpf(0)]*(1+max(sy.Poly(p, s).degree() for p in rest))
        for k, p in enumerate(rest):
            for (e,), c in sy.Poly(p, s).terms():
                coeffs[e] += tau**k*mp.mpf(sy.N(c, 90))
        coeffs = coeffs[::-1]
        while coeffs and abs(coeffs[0]) < mp.mpf(10)**-60:
            coeffs = coeffs[1:]
        if len(coeffs) > 1:
            for z in mp.polyroots(coeffs, maxsteps=400, extraprec=400):
                if abs(mp.im(z)) < mp.mpf(10)**-25 and mpf(lo) <= mp.re(z) <= mpf(hi):
                    t = mp.re(z)
                    rows.append(('point', t, point(t), 'crossing'))
    out = [('overlap', F(a), F(b)) for a, b in overlaps]
    # Points inside an overlap (its knots) belong to it.
    for r in rows:
        t = mp.mpf(sy.N(r[1], 90)) if not isinstance(r[1], mp.mpf) else r[1]
        if any(mpf(a) <= t <= mpf(b) for a, b in overlaps):
            continue
        if any(abs(t-mp.mpf(sy.N(q[1], 90)) if not isinstance(q[1], mp.mpf) else abs(t-q[1])) < mp.mpf(10)**-40
               for q in out if q[0] == 'point'):
            continue
        out.append(('point', t, r[2], r[3]))
    out.sort(key=lambda r: float(r[1]))
    return out or ['empty']


def sympy_implicit(surface, p):
    """The surface's function at a sympy point, from its exact data."""
    import sympy as sy
    Q = lambda v: sy.Rational(v.numerator, v.denominator)
    o, a = surface.axes()
    o, a = [Q(v) for v in o], [Q(v) for v in a]
    w = [p[i]-o[i] for i in range(3)]
    ww = sum(c*c for c in w)
    aa = sum(c*c for c in a)
    wa = sum(w[i]*a[i] for i in range(3))
    k = surface.kind
    if k == 'plane':
        return wa
    if k == 'sphere':
        return ww-Q(F(surface.radius))**2
    if k == 'cylinder':
        return ww-wa**2/aa-Q(F(surface.radius))**2
    R2, r2 = Q(F(surface.radius))**2, Q(F(surface.minor))**2
    return (ww+R2-r2)**2-4*R2*(ww-wa**2/aa)


def number(x):
    x = mpf(x)
    return '0.0' if abs(x) < mp.mpf(10)**-60 else mp.nstr(x, 30, min_fixed=-5, max_fixed=5)


def text(row):
    if isinstance(row, str):
        return row
    if row[0] == 'overlap':
        return ' '.join(['overlap', number(row[1]), number(row[2])])
    _, s, p, kind = row
    return ' '.join(['point', number(s)]+[number(v) for v in p]+[kind])
