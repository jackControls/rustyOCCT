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

Rows, sorted by the curve parameter:

* `empty`, or `contained`;
* `point s x y z crossing|tangent` (a line's `s`, a circle's angle in
  `(-pi, pi]`).
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
    _, s, p, kind = row
    return ' '.join(['point', number(s)]+[number(v) for v in p]+[kind])
