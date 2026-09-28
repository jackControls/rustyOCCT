#!/usr/bin/env python3
"""Independent reference for S7d.1 of REVIEW_NOTES.md: pairs of lines,
circles, ellipses and hyperbolas' branches.

Each curve is its exact point set from its stored data (`stored_axes`, the
kernel's `Frame3::new` step by step): a line `p0 + s (p1 - p0)`; a circle
the points of the plane through its origin normal to its stored normal at
its radius from the origin; an ellipse or a hyperbola the points
`o + xi x + eta y` of its stored axes with `xi^2 / A^2 +- eta^2 / B^2 = 1`
(`xi > 0` on the hyperbola's branch). Each is written as polynomial
equations in `(X, Y, Z)` (two planes for a line; a plane and a sphere for a
circle; a plane and the conic's equation in exact frame coordinates for the
others), and the two curves' equations are solved together by sympy's
Groebner basis: nothing is parameterised or eliminated as the kernel does.
Coincidence: the second curve's equations reduce to zero modulo the first's
basis (and, for hyperbolas, the branches agree). A tangency: the curves'
tangent directions (a line's direction; the cross product of a conic's two
gradients) are parallel at the point, exactly.

Rows, sorted by the first curve's parameter: `empty`, `coincident`, or
`point sa sb x y z crossing|tangent` (a line's `s`, a circle's angle from
its stored `x`, an ellipse's angle, a hyperbola's `t`).
"""
from fractions import Fraction as F

import mpmath as mp

from curve_surface_reference import stored_axes
from analytic_intersection_reference import mpf

mp.mp.dps = 80


def _Q(v):
    import sympy as sy
    v = F(v)
    return sy.Rational(v.numerator, v.denominator)


def _cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])


class Curve:
    """kind in line, circle, ellipse, hyperbola; values as the fixtures'."""

    def __init__(self, kind, values):
        self.kind, self.values = kind, tuple(values)
        if kind == 'line':
            self.p0 = tuple(F(v) for v in values[:3])
            self.d = tuple(F(b)-F(a) for a, b in zip(values[:3], values[3:6]))
        else:
            o, x, y, n = stored_axes(tuple(values[:9]))
            self.o, self.x, self.y, self.n = (tuple(F(v) for v in w) for w in (o, x, y, n))
            self.r = F(values[9])
            self.minor = F(values[10]) if kind in ('ellipse', 'hyperbola') else None

    def frame_coordinates(self, p):
        """(xi, eta, zeta) with p - o = xi x + eta y + zeta (x cross y)."""
        import sympy as sy
        m = _cross(self.x, self.y)
        M = sy.Matrix([[_Q(self.x[i]), _Q(self.y[i]), _Q(m[i])] for i in range(3)])
        rel = sy.Matrix([p[i]-_Q(self.o[i]) for i in range(3)])
        return list(M.LUsolve(rel))

    def equations(self, p):
        import sympy as sy
        if self.kind == 'line':
            d = self.d
            e = next(e for e in ((1, 0, 0), (0, 1, 0), (0, 0, 1)) if any(_cross(d, e)))
            n1 = _cross(d, e)
            n2 = _cross(d, n1)
            return [sy.expand(sum(_Q(n[i])*(p[i]-_Q(self.p0[i])) for i in range(3))) for n in (n1, n2)]
        if self.kind == 'circle':
            rel = [p[i]-_Q(self.o[i]) for i in range(3)]
            return [sy.expand(sum(_Q(self.n[i])*rel[i] for i in range(3))),
                    sy.expand(sum(c*c for c in rel)-_Q(self.r)**2)]
        xi, eta, zeta = self.frame_coordinates(p)
        sign = -1 if self.kind == 'hyperbola' else 1
        return [sy.expand(zeta), sy.expand(xi**2/_Q(self.r)**2+sign*eta**2/_Q(self.minor)**2-1)]

    def tangent(self, p, point):
        """The tangent direction at a point (sympy numbers)."""
        import sympy as sy
        if self.kind == 'line':
            return [_Q(v) for v in self.d]
        X = sy.symbols('X Y Z')
        eqs = self.equations(list(X))
        grads = [[sy.diff(e, v).subs(dict(zip(X, point))) for v in X] for e in eqs]
        return list(sy.Matrix(grads[0]).cross(sy.Matrix(grads[1])))

    def parameter(self, point):
        pm = [mp.mpf(sy_n(c)) for c in point]
        if self.kind == 'line':
            d = [mpf(v) for v in self.d]
            return sum((a-mpf(b))*c for a, b, c in zip(pm, self.p0, d))/sum(c*c for c in d)
        rel = [a-mpf(b) for a, b in zip(pm, self.o)]
        if self.kind == 'circle':
            return mp.atan2(sum(a*mpf(b) for a, b in zip(rel, self.y)),
                            sum(a*mpf(b) for a, b in zip(rel, self.x)))
        xi, eta, _ = [mp.mpf(sy_n(c)) for c in self.frame_coordinates(list(point))]
        if self.kind == 'ellipse':
            return mp.atan2(eta/mpf(self.minor), xi/mpf(self.r))
        return mp.asinh(eta/mpf(self.minor))

    def on_branch(self, point):
        if self.kind != 'hyperbola':
            return True
        return self.frame_coordinates(list(point))[0] > 0


def sy_n(c):
    import sympy as sy
    return sy.N(c, 90)


def rows(a, b):
    import sympy as sy
    X = sy.symbols('X Y Z', real=True)
    ea, eb = a.equations(list(X)), b.equations(list(X))
    basis = sy.groebner(ea, *X, order='lex', domain=sy.QQ)
    if all(basis.reduce(e)[1] == 0 for e in eb):
        # The same point set, unless two hyperbolas take opposite branches
        # (a vertex of the first on the second's branch decides).
        if a.kind == 'hyperbola' and b.kind == 'hyperbola':
            vertex = [_Q(a.o[i])+_Q(a.r)*_Q(a.x[i]) for i in range(3)]
            if not b.on_branch(vertex):
                return ['empty']
        return ['coincident']
    basis = sy.groebner(ea+eb, *X, order='lex', domain=sy.QQ)
    if list(basis) == [1]:
        return ['empty']
    out = []
    for sol in sy.solve(list(basis), list(X), dict=True):
        point = [sol[v] for v in X]
        if not all(c.is_real for c in point):
            continue
        if not (a.on_branch(point) and b.on_branch(point)):
            continue
        ta, tb = sy.Matrix(a.tangent(X, point)), sy.Matrix(b.tangent(X, point))
        tangent = all(sy.simplify(c) == 0 for c in ta.cross(tb))
        pm = tuple(mp.mpf(sy_n(c)) for c in point)
        out.append(('point', a.parameter(point), b.parameter(point), pm, 'tangent' if tangent else 'crossing'))
    out.sort(key=lambda r: float(r[1]))
    return out or ['empty']


def number(x):
    x = mpf(x)
    return '0.0' if abs(x) < mp.mpf(10)**-60 else mp.nstr(x, 30, min_fixed=-5, max_fixed=5)


def text(row):
    if isinstance(row, str):
        return row
    _, sa, sb, p, kind = row
    return ' '.join(['point', number(sa), number(sb)]+[number(v) for v in p]+[kind])
