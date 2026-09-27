"""Independent reference for spline cells (R4 of REVIEW_NOTES.md).

Rational B-spline edges, pcurves and faces in `Fraction`s. A cell is C1 in
its own parameterisation when, at every knot strictly inside a nonperiodic
domain (every knot, the seam included, of a periodic one), the homogeneous
curve's one-sided first derivatives agree; a surface when every control row
across each knot line is. The kernel decides the same property by one
exact knot removal; this module never removes a knot, it differentiates
each side by de Boor's recurrence.

The data follow OCCT's (and the kernel's) conventions, which define the
curve: nonperiodic multiplicities sum to poles + degree + 1 and the domain
is [flat[p], flat[poles]]; a periodic basis repeats its first and last
multiplicities, has sum(mults) - mults[0] poles used cyclically, and its
flat sequence is extended by degree + 1 - mults[0] knots on each side from
the neighbouring period.
"""
from dataclasses import dataclass
from fractions import Fraction as F


@dataclass
class Basis:
    degree: int
    knots: list          # distinct, increasing floats
    mults: list
    periodic: bool = False

    def flat(self):
        """(flat knots as Fractions, pole count, (a, b) domain)."""
        p = self.degree
        flat = [F(k) for k, m in zip(self.knots, self.mults) for _ in range(m)]
        if self.periodic:
            poles = sum(self.mults)-self.mults[0]
            period = F(self.knots[-1])-F(self.knots[0])
            ext = p+1-self.mults[0]
            before = [k-period for k in flat[poles-ext:poles]]
            after = [k+period for k in flat[self.mults[0]:self.mults[0]+ext]]
            return before+flat+after, poles, (F(self.knots[0]), F(self.knots[-1]))
        poles = sum(self.mults)-p-1
        return flat, poles, (flat[p], flat[poles])


@dataclass
class BSpline3:
    """A rational B-spline edge over a range (its whole domain when None)."""
    basis: Basis
    poles: list          # 3-tuples of floats
    weights: list
    range: tuple = None
    reversed: bool = False


@dataclass
class BSpline2:
    """A planar rational B-spline pcurve over a range (its whole domain when
    None)."""
    basis: Basis
    poles: list          # 2-tuples
    weights: list
    range: tuple = None
    reversed: bool = False


def flip(curve):
    """The same span traversed the other way (exact: no knot moves)."""
    import dataclasses
    return dataclasses.replace(curve, reversed=not curve.reversed)


def _dir(curve, t):
    """A span fraction as a fraction of the range in the curve's direction."""
    return 1-t if curve.reversed else t


def span_of(curve):
    """The curve's range as Fractions: the fraction maps affinely onto it."""
    if curve.range is not None:
        return F(curve.range[0]), F(curve.range[1])
    return curve.basis.flat()[2]


def closed_period(curve):
    """A periodic curve over a whole period (a ring)."""
    a, e = span_of(curve)
    _, _, (da, de) = curve.basis.flat()
    return curve.basis.periodic and e-a == de-da


def _reduce(basis, u):
    """A periodic parameter moved into the base period [a, e] (e itself
    kept, so its left side is the period's end)."""
    if not basis.periodic:
        return u
    _, _, (a, e) = basis.flat()
    period = e-a
    while u > e:
        u -= period
    while u < a:
        u += period
    return u


@dataclass
class BSplineSurface:
    """A tensor-product rational B-spline surface; poles U-major."""
    u: Basis
    v: Basis
    poles: list          # 3-tuples, index u * v_count + v
    weights: list


def homogeneous(poles, weights):
    return [tuple(F(w)*F(x) for x in p)+(F(w),) for p, w in zip(poles, weights)]


def _de_boor(flat, ctrl, p, k, u):
    """Value on span k (flat[k] <= u <= flat[k+1]) of the degree-p spline
    whose controls for basis j are ctrl[j], by de Boor's recurrence."""
    d = {j: ctrl[j] for j in range(k-p, k+1)}
    for r in range(1, p+1):
        for j in range(k, k-p+r-1, -1):
            a = (u-flat[j])/(flat[j+p-r+1]-flat[j])
            d[j] = tuple((1-a)*x+a*y for x, y in zip(d[j-1], d[j]))
    return d[k]


def _derivative(flat, ctrl, p, k, u):
    """First derivative on span k: the degree p-1 spline of the controls
    p (P_j - P_{j-1}) / (t_{j+p} - t_j)."""
    q = {j: tuple(p*(x-y)/(flat[j+p]-flat[j]) for x, y in zip(ctrl[j], ctrl[j-1]))
         for j in range(k-p+1, k+1)}
    if p == 1:
        return q[k]
    return _de_boor(flat, q, p-1, k, u)


def _controls(basis, hpoles):
    """Controls by basis index over the (extended) flat sequence."""
    flat, n, _ = basis.flat()
    count = len(flat)-basis.degree-1
    return {j: hpoles[j % n if basis.periodic else j] for j in range(count)}


def _tested(basis, lo=None, hi=None, closed=None):
    """(u for the left side, u for the right side) of each knot to test:
    those strictly inside the range [lo, hi] (the domain by default), and on
    a closed period the one at its ends. Periodic knots count at every
    translate, evaluated in the base period (the seam's left side is the
    period's end)."""
    _, _, (a, b) = basis.flat()
    lo = a if lo is None else lo
    hi = b if hi is None else hi
    closed = basis.periodic and hi-lo == b-a if closed is None else closed
    ks = [F(k) for k in basis.knots]
    if not basis.periodic:
        return [(k, k) for k in ks if lo < k < hi]
    out = []
    period = b-a
    for k in ks[:-1]:
        t = k-((k-lo)//period)*period
        while t <= hi:
            if lo < t < hi or (closed and t in (lo, hi)):
                out.append((b, a) if k == a else (k, k))
                break
            t += period
    return out


def _one_sided(basis, ctrl, u, left):
    flat, _, _ = basis.flat()
    p = basis.degree
    if left:
        k = max(j for j in range(p, len(flat)-p-1) if flat[j] < u <= flat[j+1])
    else:
        k = min(j for j in range(p, len(flat)-p-1) if flat[j] <= u < flat[j+1])
    return _de_boor(flat, ctrl, p, k, u), _derivative(flat, ctrl, p, k, u)


def _row_c1(basis, hpoles, lo=None, hi=None):
    ctrl = _controls(basis, hpoles)
    for ul, ur in _tested(basis, lo, hi):
        if _one_sided(basis, ctrl, ul, True) != _one_sided(basis, ctrl, ur, False):
            return False
    return True


def curve_c1(curve):
    """C1 of the homogeneous curve at every tested knot."""
    lo, hi = span_of(curve)
    return _row_c1(curve.basis, homogeneous(curve.poles, curve.weights), lo, hi)


def surface_c1(s):
    """C1 across every knot line: each control row across it is C1."""
    h = homogeneous(s.poles, s.weights)
    nu = len(s.poles)//_pole_count(s.v)
    nv = _pole_count(s.v)
    rows_u = [[h[i*nv+j] for i in range(nu)] for j in range(nv)]
    rows_v = [[h[i*nv+j] for j in range(nv)] for i in range(nu)]
    return all(_row_c1(s.u, r) for r in rows_u) and all(_row_c1(s.v, r) for r in rows_v)


def _pole_count(basis):
    return basis.flat()[1]


def end_point(curve, t):
    """The exact point at fraction 0 or 1 (the domain's ends)."""
    basis = curve.basis
    ctrl = _controls(basis, homogeneous(curve.poles, curve.weights))
    a, b = span_of(curve)
    t = _dir(curve, t)
    u = _reduce(basis, a if t == 0 else b)
    _, _, (_, e) = basis.flat()
    value, _ = _one_sided(basis, ctrl, u, left=(u == e or (t != 0 and not basis.periodic)))
    w = value[-1]
    return tuple(x/w for x in value[:-1])


def curve_valid(curve, tol):
    """Degenerate when every pole lies within tolerance of the first."""
    first = [F(x) for x in curve.poles[0]]
    t2 = F(float(tol))**2  # the binary64 tolerance, exactly
    return any(sum((F(x)-y)**2 for x, y in zip(p, first)) > t2 for p in curve.poles)


def pcurve_valid(curve):
    """Degenerate when every pole is the first."""
    return any(tuple(p) != tuple(curve.poles[0]) for p in curve.poles)


def reversed_curve(c):
    """The same curve traversed backwards over the mirrored domain: knots
    k -> a + b - k (exact in binary64 only when the sums are), poles and
    weights reversed, the direction flag kept. Returns None when a mirrored
    knot is not exact (native rows need a curve; the cells use `flip`)."""
    b = c.basis
    a, e = F(b.knots[0]), F(b.knots[-1])
    knots = [a+e-F(k) for k in reversed(b.knots)]
    if any(float(k) != k for k in knots):
        return None
    basis = Basis(b.degree, [float(k) for k in knots], list(reversed(b.mults)), b.periodic)
    rng = None
    if c.range is not None:
        lo, hi = a+e-F(c.range[1]), a+e-F(c.range[0])
        if float(lo) != lo or float(hi) != hi:
            return None
        rng = (float(lo), float(hi))
    return type(c)(basis, list(reversed(c.poles)), list(reversed(c.weights)), rng, c.reversed)


def unflagged(c):
    """The same traversal as an unflagged curve: a flagged span is its
    mirrored curve run forward. None when a mirrored knot is not exact."""
    import dataclasses
    if not c.reversed:
        return c
    m = reversed_curve(c)
    return None if m is None else dataclasses.replace(m, reversed=False)


def reparameterized(c, first, last):
    """The same curve with its domain moved affinely onto [first, last]; the
    second value says whether every knot moved exactly."""
    b = c.basis
    a, e = span_of(c)
    lo, hi = F(first), F(last)
    exact_knots = [lo+(F(k)-a)*(hi-lo)/(e-a) for k in b.knots]
    knots = [float(k) for k in exact_knots]
    exact = all(F(k) == x for k, x in zip(knots, exact_knots))
    return type(c)(Basis(b.degree, knots, list(b.mults), b.periodic), list(c.poles), list(c.weights)), exact


def native_basis(b, number):
    """The native rows' basis: `DEG PERIODIC NK knots... mults...`."""
    return encode_basis(b, number)


def encode_basis(b, number):
    return (f'{b.degree} {int(b.periodic)} {len(b.knots)} '
            + ' '.join(map(number, b.knots))+' '+' '.join(map(str, b.mults)))


def encode_curve(c, number, with_range=True):
    """`bspline DEG PERIODIC NK knots... mults... N poles... weights...`,
    then `range FIRST LAST` for a range other than the whole domain."""
    poles = ' '.join(number(x) for p in c.poles for x in p)
    tail = f' range {number(c.range[0])} {number(c.range[1])}' if with_range and c.range is not None else ''
    tail += ' reversed' if with_range and c.reversed else ''
    return (f'bspline {encode_basis(c.basis, number)} {len(c.poles)} {poles} '
            + ' '.join(map(number, c.weights))+tail)


def encode_surface(s, number):
    poles = ' '.join(number(x) for p in s.poles for x in p)
    return (f'bspline {encode_basis(s.u, number)} {encode_basis(s.v, number)} {len(s.poles)} {poles} '
            + ' '.join(map(number, s.weights)))


# ---------------------------------------------------------------- evaluation (S4b-d)

def fraction(x):
    """An exact Fraction of a float, Fraction or mpmath number."""
    if isinstance(x, F):
        return x
    if isinstance(x, (int, float)):
        return F(x)
    import mpmath as mp
    man, exp = mp.mpf(x).man_exp
    return F(int(man))*F(2)**exp if exp >= 0 else F(int(man), 2**(-exp))


def _span(flat, p, u, left):
    if left:
        return max(j for j in range(p, len(flat)-p-1) if flat[j] < u <= flat[j+1])
    return min(j for j in range(p, len(flat)-p-1) if flat[j] <= u < flat[j+1])


def curve_jet(curve, t):
    """(point, derivative) in the curve's fraction t, exactly: the parameter
    is a + t (b - a), the derivative by t."""
    b = curve.basis
    flat, _, (_, de) = b.flat()
    a, e = span_of(curve)
    t = _dir(curve, fraction(t))
    sign = -1 if curve.reversed else 1
    u = _reduce(b, a+t*(e-a))
    ctrl = _controls(b, homogeneous(curve.poles, curve.weights))
    k = _span(flat, b.degree, u, left=(u == de))
    h = _de_boor(flat, ctrl, b.degree, k, u)
    dh = _derivative(flat, ctrl, b.degree, k, u)
    w, dw = h[-1], dh[-1]
    point = tuple(x/w for x in h[:-1])
    speed = tuple(sign*(dx*w-x*dw)/(w*w)*(e-a) for x, dx in zip(h[:-1], dh[:-1]))
    return point, speed


def surface_jet(s, u, v):
    """(point, S_u, S_v) exactly at parameters (u, v) of a spline surface."""
    u, v = fraction(u), fraction(v)
    h = homogeneous(s.poles, s.weights)
    nv = _pole_count(s.v)
    nu = len(s.poles)//nv
    fu, _, (ua, ub) = s.u.flat()
    fv, _, (va, vb) = s.v.flat()
    ku = _span(fu, s.u.degree, u, left=(u == ub))
    kv = _span(fv, s.v.degree, v, left=(v == vb))
    cu = _controls(s.u, list(range(nu)))
    cv = _controls(s.v, list(range(nv)))
    # Along v for every active u control, then along u.
    rows, drows = {}, {}
    for j in range(ku-s.u.degree, ku+1):
        i = cu[j]
        ctrl = {m: h[i*nv+cv[m]] for m in range(kv-s.v.degree, kv+1)}
        rows[j] = _de_boor(fv, ctrl, s.v.degree, kv, v)
        drows[j] = _derivative(fv, ctrl, s.v.degree, kv, v)
    x = _de_boor(fu, rows, s.u.degree, ku, u)
    xu = _derivative(fu, rows, s.u.degree, ku, u)
    xv = _de_boor(fu, drows, s.u.degree, ku, u)
    w = x[-1]
    point = tuple(c/w for c in x[:-1])
    su = tuple((d*w-c*xu[-1])/(w*w) for c, d in zip(x[:-1], xu[:-1]))
    sv = tuple((d*w-c*xv[-1])/(w*w) for c, d in zip(x[:-1], xv[:-1]))
    return point, su, sv


def blossom_piece(curve, t0, t1):
    """The homogeneous Bézier control points of the curve on fractions
    [t0, t1] inside one span, by blossoming: control i is the blossom at
    p - i copies of u0 and i copies of u1."""
    b = curve.basis
    flat, _, (_, de) = b.flat()
    a, e = span_of(curve)
    p = b.degree
    if curve.reversed:
        t0, t1 = 1-fraction(t1), 1-fraction(t0)
    u0, u1 = a+fraction(t0)*(e-a), a+fraction(t1)*(e-a)
    shift = _reduce(b, u0)-u0
    u0, u1 = u0+shift, u1+shift
    k = _span(flat, p, u0, left=False)
    ctrl = _controls(b, homogeneous(curve.poles, curve.weights))
    out = []
    for i in range(p+1):
        args = [u0]*(p-i)+[u1]*i
        d = {j: ctrl[j] for j in range(k-p, k+1)}
        for r in range(1, p+1):
            for j in range(k, k-p+r-1, -1):
                al = (args[r-1]-flat[j])/(flat[j+p-r+1]-flat[j])
                d[j] = tuple((1-al)*x+al*y for x, y in zip(d[j-1], d[j]))
        out.append(d[k])
    return out[::-1] if curve.reversed else out


def knot_fractions(curve):
    """The curve's distinct knots strictly inside its range, as fractions
    (periodic knots at every translate)."""
    a, e = span_of(curve)
    _, _, (da, de) = curve.basis.flat()
    out = set()
    for k in curve.basis.knots:
        k = F(k)
        if curve.basis.periodic:
            period = de-da
            k -= ((k-a)//period)*period
            while k < e:
                if a < k:
                    out.add((k-a)/(e-a))
                k += period
        elif a < k < e:
            out.add((k-a)/(e-a))
    return sorted(1-f for f in out) if curve.reversed else sorted(out)


def patch_boxes(s):
    """The surface's knot spans as boxes ((u0, u1), (v0, v1))."""
    def spans(b):
        _, _, (a, e) = b.flat()
        ks = [a]+[F(k) for k in b.knots if a < F(k) < e]+[e]
        return list(zip(ks, ks[1:]))
    return [(su, sv) for su in spans(s.u) for sv in spans(s.v)]


def curve_jet_mp(curve, t):
    """(point, derivative in t) in mpmath arithmetic, for quadrature."""
    import mpmath as mp
    b = curve.basis
    flat, _, (da, de) = b.flat()
    flat = [mp.mpf(k.numerator)/k.denominator for k in flat]
    a, e = span_of(curve)
    a, e = mp.mpf(a.numerator)/a.denominator, mp.mpf(e.numerator)/e.denominator
    da, de = mp.mpf(da.numerator)/da.denominator, mp.mpf(de.numerator)/de.denominator
    sign = -1 if curve.reversed else 1
    u = a+(1-mp.mpf(t) if curve.reversed else mp.mpf(t))*(e-a)
    if b.periodic:
        while u > de:
            u -= de-da
        while u < da:
            u += de-da
    e_domain = de
    ctrl = _controls(b, [tuple(mp.mpf(w)*mp.mpf(x) for x in p)+(mp.mpf(w),)
                         for p, w in zip(curve.poles, curve.weights)])
    p = b.degree
    k = max(j for j in range(p, len(flat)-p-1) if flat[j] <= u) if u < e_domain else \
        max(j for j in range(p, len(flat)-p-1) if flat[j] < u)
    k = min(k, len(flat)-p-2)
    h = _de_boor(flat, ctrl, p, k, u)
    dh = _derivative(flat, ctrl, p, k, u)
    w, dw = h[-1], dh[-1]
    return (tuple(x/w for x in h[:-1]),
            tuple(sign*(dx*w-x*dw)/(w*w)*(e-a) for x, dx in zip(h[:-1], dh[:-1])))


def surface_jet_mp(s, u, v):
    """(point, S_u, S_v) in mpmath arithmetic, for quadrature."""
    import mpmath as mp
    m = lambda x: mp.mpf(x.numerator)/x.denominator if isinstance(x, F) else mp.mpf(x)
    u, v = m(u), m(v)
    nv = _pole_count(s.v)
    nu = len(s.poles)//nv
    h = [tuple(mp.mpf(w)*mp.mpf(x) for x in p)+(mp.mpf(w),) for p, w in zip(s.poles, s.weights)]
    fu, _, (ua, ub) = s.u.flat()
    fv, _, (va, vb) = s.v.flat()
    fu, fv = [m(k) for k in fu], [m(k) for k in fv]
    ua, ub, va, vb = m(ua), m(ub), m(va), m(vb)

    def span(flat, p, x, end):
        k = max(j for j in range(p, len(flat)-p-1) if (flat[j] <= x if x < end else flat[j] < x))
        return min(k, len(flat)-p-2)
    ku, kv = span(fu, s.u.degree, u, ub), span(fv, s.v.degree, v, vb)
    cu, cv = _controls(s.u, list(range(nu))), _controls(s.v, list(range(nv)))
    rows, drows = {}, {}
    for j in range(ku-s.u.degree, ku+1):
        ctrl = {q: h[cu[j]*nv+cv[q]] for q in range(kv-s.v.degree, kv+1)}
        rows[j] = _de_boor(fv, ctrl, s.v.degree, kv, v)
        drows[j] = _derivative(fv, ctrl, s.v.degree, kv, v)
    x = _de_boor(fu, rows, s.u.degree, ku, u)
    xu = _derivative(fu, rows, s.u.degree, ku, u)
    xv = _de_boor(fu, drows, s.u.degree, ku, u)
    w = x[-1]
    return (tuple(c/w for c in x[:-1]),
            tuple((d*w-c*xu[-1])/(w*w) for c, d in zip(x[:-1], xu[:-1])),
            tuple((d*w-c*xv[-1])/(w*w) for c, d in zip(x[:-1], xv[:-1])))


# ---------------------------------------------------------------- interval jets (S4b Taylor path)

def _iv(x):
    import mpmath as mp
    if isinstance(x, F):
        return mp.iv.mpf(x.numerator)/x.denominator
    return mp.iv.mpf(x)


def spline_jet2_iv(curve, t):
    """Value, first and second derivative in the fraction t (an mpmath
    interval inside one span) of a spline, by interval de Boor on the
    homogeneous controls and the quotient rule."""
    b = curve.basis
    flat, _, (_, de) = b.flat()
    a, e = span_of(curve)
    p = b.degree
    fl = [_iv(k) for k in flat]
    sign = -1 if curve.reversed else 1
    if curve.reversed:
        import mpmath as _mp
        t = 1-t
    mid = a+fraction(t.mid)*(e-a)
    shift = _reduce(b, mid)-mid
    a, e = a+shift, e+shift
    a_, e_ = _iv(a), _iv(e)
    u = a_+t*(e_-a_)
    mid = mid+shift
    k = _span(flat, p, mid, left=(mid == de))
    ctrl = _controls(b, [tuple(_iv(w)*_iv(x) for x in q)+(_iv(w),) for q, w in zip(curve.poles, curve.weights)])
    h = _de_boor(fl, ctrl, p, k, u)
    if p >= 1:
        q1 = {j: tuple(p*(x-y)/(fl[j+p]-fl[j]) for x, y in zip(ctrl[j], ctrl[j-1])) for j in range(k-p+1, k+1)}
        h1 = q1[k] if p == 1 else _de_boor(fl, q1, p-1, k, u)
    if p >= 2:
        q2 = {j: tuple((p-1)*(x-y)/(fl[j+p-1]-fl[j]) for x, y in zip(q1[j], q1[j-1])) for j in range(k-p+2, k+1)}
        h2 = q2[k] if p == 2 else _de_boor(fl, q2, p-2, k, u)
    else:
        h2 = tuple(0*x for x in h)
    s = (e_-a_)*sign
    w, w1, w2 = h[-1], h1[-1]*s, h2[-1]*s*s
    v = [x/w for x in h[:-1]]
    v1 = [(x*s-y*w1)/w for x, y in zip(h1[:-1], v)]
    v2 = [(x*s*s-2*y*w1-z*w2)/w for x, y, z in zip(h2[:-1], v1, v)]
    return v, v1, v2


def surface_jet2_iv(s, u, v):
    """S, S_u, S_v, S_uu, S_uv, S_vv of a spline surface over an interval box
    (u, v), the hull over the knot spans the box meets."""
    import mpmath as mp
    nv = _pole_count(s.v)
    nu = len(s.poles)//nv
    fu, _, (ua, ub) = s.u.flat()
    fv, _, (va, vb) = s.v.flat()
    cu, cv = _controls(s.u, list(range(nu))), _controls(s.v, list(range(nv)))
    h = [tuple(_iv(w)*_iv(x) for x in q)+(_iv(w),) for q, w in zip(s.poles, s.weights)]
    flu, flv = [_iv(k) for k in fu], [_iv(k) for k in fv]

    def spans(flat, p, lo, hi):
        return [j for j in range(p, len(flat)-p-1) if flat[j] < flat[j+1] and flat[j] <= hi and flat[j+1] >= lo]

    def derivs(flat, fl, p, ctrl, k, x):
        d0 = _de_boor(fl, ctrl, p, k, x)
        q1 = {j: tuple(p*(a-b)/(fl[j+p]-fl[j]) for a, b in zip(ctrl[j], ctrl[j-1])) for j in range(k-p+1, k+1)}
        d1 = q1[k] if p == 1 else _de_boor(fl, q1, p-1, k, x)
        if p >= 2:
            q2 = {j: tuple((p-1)*(a-b)/(fl[j+p-1]-fl[j]) for a, b in zip(q1[j], q1[j-1])) for j in range(k-p+2, k+1)}
            d2 = q2[k] if p == 2 else _de_boor(fl, q2, p-2, k, x)
        else:
            d2 = tuple(0*a for a in d0)
        return d0, d1, d2
    lo_u, hi_u = fraction(u.a), fraction(u.b)
    lo_v, hi_v = fraction(v.a), fraction(v.b)
    out = None
    for ku in spans(fu, s.u.degree, lo_u, hi_u):
        for kv in spans(fv, s.v.degree, lo_v, hi_v):
            uu = mp.iv.mpf([max(fu[ku], lo_u), min(fu[ku+1], hi_u)]) if False else \
                mp.iv.mpf([_iv(max(fu[ku], lo_u)).a, _iv(min(fu[ku+1], hi_u)).b])
            vv = mp.iv.mpf([_iv(max(fv[kv], lo_v)).a, _iv(min(fv[kv+1], hi_v)).b])
            rows = {}
            for j in range(ku-s.u.degree, ku+1):
                ctrl = {q: h[cu[j]*nv+cv[q]] for q in range(kv-s.v.degree, kv+1)}
                rows[j] = derivs(fv, flv, s.v.degree, ctrl, kv, vv)
            # Along u for each v-derivative order.
            byv = [{j: rows[j][o] for j in rows} for o in range(3)]
            x0 = derivs(fu, flu, s.u.degree, byv[0], ku, uu)   # H, H_u, H_uu
            x1 = derivs(fu, flu, s.u.degree, byv[1], ku, uu)   # H_v, H_uv
            x2 = _de_boor(flu, byv[2], s.u.degree, ku, uu)     # H_vv
            H, Hu, Huu = x0
            Hv, Huv = x1[0], x1[1]
            Hvv = x2
            w, wu, wv, wuu, wuv, wvv = H[-1], Hu[-1], Hv[-1], Huu[-1], Huv[-1], Hvv[-1]
            S = [c/w for c in H[:-1]]
            Su = [(c-a*wu)/w for c, a in zip(Hu[:-1], S)]
            Sv = [(c-a*wv)/w for c, a in zip(Hv[:-1], S)]
            Suu = [(c-2*b*wu-a*wuu)/w for c, a, b in zip(Huu[:-1], S, Su)]
            Suv = [(c-b*wv-d*wu-a*wuv)/w for c, a, b, d in zip(Huv[:-1], S, Su, Sv)]
            Svv = [(c-2*d*wv-a*wvv)/w for c, a, d in zip(Hvv[:-1], S, Sv)]
            jet = [S, Su, Sv, Suu, Suv, Svv]
            if out is None:
                out = jet
            else:
                out = [[mp.iv.mpf([min(x.a, y.a), max(x.b, y.b)]) for x, y in zip(p, q)] for p, q in zip(out, jet)]
    return out
