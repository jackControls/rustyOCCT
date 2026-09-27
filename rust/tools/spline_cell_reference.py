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
    """A rational B-spline edge over its whole domain."""
    basis: Basis
    poles: list          # 3-tuples of floats
    weights: list


@dataclass
class BSpline2:
    """A planar rational B-spline pcurve over its whole domain."""
    basis: Basis
    poles: list          # 2-tuples
    weights: list


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


def _tested(basis):
    """(u for the left side, u for the right side) of each knot to test."""
    _, _, (a, b) = basis.flat()
    ks = [F(k) for k in basis.knots]
    if basis.periodic:
        # The seam: its left side is the domain's end.
        return [(b, a)]+[(k, k) for k in ks[1:-1]]
    return [(k, k) for k in ks if a < k < b]


def _one_sided(basis, ctrl, u, left):
    flat, _, _ = basis.flat()
    p = basis.degree
    if left:
        k = max(j for j in range(p, len(flat)-p-1) if flat[j] < u <= flat[j+1])
    else:
        k = min(j for j in range(p, len(flat)-p-1) if flat[j] <= u < flat[j+1])
    return _de_boor(flat, ctrl, p, k, u), _derivative(flat, ctrl, p, k, u)


def _row_c1(basis, hpoles):
    ctrl = _controls(basis, hpoles)
    for ul, ur in _tested(basis):
        if _one_sided(basis, ctrl, ul, True) != _one_sided(basis, ctrl, ur, False):
            return False
    return True


def curve_c1(curve):
    """C1 of the homogeneous curve at every tested knot."""
    return _row_c1(curve.basis, homogeneous(curve.poles, curve.weights))


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
    _, _, (a, b) = basis.flat()
    u = a if t == 0 else b
    value, _ = _one_sided(basis, ctrl, u, left=(t != 0))
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


def encode_basis(b, number):
    return (f'{b.degree} {int(b.periodic)} {len(b.knots)} '
            + ' '.join(map(number, b.knots))+' '+' '.join(map(str, b.mults)))


def encode_curve(c, number):
    """`bspline DEG PERIODIC NK knots... mults... N poles... weights...`."""
    poles = ' '.join(number(x) for p in c.poles for x in p)
    return (f'bspline {encode_basis(c.basis, number)} {len(c.poles)} {poles} '
            + ' '.join(map(number, c.weights)))


def encode_surface(s, number):
    poles = ' '.join(number(x) for p in s.poles for x in p)
    return (f'bspline {encode_basis(s.u, number)} {encode_basis(s.v, number)} {len(s.poles)} {poles} '
            + ' '.join(map(number, s.weights)))
