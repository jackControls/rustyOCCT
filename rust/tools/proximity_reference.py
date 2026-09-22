"""Independent Fraction distance formulas and convex optimality certificates.

The reference uses analytic projections, cross products, and boundary candidates;
it does not use the kernel's face enumeration or Gram-system elimination.
"""
from fractions import Fraction as Q
import math
import struct

COUNTS = {'P': 1, 'L': 2, 'S': 2, 'F': 3, 'T': 3}


def parse(row, hexadecimal=False):
    words = iter(row.split())
    name = next(words)
    shapes = []
    for _ in range(2):
        kind = next(words)
        def number():
            token = next(words)
            return struct.unpack('>d', bytes.fromhex(token))[0] if hexadecimal else float(token)
        shapes.append((kind, [tuple(number() for _ in range(3)) for _ in range(COUNTS[kind])]))
    return name, shapes, list(words)


def exact(shape):
    kind, points = shape
    return kind, [tuple(Q.from_float(x) for x in p) for p in points]


def sub(a, b): return tuple(x-y for x, y in zip(a, b))
def add(a, b): return tuple(x+y for x, y in zip(a, b))
def mul(a, t): return tuple(x*t for x in a)
def dot(a, b): return sum(x*y for x, y in zip(a, b))
def cross(a, b): return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])
def norm2(a): return dot(a, a)
def normal(p): return cross(sub(p[1], p[0]), sub(p[2], p[0]))
def edges(p): return [(p[i], p[(i+1) % 3]) for i in range(3)]


def validate(shape):
    kind, points = shape
    if not all(math.isfinite(x) for p in points for x in p): return 'N'
    _, p = exact(shape)
    if kind == 'L' and p[0] == p[1]: return 'D'
    if kind in 'FT' and norm2(normal(p)) == 0: return 'D'
    return 'R'


def on_line(p, line, segment):
    a, b = line
    d = sub(b, a)
    t = dot(sub(p, a), d) / norm2(d) if norm2(d) else Q(0)
    if segment: t = max(Q(0), min(Q(1), t))
    return add(a, mul(d, t))


def inside(p, tri):
    # Oriented edge halfspaces, independent of affine-coordinate elimination.
    n = normal(tri)
    return all(dot(cross(sub(b, a), sub(p, a)), n) >= 0 for a, b in edges(tri))


def point_triangle(p, tri):
    n = normal(tri)
    q = sub(p, mul(n, dot(sub(p, tri[0]), n)/norm2(n)))
    if inside(q, tri): return norm2(sub(p, q))
    return min(norm2(sub(p, on_line(p, edge, True))) for edge in edges(tri))


def line_line(a, b, bounded_a, bounded_b):
    u, v, w = sub(a[1], a[0]), sub(b[1], b[0]), sub(a[0], b[0])
    if not norm2(u): return norm2(sub(a[0], on_line(a[0], b, bounded_b)))
    if not norm2(v): return norm2(sub(b[0], on_line(b[0], a, bounded_a)))
    uv, uu, vv, uw, vw = dot(u, v), norm2(u), norm2(v), dot(u, w), dot(v, w)
    det = norm2(cross(u, v))
    candidates = []
    if det:
        t, s = (uv*vw-vv*uw)/det, (uu*vw-uv*uw)/det
        if (not bounded_a or 0 <= t <= 1) and (not bounded_b or 0 <= s <= 1):
            candidates.append(norm2(sub(add(a[0], mul(u, t)), add(b[0], mul(v, s)))))
    elif not bounded_a and not bounded_b:
        return norm2(cross(w, u))/uu
    if bounded_a:
        candidates += [norm2(sub(p, on_line(p, b, bounded_b))) for p in a]
    if bounded_b:
        candidates += [norm2(sub(p, on_line(p, a, bounded_a))) for p in b]
    return min(candidates)


def line_triangle(line, tri, bounded):
    n, d = normal(tri), sub(line[1], line[0])
    denominator = dot(n, d)
    if denominator:
        t = dot(n, sub(tri[0], line[0]))/denominator
        if (not bounded or 0 <= t <= 1) and inside(add(line[0], mul(d, t)), tri):
            return Q(0)
    candidates = [line_line(line, edge, bounded, True) for edge in edges(tri)]
    if bounded: candidates += [point_triangle(p, tri) for p in line]
    return min(candidates)


def distance2(a, b):
    """Analytic minimum for two already validated rational primitives."""
    ka, pa = a
    kb, pb = b
    if ka == 'F' or kb == 'F':
        if ka != 'F': return distance2(b, a)
        n = normal(pa)
        heights = [dot(n, sub(p, pa[0])) for p in pb]
        if kb == 'F':
            if norm2(cross(n, normal(pb))): return Q(0)
        elif kb == 'L':
            if heights[0] != heights[1]: return Q(0)
        elif min(heights) <= 0 <= max(heights):
            return Q(0)
        return min(h*h for h in heights)/norm2(n)
    if ka == 'P' or kb == 'P':
        if ka != 'P': return distance2(b, a)
        p = pa[0]
        if kb == 'P': return norm2(sub(p, pb[0]))
        if kb == 'T': return point_triangle(p, pb)
        return norm2(sub(p, on_line(p, pb, kb == 'S')))
    if ka == 'T' or kb == 'T':
        if kb != 'T': return distance2(b, a)
        if ka == 'T':
            return min([line_triangle(e, pb, True) for e in edges(pa)] +
                       [line_triangle(e, pa, True) for e in edges(pb)])
        return line_triangle(pa, pb, ka == 'S')
    return line_line(pa, pb, ka == 'S', kb == 'S')


def certify(shapes, points, parameters, squared_distance):
    """Prove global optimality from membership and supporting halfspaces.

    This does not solve a distance problem, select active faces, or invert a
    matrix. It checks the returned construction against each original vertex
    or unbounded direction. All comparisons use exact rational arithmetic.
    """
    delta = sub(points[0], points[1])
    assert squared_distance == norm2(delta)
    for operand, ((kind, vertices), p, parameters) in enumerate(zip(shapes, points, parameters)):
        assert len(parameters) == len(vertices)-1
        constructed = vertices[0]
        for t, v in zip(parameters, vertices[1:]):
            constructed = add(constructed, mul(sub(v, vertices[0]), t))
        assert constructed == p, (operand, 'membership construction')
        if kind in 'PST':
            assert all(t >= 0 for t in parameters) and sum(parameters) <= 1
            values = [dot(delta, sub(v, p)) for v in vertices]
            assert all(s >= 0 if operand == 0 else s <= 0 for s in values), (operand, 'support', values)
        else:
            assert all(dot(delta, sub(v, vertices[0])) == 0 for v in vertices[1:]), (operand, 'orthogonality')


def bits(x): return struct.pack('>d', x).hex()


def enclosure(value, sqrt=False):
    """Independent reference bracket, seeded by correctly rounded rational
    conversion (or high-precision integer square root), then tightened exactly.
    """
    maximum = Q.from_float(float.fromhex('0x1.fffffffffffffp+1023'))
    if value > (maximum*maximum if sqrt else maximum) or (not sqrt and value < -maximum): return ['U']
    if sqrt:
        assert value >= 0
        # sqrt(value) approximated on the 2^-1074 lattice, with a <1 ulp seed
        # for every binary64 magnitude, without overflowing/underflowing floats.
        root = math.isqrt((value.numerator << 2148)//value.denominator)
        guess = float(Q(root, 1 << 1074))
        compare = lambda x: (Q.from_float(x)**2 > value)-(Q.from_float(x)**2 < value)
    else:
        guess = float(value)
        compare = lambda x: (Q.from_float(x) > value)-(Q.from_float(x) < value)
    lower = upper = guess
    while compare(lower) > 0: lower = math.nextafter(lower, -math.inf)
    while compare(upper) < 0: upper = math.nextafter(upper, math.inf)
    return ['E', bits(lower), bits(upper)]
