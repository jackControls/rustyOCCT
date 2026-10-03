#!/usr/bin/env python3
"""Independent triple points for S9e.3b of REVIEW_NOTES.md: where a face of
the last solid of a chain meets the given result's edges on meetings of two
of its solids' surfaces (`chained_curved_boolean_reference.py`'s chains).

S9e.3a's reference decides a chain by its solids' memberships, so a partner
across a given meeting of two curved faces needs nothing new there; this
module checks that the fixtures reach S9e.3b's meetings at all, and that
they do so clear of tangencies, by finding those meetings on their own:

* **Tracing.** For two solids `i` and `j` of the given result, each face of
  `i` is swept by its families of curves (S9d.4b.2's, as the chained
  reference sweeps it), and along each scanned curve the roots of each
  surface of `j` (`curve_roots`) are points of the two surfaces' meeting;
  matched from one scanned curve to the next (nearest points, where their
  number does not change), they trace the meeting across the face. Each
  surface of the last solid is evaluated at the traced points; where its
  sign changes between two matched points, the triple point between them is
  refined by Newton's method on the curve's two parameters (the second
  surface's and the third's functions both zero, the family's derivatives
  in closed form) to 40 digits. Both families of every face of both solids
  trace the same meeting, so each triple point is found several times
  (`found`), from either solid's faces.
* **On the given result's edge and the partner's face.** A triple point is
  kept where the given result's set (its expression of the solids'
  memberships) holds one or three of the four quadrants about the meeting
  there (the point pushed by `PUSH` of the case's size along both unit
  normals, either way: an edge of the given result along the meeting, not a
  face continuing across it), and the last solid's membership differs across
  its surface there (its face holds the point).
* **Margins.** Each kept point's sine: the determinant of the three unit
  normals (zero where the third surface is tangent to the meeting, the
  three normals dependent). A declared tangency is checked at its declared
  point: the three functions zero and the determinant zero within 1e-30,
  the point on the given result's edge and on the partner's face.

Nothing here uses the kernel, a resultant or an eliminant: the points are
roots along the reference's own curves, refined in floating point of 40
digits.
"""
import mpmath as mp

from torus_curved_boolean_reference import Z, curve_roots, norm
from sphere_boolean_reference import cross, dot
import chained_curved_boolean_reference as ref

mp.mp.dps = 40

SCAN = 360
# The push off the meeting to read the quadrants, relative to the size.
PUSH = mp.mpf(10)**-18
# Two triple points within this of the size are one.
SAME = mp.mpf(10)**-20


def face_surf(I, face):
    """The surface a face of an input lies on."""
    if I.kind == 'prism' or face.name != 'wall':
        return next(s for s in I.surfs if s.label == face.name)
    return I.surfs[0]


def surf_kind(I, s):
    """'plane', 'cylinder', 'sphere', 'cone' or 'torus'."""
    if s.degree == 1:
        return 'plane'
    if I.kind == 'prism':
        return 'cylinder'
    return I.kind


def in_scope(k1, k2, k3):
    """Whether a meeting of surfaces of kinds `k1` and `k2` met by `k3` is
    S9e.3b's: two curved surfaces (not two spheres, whose meeting is a
    circle) met by any; a plane and a cone or a torus (a section) met by a
    curved one."""
    curved = ('cylinder', 'sphere', 'cone', 'torus')
    if k1 in curved and k2 in curved:
        return not (k1 == k2 == 'sphere')
    pair = {k1, k2}
    return 'plane' in pair and bool(pair & {'cone', 'torus'}) and k3 in curved


def unit(v):
    n = norm(v)
    return tuple(x/n for x in v)


def det3(a, b, c):
    return dot(a, cross(b, c))


class Tracer:
    """The chain's inputs, the given result's expression and its solids,
    and the partner (the last solid)."""

    def __init__(self, chain, given, partner):
        self.chain = chain
        self.inputs = chain.inputs
        self.given = given
        self.solids = sorted(ref.leaves(given))
        self.partner = partner
        self.size = chain.size

    def members(self, X):
        return [I.contains(X) for I in self.inputs]

    def on_edge(self, X, n1, n2):
        """Whether the given result has an edge along the meeting at `X`."""
        e = PUSH*self.size
        u1, u2 = unit(n1), unit(n2)
        count = 0
        for s1 in (1, -1):
            for s2 in (1, -1):
                P = tuple(X[i]+e*(s1*u1[i]+s2*u2[i]) for i in range(3))
                count += bool(ref.evaluate(self.given, self.members(P)))
        return count in (1, 3)

    def on_partner(self, X, n3):
        e = PUSH*self.size
        u = unit(n3)
        P = self.inputs[self.partner]
        f = P.contains(tuple(X[i]+e*u[i] for i in range(3)))
        b = P.contains(tuple(X[i]-e*u[i] for i in range(3)))
        return f != b

    def newton(self, fam, s2, s3, a, b):
        """The triple point near `(a, b)` on the family's curves: `s2` and
        `s3` both zero; None where it does not converge."""
        for _ in range(60):
            cur = fam.curve(a)
            X, Xa, Xb = cur.frame_at(b)
            f2, f3 = s2.f(X), s3.f(X)
            g2, g3 = s2.grad(X), s3.grad(X)
            j = ((dot(g2, Xa), dot(g2, Xb)), (dot(g3, Xa), dot(g3, Xb)))
            d = j[0][0]*j[1][1]-j[0][1]*j[1][0]
            if d == 0:
                return None
            da = (f2*j[1][1]-f3*j[0][1])/d
            db = (j[0][0]*f3-j[1][0]*f2)/d
            a, b = a-da, b-db
            if abs(da)+abs(db) < mp.mpf(10)**-36*(1+abs(a)+abs(b)):
                X = fam.curve(a).point(b)
                if abs(s2.f(X)) < mp.mpf(10)**-30*self.size**s2.degree \
                        and abs(s3.f(X)) < mp.mpf(10)**-30*self.size**s3.degree:
                    return X
        return None

    def traced(self, fam, s2, s3s, scan=SCAN):
        """The triple points found on a family's curves (the meeting of the
        face's surface with `s2`, crossed by each of `s3s`): [(X, s3)]."""
        span = fam.a1-fam.a0
        avals = [fam.a0+span*(k+mp.mpf(1)/2)/scan for k in range(scan)]
        rows = []
        for a in avals:
            cur = fam.curve(a)
            if cur.degenerate():
                rows.append((a, cur, []))
                continue
            bs = curve_roots(cur, s2)
            rows.append((a, cur, [(b, cur.point(b), [s3.f(cur.point(b)) for s3 in s3s]) for b in bs]))
        pairs = list(zip(rows, rows[1:]))
        if fam.periodic:
            last, first = rows[-1], rows[0]
            pairs.append((last, (first[0]+span, first[1], first[2])))
        out = []
        for (a0, cur, r0), (a1, _, r1) in pairs:
            if not r0 or len(r0) != len(r1):
                continue
            used = set()
            for b0, X0, v0 in r0:
                best, k = None, None
                for kk, (_, X1, _) in enumerate(r1):
                    if kk in used:
                        continue
                    d = norm(tuple(X1[i]-X0[i] for i in range(3)))
                    if best is None or d < best:
                        best, k = d, kk
                used.add(k)
                b1, X1, v1 = r1[k]
                if cur.kind == 'circle':
                    # The nearer image of b1 about b0.
                    b1 = b0+((b1-b0+mp.pi) % (2*mp.pi))-mp.pi
                for s3, x0, x1 in zip(s3s, v0, v1):
                    if x0*x1 >= 0:
                        continue
                    X = self.newton(fam, s2, s3, (a0+a1)/2, (b0+b1)/2)
                    if X is not None:
                        out.append((X, s3))
        return out

    def points(self):
        """[(X, sine, kinds, found)] of the triple points on the given
        result's edges inside the partner's faces."""
        found = []
        for i in self.solids:
            for j in self.solids:
                if j == i:
                    continue
                I, J = self.inputs[i], self.inputs[j]
                s3s = self.inputs[self.partner].surfs
                for face in I.faces:
                    s1 = face_surf(I, face)
                    for s2 in J.surfs:
                        for fam in face.families:
                            for X, s3 in self.traced(fam, s2, s3s):
                                found.append((X, (I, s1), (J, s2), s3))
        merged = []
        for X, (I, s1), (J, s2), s3 in found:
            for m in merged:
                if norm(tuple(X[i]-m[0][i] for i in range(3))) < SAME*self.size:
                    m[4] += 1
                    break
            else:
                merged.append([X, (I, s1), (J, s2), s3, 1])
        out = []
        P = self.inputs[self.partner]
        for X, (I, s1), (J, s2), s3, n in merged:
            n1, n2, n3 = s1.grad(X), s2.grad(X), s3.grad(X)
            if not self.on_edge(X, n1, n2) or not self.on_partner(X, n3):
                continue
            sine = abs(det3(unit(n1), unit(n2), unit(n3)))
            kinds = (surf_kind(I, s1), surf_kind(J, s2), surf_kind(P, s3))
            out.append((X, sine, kinds, n))
        return out

    def declared(self, X):
        """The deviations at a declared tangency `X`: the largest of the
        three surfaces' values (relative to the size to their degrees) and
        the unit normals' determinant, with whether it lies on the given
        result's edge and the partner's face (the surfaces vanishing there:
        one of each given solid's, one of the partner's)."""
        X = tuple(mp.mpf(x) for x in X)
        zero = mp.mpf(10)**-30

        def vanishing(I):
            return [s for s in I.surfs if abs(s.f(X)) <= zero*self.size**s.degree]
        given = [(k, s) for k in self.solids for s in vanishing(self.inputs[k])]
        third = vanishing(self.inputs[self.partner])
        assert len(given) == 2 and len(third) == 1, 'a declared tangency on two given surfaces and the partner\'s'
        (_, s1), (_, s2) = given
        s3 = third[0]
        n1, n2, n3 = s1.grad(X), s2.grad(X), s3.grad(X)
        values = max(abs(s.f(X))/self.size**s.degree for s in (s1, s2, s3))
        d = abs(det3(unit(n1), unit(n2), unit(n3)))
        return values, d, self.on_edge(X, n1, n2), self.on_partner(X, n3)


__all__ = ['Tracer', 'in_scope', 'SCAN']
