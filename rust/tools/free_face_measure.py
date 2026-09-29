#!/usr/bin/env python3
"""Independent area and centre of a free spline face of a .brep file, in
high precision (mpmath), for reviewing compare_brep_io.py's measure
differences without Rust or OCCT.

The face's surface and pcurves are read as the binary64 values OCCT and the
kernel read. By Green's theorem on the parameter domain,
`∫∫_D f |S_u × S_v| du dv = ∮ Q dv` with `Q(u, v) = ∫_{u*}^{u} f(s, v)
|S_u × S_v|(s, v) ds` for `f = 1, x, y, z`, along every pcurve in its wire's
order and orientation, split at its knots, and along the straight UV chord
closing each gap between consecutive fins of a wire, as the kernel closes
its loops. Every integral is Gauss–Legendre at 40 digits; the result is
computed with `n` and `2n` nodes and both are printed, so the digits a
review records are those that agree. `--open` drops the chords and takes
`u*` from `--baseline`: the open sum is what BRepGProp's Green integration
computes from the face's UMin, which is how a UV gap in a wire moves OCCT's
measure.

Scope: nonperiodic B-spline surfaces (rational or not) and B-spline or line
pcurves, without locations. Usage:

    free_face_measure.py FILE.brep RECORD [--nodes 12] [--open --baseline U]
"""
import argparse
import sys

import mpmath as mp

import brep_io_reference as reference

mp.mp.dps = 40


def tables(text):
    """{'Curve2ds', 'Curves', 'Surfaces': [record words]} as written."""
    r = reference.Reader(text)
    while r.word() != 'Curve2ds':
        pass
    r.at -= 1
    out = {}
    for name, parse in (('Curve2ds', reference.curve2), ('Curves', reference.curve3)):
        assert r.word() == name
        out[name] = []
        for _ in range(r.int()):
            start = r.at
            parse(r)
            out[name].append(r.words[start:r.at])
    while r.word() != 'Surfaces':
        pass
    out['Surfaces'] = []
    for _ in range(r.int()):
        start = r.at
        reference.surface(r)
        out['Surfaces'].append(r.words[start:r.at])
    return out


def D(word):
    return mp.mpf(float(word))


def knots(words, count):
    flat = []
    for _ in range(count):
        x, m = D(next(words)), int(next(words))
        flat += [x]*m
    return flat


def pcurve(words):
    """A callable t -> ((u, v), (u', v')) and its breakpoints."""
    it = iter(words)
    kind = int(next(it))
    if kind == 1:
        p = [D(next(it)) for _ in range(2)]
        d = [D(next(it)) for _ in range(2)]
        return lambda t: ([p[0]+t*d[0], p[1]+t*d[1]], d), []
    if kind != 7:
        raise ValueError(f'pcurve type {kind} is out of scope')
    rational, periodic, degree, n, k = (int(next(it)) for _ in range(5))
    if periodic:
        raise ValueError('periodic pcurve is out of scope')
    poles, weights = [], []
    for _ in range(n):
        poles.append([D(next(it)) for _ in range(2)])
        weights.append(D(next(it)) if rational else mp.mpf(1))
    flat = knots(it, k)

    def point(t):
        span, N = basis(degree, flat, t)
        num, den = [mp.mpf(0)]*2, mp.mpf(0)
        for j in range(degree+1):
            i = span-degree+j
            den += N[j]*weights[i]
            for d in range(2):
                num[d] += N[j]*weights[i]*poles[i][d]
        return [x/den for x in num]

    h = mp.mpf('1e-20')
    return (lambda t: (point(t), [(b-a)/(2*h) for a, b in zip(point(t-h), point(t+h))]),
            sorted(set(flat)))


def basis(p, flat, t):
    """(span, the p + 1 nonzero basis values) at t; outside the knots the
    end spans' polynomials continue."""
    n = len(flat)-p-1
    span = p
    while span < n-1 and t >= flat[span+1]:
        span += 1
    N = [mp.mpf(1)]+[mp.mpf(0)]*p
    left, right = [0]*(p+1), [0]*(p+1)
    for j in range(1, p+1):
        left[j], right[j] = t-flat[span+1-j], flat[span+j]-t
        saved = mp.mpf(0)
        for r in range(j):
            tmp = N[r]/(right[r+1]+left[j-r])
            N[r] = saved+right[r+1]*tmp
            saved = left[j-r]*tmp
        N[j] = saved
    return span, N


class Surface:
    def __init__(self, words):
        it = iter(words)
        if int(next(it)) != 9:
            raise ValueError('only B-spline surfaces are in scope')
        ru, rv, pu, pv, self.du, self.dv, nu, nv, ku, kv = (int(next(it)) for _ in range(10))
        if pu or pv:
            raise ValueError('periodic surface is out of scope')
        self.poles = [[None]*nv for _ in range(nu)]
        self.weights = [[mp.mpf(1)]*nv for _ in range(nu)]
        for i in range(nu):
            for j in range(nv):
                self.poles[i][j] = [D(next(it)) for _ in range(3)]
                if ru or rv:
                    self.weights[i][j] = D(next(it))
        self.fu, self.fv = knots(it, ku), knots(it, kv)

    def point(self, u, v):
        su, Nu = basis(self.du, self.fu, u)
        sv, Nv = basis(self.dv, self.fv, v)
        num, den = [mp.mpf(0)]*3, mp.mpf(0)
        for a in range(self.du+1):
            for c in range(self.dv+1):
                i, j = su-self.du+a, sv-self.dv+c
                w = Nu[a]*Nv[c]*self.weights[i][j]
                den += w
                for d in range(3):
                    num[d] += w*self.poles[i][j][d]
        return [x/den for x in num]

    def density(self, u, v, h=mp.mpf('1e-20')):
        """|S_u × S_v| by central differences (error near h^2, 1e-40)."""
        su = [(b-a)/(2*h) for a, b in zip(self.point(u-h, v), self.point(u+h, v))]
        sv = [(b-a)/(2*h) for a, b in zip(self.point(u, v-h), self.point(u, v+h))]
        n = [su[1]*sv[2]-su[2]*sv[1], su[2]*sv[0]-su[0]*sv[2], su[0]*sv[1]-su[1]*sv[0]]
        return mp.sqrt(sum(x*x for x in n))


def fins(text, record):
    """(surface index, {wire: [(pcurve index, (first, last), orientation)]})."""
    _, _, shapes, _ = reference.read(text)
    kind, data, wires = shapes[record]
    if kind != 'Fa' or data is None or data[1] != 0:
        raise ValueError('not a face without location')
    out = {}
    for k, (_, w, wl) in enumerate(wires):
        if wl != 0:
            raise ValueError('wire location is out of scope')
        for eo, e, el in shapes[w][2]:
            reps = [x for x in shapes[e][1][1] if x[0] == 'pcurve' and x[2] == data[0] and x[3] == 0]
            if el != 0 or len(reps) != 1 or len(reps[0][1]) != 1:
                raise ValueError('edge location or seam is out of scope')
            out.setdefault(k, []).append((reps[0][1][0], reps[0][4], eo))
    return data[0], out


def measure(text, record, nodes, chords=True, baseline=None):
    """The signed [area, x, y, z moments] of the face."""
    table = tables(text)
    surf, wires = fins(text, record)
    S = Surface(table['Surfaces'][surf-1])
    ustar = S.fu[0] if baseline is None else mp.mpf(baseline)
    xs, ws = zip(*legendre(nodes))

    def Q(u, v):
        acc = [mp.mpf(0)]*4
        half, mid = (u-ustar)/2, (u+ustar)/2
        for x, w in zip(xs, ws):
            s = mid+half*x
            p, j = S.point(s, v), S.density(s, v)
            acc[0] += w*j
            for d in range(3):
                acc[d+1] += w*j*p[d]
        return [a*half for a in acc]

    total = [mp.mpf(0)]*4

    def along(at, a, b, sign):
        half, mid = (b-a)/2, (b+a)/2
        for x, w in zip(xs, ws):
            (u, v), (_, dv) = at(mid+half*x)
            q = Q(u, v)
            for d in range(4):
                total[d] += sign*w*half*q[d]*dv

    for wire in wires.values():
        ends = []
        for index, (first, last), orientation in wire:
            at, breaks = pcurve(table['Curve2ds'][index-1])
            cuts = [mp.mpf(first)]+[k for k in breaks if first < k < last]+[mp.mpf(last)]
            sign = 1 if orientation == '+' else -1
            for a, b in zip(cuts, cuts[1:]):
                along(at, a, b, sign)
            e = [at(mp.mpf(first))[0], at(mp.mpf(last))[0]]
            ends.append(e if orientation == '+' else e[::-1])
        for (_, p), (q, _) in zip(ends, ends[1:]+ends[:1]) if chords else []:
            if p != q:
                d = [b-a for a, b in zip(p, q)]
                along(lambda t: ([a+t*da for a, da in zip(p, d)], d), mp.mpf(0), mp.mpf(1), 1)
    return total


def legendre(n):
    """n-point Gauss–Legendre nodes and weights on [-1, 1] at the working
    precision, by Newton's method on P_n."""
    out = []
    for i in range(1, n+1):
        x = mp.cos(mp.pi*(i-mp.mpf(1)/4)/(n+mp.mpf(1)/2))
        for _ in range(100):
            p0, p1 = mp.mpf(1), x
            for k in range(2, n+1):
                p0, p1 = p1, ((2*k-1)*x*p1-(k-1)*p0)/k
            dp = n*(x*p1-p0)/(x*x-1)
            step = p1/dp
            x -= step
            if abs(step) < mp.mpf(10)**(-mp.mp.dps+2):
                break
        out.append((x, 2/((1-x*x)*dp*dp)))
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('file')
    parser.add_argument('record', type=int)
    parser.add_argument('--nodes', type=int, default=12)
    parser.add_argument('--open', action='store_true', help='no chords (BRepGProp\'s open sum)')
    parser.add_argument('--baseline', help='u* of the open sum (BRepGProp: the face\'s UMin)')
    args = parser.parse_args()
    text = open(args.file, errors='replace').read()
    for n in (args.nodes, 2*args.nodes):
        t = measure(text, args.record, n, not args.open, args.baseline)
        print(n, 'nodes: area', mp.nstr(abs(t[0]), 20), 'centre',
              ' '.join(mp.nstr(x/t[0], 20) for x in t[1:]))


if __name__ == '__main__':
    sys.exit(main())
