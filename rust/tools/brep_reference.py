"""Generic B-rep cases, the shared text protocol, and an independent validator.

The validator is deliberately different from production: geometric decisions
use high-precision mpmath evaluation with explicit Lipschitz margins, Green's
theorem quadrature for areas and volumes, and ray casting for containment.
Production instead uses certified interval branch-and-bound. Cases are built
far from every tolerance boundary, so a margin failure raises instead of
guessing. No Rust or OCCT result supplies an expected issue.

Conventions (shared with Rust): every 3D curve and pcurve is parameterized by a
normalized fraction t in [0,1]. A use's pcurve follows the oriented face's
traversal, so a reversed use pairs pcurve t with edge fraction 1-t. Loops follow
the oriented face: outer loops are counter-clockwise about the oriented normal.
"""
from dataclasses import dataclass, field, replace
from fractions import Fraction as F
import copy
import math

import mpmath as mp

mp.mp.dps = 40


# Correctly rounded binary64 trigonometry and norms. Platform libm results can
# differ in the last bit (macOS and glibc do), which would make the generated
# fixtures and native inputs depend on the host.
def cos_rn(x):
    with mp.workprec(256):
        return float(mp.cos(mp.mpf(x)))


def sin_rn(x):
    with mp.workprec(256):
        return float(mp.sin(mp.mpf(x)))


def hypot_rn(*xs):
    """Correctly rounded Euclidean norm; math.hypot changed in CPython 3.10."""
    with mp.workprec(256):
        return float(mp.sqrt(mp.fsum(mp.mpf(x)**2 for x in xs)))


def atan2_rn(y, x):
    with mp.workprec(256):
        return float(mp.atan2(mp.mpf(y), mp.mpf(x)))
TAU = 2*math.pi


# ---------------------------------------------------------------- model

@dataclass
class Frame:
    origin: tuple
    normal: tuple
    x: tuple  # hint; the effective frame is orthonormalized like Frame3::new


@dataclass
class Line3:
    start: tuple
    end: tuple


@dataclass
class Arc3:
    frame: Frame
    radius: float
    start: float
    sweep: float


@dataclass
class Edge:
    start: int
    end: int
    curve: object


@dataclass
class Line2:
    start: tuple
    end: tuple


@dataclass
class Arc2:
    center: tuple
    radius: float
    start: float
    sweep: float


@dataclass
class Use:
    edge: int
    forward: bool
    pcurve: object


@dataclass
class Plane:
    frame: Frame


@dataclass
class Cylinder:
    frame: Frame
    radius: float


@dataclass
class Face:
    surface: object
    forward: bool
    loops: list = field(default_factory=list)


@dataclass
class Model:
    name: str
    tolerance: float
    vertices: list = field(default_factory=list)
    edges: list = field(default_factory=list)
    faces: list = field(default_factory=list)
    shells: list = field(default_factory=list)


# ---------------------------------------------------------------- protocol

def number(x):
    return repr(float(x))


def encode(m):
    """Line protocol read by the Rust fixture test."""
    out = [f'case {m.name}', f'tolerance {number(m.tolerance)}']
    frame = lambda f: ' '.join(map(number, (*f.origin, *f.normal, *f.x)))
    for v in m.vertices:
        out.append('v '+' '.join(map(number, v)))
    for e in m.edges:
        c = e.curve
        if isinstance(c, Line3):
            out.append(f'e {e.start} {e.end} line '+' '.join(map(number, (*c.start, *c.end))))
        else:
            out.append(f'e {e.start} {e.end} arc {frame(c.frame)} '+' '.join(map(number, (c.radius, c.start, c.sweep))))
    for f in m.faces:
        s = f.surface
        o = 'F' if f.forward else 'R'
        if isinstance(s, Plane):
            out.append(f'f plane {frame(s.frame)} {o}')
        else:
            out.append(f'f cylinder {frame(s.frame)} {number(s.radius)} {o}')
        for loop in f.loops:
            out.append('l')
            for u in loop:
                p = u.pcurve
                o = 'F' if u.forward else 'R'
                if isinstance(p, Line2):
                    out.append(f'u {u.edge} {o} line '+' '.join(map(number, (*p.start, *p.end))))
                else:
                    out.append(f'u {u.edge} {o} arc '+' '.join(map(number, (*p.center, p.radius, p.start, p.sweep))))
    for s in m.shells:
        out.append('s '+' '.join(map(str, s)))
    out.append('end')
    return '\n'.join(out)


# ---------------------------------------------------------------- geometry (exact frames, mp evaluation)

def vec(a):
    return [mp.mpf(x) for x in a]


def sub(a, b):
    return [x-y for x, y in zip(a, b)]


def add(a, b):
    return [x+y for x, y in zip(a, b)]


def mul(a, s):
    return [x*s for x in a]


def ltr_sum(xs):
    """Left-to-right binary64 sum. Python 3.12's sum() compensates float
    rounding (Neumaier), so it differs from 3.9 in the last bits; fixture and
    native-input bytes must not depend on the Python version."""
    total = 0.0
    for x in xs:
        total += x
    return total


def dot(a, b):
    return ltr_sum(x*y for x, y in zip(a, b))


def cross(a, b):
    return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]


def norm(a):
    return mp.sqrt(dot(a, a))


def axes(f):
    """Orthonormal (origin, x, y, n) exactly as Frame3::new, in high precision."""
    n = vec(f.normal)
    h = vec(f.x)
    if norm(n) == 0 or norm(h) == 0:
        return None
    n = mul(n, 1/norm(n))
    h = mul(h, 1/norm(h))
    y = cross(n, h)
    if norm(y) <= mp.mpf('1e-12'):
        return None
    y = mul(y, 1/norm(y))
    x = cross(y, n)
    return vec(f.origin), x, y, n


def finite(values):
    return all(math.isfinite(v) for v in values)


def curve_valid(c, tol):
    if isinstance(c, Line3):
        return finite((*c.start, *c.end)) and norm(sub(vec(c.end), vec(c.start))) > tol
    return (finite((*c.frame.origin, *c.frame.normal, *c.frame.x, c.radius, c.start, c.sweep))
            and axes(c.frame) is not None and c.radius > tol and 0 < abs(c.sweep) <= TAU)


def surface_valid(s, tol):
    frame_ok = finite((*s.frame.origin, *s.frame.normal, *s.frame.x)) and axes(s.frame) is not None
    return frame_ok and (isinstance(s, Plane) or (math.isfinite(s.radius) and s.radius > tol))


def pcurve_valid(p):
    if isinstance(p, Line2):
        return finite((*p.start, *p.end)) and p.start != p.end
    return finite((*p.center, p.radius, p.start, p.sweep)) and p.radius > 0 and 0 < abs(p.sweep) <= TAU


def curve_point(c, t):
    t = mp.mpf(t)
    if isinstance(c, Line3):
        a, b = vec(c.start), vec(c.end)
        return add(a, mul(sub(b, a), t))
    o, x, y, _ = axes(c.frame)
    a = mp.mpf(c.start)+mp.mpf(c.sweep)*t
    return add(o, add(mul(x, c.radius*mp.cos(a)), mul(y, c.radius*mp.sin(a))))


def pcurve_point(p, t):
    t = mp.mpf(t)
    if isinstance(p, Line2):
        a, b = vec(p.start), vec(p.end)
        return add(a, mul(sub(b, a), t))
    a = mp.mpf(p.start)+mp.mpf(p.sweep)*t
    c = vec(p.center)
    return [c[0]+p.radius*mp.cos(a), c[1]+p.radius*mp.sin(a)]


def surface_point(s, uv):
    o, x, y, n = axes(s.frame)
    u, v = uv
    if isinstance(s, Plane):
        return add(o, add(mul(x, u), mul(y, v)))
    radial = add(mul(x, s.radius*mp.cos(u)), mul(y, s.radius*mp.sin(u)))
    return add(add(o, radial), mul(n, v))


class Harmonic:
    """A(t) = a0 + a1 t + sum over frequencies w of c_w cos(w t) + s_w sin(w t).

    Frequencies are exact binary64 magnitudes. Terms of equal frequency merge.
    Two terms at nearby frequencies w1 < w2 are Re(z1 e^{i w1 t}) and
    Re(z2 e^{i w2 t}) with z = c - i s; on [0, 1] their sum is at most the
    ellipse bound of z1 + z2 plus |z2| (w2 - w1), because |e^{i d t} - 1| <=
    d t. Adjacent terms by frequency are bounded that way when it is smaller.
    """
    def __init__(self, dim):
        self.a0 = [mp.mpf(0)]*dim
        self.a1 = [mp.mpf(0)]*dim
        self.terms = {}

    def affine(self, a0, a1, sign=1):
        self.a0 = add(self.a0, mul(a0, sign))
        self.a1 = add(self.a1, mul(a1, sign))

    def rotating(self, alpha, omega, cos_vector, sin_vector, sign=1):
        """sign*(cos(alpha+omega t) cos_vector + sin(alpha+omega t) sin_vector)."""
        if omega == 0:
            self.affine(add(mul(cos_vector, mp.cos(alpha)), mul(sin_vector, mp.sin(alpha))), [0]*len(cos_vector), sign)
            return
        w = abs(omega)
        # cos(alpha + omega t) = cos(alpha) cos(wt) - sgn sin(alpha) sin(wt)
        g = 1 if omega > 0 else -1
        ca, sa = mp.cos(alpha), mp.sin(alpha)
        cterm = add(mul(cos_vector, ca), mul(sin_vector, sa))
        sterm = add(mul(cos_vector, -g*sa), mul(sin_vector, g*ca))
        c, s = self.terms.get(w, ([mp.mpf(0)]*len(cos_vector),)*2)
        self.terms[w] = (add(c, mul(cterm, sign)), add(s, mul(sterm, sign)))

    def upper(self):
        bound = max(norm(self.a0), norm(add(self.a0, self.a1)))

        def ellipse(c, s):
            cc, ss, cs = dot(c, c), dot(s, s), dot(c, s)
            return mp.sqrt((cc+ss+mp.sqrt((cc-ss)**2+4*cs*cs))/2)
        terms = sorted(self.terms.items())
        k = 0
        while k < len(terms):
            w1, (c1, s1) = terms[k]
            single = ellipse(c1, s1)
            if k+1 < len(terms):
                w2, (c2, s2) = terms[k+1]
                apart = single+ellipse(c2, s2)
                paired = ellipse(add(c1, c2), add(s1, s2))+mp.sqrt(dot(c2, c2)+dot(s2, s2))*(w2-w1)
                if paired < apart:
                    bound += paired
                    k += 2
                    continue
            bound += single
            k += 1
        return bound


def curve_harmonic(c, forward, h, sign):
    """Add sign*C(tau(t)) with tau(t) = t or 1-t."""
    if isinstance(c, Line3):
        a, b = vec(c.start), vec(c.end)
        if forward:
            h.affine(a, sub(b, a), sign)
        else:
            h.affine(b, sub(a, b), sign)
        return
    o, x, y, _ = axes(c.frame)
    alpha, omega = mp.mpf(c.start), mp.mpf(c.sweep)
    if not forward:
        alpha, omega = alpha+omega, -omega
    h.affine(o, [0, 0, 0], sign)
    h.rotating(alpha, omega, mul(x, c.radius), mul(y, c.radius), sign)


def use_harmonic(s, p, h, sign):
    """Add sign*S(P(t)); returns False when not a harmonic expression."""
    o, x, y, n = axes(s.frame)
    if isinstance(s, Plane):
        if isinstance(p, Line2):
            a, b = vec(p.start), vec(p.end)
            h.affine(add(o, add(mul(x, a[0]), mul(y, a[1]))), add(mul(x, b[0]-a[0]), mul(y, b[1]-a[1])), sign)
        else:
            c = vec(p.center)
            h.affine(add(o, add(mul(x, c[0]), mul(y, c[1]))), [0, 0, 0], sign)
            h.rotating(mp.mpf(p.start), mp.mpf(p.sweep), mul(x, p.radius), mul(y, p.radius), sign)
        return True
    if isinstance(p, Line2):
        u0, v0 = mp.mpf(p.start[0]), mp.mpf(p.start[1])
        du, dv = mp.mpf(p.end[0])-u0, mp.mpf(p.end[1])-v0
        h.affine(add(o, mul(n, v0)), mul(n, dv), sign)
        h.rotating(u0, du, mul(x, s.radius), mul(y, s.radius), sign)
        return True
    return False


def deviation_bounds(c, s, p, forward, samples=256):
    """(certain lower bound by sampling, rigorous harmonic upper bound or inf)."""
    def gap(t):
        return norm(sub(curve_point(c, t if forward else 1-t), surface_point(s, pcurve_point(p, t))))
    low = max(gap(mp.mpf(i)/samples) for i in range(samples+1))
    h = Harmonic(3)
    curve_harmonic(c, forward, h, 1)
    high = h.upper() if use_harmonic(s, p, h, -1) else mp.inf
    return low, high


# ---------------------------------------------------------------- independent validator

def issue(kind, entity):
    return (kind, entity)


def use_vertices(m, u):
    e = m.edges[u.edge]
    return (e.start, e.end) if u.forward else (e.end, e.start)


def validate(m):
    """Complete sorted issue list for the milestone-1 contract."""
    issues = []
    nv, ne, nf = len(m.vertices), len(m.edges), len(m.faces)
    for i, e in enumerate(m.edges):
        if not (0 <= e.start < nv and 0 <= e.end < nv):
            issues.append(issue('reference', f'edge {i}'))
    for fi, f in enumerate(m.faces):
        for li, loop in enumerate(f.loops):
            for ui, u in enumerate(loop):
                if not 0 <= u.edge < ne:
                    issues.append(issue('reference', f'use {fi}.{li}.{ui}'))
    for si, s in enumerate(m.shells):
        if any(not 0 <= fi < nf for fi in s):
            issues.append(issue('reference', f'shell {si}'))
    if issues:
        return sorted(issues)
    tol = mp.mpf(m.tolerance)

    used_vertices = {v for e in m.edges for v in (e.start, e.end)}
    for v in range(nv):
        if v not in used_vertices:
            issues.append(issue('unused_vertex', f'vertex {v}'))
    uses = {}
    for fi, f in enumerate(m.faces):
        for li, loop in enumerate(f.loops):
            for ui, u in enumerate(loop):
                uses.setdefault(u.edge, []).append((fi, li, ui, u.forward))
    for e in range(ne):
        if e not in uses:
            issues.append(issue('unused_edge', f'edge {e}'))
    owner = {}
    for si, s in enumerate(m.shells):
        if not s:
            issues.append(issue('empty_shell', f'shell {si}'))
        for fi in s:
            owner.setdefault(fi, []).append(si)
    for fi, f in enumerate(m.faces):
        if fi not in owner:
            issues.append(issue('face_without_shell', f'face {fi}'))
        elif len(owner[fi]) > 1:
            issues.append(issue('face_reused', f'face {fi}'))
        if not f.loops:
            issues.append(issue('empty_face', f'face {fi}'))
        for li, loop in enumerate(f.loops):
            if not loop:
                issues.append(issue('empty_loop', f'loop {fi}.{li}'))
                continue
            ends = [use_vertices(m, u) for u in loop]
            if any(ends[k][1] != ends[(k+1) % len(ends)][0] for k in range(len(ends))):
                issues.append(issue('open_loop', f'loop {fi}.{li}'))

    # Edge accounting uses only faces owned by exactly one shell.
    shell_of = {fi: s[0] for fi, s in owner.items() if len(s) == 1}
    bad_shells = set()
    for e, lst in uses.items():
        lst = [x for x in lst if x[0] in shell_of]
        shells = {shell_of[x[0]] for x in lst}
        if len(shells) > 1:
            issues.append(issue('edge_across_shells', f'edge {e}'))
            bad_shells |= shells
        elif len(lst) == 1:
            issues.append(issue('free_edge', f'edge {e}'))
            bad_shells |= shells
        elif len(lst) > 2:
            issues.append(issue('non_manifold_edge', f'edge {e}'))
            bad_shells |= shells
        elif len(lst) == 2 and lst[0][3] == lst[1][3]:
            issues.append(issue('same_sense_uses', f'edge {e}'))
            bad_shells |= shells
    for (kind, ent) in issues:
        if kind in ('open_loop', 'empty_loop', 'empty_face'):
            fi = int(ent.split()[1].split('.')[0])
            bad_shells |= set(owner.get(fi, []))
    for si, s in enumerate(m.shells):
        if si in bad_shells or not s or any(len(owner.get(fi, [])) != 1 for fi in s):
            bad_shells.add(si)
            continue
        faces = set(s)
        # Face connectivity through shared edges.
        adjacency = {fi: set() for fi in faces}
        for e, lst in uses.items():
            fs = {x[0] for x in lst if x[0] in faces}
            for a in fs:
                adjacency[a] |= fs-{a}
        seen, stack = set(), [s[0]]
        while stack:
            f = stack.pop()
            if f not in seen:
                seen.add(f)
                stack.extend(adjacency[f])
        if seen != faces:
            issues.append(issue('disconnected_shell', f'shell {si}'))
            bad_shells.add(si)
            continue
        # Vertex links: corners connect the incoming and outgoing edges at v.
        links = {}
        for fi in s:
            for loop in m.faces[fi].loops:
                for k, u in enumerate(loop):
                    w = loop[(k+1) % len(loop)]
                    v = use_vertices(m, u)[1]
                    g = links.setdefault(v, {})
                    g.setdefault(u.edge, set()).add(w.edge)
                    g.setdefault(w.edge, set()).add(u.edge)
        pinched = False
        for v in sorted(links):
            g = links[v]
            start = next(iter(g))
            seen, stack = set(), [start]
            while stack:
                x = stack.pop()
                if x not in seen:
                    seen.add(x)
                    stack.extend(g[x])
            if seen != set(g):
                issues.append(issue('non_manifold_vertex', f'vertex {v}'))
                pinched = True
        if pinched:
            bad_shells.add(si)
            continue
        vs = {x for fi in s for loop in m.faces[fi].loops for u in loop for x in use_vertices(m, u)}
        es = {u.edge for fi in s for loop in m.faces[fi].loops for u in loop}
        loops = sum(len(m.faces[fi].loops) for fi in s)
        chi = len(vs)-len(es)+2*len(s)-loops
        if chi % 2 or chi > 2:
            issues.append(issue('euler', f'shell {si}'))
            bad_shells.add(si)

    # ------------------------------------------------ certified-in-production geometry
    vertex_ok = []
    for v, p in enumerate(m.vertices):
        vertex_ok.append(finite(p))
        if not vertex_ok[-1]:
            issues.append(issue('degenerate_vertex', f'vertex {v}'))
    curve_ok = []
    for i, e in enumerate(m.edges):
        ok = curve_valid(e.curve, tol)
        curve_ok.append(ok)
        if not ok:
            issues.append(issue('degenerate_curve', f'edge {i}'))
    surface_ok = []
    for fi, f in enumerate(m.faces):
        ok = surface_valid(f.surface, tol)
        surface_ok.append(ok)
        if not ok:
            issues.append(issue('degenerate_surface', f'face {fi}'))
    geometry_bad_faces = {fi for fi, ok in enumerate(surface_ok) if not ok}
    for i, e in enumerate(m.edges):
        if not curve_ok[i]:
            continue
        for end, (v, t) in enumerate(((e.start, 0), (e.end, 1))):
            if not vertex_ok[v]:
                continue
            d = norm(sub(curve_point(e.curve, t), vec(m.vertices[v])))
            margin_check(d, tol, f'vertex on edge {i}')
            if d > tol:
                issues.append(issue('vertex_off_curve', f'edge {i} {"start" if end == 0 else "end"}'))
    for fi, f in enumerate(m.faces):
        for li, loop in enumerate(f.loops):
            for ui, u in enumerate(loop):
                ent = f'use {fi}.{li}.{ui}'
                if not pcurve_valid(u.pcurve):
                    issues.append(issue('degenerate_pcurve', ent))
                    geometry_bad_faces.add(fi)
                    continue
                if not surface_ok[fi] or not curve_ok[u.edge]:
                    continue
                low, high = deviation_bounds(m.edges[u.edge].curve, f.surface, u.pcurve, u.forward)
                assert low <= high*(1+mp.mpf('1e-30'))+mp.mpf('1e-30'), (m.name, ent, low, high)
                if low > tol:
                    issues.append(issue('pcurve_off_edge', ent))
                    geometry_bad_faces.add(fi)
                elif high > tol:
                    raise ArithmeticError(f'{m.name}: {ent} is too close to tolerance for the oracle')

    # UV continuity: each use must end where the next begins, in surface
    # length (radius*du on a cylinder).
    for fi, f in enumerate(m.faces):
        if not surface_ok[fi]:
            continue
        for li, loop in enumerate(f.loops):
            if any(not pcurve_valid(u.pcurve) for u in loop):
                continue
            for ui, u in enumerate(loop):
                w = loop[(ui+1) % len(loop)]
                a, b = pcurve_point(u.pcurve, 1), pcurve_point(w.pcurve, 0)
                du, dv = a[0]-b[0], a[1]-b[1]
                if isinstance(f.surface, Cylinder):
                    du *= f.surface.radius
                d = mp.sqrt(du*du+dv*dv)
                margin_check(d, tol, f'{m.name}: uv gap')
                if d > tol:
                    issues.append(issue('uv_gap', f'use {fi}.{li}.{ui}'))
                    geometry_bad_faces.add(fi)

    # Loop winding and imbrication, for structurally and geometrically sound faces.
    structural_faces = {int(ent.split()[1].split('.')[0]) for k, ent in issues if k in ('open_loop', 'empty_loop')}
    for fi, f in enumerate(m.faces):
        if fi in geometry_bad_faces or fi in structural_faces or not f.loops:
            continue
        areas = [loop_area(loop) for loop in f.loops]
        for li, a in enumerate(areas):
            want = (1 if li == 0 else -1)*(1 if f.forward else -1)
            if abs(a) < mp.mpf('1e-25'):
                raise ArithmeticError(f'{m.name}: loop area too small for the oracle')
            if a*want <= 0:
                issues.append(issue('loop_winding', f'loop {fi}.{li}'))
        for li in range(1, len(f.loops)):
            point = pcurve_point(f.loops[li][0].pcurve, 0)
            if winding(f.loops[0], point) == 0:
                issues.append(issue('inner_loop_outside', f'loop {fi}.{li}'))

    # Shell orientation and cavity nesting, for fully sound shells.
    bad_faces = geometry_bad_faces | structural_faces | {int(ent.split()[1].split('.')[0]) for k, ent in issues if k == 'loop_winding'}
    sound = [si for si, s in enumerate(m.shells) if si not in bad_shells and not (set(s) & bad_faces)]
    oriented = []
    for si in sound:
        volume = shell_volume(m, m.shells[si])
        if abs(volume) < mp.mpf('1e-25'):
            raise ArithmeticError(f'{m.name}: shell volume too small for the oracle')
        if (volume > 0) != (si == 0):
            issues.append(issue('shell_orientation', f'shell {si}'))
        else:
            oriented.append(si)
    if 0 in oriented:
        cavities = [si for si in oriented if si != 0]
        for si in cavities:
            first = m.faces[m.shells[si][0]].loops[0][0]
            if not vertex_ok[use_vertices(m, first)[0]]:
                continue
            point = vec(m.vertices[use_vertices(m, first)[0]])
            if not inside(m, m.shells[0], point):
                issues.append(issue('cavity_outside', f'shell {si}'))
            elif any(inside(m, m.shells[sj], point) for sj in cavities if sj != si):
                issues.append(issue('nested_cavity', f'shell {si}'))
    return sorted(set(issues))


def margin_check(value, tol, what):
    """Refuse a decision within 0.1% of the tolerance."""
    if abs(value-tol) < tol*mp.mpf('1e-3'):
        raise ArithmeticError(f'{what}: too close to tolerance for the oracle')


def loop_area(loop, steps=2048):
    """Signed area 1/2 * integral(u dv - v du) by Gauss-Legendre per use, with
    the loop closed by straight chords between consecutive uses."""
    total = mp.mpf(0)
    for k, u in enumerate(loop):
        p = u.pcurve
        f = lambda t: area_integrand(p, t)
        total += mp.quad(f, [0, 1])
        a, b = pcurve_point(p, 1), pcurve_point(loop[(k+1) % len(loop)].pcurve, 0)
        total += a[0]*b[1]-b[0]*a[1]
    return total/2


def area_integrand(p, t):
    if isinstance(p, Line2):
        a, b = vec(p.start), vec(p.end)
        q = add(a, mul(sub(b, a), t))
        d = sub(b, a)
        return q[0]*d[1]-q[1]*d[0]
    a = mp.mpf(p.start)+mp.mpf(p.sweep)*t
    c = vec(p.center)
    q = [c[0]+p.radius*mp.cos(a), c[1]+p.radius*mp.sin(a)]
    d = [-p.radius*mp.sin(a)*p.sweep, p.radius*mp.cos(a)*p.sweep]
    return q[0]*d[1]-q[1]*d[0]


def winding(loop, point):
    """Winding number of a closed UV loop around point, by angle integration."""
    total = mp.mpf(0)
    for u in loop:
        def dtheta(t, p=u.pcurve):
            q = pcurve_point(p, t)
            if isinstance(p, Line2):
                d = sub(vec(p.end), vec(p.start))
            else:
                a = mp.mpf(p.start)+mp.mpf(p.sweep)*t
                d = [-p.radius*mp.sin(a)*p.sweep, p.radius*mp.cos(a)*p.sweep]
            r = sub(q, point)
            return (r[0]*d[1]-r[1]*d[0])/dot(r, r)
        total += mp.quad(dtheta, [0, mp.mpf(1)/2, 1])
    return int(mp.nint(total/(2*mp.pi)))


def face_normal(face, uv):
    s = face.surface
    o, x, y, n = axes(s.frame)
    if isinstance(s, Plane):
        base = n
    else:
        base = add(mul(x, mp.cos(uv[0])), mul(y, mp.sin(uv[0])))
    return base if face.forward else mul(base, -1)


def shell_volume(m, shell):
    """(1/3) sum over faces of integral x.n dA, via Green's theorem in UV."""
    total = mp.mpf(0)
    for fi in shell:
        face = m.faces[fi]
        s = face.surface
        o, x, y, n = axes(s.frame)
        # A reversed face's loops wind clockwise in UV, so the boundary
        # integral already carries the reversed normal's sign.
        # integral over the UV region of g(u,v) = (S(u,v).N(u)) * |S_u x S_v|
        # = boundary integral of G(u,v) dv with dG/du = g (region counted by
        # loop winding, positive for the forward-face outer loop).
        if isinstance(s, Plane):
            h = dot(o, n)
            G = lambda u, v: h*u  # g = h, |S_u x S_v| = 1
        else:
            r = mp.mpf(s.radius)
            ox, oy = dot(o, x), dot(o, y)
            # g = (ox cos u + oy sin u + r) * r
            G = lambda u, v, r=r, ox=ox, oy=oy: r*(ox*mp.sin(u)-oy*mp.cos(u)+r*u)
        for loop in face.loops:
            for u in loop:
                p = u.pcurve
                def integrand(t, p=p, G=G):
                    q = pcurve_point(p, t)
                    if isinstance(p, Line2):
                        dv = p.end[1]-p.start[1]
                    else:
                        a = mp.mpf(p.start)+mp.mpf(p.sweep)*t
                        dv = p.radius*mp.cos(a)*p.sweep
                    return G(q[0], q[1])*dv
                total += mp.quad(integrand, [0, mp.mpf(1)/2, 1])
    return total/3


def inside(m, shell, point):
    """Ray parity against the shell's faces, in high precision."""
    direction = [mp.mpf('0.5773502691896257'), mp.mpf('0.6123724356957945'), mp.mpf('0.5400617248673217')]
    hits = 0
    for fi in shell:
        face = m.faces[fi]
        s = face.surface
        o, x, y, n = axes(s.frame)
        rel = sub(point, o)
        # One entry per ray parameter; each lists the UV aliases of that hit.
        crossings = []
        if isinstance(s, Plane):
            den = dot(direction, n)
            if den != 0:
                t = -dot(rel, n)/den
                if t > 0:
                    q = add(rel, mul(direction, t))
                    crossings.append([[dot(q, x), dot(q, y)]])
        else:
            pr = sub(rel, mul(n, dot(rel, n)))
            pd = sub(direction, mul(n, dot(direction, n)))
            a, b, c = dot(pd, pd), 2*dot(pr, pd), dot(pr, pr)-mp.mpf(s.radius)**2
            disc = b*b-4*a*c
            if a != 0 and disc > 0:
                for t in ((-b-mp.sqrt(disc))/(2*a), (-b+mp.sqrt(disc))/(2*a)):
                    if t > 0:
                        q = add(rel, mul(direction, t))
                        angle, v = mp.atan2(dot(q, y), dot(q, x)), dot(q, n)
                        crossings.append([[angle+2*mp.pi*k, v] for k in (-1, 0, 1)])
        for aliases in crossings:
            if any(sum(winding(loop, uv) for loop in face.loops) for uv in aliases):
                hits += 1
    return hits % 2 == 1


# ---------------------------------------------------------------- OCCT encoding

def representable(m):
    """Reference errors cannot be built as OCCT shapes."""
    return not any(k == 'reference' for k, _ in validate_references(m))


def validate_references(m):
    out = []
    nv, ne, nf = len(m.vertices), len(m.edges), len(m.faces)
    for i, e in enumerate(m.edges):
        if not (0 <= e.start < nv and 0 <= e.end < nv):
            out.append(('reference', f'edge {i}'))
    for f in m.faces:
        for loop in f.loops:
            for u in loop:
                if not 0 <= u.edge < ne:
                    out.append(('reference', 'use'))
    for s in m.shells:
        if any(not 0 <= fi < nf for fi in s):
            out.append(('reference', 'shell'))
    return out


def edge_range(c):
    """OCCT curve, parameter range, and our fraction -> OCCT parameter."""
    if isinstance(c, Line3):
        d = [b-a for a, b in zip(c.start, c.end)]
        length = math.sqrt(ltr_sum(x*x for x in d))
        unit = [x/length for x in d] if length else [1.0, 0.0, 0.0]
        return ('line', (*c.start, *unit)), 0.0, length
    f = c.frame
    if c.sweep >= 0:
        return ('circle', (*f.origin, *f.normal, *f.x, c.radius)), c.start, c.start+c.sweep
    flipped = tuple(-x for x in f.normal)
    return ('circle', (*f.origin, *flipped, *f.x, c.radius)), -c.start, -c.start-c.sweep


def pcurve_encoding(p, forward, first, last):
    """Geom2d curve whose parameter equals the edge parameter; exact flag."""
    span = last-first
    if isinstance(p, Line2):
        a, b = (p.start, p.end) if forward else (p.end, p.start)
        d = (b[0]-a[0], b[1]-a[1])
        length = hypot_rn(*d)
        unit = (d[0]/length, d[1]/length) if length else (1.0, 0.0)
        origin = (a[0]-first*unit[0], a[1]-first*unit[1])
        exact = span > 0 and abs(length-span) <= 1e-12*max(1.0, span)
        return ('line', (*origin, *unit)), exact
    speed = p.sweep/span if span else 0.0
    if forward:
        sense, beta = (1 if speed > 0 else -1), p.start
    else:
        sense, beta = (-1 if speed > 0 else 1), p.start+p.sweep
    beta -= sense*first
    exact = span > 0 and abs(abs(speed)-1) <= 1e-12
    return ('circle', (*p.center, cos_rn(beta), sin_rn(beta), p.radius, sense)), exact


def native(m):
    """Explicit OCCT construction rows for occt_brep_check_oracle.cpp."""
    out = [f'case {m.name} {number(m.tolerance)}']
    ranges = []
    for v in m.vertices:
        out.append('v '+' '.join(map(number, v)))
    for e in m.edges:
        (kind, values), first, last = edge_range(e.curve)
        ranges.append((first, last))
        out.append(f'e {e.start} {e.end} {kind} '+' '.join(map(number, values))+f' {number(first)} {number(last)}')
    inexact = []
    for fi, f in enumerate(m.faces):
        s = f.surface
        values = (*s.frame.origin, *s.frame.normal, *s.frame.x) + (() if isinstance(s, Plane) else (s.radius,))
        kind = 'plane' if isinstance(s, Plane) else 'cylinder'
        out.append(f'f {kind} '+' '.join(map(number, values))+(' F' if f.forward else ' R'))
        for li, loop in enumerate(f.loops):
            out.append('w')
            for ui, u in enumerate(loop):
                first, last = ranges[u.edge]
                (pk, pv), exact = pcurve_encoding(u.pcurve, u.forward, first, last)
                if not exact:
                    inexact.append(f'{fi}.{li}.{ui}')
                # OCCT orientation inside the underlying FORWARD face.
                o = 'F' if u.forward == f.forward else 'R'
                out.append(f'u {u.edge} {o} {pk} '+' '.join(map(number, pv)))
    for s in m.shells:
        out.append('s '+' '.join(map(str, s)))
    out.append('end')
    return '\n'.join(out), inexact


# ------------------------------------------------------------------ native comparison

# BRepCheck_Status values (BRepCheck_Status.hxx order) that correspond to each
# issue class of the Rust contract. An issue class passes when the native
# analyzer reports any of its statuses anywhere in the solid. Classes with no
# BRepCheck counterpart compare by verdict only.
CORRESPONDING = {
    'vertex_off_curve': {1},
    'pcurve_off_edge': {8, 9},
    'uv_gap': {28},
    'free_edge': {28},
    'open_loop': {28, 32},
    'same_sense_uses': {32},
    'loop_winding': {32},
    'non_manifold_edge': {14},
    'disconnected_shell': {29},
    'face_reused': {26, 29},
    'empty_shell': {24},
    'empty_loop': {16},
    'inner_loop_outside': {23},
    'shell_orientation': {30, 35},
    'cavity_outside': {30},
    'nested_cavity': {30},
}


def decode_native(line, name):
    """(valid, {label: statuses}, absent labels) from one probe output row."""
    words = line.split()
    if len(words) < 3 or words[0] != name or words[1] != 'R' or words[2] not in ('0', '1'):
        raise ValueError(f'malformed native row for {name}: {line!r}')
    statuses, absent = {}, set()
    for token in words[3:]:
        label, *codes = token.split(':')
        if codes == ['absent']:
            absent.add(label)
        elif codes and all(c.isdigit() and 0 < int(c) < 37 for c in codes):
            statuses[label] = {int(c) for c in codes}
        else:
            raise ValueError(f'malformed native status {token!r}')
    return words[2] == '1', statuses, absent


def structure_only(m):
    """Native labels of OCCT structure the seamless model does not have: seam
    edges (used twice by one loop of a face) and vertices used only by seams
    and by edges closed on them."""
    seams = set()
    for f in m.faces:
        for loop in f.loops:
            used = [u.edge for u in loop]
            seams |= {e for e in used if used.count(e) == 2}
    labels = {f'e{e}' for e in seams}
    for v in range(len(m.vertices)):
        incident = [k for k, e in enumerate(m.edges) if v in (e.start, e.end)]
        if (incident and any(k in seams for k in incident)
                and all(k in seams or m.edges[k].start == m.edges[k].end for k in incident)):
            labels.add(f'v{v}')
    return labels


def compare_native(issues, native, inexact, structure=frozenset()):
    """Sorted difference classes between a Rust issue list and a native row.
    Statuses on structure-only labels cannot stand for a Rust issue class."""
    valid, statuses, _ = native
    located = [codes for label, codes in statuses.items() if label not in structure]
    reported = set().union(*located) if located else set()
    differences = set()
    if valid != (not issues):
        differences.add('verdict')
    for kind in {i.split(':', 1)[0] for i in issues}:
        if kind in CORRESPONDING and not CORRESPONDING[kind] & reported:
            differences.add('unmatched_'+kind)
    if inexact:
        # A pcurve whose parameter speed OCCT cannot match is not a verdict.
        differences.add('inexact_pcurve_encoding')
    return sorted(differences)
