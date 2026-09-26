"""Independent reference for revolved primitives (S3 of REVIEW_NOTES.md): the
right circular cone and frustum, as `BRepPrimAPI_MakeCone(gp_Ax2, r1, r2, h)`
and the kernel's cone builder make them, the sphere and spherical zone,
as `BRepPrimAPI_MakeSphere(gp_Ax2, R, angle1, angle2)` makes them, and the
torus, as `BRepPrimAPI_MakeTorus(gp_Ax2, R1, R2, angle1, angle2, angle)`.

Everything is exact geometry in mpmath from the specification, never from
either implementation:

* mass properties: volume, surface area, centre of mass and the matrix of
  inertia about it, by integrating the disc moments along the axis
  (polynomials in the height, integrated exactly);
* OCCT's structure, from the source review of `BRepPrim_OneAxis` and
  `BRepPrim_Cone`: the lateral face has one seam along the frame's x axis
  (angle 0), a zero radius end is an apex with a degenerated edge, each end
  circle is closed at its seam vertex. So an apex cone has 2 vertices, 3
  edges (base circle, seam, degenerated), 2 wires, 2 faces; a frustum 2
  vertices, 3 edges, 3 wires, 3 faces; one shell and one solid;
* every face (type, area, centre of its surface), edge (degenerated or not,
  closed or not, length, the point at its middle parameter) and vertex,
  unordered, to be matched by geometry.

For the sphere, from the source review of `BRepPrim_Sphere` and
`BRepPrim_OneAxis`: the lateral face is a `Geom_SphericalSurface` on the
frame, latitudes `angle1 <= v <= angle2` (radians); its meridian is the
circle about `-y` through `x` and the axis, offset by `2 pi`, so the seam
edge runs along `x` from `angle1` to `angle2` with its middle at the mean
latitude. An end at latitude `±pi/2` (the binary64 value, as DRAW converts
`±90` degrees) is a pole: one vertex and a degenerated edge; any other end is
a circle closed at its seam vertex, bounding a planar disc. So every sphere
or zone has 2 vertices and 3 edges, `1 + caps` wires and faces, one shell
and one solid.
"""
from dataclasses import dataclass

import mpmath as mp

mp.mp.dps = 40


@dataclass
class Cone:
    name: str
    origin: tuple
    normal: tuple
    x: tuple
    r1: float
    r2: float
    height: float


def vec(a):
    return [mp.mpf(x) for x in a]


def add(a, b):
    return [x+y for x, y in zip(a, b)]


def sub(a, b):
    return [x-y for x, y in zip(a, b)]


def mul(a, s):
    return [x*s for x in a]


def dot(a, b):
    return sum((x*y for x, y in zip(a, b)), mp.mpf(0))


def cross(a, b):
    return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]


def norm(a):
    return mp.sqrt(dot(a, a))


def axes(c):
    """gp_Ax2 as OCCT builds it: the normal, and x made orthogonal to it."""
    n = vec(c.normal)
    n = mul(n, 1/norm(n))
    x = vec(c.x)
    x = sub(x, mul(n, dot(x, n)))
    x = mul(x, 1/norm(x))
    return vec(c.origin), x, cross(n, x), n


def radius(c, z):
    return mp.mpf(c.r1)+(mp.mpf(c.r2)-mp.mpf(c.r1))*z/mp.mpf(c.height)


def integrate(f, h):
    """Exact for the polynomials in z used here (degree at most 6)."""
    return mp.quad(f, [0, h])


def mass(c):
    """(volume, area, centroid, inertia matrix about the centroid) in world axes."""
    h = mp.mpf(c.height)
    o, x, y, n = axes(c)
    r = lambda z: radius(c, z)
    volume = integrate(lambda z: mp.pi*r(z)**2, h)
    zc = integrate(lambda z: z*mp.pi*r(z)**2, h)/volume
    slant = mp.sqrt(h**2+(mp.mpf(c.r2)-mp.mpf(c.r1))**2)
    area = mp.pi*(mp.mpf(c.r1)**2+mp.mpf(c.r2)**2)+mp.pi*(mp.mpf(c.r1)+mp.mpf(c.r2))*slant
    # Local frame: the axis moment and the transverse moment about the
    # centroid (discs: pi r^4/2 about their axis, pi r^4/4 about a diameter).
    axial = integrate(lambda z: mp.pi*r(z)**4/2, h)
    transverse = integrate(lambda z: mp.pi*r(z)**4/4+(z-zc)**2*mp.pi*r(z)**2, h)
    centre = add(o, mul(n, zc))
    local = [[transverse, 0, 0], [0, transverse, 0], [0, 0, axial]]
    basis = [x, y, n]
    inertia = [[sum(basis[a][i]*local[a][b]*basis[b][j] for a in range(3) for b in range(3))
                for j in range(3)] for i in range(3)]
    return volume, area, centre, inertia


def counts(c):
    """(vertices, edges, wires, faces, shells, solids) of OCCT's cone."""
    ends = [c.r1 > 0, c.r2 > 0]
    faces = 1+sum(ends)
    return (2, 3, faces, faces, 1, 1)


def faces(c):
    """[(type, area, centre)] of every face."""
    h = mp.mpf(c.height)
    o, x, y, n = axes(c)
    r1, r2 = mp.mpf(c.r1), mp.mpf(c.r2)
    slant = mp.sqrt(h**2+(r2-r1)**2)
    # Lateral surface centroid height: (r1 + 2 r2) h / (3 (r1 + r2)).
    out = [('cone', mp.pi*(r1+r2)*slant, add(o, mul(n, h*(r1+2*r2)/(3*(r1+r2)))))]
    if c.r1 > 0:
        out.append(('plane', mp.pi*r1**2, o))
    if c.r2 > 0:
        out.append(('plane', mp.pi*r2**2, add(o, mul(n, h))))
    return out


def edges(c):
    """[(degenerated, closed, length, point at the middle parameter)]."""
    h = mp.mpf(c.height)
    o, x, y, n = axes(c)
    r1, r2 = mp.mpf(c.r1), mp.mpf(c.r2)
    out = []
    for r, z in ((r1, 0), (r2, h)):
        centre = add(o, mul(n, z))
        if r > 0:
            # Parameter from 0 to 2 pi, starting at the seam on x.
            out.append(('regular', 'closed', 2*mp.pi*r, sub(centre, mul(x, r))))
        else:
            out.append(('degenerated', 'closed', mp.mpf(0), centre))
    a, b = add(o, mul(x, r1)), add(add(o, mul(n, h)), mul(x, r2))
    out.append(('regular', 'open', norm(sub(b, a)), mul(add(a, b), mp.mpf(1)/2)))
    return out


def vertices(c):
    h = mp.mpf(c.height)
    o, x, y, n = axes(c)
    return [add(o, mul(x, mp.mpf(c.r1))), add(add(o, mul(n, h)), mul(x, mp.mpf(c.r2)))]


# ------------------------------------------------------------------ spheres

HALF_PI = 1.5707963267948966


@dataclass
class Sphere:
    name: str
    origin: tuple
    normal: tuple
    x: tuple
    radius: float
    a1: float
    a2: float


def sphere_ends(c):
    """[(latitude, is a pole)] for the lower and upper end."""
    return [(c.a1, c.a1 == -HALF_PI), (c.a2, c.a2 == HALF_PI)]


def sphere_mass(c):
    """(volume, area, centroid, inertia matrix about the centroid)."""
    o, x, y, n = axes(c)
    R = mp.mpf(c.radius)
    z1, z2 = R*mp.sin(mp.mpf(c.a1)), R*mp.sin(mp.mpf(c.a2))
    r2 = lambda z: R**2-z**2
    span = lambda f: mp.quad(f, [z1, z2])
    volume = span(lambda z: mp.pi*r2(z))
    zc = span(lambda z: z*mp.pi*r2(z))/volume
    caps = sum(mp.pi*r2(R*mp.sin(mp.mpf(a))) for a, pole in sphere_ends(c) if not pole)
    area = 2*mp.pi*R*(z2-z1)+caps
    axial = span(lambda z: mp.pi*r2(z)**2/2)
    transverse = span(lambda z: mp.pi*r2(z)**2/4+(z-zc)**2*mp.pi*r2(z))
    centre = add(o, mul(n, zc))
    local = [[transverse, 0, 0], [0, transverse, 0], [0, 0, axial]]
    basis = [x, y, n]
    inertia = [[sum(basis[a][i]*local[a][b]*basis[b][j] for a in range(3) for b in range(3))
                for j in range(3)] for i in range(3)]
    return volume, area, centre, inertia


def sphere_counts(c):
    caps = sum(1 for _, pole in sphere_ends(c) if not pole)
    return (2, 3, 1+caps, 1+caps, 1, 1)


def sphere_faces(c):
    o, x, y, n = axes(c)
    R = mp.mpf(c.radius)
    z1, z2 = R*mp.sin(mp.mpf(c.a1)), R*mp.sin(mp.mpf(c.a2))
    # A zone's area is uniform in height (Archimedes).
    out = [('sphere', 2*mp.pi*R*(z2-z1), add(o, mul(n, (z1+z2)/2)))]
    for a, pole in sphere_ends(c):
        if not pole:
            z = R*mp.sin(mp.mpf(a))
            out.append(('plane', mp.pi*(R**2-z**2), add(o, mul(n, z))))
    return out


def sphere_edges(c):
    o, x, y, n = axes(c)
    R = mp.mpf(c.radius)
    out = []
    for a, pole in sphere_ends(c):
        a = mp.mpf(a)
        centre = add(o, mul(n, R*mp.sin(a)))
        if pole:
            out.append(('degenerated', 'closed', mp.mpf(0), centre))
        else:
            r = R*mp.cos(a)
            out.append(('regular', 'closed', 2*mp.pi*r, sub(centre, mul(x, r))))
    a1, a2 = mp.mpf(c.a1), mp.mpf(c.a2)
    m = (a1+a2)/2
    out.append(('regular', 'open', R*(a2-a1), add(o, add(mul(x, R*mp.cos(m)), mul(n, R*mp.sin(m))))))
    return out


def sphere_vertices(c):
    o, x, y, n = axes(c)
    R = mp.mpf(c.radius)
    return [add(o, add(mul(x, R*mp.cos(mp.mpf(a))), mul(n, R*mp.sin(mp.mpf(a))))) for a, _ in sphere_ends(c)]


# ------------------------------------------------------------------ tori

TWO_PI = 6.283185307179586


@dataclass
class Torus:
    """BRepPrimAPI_MakeTorus(gp_Ax2, major, minor, a1, a2, angle): the minor
    circle about O + major x in the half-plane of x and the axis, latitudes
    a1..a2 of it, revolved by `angle` from x. From the source review of
    `BRepPrim_Torus` and `BRepPrim_OneAxis`: a closed meridian (a2 - a1 a
    full turn) has no top or bottom; otherwise the ends revolve into circles
    bounding planar discs. A partial `angle` adds planar start and end faces
    in the meridian half-planes. So the whole torus has 1 vertex (at v = a1,
    u = 0) and 2 edges (the meridian circle at u = 0 and the latitude circle
    at v = a1, each a seam), 1 wire, 1 face; a v-segment 2 vertices, 3 edges
    (two latitude circles, the meridian arc seam), 3 wires, 3 faces; a wedge 2
    vertices, 3 edges (the meridian circles at u = 0 and u = angle, closed at
    v = a1, and the latitude arc at v = a1 between them), 3 wires, 3 faces."""
    name: str
    origin: tuple
    normal: tuple
    x: tuple
    major: float
    minor: float
    a1: float
    a2: float
    angle: float


def torus_closed(c):
    return c.a2-c.a1 == TWO_PI


def meridian_moments(c, k, m):
    """The integral of rho^k z^m over the meridian region: the minor disc, or
    the region between the arc a1..a2, the segments from its ends to the axis
    and the axis, by Green (the integral of rho^(k+1) z^m / (k+1) dz round
    the boundary, counterclockwise in (rho, z))."""
    R, r = mp.mpf(c.major), mp.mpf(c.minor)
    a1, a2 = mp.mpf(c.a1), mp.mpf(c.a2)

    def arc(k, m):
        # rho = R + r cos v, z = r sin v, dz = r cos v dv, v from a1 to a2.
        # The segments from the arc's ends to the axis have dz = 0 and the
        # axis rho = 0: only the arc contributes.
        G = lambda rho, z: rho**(k+1)*z**m/(k+1)
        return mp.quad(lambda v: G(R+r*mp.cos(v), r*mp.sin(v))*r*mp.cos(v), [a1, (a1+a2)/2, a2])
    # The arc runs counterclockwise about the region when it bulges away from
    # the axis and clockwise when it bulges toward it: the signed area says.
    sign = 1 if arc(0, 0) > 0 else -1
    return sign*arc(k, m)


def torus_mass(c):
    """(volume, area, centroid, inertia about the centroid): a partial
    revolution by alpha of the meridian region, from its moments."""
    o, x, y, n = axes(c)
    al = mp.mpf(c.angle)
    M = lambda k, m: meridian_moments(c, k, m)
    volume = al*M(1, 0)
    sx, sy, sz = mp.sin(al)*M(2, 0), (1-mp.cos(al))*M(2, 0), al*M(1, 1)
    cx, cy, cz = sx/volume, sy/volume, sz/volume
    J = {
        'xx': M(3, 0)*(al/2+mp.sin(2*al)/4), 'yy': M(3, 0)*(al/2-mp.sin(2*al)/4),
        'zz': al*M(1, 2), 'xy': M(3, 0)*mp.sin(al)**2/2,
        'xz': mp.sin(al)*M(2, 1), 'yz': (1-mp.cos(al))*M(2, 1),
    }
    # Second moments about the centroid, then the inertia tensor.
    c_ = [cx, cy, cz]
    names = ['x', 'y', 'z']
    second = [[J[''.join(sorted(names[i]+names[j]))] - volume*c_[i]*c_[j] for j in range(3)] for i in range(3)]
    trace = second[0][0]+second[1][1]+second[2][2]
    local = [[(trace if i == j else 0)-second[i][j] for j in range(3)] for i in range(3)]
    basis = [x, y, n]
    inertia = [[sum(basis[a][i]*local[a][b]*basis[b][j] for a in range(3) for b in range(3))
                for j in range(3)] for i in range(3)]
    centre = add(o, add(add(mul(x, cx), mul(y, cy)), mul(n, cz)))
    area = sum(a for _, a, _ in torus_faces(c))
    return volume, area, centre, inertia


def torus_counts(c):
    if torus_closed(c) and c.angle == TWO_PI:
        return (1, 2, 1, 1, 1, 1)
    return (2, 3, 3, 3, 1, 1)


def torus_faces(c):
    """[(type, area, centre)]: the lateral face by Pappus in the partial
    revolution, the discs of a v-segment's ends, a wedge's end faces."""
    o, x, y, n = axes(c)
    R, r, al = mp.mpf(c.major), mp.mpf(c.minor), mp.mpf(c.angle)
    a1, a2 = mp.mpf(c.a1), mp.mpf(c.a2)
    rho = lambda v: R+r*mp.cos(v)
    z = lambda v: r*mp.sin(v)
    ds = lambda f: mp.quad(lambda v: f(v)*r, [a1, (a1+a2)/2, a2])
    L1, L2, Lz = ds(rho), ds(lambda v: rho(v)**2), ds(lambda v: z(v)*rho(v))
    area = al*L1
    local = [mp.sin(al)*L2/area, (1-mp.cos(al))*L2/area, al*Lz/area]
    out = [('torus', area, add(o, add(add(mul(x, local[0]), mul(y, local[1])), mul(n, local[2]))))]
    if not torus_closed(c):
        for a in (a1, a2):
            out.append(('plane', mp.pi*rho(a)**2, add(o, mul(n, z(a)))))
    if c.angle != TWO_PI:
        e = add(mul(x, mp.cos(al)), mul(y, mp.sin(al)))
        for d in (x, e):
            out.append(('plane', mp.pi*r**2, add(o, mul(d, R))))
    return out


def torus_edges(c):
    """[(degenerated, closed, length, point at the middle parameter)]."""
    o, x, y, n = axes(c)
    R, r, al = mp.mpf(c.major), mp.mpf(c.minor), mp.mpf(c.angle)
    a1, a2 = mp.mpf(c.a1), mp.mpf(c.a2)
    e = lambda u: add(mul(x, mp.cos(u)), mul(y, mp.sin(u)))
    point = lambda u, v: add(o, add(mul(e(u), R+r*mp.cos(v)), mul(n, r*mp.sin(v))))
    out = []
    if c.angle == TWO_PI:
        # The meridian at u = 0: a seam, closed when the meridian is.
        m = (a1+a2)/2
        out.append(('regular', 'closed' if torus_closed(c) else 'open', r*(a2-a1), point(0, m)))
        ends = [a1] if torus_closed(c) else [a1, a2]
        for a in ends:
            # A latitude circle from u = 0, its middle at u = pi.
            out.append(('regular', 'closed', 2*mp.pi*(R+r*mp.cos(a)), point(mp.pi, a)))
    else:
        # The meridian circles at u = 0 and u = angle, closed at v = a1, and
        # the latitude arc at v = a1 between them.
        for u in (mp.mpf(0), al):
            out.append(('regular', 'closed', 2*mp.pi*r, point(u, a1+mp.pi)))
        out.append(('regular', 'open', al*(R+r*mp.cos(a1)), point(al/2, a1)))
    return out


def torus_vertices(c):
    o, x, y, n = axes(c)
    R, r, al = mp.mpf(c.major), mp.mpf(c.minor), mp.mpf(c.angle)
    e = lambda u: add(mul(x, mp.cos(u)), mul(y, mp.sin(u)))
    point = lambda u, v: add(o, add(mul(e(u), R+r*mp.cos(mp.mpf(v))), mul(n, r*mp.sin(mp.mpf(v)))))
    if c.angle == TWO_PI:
        return [point(0, c.a1)] if torus_closed(c) else [point(0, c.a1), point(0, c.a2)]
    return [point(0, c.a1), point(al, c.a1)]
