"""Independent reference for revolved primitives (S3 of REVIEW_NOTES.md): the
right circular cone and frustum, as `BRepPrimAPI_MakeCone(gp_Ax2, r1, r2, h)`
and the kernel's cone builder make them.

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
