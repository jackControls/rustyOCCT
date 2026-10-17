"""The exact volume and areas of a ball crossing a rod's wall (HANDOFF item 5
(d), "The shallow lens refined" in REVIEW_NOTES.md): the lens the two share,
from the bodies' stored binary64 data, computed in mpmath to 30 digits.

The rod is `Solid::cylinder_with(frame, R, -3, 3)`: its wall the affine
cylinder `o' + l0 x + l1 y + l2 n` with `l0^2 + l1^2 = R^2` on the frame's
stored axes (`Frame3::new`'s, emulated in `frame_axes`: each direction
scaled by its largest coordinate and divided by its nested correctly
rounded `hypot`), `o'` the stored bottom cap's origin `o - 3 n` (binary64
operations, as the kernel rounds them); its caps the planes through `o'`
and `o + 3 n` normal to `n`. The ball is `Solid::sphere_with` about a
binary64 centre `c` in a level frame, radius `r`: `|p - c| = r` exactly.

The lens's volume, in the cylinder's affine coordinates `(s cos t, s sin t,
l2)`: for each `t` the ball's chord along `n` through the point at radius
`s` has length `2 sqrt(D(s)) / (n . n)`, `D(s) = A s^2 + B s + C` the
quadratic's discriminant (`A < 0`), so the volume is `|det(x, y, n)|`
times the integral over `t` of `integral from s1 to R of 2 sqrt(D) s / (n .
n) ds` (`s1` the smaller root of `D`), the inner integral in closed form
(`(m^2 - k^2 w^2)^(3/2)` and `asin`), the outer by tanh-sinh quadrature
between the angles where `D(R) = 0` (the lens's ends on the wall), where
the inner integral vanishes like a power three halves. The wall's part of
the lens's area: `R |e'(t) x n| 2 sqrt(D(R)) / (n . n)` over the same
angles. The sphere's part: `r^2` times the integral over `phi` of `1 - cos
psi(phi)` (written `2 sin^2(psi / 2)`), `psi(phi)` the polar angle about
the direction from `c` to the axis where the sphere's point leaves the
cylinder, found to the working precision (bisection, then a bracketing
secant), integrated by the periodic trapezoid rule; the region is
star-shaped about that direction for every case here (a crossing shallower
than the ball's radius). No `math.cos`, `sin` or `hypot`: binary64 values
are rounded once from mpmath.
"""
import mpmath as mp

from brep_reference import hypot_rn

DPS = 60
DIGITS = 30
# The sphere's part's periodic trapezoid rule: its integrand is analytic in
# `phi`, so the rule converges geometrically (128 points agree with
# tanh-sinh quadrature to every digit kept, `test_shallow_lens_reference`).
POINTS = 128


def unit(v):
    """Vec3::normalized: scaled by the largest coordinate, then divided by
    the nested correctly rounded hypot (the emulated glibc's)."""
    s = max(abs(x) for x in v)
    w = tuple(x/s for x in v)
    length = hypot_rn(hypot_rn(w[0], w[1]), w[2])
    return tuple(x/length for x in w)


def cross64(a, b):
    """A binary64 cross product, as Vec3::cross rounds it."""
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])


def frame_axes(origin, normal, hint):
    """Frame3::new's stored axes: unit normal, y = normal x hint, x = y x
    normal, each normalized."""
    n = unit(normal)
    h = unit(hint)
    y = unit(cross64(n, h))
    x = unit(cross64(y, n))
    return tuple(origin), x, y, n


def at_height(origin, n, w):
    """`Frame3::point(0, 0, w)`'s binary64 origin `o + w n`."""
    return tuple(origin[i]+w*n[i] for i in range(3))


def mpv(v):
    return [mp.mpf(x) for x in v]


def dot(a, b):
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]


def sub(a, b):
    return [a[i]-b[i] for i in range(3)]


def add(a, b):
    return [a[i]+b[i] for i in range(3)]


def scale(a, k):
    return [x*k for x in a]


def cross(a, b):
    return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]


def norm(a):
    return mp.sqrt(dot(a, a))


def det3(a, b, c):
    return dot(a, cross(b, c))


def world64(origin, x, y, n, u, v, w):
    """A point `o + u x + v y + w n` of a frame, exact, rounded once."""
    o, x, y, n = mpv(origin), mpv(x), mpv(y), mpv(n)
    p = add(add(o, scale(x, mp.mpf(u))), add(scale(y, mp.mpf(v)), scale(n, mp.mpf(w))))
    return tuple(float(t) for t in p)


class Lens:
    """The rod (stored frame, radius `R`, heights `-3` to `3`) and the ball
    (centre `c`, radius `r`): the lens's volume and areas, and the bodies'
    own."""

    def __init__(self, frame, radius, centre, ball_radius, low=-3.0, high=3.0):
        origin, x, y, n = frame
        self.bottom = at_height(origin, n, low)
        self.top = at_height(origin, n, high)
        self.o = mpv(self.bottom)
        self.x, self.y, self.n = mpv(x), mpv(y), mpv(n)
        self.big = mp.mpf(radius)
        self.c = mpv(centre)
        self.r = mp.mpf(ball_radius)
        self.w0 = sub(self.o, self.c)
        self.a = dot(self.n, self.n)
        # Local coordinates on the stored axes: the inverse's rows.
        det = det3(self.x, self.y, self.n)
        self.rows = [scale(cross(self.y, self.n), 1/det), scale(cross(self.n, self.x), 1/det),
                     scale(cross(self.x, self.y), 1/det)]

    def e(self, t):
        return add(scale(self.x, mp.cos(t)), scale(self.y, mp.sin(t)))

    def de(self, t):
        return add(scale(self.x, -mp.sin(t)), scale(self.y, mp.cos(t)))

    def quadratic(self, t):
        """`D(s) = A s^2 + B s + C` at angle `t`."""
        e, n, w0, a = self.e(t), self.n, self.w0, self.a
        ne, nw = dot(n, e), dot(n, w0)
        big_a = ne*ne-a*dot(e, e)
        big_b = 2*nw*ne-2*a*dot(w0, e)
        big_c = nw*nw-a*(dot(w0, w0)-self.r*self.r)
        return big_a, big_b, big_c

    def disc_wall(self, t):
        a2, b2, c2 = self.quadratic(t)
        return a2*self.big*self.big+b2*self.big+c2

    def local(self, p):
        d = sub(p, self.o)
        return [dot(row, d) for row in self.rows]

    def centre_angle(self):
        """The ball centre's angle in the cylinder's coordinates."""
        local = self.local(self.c)
        return mp.atan2(local[1], local[0])

    def ends(self):
        """The angles of the lens's ends on the wall: `D(R) = 0` either side
        of the centre's angle."""
        tc = self.centre_angle()
        if self.disc_wall(tc) <= 0:
            raise ValueError('the ball does not cross the wall')
        out = []
        for side in (-1, 1):
            lo, hi = tc, tc+side*mp.mpf(1)
            if self.disc_wall(hi) >= 0:
                raise ValueError('no end within a radian')
            out.append(bisect(self.disc_wall, lo, hi))
        return out[0], out[1]

    def inner(self, t):
        """`integral from s1 to R of 2 sqrt(D(s)) s / (n . n) ds`, closed
        form: with `D = m^2 - k^2 (s - s0)^2` and `z = k (s - s0) / m`, the
        integral from the lower root (`z = -1`) is `-m^3 (1 - z^2)^(3/2) /
        (3 k^2) + s0 m^2 / (2 k) (z sqrt(1 - z^2) + acos(-z))`, `1 - z^2`
        as `(1 - z) (1 + z)` and `acos(-z)` as `2 asin(sqrt((1 + z) / 2))`,
        so that a lens's end, `z` near `-1`, keeps every digit."""
        a2, b2, c2 = self.quadratic(t)
        k2 = -a2
        k = mp.sqrt(k2)
        s0 = b2/(2*k2)
        m2 = c2+b2*b2/(4*k2)
        if m2 <= 0:
            return mp.mpf(0)
        m = mp.sqrt(m2)
        if self.big >= s0+m/k:
            one_plus, one_minus = mp.mpf(2), mp.mpf(0)
        else:
            one_plus = (m+k*(self.big-s0))/m
            one_minus = (m-k*(self.big-s0))/m
        if one_plus <= 0:
            return mp.mpf(0)
        z = one_plus-1
        root = mp.sqrt(one_plus*one_minus)
        from_end = 2*mp.asin(mp.sqrt(one_plus/2))
        value = -m**3*root**3/(3*k2)+s0*m2/(2*k)*(z*root+from_end)
        return 2/self.a*value

    def lens_volume(self):
        ta, tb = self.ends()
        jac = abs(det3(self.x, self.y, self.n))
        return jac*mp.quad(self.inner, [ta, tb])

    def wall_area(self):
        ta, tb = self.ends()

        def f(t):
            d = self.disc_wall(t)
            if d <= 0:
                return mp.mpf(0)
            return self.big*norm(cross(self.de(t), self.n))*2*mp.sqrt(d)/self.a
        return mp.quad(f, [ta, tb])

    def inside(self, p):
        """`l0^2 + l1^2 - R^2` of a point (negative inside the cylinder)."""
        local = self.local(p)
        return local[0]**2+local[1]**2-self.big**2

    def polar_axes(self):
        """The unit direction `q` from `c` toward the axis (perpendicular to
        `n`), and two unit directions across it."""
        nn = scale(self.n, 1/norm(self.n))
        d = sub(self.c, self.o)
        perp = sub(d, scale(nn, dot(d, nn)))
        q = scale(perp, -1/norm(perp))
        e1 = nn
        e2 = cross(q, e1)
        return q, e1, e2

    def sphere_area(self, points=None):
        q, e1, e2 = self.polar_axes()

        def boundary(phi):
            direction = add(scale(e1, mp.cos(phi)), scale(e2, mp.sin(phi)))

            def g(psi):
                p = add(self.c, scale(add(scale(q, mp.cos(psi)), scale(direction, mp.sin(psi))), self.r))
                return self.inside(p)
            hi = mp.mpf(2)**-30
            while g(hi) < 0:
                hi *= 2
                if hi > mp.pi/2:
                    raise ValueError('the sphere\'s part is not a cap')
            return bisect(g, mp.mpf(0), hi)

        def f(phi):
            return 2*mp.sin(boundary(phi)/2)**2
        if points is None:
            total = mp.quad(f, mp.linspace(0, 2*mp.pi, 9))
        else:
            # The periodic trapezoid rule (the two-way check).
            h = 2*mp.pi/points
            total = h*mp.fsum(f(h*j) for j in range(points))
        return self.r*self.r*total

    def rod_volume(self):
        """`|det| pi R^2` times the caps' distance along `n` over `n . n`."""
        jac = abs(det3(self.x, self.y, self.n))
        height = dot(sub(mpv(self.top), self.o), self.n)/self.a
        return jac*mp.pi*self.big*self.big*height

    def rod_area(self):
        """The wall (`v` from the bottom cap to the top one at each angle)
        and the two caps (`pi R^2 |x' x y'|`, `x'`, `y'` the axes' parts
        across `n`)."""
        a = self.a
        top = dot(sub(mpv(self.top), self.o), self.n)

        def wall(t):
            e = self.e(t)
            v_top = (top-self.big*dot(e, self.n))/a
            v_bottom = (-self.big*dot(e, self.n))/a
            return self.big*norm(cross(self.de(t), self.n))*(v_top-v_bottom)
        side = mp.quad(wall, mp.linspace(0, 2*mp.pi, 5))
        xp = sub(self.x, scale(self.n, dot(self.x, self.n)/a))
        yp = sub(self.y, scale(self.n, dot(self.y, self.n)/a))
        cap = mp.pi*self.big*self.big*norm(cross(xp, yp))
        return side+2*cap

    def ball_volume(self):
        return 4*mp.pi*self.r**3/3

    def ball_area(self):
        return 4*mp.pi*self.r**2


def bisect(f, lo, hi):
    """A root of `f` between `lo` and `hi` (opposite signs) to the working
    precision: forty halvings, then Anderson-Bjorck's bracketing secant."""
    flo = f(lo)
    for _ in range(40):
        mid = (lo+hi)/2
        fm = f(mid)
        if fm == 0:
            return mid
        if (fm < 0) == (flo < 0):
            lo, flo = mid, fm
        else:
            hi = mid
    return mp.findroot(f, (lo, hi), solver='anderson')


def text(x):
    """A reference value's 30 significant digits."""
    return mp.nstr(x, DIGITS, min_fixed=1, max_fixed=0)


def measures(lens):
    """Every reference of a case, in the fixture's order."""
    with mp.workdps(DPS):
        v = lens.lens_volume()
        wall = lens.wall_area()
        sph = lens.sphere_area(POINTS)
        return {
            'lens_volume': v,
            'lens_area': wall+sph,
            'lens_wall_area': wall,
            'lens_sphere_area': sph,
            'rod_volume': lens.rod_volume(),
            'rod_area': lens.rod_area(),
            'ball_volume': lens.ball_volume(),
            'ball_area': lens.ball_area(),
        }
