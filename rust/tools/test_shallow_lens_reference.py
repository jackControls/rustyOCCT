"""Two-way checks of shallow_lens_reference.py (the shallow lens, HANDOFF item
5 (d)): its closed-form inner integral against quadrature, its sphere's part
by a second rule and, in the level frame, by an independent formula (the
sphere as two graphs over the axis's plane), the shallow lenses against
their paraboloid asymptotics, and the frames against `Frame3::new`'s stored
axes the Rust kernel reports (`tests/shallow_lens_booleans.rs` checks the
fixture's bits)."""
import unittest

import mpmath as mp

import generate_shallow_lens_fixtures as gen
import shallow_lens_reference as ref

DPS = 60


def lens(frame_name, d):
    for name, origin, normal, hint in gen.FRAMES:
        if name == frame_name:
            frame = ref.frame_axes(origin, normal, hint)
            centre = ref.world64(*frame, gen.RADIUS+gen.BALL*(1.0-d), 0.0, gen.HEIGHT)
            return ref.Lens(frame, gen.RADIUS, centre, gen.BALL)
    raise KeyError(frame_name)


def close(a, b, relative):
    return abs(a-b) <= relative*max(abs(a), abs(b))


class ShallowLensReference(unittest.TestCase):
    def test_frames_are_the_kernels(self):
        o, x, y, n = ref.frame_axes((1.0, -2.0, 0.5), (0.0, 3.0, 4.0), (1.0, 0.0, 0.0))
        # Rust's Frame3::new on macOS and under the emulated glibc hypot.
        self.assertEqual(n, (0.0, 0.6, 0.8))
        self.assertEqual(y, (0.0, 0.8, -0.5999999999999999))
        self.assertEqual(x, (1.0, -0.0, 0.0))
        self.assertEqual(ref.at_height(o, n, -3.0), (1.0, -3.8, -1.9000000000000004))

    def test_inner_integral_is_its_quadrature(self):
        with mp.workdps(DPS):
            for frame_name, d in [('level', 1e-6), ('turned', 1.4e-7), ('turned', 0.25)]:
                lens_ = lens(frame_name, d)
                ta, tb = lens_.ends()
                for k in (1, 3, 5):
                    t = ta+(tb-ta)*k/6

                    def f(s):
                        a2, b2, c2 = lens_.quadratic(t)
                        q = a2*s*s+b2*s+c2
                        return 2*mp.sqrt(max(q, 0))*s/lens_.a
                    a2, b2, c2 = lens_.quadratic(t)
                    s1 = (-b2+mp.sqrt(b2*b2-4*a2*c2))/(2*a2)
                    numeric = mp.quad(f, [s1, lens_.big])
                    self.assertTrue(close(lens_.inner(t), numeric, mp.mpf(10)**-30),
                                    (frame_name, d, k))

    def test_volume_and_wall_two_rules(self):
        """The outer integrals again over four pieces (each end's power
        singularity alone in its piece)."""
        with mp.workdps(DPS):
            for frame_name, d in [('turned', 1.4e-7), ('level', 1e-5), ('turned', 0.05)]:
                lens_ = lens(frame_name, d)
                ta, tb = lens_.ends()
                jac = abs(ref.det3(lens_.x, lens_.y, lens_.n))
                pieces = mp.linspace(ta, tb, 5)
                volume = jac*mp.quad(lens_.inner, pieces)
                self.assertTrue(close(lens_.lens_volume(), volume, mp.mpf(10)**-32), (frame_name, d))

                def wall(t):
                    q = lens_.disc_wall(t)
                    return lens_.big*ref.norm(ref.cross(lens_.de(t), lens_.n))*2*mp.sqrt(max(q, 0))/lens_.a
                self.assertTrue(close(lens_.wall_area(), mp.quad(wall, pieces), mp.mpf(10)**-32),
                                (frame_name, d))

    def test_sphere_part_two_rules(self):
        with mp.workdps(DPS):
            for frame_name, d in [('turned', 1e-6), ('level', 0.25)]:
                lens_ = lens(frame_name, d)
                trapezoid = lens_.sphere_area(ref.POINTS)
                finer = lens_.sphere_area(2*ref.POINTS)
                quad = lens_.sphere_area()
                self.assertTrue(close(trapezoid, finer, mp.mpf(10)**-32), (frame_name, d))
                self.assertTrue(close(trapezoid, quad, mp.mpf(10)**-32), (frame_name, d))

    def test_level_sphere_part_as_two_graphs(self):
        """In the level frame the sphere's two halves are graphs over the
        axis's plane: `2 r / sqrt(r^2 - rho^2)` over the disc of the rod and
        the ball's shadow, its inner integral in `s` in closed form."""
        with mp.workdps(DPS):
            for d in (1e-5, 1.4e-7, 0.25):
                lens_ = lens('level', d)
                cx, cy = lens_.c[0], lens_.c[1]
                dist = mp.sqrt(cx*cx+cy*cy)
                r, big = lens_.r, lens_.big

                def inner(t):
                    a2 = r*r-dist*dist*mp.sin(t)**2
                    if a2 <= 0:
                        return mp.mpf(0)
                    a = mp.sqrt(a2)
                    lo = dist*mp.cos(t)-a
                    hi = min(big, dist*mp.cos(t)+a)
                    if hi <= lo:
                        return mp.mpf(0)

                    def g(s):
                        w = s-dist*mp.cos(t)
                        arg = max(min(w/a, mp.mpf(1)), mp.mpf(-1))
                        return 2*r*(-mp.sqrt(max(a2-w*w, 0))+dist*mp.cos(t)*mp.asin(arg))
                    return g(hi)-g(lo)
                ta, tb = lens_.ends()
                graphs = mp.quad(inner, [ta, tb])
                self.assertTrue(close(graphs, lens_.sphere_area(ref.POINTS), mp.mpf(10)**-30), d)

    def test_shallow_lenses_follow_the_paraboloids(self):
        """A lens `delta` deep: `pi delta^2 sqrt(a b)` and area `4 pi delta
        sqrt(a b)`, `a = r` along the axis, `b = r R / (r + R)` across it,
        within `10 delta / r` relatively."""
        with mp.workdps(DPS):
            for frame_name in ('level', 'turned'):
                for d in (1e-5, 1e-6, 3e-7, 1.4e-7):
                    lens_ = lens(frame_name, d)
                    m = ref.measures(lens_)
                    delta = mp.mpf(gen.BALL)*mp.mpf(d)
                    a, b = lens_.r, lens_.r*lens_.big/(lens_.r+lens_.big)
                    volume = mp.pi*delta**2*mp.sqrt(a*b)
                    area = 4*mp.pi*delta*mp.sqrt(a*b)
                    bound = 10*d
                    self.assertTrue(close(m['lens_volume'], volume, bound), (frame_name, d))
                    self.assertTrue(close(m['lens_area'], area, bound), (frame_name, d))

    def test_bodies(self):
        with mp.workdps(DPS):
            level = lens('level', 1e-6)
            self.assertTrue(close(level.rod_volume(), 6*mp.pi, mp.mpf(10)**-35))
            self.assertTrue(close(level.rod_area(), 14*mp.pi, mp.mpf(10)**-35))
            turned = lens('turned', 1e-6)
            # The turned axes are orthonormal but for rounding.
            self.assertTrue(close(turned.rod_volume(), 6*mp.pi, mp.mpf(10)**-15))
            self.assertTrue(close(turned.rod_area(), 14*mp.pi, mp.mpf(10)**-15))


if __name__ == '__main__':
    unittest.main()
