"""Independent checks of S9c.2's closed forms and fixture list
(`generate_procedural_boolean_fixtures.py`), without the slicing reference."""
import unittest

import mpmath as mp

import generate_procedural_boolean_fixtures as fixtures
from curve_surface_reference import stored_axes


class ProceduralBooleanReferenceTests(unittest.TestCase):
    def setUp(self):
        self.dps = mp.mp.dps
        mp.mp.dps = 40

    def tearDown(self):
        mp.mp.dps = self.dps

    def test_legendre_at_equal_radii_is_steinmetz(self):
        for r in (1, 2, 0.75):
            want = 16*mp.mpf(r)**3/3
            self.assertLess(abs(fixtures.legendre(r, r)-want), mp.mpf(10)**-36*want)

    def test_legendre_agrees_with_the_quadrature(self):
        for rA, rB in ((2.5, 1.5), (2, 0.5), (3, 2.9)):
            P = fixtures.Perpendicular(rA, 0, 0, -50, 50, rB, 0, 0, -50, 50, fixtures.WORLD)
            volume, moments, _ = P.volume()
            want = fixtures.legendre(rA, rB)
            self.assertLess(abs(volume-want), mp.mpf(10)**-35*want)
            self.assertLess(max(abs(m) for m in moments), mp.mpf(10)**-35*want)

    def test_areas_at_equal_radii(self):
        # Each wall's part inside the other cylinder: 8 r^2 (16 r^2 in all).
        P = fixtures.Perpendicular(2, 0, 0, -50, 50, 2, 0, 0, -50, 50, fixtures.WORLD)
        for area in P.areas():
            self.assertLess(abs(area-32), mp.mpf(10)**-35)

    def test_a_pipe_ending_on_the_axis_takes_half(self):
        full = fixtures.Perpendicular(2, 0, 0, -50, 50, 1, 0.5, 0, -50, 50, fixtures.WORLD).volume()[0]
        half = fixtures.Perpendicular(2, 0, 0, -50, 50, 1, 0.5, 0, -50, 0, fixtures.WORLD).volume()[0]
        self.assertLess(abs(2*half-full), mp.mpf(10)**-35*full)

    def test_cases_cover_both_steps_in_the_curved_frames(self):
        listed = fixtures.cases()
        self.assertEqual(len({c.name for c in listed}), len(listed))
        steps = {fixtures.step(c) for c in listed}
        self.assertEqual(steps, {'S9c.2a', 'S9c.2b'})
        for c in listed:
            self.assertNotEqual(stored_axes(c.obj.frame)[1:], stored_axes(c.tool.frame)[1:])
            if c.kind == 'degenerate':
                self.assertTrue(c.reason)
        # The thin pipe is the object of at least one pair.
        self.assertTrue(any(c.obj.boundaries[0].circle[2] < c.tool.boundaries[0].circle[2] for c in listed))


if __name__ == '__main__':
    unittest.main()
