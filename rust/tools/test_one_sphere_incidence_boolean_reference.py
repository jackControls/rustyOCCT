"""Independent checks of S9e.4b.3c.2's fixtures
(`generate_one_sphere_incidence_boolean_fixtures.py`): the cells' closed
forms against direct integrations, the cell reference on single pieces
against their textbook measures and on pieces of one sphere against their
set identities (the hemisphere common the wedge is the wedge, their fuse the
hemisphere), the touching of two octants about one axis, the bodies' rows
read back as axial pieces, and the case list and its protocol rows."""
import unittest

import mpmath as mp

import generate_one_sphere_incidence_boolean_fixtures as fixtures


def piece(name):
    return fixtures.body_input(fixtures.BODIES[name])


class OneSphereIncidenceReferenceTests(unittest.TestCase):
    def setUp(self):
        self.dps = mp.mp.dps
        mp.mp.dps = 40

    def tearDown(self):
        mp.mp.dps = self.dps

    def test_closed_forms(self):
        # The cells' antiderivatives against quadrature of their integrands.
        r = mp.mpf(5)
        for a, b in ((-5, 5), (0, 5), (mp.mpf(5)/2, 5), (-3, 1)):
            a, b = mp.mpf(a), mp.mpf(b)
            self.assertLess(abs(fixtures.G(r, b)-fixtures.G(r, a)-mp.quad(lambda h: r*r-h*h, [a, b])),
                            mp.mpf(10)**-30)
            self.assertLess(abs(fixtures.K(r, b)-fixtures.K(r, a)-mp.quad(lambda h: (r*r-h*h)**1.5/3, [a, b])),
                            mp.mpf(10)**-25)
            self.assertLess(abs(fixtures.L(r, b)-fixtures.L(r, a)-mp.quad(lambda h: mp.sqrt(r*r-h*h), [a, b])),
                            mp.mpf(10)**-25)

    def test_a_wedge_by_its_cells(self):
        # The wedge between 45 and 135 degrees above the centre: a quarter
        # of the hemisphere, its area a quarter of the hemisphere's sphere
        # and three quarter discs; its centre on the bisector.
        w = piece('incidence_wedge45')
        cells = fixtures.Cells([w], 0)
        V, moments, A, div = cells.measures()
        r = mp.mpf(5)
        self.assertLess(abs(V-mp.pi*r**3/6), mp.mpf(10)**-30)
        self.assertLess(abs(div-V), mp.mpf(10)**-30)
        self.assertLess(abs(A-(mp.pi*r*r/2+3*mp.pi*r*r/4)), mp.mpf(10)**-30)
        c = [m/V for m in moments]
        self.assertLess(abs(c[0]-5), mp.mpf(10)**-30)
        self.assertGreater(c[1], 5)
        # The hemisphere's centre at 3 r / 8 above its base.
        self.assertLess(abs(c[2]-4-3*r/8), mp.mpf(10)**-30)
        self.assertEqual(cells.components(), (1, False))

    def test_pieces_of_one_sphere_by_their_cells(self):
        hemi, wedge = piece('incidence_hemi'), piece('incidence_wedge45')
        common = fixtures.Cells([hemi, wedge], ('common', 0, 1)).measures()
        fuse = fixtures.Cells([hemi, wedge], ('fuse', 0, 1)).measures()
        alone = fixtures.Cells([hemi, wedge], 1).measures()
        whole = fixtures.Cells([hemi, wedge], 0).measures()
        for got, want in ((common, alone), (fuse, whole)):
            self.assertLess(abs(got[0]-want[0]), mp.mpf(10)**-30)
            self.assertLess(abs(got[2]-want[2]), mp.mpf(10)**-30)
        # The hemisphere less the wedge: three quarters of it, one solid.
        cut = fixtures.Cells([hemi, wedge], ('cut', 0, 1))
        self.assertLess(abs(cut.measures()[0]-3*whole[0]/4), mp.mpf(10)**-30)
        self.assertEqual(cut.components(), (1, False))

    def test_octants_about_one_axis_touch(self):
        octant, back = piece('incidence_octant'), piece('incidence_back')
        fuse = fixtures.Cells([octant, back], ('fuse', 0, 1))
        self.assertEqual(fuse.components(), (2, True))
        self.assertEqual(fixtures.Cells([octant, back], ('common', 0, 1)).measures()[0], 0)

    def test_bodies_read_as_axial_pieces(self):
        # The wedge of `so3`: its half-planes' outward normals at the angle
        # of (12, 5) less a quarter and half a turn more.
        w = piece('incidence_wedge23')
        phi = mp.atan2(5, 12)
        want = sorted(((phi+mp.pi) % (2*mp.pi), (phi-mp.pi/2) % (2*mp.pi)))
        got = sorted(p % (2*mp.pi) for p in w.phis)
        for a, b in zip(got, want):
            self.assertLess(abs(a-b), mp.mpf(10)**-35)
        self.assertEqual(w.lo, 0)
        self.assertEqual(piece('incidence_high23').lo, mp.mpf(5)/2)
        # The hemisphere and the half: no half-plane, one.
        self.assertEqual(piece('incidence_hemi').phis, [])
        self.assertEqual(len(piece('incidence_half').phis), 1)
        # A box whose near face crosses the ball off the axis is refused.
        bad = fixtures.Body('bad', 'z', fixtures.on_axis('z', 0.0, 1.0), fixtures.along('z', 1, 0), (8.0, 8.0, 8.0))
        with self.assertRaises(AssertionError):
            fixtures.body_input(bad)

    def test_cases(self):
        listed = fixtures.all_cases()
        fixtures.validate(listed)
        self.assertEqual(len(listed), 45)
        kinds = {k: sum(1 for c in listed if c.kind == k) for k in ('solid', 'empty', 'degenerate')}
        self.assertEqual(kinds, {'solid': 29, 'empty': 4, 'degenerate': 12})
        for klass in fixtures.CLASSES:
            self.assertTrue(any(c.klass == klass for c in listed), klass)
        for name in ('hemi_wedge_fuse', 'cap_wedge_cut', 'chain_zone_cut'):
            c = next(c for c in listed if c.name == name)
            want = sum(s[0] == 'imported' for s in c.items)
            self.assertEqual(c.encode().count('\nbrep imported/'), want, name)
            self.assertEqual(c.native().count('\nbrep imported/'), want, name)
        chain = next(c for c in listed if c.name == 'chain_zone_cut')
        self.assertEqual(chain.encode().count('\nthen '), 1)
        self.assertEqual(chain.expr(), ('cut', ('fuse', 0, 1), 2))
        # Every body written once, its box's row, the hemisphere divided.
        text = fixtures.bodies_text()
        self.assertEqual(text.count('write '), len(fixtures.BODIES))
        self.assertEqual(text.count('\nbox '), len(fixtures.BODIES))
        self.assertEqual(text.count('\ndivide\n'), 1)
        used = {s[1] for c in listed for s in c.items if s[0] == 'imported'}
        self.assertEqual(used, set(fixtures.BODIES))
        # None of the earlier steps' bodies' names.
        import generate_imported_boolean_fixtures as a
        import generate_imported_arcs_boolean_fixtures as b
        import generate_imported_polyhedra_boolean_fixtures as c
        import generate_imported_pieces_boolean_fixtures as d
        import generate_one_sphere_boolean_fixtures as e
        for earlier in (a.BODIES, b.BODIES, c.BODIES, d.BODIES, e.BODIES):
            self.assertFalse(set(fixtures.BODIES) & set(earlier))


if __name__ == '__main__':
    unittest.main()
