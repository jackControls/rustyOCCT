"""Independent checks of S9e.4b.3c.1's fixtures
(`generate_one_sphere_boolean_fixtures.py`): the cap's closed form against a
direct integration of its slices, the chained reference on the cap and on
two pieces of one sphere against their closed forms (the hemisphere common
the cap is the cap, their fuse the hemisphere), the frames' normals,
and the case list and its protocol rows."""
import unittest

import mpmath as mp

import chained_curved_boolean_reference as ref
import generate_one_sphere_boolean_fixtures as fixtures


class OneSphereReferenceTests(unittest.TestCase):
    def setUp(self):
        self.dps = mp.mp.dps
        mp.mp.dps = 40

    def tearDown(self):
        mp.mp.dps = self.dps

    def test_cap_closed_form(self):
        # The cap above height h: its volume the integral of the slices'
        # discs, its area the zone and the disc (the hemisphere's at h = 0).
        for r, h in ((5, 0), (5, 2.5), (3, -1)):
            R, H = mp.mpf(r), mp.mpf(h)
            V, A = fixtures.cap_closed(r, h)
            self.assertLess(abs(V-mp.quad(lambda w: mp.pi*(R**2-w**2), [H, R])), mp.mpf(10)**-30)
            self.assertLess(abs(A-(2*mp.pi*R*(R-H)+mp.pi*(R**2-H**2))), mp.mpf(10)**-30)
        V, A = fixtures.cap_closed(5, 0)
        self.assertLess(abs(V-2*mp.pi*125/3), mp.mpf(10)**-30)
        self.assertLess(abs(A-3*mp.pi*25), mp.mpf(10)**-30)

    def test_frames(self):
        # No frame's normal lies in a split plane of the kernel's whole
        # sphere's rotations (rows (2, 3, 6), (3, -6, 2); (1, 4, 8), (4, 7,
        # -4); (2, 6, 9), (6, 7, -6); (1, 2, 2), (2, 1, -2)).
        rows = [((2, 3, 6), (3, -6, 2)), ((1, 4, 8), (4, 7, -4)), ((2, 6, 9), (6, 7, -6)), ((1, 2, 2), (2, 1, -2))]
        cross = lambda a, b: (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])
        for name in ('SKEW', 'SKEW2', 'SKEW4'):
            n = fixtures.FRAMES[name][:3]
            for a, b in rows:
                self.assertNotEqual(sum(p*q for p, q in zip(n, cross(a, b))), 0, name)

    def test_pieces_of_one_sphere_on_the_chained_reference(self):
        # The hemisphere and the cap of one sphere and frame: their common
        # the cap, their fuse the hemisphere, by the chained reference with
        # both pieces' spheres one surface.
        scan, ref.SCAN = ref.SCAN, 48
        try:
            hemi, cap = fixtures.BODIES['sphere_hemi'], fixtures.BODIES['sphere_cap']
            chain = ref.Chain([fixtures.construction(s, 0) for s in hemi.specs()+cap.specs()])
            common = chain.measures(('common', ('common', 0, 1), ('common', 2, 3)))
            fuse = chain.measures(('fuse', ('common', 0, 1), ('common', 2, 3)))
        finally:
            ref.SCAN = scan
        Vc, Ac = cap.closed()
        Vh, Ah = hemi.closed()
        self.assertLess(abs(common[0]-Vc), mp.mpf(10)**-25)
        self.assertLess(abs(common[2]-Ac), mp.mpf(10)**-25)
        self.assertLess(abs(fuse[0]-Vh), mp.mpf(10)**-25)
        self.assertLess(abs(fuse[2]-Ah), mp.mpf(10)**-25)

    def test_cases(self):
        listed = fixtures.all_cases()
        fixtures.validate(listed)
        self.assertEqual(len(listed), 22)
        kinds = {k: sum(1 for c in listed if c.kind == k) for k in ('solid', 'empty', 'degenerate', 'unsupported')}
        self.assertEqual(kinds, {'solid': 19, 'empty': 1, 'degenerate': 2, 'unsupported': 0})
        for klass in fixtures.CLASSES:
            self.assertTrue(any(c.klass == klass for c in listed), klass)
        for name in ('hemi_cap_fuse', 'cap_octant_cut', 'ball_cap_common', 'chain_cap_cut'):
            c = next(c for c in listed if c.name == name)
            want = sum(s[0] == 'imported' for s in c.items)
            self.assertEqual(c.encode().count('\nbrep imported/'), want, name)
            self.assertEqual(c.native().count('\nbrep imported/'), want, name)
        chain = next(c for c in listed if c.name == 'chain_cap_cut')
        self.assertEqual(chain.encode().count('\nthen '), 1)
        self.assertEqual(chain.expr(), ('cut', ('cut', ('common', 0, 1), ('common', 2, 3)), 4))
        # Every body written once, the divided ones with their `divide` row.
        text = fixtures.bodies_text()
        self.assertEqual(text.count('write '), len(fixtures.BODIES))
        self.assertEqual(text.count('\ndivide\n'), sum(b.divide for b in fixtures.BODIES.values()))
        self.assertEqual(sum(b.divide for b in fixtures.BODIES.values()), 2)
        used = {s[1] for c in listed for s in c.items if s[0] == 'imported'}
        self.assertEqual(used, set(fixtures.BODIES))
        # None of the earlier steps' bodies' names.
        import generate_imported_boolean_fixtures as a
        import generate_imported_arcs_boolean_fixtures as b
        import generate_imported_polyhedra_boolean_fixtures as c
        import generate_imported_pieces_boolean_fixtures as d
        for earlier in (a.BODIES, b.BODIES, c.BODIES, d.BODIES):
            self.assertFalse(set(fixtures.BODIES) & set(earlier))


if __name__ == '__main__':
    unittest.main()
