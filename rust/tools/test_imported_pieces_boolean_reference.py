"""Independent checks of S9e.4b.3a's fixtures
(`generate_imported_pieces_boolean_fixtures.py`): the pieces' closed forms
against a direct quadrature of their sections, the chained reference on the
octant against its closed form, the test that a body is no S9e.4a
construction, and the case list and its protocol rows."""
import unittest

import mpmath as mp

import chained_curved_boolean_reference as ref
import generate_imported_pieces_boolean_fixtures as fixtures


class ImportedPiecesReferenceTests(unittest.TestCase):
    def setUp(self):
        self.dps = mp.mp.dps
        mp.mp.dps = 40

    def tearDown(self):
        mp.mp.dps = self.dps

    def test_closed_forms(self):
        # On exact axes the corner is an eighth of the ball, the lune a
        # quarter and the half a half; on turned axes the corner's solid
        # angle by Van Oosterom and Strackee against the sum of its faces'
        # spherical excess (`Omega = A + B + C - pi`, its dihedral angles
        # between the planes of two axes each).
        R = mp.mpf(4)
        xy = fixtures.at('XY', (0, 0, 0))
        V, A = fixtures.corner_closed(4, xy)
        self.assertLess(abs(V-mp.pi*R**3/6), mp.mpf(10)**-35)
        self.assertLess(abs(A-(mp.pi*R**2/2+3*mp.pi*R**2/4)), mp.mpf(10)**-35)
        V, A = fixtures.lune_closed(4, xy)
        self.assertLess(abs(V-mp.pi*R**3/3), mp.mpf(10)**-35)
        self.assertLess(abs(A-2*mp.pi*R**2), mp.mpf(10)**-35)
        for name in ('SKEW', 'SKEW2', 'SKEW3'):
            frame = fixtures.at(name, (0, 0, 0))
            x, y, n = fixtures.axes(frame)
            cross = lambda a, b: [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
            unit = lambda v: [c/mp.sqrt(sum(d*d for d in v)) for c in v]
            # Each edge's dihedral angle: between its two faces' inward
            # normals, pi less their angle.
            faces = [unit(cross(y, n)), unit(cross(n, x)), unit(cross(x, y))]
            inward = [f if sum(p*q for p, q in zip(f, a)) > 0 else [-c for c in f] for f, a in zip(faces, (x, y, n))]
            dihedral = [mp.pi-fixtures.angle(inward[i], inward[j]) for i, j in ((0, 1), (1, 2), (2, 0))]
            omega = sum(dihedral)-mp.pi
            V, _ = fixtures.corner_closed(4, frame)
            self.assertLess(abs(V-omega*R**3/3), mp.mpf(10)**-30, name)

    def test_octant_on_the_chained_reference(self):
        # The octant (the ball's corner through its centre on SKEW's axes):
        # the chained reference's measures of `P common B` its closed form.
        scan, ref.SCAN = ref.SCAN, 48
        try:
            b = fixtures.BODIES['octant']
            chain = ref.Chain([fixtures.construction(s, 0) for s in b.specs()])
            V, _, A = chain.measures(('common', 0, 1))
        finally:
            ref.SCAN = scan
        V0, A0 = b.closed()
        self.assertLess(abs(V-V0), mp.mpf(10)**-28)
        self.assertLess(abs(A-A0), mp.mpf(10)**-28)

    def test_bodies_are_this_steps(self):
        # Every piece has a plane of its box not normal to its sphere's
        # axis; a cap's box only planes normal or along it is S9e.4a's.
        for b in fixtures.BODIES.values():
            self.assertEqual(fixtures.oblique(b) or b.earlier, True, b.name)
        # A box whose planes are all normal to a sphere's axis or along it
        # leaves no plane normal only; a cylinder's needs one oblique.
        cyl = fixtures.Body('c', 'cylinder', fixtures.prism([fixtures.disc(0.0, 0.0, 1.0)], fixtures.at('XY', (0, 0, 0)),
                                                            0.0, 1.0),
                            fixtures.box(fixtures.at('XY', (-5, -5, 1)), (10.0, 10.0, 10.0)))
        self.assertFalse(fixtures.oblique(cyl))
        tilted = fixtures.Body('t', 'cylinder', cyl.primitive, fixtures.box(fixtures.at('TILT', (-5, -5, 1)),
                                                                            (10.0, 10.0, 10.0)))
        self.assertTrue(fixtures.oblique(tilted))
        # The half ball's sphere seam (its frame's x) points away from it.
        half = fixtures.BODIES['half']
        self.assertEqual(half.primitive[2][6:9], tuple(-c for c in fixtures.FRAMES['SKEW2'][3:]))

    def test_cases(self):
        listed = fixtures.all_cases()
        fixtures.validate(listed)
        self.assertEqual(len(listed), 45)
        kinds = {k: sum(1 for c in listed if c.kind == k) for k in ('solid', 'empty', 'degenerate', 'unsupported')}
        self.assertEqual(kinds, {'solid': 33, 'empty': 0, 'degenerate': 6, 'unsupported': 6})
        for klass in fixtures.CLASSES:
            self.assertTrue(any(c.klass == klass for c in listed), klass)
        for name in ('octant_box_fuse', 'box_octant_cut', 'pieces_common', 'half_box_fuse', 'chain_octant_cut'):
            c = next(c for c in listed if c.name == name)
            enc, nat = c.encode(), c.native()
            want = sum(s[0] in ('imported', 'base') for s in c.items)
            self.assertEqual(enc.count('\nbrep imported/'), want, name)
            self.assertEqual(nat.count('\nbrep imported/'), want, name)
        chain = next(c for c in listed if c.name == 'chain_octant_cut')
        self.assertEqual(chain.encode().count('\nthen '), 1)
        self.assertEqual(chain.expr(), ('cut', ('cut', ('common', 0, 1), 2), 3))
        both = next(c for c in listed if c.name == 'pieces_fuse')
        self.assertEqual(both.expr(), ('fuse', ('common', 0, 1), ('common', 2, 3)))
        # Every body written once: a sphere's row, a Boolean, a box's prism.
        text = fixtures.bodies_text()
        self.assertEqual(text.count('write '), len(fixtures.BODIES))
        self.assertEqual(text.count('\nsphere '), len(fixtures.BODIES))
        self.assertEqual(text.count('\nboolean '), len(fixtures.BODIES))
        used = {s[1] for c in listed for s in c.items if s[0] == 'imported'}
        self.assertEqual(used, set(fixtures.BODIES))
        # None of the earlier steps' bodies' names.
        import generate_imported_boolean_fixtures as a
        import generate_imported_arcs_boolean_fixtures as b
        import generate_imported_polyhedra_boolean_fixtures as c
        for earlier in (a.BODIES, b.BODIES, c.BODIES):
            self.assertFalse(set(fixtures.BODIES) & set(earlier))


if __name__ == '__main__':
    unittest.main()
