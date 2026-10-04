"""Independent checks of S9e.4b.2's reference
(`imported_polyhedra_boolean_reference.py`) and fixtures
(`generate_imported_polyhedra_boolean_fixtures.py`): convex hulls against
known solids, the wedge's corners and prismatoid volume against direct
counts, Booleans of cubes and of an L-shaped union known in closed form,
the margins on contacts made on purpose, the prism test, and the case list
and its protocol rows."""
from fractions import Fraction as F
import unittest

import mpmath as mp

import imported_polyhedra_boolean_reference as ref
import generate_imported_polyhedra_boolean_fixtures as fixtures
import polyhedral_reference as pr


def cube(x, y, z, s):
    return [(x+i*s, y+j*s, z+k*s) for i in (0, 1) for j in (0, 1) for k in (0, 1)]


class ImportedPolyhedraReferenceTests(unittest.TestCase):
    def setUp(self):
        self.dps = mp.mp.dps
        mp.mp.dps = 40

    def tearDown(self):
        mp.mp.dps = self.dps

    def test_hulls(self):
        # A cube's hull: six square faces, its volume and centre exact; a
        # tetrahedron's four triangles, its volume the determinant's sixth;
        # interior points ignored.
        faces = ref.hull(cube(1, 2, 3, 2)+[(2, 3, 4)])
        self.assertEqual(len(faces), 6)
        self.assertTrue(all(len(f) == 4 for f in faces))
        v, m = pr.volume(faces)
        self.assertEqual(v, 8)
        self.assertEqual(tuple(x/v for x in m), (2, 3, 4))
        tet = ref.hull([(0, 0, 0), (3, 0, 0), (0, 4, 0), (0, 0, 5)])
        self.assertEqual(len(tet), 4)
        self.assertEqual(pr.volume(tet)[0], 10)

    def test_wedge_corners_and_prismatoid(self):
        # In the XY frame (exact axes: x, y = n x x, n) a wedge's corners are
        # its local coordinates: a pyramid of base 8 by 4 and height 6 with
        # its apex at (4, 6, 2). Its volume the hull's, the prismatoid
        # formula's and `A h / 3`; a frustum's `h (A0 + A1 + sqrt(A0 A1)) /
        # 3`.
        frame = fixtures.at('XY', (0, 0, 0))
        pts = ref.wedge_points(frame, 8.0, 6.0, 4.0, 4.0, 2.0, 4.0, 2.0)
        self.assertIn((F(4), F(6), F(2)), pts)
        self.assertIn((F(8), F(0), F(4)), pts)
        hull = pr.volume(ref.hull(pts))[0]
        self.assertEqual(hull, ref.wedge_closed(frame, 8.0, 6.0, 4.0, 4.0, 2.0, 4.0, 2.0))
        self.assertEqual(hull, F(8*4*6, 3))
        frustum = ref.wedge_points(frame, 8.0, 3.0, 8.0, 2.0, 2.0, 6.0, 6.0)
        self.assertEqual(pr.volume(ref.hull(frustum))[0], F(3, 3)*(64+16+32))

    def test_cubes(self):
        # Two cubes of side 2 overlapping in a 1 x 2 x 2 slab: fuse 12, cut
        # 4, common 4; areas 32, 16, 16; one solid each. Corner to corner
        # they meet at a point (two solids, the common empty).
        a, b = ref.Body([ref.hull(cube(0, 0, 0, 2))]), ref.Body([ref.hull(cube(1, 0, 0, 2))])
        for op, (n, v, area) in {'fuse': (1, 12, 32), 'cut': (1, 4, 16), 'common': (1, 4, 16)}.items():
            got = ref.result(a, b, op)
            self.assertEqual((got[0], got[1]), (n, v), op)
            self.assertLess(abs(got[2]-area), mp.mpf(10)**-35, op)
            self.assertEqual(got[5]['volume_three_ways'], 0, op)
            self.assertLess(got[5]['area_two_ways'], mp.mpf(10)**-35, op)
        c = ref.Body([ref.hull(cube(2, 2, 2, 2))])
        self.assertEqual(ref.result(a, c, 'fuse')[0], 2)
        self.assertEqual(ref.result(a, c, 'common')[0], 0)

    def test_union_boundary(self):
        # Two cells sharing a face (a cube of side 2 and a 1 by 2 by 2 box
        # beside it): the union's boundary leaves the shared square out (the
        # 3 by 2 by 2 box's area, 32, its 8 corners), its divergence volume
        # its cells'.
        body = ref.Body([ref.hull(cube(0, 0, 0, 2)), ref.hull([(2, 0, 0), (3, 0, 0), (2, 2, 0), (3, 2, 0),
                                                                (2, 0, 2), (3, 0, 2), (2, 2, 2), (3, 2, 2)])])
        area = sum((pr.area(f) for f in body.faces), mp.mpf(0))
        self.assertLess(abs(area-32), mp.mpf(10)**-35)
        self.assertEqual(ref.divergence_volume(body.faces), 12)
        self.assertEqual(len(body.vertices()), 8)

    def test_margins(self):
        # A cube on another's top, shifted by half: vertices on the other's
        # face counted as contacts; a unit cube a quarter above, the least
        # clearance a quarter over the size, 13/4.
        a = ref.Body([ref.hull(cube(0, 0, 0, 2))])
        b = ref.Body([ref.hull(cube(1, 1, 2, 2))])
        clearance, contacts, _ = ref.margins(a, b)
        self.assertGreater(contacts, 0)
        c = ref.Body([ref.hull(cube(F(1, 2), F(1, 2), F(9, 4), 1))])
        clearance, contacts, _ = ref.margins(a, c)
        self.assertEqual(contacts, 0)
        self.assertLess(abs(clearance-mp.mpf(1)/13), mp.mpf(10)**-35)

    def test_prism_test(self):
        box = ref.Body([ref.hull(cube(0, 0, 0, 2))])
        self.assertTrue(fixtures.is_prism(box))
        self.assertFalse(fixtures.is_prism(fixtures.BODIES['pyramid'].body()))
        self.assertTrue(fixtures.is_prism(fixtures.BODIES['ridge'].body()))
        self.assertFalse(fixtures.is_prism(fixtures.BODIES['ell'].body()))

    def test_cases_and_rows(self):
        cases = fixtures.all_cases()
        fixtures.validate(cases)
        self.assertEqual(len(cases), 48)
        steps = next(c for c in cases if c.name == 'steps_fuse')
        text = steps.encode()
        self.assertIn('brep imported/steps_low.brep', text)
        self.assertIn('brep imported/steps_high.brep', text)
        chain = next(c for c in cases if c.name == 'chain_vanes_cut')
        self.assertIn('then cut 95 swapped', chain.encode())
        self.assertIn('then cut swapped', chain.native())
        flush = [c.kind for c in cases if c.group == 'pyramid_flush']
        self.assertEqual(flush, ['degenerate', 'solid', 'empty'])


if __name__ == '__main__':
    unittest.main()
