"""Independent checks of S9e.4a's fixtures (`generate_imported_boolean_fixtures.py`):
the references on results known without them (two boxes by their grid
cells, exactly; the coaxial sections against the chained reference; the
hemispheres' closed forms against a direct quadrature of their sections),
the reader of the bodies' files on a file of known content, and the case
list and its protocol rows (an imported input's `brep` row)."""
from fractions import Fraction as F
import unittest

import mpmath as mp

import chained_curved_boolean_reference as ref
import generate_imported_boolean_fixtures as fixtures
from generate_curved_boolean_fixtures import square
from test_chained_boolean_reference import cells


def mpq(x):
    return mp.mpf(x.numerator)/x.denominator


class ImportedReferenceTests(unittest.TestCase):
    def setUp(self):
        self.dps = mp.mp.dps
        mp.mp.dps = 40

    def tearDown(self):
        mp.mp.dps = self.dps

    def test_boxes_by_their_cells(self):
        # The imported box's construction (`MakeBox` of [0, 10]^2 x [0, 4])
        # against a box overlapping it: the chained reference's measures are
        # the grid cells' exactly.
        scan, ref.SCAN = ref.SCAN, 48
        try:
            a = (0, 10, 0, 10, 0, 4)
            b = (3, 12, -1, 6, 1, 7)
            box = lambda x, op: fixtures.construction(
                fixtures.prism([square(float(x[0]), float(x[2]), float(x[1]), float(x[3]))],
                               fixtures.at('XY', (0, 0, 0)), x[4], x[5]), op)
            chain = ref.Chain([box(a, 91), box(b, 92)])
            grid = [tuple(F(x) for x in a), tuple(F(x) for x in b)]
            for op in fixtures.OPS:
                vol, area, centre = cells(grid, lambda *m: ref.evaluate((op, 0, 1), list(m)))
                v, mom, s = chain.measures((op, 0, 1))
                self.assertLess(abs(v-mpq(vol)), mp.mpf(10)**-30, op)
                self.assertLess(abs(s-mpq(area)), mp.mpf(10)**-30, op)
                for x, y in zip(mom, centre):
                    self.assertLess(abs(x/v-mpq(y)), mp.mpf(10)**-30, op)
        finally:
            ref.SCAN = scan

    def test_coaxial_sections_against_the_chained_reference(self):
        # The imported cylinder (radius 3 about (5, 5), z in [-1, 5]) and the
        # box [0, 10]^2 x [0, 4]: the coaxial sections' closed forms are the
        # chained reference's sweeps'.
        scan, ref.SCAN = ref.SCAN, 64
        try:
            cases = [fixtures.construction(fixtures.resolve(fixtures.imported('cyl')), 91),
                     fixtures.construction(fixtures.B, 92)]
            chain = ref.Chain(cases)
            forms = fixtures.coaxial_forms()
            for op in fixtures.OPS:
                V, A, c = forms[('cyl_box', op)]
                v, mom, a = chain.measures((op, 0, 1))
                self.assertLess(abs(v-V), mp.mpf(10)**-30, op)
                self.assertLess(abs(a-A), mp.mpf(10)**-30, op)
                for i in range(3):
                    self.assertLess(abs(mom[i]/v-c[i]), mp.mpf(10)**-30, op)
        finally:
            ref.SCAN = scan

    def test_hemispheres_by_their_sections(self):
        # The ball of radius 3 about (5, 5, 4) against [0, 10]^2 x [0, 4]:
        # volumes and heights of the centres by quadrature of the horizontal
        # sections, areas by the zones' (2 pi R dz) and the faces'.
        forms = fixtures.hemisphere_forms()
        R = mp.mpf(3)
        sec = lambda z: mp.pi*(R*R-(z-4)**2)
        lower = mp.quad(sec, [1, 4])
        lower_z = mp.quad(lambda z: z*sec(z), [1, 4])/lower
        upper_z = mp.quad(lambda z: z*sec(z), [4, 7])/lower
        V, A, c = forms[('ball_box', 'common')]
        self.assertLess(abs(V-lower), mp.mpf(10)**-35)
        self.assertLess(abs(c[2]-lower_z), mp.mpf(10)**-35)
        self.assertLess(abs(A-(mp.quad(lambda z: 2*mp.pi*R, [1, 4])+mp.pi*R*R)), mp.mpf(10)**-35)
        V, A, c = forms[('ball_box', 'cut')]
        self.assertLess(abs(c[2]-upper_z), mp.mpf(10)**-35)
        V, A, c = forms[('ball_box', 'fuse')]
        self.assertLess(abs(V-(400+lower)), mp.mpf(10)**-35)
        self.assertLess(abs(c[2]-(400*2+lower*upper_z)/(400+lower)), mp.mpf(10)**-35)
        self.assertLess(abs(A-(360-mp.pi*R*R+2*mp.pi*R*R)), mp.mpf(10)**-35)

    def test_the_files_reader(self):
        # The hemisphere OCCT wrote (`MakeSphere` of radius 4, latitudes 0 to
        # pi / 2): a spherical face and its equator's disc, two vertices (the
        # pole and the equator's seam vertex, OCCT's 15 digits off the axis
        # and the plane by 2e-16), both on the construction's surfaces within
        # 1e-15 and off every surface where moved by 1e-9.
        faces, vertices = fixtures.stored_records((fixtures.ROOT/'fixtures/imported/dome.brep').read_text())
        self.assertEqual(sorted(faces), ['plane', 'sphere'])
        self.assertEqual(len(vertices), 2)
        want, surfs = fixtures.construction_surfaces(fixtures.BODIES['dome'].spec)
        self.assertEqual(sorted(want), sorted(faces))
        for v in vertices:
            X = tuple(mp.mpf(c) for c in v)
            self.assertTrue(any(abs(s.f(X))/fixtures.tc.norm(s.grad(X)) < mp.mpf(10)**-15 for s in surfs), v)
            off = (X[0]+mp.mpf(10)**-9, X[1], X[2]+mp.mpf(10)**-9)
            self.assertTrue(all(abs(s.f(off))/fixtures.tc.norm(s.grad(off)) > mp.mpf(10)**-11 for s in surfs), v)
        # A prism's top shares its bottom's records under a translation: the
        # tangent profile's eight vertices from four records.
        faces, vertices = fixtures.stored_records((fixtures.ROOT/'fixtures/imported/dee.brep').read_text())
        self.assertEqual(len(faces), 6)
        self.assertEqual(sorted({v[2] for v in vertices}), [0.0, 6.0])
        self.assertEqual(len(set(vertices)), 8)

    def test_cases(self):
        listed = fixtures.all_cases()
        fixtures.validate(listed)
        self.assertEqual(len(listed), 69)
        kinds = {k: sum(1 for c in listed if c.kind == k) for k in ('solid', 'empty', 'degenerate', 'unsupported')}
        self.assertEqual(kinds, {'solid': 57, 'empty': 0, 'degenerate': 9, 'unsupported': 3})
        for klass in fixtures.CLASSES:
            self.assertTrue(any(c.klass == klass for c in listed), klass)
        # An imported input's one row in both protocols, as object, as tool
        # and both.
        for name in ('box_slab_fuse', 'box_cyl_cut', 'both_common', 'chain_drill_cut'):
            c = next(c for c in listed if c.name == name)
            enc, nat = c.encode(), c.native()
            want = sum(s[0] == 'imported' for s in c.specs)
            self.assertEqual(enc.count('\nbrep imported/'), want, name)
            self.assertEqual(nat.count('\nbrep imported/'), want, name)
        chain = next(c for c in listed if c.name == 'chain_drill_cut')
        self.assertEqual(chain.encode().count('\nthen '), 1)
        # Every body written once, by its rows.
        text = fixtures.bodies_text()
        self.assertEqual(text.count('write '), len(fixtures.BODIES))
        self.assertIn('\nbox 0.0 0.0 0.0 0.0 0.0 1.0 1.0 0.0 0.0 10.0 10.0 4.0\n', text)
        self.assertIn('\ncylinder 5.0 5.0 -1.0 0.0 0.0 1.0 1.0 0.0 0.0 3.0 6.0\n', text)
        used = {s[1] for c in listed for s in c.specs if s[0] == 'imported'}
        self.assertEqual(used, set(fixtures.BODIES))


if __name__ == '__main__':
    unittest.main()
