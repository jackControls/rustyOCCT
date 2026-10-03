"""Independent checks of S9e.3a's reference
(`chained_curved_boolean_reference.py`) on chains whose results are known
without it: boxes (volumes, areas and centres by their grid cells, exactly),
with faces of several boxes on one plane of both orientations and a deeper
chain; a hemisphere less a coaxial cylinder (closed forms); a chain of three
prisms against S9e.1's slicing reference (`chained_boolean_reference.py`,
another engine); and of the fixture list
(`generate_given_curved_boolean_fixtures.py`)."""
from fractions import Fraction as F
import unittest

import mpmath as mp

import chained_curved_boolean_reference as ref
import chained_boolean_reference as sliced
import generate_given_curved_boolean_fixtures as fixtures
from generate_curved_boolean_fixtures import at, disc, square
from test_chained_boolean_reference import cells


def box(x0, x1, y0, y1, z0, z1, op):
    return fixtures.solid(fixtures.prism([square(float(x0), float(y0), float(x1), float(y1))], at('XY', (0, 0, 0)),
                                         z0, z1), op)


def side_box(y0, z0, y1, z1, x0, x1, op):
    """A box in the SIDE frame (u = y, v = z, w = x)."""
    return fixtures.solid(fixtures.prism([square(float(y0), float(z0), float(y1), float(z1))],
                                         at('SIDE', (0, 0, 0)), x0, x1), op)


def mpq(x):
    return mp.mpf(x.numerator)/x.denominator


class GivenCurvedReferenceTests(unittest.TestCase):
    def setUp(self):
        self.dps = mp.mp.dps
        mp.mp.dps = 40

    def tearDown(self):
        mp.mp.dps = self.dps

    def check_cells(self, boxes, solids, exprs):
        # Boxes' faces meet in few events: a coarser scan finds them all
        # (a missed one would be refined, or fail the cells' equality).
        scan, ref.SCAN = ref.SCAN, 48
        try:
            self._check_cells(boxes, solids, exprs)
        finally:
            ref.SCAN = scan

    def _check_cells(self, boxes, solids, exprs):
        chain = ref.Chain(solids)
        grid = [tuple(F(x) for x in b) for b in boxes]
        for expr in exprs:
            vol, area, centre = cells(grid, lambda *m: ref.evaluate(expr, list(m)))
            v, mom, a = chain.measures(expr)
            self.assertLess(abs(v-mpq(vol)), mp.mpf(10)**-30, expr)
            self.assertLess(abs(a-mpq(area)), mp.mpf(10)**-30, expr)
            if vol:
                for x, y in zip(mom, centre):
                    self.assertLess(abs(x/v-mpq(y)), mp.mpf(10)**-30, expr)
            v1, _, a1 = chain.measures(expr, 1)
            self.assertLess(abs(v1-v)+abs(a1-a), mp.mpf(10)**-30, expr)

    def test_boxes_by_their_cells(self):
        # The first Boolean in frames of different axes (XY and SIDE); the
        # third crossing both, sharing none of their planes.
        a = (0, 4, 0, 3, 0, 2)
        b = (1, 3, -1, 2, 1, 3)
        c = (2, 5, 1, 4, -1, 1.5)
        solids = [box(*a, 91), side_box(b[2], b[4], b[3], b[5], b[0], b[1], 92), box(*c, 94)]
        exprs = [ref.chain_expr(op1, op2, sw) for op1 in ('fuse', 'cut') for op2 in ('fuse', 'cut', 'common')
                 for sw in (False, True)]
        self.check_cells([a, b, c], solids, exprs)

    def test_boxes_sharing_planes(self):
        # Faces on one plane in both orientations (a face of several boxes
        # counted once), and a deeper chain with a fourth box.
        a = (0, 4, 0, 3, 0, 2)
        b = (0, 2, 0, 3, 1, 3)
        c = (2, 4, 0, 3, 0, 2)
        d = (1, 3, 1, 2, 2, 4)
        solids = [box(*x, 91+k) for k, x in enumerate((a, b, c, d))]
        exprs = [('cut', ('fuse', 0, 1), 2), ('fuse', ('cut', 0, 1), 2), ('common', 3, ('fuse', ('cut', 0, 1), 2)),
                 ('fuse', ('fuse', ('fuse', 0, 1), 2), 3)]
        self.check_cells([a, b, c, d], solids, exprs)

    def test_hemisphere_less_a_cylinder(self):
        # A sphere of radius 3 above its centre's plane (a box below cut
        # away), less a coaxial cylinder of radius 1: closed forms.
        R, r = mp.mpf(3), mp.mpf(1)
        solids = [fixtures.solid(fixtures.sphere(3.0, at('XY', (0, 0, 0))), 91),
                  box(-4, 4, -4, 4, -4, 0, 92),
                  fixtures.solid(fixtures.prism([disc(0.0, 0.0, 1.0)], at('XY', (0, 0, 0)), -5.0, 5.0), 94)]
        chain = ref.Chain(solids)
        h = mp.sqrt(R*R-r*r)
        # The cylinder's part in the hemisphere: 2 pi/3 (R^3 - h^3).
        inner = 2*mp.pi/3*(R**3-h**3)
        V = 2*mp.pi*R**3/3-inner
        # Its walls: the sphere's zone above h (2 pi R (R - h)) removed, the
        # cylinder's wall up to h and the base's disc of radius 1 removed.
        A = 3*mp.pi*R*R-2*mp.pi*R*(R-h)+2*mp.pi*r*h-mp.pi*r*r
        # The centre's height: the hemisphere's moment less the core's.
        mz_hemi = mp.pi*R**4/4
        mz_core = mp.quad(lambda z: z*mp.pi*min(r*r, R*R-z*z), [0, h, R])
        v, mom, a = chain.measures(('cut', ('cut', 0, 1), 2))
        self.assertLess(abs(v-V), mp.mpf(10)**-30)
        self.assertLess(abs(a-A), mp.mpf(10)**-30)
        self.assertLess(abs(mom[2]-(mz_hemi-mz_core)), mp.mpf(10)**-30)
        self.assertLess(abs(mom[0])+abs(mom[1]), mp.mpf(10)**-30)

    def test_against_the_slicing_reference(self):
        # S9e.1's quarter cylinder bored coaxially (`QUARTER`): three
        # prisms, by both engines.
        import generate_chained_boolean_fixtures as chained
        obj, tool = [chained.prism(*x, op=91+k) for k, x in enumerate(chained.QUARTER)]
        third = chained.prism([disc(0.0, 0.0, 1.5)], at('XY', (0, 0, 0)), 1.0, 7.0, 94)
        mine = ref.Chain([obj, tool, third])
        for op2 in ('fuse', 'cut'):
            other = sliced.Chain(obj, tool, third, 'common', False).result(op2)
            v, mom, a = mine.measures(ref.chain_expr('common', op2, False))
            self.assertLess(abs(v-other[1]), mp.mpf(10)**-30)
            self.assertLess(abs(a-other[2]), mp.mpf(10)**-30)
            for i in range(3):
                self.assertLess(abs(mom[i]/v-other[3][i]), mp.mpf(10)**-30)

    def test_cases(self):
        listed = fixtures.cases()
        fixtures.validate(listed)
        self.assertEqual(len(listed), 48)
        self.assertEqual(sum(1 for c in listed if c.kind == 'degenerate'), 6)
        self.assertTrue(any(c.stages[-1][1] for c in listed))
        self.assertTrue(any(len(c.stages) > 1 for c in listed))
        for klass in fixtures.CLASSES:
            self.assertTrue(any(c.klass == klass for c in listed), klass)
        # The protocol's further then rows.
        deep = next(c for c in listed if len(c.stages) > 1)
        self.assertEqual(deep.encode().count('\nthen '), len(deep.stages))


if __name__ == '__main__':
    unittest.main()
