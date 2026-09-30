"""Independent checks of S9e.2's extension of the chained reference
(`chained_boolean_reference.py`: a selector picking one solid of the first
result) on chains whose results are known without it: four boxes (a box
severed by a slab, the selector holding one part, a third crossing it) by
the grid of their cells, the selector's separation measure, a selector
holding the whole first result against the three-prism chain, and the
fixture list (`generate_given_boolean_fixtures.py`)."""
import unittest

import mpmath as mp

import chained_boolean_reference as ref
import generate_given_boolean_fixtures as fixtures
from generate_curved_boolean_fixtures import at, square
from test_chained_boolean_reference import box, cells, side_box


class GivenBooleanReferenceTests(unittest.TestCase):
    def setUp(self):
        self.dps = mp.mp.dps
        mp.mp.dps = 40

    def tearDown(self):
        mp.mp.dps = self.dps

    def chain(self, a, b, c, d, op1, swapped):
        obj = box(*a, op=91)
        tool = side_box(b[2], b[4], b[3], b[5], b[0], b[1], op=92)
        third = box(*c, op=94)
        selector = box(*d, op=96)
        return ref.Chain(obj, tool, third, op1, swapped, selector=selector)

    def check_boxes(self, a, b, c, d, op1, swapped):
        """The four boxes' chain against the grid of their cells."""
        from fractions import Fraction as F
        chain = self.chain(a, b, c, d, op1, swapped)
        grid = [tuple(F(x) for x in bx) for bx in (a, b, c, d)]
        for op2 in ref.OPS:
            vol, area, centre = cells(grid, ref.chained(op1, op2, swapped, True))
            n, v, s, cen = chain.result(op2)
            if vol == 0:
                self.assertEqual(n, 0)
                continue
            self.assertLess(abs(v-mp.mpf(vol.numerator)/vol.denominator), mp.mpf(10)**-30)
            self.assertLess(abs(s-mp.mpf(area.numerator)/area.denominator), mp.mpf(10)**-30)
            for x, y in zip(cen, centre):
                self.assertLess(abs(x-mp.mpf(y.numerator)/y.denominator), mp.mpf(10)**-30)
        return chain

    def test_four_boxes_by_their_cells(self):
        # A box severed by a slab of the SIDE frame (two solids, x < 2 and
        # x > 3), the selector holding the first, a third crossing both.
        a = (0, 6, 0, 3, 0, 2)
        b = (2, 3, -1, 4, -1, 3)
        c = (1, 5, 1, 4, 1, 3)
        d = (-1, 2.5, -1, 4, -1, 3)
        for op1 in ('cut',):
            for swapped in (False, True):
                chain = self.check_boxes(a, b, c, d, op1, swapped)
                self.assertEqual(ref.separation(chain), 0)

    def test_a_selector_across_the_result_is_measured(self):
        # The selector's face x = 1 crosses the first solid: its area there
        # (3 x 2) is what separation reports.
        chain = self.chain((0, 6, 0, 3, 0, 2), (2, 3, -1, 4, -1, 3), (1, 5, 1, 4, 1, 3),
                           (-1, 1, -1, 4, -1, 3), 'cut', False)
        self.assertLess(abs(ref.separation(chain)-6), mp.mpf(10)**-30)

    def test_a_selector_holding_everything_is_the_chain(self):
        # A selector around the whole first result leaves S9e.1's chain.
        obj, tool = [fixtures.prism(*x, op=91+k) for k, x in enumerate(fixtures.SEVERED)]
        third = fixtures.prism([square(3.0, -1.0, 7.0, 4.0)], at('XY', (0, 0, 0)), 1.0, 7.0, 94)
        whole = fixtures.prism([square(-2.0, -2.0, 12.0, 12.0)], at('XY', (0, 0, 0)), -2.0, 7.0, 96)
        a = ref.Chain(obj, tool, third, 'cut', True, selector=whole)
        b = ref.Chain(obj, tool, third, 'cut', True)
        self.assertEqual(ref.separation(a), 0)
        for op2 in ref.OPS:
            x, y = a.result(op2), b.result(op2)
            self.assertEqual(x[0], y[0])
            self.assertLess(abs(x[1]-y[1])+abs(x[2]-y[2]), mp.mpf(10)**-30)

    def test_cases(self):
        listed = fixtures.cases()
        fixtures.validate(listed)
        self.assertEqual(len(listed), 36)
        self.assertEqual(sum(1 for c in listed if c.kind == 'degenerate'), 6)
        self.assertTrue(any(c.swapped for c in listed))
        for klass in fixtures.CLASSES:
            self.assertTrue(any(c.klass == klass for c in listed), klass)
        # At most two axis directions per chain, the selector's among them
        # (the reference's slicing).
        for c in listed:
            prisms = [c.obj, c.tool, c.third]+([c.selector] if c.selector is not None else [])
            ref.Slicing3([ref.Prism(x) for x in prisms])
        # DRAW's rollex: bcut_simple/L3 to L6 cut the stack by the cylinder.
        names = {c.name for c in listed}
        self.assertTrue({'rollex_turned_cut', 'rollex_flat_cut'} <= names)


if __name__ == '__main__':
    unittest.main()
