"""Independent checks of S9e.1's chained reference
(`chained_boolean_reference.py`) on chains whose results are known without
it: three boxes (volumes, areas and centres by inclusion and exclusion over
axis-aligned boxes), a third prism far from the first result, the common's
symmetry in its arguments, and of the fixture list
(`generate_chained_boolean_fixtures.py`)."""
import unittest

import mpmath as mp

import chained_boolean_reference as ref
import generate_chained_boolean_fixtures as fixtures
from generate_curved_boolean_fixtures import at, disc, square


def prism(boundaries, frame, start, end, op):
    return fixtures.prism(boundaries, frame, start, end, op)


def box(x0, x1, y0, y1, z0, z1, op):
    return prism([square(float(x0), float(y0), float(x1), float(y1))], at('XY', (0, 0, 0)), float(z0), float(z1), op)


def side_box(y0, z0, y1, z1, x0, x1, op):
    """A box in the SIDE frame (u = y, v = z, w = x)."""
    return prism([square(float(y0), float(z0), float(y1), float(z1))], at('SIDE', (0, 0, 0)), float(x0), float(x1), op)


def cells(boxes, fn):
    """Volume, area and moments of the region `fn(in each box)` over the grid
    of the boxes' coordinates (exact Fractions)."""
    from fractions import Fraction as F
    xs = sorted({b[k] for b in boxes for k in (0, 1)})
    ys = sorted({b[k] for b in boxes for k in (2, 3)})
    zs = sorted({b[k] for b in boxes for k in (4, 5)})

    def inside(b, p):
        return b[0] < p[0] < b[1] and b[2] < p[1] < b[3] and b[4] < p[2] < b[5]

    def held(i, j, k):
        if not (0 <= i < len(xs)-1 and 0 <= j < len(ys)-1 and 0 <= k < len(zs)-1):
            return False
        p = ((xs[i]+xs[i+1])/2, (ys[j]+ys[j+1])/2, (zs[k]+zs[k+1])/2)
        return fn(*[inside(b, p) for b in boxes])

    vol, mom, area = F(0), [F(0)]*3, F(0)
    for i in range(len(xs)-1):
        for j in range(len(ys)-1):
            for k in range(len(zs)-1):
                if not held(i, j, k):
                    continue
                dx, dy, dz = xs[i+1]-xs[i], ys[j+1]-ys[j], zs[k+1]-zs[k]
                v = dx*dy*dz
                vol += v
                c = ((xs[i]+xs[i+1])/2, (ys[j]+ys[j+1])/2, (zs[k]+zs[k+1])/2)
                mom = [m+v*x for m, x in zip(mom, c)]
                for (di, dj, dk), a in (((1, 0, 0), dy*dz), ((-1, 0, 0), dy*dz), ((0, 1, 0), dx*dz),
                                        ((0, -1, 0), dx*dz), ((0, 0, 1), dx*dy), ((0, 0, -1), dx*dy)):
                    if not held(i+di, j+dj, k+dk):
                        area += a
    return vol, area, [m/vol for m in mom] if vol else None


class ChainedBooleanReferenceTests(unittest.TestCase):
    def setUp(self):
        self.dps = mp.mp.dps
        mp.mp.dps = 40

    def tearDown(self):
        mp.mp.dps = self.dps

    def check_boxes(self, a, b, c, op1, swapped):
        """The three boxes' chain against the grid of their cells."""
        from fractions import Fraction as F
        obj = box(*a, op=91)
        tool = side_box(b[2], b[4], b[3], b[5], b[0], b[1], op=92)
        third = box(*c, op=94)
        chain = ref.Chain(obj, tool, third, op1, swapped)
        grid = [tuple(F(x) for x in bx) for bx in (a, b, c)]
        for op2 in ref.OPS:
            vol, area, centre = cells(grid, ref.chained(op1, op2, swapped))
            n, v, s, cen = chain.result(op2)
            if vol == 0:
                self.assertEqual(n, 0)
                continue
            self.assertLess(abs(v-mp.mpf(vol.numerator)/vol.denominator), mp.mpf(10)**-30)
            self.assertLess(abs(s-mp.mpf(area.numerator)/area.denominator), mp.mpf(10)**-30)
            for x, y in zip(cen, centre):
                self.assertLess(abs(x-mp.mpf(y.numerator)/y.denominator), mp.mpf(10)**-30)

    def test_three_boxes_by_their_cells(self):
        # The first Boolean in frames of different axes (XY and SIDE); the
        # third crossing both, sharing none of their planes.
        a = (0, 4, 0, 3, 0, 2)
        b = (1, 3, -1, 2, 1, 3)
        c = (2, 5, 1, 4, -1, 1.5)
        for op1 in ('fuse', 'cut', 'common'):
            for swapped in (False, True):
                self.check_boxes(a, b, c, op1, swapped)

    def test_boxes_sharing_planes(self):
        # The third's faces on the first's planes (both orientations): a
        # face on several prisms counted once.
        a = (0, 4, 0, 3, 0, 2)
        b = (0, 2, 0, 3, 1, 3)
        c = (2, 4, 0, 3, 0, 2)
        for op1 in ('fuse', 'cut'):
            self.check_boxes(a, b, c, op1, False)

    def test_far_third_prism(self):
        # A cylinder far from a box less a tilted hole: the fuse adds its
        # closed forms, the cut and the swapped common leave nothing of it.
        obj, tool = [fixtures.prism(*x, op=91+k) for k, x in enumerate(fixtures.HOLED)]
        third = prism([disc(0.0, 0.0, 1.0)], at('XY', (30, 0, 0)), 0.0, 2.0, 94)
        chain = ref.Chain(obj, tool, third, 'cut', False)
        pair = ref.cref.Pair(obj, tool)
        vx, _ = pair.volumes()['cut']
        ax = pair.area('cut')
        vc, _, ac = chain.prisms[2].measures()
        n, v, s, _ = chain.result('fuse')
        self.assertEqual(n, 2)
        self.assertLess(abs(v-vx-vc), mp.mpf(10)**-30)
        self.assertLess(abs(s-ax-ac), mp.mpf(10)**-30)
        n, v, s, _ = chain.result('cut')
        self.assertLess(abs(v-vx)+abs(s-ax), mp.mpf(10)**-30)
        self.assertEqual(chain.result('common')[0], 0)

    def test_common_is_symmetric(self):
        first = fixtures.QUARTER
        obj, tool = [fixtures.prism(*x, op=91+k) for k, x in enumerate(first)]
        third = prism([disc(0.0, 0.0, 1.5)], at('XY', (0, 0, 0)), 1.0, 7.0, 94)
        a = ref.Chain(obj, tool, third, 'common', False).result('common')
        b = ref.Chain(obj, tool, third, 'common', True).result('common')
        self.assertEqual(a[0], b[0])
        self.assertLess(abs(a[1]-b[1])+abs(a[2]-b[2]), mp.mpf(10)**-30)

    def test_cases(self):
        listed = fixtures.cases()
        fixtures.validate(listed)
        self.assertEqual(len(listed), 30)
        self.assertEqual(sum(1 for c in listed if c.kind == 'degenerate'), 6)
        self.assertTrue(any(c.swapped for c in listed))
        # The first Boolean's frames always differ in their axes; at most two
        # axis directions per chain (the reference's slicing).
        for c in listed:
            ref.Slicing3([ref.Prism(x) for x in (c.obj, c.tool, c.third)])


if __name__ == '__main__':
    unittest.main()
