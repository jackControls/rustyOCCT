"""Independent checks of the degree-elevation coefficient oracle itself."""
from fractions import Fraction as F
from math import comb
import hashlib
import json
from pathlib import Path
import unittest
from degree_elevation_reference import elevate, elevated_axis, elevate_surface, count
from surface_knot_reference import domain


class DegreeElevationReferenceTests(unittest.TestCase):
    def test_bezier_coefficients_agree_with_direct_binomial_identity(self):
        for p in (1, 2, 5):
            controls = tuple((F(i*(i+1), 3), F(-i*i, 7), F(2*i-3), F(i+1, 2))
                             for i in range(p+1))
            original = (p, False, (F(-3), F(2)), (p+1, p+1), controls)
            for q in (p+1, 25):
                expected = tuple(tuple(sum(
                    F(comb(p,i)*comb(q-p,j-i), comb(q,j))*controls[i][c]
                    for i in range(max(0,j-(q-p)), min(p,j)+1))
                    for c in range(4)) for j in range(q+1))
                self.assertEqual(elevate(original, q)[4], expected)

    def test_unclamped_axis_preserves_domain_and_complete_control_count(self):
        axis = (2, False, tuple(map(F, range(8))), (1,)*8)
        result = elevated_axis(axis, 3)
        self.assertEqual(result, (3, False, tuple(map(F, range(1,7))), (2,)*6))
        self.assertEqual(domain(result), domain(axis))
        self.assertEqual(count(result), 8)
        for target in (0, 1, 26):
            with self.assertRaises(ValueError):
                elevated_axis(axis, target)

    def test_tensor_orders_and_identity_preserve_all_fields(self):
        axes = ((1, False, (F(0),F(1)), (2,2)),
                (2, True, (F(0),F(1),F(3)), (1,2,1)))
        controls = tuple((F(i),F(i*i,3),F(5-i),F(i+1)) for i in range(6))
        surface = axes, controls
        self.assertEqual(elevate_surface(surface, 1, 2), surface)
        self.assertEqual(elevate_surface(surface, 3, 4),
                         elevate_surface(surface, 3, 4, order=(1,0)))

    def test_fixture_inputs_are_the_captured_preimplementation_inputs(self):
        fixtures = Path(__file__).resolve().parents[1]/'fixtures'
        record = json.loads((fixtures/'occt-degree-elevation-capture.json').read_text())
        for kind in ('curves', 'surfaces'):
            data = (fixtures/f'degree-elevation-{kind}.txt').read_bytes()
            capture = record[kind]
            self.assertEqual(hashlib.sha256(data).hexdigest(), capture['input_sha256'])
            self.assertFalse(capture['implementation_present'])
            self.assertEqual(capture['source_reference'], record['source_reference'])
            self.assertEqual(capture['head'], record['before_implementation_head'])


if __name__ == '__main__':
    unittest.main()
