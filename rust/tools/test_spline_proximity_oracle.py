#!/usr/bin/env python3
"""Deliberate failures for independent minima certification and native reporting."""
from copy import deepcopy
from fractions import Fraction as F
import json
import math
from pathlib import Path
import unittest

from compare_spline_proximity import review_for, original_capture, verify_family_libraries
from spline_proximity_reference import verify_rust, decode_rust, decode_native, compare_native

ROOT = Path(__file__).resolve().parents[1]
CASES = {c['name']: c for c in json.loads((ROOT/'fixtures/spline-proximity-inputs.json').read_text())['cases']}
LINE = 'line_interior R 1 0 4 4 P .5 .5 .5 .5 0 0 0 0'


class SplineProximityOracle(unittest.TestCase):
    def test_newer_api_loads_its_actual_toolkit_without_unused_legacy_wrapper(self):
        for suffix in ['.so.8.1', '.8.1.dylib']:
            base = [{'path': '/pinned/lib/libTKGeomBase'+suffix}]
            legacy = base+[{'path': '/pinned/lib/libTKGeomAlgo'+suffix}]
            verify_family_libraries('extremapc', base)
            verify_family_libraries('legacy', legacy)
            with self.assertRaises(ValueError):
                verify_family_libraries('legacy', base)
            for family in ['legacy', 'extremapc']:
                with self.assertRaises(ValueError):
                    verify_family_libraries(family, legacy[1:])
        with self.assertRaises(ValueError):
            verify_family_libraries('unknown', legacy)

    def test_original_observations_and_inputs_remain_unchanged(self):
        original_capture()

    def test_exact_line_witness_and_deliberately_wrong_coordinate(self):
        verify_rust(CASES['line_interior'], LINE)
        bad = LINE.replace('P .5 .5 .5 .5', 'P .5 .5 .25 .25')
        with self.assertRaises(ValueError):
            verify_rust(CASES['line_interior'], bad)

    def test_complete_irrational_tie_cannot_lose_one_parameter(self):
        hi = math.sqrt(.5)
        lo = math.nextafter(hi, -math.inf)
        points = [f'P {-hi:.17g} {-lo:.17g} {-hi:.17g} {-lo:.17g} .5 .5 0 0',
                  f'P {lo:.17g} {hi:.17g} {lo:.17g} {hi:.17g} .5 .5 0 0']
        name = 'parabola_two_equal_minima'
        verify_rust(CASES[name], f'{name} R 2 0 .75 .75 '+' '.join(points))
        with self.assertRaises(ValueError):
            verify_rust(CASES[name], f'{name} R 1 0 .75 .75 '+points[0])

    def test_wrong_distance_and_wide_bounds_fail(self):
        for bad in [LINE.replace('4 4 P', '5 5 P'), LINE.replace('.5 .5', '.49 .51', 1)]:
            with self.assertRaises(ValueError):
                verify_rust(CASES['line_interior'], bad)

    def test_malformed_nonfinite_and_trailing_output_fail(self):
        for bad in [LINE+' trailing', LINE.replace('4 4', 'nan nan'), LINE.replace('R 1', 'R -1'), LINE[:-2]]:
            with self.assertRaises((ValueError, StopIteration)):
                decode_rust(bad)

    def test_native_whole_interval_is_separate_from_two_endpoints(self):
        actual = {'distance': (1., 1.), 'points': [], 'intervals': [[F(0), F(1)]], 'domain': [F(0), F(1)]}
        native = {'points': [{'parameter': 0., 'xyz': [1., 0., 0.], 'distance': 1.},
                             {'parameter': 1., 'xyz': [0., 1., 0.], 'distance': 1.}],
                  'infinite': False, 'infinite_distance': None}
        self.assertEqual(compare_native(actual, native), ['whole_minimum_interval_not_represented'])
        native.update(infinite=True, infinite_distance=1.)
        self.assertEqual(compare_native(actual, native), [])
        actual['intervals'] = [[F(0), F(1, 2)]]
        self.assertIn('minimum_interval_ranges_not_represented', compare_native(actual, native))

    def test_native_duplicate_cannot_cover_distinct_periodic_parameters(self):
        point = {'parameter': (0., 0.), 'xyz': [(0., 0.)]*3}
        actual = {'distance': (0., 0.), 'points': [point, dict(point, parameter=(1e-10, 1e-10))], 'intervals': []}
        point = {'parameter': 0., 'xyz': [0., 0., 0.], 'distance': 0.}
        native = {'points': [point, point], 'infinite': False}
        self.assertEqual(compare_native(actual, native), ['minimum_witness_coverage_1_of_2'])
        native['points'][1] = dict(point, parameter=1e-10)
        self.assertEqual(compare_native(actual, native), [])

    def test_native_status_identity_and_selected_index_are_checked(self):
        good = 'x R 0 0 1 0 1 .5 1 0 4 .5 0 0 M 0 4\nx R 1 0 1 0 1 .5 1 0 4 .5 0 0 M 0 4\n'
        decode_native('extremapc', good, 0., 1.)
        no_solution = 'x R 0 3 1 0 0\nx R 1 3 1 0 0\n'
        decode_native('extremapc', no_solution, 0., 1.)
        for bad in [good.replace('M 0 4', 'M 1 4'), good.replace('R 1 0 1', 'R 1 1 1'), good.replace('x R 1', 'y R 1')]:
            with self.assertRaises(ValueError):
                decode_native('extremapc', bad, 0., 1.)

    def test_review_requires_every_fingerprint_and_evidence(self):
        evidence = dict(family='legacy', case='example', source_reference='pin', input_sha256='input',
                        exact_sha256='exact', native_stdout_sha256='stdout', native_stderr_sha256='stderr',
                        differences=['whole_minimum_interval_not_represented'], oracle='OCCT 8.1.0')
        review = dict(evidence, reason='contract difference', independent_evidence='exact polynomial identity')
        self.assertEqual(review_for(evidence, [review]), review)
        for key in evidence:
            changed = deepcopy(review)
            changed[key] = 'changed'
            self.assertIsNone(review_for(evidence, [changed]))
        for key in ['reason', 'independent_evidence']:
            changed = deepcopy(review)
            del changed[key]
            self.assertIsNone(review_for(evidence, [changed]))


if __name__ == '__main__':
    unittest.main()
