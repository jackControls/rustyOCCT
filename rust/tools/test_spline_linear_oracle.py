#!/usr/bin/env python3
"""Deliberate failures for independent preimage certification and native reporting."""
from fractions import Fraction as F
import unittest

from compare_spline_linear import original_capture, queries, review_for, segmented
from spline_linear_reference import compare_native, decode_native, native_cases, verify_rust

CASES = {c['name']: c for c in native_cases()}
SECANT = ('parabola_secant R 2 0 P -1.0 -1.0 0.25 0.25 -1.0 -1.0 1.0 1.0 0.0 0.0 '
          'P 1.0 1.0 0.75 0.75 1.0 1.0 1.0 1.0 0.0 0.0')
RETRACE = ('collinear_disconnected_parameters R 0 2 '
           'I -1.4142135623730951 -1.414213562373095 1.0 1.0 2.0 2.0 0.0 0.0 0.0 0.0 '
           '-1.0 -1.0 0.0 0.0 1.0 1.0 0.0 0.0 0.0 0.0 '
           'I 1.0 1.0 0.0 0.0 1.0 1.0 0.0 0.0 0.0 0.0 '
           '1.414213562373095 1.4142135623730951 1.0 1.0 2.0 2.0 0.0 0.0 0.0 0.0')


def native(case, *parts):
    return decode_native(f'{case} R 1 {len(parts)} '+' '.join(parts), CASES[case])


class SplineLinearOracle(unittest.TestCase):
    def test_original_observations_and_inputs_remain_unchanged(self):
        original_capture()

    def test_certified_secant_and_deliberately_wrong_rows(self):
        verify_rust(CASES['parabola_secant'], SECANT)
        for bad in [SECANT.replace('P 1.0 1.0 0.75', 'P 0.5 0.5 0.75'),
                    SECANT.replace('0.25 0.25', '0.5 0.5'),
                    SECANT.replace('P -1.0 -1.0', 'P -1.0 -0.5'),
                    'parabola_secant R 1 0 P 1.0 1.0 0.75 0.75 1.0 1.0 1.0 1.0 0.0 0.0']:
            with self.assertRaises(ValueError):
                verify_rust(CASES['parabola_secant'], bad)

    def test_irrational_interval_endpoints_cannot_move_one_float(self):
        verify_rust(CASES['collinear_disconnected_parameters'], RETRACE)
        with self.assertRaises(ValueError):
            verify_rust(CASES['collinear_disconnected_parameters'],
                        RETRACE.replace('I -1.4142135623730951 -1.414213562373095',
                                        'I -1.414213562373095 -1.4142135623730949'))

    def test_native_differences_are_classified(self):
        actual = verify_rust(CASES['parabola_secant'], SECANT)
        # Native segment edge length is 4, so v = 4s.
        both = ['P 1 -1 1 -1 1 0 -1 1 0', 'P 1 1 3 1 1 0 1 1 0']
        self.assertEqual(compare_native(actual, native('parabola_secant', *both)), [])
        self.assertEqual(compare_native(actual, native('parabola_secant', both[0])),
                         ['missing_isolated_point'])
        self.assertEqual(compare_native(actual, native('parabola_secant', *both, 'P 1 0 2 0 0 0 0 0 0')),
                         ['extra_native_vertex'])
        # A native vertex cannot match a correct parameter on the wrong line parameter.
        self.assertEqual(compare_native(actual, native('parabola_secant', both[0], 'P 1 1 2 1 1 0 1 1 0')),
                         ['extra_native_vertex', 'missing_isolated_point'])
        retrace = verify_rust(CASES['collinear_disconnected_parameters'], RETRACE)
        first = 'I 1 -1.4142136 -0.99999985 1 0 1'
        second = 'I 1 1 1.4142136 1 0 1'
        name = 'collinear_disconnected_parameters'
        self.assertEqual(compare_native(retrace, native(name, first, second)), [])
        self.assertEqual(compare_native(retrace, native(name, first)), ['missing_overlap'])
        self.assertEqual(compare_native(retrace, native(name, first, 'P 1 1.2 0.44 1.44 0 0 1.44 0 0')),
                         ['overlap_reported_as_vertex'])
        self.assertEqual(compare_native(retrace, native(name, first, 'I 1 1.2 1.4142136 1 0 1')),
                         ['overlap_range'])

    def test_periodic_turn_windows_keep_integer_offsets(self):
        windows = queries(CASES['periodic_multiple_turns'])
        self.assertEqual([(w['range'], offset) for w, offset in windows],
                         [([F(0), F(3)], F(-3)), ([F(0), F(3)], F(0)), ([F(0), F(3)], F(3))])
        self.assertEqual(len(queries(CASES['periodic_one_turn'])), 1)
        parts = segmented(CASES['chebyshev_25'])
        self.assertEqual(len(parts), 32)
        self.assertEqual((parts[0][0]['range'][0], parts[-1][0]['range'][1]), (F(-1), F(1)))

    def test_review_requires_rationale_and_exact_fingerprint(self):
        evidence = {'case': 'x', 'native_stdout_sha256': 'a', 'differences': ['missing_overlap']}
        self.assertIsNone(review_for(evidence, [dict(evidence)]))
        review = dict(evidence, reason='r', independent_evidence='e')
        self.assertIs(review_for(evidence, [review]), review)
        self.assertIsNone(review_for(dict(evidence, native_stdout_sha256='b'), [review]))


if __name__ == '__main__':
    unittest.main()
