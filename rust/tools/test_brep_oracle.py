#!/usr/bin/env python3
"""Deliberate failures for native BRepCheck decoding, classification and review."""
import unittest

from brep_reference import CORRESPONDING, compare_native, decode_native
from compare_brep import original_capture, review_for, same_inputs


class BrepOracle(unittest.TestCase):
    def test_original_observations_and_inputs_remain_unchanged(self):
        original_capture()

    def test_regenerated_inputs_may_only_drift_in_the_last_bits(self):
        same_inputs('v 0.5 1 2\n', 'v 0.49999999999999994 1 2\n')
        for bad in ['v 0.5 1 3\n', 'v 0.5 1\n', 'e 0.5 1 2\n', 'v 0.5 1 2 3\n', 'v 0.5 1 nan\n']:
            with self.assertRaises(ValueError):
                same_inputs('v 0.5 1 2\n', bad)

    def test_native_rows_decode_or_fail(self):
        valid, statuses, absent = decode_native('box_extra_vertex R 1 v8:absent', 'box_extra_vertex')
        self.assertEqual((valid, statuses, absent), (True, {}, {'v8'}))
        valid, statuses, _ = decode_native('x R 0 e1:8:11 w2.0:32', 'x')
        self.assertEqual((valid, statuses), (False, {'e1': {8, 11}, 'w2.0': {32}}))
        for bad in ['y R 0', 'x E native_exception', 'x R 2', 'x R 0 e1', 'x R 0 e1:0',
                    'x R 0 e1:37', 'x R 0 e1:eight', '']:
            with self.assertRaises(ValueError):
                decode_native(bad, 'x')

    def test_differences_are_classified(self):
        row = lambda text: decode_native('x '+text, 'x')
        free = ['free_edge:edge 1', 'free_edge:edge 2']
        self.assertEqual(compare_native(free, row('R 0 s0:28'), []), [])
        self.assertEqual(compare_native([], row('R 1'), []), [])
        self.assertEqual(compare_native(free, row('R 1'), []), ['unmatched_free_edge', 'verdict'])
        self.assertEqual(compare_native([], row('R 0 s0:28'), []), ['verdict'])
        self.assertEqual(compare_native(free, row('R 0 s0:32'), []), ['unmatched_free_edge'])
        # Classes without a BRepCheck counterpart compare by verdict only.
        self.assertEqual(compare_native(['euler:shell 0'], row('R 0 s0:29'), []), [])
        self.assertEqual(compare_native(free, row('R 0 s0:28'), ['2.0.1']), ['inexact_pcurve_encoding'])

    def test_correspondence_uses_real_statuses(self):
        # BRepCheck_Status has 37 values; 0 is NoError.
        for kind, codes in CORRESPONDING.items():
            self.assertTrue(codes and all(0 < c < 37 for c in codes), kind)

    def test_review_requires_rationale_and_exact_fingerprint(self):
        evidence = {'case': 'x', 'native_stdout_sha256': 'a', 'differences': ['verdict']}
        self.assertIsNone(review_for(evidence, [dict(evidence)]))
        review = dict(evidence, reason='r', independent_evidence='e')
        self.assertIs(review_for(evidence, [review]), review)
        self.assertIsNone(review_for(dict(evidence, native_stdout_sha256='b'), [review]))
        self.assertIsNone(review_for(dict(evidence, differences=['verdict', 'unmatched_uv_gap']), [review]))


if __name__ == '__main__':
    unittest.main()
