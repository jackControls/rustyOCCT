#!/usr/bin/env python3
"""Deliberate failures for native BRepCheck decoding, classification and review."""
import unittest

from brep_reference import CORRESPONDING, compare_native, decode_native, structure_only
from compare_brep import (decode_tolerances, enclosure_capture, generate, original_capture, review_for,
                          same_inputs, same_tolerances)


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

    def test_structure_only_statuses_cannot_stand_for_an_issue(self):
        row = lambda text: decode_native('x '+text, 'x')
        gap = ['uv_gap:fin 3']
        self.assertEqual(compare_native(gap, row('R 0 e2:28'), []), [])
        self.assertEqual(compare_native(gap, row('R 0 e2:28'), [], {'e2'}), ['unmatched_uv_gap'])
        self.assertEqual(compare_native(gap, row('R 0 e2:28 w2.0:28'), [], {'e2'}), [])
        # A verdict never depends on where the statuses are.
        self.assertEqual(compare_native([], row('R 0 e2:28'), [], {'e2'}), ['verdict'])

    def test_seams_and_their_vertices_are_structure_only(self):
        cylinder = next(m for m in generate()[0] if m.name == 'cylinder')
        self.assertEqual(structure_only(cylinder), {'e2', 'v0', 'v1'})
        box = next(m for m in generate()[0] if m.name == 'box')
        self.assertEqual(structure_only(box), set())

    def test_correspondence_uses_real_statuses(self):
        # BRepCheck_Status has 37 values; 0 is NoError.
        for kind, codes in CORRESPONDING.items():
            self.assertTrue(codes and all(0 < c < 37 for c in codes), kind)

    def test_tolerance_rows_decode_or_fail(self):
        row = 'x T v0:1e-07:0 e0:1e-07 u0.0.0:3.7e-16 f0:1e-07'
        self.assertEqual(decode_tolerances(row, 'x'),
                         {'v0': (1e-7, 0.0), 'e0': (1e-7, None), 'u0.0.0': (None, 3.7e-16), 'f0': (1e-7, None)})
        for bad in ['y T v0:1:0', 'x N v0:1:0', 'x T v0:1', 'x T e0:1:2', 'x T u0:a', 'x T q0:1']:
            with self.assertRaises(ValueError):
                decode_tolerances(bad, 'x')

    def test_tolerance_observations_are_reproduced_or_refused(self):
        base = {'v0': (1e-7, 0.0), 'u0.0.0': (None, 3.7e-16), 'u0.0.1': (None, 2.0)}
        self.assertTrue(same_tolerances(base, dict(base, **{'u0.0.0': (None, 4.5e-16)})))
        self.assertTrue(same_tolerances(base, dict(base, **{'u0.0.1': (None, 2.0*(1+5e-10))})))
        for changed in [{'v0': (2e-7, 0.0)}, {'u0.0.1': (None, 2.0*(1+1e-8))}, {'u0.0.0': (None, 2e-15)}]:
            self.assertFalse(same_tolerances(base, dict(base, **changed)))
        self.assertFalse(same_tolerances(base, {'v0': (1e-7, 0.0)}))
        with self.assertRaises(ValueError):
            enclosure_capture({'box': {'v0': (1.0, 0.0)}})

    def test_review_requires_rationale_and_exact_fingerprint(self):
        evidence = {'case': 'x', 'native_stdout_sha256': 'a', 'differences': ['verdict']}
        self.assertIsNone(review_for(evidence, [dict(evidence)]))
        review = dict(evidence, reason='r', independent_evidence='e')
        self.assertIs(review_for(evidence, [review]), review)
        self.assertIsNone(review_for(dict(evidence, native_stdout_sha256='b'), [review]))
        self.assertIsNone(review_for(dict(evidence, differences=['verdict', 'unmatched_uv_gap']), [review]))


if __name__ == '__main__':
    unittest.main()
