#!/usr/bin/env python3
"""Deliberate failures for the native history comparison."""
import copy
import json
import unittest

from compare_history import ORIGINAL, compare, original_capture, parse_native, parse_rust, review_for

# The Rust history probe's rows for the captured 'circle' case.
RUST = """R circle
G end_cap 0 face 0 0 | F plane 28.274333882308138 0.9999999999999999 2.0 4.0
G start_cap 0 face 0 0 | F plane 28.274333882308138 0.9999999999999999 2.0 0.0
G wall 0 0 edge 0 | F cylinder 75.39822368615503 1.0 2.0 2.0
G bottom_edge 0 0 edge 0 | E circle 4.0 2.0 0.0 -2.0 2.0000000000000004 0.0 4.0 1.9999999999999993 0.0
G top_edge 0 0 edge 0 | E circle 4.0 2.0 4.0 -2.0 2.0000000000000004 4.0 4.0 1.9999999999999993 4.0
G seam_vertex 1 0 vertex 0 | V 4.0 2.0 4.0
G seam_vertex 0 0 vertex 0 | V 4.0 2.0 0.0
G seam 0 0 vertex 0 | E line 4.0 2.0 0.0 4.0 2.0 2.0 4.0 2.0 4.0
B S 113.09733552923255 1.0 2.0 2.0
end
"""


def native():
    records = json.loads((ORIGINAL/'native.json').read_text())
    record = next(r for r in records if r['case'] == 'circle')
    return parse_native(record['stdout'], 'circle')


class HistoryOracle(unittest.TestCase):
    def test_original_observations_and_inputs_remain_unchanged(self):
        original_capture()

    def test_matching_case_and_deliberate_failures(self):
        rust = parse_rust(RUST)['circle']
        self.assertEqual(compare(native(), rust), [])
        # Swapping the start and end seam vertices moves both.
        swapped = copy.deepcopy(rust)
        rel = swapped['relations']
        rel[5], rel[6] = (rel[5][0], 0, *rel[5][2:]), (rel[6][0], 1, *rel[6][2:])
        self.assertEqual(compare(native(), swapped), ['signature_vertex_first', 'signature_vertex_last'])
        # A missing relation, a moved face, an extra relation.
        missing = copy.deepcopy(rust)
        del missing['relations'][2]
        self.assertEqual(compare(native(), missing), ['count_F', 'missing_rust_edge_gen'])
        moved = copy.deepcopy(rust)
        moved['relations'][2] = moved['relations'][2][:3]+('F cylinder 75.39822368615503 1.0 2.0 2.5',)
        self.assertEqual(compare(native(), moved), ['signature_edge_gen'])
        extra = copy.deepcopy(rust)
        extra['relations'].append(('wall', 0, '0 edge 7', 'F plane 1 0 0 0'))
        self.assertEqual(compare(native(), extra), ['count_F', 'unreached_rust_relation'])
        body = copy.deepcopy(rust)
        body['bodies'][0] = 'S 113.2 1.0 2.0 2.0'
        self.assertEqual(compare(native(), body), ['body_signature'])

    def test_native_completeness_and_transforms_are_checked(self):
        rust = parse_rust(RUST)['circle']
        incomplete = native()
        incomplete['outputs'][0] = (incomplete['outputs'][0][0], 'none')
        self.assertEqual(compare(incomplete, rust), ['native_history_incomplete'])
        invalid = native()
        invalid['valid'] = False
        self.assertEqual(compare(invalid, rust), ['native_invalid'])
        a, b = 'V 1 2 3', 'V 4 5 6'
        n, r = native(), copy.deepcopy(rust)
        n['transforms'][0] = [(a, [b])]
        r['transforms'][0] = [(a, b)]
        self.assertEqual(compare(n, r), [])
        r['transforms'][0] = [(a, 'V 4 5 7')]
        self.assertEqual(compare(n, r), ['transform0_pair'])
        n['transforms'][0] = [(a, [b, b])]
        self.assertEqual(compare(n, r), ['transform0_arity'])
        r['transforms'][1] = []
        self.assertIn('transform_steps', compare(n, r))

    def test_malformed_native_output_fails(self):
        for bad in ['', 'R circle 1\n', 'R other 1\nend', 'R circle 1\nX 1\nend']:
            with self.assertRaises((ValueError, IndexError)):
                parse_native(bad, 'circle')

    def test_review_requires_rationale_and_exact_fingerprint(self):
        evidence = {'case': 'x', 'native_stdout_sha256': 'a', 'differences': ['count_F']}
        self.assertIsNone(review_for(evidence, [dict(evidence)]))
        review = dict(evidence, reason='r', independent_evidence='e')
        self.assertIs(review_for(evidence, [review]), review)
        self.assertIsNone(review_for(dict(evidence, native_stdout_sha256='b'), [review]))


if __name__ == '__main__':
    unittest.main()
