#!/usr/bin/env python3
"""Deliberate failures for the native history comparison."""
import copy
import json
import unittest

from compare_history import ORIGINAL, compare, original_capture, parse_native, parse_rust, review_for

# The Rust history probe's rows for the captured 'circle' case: seamless.
RUST = """R circle
G end_cap 0 face 0 0 | F plane 28.274333882308138 0.9999999999999999 2.0 4.0
G start_cap 0 face 0 0 | F plane 28.274333882308138 0.9999999999999999 2.0 0.0
G region 0 face 0 0 | S 113.09733552923255 1.0 2.0 2.0
G wall 0 0 edge 0 | F cylinder 75.39822368615503 1.0 2.0 2.0
G bottom_edge 0 0 edge 0 | E circle 4.0 2.0 0.0 -2.0 2.0000000000000004 0.0 4.0 1.9999999999999993 0.0
G top_edge 0 0 edge 0 | E circle 4.0 2.0 4.0 -2.0 2.0000000000000004 4.0 4.0 1.9999999999999993 4.0
B S 113.09733552923255 1.0 2.0 2.0
C 2 3 3 3 1 1
end
"""

# The rows the current oracle appends to the original capture for 'circle':
# its counts and its structure-only seam and seam vertices.
STRUCTURE = """N 2 3 3 3 1 1
S - E line 4 2 0 4 2 2 4 2 4
S - V 4 2 0
S - V 4 2 4
"""


def native(structure=STRUCTURE):
    records = json.loads((ORIGINAL/'native.json').read_text())
    record = next(r for r in records if r['case'] == 'circle')
    stdout = record['stdout']
    assert stdout.endswith('end\n')
    return parse_native(stdout[:-len('end\n')]+structure+'end\n', 'circle')


class HistoryOracle(unittest.TestCase):
    def test_original_observations_and_inputs_remain_unchanged(self):
        original_capture()

    def test_matching_case_and_deliberate_failures(self):
        rust = parse_rust(RUST)['circle']
        self.assertEqual(compare(native(), rust), [])
        # A missing relation, a moved face, an extra relation.
        missing = copy.deepcopy(rust)
        del missing['relations'][3]
        self.assertEqual(compare(native(), missing), ['count_F', 'missing_rust_edge_gen'])
        moved = copy.deepcopy(rust)
        moved['relations'][3] = moved['relations'][3][:3]+('F cylinder 75.39822368615503 1.0 2.0 2.5',)
        self.assertEqual(compare(native(), moved), ['signature_edge_gen'])
        extra = copy.deepcopy(rust)
        extra['relations'].append(('wall', 0, '0 edge 7', 'F plane 1 0 0 0'))
        self.assertEqual(compare(native(), extra), ['count_F', 'unreached_rust_relation'])
        body = copy.deepcopy(rust)
        body['bodies'][0] = 'S 113.2 1.0 2.0 2.0'
        self.assertEqual(compare(native(), body), ['body_signature'])
        # The solid is generated from the face as the region.
        region = copy.deepcopy(rust)
        del region['relations'][2]
        self.assertEqual(compare(native(), region), ['region_signature'])
        region = copy.deepcopy(rust)
        region['relations'][2] = region['relations'][2][:3]+('S 113.2 1.0 2.0 2.0',)
        self.assertEqual(compare(native(), region), ['region_signature', 'unreached_rust_relation'])

    def test_structure_only_rows_are_verified_by_count_synthesis(self):
        rust = parse_rust(RUST)['circle']
        # Without the classification the seam and its vertices are unmapped.
        self.assertEqual(compare(native('N 2 3 3 3 1 1\n'), rust),
                         ['count_E', 'count_V', 'missing_rust_vertex_first', 'missing_rust_vertex_gen',
                          'missing_rust_vertex_last'])
        # Classifying a real edge as structure hides it, and the counts catch
        # a synthesis that disagrees with native.
        self.assertEqual(compare(native(STRUCTURE+'S - E circle 4 2 0 -2 2.0000000000000004 0 4 '
                                                   '1.9999999999999993 0\n'), rust),
                         ['count_E', 'unreached_rust_relation'])
        wrong = copy.deepcopy(rust)
        wrong['counts'] = '2 2 3 3 1 1'
        self.assertEqual(compare(native(), wrong), ['synthesized_counts'])
        self.assertEqual(compare(native(STRUCTURE.replace('N 2 3 3 3 1 1\n', '')), rust),
                         ['synthesized_counts'])

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
        # Structure-only pairs are set aside; the region moves with the body.
        n, r = native(), copy.deepcopy(rust)
        seam = 'E line 4 2 0 4 2 2 4 2 4'
        n['transforms'][0] = [(a, [b]), (seam, ['E line 5 2 0 5 2 2 5 2 4'])]
        n['structure']['0'] = {seam}
        body = n['bodies'][0][0]
        moved = 'S 113.09733552923254 2 2 2'
        n['bodies'].append((moved, True))
        r['bodies'].append(moved)
        r['transforms'][0] = [(a, b), (r['relations'][2][3], moved)]
        self.assertEqual(compare(n, r), [])
        r['transforms'][0][1] = (body, 'S 113.09733552923254 2 2 2.5')
        self.assertEqual(compare(n, r), ['transform0_region'])
        n['structure']['0'] = set()
        r['transforms'][0][1] = (body, moved)
        self.assertEqual(compare(n, r), ['transform0_arity'])

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
