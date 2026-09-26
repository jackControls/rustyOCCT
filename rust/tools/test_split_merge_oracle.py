#!/usr/bin/env python3
"""Deliberate failures for the native split/fuse comparison."""
import copy
import functools
import json
import subprocess
import unittest

from compare_occt import ROOT
from compare_split_merge import (ORIGINAL, compare, original_capture, parse_native, parse_rust,
                                 review_for)
import generate_split_merge_fixtures as generator

CASE = 'split_polygon_4'


@functools.lru_cache(maxsize=None)
def rust():
    """The Rust probe's rows for the captured scenario's native steps."""
    scenarios, _ = generator.generate()
    s = next(s for s in scenarios if s.name == CASE)
    subprocess.run(['cargo', '+stable', 'build', '--release', '--locked', '--example', 'split_merge_probe'],
                   cwd=ROOT, check=True, capture_output=True)
    text = subprocess.run([str(ROOT/'target/release/examples/split_merge_probe')],
                          input=generator.encode(generator.native_scenario(s))+'\n',
                          text=True, capture_output=True, check=True).stdout
    return parse_rust(text)[CASE]


def native():
    records = json.loads((ORIGINAL/'native.json').read_text())
    record = next(r for r in records if r['case'] == CASE)
    return parse_native(record['stdout'], CASE)


class SplitMergeOracle(unittest.TestCase):
    def test_original_observations_and_inputs_remain_unchanged(self):
        original_capture()

    def test_matching_case(self):
        self.assertEqual(compare(native(), copy.deepcopy(rust())), [])

    def test_relations_are_checked(self):
        r = copy.deepcopy(rust())
        # A split child moved off its parent.
        q = next(q for q in r['split']['queries'] if q['kind'] == 'split' and q['sig'].startswith('F'))
        q['targets'][1] = q['targets'][1].replace(' 2.75', ' 2.5')
        self.assertIn('split_split_F', compare(native(), r))
        # A deleted cap entity reported unchanged.
        r = copy.deepcopy(rust())
        q = next(q for q in r['fuse']['queries'] if q['kind'] == 'deleted')
        q['kind'] = 'unchanged'
        q['targets'] = [q['sig']]
        self.assertTrue(any(d.startswith('fuse_unchanged_') for d in compare(native(), r)))
        # A merge into the wrong face.
        r = copy.deepcopy(rust())
        merged = [q for q in r['fuse']['queries'] if q['kind'] == 'merged' and q['sig'].startswith('F')]
        merged[0]['targets'], merged[-1]['targets'] = merged[-1]['targets'], merged[0]['targets']
        self.assertIn('fuse_merged_F', compare(native(), r))

    def test_generated_outputs_and_counts_are_checked(self):
        r = copy.deepcopy(rust())
        r['split']['generated'] = [g for g in r['split']['generated'] if g['role'] != 'CutEdge']
        self.assertIn('split_native_generated_unmatched', compare(native(), r))
        r = copy.deepcopy(rust())
        g = next(g for g in r['split']['generated'] if g['role'] == 'CutFace')
        g['target'] = g['target'].replace(' 1.5', ' 1.25')
        self.assertIn('split_generated_CutFace', compare(native(), r))
        r = copy.deepcopy(rust())
        r['split']['counts'][0] = '8 12 6 6 1 2'
        self.assertIn('split_counts', compare(native(), r))
        r = copy.deepcopy(rust())
        r['split']['outputs'].append('V 1 2 3')
        self.assertIn('split_rust_output_unmatched', compare(native(), r))
        r = copy.deepcopy(rust())
        del r['composed']
        self.assertEqual(compare(native(), r), ['stages'])

    def test_native_side_is_checked(self):
        n = native()
        q = next(q for q in n['stages']['split']['queries'] if len(q['modified']) == 2)
        q['modified'] = q['modified'][:1]
        self.assertTrue(any(d.startswith('split_split_') for d in compare(n, copy.deepcopy(rust()))))
        n = native()
        n['stages']['split']['bodies'][0] = (n['stages']['split']['bodies'][0][0], False)
        self.assertIn('split_native_invalid', compare(n, copy.deepcopy(rust())))
        n = native()
        n['inputs_valid'] = False
        self.assertIn('native_invalid_input', compare(n, copy.deepcopy(rust())))
        # A structure-only row set aside must still be accounted for by counts.
        n = native()
        n['stages']['split']['structure'].add(n['stages']['split']['queries'][0]['sig'])
        self.assertIn('split_unchanged_V', compare(n, copy.deepcopy(rust())))

    def test_malformed_native_output_fails(self):
        for bad in ['', f'R {CASE}\n', f'R other\nend', f'R {CASE}\nX 1\nend',
                    f'R {CASE}\nQ split P v 0 mod 1 gen 0 del 0 | V 0 0 0\nend']:
            with self.assertRaises((ValueError, IndexError)):
                parse_native(bad, CASE)

    def test_review_requires_rationale_and_exact_fingerprint(self):
        evidence = {'case': 'x', 'native_stdout_sha256': 'a', 'differences': ['fuse_counts']}
        self.assertIsNone(review_for(evidence, [dict(evidence)]))
        review = dict(evidence, reason='r', independent_evidence='e')
        self.assertIs(review_for(evidence, [review]), review)
        self.assertIsNone(review_for(dict(evidence, native_stdout_sha256='b'), [review]))


if __name__ == '__main__':
    unittest.main()
