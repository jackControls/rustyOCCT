"""Reviewed native differences must never weaken the independent Rust check."""
from copy import deepcopy
import json
import math
import struct
import unittest
from compare_splines import ROOT, cases, compare_observations, reviewed_divergence


class SplineOracleTests(unittest.TestCase):
    def setUp(self):
        self.review=json.loads((ROOT/'rust/fixtures/occt-spline-divergences.json').read_text())[0]
        self.row=next(row for row in cases().splitlines() if row.split()[0]==self.review['case'])
        self.native=struct.unpack('>d',bytes.fromhex(self.review['native_bits']))[0]
        self.rust=-0.47077685157054794

    def check(self, **override):
        args=dict(oracle=self.review['oracle'],label=self.review['case'],component=6,native=self.native,rust=self.rust,input_row=self.row)
        args.update(override)
        return reviewed_divergence(**args)

    def test_review_recomputes_independent_exact_value(self):
        result=self.check()
        self.assertIsNotNone(result)
        self.assertLess(result['exact_upper'],self.native)

    def test_changed_input_cannot_use_review(self):
        self.assertIsNone(self.check(input_row=self.row+' '))

    def test_changed_native_observation_cannot_use_review(self):
        self.assertIsNone(self.check(native=math.nextafter(self.native,math.inf)))
        self.assertIsNone(self.check(oracle='OCCT 7.9.3'))
        self.assertIsNone(self.check(component=7))

    def test_wrong_rust_value_cannot_use_review(self):
        self.assertIsNone(self.check(rust=self.native))
        self.assertIsNone(self.check(rust=math.nextafter(self.rust,-math.inf)))

    def test_reviewed_case_is_not_counted_as_matching(self):
        name=self.review['case']
        expected={name:[0.]*6+[self.native,0.,0.]}
        actual=deepcopy(expected)
        actual[name][6]=self.rust
        result=compare_observations(self.review['oracle'],{name:self.row},expected,actual)
        self.assertEqual((result['matched_cases'],result['reviewed_divergence_cases'],result['unexpected_mismatch_cases']),(0,1,0))
        actual[name][6]=1.
        result=compare_observations(self.review['oracle'],{name:self.row},expected,actual)
        self.assertEqual(result['unexpected_mismatch_cases'],1)


if __name__=='__main__': unittest.main()
