#!/usr/bin/env python3
"""Deliberate failures of the STEP bridge (compare_step.py): wrong counts,
verdicts, measures, centres or bodies on either side are differences, and a
kernel enclosure that misses the reference is caught."""
import unittest

import compare_step as bridge


ROWS = [(70, 'solid', [2, 3, 3, 3, 1, 1], 942.4777960769379, 534.0707511102648, [0.0, 0.0, 6.0])]
NATIVE = ('done', [('solid', [2, 3, 3, 3, 1, 1], True, 942.4777960769379, 534.0707511102648,
                    [0.0, 0.0, 6.0], 1e-7)])


def rust(volume=(942.47779607, 942.47779608), counts='2 3 3 3 1 1', status='ok'):
    enclosure = [*volume, 534.07075111, 534.07075112, -1e-15, 1e-15, -1e-15, 1e-15, 5.99999999, 6.00000001]
    return [['70', 'solid', status, *counts.split(), *map(repr, enclosure)]]


class StepBridge(unittest.TestCase):
    def test_agreement(self):
        self.assertEqual(bridge.differences(NATIVE, ROWS), [])
        self.assertEqual(bridge.rust_differences(rust(), ROWS, NATIVE), [])

    def test_native_differences(self):
        body = NATIVE[1][0]
        wrong = lambda **kw: ('done', [tuple(kw.get(k, v) for k, v in zip(
            ['kind', 'counts', 'valid', 'volume', 'area', 'centre', 'tol'], body))])
        self.assertEqual(bridge.differences(wrong(counts=[2, 3, 3, 3, 1, 0]), ROWS), ['counts'])
        self.assertEqual(bridge.differences(wrong(valid=False), ROWS), ['invalid'])
        self.assertEqual(bridge.differences(wrong(volume=942.48), ROWS), ['measure'])
        self.assertEqual(bridge.differences(wrong(centre=[0.0, 0.0, 6.001]), ROWS), ['centre'])
        self.assertEqual(bridge.differences(wrong(kind='sheet'), ROWS), ['bodies'])
        self.assertEqual(bridge.differences(('read_status_3', []), ROWS), ['not_done'])

    def test_kernel_differences(self):
        self.assertEqual(bridge.rust_differences(rust(volume=(943.0, 944.0)), ROWS, NATIVE),
                         ['rust_outside_native', 'rust_outside_reference'])
        self.assertEqual(bridge.rust_differences(rust(counts='2 3 3 3 1 0'), ROWS, NATIVE), ['rust_counts'])
        self.assertEqual(bridge.rust_differences([['70', 'solid', 'unsupported', 'ELLIPSE']], ROWS, NATIVE),
                         ['rust_rejected'])
        self.assertEqual(bridge.rust_differences([['70', 'solid', 'invalid', 'edge_not_c1:edge', '0',
                                                   'face_not_c1:face', '1', 'edge_not_c1:edge', '2']],
                                                 ROWS, NATIVE), ['rust_invalid:edge_not_c1,face_not_c1'])
        self.assertEqual(bridge.rust_differences([['error', 'line', '1']], ROWS, NATIVE), ['rust_error'])
        self.assertEqual(bridge.rust_differences([], ROWS, NATIVE), ['rust_bodies'])


if __name__ == '__main__':
    unittest.main()
