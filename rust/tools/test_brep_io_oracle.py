#!/usr/bin/env python3
"""Deliberate failures for the .brep interop bridge and the independent reader."""
import json
import unittest

from brep_io_reference import IDENTITY, read, summary
from compare_brep_io import CAPTURE, PROPERTY_BOUND, close, original_capture, parse


def captured():
    out = {}
    for line in (CAPTURE/'native.txt').read_text().splitlines():
        w = line.split()
        if w[0] == 'F':
            current = out.setdefault(w[1], [])
        else:
            current.append((w[1], tuple(int(x) for x in w[2:8]), [float(x) for x in w[8:]]))
    return out


class BrepIoOracle(unittest.TestCase):
    def test_capture_is_unchanged_and_detects_drift(self):
        observed = captured()
        original_capture(observed)
        name = next(n for n, s in observed.items() if s)
        verdict, counts, props = observed[name][0]
        for changed in [('invalid', counts, props), (verdict, (0,)+counts[1:], props),
                        (verdict, counts, [props[0]*(1+1e-9)]+props[1:])]:
            bad = dict(observed, **{name: [changed]+observed[name][1:]})
            with self.assertRaises(ValueError):
                original_capture(bad)
        with self.assertRaises(ValueError):
            original_capture({k: v for k, v in observed.items() if k != name})

    def test_native_output_parses_or_fails(self):
        good = 'F a 1\nS valid 1 2 3 4 5 6 1 2 3 4 5\nF b 0\n'
        self.assertEqual(parse(good, ['a', 'b'])['a'][0][1], (1, 2, 3, 4, 5, 6))
        for bad in ['F a unreadable\nF b 0\n', 'F a 1\nS valid 1 2\nF b 0\n', 'F b 0\nF a 0\n',
                    'S valid 1 2 3 4 5 6 1 2 3 4 5\n', 'F a 0\n', 'F a failure x\nF b 0\n']:
            with self.assertRaises(ValueError):
                parse(bad, ['a', 'b'])

    def test_properties_compare_relatively(self):
        base = [1000.0, 600.0, 1.0, 2.0, 3.0]
        self.assertEqual(close(base, base), 0.0)
        self.assertLessEqual(close([1000.0*(1+1e-13)]+base[1:], base), PROPERTY_BOUND)
        self.assertGreater(close(base[:2]+[1.0+1e-9]+base[3:], base), PROPERTY_BOUND)
        self.assertFalse(close([float('nan')]+base[1:], base) <= PROPERTY_BOUND)

    def test_reader_follows_the_format(self):
        # A matrix record is always numbered; composite 2 reuses it squared.
        text = '''CASCADE Topology V1, (c) Matra-Datavision
Locations 2
1
1 0 0 5
0 1 0 0
0 0 1 0
2 1 2 0
Curve2ds 1
8 0 1 1 0 0 1 0
Curves 0
Polygon3D 0
PolygonOnTriangulations 0
Surfaces 1
3 0 0 0 0 0 1 1 0 0 0 1 0 1 0.5
Triangulations 0

TShapes 1
Co

1100000
*

+1 0
'''
        locations, tables, shapes, root = read(text)
        self.assertEqual(locations[1][3], 10.0)
        self.assertEqual(locations[0][:3], IDENTITY[:3])
        self.assertEqual(tables['Curve2ds'], ['line'])
        self.assertEqual(tables['Surfaces'], ['ConicalSurface'])
        self.assertEqual(summary(text), ({'ConicalSurface': 1}, []))
        with self.assertRaises(Exception):
            read(text.replace('Surfaces 1', 'Surfaces 2'))


if __name__ == '__main__':
    unittest.main()
