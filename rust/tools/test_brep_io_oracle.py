#!/usr/bin/env python3
"""Deliberate failures for the .brep interop bridge and the independent reader."""
import json
from pathlib import Path
import tempfile
import unittest

from brep_io_reference import IDENTITY, read, summary
from compare_brep_io import (ADAPTIVE_SOURCE, CAPTURE, FREE_SOURCE, PROPERTY_BOUND, REVIEWS, V3_LINE,
                             close, measure_differences, original_capture, parse, review_measure,
                             v3_capture, v3_files)


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

    def test_v3_capture_is_a_clean_pre_implementation_record(self):
        solids, free, metadata = v3_capture()
        self.assertEqual(len(metadata['inputs_sha256']), 26)
        self.assertEqual(sum(map(len, solids.values())), 11)
        self.assertEqual(sum(map(len, free.values())), 28)
        self.assertTrue(all(v.startswith('refused line ') for v in metadata['rust_refusals'].values()))
        self.assertIs(metadata['rust_v3_header_read'], False)

    def test_v3_files_match_the_whole_header_line(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            (root/'sub').mkdir()
            texts = {'a.brep': b'DBRep_DrawableShape\n\n'+V3_LINE+b'\nLocations 0\n',
                     'sub/b.brep': V3_LINE+b'\r\nLocations 0\r\n',
                     'c.brep': b'DBRep_DrawableShape\n\nCASCADE Topology V1, (c) Matra-Datavision\n',
                     'd.brep': b'Open CASCADE Topology V3 (c)\n\x00\x01',
                     'e.brep': V3_LINE+b' and more\n',
                     'f.txt': V3_LINE+b'\n'}
            for name, text in texts.items():
                (root/name).write_bytes(text)
            self.assertEqual([str(p.relative_to(root)) for p in v3_files(root)], ['a.brep', 'sub/b.brep'])

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

    def test_adaptive_probe_changes_only_the_integration(self):
        # The S6 probe stays pinned; the adaptive one differs in its header,
        # its banner and SurfaceProperties' Eps alone.
        body = lambda text: text[text.index('#include'):]
        adaptive = body(ADAPTIVE_SOURCE.read_text())
        self.assertIn('BRepGProp::SurfaceProperties(s, g, 1e-12);', adaptive)
        restored = adaptive.replace('SurfaceProperties(s, g, 1e-12)', 'SurfaceProperties(s, g)').replace(
            ' BRepGProp adaptive 1e-12"', ' BRepGProp"')
        self.assertEqual(restored, body(FREE_SOURCE.read_text()))

    def test_measure_reviews_need_an_independent_measure_inside(self):
        enclosure = [(100.0, 100.0+1e-10), (1.0, 1.0+1e-12), (2.0, 2.0), (-3.0-1e-12, -3.0)]
        default, adaptive = [100.001, 1.0, 2.0, -3.0], [100.0, 1.0, 2.0, -3.0]
        self.assertEqual(measure_differences(enclosure, default, adaptive), ['default_integration'])
        self.assertEqual(measure_differences(enclosure, default, default), ['adaptive_outside:measure'])
        base = dict(case='a.brep-free-1', source_reference=None, oracle='o', input_sha256=None,
                    native_sha256=None, differences=None, reason='r', independent_evidence='e',
                    independent_measure=['100.00000000005', '1.0000000000005', '2', '-3.0000000000005'])
        args = ('a.brep-free-1', 'text', 'o', enclosure, default, adaptive)
        # The fingerprint comes from the run; a review matches only it.
        _, failure = review_measure(*args, [])
        self.assertEqual(failure['reason'], 'unreviewed measure difference')
        review = dict(base, **{k: failure[k] for k in ['source_reference', 'input_sha256', 'native_sha256',
                                                        'differences']})
        self.assertEqual(review_measure(*args, [review]), (review, None))
        for changed, reason in [({'independent_measure': None}, 'the review records no independent measure'),
                                ({'independent_measure': ['100.0000000002', '1', '2', '-3']},
                                 'independent measure outside the enclosure: measure'),
                                ({'independent_measure': ['100', '1', '2.0000000000001', '-3']},
                                 'independent measure outside the enclosure: cy')]:
            _, failure = review_measure(*args, [dict(review, **changed)])
            self.assertEqual(failure['reason'], reason)
        # Another default row, file or adaptive verdict is unreviewed.
        for other in [('a.brep-free-1', 'text', 'o', enclosure, [100.002]+default[1:], adaptive),
                      ('a.brep-free-1', 'other', 'o', enclosure, default, adaptive),
                      ('a.brep-free-1', 'text', 'o', enclosure, default, default)]:
            self.assertEqual(review_measure(*other, [review])[1]['reason'], 'unreviewed measure difference')

    def test_committed_measure_reviews_are_complete(self):
        for review in json.loads(REVIEWS.read_text())['reviews']:
            self.assertTrue(review['reason'] and review['independent_evidence'])
            if review.get('native_sha256'):
                self.assertEqual(len(review['independent_measure']), 4)
                self.assertTrue(all(isinstance(x, str) for x in review['independent_measure']))

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
        self.assertEqual(tables['Surfaces'], ['cone'])
        self.assertEqual(summary(text), ({}, []))
        with self.assertRaises(Exception):
            read(text.replace('Surfaces 1', 'Surfaces 2'))
        # Version 3 (OCCT 7.6 on): its own line, and a normals flag and
        # normals in each triangulation record.
        v1 = 'CASCADE Topology V1, (c) Matra-Datavision'
        mesh = text.replace('\nTriangulations 0', '\nTriangulations 1\n3 1 0 0.1\n0 0 0 1 0 0 0 1 0\n1 2 3')
        v3 = mesh.replace(v1, 'CASCADE Topology V3, (c) Open Cascade').replace(
            '3 1 0 0.1', '3 1 0 1 0.1').replace('1 2 3\n', '1 2 3\n0 0 1 0 0 1 0 0 1\n')
        self.assertEqual(read(v3), read(mesh))
        self.assertEqual(read(text.replace(v1, 'CASCADE Topology V3, (c) Matra-Datavision')), read(text))
        for bad in ['CASCADE Topology V1, (c) Open Cascade', 'CASCADE Topology V3, (c) Open Cascad']:
            with self.assertRaises(AssertionError):
                read(text.replace(v1, bad))
        with self.assertRaises(Exception):
            read(v3.replace('3 1 0 1 0.1', '3 1 0 0.1'))


if __name__ == '__main__':
    unittest.main()
