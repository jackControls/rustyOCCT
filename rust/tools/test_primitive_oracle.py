#!/usr/bin/env python3
"""Deliberate failures for the primitive bridge: every difference class is caught."""
import copy
from pathlib import Path
import unittest

from compare_primitives import CAPTURE, differences, parse_expected, parse_observations, sizes
import generate_primitive_fixtures


class PrimitiveOracle(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.expected = parse_expected(generate_primitive_fixtures.generate()['primitive-expected.tsv'])
        cls.native = parse_observations((CAPTURE/'native.txt').read_text())
        cls.size = sizes()

    def test_the_capture_matches_the_reference(self):
        for name, want in self.expected.items():
            self.assertEqual(differences(self.native[name], want, self.size[name]), [], name)

    def test_every_difference_class_is_caught(self):
        name = 'frustum'
        base = self.native[name]
        cases = {
            'verdict': lambda o: o.update(verdict='invalid'),
            'counts': lambda o: o.update(counts=(3,)+o['counts'][1:]),
            'mass_properties': lambda o: o['props'].__setitem__(0, o['props'][0]*(1+1e-7)),
            'faces': lambda o: o['faces'].pop(),
            'edges': lambda o: o['edges'][0].__setitem__(3, [x+1e-3 for x in o['edges'][0][3]]),
            'vertices': lambda o: o['vertices'].append([0.0, 0.0, 0.0]),
        }
        for kind, change in cases.items():
            changed = copy.deepcopy(base)
            if kind == 'edges':
                changed['edges'] = [list(e) for e in changed['edges']]
            change(changed)
            self.assertIn(kind, differences(changed, self.expected[name], self.size[name]), kind)
        # A face of the wrong type with the right area is still a difference.
        changed = copy.deepcopy(base)
        changed['faces'] = [('cylinder',)+f[1:] if f[0] == 'cone' else f for f in changed['faces']]
        self.assertIn('faces', differences(changed, self.expected[name], self.size[name]))


if __name__ == '__main__':
    unittest.main()
