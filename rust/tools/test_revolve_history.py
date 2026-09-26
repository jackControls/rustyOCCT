#!/usr/bin/env python3
"""Deliberate failures for the revolve history bridge: the captured MakeRevol
observations against hand-written kernel histories of two cones, which must
match, and every failure class caught when a relation is wrong."""
import copy
import unittest

from compare_revolve_history import CAPTURE, compare, parse_native, parse_rust
import compare_primitives
import generate_primitive_fixtures

# The kernel's histories of `frustum` and `apex`, as examples/cone_probe.rs
# prints them.
RUST = """case frustum
C 2 3 3 3 1 1
G P0.b0 region S 21.991148575128555 0.9999999999999999 2.0 4.178571428571429
G P0.s0 start_cap F plane 12.566370614359172 0.9999999999999999 2.0 3.0
G P0.s1 wall F cone 29.80376479738831 0.9999999999999999 2.0 4.333333333333334
G P0.s2 end_cap F plane 3.141592653589793 1.0 2.0 6.0
G P0.v1 bottom_edge E circle 3.0 2.0 3.0 -1.0 2.0000000000000004 3.0 3.0 1.9999999999999996 3.0
G P0.v2 top_edge E circle 2.0 2.0 6.0 0.0 2.0 6.0 2.0 1.9999999999999998 6.0
B S 21.991148575128555 0.9999999999999999 2.0 4.178571428571429
end
case apex
C 2 3 2 2 1 1
G P0.b0 region S 12.566370614359181 -1.1485336349731071e-16 -7.067899292142185e-17 0.7499999999999996
G P0.s0 start_cap F plane 12.566370614359172 -7.796343665038751e-17 0.0 0.0
G P0.s1 wall F cone 22.654346798277956 -9.423865722855219e-17 -9.423865722855199e-17 0.9999999999999996
G P0.v1 bottom_edge E circle 2.0 0.0 0.0 -2.0 2.4492935982947064e-16 0.0 2.0 -4.898587196589413e-16 0.0
G P0.v2 apex V 0.0 0.0 3.0
B S 12.566370614359181 -1.1485336349731071e-16 -7.067899292142185e-17 0.7499999999999996
end
"""


class RevolveHistory(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.native = parse_native((CAPTURE/'native.txt').read_text())
        cls.rust = parse_rust(RUST)
        cls.size = compare_primitives.sizes()
        cls.cases = {c.name: c for c in generate_primitive_fixtures.corpus()}

    def run_case(self, name, rust=None, native=None):
        return compare(self.cases[name], native or self.native[name], rust or self.rust[name],
                       self.size[name])[0]

    def test_the_capture_matches_the_kernel_histories(self):
        for name in ('frustum', 'apex'):
            self.assertEqual(self.run_case(name), [], name)

    def test_every_failure_class_is_caught(self):
        def relation(rust, role, change):
            k = next(i for i, r in enumerate(rust['relations']) if r[1] == role)
            parent, role, (kind, numbers) = rust['relations'][k]
            rust['relations'][k] = change(parent, role, kind, list(numbers))

        cases = {
            # A ring from the wrong rim point.
            'unmapped_generated': ('frustum', lambda r: relation(
                r, 'bottom_edge', lambda p, role, k, n: ('P0.v2', role, (k, n)))),
            # An apex elsewhere.
            'unmapped_generated': ('apex', lambda r: relation(
                r, 'apex', lambda p, role, k, n: (p, role, (k, [n[0], n[1], n[2]+1e-3])))),
            # A cap of the wrong area.
            'kernel_relation_unreached': ('frustum', lambda r: relation(
                r, 'end_cap', lambda p, role, k, n: (p, role, (k, [n[0]*(1+1e-6)]+n[1:])))),
            'counts': ('frustum', lambda r: r.update(counts=[3, 3, 3, 3, 1, 1])),
            'body': ('apex', lambda r: r.update(body=('S', [r['body'][1][0]*(1+1e-6)]+r['body'][1][1:]))),
        }
        for kind, (name, change) in cases.items():
            rust = copy.deepcopy(self.rust[name])
            change(rust)
            self.assertIn(kind, self.run_case(name, rust=rust), kind)
        # A missing cap leaves its native disc uncovered.
        rust = copy.deepcopy(self.rust['frustum'])
        rust['relations'] = [r for r in rust['relations'] if r[1] != 'start_cap']
        self.assertIn('native_output_uncovered', self.run_case('frustum', rust=rust))
        # OCCT reporting a radial edge as generating would change the
        # deleted set.
        native = copy.deepcopy(self.native['frustum'])
        native['deleted'] = [d for d in native['deleted'] if d != ('edge', 0)]
        self.assertIn('deleted_set', self.run_case('frustum', native=native))
        # A FirstShape reaching the wall: no rule maps it.
        native = copy.deepcopy(self.native['frustum'])
        wall = next(q[2] for q in native['queries'] if q[1] == 'gen' and q[2][0] == 'F cone')
        k = next(i for i, q in enumerate(native['queries']) if q[1] == 'first' and q[0] == ('vertex', 1))
        native['queries'][k] = (('vertex', 1), 'first', wall)
        self.assertIn('unmapped_first_or_last', self.run_case('frustum', native=native))
        native = copy.deepcopy(self.native['frustum'])
        native['valid'] = False
        self.assertIn('native_invalid', self.run_case('frustum', native=native))


if __name__ == '__main__':
    unittest.main()
