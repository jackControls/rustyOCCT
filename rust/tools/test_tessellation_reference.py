#!/usr/bin/env python3
"""Deliberate failures the tessellation reference must catch.

A hand-made mesh of the fixture box (two triangles per face, outward) and of
a coarse octahedral sphere passes; each mutation (a flipped triangle, a
missing one, a node moved off the boundary or inward past the request, a
bound understated, a mesh reported over the request) must be reported by
`tessellation_reference.check` with its own failure name.
"""
import copy
import math
import unittest

import generate_tessellation_fixtures as fixtures
import tessellation_reference as ref


def box_mesh():
    case = next(c for c in fixtures.cases() if c.name == 'box')
    (x0, y0, z0), (w, d, h) = case.box
    nodes = [(x0+i*w, y0+j*d, z0+k*h) for k in (0, 1) for j in (0, 1) for i in (0, 1)]
    quads = [(0, 2, 3, 1), (4, 5, 7, 6), (0, 1, 5, 4), (2, 6, 7, 3), (0, 4, 6, 2), (1, 3, 7, 5)]
    triangles = []
    for face, (a, b, c, e) in enumerate(quads):
        triangles += [(a, b, c, face, 1e-12, 0.0), (a, c, e, face, 1e-12, 0.0)]
    edges = []
    return case, ref.Mesh(nodes, triangles, edges, 1e-12, 0.0)


def octahedron(case):
    R = case.sphere[0]
    nodes = [(R, 0.0, 0.0), (-R, 0.0, 0.0), (0.0, R, 0.0), (0.0, -R, 0.0), (0.0, 0.0, R), (0.0, 0.0, -R)]
    faces = [(0, 2, 4), (2, 1, 4), (1, 3, 4), (3, 0, 4), (2, 0, 5), (1, 2, 5), (3, 1, 5), (0, 3, 5)]
    # The deviation of a face's centre: R - R / sqrt(3).
    bound = R*(1-1/math.sqrt(3))*1.000001
    return ref.Mesh(nodes, [(a, b, c, 0, bound, 1.6) for a, b, c in faces], [], bound, 1.6)


class ReferenceCatchesMutations(unittest.TestCase):
    def test_box_passes_and_every_mutation_is_named(self):
        case, mesh = box_mesh()
        _, failures = ref.check(case, 0.4, 0.5, mesh)
        self.assertEqual(failures, [])
        flipped = copy.deepcopy(mesh)
        a, b, c, f, d, t = flipped.triangles[0]
        flipped.triangles[0] = (a, c, b, f, d, t)
        self.assertIn('misoriented', ref.check(case, 0.4, 0.5, flipped)[1])
        missing = copy.deepcopy(mesh)
        missing.triangles.pop()
        failures = ref.check(case, 0.4, 0.5, missing)[1]
        self.assertIn('not_closed', failures)
        self.assertIn('euler', failures)
        moved = copy.deepcopy(mesh)
        x, y, z = moved.nodes[7]
        moved.nodes[7] = (x+1.0, y, z)
        self.assertIn('node_off_boundary', ref.check(case, 0.4, 0.5, moved)[1])
        over = copy.deepcopy(mesh)
        over.deflection = 0.5
        self.assertIn('reported_over_request', ref.check(case, 0.4, 0.5, over)[1])

    def test_sphere_deflection_and_bounds(self):
        case = next(c for c in fixtures.cases() if c.name == 'sphere')
        mesh = octahedron(case)
        R = case.sphere[0]
        deviation = R*(1-1/math.sqrt(3))
        measured, failures = ref.check(case, deviation*1.01, 1.6, mesh)
        self.assertEqual(failures, [])
        # The samples miss the centre, the deepest point: within 10% of it.
        self.assertTrue(0.9*deviation <= measured['deflection'] <= deviation)
        self.assertIn('deflection_exceeded', ref.check(case, deviation*0.9, 1.6, mesh)[1])
        understated = copy.deepcopy(mesh)
        understated.triangles = [(a, b, c, f, d*0.5, t) for a, b, c, f, d, t in mesh.triangles]
        self.assertIn('bound_unsound', ref.check(case, deviation*1.01, 1.6, understated)[1])
        inward = copy.deepcopy(mesh)
        inward.triangles = [(a, c, b, f, d, t) for a, b, c, f, d, t in mesh.triangles]
        failures = ref.check(case, deviation*1.01, 1.6, inward)[1]
        self.assertIn('normal_inward', failures)
        self.assertIn('volume', failures)


if __name__ == '__main__':
    unittest.main()
