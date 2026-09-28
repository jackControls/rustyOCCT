#!/usr/bin/env python3
"""Deliberate failures the tessellation reference must catch.

A hand-made mesh of the fixture box (two triangles per face, outward) and of
a coarse octahedral sphere passes; each mutation (a flipped triangle, a
missing one, a node moved off the boundary or inward past the request, a
bound understated, a mesh reported over the request) must be reported by
`tessellation_reference.check` with its own failure name. For the spline
bodies of T-b, a hand-made mesh of the dome (its top a grid of points of
the spline graph) passes and its mutations are named, and a spline profile
decides membership and distance around its extreme point.
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

    def test_spline_profile_membership_and_distance(self):
        # The bulge's quadratic from (3, 0) over (4, 1) to (3, 2) reaches
        # x = 3.5 at y = 1, where its tangent is vertical.
        case = next(c for c in fixtures.spline_cases() if c.name == 'spline_bulge')
        path = ref.body(case).shape.paths[0]
        self.assertTrue(path.inside((3.49, 1.0)))
        self.assertFalse(path.inside((3.51, 1.0)))
        self.assertTrue(path.inside((1.0, 1.0)))
        self.assertAlmostEqual(path.distance((3.6, 1.0)), 0.1, delta=1e-12)
        self.assertAlmostEqual(path.distance((3.4, 1.0)), 0.1, delta=1e-12)

    def test_dome_passes_and_every_mutation_is_named(self):
        case, mesh = dome_mesh()
        measured, failures = ref.check(case, 0.2, 1.0, mesh)
        self.assertEqual(failures, [])
        self.assertTrue(0 < measured['deflection'] < 0.2)
        moved = copy.deepcopy(mesh)
        x, y, z = moved.nodes[-1]
        moved.nodes[-1] = (x, y, z+0.3)
        self.assertIn('node_off_boundary', ref.check(case, 0.2, 1.0, moved)[1])
        understated = copy.deepcopy(mesh)
        understated.triangles = [(a, b, c, f, 1e-9 if f == 5 else d, t)
                                 for a, b, c, f, d, t in mesh.triangles]
        self.assertIn('bound_unsound', ref.check(case, 0.2, 1.0, understated)[1])
        flipped = copy.deepcopy(mesh)
        a, b, c, f, d, t = flipped.triangles[-1]
        flipped.triangles[-1] = (a, c, b, f, d, t)
        self.assertIn('misoriented', ref.check(case, 0.2, 1.0, flipped)[1])
        self.assertIn('deflection_exceeded', ref.check(case, measured['deflection']*0.5, 1.0, mesh)[1])


def dome_mesh():
    """The spline dome by hand: its bottom and four walls fanned from a
    corner, its top a 6 x 4 grid of surface points, every triangle
    outward; each bound 0.2."""
    case = next(c for c in fixtures.spline_cases() if c.name == 'spline_dome')
    w, d, h = case.size
    surface = ref.PatchSurface(case.surface)
    nu, nv = 6, 4
    nodes = [(0.0, 0.0, 0.0), (w, 0.0, 0.0), (w, d, 0.0), (0.0, d, 0.0)]
    top = {}
    for j in range(nv+1):
        for i in range(nu+1):
            u, v = w*i/nu, d*j/nv
            p = surface.jet(u, v)[0]
            # The boundary rows are the top edges at height h exactly.
            if i in (0, nu) or j in (0, nv):
                p = (u, v, h)
            top[i, j] = len(nodes)
            nodes.append(p)
    triangles = []

    def fan(face, ring, outward):
        for k in range(1, len(ring)-1):
            a, b, c = ring[0], ring[k], ring[k+1]
            n = ref.cross(ref.sub(nodes[b], nodes[a]), ref.sub(nodes[c], nodes[a]))
            triangles.append((a, b, c, face, 0.2, 0.0) if ref.dot(n, outward) > 0 else (a, c, b, face, 0.2, 0.0))
    fan(0, [0, 1, 2, 3], (0.0, 0.0, -1.0))
    fan(1, [0, 1]+[top[i, 0] for i in range(nu, -1, -1)], (0.0, -1.0, 0.0))
    fan(2, [1, 2]+[top[nu, j] for j in range(nv, -1, -1)], (1.0, 0.0, 0.0))
    fan(3, [2, 3]+[top[i, nv] for i in range(nu+1)], (0.0, 1.0, 0.0))
    fan(4, [3, 0]+[top[0, j] for j in range(nv+1)], (-1.0, 0.0, 0.0))
    for j in range(nv):
        for i in range(nu):
            a, b, c, e = top[i, j], top[i+1, j], top[i+1, j+1], top[i, j+1]
            triangles += [(a, b, c, 5, 0.2, 0.0), (a, c, e, 5, 0.2, 0.0)]
    # The last node, an interior top node, stays last for the mutations.
    last = top[nu-1, nv-1]
    nodes.append(nodes[last])
    triangles = [tuple(len(nodes)-1 if x == last and k < 3 else x for k, x in enumerate(t)) for t in triangles]
    return case, ref.Mesh(nodes, triangles, [], 0.2, 0.0)


if __name__ == '__main__':
    unittest.main()
