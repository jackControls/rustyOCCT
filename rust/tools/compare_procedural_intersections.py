#!/usr/bin/env python3
"""Source-pinned GeomInt_IntSS (GeomAPI_IntSS's engine) observations beside the S7b reference.

The independent reference (procedural_intersection_reference.py) classifies
every case of procedural-intersection-cases.txt (the quadric pairs of S7b.1
and S7b.2, a torus with a plane, a sphere or a coaxial surface in S7b.3a)
and gives its canonical rows;
the native probe (occt_procedural_intersection_oracle.cpp) runs
GeomInt_IntSS on the kernel's stored frames and samples its lines. Every
native sample must lie within 1e-6 of the case's size from the reference
curve, and every component of the reference curve (a loop, each ring) must
carry native samples; an empty or single-point reference must be matched by
no line or by a point at it, circles (a torus's special and coaxial pairs)
by samples on them. `--capture` records the native observations
before any kernel code for these pairs exists; later runs must reproduce
them. Differences need a fingerprinted review. The kernel's certified curves
(`procedural_intersection_probe`) must have the reference's class and
components, contain its loop ranges, ends and ring and node points, and put
the middle of each loop within 1e-12 of the reference's.
"""
import argparse
import json
import math
from pathlib import Path
import subprocess
import sys

import mpmath as mp

from build_pinned_occt import SOURCE, digest
from compare_brep import review_for, run, sha, write
from compare_brep_io import build
from compare_degree_elevation import verify_sdk
from compare_occt import ROOT
import analytic_intersection_reference as ana
import generate_procedural_intersection_fixtures as fixtures
import procedural_intersection_reference as ref
from identity_reference import frame_axes

SOURCE_FILE = ROOT/'rust/tools/occt_procedural_intersection_oracle.cpp'
REVIEWS = ROOT/'rust/fixtures/occt-procedural-intersection-divergences.json'
KERNEL_FILE = ROOT/'rust/kernel/src/intersection/procedural.rs'
TORUS_FILE = ROOT/'rust/kernel/src/intersection/toroidal.rs'
BOUND = 1e-6
# Each sub-step's native observations were captured before its kernel code:
# (directory, case-name prefixes, whether its kernel code exists).
CAPTURES = {
    's7b1': (ROOT/'rust/fixtures/occt-procedural-intersection-preimplementation', ('cc_', 'cs_'),
             lambda: KERNEL_FILE.exists()),
    's7b2': (ROOT/'rust/fixtures/occt-procedural-cone-preimplementation', ('ck_', 'ks_'),
             lambda: KERNEL_FILE.exists() and 'Surface::Cone' in KERNEL_FILE.read_text()),
    's7b3a': (ROOT/'rust/fixtures/occt-procedural-torus-preimplementation', ('tp_', 'ts_', 'tx_'),
              lambda: TORUS_FILE.exists()),
}


def native_input(prefixes=None):
    rows = []
    for name, a, b in fixtures.cases():
        if prefixes and not name.startswith(prefixes):
            continue
        rows.append(f'case {name}')
        for kind, values in (a, b):
            o, x, _, n = frame_axes(tuple(values[:9]))
            rows.append(' '.join(['surface', kind, *map(repr, (*o, *n, *x)),
                                  *(repr(float(v)) for v in values[9:])]))
        rows.append('end')
    return '\n'.join(rows)+'\n'


def parse_native(stdout):
    """{case: (status, [[line points]], [points])}."""
    out, current = {}, None
    for line in stdout.splitlines():
        w = line.split()
        if w[0] == 'L':
            current[1].append([])
        elif w[0] == 'p':
            current[1][-1].append([float(x) for x in w[1:4]])
        elif w[0] == 'P':
            current[2].append([float(x) for x in w[1:4]])
        else:
            current = (w[1], [], [])
            out[w[0]] = current
    return out


def surfaces_of(name, a, b):
    _, surfaces = ana.parse(fixtures.encode(name, a, b))
    return surfaces


class FloatTorus:
    """A torus curve of the reference evaluated in binary64 for distances
    (at the 1e-6 bound, 80 digits only cost time): the same frame, `D` and
    points as `procedural_intersection_reference.TorusCurve`."""
    def __init__(self, cv, other):
        self.o, self.a, self.x, self.y = ([float(v) for v in w] for w in (cv.o, cv.a, cv.x, cv.y))
        self.R, self.r = float(cv.R), float(cv.r)
        oo, ao = other.axes()
        o = [float(v) for v in oo]
        if other.kind == 'plane':
            n = [float(v) for v in ao]
            self.f = lambda q: sum((a-b)*c for a, b, c in zip(q, o, n))
        else:
            r2 = float(other.radius)**2
            self.f = lambda q: sum((a-b)**2 for a, b in zip(q, o))-r2

    def meridian(self, phi, t):
        c, s = math.cos(phi), math.sin(phi)
        k = self.R+self.r*math.cos(t)
        return [o+k*(c*x+s*y)+self.r*math.sin(t)*a for o, x, y, a in zip(self.o, self.x, self.y, self.a)]

    def coefficients(self, phi):
        g0, g1, g2 = (self.f(self.meridian(phi, t)) for t in (0.0, math.pi/2, math.pi))
        f0 = (g0+g2)/2
        return f0, (g0-g2)/2, g1-f0

    def D(self, phi):
        f0, alpha, beta = self.coefficients(phi)
        return alpha*alpha+beta*beta-f0*f0

    def point(self, phi, sign):
        f0, alpha, beta = self.coefficients(phi)
        k = max(-1.0, min(1.0, -f0/math.hypot(alpha, beta)))
        return self.meridian(phi, math.atan2(beta, alpha)+sign*math.acos(k))


def float_distance(cv, p, branch=None, span=None):
    """distance_at for a FloatTorus, in binary64."""
    rel = [a-b for a, b in zip(p, cv.o)]
    u = math.atan2(sum(a*b for a, b in zip(rel, cv.y)), sum(a*b for a, b in zip(rel, cv.x)))
    lo, hi = u-0.02, u+0.02
    if span is not None:
        u0, u1 = float(span[0]), float(span[1])
        u -= 2*math.pi*math.floor((u-u0)/(2*math.pi))
        lo, hi = max(u-0.02, u0), min(u+0.02, u1)
        if lo > hi:
            return min(math.dist(p, cv.point(t, 1)) for t in (u0, u1))

    def dist(t, sign):
        return math.dist(p, cv.point(t, sign)) if cv.D(t) >= 0 else math.inf
    best = math.inf
    ts = [lo+(hi-lo)*k/200 for k in range(201)]
    # A loop's end inside the window, by bisection on the sign of D: the
    # curve is vertical there, so the grid alone misses its nearest points.
    for t0, t1 in zip(ts, ts[1:]):
        if (cv.D(t0) >= 0) != (cv.D(t1) >= 0):
            a, b = (t0, t1) if cv.D(t0) >= 0 else (t1, t0)
            for _ in range(60):
                m = (a+b)/2
                a, b = (m, b) if cv.D(m) >= 0 else (a, m)
            best = min(best, dist(a, 1))
    for sign in ([branch] if branch else [1, -1]):
        ds = [dist(t, sign) for t in ts]
        k = min(range(len(ts)), key=lambda j: ds[j])
        a, b = ts[max(k-1, 0)], ts[min(k+1, 200)]
        g = (math.sqrt(5)-1)/2
        for _ in range(60):
            c1, c2 = b-g*(b-a), a+g*(b-a)
            if dist(c1, sign) < dist(c2, sign):
                b = c2
            else:
                a = c1
        best = min(best, ds[k], dist((a+b)/2, sign))
    return best


def distance(cv, p, branch=None, span=None):
    """Distance from a point to the reference curve (or one branch, or one
    loop's parameter range): the point's angle on the ruled surface, then
    the nearer branch there."""
    if isinstance(cv, FloatTorus):
        return float_distance(cv, p, branch, span)
    # 24 digits: distances are compared with 1e-6, and 80 digits cost time.
    with mp.workdps(24):
        return mp_distance(cv, p, branch, span)


def mp_distance(cv, p, branch, span):
    pm = [mp.mpf(x) for x in p]
    if isinstance(cv, ref.ConeCurve):
        # From the apex a point of the lower nappe (v < 0) lies along
        # -d(u): its angle about the axis is u + pi. Try both.
        rel = [x-y for x, y in zip(pm, cv.V)]
        u = mp.atan2(sum(a*b for a, b in zip(rel, cv.y)), sum(a*b for a, b in zip(rel, cv.x)))
        return min(distance_at(cv, pm, t, branch, span) for t in (u, u+mp.pi))
    rel = [x-y for x, y in zip(pm, cv.o)]
    u = mp.atan2(sum(a*b for a, b in zip(rel, cv.y)), sum(a*b for a, b in zip(rel, cv.x)))
    return distance_at(cv, pm, u, branch, span)


def distance_at(cv, pm, u, branch, span):
    """The Euclidean distance from a point to the curve near parameter u:
    near a loop's end the branches are vertical in (u, v), so the point's
    own angle is not its nearest parameter; search a neighbourhood of u on
    both branches (within a loop's range when given) and refine by golden
    section."""
    lo, hi = u-mp.mpf('0.02'), u+mp.mpf('0.02')
    if span is not None:
        u0, u1 = span
        k = mp.floor((u-u0)/(2*mp.pi))
        u = u-2*mp.pi*k
        lo, hi = max(u-mp.mpf('0.02'), u0), min(u+mp.mpf('0.02'), u1)
        if lo > hi:
            ends = [cv.point(u0, 1), cv.point(u1, 1)]
            return float(min(mp.sqrt(sum((a-b)**2 for a, b in zip(pm, q))) for q in ends))

    def dist(t, sign):
        if cv.D(t) < 0:
            return mp.inf
        q = cv.point(t, sign)
        return mp.sqrt(sum((a-b)**2 for a, b in zip(pm, q)))
    best = mp.inf
    for sign in ([branch] if branch else [1, -1]):
        ts = [lo+(hi-lo)*k/200 for k in range(201)]
        ds = [dist(t, sign) for t in ts]
        k = min(range(len(ts)), key=lambda j: ds[j])
        a, b = ts[max(k-1, 0)], ts[min(k+1, 200)]
        g = (mp.sqrt(5)-1)/2
        for _ in range(60):
            c1, c2 = b-g*(b-a), a+g*(b-a)
            if dist(c1, sign) < dist(c2, sign):
                b = c2
            else:
                a = c1
        best = min(best, ds[k], dist((a+b)/2, sign))
    return float(best)


def circle_distance(row, p):
    """The distance from a point to a reference circle row."""
    c, n, r = [mp.mpf(x) for x in row[1]], [mp.mpf(x) for x in row[2]], mp.mpf(row[3])
    w = [mp.mpf(a)-b for a, b in zip(p, c)]
    h = sum(a*b for a, b in zip(w, n))
    radial = mp.sqrt(max(mp.mpf(0), sum(a*a for a in w)-h*h))
    return float(mp.sqrt((radial-r)**2+h*h))


def surface_distance(s, p):
    """The distance from a point to a plane, sphere or torus of a case."""
    o, a = s.axes()
    rel = [mp.mpf(x)-ana.mpf(y) for x, y in zip(p, o)]
    am = [ana.mpf(x) for x in a]
    la = mp.sqrt(sum(x*x for x in am))
    h = sum(x*y for x, y in zip(rel, am))/la
    if s.kind == 'plane':
        return float(abs(h))
    if s.kind == 'sphere':
        return float(abs(mp.sqrt(sum(x*x for x in rel))-ana.mpf(s.radius)))
    rho = mp.sqrt(max(mp.mpf(0), sum(x*x for x in rel)-h*h))
    return float(abs(mp.sqrt((rho-ana.mpf(s.radius))**2+h*h)-ana.mpf(s.minor)))


def differences(name, surfaces, native, rows):
    """What separates the native result from the reference's rows."""
    status, lines, points = native
    if status != 'done':
        return ['not_done']
    samples = [p for line in lines for p in line]
    scale = max([1.0]+[abs(x) for s in surfaces for x in s.frame[:3]]+[s.radius for s in surfaces])
    tol = BOUND*scale
    kind = rows[0] if isinstance(rows[0], str) else rows[0][0]
    if kind in ('empty', 'same'):
        return [] if not samples and not points else ['spurious']
    if kind == 'not_conic':
        # Not parameterised yet (a sphere containing a meridian circle):
        # native samples must lie on both surfaces.
        if not samples:
            return ['missed_curve']
        return [] if all(surface_distance(s, p) <= tol for s in surfaces for p in samples+points) \
            else ['off_surface']
    if kind == 'point':
        targets = [[float(x) for x in row[1]] for row in rows]
        near = lambda t: [p for p in samples+points if max(abs(a-b) for a, b in zip(p, t)) <= 1e-4*scale]
        if not all(near(t) for t in targets):
            return ['missed_tangent_point']
        return [] if sum(len(near(t)) for t in targets) == len(samples+points) else ['spurious']
    if kind == 'circle':
        out = []
        d = [[circle_distance(row, p) for row in rows] for p in samples+points]
        if any(min(x) > tol for x in d):
            out.append('off_curve')
        if {min(range(len(rows)), key=x.__getitem__) for x in d if min(x) <= tol} != set(range(len(rows))):
            out.append('missed_component')
        return out
    _, cv = ref.curve(*surfaces)
    if isinstance(cv, ref.TorusCurve):
        cv = FloatTorus(cv, next(s for s in surfaces if s.kind != 'torus'))
    # Components: each loop by its range, each ring by its branch, a
    # figure-eight whole.
    if kind == 'loop':
        comps = [dict(span=(r[1], r[2])) for r in rows]
    elif kind == 'rings':
        comps = [dict(branch=1), dict(branch=-1)]
    else:
        comps = [dict()]
    out = []
    if any(distance(cv, p) > tol for p in samples+points):
        out.append('off_curve')
    covered = set()
    for p in samples:
        d = [distance(cv, p, **c) for c in comps]
        k = min(range(len(d)), key=d.__getitem__)
        if d[k] <= tol:
            covered.add(k)
    if covered != set(range(len(comps))):
        out.append('missed_component')
    return out


def rust_rows():
    """{case: [(kind, [(lo, hi)])]} from the kernel's probe, one per row."""
    subprocess.run(['cargo', '+stable', 'build', '--release', '--locked', '--example',
                    'procedural_intersection_probe'], cwd=ROOT, check=True)
    text = (ROOT/'rust/fixtures/procedural-intersection-cases.txt').read_text()
    rows = subprocess.run([str(ROOT/'target/release/examples/procedural_intersection_probe')], input=text,
                          text=True, capture_output=True, timeout=600, check=True).stdout
    out = {}
    for line in rows.splitlines():
        w = line.split()
        v = [float(x) for x in w[2:]] if w[1] != 'error' else []
        out.setdefault(w[0], []).append((w[1], list(zip(v[::2], v[1::2]))))
    return out


def rust_differences(rust, expected):
    """The kernel against the reference, row by row: its kinds, every exact
    reference number inside the kernel's enclosure (up to 1e-25 relative,
    the printing), a loop's middle within 1e-12."""
    want = [e if isinstance(e, str) else e[0] for e in expected]
    if [k for k, _ in rust] != want:
        return ['rust_class']
    for (kind, got), row in zip(rust, expected):
        if isinstance(row, str):
            continue
        # The reference's numbers as it prints them (exact zeros as zero).
        numbers = []
        for part in row[1:]:
            numbers += [float(ref.number(x)) for x in (part if isinstance(part, (list, tuple)) else [part])]
        if len(numbers) != len(got):
            return ['rust_values']
        inside = lambda k: got[k][0]-1e-25*abs(numbers[k]) <= numbers[k] <= got[k][1]+1e-25*abs(numbers[k])
        exact = range(len(numbers)) if kind != 'loop' else range(8)
        if not all(inside(k) for k in exact):
            return ['rust_outside_reference']
        if kind == 'loop':
            scale = max([1.0]+[abs(x) for x in numbers])
            if any(abs((got[k][0]+got[k][1])/2-numbers[k]) > 1e-12*scale for k in range(8, 14)):
                return ['rust_loop_middle']
    return []


def capture(executable, env, key, sdk_manifest):
    CAPTURE, prefixes, exists = CAPTURES[key]
    text = native_input(prefixes)
    record = run(executable, text, env)
    if record['exit_code'] != 0:
        raise SystemExit('native procedural intersection run failed: '+json.dumps(record)[:2000])
    CAPTURE.mkdir(parents=True, exist_ok=True)
    (CAPTURE/'inputs.txt').write_text(text)
    (CAPTURE/'native.txt').write_text(record['stdout'])
    (CAPTURE/'oracle.cpp').write_text(SOURCE_FILE.read_text())
    status = subprocess.run(['git', 'status', '--porcelain', '--', 'rust'], cwd=ROOT, text=True,
                            capture_output=True, check=True).stdout.splitlines()
    revision = subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True, capture_output=True,
                              check=True).stdout.strip()
    write(CAPTURE/'capture.json', {
        'source_reference': SOURCE, 'oracle': next(iter(record['stderr'].splitlines()), None),
        'platform': sys.platform, 'rust_revision': revision,
        'rust_procedural_intersection_exists': exists(),
        'rust_worktree_uncommitted': status, 'sdk_manifest_sha256': digest(sdk_manifest),
        'input_sha256': digest(CAPTURE/'inputs.txt'), 'probe_source_sha256': digest(CAPTURE/'oracle.cpp'),
        'observations_sha256': digest(CAPTURE/'native.txt')})


def captured(observed):
    """True when every capture was taken on this platform and is reproduced
    exactly. On another platform IntPatch's walking lines differ in their
    last digits (and near degeneracies in their pieces), so there the
    captures' integrity and cases are checked and the native observations
    are held to the reference instead."""
    return all([captured_one(observed, CAPTURE, prefixes) for CAPTURE, prefixes, _ in CAPTURES.values()])


def captured_one(observed, CAPTURE, prefixes):
    metadata = json.loads((CAPTURE/'capture.json').read_text())
    if metadata['source_reference'] != SOURCE or metadata['rust_procedural_intersection_exists']:
        raise ValueError('procedural intersection capture was not a clean pre-implementation reference')
    for key, name in [('input_sha256', 'inputs.txt'), ('probe_source_sha256', 'oracle.cpp'),
                      ('observations_sha256', 'native.txt')]:
        if metadata[key] != digest(CAPTURE/name):
            raise ValueError('procedural intersection evidence changed: '+name)
    if (CAPTURE/'inputs.txt').read_text() != native_input(prefixes):
        raise ValueError('the native inputs differ from the captured ones')
    was = parse_native((CAPTURE/'native.txt').read_text())
    if set(was) != {n for n in observed if n.startswith(prefixes)}:
        raise ValueError('native cases differ from the capture')
    if metadata['platform'] != sys.platform:
        return False
    for name, (status, lines, points) in was.items():
        now = observed[name]
        flat = lambda x: [v for p in [q for line in x[1] for q in line]+x[2] for v in p]
        if now[0] != status or [len(l) for l in now[1]] != [len(l) for l in lines] or len(now[2]) != len(points) \
                or any(abs(a-b) > 1e-9*max(1.0, abs(a)) for a, b in zip(flat(now), flat((status, lines, points)))):
            raise ValueError('native observation of '+name+' differs from the capture')
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--occt-root', type=Path, required=True)
    parser.add_argument('--sdk-manifest', type=Path, required=True)
    parser.add_argument('--output', type=Path, default=ROOT/'target/procedural-intersection-oracle')
    parser.add_argument('--strict-native', action='store_true')
    parser.add_argument('--capture', choices=sorted(CAPTURES),
                        help='record one sub-step\'s native observations (before its implementation only)')
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    prefix = args.occt_root.resolve()
    verify_sdk(prefix, args.sdk_manifest)
    for name, text in fixtures.generate().items():
        if text != (ROOT/'rust/fixtures'/name).read_text():
            raise ValueError('independent fixture regeneration changed: '+name)
    executable, env, loaded, command = build(prefix, output, SOURCE_FILE, 'procedural-oracle', ('TKGeomAlgo',))
    text = native_input()
    if args.capture:
        capture(executable, env, args.capture, args.sdk_manifest)
        print('captured', args.capture)
        return
    record = run(executable, text, env)
    if record['exit_code'] != 0:
        raise SystemExit('native procedural intersection run failed: '+json.dumps(record)[:2000])
    observed = parse_native(record['stdout'])
    # Off the capture's platform a review holds for the same case and kinds
    # of difference, its fingerprint being the capture platform's output.
    strict = captured(observed)
    oracle = next(iter(record['stderr'].splitlines()), None)
    reviews = [] if args.strict_native or not REVIEWS.exists() else json.loads(REVIEWS.read_text())['reviews']
    report = {'source_reference': SOURCE, 'oracle': oracle, 'capture_reproduced': strict,
              'cases': 0, 'native_samples': 0,
              'rust_within_reference': 0, 'matches': [], 'reviewed_differences': [], 'failures': []}
    rust = rust_rows()
    for name, a, b in fixtures.cases():
        report['cases'] += 1
        rows = ref.rows(*surfaces_of(name, a, b))
        wrong = rust_differences(rust[name], rows)
        if wrong:
            report['failures'].append({'case': name, 'reason': ' '.join(wrong), 'rust': rust[name]})
            continue
        report['rust_within_reference'] += 1
        native = observed[name]
        report['native_samples'] += sum(len(l) for l in native[1])+len(native[2])
        found = differences(name, surfaces_of(name, a, b), native, rows)
        if not found:
            report['matches'].append(name)
            continue
        evidence = {'case': name, 'source_reference': SOURCE, 'oracle': oracle,
                    'native_sha256': sha(json.dumps(native)), 'differences': found}
        review = review_for(evidence if strict else {k: v for k, v in evidence.items() if k != 'native_sha256'},
                            reviews)
        report['reviewed_differences' if review else 'failures'].append(
            review or dict(evidence, lines=len(native[1]), points=len(native[2])))
    write(output/'capture.json', {'source_reference': SOURCE, 'oracle': oracle,
                                  'sdk_manifest_sha256': digest(args.sdk_manifest),
                                  'source_sha256': digest(SOURCE_FILE), 'probe_sha256': digest(executable),
                                  'loaded_libraries': loaded, 'build_command': command})
    write(output/'report.json', report)
    print(json.dumps({k: len(v) if isinstance(v, list) else v for k, v in report.items()}, indent=2))
    if report['failures']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
