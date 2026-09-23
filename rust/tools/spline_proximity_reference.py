"""Certify the complete Rust minimum set and all reported binary64 views.

Independent Cox equations, QQ irreducible factors, VAS refinement and resultant
images supply the expected answer. Native comparisons run only after these
checks; an OCCT observation is never a mathematical expected value.
"""
from fractions import Fraction as F
import math

from exact_polynomial_oracle import Root, poly
from generate_spline_proximity_fixtures import coefficients, equations, minima


class Words:
    def __init__(self, line):
        self.words = iter(line.split())

    def take(self, convert=str):
        return convert(next(self.words))

    def end(self):
        if next(self.words, None) is not None:
            raise ValueError('trailing output')


def finite(word):
    value = float(word)
    if not math.isfinite(value):
        raise ValueError('non-finite output')
    return value


def bounds(words):
    a, b = words.take(finite), words.take(finite)
    if a > b or (a != b and math.nextafter(a, math.inf) != b):
        raise ValueError('output is not a tight finite binary64 enclosure')
    return a, b


def certify(view, compare):
    a, b = map(F, view)
    if a == b:
        if compare(a) != 0:
            raise ValueError('incorrect exact binary64 value')
    elif compare(a) != 1 or compare(b) != -1:
        raise ValueError('incorrect binary64 enclosure')


def decode_rust(line):
    words = Words(line)
    name, status = words.take(), words.take()
    if status != 'R':
        raise ValueError('Rust did not return a complete minimum set')
    np, ni = words.take(int), words.take(int)
    if not 0 <= np <= 65536 or not 0 <= ni <= 4096:
        raise ValueError('invalid minimum counts')
    distance = bounds(words)
    points, intervals = [], []
    for _ in range(np):
        if words.take() != 'P':
            raise ValueError('missing point marker')
        points.append({'parameter': bounds(words), 'xyz': [bounds(words) for _ in range(3)]})
    for _ in range(ni):
        if words.take() != 'I':
            raise ValueError('missing interval marker')
        a, b = words.take(F), words.take(F)
        if a >= b:
            raise ValueError('invalid minimum interval')
        intervals.append([a, b])
    words.end()
    return name, {'distance': distance, 'points': points, 'intervals': intervals}


def verify_rust(case, line):
    name, actual = decode_rust(line)
    if name != case['name']:
        raise ValueError('incorrect Rust case identity')
    cells = equations(case)
    points, spans = minima(cells)
    if len(actual['points']) != len(points) or actual['intervals'] != spans:
        raise ValueError('incomplete or incorrect minimum parameter set')
    query = list(map(F, case['query']))
    for (a, b, p, lo, hi), observed in zip(points, actual['points']):
        root = Root(poly(p), 1, (lo, hi))
        certify(observed['parameter'], lambda x: root.compare((x-a)/(b-a)))
        cell = next(c for c in cells if root.compare((c['lo']-a)/(b-a)) >= 0
                    and root.compare((c['hi']-a)/(b-a)) <= 0)
        if cell['lo'] == cell['hi']:
            h = cell['h']
        else:
            change = poly([(a-cell['lo'])/(cell['hi']-cell['lo']), (b-a)/(cell['hi']-cell['lo'])])
            h = [v.compose(change) for v in cell['h']]
        for coordinate, view in zip(h[:3], observed['xyz']):
            certify(view, lambda x: root.sign_at(coefficients(coordinate-x*h[3])))
        n = sum(((h[i]-query[i]*h[3])**2 for i in range(3)), poly([]))
        # Verify the shared minimum distance against EVERY returned minimizer.
        minimum_sign = lambda x: root.sign_at(coefficients(n-x*h[3]*h[3]))
        certify(actual['distance'], minimum_sign)
    if spans:
        cell = next(c for c in cells if c['lo'] == spans[0][0])
        distance = F(cell['n'].eval(0)/(cell['w'].eval(0)**2))
        certify(actual['distance'], lambda x: (distance > x)-(distance < x))
    if not points and not spans:
        raise ValueError('empty global minimum on a nonempty compact domain')
    actual['domain'] = list(map(F, case['range']))
    return actual


def flag(words):
    value = words.take(int)
    if value not in (0, 1):
        raise ValueError('invalid native flag')
    return bool(value)


def decode_native(family, text, first, last):
    """Keep all stationary candidates. Legacy receives its reported endpoints;
    the newer API already has an explicit PerformWithEndpoints variant.
    """
    rows = text.splitlines()
    if len(rows) != (1 if family == 'legacy' else 2):
        raise ValueError('incorrect native row count')
    results = []
    names = []
    for mode, row in enumerate(rows):
        w = Words(row)
        names.append(w.take())
        if w.take() != 'R':
            raise ValueError('native exception or invalid output')
        if family == 'legacy':
            done, infinite = flag(w), False
        else:
            if w.take(int) != mode:
                raise ValueError('incorrect native variant')
            status = w.take(int)
            if status not in range(6):
                raise ValueError('unknown native status')
            done, infinite = flag(w), flag(w)
            # OCCT considers both OK and NoSolution completed finite searches.
            if done != (status in (0, 3)) or infinite != (status == 2):
                raise ValueError('inconsistent native status flags')
        count = w.take(int)
        if not 0 <= count <= 100000:
            raise ValueError('invalid native candidate count')
        points = []
        for _ in range(count):
            parameter = w.take(finite)
            minimum = flag(w)
            if family != 'legacy':
                flag(w)
            d = w.take(finite)
            xyz = [w.take(finite) for _ in range(3)]
            if d < 0:
                raise ValueError('negative squared distance')
            points.append({'parameter': parameter, 'xyz': xyz, 'distance': d, 'minimum': minimum})
        if family == 'legacy':
            if w.take() != 'E':
                raise ValueError('missing legacy endpoints')
            da, db = w.take(finite), w.take(finite)
            if da < 0 or db < 0:
                raise ValueError('negative endpoint squared distance')
            a, b = [w.take(finite) for _ in range(3)], [w.take(finite) for _ in range(3)]
            points += [{'parameter': first, 'xyz': a, 'distance': da, 'minimum': False},
                       {'parameter': last, 'xyz': b, 'distance': db, 'minimum': False}]
        infinite_distance = None
        if infinite:
            if w.take() != 'I':
                raise ValueError('missing infinite-distance value')
            infinite_distance = w.take(finite)
            if infinite_distance < 0:
                raise ValueError('negative infinite squared distance')
        if count:
            if w.take() != 'M':
                raise ValueError('missing selected native minimum')
            if family == 'legacy':
                w.take(finite)
                w.take(finite)
            else:
                index = w.take(int)
                d = w.take(finite)
                if not 0 <= index < count or d != points[index]['distance']:
                    raise ValueError('incorrect native selected minimum index or value')
        w.end()
        results.append({'done': done, 'infinite': infinite, 'infinite_distance': infinite_distance,
                        'points': points, 'stationary_candidates': count})
    if len(set(names)) != 1:
        raise ValueError('native variant names differ')
    return names[0], results[-1]


# Native APIs return floating approximations. These fixed budgets do not decide
# Rust identity, tie equality, or which candidates are mathematical minima.
ABSOLUTE = 1e-6
RELATIVE = 1e-10


def nearby(value, exact_bounds):
    a, b = exact_bounds
    budget = ABSOLUTE+RELATIVE*max(abs(a), abs(b))
    return a-budget <= value <= b+budget


def compare_native(actual, native):
    differences = []
    if actual['intervals']:
        if actual['intervals'] != [actual['domain']]:
            differences.append('minimum_interval_ranges_not_represented')
        if not native['infinite']:
            differences.append('whole_minimum_interval_not_represented')
        elif not nearby(native['infinite_distance'], actual['distance']):
            differences.append('infinite_squared_distance')
    elif native['infinite']:
        differences.append('unexpected_infinite_minima')
    candidates = native['points']
    if candidates:
        if not nearby(min(p['distance'] for p in candidates), actual['distance']):
            differences.append('global_squared_distance')
    elif not native['infinite']:
        differences.append('no_native_candidates')
    # Exact duplicate observations cannot cover two distinct parameter aliases.
    unique = { (p['parameter'], *p['xyz']): p for p in candidates }
    candidates = list(unique.values())
    edges = [[j for j, p in enumerate(candidates)
              if nearby(p['parameter'], want['parameter'])
              and nearby(p['distance'], actual['distance'])
              and all(nearby(x,b) for x,b in zip(p['xyz'],want['xyz']))]
             for want in actual['points']]
    assigned = {}

    def match(i, visited):
        for j in edges[i]:
            if j in visited:
                continue
            visited.add(j)
            if j not in assigned or match(assigned[j], visited):
                assigned[j] = i
                return True
        return False

    covered = sum(match(i, set()) for i in range(len(edges)))
    if covered != len(edges):
        differences.append(f'minimum_witness_coverage_{covered}_of_{len(edges)}')
    return differences
