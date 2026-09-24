"""Independent certification of Rust spline/linear rows and native comparison.

Expected parameters come from generate_spline_linear_fixtures.py (Cox, QQ
factorization, VAS, rational sample-point clipping). A Rust row passes only if
it has exactly the independent point/interval counts, and every tight binary64
parameter bound contains its irreducible root. Line-parameter bounds must meet
the independent rigorous bracket. Native results are compared afterwards and
never certify completeness.
"""
from fractions import Fraction as F
import math

from exact_polynomial_oracle import Root, poly
from generate_spline_linear_fixtures import make_cases, preimage

ABSOLUTE = 1e-6
RELATIVE = 1e-10
NATIVE_CASES = 31


def native_line(case):
    words = [case['name'], case['kind'], case['degree'], int(case['periodic']),
             len(case['controls']), len(case['knots']), *case['range'],
             *case['endpoints'][0], *case['endpoints'][1],
             *(x for p in case['controls'] for x in p),
             *(x for pair in zip(case['knots'], case['multiplicities']) for x in pair)]
    return ' '.join(w if isinstance(w, str) else str(w) if isinstance(w, int)
                    else format(float(w), '.17g') for w in words)


def native_cases():
    return make_cases()[:NATIVE_CASES]


def native_inputs():
    return '\n'.join(map(native_line, native_cases()))+'\n'


def expected(case):
    """Independent points and intervals as (first, last, P, lo, hi, s_lo, s_hi)."""
    def parse(words):
        first, last, n = words[0], words[1], int(words[2])
        p = words[3:3+n]
        lo, hi, s_lo, s_hi = words[3+n:7+n]
        return dict(first=F(first), last=F(last), p=[F(x) for x in p],
                    lo=F(lo), hi=F(hi), s=(F(s_lo), F(s_hi)))
    points, intervals = preimage(case)
    return ([parse(p) for p in points], [[parse(a), parse(b)] for a, b in intervals])


class Words:
    def __init__(self, text):
        self.words = text.split()
        self.i = 0

    def take(self, kind=str):
        if self.i >= len(self.words):
            raise ValueError('truncated row')
        self.i += 1
        return kind(self.words[self.i-1])

    def end(self):
        if self.i != len(self.words):
            raise ValueError('trailing row data')


def tight(lo, hi):
    if not (math.isfinite(lo) and math.isfinite(hi) and (lo == hi or math.nextafter(lo, math.inf) == hi)):
        raise ValueError('bound is not a tight finite enclosure')


def decode_rust(line):
    w = Words(line)
    name = w.take()
    if w.take() != 'R':
        raise ValueError('Rust error row')
    np, ni = w.take(int), w.take(int)

    def parameter():
        values = [w.take(float) for _ in range(10)]
        return dict(u=values[0:2], s=values[2:4], xyz=[values[4:6], values[6:8], values[8:10]])
    points = []
    for _ in range(np):
        if w.take() != 'P':
            raise ValueError('expected Rust point')
        points.append(parameter())
    intervals = []
    for _ in range(ni):
        if w.take() != 'I':
            raise ValueError('expected Rust interval')
        intervals.append([parameter(), parameter()])
    w.end()
    return name, dict(points=points, intervals=intervals)


def certify(observed, want):
    """observed Rust bounds contain the independent algebraic parameter."""
    lo, hi = observed['u']
    tight(lo, hi)
    root = Root(poly(want['p']), 1, (want['lo'], want['hi']))
    width = want['last']-want['first']
    local = lambda x: (F(x)-want['first'])/width
    if root.compare(local(lo)) < 0 or root.compare(local(hi)) > 0:
        raise ValueError('Rust parameter bound misses the independent root')
    s_lo, s_hi = observed['s']
    tight(s_lo, s_hi)
    if F(s_lo) > want['s'][1] or F(s_hi) < want['s'][0]:
        raise ValueError('Rust line parameter misses the independent bracket')


def verify_rust(case, line):
    name, actual = decode_rust(line)
    if name != case['name']:
        raise ValueError('incorrect Rust case identity')
    points, intervals = expected(case)
    if len(actual['points']) != len(points) or len(actual['intervals']) != len(intervals):
        raise ValueError('incomplete or extra Rust preimage components')
    for observed, want in zip(actual['points'], points):
        certify(observed, want)
    for observed, want in zip(actual['intervals'], intervals):
        for o, w in zip(observed, want):
            certify(o, w)
    return actual


def near(value, bounds):
    a, b = bounds
    budget = ABSOLUTE+RELATIVE*max(abs(a), abs(b))
    return a-budget <= value <= b+budget


def decode_native(text, case, offset=0.):
    """Common parts of one native query row, shifted by an integer period offset.
    Line parameters are arc lengths from A along a unit direction; they are
    converted to A + s(B-A) units."""
    rows = text.splitlines()
    if len(rows) != 1:
        raise ValueError('incorrect native row count')
    w = Words(rows[0])
    if w.take() != case['name']:
        raise ValueError('wrong native identity')
    if w.take() != 'R' or w.take(int) != 1:
        raise ValueError('native exception or incomplete intersection')
    a, b = [[float(x) for x in p] for p in case['endpoints']]
    length = math.dist(a, b)
    points, edges = [], []
    for _ in range(w.take(int)):
        tag, curve_first = w.take(), w.take(int)
        if tag == 'P':
            u, v = w.take(float), w.take(float)
            first, second = [w.take(float) for _ in range(3)], [w.take(float) for _ in range(3)]
            if not curve_first:
                u, v, first, second = v, u, second, first
            points.append(dict(u=u+offset, s=v/length, curve=first, line=second))
        elif tag == 'I':
            r = [w.take(float), w.take(float)]
            lines = [[w.take(float), w.take(float)] for _ in range(w.take(int))]
            if not curve_first:
                raise ValueError('curve range reported second')
            edges.append(dict(u=[r[0]+offset, r[1]+offset], s=[[x/length for x in l] for l in lines]))
        else:
            raise ValueError('unknown native common part')
    w.end()
    return dict(points=points, edges=edges)


def merge_native(parts):
    """Combine adapted windows: identical-within-budget vertices and edges that
    touch within budget are one observation."""
    points, edges = [], []
    for part in parts:
        for p in part['points']:
            if not any(near(p['u'], (q['u'], q['u'])) and near(p['s'], (q['s'], q['s'])) for q in points):
                points.append(p)
        edges.extend(part['edges'])
    edges.sort(key=lambda e: e['u'][0])
    merged = []
    for e in edges:
        if merged and near(e['u'][0], (merged[-1]['u'][1], merged[-1]['u'][1])):
            merged[-1] = dict(u=[merged[-1]['u'][0], max(merged[-1]['u'][1], e['u'][1])],
                              s=merged[-1]['s']+e['s'])
        else:
            merged.append(dict(e))
    return dict(points=sorted(points, key=lambda p: p['u']), edges=merged)


def compare_native(actual, native):
    """Sorted, deduplicated difference labels; empty means a complete match."""
    differences = set()
    points, edges = native['points'], native['edges']
    # Maximum bipartite matching of exact points to native vertices.
    options = [[j for j, p in enumerate(points)
                if near(p['u'], want['u']) and near(p['s'], want['s'])]
               for want in actual['points']]
    assigned = {}

    def match(i, seen):
        for j in options[i]:
            if j not in seen:
                seen.add(j)
                if j not in assigned or match(assigned[j], seen):
                    assigned[j] = i
                    return True
        return False
    for i in range(len(options)):
        if not match(i, set()):
            differences.add('missing_isolated_point')
    intervals = [[x['u'][0], y['u'][1]] for x, y in actual['intervals']]
    for j, p in enumerate(points):
        if j in assigned:
            continue
        if any(near(p['u'], iv) for iv in intervals):
            differences.add('overlap_reported_as_vertex')
        else:
            differences.add('extra_native_vertex')
    used = set()
    for iv in intervals:
        exact = [k for k, e in enumerate(edges) if near(e['u'][0], (iv[0], iv[0])) and near(e['u'][1], (iv[1], iv[1]))]
        if exact:
            used.add(exact[0])
            continue
        overlapping = [k for k, e in enumerate(edges)
                       if e['u'][0] <= iv[1]+ABSOLUTE and e['u'][1] >= iv[0]-ABSOLUTE]
        if overlapping:
            used.update(overlapping)
            differences.add('overlap_range')
        elif not any(near(p['u'], iv) for p in points):
            differences.add('missing_overlap')
    for k, e in enumerate(edges):
        if k not in used:
            differences.add('extra_native_edge')
    return sorted(differences)
