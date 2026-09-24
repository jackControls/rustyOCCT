"""Independent prototype oracle: transform equations first, factor over QQ,
then identify/order global irreducible roots by canonical indices and VAS.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import random
import sys
import sympy as s

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tools'))
from generate_root_comparison_fixtures import roots, compare

X = s.Symbol('x')
BASE = Path(__file__).resolve().parents[1] / 'fixtures'
rng = random.Random(0xAF_2026)


def side(global_equation, offset, width, domain=(-3, 3), extra=False):
    offset, width = map(s.Rational, (offset, width))
    local = s.Poly(global_equation.subs(X, offset + width * X), X, domain=s.QQ)
    if extra:
        local *= s.Poly((X * X + 1) ** 2, X, domain=s.QQ)
    lo, hi = [(s.Rational(v) - offset) / width for v in domain]
    # The expected equation is expanded in global parameter units and factored
    # independently; the production comparator composes only a sign query.
    global_poly = s.Poly(local.as_expr().subs(X, (X - offset) / width), X, domain=s.QQ)
    expected = roots([F(global_poly.nth(i)) for i in range(global_poly.degree()+1)], tuple(map(F, domain)))
    return [str(local.degree()+1), str(lo), str(hi), str(offset), str(width),
            *(str(local.nth(i)) for i in range(local.degree()+1)), str(len(expected))], expected


equations = [X**2-2, X**3-3, (X+2)*(X-1)*(X-3), X**5-X-1]
pairs = []
for i in range(48):
    g = equations[i % len(equations)]
    h = g if i % 2 else equations[(i+1) % len(equations)]
    pairs.append((f'random_{i}',
                  side(g, s.Rational(rng.randrange(-20,21),7), s.Rational(rng.randrange(1,10),3), extra=i%3==0),
                  side(h, s.Rational(rng.randrange(-20,21),5), s.Rational(rng.randrange(1,10),11), extra=i%4==0)))
for bits in (70, 180, 600):
    epsilon = s.Rational(1, 2**bits)
    for sign in (-1, 1):
        pairs.append((f'sub_float_{bits}_{sign}', side(X**2-2, 1, 2),
                      side((X-sign*epsilon)**2-2, -3, s.Rational(1,7))))
for bits in (1024, 2048):
    huge = s.Integer(2)**bits
    for sign in (-1, 1):
        for width in (huge, 1/huge):
            # Equal global values whose local parameter representations differ
            # by values beyond binary64 range.
            pairs.append((f'extreme_{bits}_{sign}_{width>1}',
                          side(X**2-2, sign*huge, width), side(X**2-2, 0, 1, extra=True)))
for endpoint in (-3, -1, 0, 1, 3):
    pairs.append((f'rational_{endpoint}', side(X-endpoint, 9, 2),
                  side((X-endpoint)*(X**2+1), -5, s.Rational(1,3))))
for degree in (9,25):
    pairs.append((f'chebyshev_{degree}', side(s.chebyshevt(degree,X), 0, 1),
                  side(s.chebyshevt(degree,X), s.Rational(-7,3), s.Rational(2,5))))

rows=[];counts={-1:0,0:0,1:0}
for name,(aw,ar),(bw,br) in pairs:
    orderings=[compare(a,b) for a in ar for b in br]
    for value in orderings: counts[value]+=1
    rows.append(' '.join([name,*aw,*bw,*map(str,orderings)]))
path=BASE/'affine-parameters.tsv'
content=('# name [n lower upper offset width coefficients root_count]x2 comparisons\n'+'\n'.join(rows)+'\n')
import argparse
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--check',action='store_true')
args=parser.parse_args()
if args.check:
    if path.read_text()!=content:
        parser.error('independent fixture changed; investigate before updating')
else:
    path.write_text(content)
report=dict(pairs=len(pairs),comparisons=sum(counts.values()),orderings=counts,
            oracle='QQ irreducible factor identity and canonical root index; VAS interval separation',
            fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),)

print(json.dumps(report))
