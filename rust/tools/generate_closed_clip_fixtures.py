"""Independent closed clipping sets from SymPy inequality solving over QQ."""
from pathlib import Path
import hashlib
import json
import random
import sympy as s

BASE=Path(__file__).resolve().parents[1]/'fixtures'
X=s.Symbol('x',real=True)
rng=random.Random(0xC11_2206)
factors=[X+2,X+1,X,X-1,X-2,X*X-2,X*X-3,X*X+1]
cases=[]
for i in range(64):
    def polynomial():
        p=s.Integer(rng.choice([-1,1]))
        for _ in range(1+i%5):
            p*=rng.choice(factors)**rng.randrange(1,4)
        return s.expand(p)
    cases.append((f'factored_{i}',polynomial(),polynomial(),-3,3))
cases.append(('degree25', (X*X-2)**6*(X+1)**7*(X-2)**6, 9-X*X, -3, 3))
cases.extend([
    ('all',s.Integer(0),s.Integer(0),-3,3),
    ('none',s.Integer(0),s.Integer(-1),-3,3),
    ('touch',-X**2,s.Integer(1),-3,3),
    ('irrational_touch',-(X**2-2)**2,s.Integer(1),-3,3),
    ('singleton_hit',-X**2,s.Integer(1),0,0),
    ('singleton_miss',-X**2,s.Integer(1),1,1),
])
for bits in (70,180,600):
    e=s.Rational(1,2**bits)
    cases.append((f'gap_{bits}',(X-1)*(X-1-e),s.Integer(0),0,2))
    cases.append((f'interval_{bits}',-(X-1)*(X-1-e),s.Integer(0),0,2))
for bits in (1024,2048):
    huge=s.Integer(2)**bits
    cases.append((f'huge_{bits}',-(X-huge)*(X-huge-1),s.Integer(0),huge-1,huge+2))


def equation(expr):
    p=s.Poly(expr,X,domain=s.QQ)
    coefficients=list(reversed(p.all_coeffs()))
    return [str(len(coefficients)),*map(str,coefficients)]


def endpoint(expr):
    if expr.is_Rational:
        return [*equation(X-expr),str(expr),str(expr)]
    polynomial=s.Poly(s.minpoly(expr,X),X,domain=s.QQ)
    selected=[]
    for (lo,hi),_ in polynomial.intervals(eps=s.Rational(1,2**20)):
        if bool(expr>lo) and bool(expr<hi):selected.append((lo,hi))
    assert len(selected)==1,(expr,selected)
    lo,hi=selected[0]
    return [*equation(polynomial.as_expr()),str(lo),str(hi)]


rows=[];count_points=0;count_intervals=0;max_degree=0
for name,p,q,lo,hi in cases:
    allowed=s.Interval(lo,hi)
    for e in (p,q):
        allowed=allowed.intersect(s.solve_univariate_inequality(e>=0,X,relational=False))
    parts=list(allowed.args) if isinstance(allowed,s.Union) else [allowed]
    points=[];intervals=[]
    for part in parts:
        if part==s.S.EmptySet:continue
        if isinstance(part,s.FiniteSet):points.extend(part)
        else:
            assert isinstance(part,s.Interval) and not part.left_open and not part.right_open
            intervals.append((part.start,part.end))
    points.sort();intervals.sort(key=lambda v:v[0])
    row=[name,*equation(p),*equation(q),str(lo),str(hi),str(len(points)),str(len(intervals))]
    for point in points:row+=endpoint(point)
    for a,b in intervals:row+=endpoint(a)+endpoint(b)
    rows.append(' '.join(row));count_points+=len(points);count_intervals+=len(intervals)
    max_degree=max(max_degree,s.Poly(p,X).degree(),s.Poly(q,X).degree())
path=BASE/'closed-clip.tsv';content=('# name [n coefficients]x2 lo hi points intervals; each endpoint has n coefficients lo hi\n'+'\n'.join(rows)+'\n')
import argparse
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--check',action='store_true')
args=parser.parse_args()
if args.check:
    if path.read_text()!=content:
        parser.error('independent fixture changed; investigate before updating')
else:
    path.write_text(content)
report=dict(cases=len(cases),isolated_points=count_points,maximal_closed_intervals=count_intervals,
            max_equation_degree=int(max_degree),oracle='SymPy exact univariate inequality intersection; minimal polynomial and certified real isolator at every endpoint',
            fixture_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),)
print(json.dumps(report))
