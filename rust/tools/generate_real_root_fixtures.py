#!/usr/bin/env python3
"""Exact polynomial root bounds, multiplicities, and algebraic sign queries."""
import argparse
from fractions import Fraction as F
from pathlib import Path
import random
from exact_polynomial_oracle import roots, multiply
from generate_curved_fixtures import enclosure
from generate_spatial_fixtures import bits, value

ROOT = Path(__file__).resolve().parents[1]


def generate():
    tiny, maximum = value(1), value(0x7fefffffffffffff)
    cases = [
        ('zero',[],[1.],None), ('constant',[2.],[-1.],None),
        ('large_positive',[-maximum,tiny],[1.,1.],None),
        ('large_negative',[maximum,tiny],[1.,1.],None),
        ('small',[-tiny,0.,1.],[0.,1.],None),
        ('tiny_coefficients',[-tiny,2*tiny,3*tiny],[0.,1.],None),
        ('negative_lead',[4.,-12.,8.,8.,-11.,1.,3.,-1.],[-2.,0.,1.],None),
    ]
    for n in [9,25]:
        a=3*2**18
        p=[-2.,float(4*a),float(-2*a*a)]+[0.]*(n-3)+[1.]
        cases.append((f'cluster_{n}',p,[-1.,float(a)],(0.,1.)))
    rng=random.Random(0x51_9A_23)
    for i in range(55):
        p=[F(rng.choice([-1,1]))]
        for _ in range(1+i%7):
            factor=[F(rng.randint(-5,5)),F(rng.randint(1,4))] if i%3 else [F(1),F(0),F(1)]
            for _ in range(1+i%3): p=multiply(p,factor)
        if len(p)>26: continue
        q=p if i%4==0 else [F(rng.randint(-4,4)) for _ in range(1+i%26)]
        cases.append((f'factored_{i}',list(map(float,p)),list(map(float,q)),None if i%2 else (-2.,2.)))
    for i in range(30):
        p=[float(rng.randint(-5,5)) for _ in range(2+i%25)]
        q=[float(rng.randint(-5,5)) for _ in range(1+i%26)]
        cases.append((f'dense_{i}',p,q,None))
    rows=[]
    for name,p,q,domain in cases:
        found=roots(p,*(domain or (None,None)))
        record=['A'] if found is None else ['N',str(len(found))]
        for root in found or []:
            record.append(str(root.multiplicity))
            try: record += ['E',*enclosure(root.compare)]
            except OverflowError: record += ['U']
            record.append(str(root.sign_at(q)))
        header=[name,str(len(p)),str(len(q)),*(map(bits,domain) if domain else ['*','*']),*map(bits,p),*map(bits,q)]
        rows.append(' '.join(header+record))
    return '# name np nq lower upper coefficient_bits query_bits status count [multiplicity E lo hi / U; query_sign]\n'+'\n'.join(rows)+'\n'


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--check',action='store_true');args=parser.parse_args()
    text=generate();path=ROOT/'fixtures/real-roots.tsv'
    if args.check:
        if path.read_text()!=text: parser.error('real-root fixtures changed; investigate before updating')
    else:path.write_text(text)
    print(f'real-roots.tsv: {len(text.splitlines())-1} independent exact cases')


if __name__=='__main__':main()
