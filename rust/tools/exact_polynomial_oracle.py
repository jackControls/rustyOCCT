"""Independent test oracle: SymPy factorization and continued-fraction isolation.

Production uses subresultant integer Sturm/Sturm-Tarski sequences. This oracle
uses irreducible QQ factors, VAS isolation/refinement, and rational interval
evaluation. No Rust observation or approximate root supplies an expected value.
Install the hash-pinned math-oracle-requirements.txt in a test-only environment.
"""
from fractions import Fraction as F
from functools import cmp_to_key
import sympy as sp

X = sp.Symbol('x')


def trim(p):
    p = list(map(F, p))
    while p and not p[-1]: p.pop()
    return p


def add(a, b):
    result = [F(0)] * max(len(a), len(b))
    for i, x in enumerate(a): result[i] += x
    for i, x in enumerate(b): result[i] += x
    return trim(result)


def multiply(a, b):
    if not a or not b: return []
    result = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): result[i+j] += x*y
    return trim(result)


def evaluate(p, x):
    value = F(0)
    for c in reversed(p): value = value*x + c
    return value


def sign(x): return (x > 0) - (x < 0)


def poly(coefficients):
    return sp.Poly.from_list(list(reversed(trim(coefficients))) or [0], X, domain=sp.QQ)


class Root:
    def __init__(self, factor, multiplicity, interval):
        self.factor = factor
        self.coefficients = list(map(F, reversed(factor.all_coeffs())))
        self.multiplicity = multiplicity
        self.low, self.high = map(F, interval)

    def refine(self):
        if self.low != self.high:
            before = self.low, self.high
            self.low, self.high = map(F, self.factor.refine_root(self.low, self.high, eps=(self.high-self.low)/16, fast=True))
            if before == (self.low,self.high):
                raise ArithmeticError(f'no refinement: {self.factor}, {before}')

    def sign_at(self, coefficients):
        # Irreducibility makes zero detection an exact polynomial remainder
        # check. Nonzero remainders eventually separate from zero by continuity.
        remainder = trim(coefficients)
        divisor = self.coefficients
        while len(remainder) >= len(divisor):
            shift = len(remainder)-len(divisor)
            scale = remainder[-1]/divisor[-1]
            for i, c in enumerate(divisor): remainder[i+shift] -= scale*c
            remainder = trim(remainder)
        if not remainder: return 0
        for _ in range(4096):
            lo = hi = F(0)
            for c in reversed(remainder):
                products = [lo*self.low, lo*self.high, hi*self.low, hi*self.high]
                lo, hi = min(products)+c, max(products)+c
            if lo > 0: return 1
            if hi < 0: return -1
            if self.low == self.high: return sign(lo)
            self.refine()
        raise ArithmeticError(f'independent root sign oracle exhausted refinement: {self.factor}; {self.low} {self.high}; remainder {remainder}')

    def compare(self, value): return self.sign_at([-F(value), F(1)])


def roots(coefficients, lower=None, upper=None):
    p = poly(coefficients)
    if p.is_zero: return None
    result = []
    for factor, multiplicity in p.factor_list()[1]:
        for interval, _ in factor.intervals():
            root = Root(factor, multiplicity, interval)
            if lower is not None and root.compare(lower) < 0: continue
            if upper is not None and root.compare(upper) > 0: continue
            result.append(root)

    def compare(a, b):
        for _ in range(4096):
            if a.high < b.low: return -1
            if b.high < a.low: return 1
            a.refine(); b.refine()
        raise ArithmeticError(f'distinct factors did not separate: {a.factor}, {(a.low,a.high)}; {b.factor}, {(b.low,b.high)}')
    return sorted(result, key=cmp_to_key(compare))
