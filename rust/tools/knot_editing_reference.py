"""Independent exact knot editing oracle.

Recover controls by solving EVERY homogeneous Cox-basis power coefficient
equation on the common span partition. Inconsistent equations prove that
removal is impossible. A separate rational Greville collocation reconstruction
is available for checking selected cases. No Boehm/de Boor insertion or
inverse-insertion recurrence is used here.
"""
from fractions import Fraction as F
from functools import lru_cache
from bisect import bisect_right
from generate_spline_fixtures import axis
from bezier_editing_reference import add, times_linear


def parse(row):
    w = row.split(); p,n,k,op = map(int,w[2:6]); a = 6
    data = list(map(lambda x:F(float(x)),w[a:a+4*n])); a += 4*n
    controls = tuple(tuple(data[i+c]*data[i+3] if c<3 else data[i+3] for c in range(4)) for i in range(0,4*n,4))
    knots = tuple(F(float(x)) for x in w[a:a+2*k:2]); mults = tuple(map(int,w[a+1:a+2*k:2])); a += 2*k
    operations = [(w[a+3*i],F(float(w[a+3*i+1])),int(w[a+3*i+2])) for i in range(op)]
    if len(w)!=a+3*op: raise ValueError('input token count')
    return w[0], (p,w[1]=='P',knots,mults,controls), operations


def data(curve):
    p,periodic,knots,mults,controls = curve
    flat,a,b = axis(p,knots,mults,periodic)
    return p,periodic,flat,a,b,len(controls)


@lru_cache(maxsize=2048)
def values(p,periodic,knots,mults,u):
    flat,a,b = axis(p,knots,mults,periodic)
    n = sum(mults)-(mults[0] if periodic else p+1)
    if periodic: u = a+(u-a)%(b-a)
    left = not periodic and u==flat[-1]
    row = [F(x<u<=y if left else x<=u<y) for x,y in zip(flat,flat[1:])]
    for d in range(1,p+1):
        row = [((u-flat[i])*row[i]/(flat[i+d]-flat[i]) if flat[i+d]!=flat[i] else F(0))
               +((flat[i+d+1]-u)*row[i+1]/(flat[i+d+1]-flat[i+1]) if flat[i+d+1]!=flat[i+1] else F(0))
               for i in range(len(row)-1)]
    result = [F(0)]*n
    for i,v in enumerate(row): result[i%n] += v
    return tuple(result)


def value(curve,u):
    p,periodic,knots,mults,controls = curve
    row = values(p,periodic,knots,mults,u)
    return [sum(a*c[j] for a,c in zip(row,controls)) for j in range(4)]


@lru_cache(maxsize=256)
def inverse(p,periodic,knots,mults):
    flat,a,b = axis(p,knots,mults,periodic)
    n = sum(mults)-(mults[0] if periodic else p+1)
    nodes = tuple(sum(flat[i+1:i+p+1])/p for i in range(n))
    # Exact Gauss-Jordan on sparse rows. Store inverse once per candidate basis;
    # input weights and control coordinates never affect this factorization.
    rows = []
    for i,u in enumerate(nodes):
        row = {j:x for j,x in enumerate(values(p,periodic,knots,mults,u)) if x}
        row[n+i] = F(1); rows.append(row)
    for j in range(n):
        pivot = next((i for i in range(j,n) if rows[i].get(j,0)),None)
        if pivot is None: raise ArithmeticError('singular independent collocation basis')
        rows[j],rows[pivot] = rows[pivot],rows[j]
        scale = rows[j][j]; rows[j] = {i:x/scale for i,x in rows[j].items()}
        for k in range(n):
            if k==j or not rows[k].get(j,0): continue
            scale = rows[k][j]
            for i,x in rows[j].items():
                y = rows[k].get(i,F(0))-scale*x
                if y: rows[k][i] = y
                else: rows[k].pop(i,None)
    return nodes,tuple(tuple((j-n,x) for j,x in row.items() if j>=n) for row in rows)


def recover_by_collocation(original,knots,mults):
    p,periodic = original[:2]
    nodes,matrix = inverse(p,periodic,knots,mults)
    samples = [value(original,u) for u in nodes]
    controls = tuple(tuple(sum(x*samples[j][c] for j,x in row) for c in range(4)) for row in matrix)
    return p,periodic,knots,mults,controls


@lru_cache(maxsize=512)
def basis_polynomials(p,periodic,knots,mults,lo,hi):
    flat,a,b = axis(p,knots,mults,periodic)
    n = sum(mults)-(mults[0] if periodic else p+1)
    if periodic:
        offset = ((lo-a)//(b-a))*(b-a)
        lo,hi = lo-offset,hi-offset
    # The query is a cell of the common knot partition. Exterior spans of
    # unclamped curves are included; finite basis functions vanish off support.
    midpoint = (lo+hi)/2
    span = bisect_right(flat,midpoint)-1
    row = [[F(1)] if i==span else [] for i in range(len(flat)-1)]
    for d in range(1,p+1):
        result = []
        for i in range(len(row)-1):
            left,right = flat[i+d]-flat[i],flat[i+d+1]-flat[i+1]
            x = times_linear(row[i],(lo-flat[i])/left,(hi-lo)/left) if left and row[i] else []
            y = times_linear(row[i+1],(flat[i+d+1]-lo)/right,(lo-hi)/right) if right and row[i+1] else []
            coefficients = add(x,y)
            while coefficients and not coefficients[-1]: coefficients.pop()
            result.append(coefficients)
        row = result
    result = [[F(0)]*(p+1) for _ in range(n)]
    for i,coefficients in enumerate(row):
        for k,x in enumerate(coefficients): result[i%n][k] += x
    return tuple(map(tuple,result))


def polynomial(curve,lo,hi):
    p,periodic,knots,mults,controls = curve
    matrix = basis_polynomials(p,periodic,knots,mults,lo,hi)
    return tuple(tuple(sum(row[k]*c[j] for row,c in zip(matrix,controls) if row[k] and c[j]) for k in range(p+1)) for j in range(4))


def partition(original,candidate):
    p,periodic,flat,a,b,n = data(original)
    cuts = set()
    if periodic:
        for curve in [original,candidate]:
            cuts.update(a+(k-a)%(b-a) for k in curve[2])
        cuts.add(b)
    else:
        cuts.update(original[2]); cuts.update(candidate[2])
    cuts = sorted(cuts)
    return list(zip(cuts,cuts[1:]))


def equal(original,candidate):
    return all(polynomial(original,lo,hi)==polynomial(candidate,lo,hi) for lo,hi in partition(original,candidate))


def recover(original,knots,mults):
    p,periodic = original[:2]
    n = sum(mults)-(mults[0] if periodic else p+1)
    candidate = p,periodic,knots,mults,((F(0),)*4,)*n
    pivots = {}
    for lo,hi in partition(original,candidate):
        matrix = basis_polynomials(p,periodic,knots,mults,lo,hi)
        rhs = polynomial(original,lo,hi)
        for k in range(p+1):
            row = {i:coeffs[k] for i,coeffs in enumerate(matrix) if coeffs[k]}
            values = [rhs[c][k] for c in range(4)]
            while row:
                j = min(row)
                if j not in pivots:
                    factor = row[j]
                    pivots[j] = {i:x/factor for i,x in row.items()},[x/factor for x in values]
                    break
                basis,answer = pivots[j]; factor = row[j]
                for i,x in basis.items():
                    y = row.get(i,F(0))-factor*x
                    if y: row[i] = y
                    else: row.pop(i,None)
                values = [v-factor*x for v,x in zip(values,answer)]
            else:
                if any(values): return None
    if len(pivots)!=n: raise ArithmeticError('incomplete independent coefficient system')
    controls = [None]*n
    for j in reversed(range(n)):
        row,answer = pivots[j]
        controls[j] = tuple(answer[c]-sum(x*controls[i][c] for i,x in row.items() if i!=j) for c in range(4))
    return p,periodic,knots,mults,tuple(controls)


@lru_cache(maxsize=1024)
def edited(curve,op,u,target):
    p,periodic,knots,mults,controls = curve
    _,_,_,a,b,n = data(curve)
    if periodic and u==b: u=a
    mapping = dict(zip(knots,mults)); old = mapping.get(u,0)
    if (op=='I' and target<=old) or (op=='R' and target>=old): return True,curve
    if op=='I':
        mapping[u] = target
        if periodic and u==a: mapping[b] = target
    elif op=='R':
        if not old: raise ValueError('missing removal knot')
        if target:
            mapping[u] = target
            if periodic and u==a: mapping[b] = target
        else:
            del mapping[u]
            if periodic and u==a:
                del mapping[b]
                start = min(mapping)
                mapping[start+b-a] = mapping[start]
    else: raise ValueError('unknown operation')
    knots,mults = map(tuple,zip(*sorted(mapping.items())))
    candidate = recover(curve,knots,mults)
    valid = candidate is not None and all(c[3]>0 for c in candidate[4])
    if op=='I' and not valid: raise AssertionError('independent refinement reconstruction failed')
    return (True,candidate) if valid else (False,curve)


def expected(row):
    name,curve,operations = parse(row); flags = []
    for op,u,m in operations:
        flag,curve = edited(curve,op,u,m); flags.append(flag)
    return flags,curve


def encode(name,flags,curve):
    p,periodic,knots,mults,controls = curve
    _,_,_,a,b,n = data(curve)
    w = [name,'R',str(len(flags)),*[str(int(f)) for f in flags],str(p),str(int(periodic)),str(n),str(len(knots)),str(a),str(b)]
    w += [str(x) for c in controls for x in c]
    w += [s for k,m in zip(knots,mults) for s in [str(k),str(m)]]
    return ' '.join(w)


def decode(text,native=False):
    import math
    result = {}
    for row in text.splitlines():
        w = iter(row.split()); name = next(w)
        if name in result or next(w)!='R': raise ValueError('duplicate or unsuccessful observation')
        flag_values = [int(next(w)) for _ in range(int(next(w)))]
        if any(f not in [0,1] for f in flag_values): raise ValueError('invalid operation flag')
        flags = list(map(bool,flag_values))
        p,periodic,n,k = (int(next(w)) for _ in range(4))
        number = (lambda:float(next(w))) if native else (lambda:F(next(w)))
        a,b = number(),number()
        controls = tuple(tuple(number() for _ in range(4)) for _ in range(n))
        knots = []; mults = []
        for _ in range(k): knots.append(number()); mults.append(int(next(w)))
        if list(w) or not 1<=p<=25 or periodic not in [0,1] or n<=p or n>4096 or not a<b:
            raise ValueError('invalid observation')
        if any(c[3]<=0 for c in controls): raise ValueError('nonpositive control weight')
        if native and not all(math.isfinite(x) for x in [a,b,*knots,*[v for c in controls for v in c]]):
            raise ValueError('nonfinite native observation')
        curve = p,bool(periodic),tuple(knots),tuple(mults),controls
        if not native and (a,b)!=data(curve)[3:5]: raise ValueError('wrong parameter domain')
        result[name] = flags,curve,(a,b)
    return result
