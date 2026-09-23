"""Complete tensor knot oracle, independent of local knot insertion/removal.

Each new basis is reconstructed from ALL Cox power-coefficient equations on
its full raw support. Factorization depends only on axes and is shared across
transverse controls. Residual equations prove non-removability; no sample
points, de Boor recurrence, production output or rounded values are used.
"""
from fractions import Fraction as F
from functools import lru_cache
import math
from generate_spline_fixtures import axis as expand_axis
from knot_editing_reference import basis_polynomials, partition


def count(axis):
    p,periodic,knots,mults=axis
    return sum(mults)-(mults[0] if periodic else p+1)


def domain(axis):
    p,periodic,knots,mults=axis
    return expand_axis(p,knots,mults,periodic)[1:]


def parse(row):
    w=iter(row.split()); name=next(w)
    du,dv,nu,nv,ku,kv,pu,pv,nop=(int(next(w)) for _ in range(9))
    controls=[]
    for _ in range(nu*nv):
        x,y,z,q=(F(float(next(w))) for _ in range(4))
        controls.append((x*q,y*q,z*q,q))
    axes=[]
    for p,periodic,k in [(du,pu,ku),(dv,pv,kv)]:
        pairs=[(F(float(next(w))),int(next(w))) for _ in range(k)]
        knots,mults=map(tuple,zip(*pairs)); axes.append((p,bool(periodic),knots,mults))
    ops=[(next(w),'UV'.index(next(w)),F(float(next(w))),int(next(w))) for _ in range(nop)]
    if list(w) or [count(a) for a in axes]!=[nu,nv]: raise ValueError('input dimensions')
    return name,(tuple(axes),tuple(controls)),ops


def add_scaled(to,source,scale):
    for i,x in source.items():
        y=to.get(i,F(0))+scale*x
        if y: to[i]=y
        else: to.pop(i,None)


@lru_cache(maxsize=1024)
def reconstruction(old,new):
    """New controls as linear forms in old controls, plus ALL constraints."""
    assert old[:2]==new[:2]
    pivots={}; constraints=[]
    original=(*old,((F(0),)*4,)*count(old))
    candidate=(*new,((F(0),)*4,)*count(new))
    for lo,hi in partition(original,candidate):
        before=basis_polynomials(*old,lo,hi); after=basis_polynomials(*new,lo,hi)
        for k in range(old[0]+1):
            row={i:c[k] for i,c in enumerate(after) if c[k]}
            rhs={i:c[k] for i,c in enumerate(before) if c[k]}
            while row:
                j=min(row)
                if j not in pivots:
                    scale=row[j]
                    pivots[j]=({i:x/scale for i,x in row.items()},{i:x/scale for i,x in rhs.items()})
                    break
                basis,answer=pivots[j]; scale=-row[j]
                add_scaled(row,basis,scale); add_scaled(rhs,answer,scale)
            else:
                if rhs: constraints.append(tuple(sorted(rhs.items())))
    if len(pivots)!=count(new): raise ArithmeticError('incomplete coefficient rank')
    controls=[None]*count(new)
    for j in reversed(range(count(new))):
        row,answer=pivots[j]; answer=answer.copy()
        for i,x in row.items():
            if i!=j: add_scaled(answer,controls[i],-x)
        controls[j]=answer
    return tuple(tuple(sorted(r.items())) for r in controls),tuple(constraints)


@lru_cache(maxsize=256)
def edit(surface,op,which,u,target):
    axes,grid=surface; old=axes[which]; p,periodic,knots,mults=old
    a,b=domain(old)
    if periodic and u==b: u=a
    mapping=dict(zip(knots,mults)); previous=mapping.get(u,0)
    if (op=='I' and target<=previous) or (op=='R' and target>=previous): return True,surface
    if op=='I':
        mapping[u]=target
        if periodic and u==a: mapping[b]=target
    elif op=='R':
        if not previous: raise ValueError('missing knot')
        if target:
            mapping[u]=target
            if periodic and u==a: mapping[b]=target
        else:
            del mapping[u]
            if periodic and u==a:
                del mapping[b]; start=min(mapping); mapping[start+b-a]=mapping[start]
    else: raise ValueError('unknown operation')
    nk,nm=map(tuple,zip(*sorted(mapping.items()))); new=(p,periodic,nk,nm)
    matrix,constraints=reconstruction(old,new)
    nu,nv=map(count,axes); transverse=[nv,nu][which]
    rows=[]
    for fixed in range(transverse):
        controls=[grid[i*nv+fixed] if which==0 else grid[fixed*nv+i] for i in range(count(old))]
        if any(sum(x*controls[i][c] for i,x in equation) for equation in constraints for c in range(4)):
            if op=='I': raise AssertionError('inconsistent refinement equations')
            return False,surface
        row=tuple(tuple(sum(x*controls[i][c] for i,x in equation) for c in range(4)) for equation in matrix)
        if any(c[3]<=0 for c in row):
            if op=='I': raise AssertionError('nonpositive refinement weight')
            return False,surface
        rows.append(row)
    axes=list(axes); axes[which]=new; nu,nv=map(count,axes)
    grid=tuple(rows[j][i] if which==0 else rows[i][j] for i in range(nu) for j in range(nv))
    return True,(tuple(axes),grid)


def expected(row):
    name,surface,ops=parse(row); flags=[]
    for op,which,u,target in ops:
        flag,surface=edit(surface,op,which,u,target); flags.append(flag)
    return flags,surface


def encode(name,flags,surface,compact=False):
    axes,controls=surface
    w=[name,'D' if compact else 'R',len(flags),*map(int,flags),*[a[0] for a in axes],*[int(a[1]) for a in axes],*map(count,axes),*[len(a[2]) for a in axes],*[x for a in axes for x in domain(a)]]
    if compact:
        denominator=1
        for c in controls:
            for x in c: denominator=denominator//math.gcd(denominator,x.denominator)*x.denominator
        w.append(denominator); w.extend(x.numerator*(denominator//x.denominator) for c in controls for x in c)
    else: w.extend(x for c in controls for x in c)
    w.extend(x for a in axes for k,m in zip(a[2],a[3]) for x in [k,m])
    return ' '.join(map(str,w))


def decode(text,native=False):
    results={}
    for line in text.splitlines():
        w=iter(line.split()); name=next(w); kind=next(w)
        if name in results or kind not in ('R','D') or (native and kind!='R'): raise ValueError('duplicate or unsuccessful observation')
        flags=[int(next(w)) for _ in range(int(next(w)))]
        if any(x not in (0,1) for x in flags): raise ValueError('invalid operation flags')
        du,dv,pu,pv,nu,nv,ku,kv=(int(next(w)) for _ in range(8))
        number=(lambda:float(next(w))) if native else (lambda:F(next(w)))
        bounds=tuple(tuple(number() for _ in range(2)) for _ in range(2))
        denominator=int(next(w)) if kind=='D' else 1
        if denominator<=0: raise ValueError('invalid denominator')
        controls=tuple(tuple(number()/denominator for _ in range(4)) for _ in range(nu*nv))
        axes=[]
        for p,periodic,k in [(du,pu,ku),(dv,pv,kv)]:
            pairs=[(number(),int(next(w))) for _ in range(k)]
            knots,mults=map(tuple,zip(*pairs)); axes.append((p,bool(periodic),knots,mults))
            if not 1<=p<=25 or periodic not in (0,1) or k<2 or any(m<1 or m>p+(not periodic and i in (0,k-1)) for i,m in enumerate(mults)) or any(x>=y for x,y in zip(knots,knots[1:])):
                raise ValueError('invalid axis')
            if periodic and mults[0]!=mults[-1]: raise ValueError('inconsistent periodic seam')
        if list(w) or [count(a) for a in axes]!=[nu,nv] or nu<=du or nv<=dv or nu*nv>4096: raise ValueError('invalid grid dimensions')
        if any(c[3]<=0 for c in controls) or any(a>=b for a,b in bounds): raise ValueError('invalid weight or domain')
        if native and not all(math.isfinite(x) for x in [*[v for b in bounds for v in b],*[k for a in axes for k in a[2]],*[v for c in controls for v in c]]): raise ValueError('nonfinite native observation')
        if not native and tuple(domain(a) for a in axes)!=bounds: raise ValueError('wrong parameter domain')
        results[name]=(list(map(bool,flags)),(tuple(axes),controls),bounds)
    return results
