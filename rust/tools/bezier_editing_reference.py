"""Independent exact Cox basis/power coefficient oracle, not spline blossoming.

Each edit is an affine substitution in homogeneous power polynomials. Convert
those coefficients to Bernstein controls only at the end. No de Casteljau,
degree-elevation recurrence, native outputs or Rust samples supply expectations.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import comb, lcm
from generate_spline_fixtures import axis


def parse(row):
    w=row.split(); degree,np,nk=map(int,w[2:5])
    data=list(map(lambda x:F(float(x)),w[9:9+4*np])); tail=w[9+4*np:]
    if len(tail)!=2*nk or len(data)!=4*np: raise ValueError('invalid input count')
    return dict(name=w[0],periodic=w[1]=='P',degree=degree,
                first=F(float(w[5])),last=F(float(w[6])),op=int(w[7]),elevation=int(w[8]),
                poles=[data[i:i+3] for i in range(0,len(data),4)],weights=data[3::4],
                knots=list(map(float,tail[::2])),mults=list(map(int,tail[1::2])))


def add(a,b):
    result=[F(0)]*max(len(a),len(b))
    for p in [a,b]:
        for i,x in enumerate(p): result[i]+=x
    return result


def times_linear(p,a,b):
    result=[F(0)]*(len(p)+1)
    for i,x in enumerate(p): result[i]+=a*x; result[i+1]+=b*x
    return result


def integer_matrix(rows):
    result=[]
    for row in rows:
        denominator=lcm(*(x.denominator for x in row))
        result.append(([x.numerator*(denominator//x.denominator) for x in row],denominator))
    return result


def apply_matrix(p,matrix):
    denominator=lcm(*(x.denominator for x in p))
    integers=[x.numerator*(denominator//x.denominator) for x in p]
    return [F(sum(x*y for x,y in zip(integers,row)),denominator*d) for row,d in matrix]


@lru_cache(maxsize=256)
def affine_matrix(n,a,b):
    return integer_matrix([[F(comb(j,i))*a**(j-i)*b**i if j>=i else F(0) for j in range(n)] for i in range(n)])


def substitute(p,a,b):
    # The same direct binomial substitution, with denominators cleared before
    # matrix multiplication. Only final coefficients require rational reduction.
    return apply_matrix(p,affine_matrix(len(p),a,b))


@lru_cache(maxsize=256)
def homogeneous(degree,poles,weights,flat,span):
    a,length=flat[span],flat[span+1]-flat[span]
    @lru_cache(None)
    def basis(i,p):
        if not p: return [F(i==span)]
        left,right=flat[i+p]-flat[i],flat[i+p+1]-flat[i+1]
        lo=times_linear(basis(i,p-1),(a-flat[i])/left,length/left) if left else []
        hi=times_linear(basis(i+1,p-1),(flat[i+p+1]-a)/right,-length/right) if right else []
        return add(lo,hi)
    result=[[F(0)]*(degree+1) for _ in range(4)]
    for i in range(span-degree,span+1):
        j=i%len(poles); coefficients=basis(i,degree)
        for c in range(4):
            scale=weights[j]*(poles[j][c] if c<3 else 1)
            for k,x in enumerate(coefficients): result[c][k]+=scale*x
    return result


def expected(row):
    v=parse(row); p=v['degree']; flat,start,end=axis(p,v['knots'],v['mults'],v['periodic'])
    first,last=v['first'],v['last']; period=end-start
    if not first<last or (not v['periodic'] and not start<=first<last<=end): raise ValueError('invalid interval')
    turns=range((first-start)//period,(last-start)//period+1) if v['periodic'] else [0]
    result=[]
    for turn in turns:
        offset=turn*period
        for span in range(p,len(flat)-p-1):
            a,b=flat[span:span+2]
            if a>=b or a<start or b>end: continue
            low,high=max(first,a+offset),min(last,b+offset)
            if low>=high: continue
            h=homogeneous(p,tuple(map(tuple,v['poles'])),tuple(v['weights']),tuple(flat),span)
            h=[substitute(c,(low-a-offset)/(b-a),(high-low)/(b-a)) for c in h]
            op=v['op']; degree=v['elevation'] if op in [4,5] else p
            if op in [1,5]:
                h=[substitute(c,F(1,4),F(1,2)) for c in h]
                low,high=low+(high-low)/4,low+3*(high-low)/4
            if op in [3,5]: h=[substitute(c,F(1),F(-1)) for c in h]
            ranges=[(F(0),F(3,8)),(F(3,8),F(1))] if op in [2,5] else [(F(0),F(1))]
            for x,y in ranges:
                coefficients=[substitute(c,x,y-x) for c in h]
                controls=[[sum(c[k]*F(comb(i,k),comb(degree,k)) for k in range(min(i,len(c)-1)+1)) for c in coefficients] for i in range(degree+1)]
                assert all(c[3]>0 for c in controls)
                result.append((degree,low+x*(high-low),low+y*(high-low),controls))
    return result


def encode(name,curves):
    return ' '.join([name,'R',str(len(curves)),*[x for d,a,b,controls in curves for x in [str(d),str(a),str(b),*map(str,[v for c in controls for v in c])]]])


def decode(text,native=False):
    import math
    result={}
    for row in text.splitlines():
        w=iter(row.split()); name=next(w)
        if name in result or next(w)!='R': raise ValueError('invalid/unsuccessful result row')
        count=int(next(w)); curves=[]
        if not 1<=count<=4096: raise ValueError('invalid arc count')
        number=(lambda:float(next(w))) if native else (lambda:F(next(w)))
        for _ in range(count):
            degree=int(next(w)); a,b=number(),number()
            if not 1<=degree<=25 or not a<b: raise ValueError('invalid arc')
            controls=[[number() for _ in range(4)] for _ in range(degree+1)]
            if not all(c[3]>0 for c in controls): raise ValueError('nonpositive weight')
            if native and not all(math.isfinite(x) for x in [a,b,*[x for c in controls for x in c]]): raise ValueError('nonfinite native value')
            curves.append((degree,a,b,controls))
        if list(w): raise ValueError('trailing result tokens')
        result[name]=curves
    return result
