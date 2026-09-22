"""Independent exact tensor Cox/power oracle. No native or Rust expectations.

Affine binomial substitutions implement edits; tensor Bernstein conversion is
performed only for output. Constant-parameter substitution supplies isocurves.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import comb, isfinite, lcm
from generate_spline_fixtures import axis
from bezier_editing_reference import add,times_linear,substitute,integer_matrix,apply_matrix


def parse(row):
    w=row.split(); du,dv,nu,nv,nuk,nvk,pu,pv=map(int,w[2:10])
    data=list(map(lambda x:F(float(x)),w[17:17+4*nu*nv])); tail=w[17+4*nu*nv:]
    if len(data)!=4*nu*nv or len(tail)!=2*(nuk+nvk): raise ValueError('invalid input count')
    return dict(name=w[0],degrees=(du,dv),counts=(nu,nv),periodic=(bool(pu),bool(pv)),
                domain=tuple(F(float(x)) for x in w[10:14]),op=int(w[14]),elevation=tuple(map(int,w[15:17])),
                controls=tuple(tuple([*(data[i+c]*data[i+3] for c in range(3)),data[i+3]]) for i in range(0,len(data),4)),
                knots=(tuple(map(float,tail[:2*nuk:2])),tuple(map(float,tail[2*nuk::2]))),
                mults=(tuple(map(int,tail[1:2*nuk:2])),tuple(map(int,tail[2*nuk+1::2]))))


@lru_cache(maxsize=128)
def spans(p,knots,mults,periodic,first,last):
    flat,start,end=axis(p,knots,mults,periodic); period=end-start
    if not first<last or (not periodic and not start<=first<last<=end): raise ValueError('invalid interval')
    turns=range((first-start)//period,(last-start)//period+1) if periodic else [0]
    result=[]
    for turn in turns:
        offset=turn*period
        for span in range(p,len(flat)-p-1):
            a,b=flat[span:span+2]
            if a>=b or a<start or b>end: continue
            low,high=max(first,a+offset),min(last,b+offset)
            if low>=high: continue
            @lru_cache(None)
            def basis(i,d):
                if not d: return [F(i==span)]
                left,right=flat[i+d]-flat[i],flat[i+d+1]-flat[i+1]
                lo=times_linear(basis(i,d-1),(a-flat[i])/left,(b-a)/left) if left else []
                hi=times_linear(basis(i+1,d-1),(flat[i+d+1]-a)/right,(a-b)/right) if right else []
                return add(lo,hi)
            coefficients=[substitute(basis(i,p),(low-a-offset)/(b-a),(high-low)/(b-a)) for i in range(span-p,span+1)]
            result.append((low,high,span-p,coefficients))
    return result


def along(grid,du,dv,direction,transform):
    """Apply a scalar polynomial transform independently to each coefficient row."""
    nu,nv=du+1,dv+1; result=None
    for fixed in range(nv if direction==0 else nu):
        indices=[i*nv+fixed for i in range(nu)] if direction==0 else [fixed*nv+j for j in range(nv)]
        rows=[transform([grid[i][c] for i in indices]) for c in range(4)]
        length=len(rows[0]); new_nu,new_nv=(length,nv) if direction==0 else (nu,length)
        if result is None: result=[[F(0)]*4 for _ in range(new_nu*new_nv)]
        for i in range(length): result[i*new_nv+fixed if direction==0 else fixed*new_nv+i]=[rows[c][i] for c in range(4)]
    return result,new_nu-1,new_nv-1


@lru_cache(maxsize=4)
def extract(degrees,counts,periodic,domain,controls,knots,mults):
    du,dv=degrees; nu,nv=counts; result=[]
    us=spans(du,knots[0],mults[0],periodic[0],*domain[:2])
    vs=spans(dv,knots[1],mults[1],periodic[1],*domain[2:])
    for ua,ub,ui,u in us:
        for va,vb,vi,v in vs:
            # Independent separable tensor sum, using Cox polynomial basis.
            local=[controls[((ui+a)%nu)*nv+(vi+b)%nv] for a in range(du+1) for b in range(dv+1)]
            um=integer_matrix([[u[a][i] for a in range(du+1)] for i in range(du+1)])
            vm=integer_matrix([[v[b][j] for b in range(dv+1)] for j in range(dv+1)])
            h,_,_=along(local,du,dv,0,lambda p:apply_matrix(p,um))
            h,_,_=along(h,du,dv,1,lambda p:apply_matrix(p,vm))
            result.append((du,dv,(ua,ub,va,vb),h))
    return result


def restrict(item,uf,ul,vf,vl):
    du,dv,domain,h=item; ua,ub,va,vb=domain
    h,du,dv=along(h,du,dv,0,lambda p:substitute(p,uf,ul-uf))
    h,du,dv=along(h,du,dv,1,lambda p:substitute(p,vf,vl-vf))
    return du,dv,(ua+uf*(ub-ua),ua+ul*(ub-ua),va+vf*(vb-va),va+vl*(vb-va)),h


@lru_cache(maxsize=128)
def bernstein_matrix(n,degree):
    return integer_matrix([[F(comb(i,k),comb(degree,k)) if k<=i else F(0) for k in range(n)] for i in range(degree+1)])


def bernstein(p,degree):
    return apply_matrix(p,bernstein_matrix(len(p),degree))


def expected(row):
    v=parse(row); result=[]; op=v['op']
    patches=extract(*(v[k] for k in ['degrees','counts','periodic','domain','controls','knots','mults']))
    for patch in patches:
        if op in [1,10]: patch=restrict(patch,F(1,4),F(3,4),F(1,4),F(3,4))
        du,dv,domain,h=patch
        eu,ev=v['elevation'] if op in [7,10] else (du,dv)
        if op in [4,10]: h,du,dv=along(h,du,dv,0,lambda p:substitute(p,F(1),F(-1)))
        if op==5: h,du,dv=along(h,du,dv,1,lambda p:substitute(p,F(1),F(-1)))
        if op in [6,10]:
            h=[h[i*(dv+1)+j] for j in range(dv+1) for i in range(du+1)]
            du,dv=dv,du; eu,ev=ev,eu; domain=domain[2:]+domain[:2]
        if op in [8,9]:
            direction=0 if op==8 else 1; t=F(3,8) if op==8 else F(5,8)
            h,du,dv=along(h,du,dv,direction,lambda p:[sum(x*t**i for i,x in enumerate(p))])
            degree=ev if op==8 else eu
            controls=[[bernstein([p[c] for p in h],degree)[i] for c in range(4)] for i in range(degree+1)]
            result.append(('C',degree,0,domain[2:] if op==8 else domain[:2],controls)); continue
        u_ranges=[(F(0),F(3,8)),(F(3,8),F(1))] if op in [2,10] else [(F(0),F(1))]
        v_ranges=[(F(0),F(5,8)),(F(5,8),F(1))] if op in [3,10] else [(F(0),F(1))]
        for a,b in u_ranges:
            for c,d in v_ranges:
                pd,pv,domain_out,g=restrict((du,dv,domain,h),a,b,c,d)
                g,pd,pv=along(g,pd,pv,0,lambda p:bernstein(p,eu))
                g,pd,pv=along(g,pd,pv,1,lambda p:bernstein(p,ev))
                assert all(p[3]>0 for p in g)
                result.append(('P',eu,ev,domain_out,g))
    return result


def encode_result(name,items):
    return ' '.join([name,'R',str(len(items)),*[str(x) for kind,du,dv,domain,controls in items for x in [kind,du,dv,*domain,*[x for p in controls for x in p]]]])


def encode_fixture(name,items):
    words=[name,'D',str(len(items))]
    for kind,du,dv,domain,controls in items:
        values=[x for p in controls for x in p]
        denominator=lcm(*(x.denominator for x in values))
        words.extend(map(str,[kind,du,dv,*domain,denominator,*[x.numerator*(denominator//x.denominator) for x in values]]))
    return ' '.join(words)


def decode(text,native=False):
    result={}
    for row in text.splitlines():
        w=iter(row.split()); name=next(w); status=next(w)
        if name in result: raise ValueError('duplicate result')
        if native and status=='E':
            if list(w): raise ValueError('malformed exception')
            result[name]=None; continue
        if status!='R': raise ValueError('unsuccessful result')
        count=int(next(w)); items=[]
        if not 1<=count<=16384: raise ValueError('invalid count')
        number=(lambda:float(next(w))) if native else (lambda:F(next(w)))
        for _ in range(count):
            kind=next(w); du,dv=int(next(w)),int(next(w))
            if kind not in ['P','C'] or not 1<=du<=25 or not (1<=dv<=25 if kind=='P' else dv==0): raise ValueError('invalid item')
            domain=tuple(number() for _ in range(4 if kind=='P' else 2))
            if not all(domain[i]<domain[i+1] for i in range(0,len(domain),2)): raise ValueError('invalid domain')
            controls=[[number() for _ in range(4)] for _ in range((du+1)*(dv+1))]
            if not native and not all(p[3]>0 for p in controls): raise ValueError('nonpositive exact weight')
            if native and not all(isfinite(x) for x in domain): raise ValueError('nonfinite domain')
            items.append((kind,du,dv,domain,controls))
        if list(w): raise ValueError('trailing result tokens')
        result[name]=items
    return result
