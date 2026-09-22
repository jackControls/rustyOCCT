"""Independent rational boundary-intersection oracle.

Closed-form cross-product line/plane formulas, edge intersections, and a
3D gift-wrapping hull. No affine RREF or halfspace-vertex enumeration used by
production. Fractions use Python integers independently of Rust bigint code.
"""
from fractions import Fraction as Q
from proximity_reference import add, sub, mul, dot, cross, norm2, normal, edges, inside


def contains(shape, point):
    kind,p=shape
    if kind=='P': return point==p[0]
    if kind in 'LS':
        d=sub(p[1],p[0]); w=sub(point,p[0])
        if norm2(d)==0: return point==p[0]
        return norm2(cross(d,w))==0 and (kind=='L' or 0<=dot(w,d)<=norm2(d))
    return dot(normal(p),sub(point,p[0]))==0 and (kind=='F' or inside(point,p))


def line(p,d):
    pivot=next(x for x in d if x)
    direction=mul(d,1/pivot)
    axis=next(i for i,x in enumerate(d) if x)
    return 'L',[sub(p,mul(direction,p[axis])),direction]


def plane(p):
    n=normal(p); k=next(x for x in n if x)
    return 'F',[mul(n,1/k),(dot(n,p[0])/k,)]


def hull(points):
    points=sorted(set(points))
    if not points: return 'E',[]
    if len(points)==1: return 'P',points
    d=sub(points[1],points[0])
    normals=[cross(d,sub(p,points[0])) for p in points[2:]]
    n=next((n for n in normals if norm2(n)),None)
    if n is None: return 'S',[points[0],points[-1]]
    result=[points[0]]
    while True:
        a=result[-1]
        b=next(p for p in points if p!=a)
        for c in points:
            turn=dot(n,cross(sub(b,a),sub(c,a)))
            if turn<0 or (turn==0 and norm2(sub(c,a))>norm2(sub(b,a))): b=c
        if b==result[0]: break
        assert b not in result
        result.append(b)
    if result[1]>result[-1]: result[1:]=reversed(result[1:])
    return 'G',result


def line_line(a,b):
    ka,pa=a; kb,pb=b
    u,v,w=sub(pa[1],pa[0]),sub(pb[1],pb[0]),sub(pb[0],pa[0])
    n=cross(u,v)
    if norm2(n):
        if dot(w,n): return 'E',[]
        t=dot(cross(w,v),n)/norm2(n)
        s=dot(cross(w,u),n)/norm2(n)
        if (ka=='S' and not 0<=t<=1) or (kb=='S' and not 0<=s<=1): return 'E',[]
        return 'P',[add(pa[0],mul(u,t))]
    if norm2(cross(w,u)): return 'E',[]
    if ka==kb=='L': return line(pa[0],u)
    return hull([p for k,ps,other in [(ka,pa,b),(kb,pb,a)] if k=='S' for p in ps if contains(other,p)])


def line_plane(a,b):
    ka,pa=a; _,pb=b
    d=sub(pa[1],pa[0]); n=normal(pb)
    den=dot(d,n); num=dot(sub(pb[0],pa[0]),n)
    if not den:
        if num: return 'E',[]
        return hull(pa) if ka=='S' else line(pa[0],d)
    t=num/den
    return ('E',[]) if ka=='S' and not 0<=t<=1 else ('P',[add(pa[0],mul(d,t))])


def intersection(a,b):
    ka,pa=a; kb,pb=b
    if ka=='S' and pa[0]==pa[1]: return intersection(('P',[pa[0]]),b)
    if kb=='S' and pb[0]==pb[1]: return intersection(a,('P',[pb[0]]))
    if ka=='P': return ('P',pa) if contains(b,pa[0]) else ('E',[])
    if kb=='P': return intersection(b,a)
    if ka in 'LS':
        if kb in 'LS': return line_line(a,b)
        hit=line_plane(a,('F',pb))
        if kb=='F' or hit[0]=='E': return hit
        if hit[0]=='P': return hit if contains(b,hit[1][0]) else ('E',[])
        # Coplanar line/segment against a triangle: collect boundary crossings.
        points=[p for p in pa if ka=='S' and contains(b,p)]
        for edge in edges(pb):
            result=intersection(a,('S',list(edge)))
            assert result[0] in 'EPS'
            points+=result[1]
        return hull(points)
    if kb in 'LS': return intersection(b,a)
    if ka==kb=='F':
        u,v=normal(pa),normal(pb); d=cross(u,v)
        if not norm2(d): return plane(pa) if contains(a,pb[0]) else ('E',[])
        p=mul(add(mul(cross(v,d),dot(u,pa[0])),mul(cross(d,u),dot(v,pb[0]))),1/norm2(d))
        return line(p,d)
    if ka=='F': return intersection(b,a)
    points=[]
    for edge in edges(pa):
        result=intersection(('S',list(edge)),b)
        assert result[0] in 'EPS'
        points+=result[1]
    if kb=='T':
        for edge in edges(pb):
            result=intersection(('S',list(edge)),a)
            assert result[0] in 'EPS'
            points+=result[1]
    return hull(points)


def encode(result):
    kind,values=result
    return ' '.join([kind,str(len(values)),*[str(x.numerator)+'/'+str(x.denominator) for v in values for x in v]])


def decode(words):
    kind=next(words); n=int(next(words))
    if kind not in 'EPSGLF' or n not in range(7): raise ValueError('invalid exact result')
    expected={'E':0,'P':1,'S':2,'L':2,'F':2}
    if (kind in expected and n!=expected[kind]) or (kind=='G' and n<3): raise ValueError('invalid result count')
    return kind,[tuple(Q(next(words)) for _ in range(1 if kind=='F' and i==1 else 3)) for i in range(n)]
