"""Independent degree-elevation oracle: all exact Cox coefficient equations.

The homogeneous control map is solved on complete support, including inactive
nonperiodic controls and a full periodic interval. It shares neither Boehm
insertion nor the production Prautzsch averaging algorithm.
"""
from fractions import Fraction as F
from functools import lru_cache
from knot_editing_reference import basis_polynomials,parse


def count(axis):
    p,periodic,k,m=axis
    return sum(m)-(m[0] if periodic else p+1)


def active_indices(axis):
    p,periodic,k,m=axis
    if periodic:return 0,len(k)-1
    total=0
    for first,mult in enumerate(m):
        total+=mult
        if total>p:break
    total=0
    for last in range(len(m)-1,-1,-1):
        total+=m[last]
        if total>p:break
    return first,last


def elevated_axis(axis,degree):
    p,periodic,knots,mults=axis
    if not p<=degree<=25:raise ValueError('invalid target degree')
    delta=degree-p
    first,last=active_indices(axis)
    pairs=[[k,m+delta] for k,m in zip(knots,mults)]
    if not periodic:
        for reverse,drop in [(False,delta*first),(True,delta*(len(knots)-1-last))]:
            if reverse:pairs.reverse()
            while drop:
                taken=min(drop,pairs[0][1]);drop-=taken;pairs[0][1]-=taken
                if not pairs[0][1]:pairs.pop(0)
            if reverse:pairs.reverse()
    knots,mults=map(tuple,zip(*pairs))
    result=(degree,periodic,knots,mults)
    assert count(result)==count(axis)+delta*(last-first)
    return result


def accumulate(target,source,scale):
    for i,x in source.items():
        y=target.get(i,F(0))+scale*x
        if y:target[i]=y
        else:target.pop(i,None)


@lru_cache(maxsize=256)
def coefficient_map(old,new):
    """Every coefficient on every common cell: controls as exact linear forms."""
    assert old[1]==new[1]
    cuts=sorted(set(old[2])|set(new[2]));pivots={};constraints=[]
    for lo,hi in zip(cuts,cuts[1:]):
        before=basis_polynomials(*old,lo,hi);after=basis_polynomials(*new,lo,hi)
        for power in range(max(old[0],new[0])+1):
            row={i:c[power] for i,c in enumerate(after) if power<len(c) and c[power]}
            rhs={i:c[power] for i,c in enumerate(before) if power<len(c) and c[power]}
            while row:
                j=min(row)
                if j not in pivots:
                    pivot=row[j]
                    pivots[j]=({i:x/pivot for i,x in row.items()},{i:x/pivot for i,x in rhs.items()})
                    break
                basis,answer=pivots[j];scale=-row[j]
                accumulate(row,basis,scale);accumulate(rhs,answer,scale)
            else:
                if rhs:constraints.append(tuple(sorted(rhs.items())))
    if constraints:raise ArithmeticError('non-nested polynomial spaces')
    if len(pivots)!=count(new):raise ArithmeticError('incomplete coefficient rank')
    controls=[None]*count(new)
    for j in reversed(range(count(new))):
        row,answer=pivots[j];answer=answer.copy()
        for i,x in row.items():
            if i!=j:accumulate(answer,controls[i],-x)
        controls[j]=answer
    return tuple(tuple(sorted(c.items())) for c in controls)


def elevate(curve,degree):
    p,periodic,knots,mults,controls=curve
    axis=curve[:4];new=elevated_axis(axis,degree)
    if degree==p:return curve
    if periodic:
        matrix=coefficient_map(axis,new)
        padded=controls;drop=0
    else:
        first,last=active_indices(axis)
        prefix,suffix=p+1-mults[0],p+1-mults[-1]
        work_mults=(p+1,*mults[1:-1],p+1)
        old_work=(p,False,knots,work_mults)
        new_work=(degree,False,knots,tuple(m+degree-p for m in work_mults))
        matrix=coefficient_map(old_work,new_work)
        zero=(F(0),)*4
        padded=(zero,)*prefix+controls+(zero,)*suffix
        drop=prefix+(degree-p)*first
        tail=suffix+(degree-p)*(len(knots)-1-last)
        work_flat=tuple(k for k,m in zip(new_work[2],new_work[3]) for _ in range(m))
        final_flat=tuple(k for k,m in zip(new[2],new[3]) for _ in range(m))
        assert work_flat[drop:len(work_flat)-tail if tail else None]==final_flat
        # Cropped basis functions vanish on the original active domain.
        assert all(work_flat[i+degree+1]<=knots[first] for i in range(drop))
        assert all(work_flat[i]>=knots[last] for i in range(len(matrix)-tail,len(matrix)))
    final=tuple(tuple(sum(x*padded[i][c] for i,x in row) for c in range(4)) for row in matrix[drop:drop+count(new)])
    if not all(c[3]>0 for c in final):raise ArithmeticError('nonpositive elevated weight')
    return (*new,final)


def elevate_surface(surface, u_degree, v_degree, order=(0, 1)):
    """Independent tensor product of the complete coefficient maps."""
    axes, grid = surface
    targets = (u_degree, v_degree)
    for which in order:
        nu, nv = map(count, axes)
        old = axes[which]
        rows = []
        for fixed in range((nv, nu)[which]):
            controls = tuple(grid[i*nv+fixed] if which == 0 else grid[fixed*nv+i]
                             for i in range(count(old)))
            curve = elevate((*old, controls), targets[which])
            rows.append(curve[4])
        axes = list(axes)
        axes[which] = curve[:4]
        axes = tuple(axes)
        nu, nv = map(count, axes)
        grid = tuple(rows[j][i] if which == 0 else rows[i][j]
                     for i in range(nu) for j in range(nv))
    return axes, grid
