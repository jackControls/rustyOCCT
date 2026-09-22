#!/usr/bin/env python3
"""OCCT spline/sphere and spline/infinite-cylinder observations.

The shared text protocol prefixes the curve kind with sphere: or cylinder:.
Its nine surface values are center/origin, axis, and (radius,0,0).
"""
from compare_spline_plane import encode, main as compare


def cases():
    shapes=[]
    def bezier(name,poles,weights=None):
        d=len(poles)-1
        shapes.append((name,'B',d,poles,weights or [1.]*len(poles),[0.,1.],[d+1]*2))
    bezier('secant',[(-2.,0.,0.),(2.,0.,0.)])
    bezier('tangent',[(-2.,1.,0.),(2.,1.,0.)])
    bezier('outside',[(-2.,2.,0.),(2.,2.,0.)])
    bezier('inside',[(-.5,0.,0.),(.5,0.,0.)])
    bezier('constant',[(1.,0.,0.)]*3)
    bezier('rational_arc',[(1.,0.,0.),(1.,1.,0.),(0.,1.,0.)],[1.,1.,2.])
    bezier('generator',[(1.,0.,-2.),(1.,0.,2.)])
    for d in [2,3,5,8,13,25]:
        bezier(f'power_{d}',[(0.,0.,0.)]*d+[(2.,0.,0.)])
    bezier('rational_space',[(-2.,0.,-1.),(2.,3.,1.),(-2.,-1.,2.),(2.,0.,-1.)],[1.,2.,.5,1.])
    shapes.append(('partial_overlap','S',2,[(1.,0.,0.),(1.,1.,0.),(0.,1.,0.),(-1.,1.,0.),(-2.,0.,0.)],[1.,1.,2.,1.,1.],[0.,1.,2.],[3,2,3]))
    shapes.append(('corner','S',1,[(-2.,0.,0.),(-1.,0.,0.),(-2.,1.,0.)],[1.]*3,[0.,1.,2.],[2,1,2]))
    shapes.append(('unclamped','S',2,[(-2.,0.,0.),(0.,1.,1.),(2.,0.,0.),(0.,-1.,-1.),(-2.,0.,0.)],[1.,2.,1.,3.,1.],[-2.,-1.,0.,1.,2.,3.,4.,5.],[1]*8))
    shapes.append(('periodic','P',1,[(-2.,0.,0.),(0.,0.,0.),(2.,0.,0.),(0.,0.,0.)],[1.]*4,[0.,1.,2.,3.,4.],[1]*5))
    rows=[]
    for primitive in ['sphere','cylinder']:
        for frame in [0,1]:
            # Exact signed-permutation/scale covariance; no rounded rotation.
            center=[0.,0.,0.] if frame==0 else [2.,-3.,5.]
            scale=1. if frame==0 else .125
            def transform(p):
                p=p if frame==0 else [p[2],-p[0],p[1]]
                return [center[c]+scale*p[c] for c in range(3)]
            axis=[0.,0.,1.] if frame==0 else [1.,0.,0.]
            for name,kind,d,poles,weights,knots,mults in shapes:
                header=[center,axis,[scale,0.,0.]]
                label=f'{primitive}_{name}_{frame}'
                rows.append(encode(label,primitive+':'+kind,d,header,list(map(transform,poles)),weights,knots,mults))
                if name in ['secant','tangent','rational_arc','partial_overlap','periodic']:
                    a,b=(3.5,8.5) if kind=='P' else (.25,.75)
                    rows.append(encode(label+'_trim',primitive+':'+kind+'T',d,header,list(map(transform,poles)),weights,knots,mults)+f' {a} {b}')
    return '\n'.join(rows)+'\n'


if __name__=='__main__':
    compare(cases,'quadric')
