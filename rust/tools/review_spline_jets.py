"""Fail-closed review of observed OCCT spline evaluation discrepancies.

Reviews pin the complete native jet, version, input, and independent rational
answer. This module never creates reviews or widens a comparison budget.
"""
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path
import struct

ROOT=Path(__file__).resolve().parents[2]


@lru_cache(maxsize=1024)
def exact_jet(family,row):
    from generate_spline_fixtures import axis, parameters, jet
    if family=='curves':
        w=row.split()
        degree,np,nk=map(int,w[2:5])
        data=list(map(lambda x:F(float(x)),w[8:8+4*np]))
        tail=w[8+4*np:]
        if len(tail)!=2*nk: raise ValueError('invalid curve data')
        knots=list(map(float,tail[::2])); mults=list(map(int,tail[1::2]))
        flat,start,end=axis(degree,knots,mults,w[1]=='P')
        at=parameters(flat,start,end,float(w[5]),w[6],w[1]=='P')[0]
        poles=[data[4*(i%np):4*(i%np)+3] for i in range(len(flat)-degree-1)]
        weights=[data[4*(i%np)+3] for i in range(len(poles))]
        values=jet(degree,poles,weights,flat,*at,int(w[7]))
    elif family=='surfaces':
        from generate_surface_fixtures import parse, jet as surface_jet
        _,_,du,dv,poles,weights,uk,um,vk,vm,pu,pv,u,v,su,sv,order=parse(row)
        uf,us,ue=axis(du,uk,um,pu); vf,vs,ve=axis(dv,vk,vm,pv)
        nu=sum(um[:-1]) if pu else sum(um)-du-1
        nv=sum(vm[:-1]) if pv else sum(vm)-dv-1
        values=surface_jet(du,dv,[list(map(F,p)) for p in poles],list(map(F,weights)),nu,nv,uf,vf,
                           parameters(uf,us,ue,u,su,pu)[0],parameters(vf,vs,ve,v,sv,pv)[0],order)
    else: raise ValueError('invalid spline family')
    return tuple(x for r in values for x in r)


def exact_hash(values):
    return hashlib.sha256(','.join(f'{x.numerator}/{x.denominator}' for x in values).encode()).hexdigest()


def minimal_bounds(value):
    # CPython's Fraction-to-float conversion supplies the nearest value; exact
    # comparison picks its adjacent endpoint. This does not call Rust.
    nearest=float(value)
    represented=F(nearest)
    return (math.nextafter(nearest,-math.inf) if represented>value else nearest,
            math.nextafter(nearest, math.inf) if represented<value else nearest)


@lru_cache(maxsize=1)
def reviews():
    path=ROOT/'rust/fixtures/occt-spline-jet-divergences.json'
    return json.loads(path.read_text())


def reviewed_jet(family,oracle,label,native,rust,row):
    for review in reviews():
        if (review['family'],review['oracle'],review['case'])!=(family,oracle,label): continue
        if review['input_sha256']!=hashlib.sha256(row.encode()).hexdigest(): continue
        if review['native_bits']!=[struct.pack('>d',x).hex() for x in native]: continue
        exact=exact_jet(family,row)
        if exact_hash(exact)!=review['exact_jet_sha256']: continue
        if len(rust)!=len(exact) or any(r not in minimal_bounds(x) for r,x in zip(rust,exact)): continue
        return {'id':review['id'],'case':label,'family':family,'input_sha256':review['input_sha256'],
                'exact_jet_sha256':review['exact_jet_sha256'],'evidence':review['evidence'],
                'reason':review['reason']}
    return None
