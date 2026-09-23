#!/usr/bin/env python3
"""Independent tensor coefficient fixtures, including native-unsupported extremes."""
import argparse
import math
from pathlib import Path
import random
import struct
from compare_surface_editing import cases,encode
from surface_editing_reference import expected,encode_fixture
ROOT=Path(__file__).resolve().parents[1]


def inputs():
    # All native families. Full degree-25 x 25 output is recomputed for every
    # operation in the native bridge; retain extraction and the combined sequence
    # for each periodicity in Cargo's complete-control fixture file.
    rows=[r for r in cases().splitlines() if r.split()[2:4]!=['25','25'] or ('_w0_' in r.split()[0] and int(r.split()[14]) in [0,10])]
    tiny=float.fromhex('0x0.0000000000001p-1022'); maximum=float.fromhex('0x1.fffffffffffffp+1023')
    basic=[(float(i),float(j),float(i*j)) for i in range(3) for j in range(3)]
    families=[
        ('tiny_rectangle',2,2,basic,[1.,2.,3.,2.,1.,2.,3.,2.,1.],[0.,tiny],[3,3],[0.,tiny],[3,3],False,False,(0.,tiny,0.,tiny)),
        ('huge_rectangle',1,1,[(0.,0.,0.),(0.,1.,0.),(1.,0.,0.),(1.,1.,1.)],[1.]*4,[-maximum,maximum],[2,2],[-maximum,maximum],[2,2],False,False,(-maximum,maximum,-maximum,maximum)),
        ('extreme_homogeneous',1,1,[(maximum,0.,-maximum),(0.,maximum,0.),(-maximum,0.,maximum),(1.,2.,3.)],[tiny,maximum,1.,tiny],[0.,1.],[2,2],[0.,1.],[2,2],False,False,(0.,1.,0.,1.)),
        ('subulp_turns',1,1,[(float(i),float(j),float(i*j)) for i in range(3) for j in range(3)],[1.+i%3 for i in range(9)],[0.,.125,.375,.5],[1]*4,[0.,.125,.375,.5],[1]*4,True,True,(2.**53,2.**53+2.,-2.**53-2.,-2.**53)),
        ('tiny_turns',1,1,[(0.,0.,0.),(1.,0.,1.),(0.,1.,2.),(1.,1.,0.)],[1.,2.,3.,4.],[0.,tiny,4*tiny],[1]*3,[0.,tiny,4*tiny],[1]*3,True,True,(tiny,5*tiny,-tiny,5*tiny)),
        ('adjacent_knots',2,2,[(float(i),float(j),float(i-j)) for i in range(4) for j in range(4)],[1.+i%3 for i in range(16)],[1.,math.nextafter(1.,math.inf),2.],[3,1,3],[1.,math.nextafter(1.,math.inf),2.],[3,1,3],False,False,(1.,2.,1.,2.)),
    ]
    for name,du,dv,poles,weights,uk,um,vk,vm,pu,pv,window in families:
        for op in range(11):
            rows.append(encode((name,'S',du,dv,poles,weights,uk,um,vk,vm,pu,pv),window,op,name=f'{name}_op{op}'))
    rng=random.Random(0x5FACE)
    for i in range(32):
        du,dv=1+i%5,1+(i//5)%5; pu,pv=bool(i&1),bool(i&2)
        uk,vk=[.1,.4,1.3],[.2,.7,1.9]
        um=[1,du,1] if pu else [du+1,1,du+1]; vm=[1,dv,1] if pv else [dv+1,1,dv+1]
        nu=sum(um[:-1]) if pu else sum(um)-du-1; nv=sum(vm[:-1]) if pv else sum(vm)-dv-1
        scale=math.ldexp(1.,rng.randint(-800,800))
        poles=[tuple(scale*rng.randint(-5,5) for _ in range(3)) for _ in range(nu*nv)]
        weights=[math.ldexp(float(rng.randint(1,5)),rng.randint(-12,12)) for _ in poles]
        shape=(f'random_{i}','S',du,dv,poles,weights,uk,um,vk,vm,pu,pv)
        window=(*((-1.,2.) if pu else (.15,1.2)),*((-2.,2.5) if pv else (.25,1.8)))
        rows.append(encode(shape,window,i%11))
    # Ordinary regressions for the retained degree-25 tensor fuzz
    # inputs. Only decoding is shared conceptually; exact expectations still
    # come from the independent Cox tensor/power oracle.
    for filename in ['unclamped-degree25-combined.bin','unclamped-degree25-extraction.bin','mutated-unclamped-degree25.bin','clamped-degree25-combined.bin','periodic-degree25-combined.bin']:
        data=(ROOT/'fuzz/regressions/surface_editing'/filename).read_bytes()
        assert list(data[:3])==[2,24,24] and data[3]==data[4] and data[5] in [0,10]
        assert data[10]%8 not in [0,1,2] and data[11]%8 not in [0,1,2]
        kind=data[3]
        if kind==2: knots,mults,periodic,n=list(map(float,range(54))),[1]*54,False,28
        elif kind==1: knots,mults,periodic,n=[0.,1.,3.],[1,25,1],True,26
        else:
            assert kind==0
            knots,mults,periodic,n=[0.,1.,3.],[26,1,26],False,27
        def byte(i): return data[i] if i<len(data) else 0
        def signed(i):
            x=byte(i); return float(x if x<128 else x-256)
        poles=[tuple(signed(16+4*i+c) for c in range(3)) for i in range(n*n)]
        weights=[float(1+byte(16+4*i+3)%8) for i in range(n*n)]
        shape=('fuzz_'+filename[:-4].replace('-','_'),'S',25,25,poles,weights,knots,mults,knots,mults,periodic,periodic)
        intervals=[(25.,26.),(26.,27.),(27.,28.)] if kind==2 else [(0.,1.),(1.,3.)]
        u,v=intervals[data[12]%len(intervals)],intervals[data[13]%len(intervals)]
        offset=-6. if periodic and data[8]&32 else 0.
        rows.append(encode(shape,tuple(x+offset for x in (*u,*v)),data[5]))
    # Complete full-exponent periodic input saved by Linux sanitizer fuzzing.
    data=(ROOT/'fuzz/regressions/surface_editing/slow-unit-3e01ced2e0056a1159d07d2a69d7405e6dddf050.bin').read_bytes()
    assert list(data[:16])==[0,1,1,1,1,10,2,2,192,255,255,255,84,72,85,170]
    values=[struct.unpack('<d',data[16+8*i:24+8*i])[0] for i in range(36)]
    poles=[tuple(values[4*i:4*i+3]) for i in range(9)]
    weights=values[3::4]
    shape=('fuzz_full_exponent_periodic','S',2,2,poles,weights,[0.,1.,3.],[1,2,1],[0.,1.,3.],[1,2,1],True,True)
    rows.append(encode(shape,(-3.,6.,-3.,6.),10))
    return rows


def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--check',action='store_true'); args=p.parse_args()
    rows=inputs(); text='# input | exact tensor controls: D uses one shared control denominator per item; see SURFACE_EDITING.md\n'+'\n'.join(row+' | '+encode_fixture(row.split()[0],expected(row)) for row in rows)+'\n'
    path=ROOT/'fixtures/surface-editing.tsv'
    if args.check:
        if path.read_text()!=text: p.error('surface editing fixtures changed; investigate before updating')
    else: path.write_text(text)
    print(f'surface-editing.tsv: {len(rows)} complete exact cases, {len(text)} bytes')

if __name__=='__main__': main()
