#!/usr/bin/env python3
"""Run bounded, instrumented libFuzzer campaigns with persistent corpus/artifacts.

Requires cargo-fuzz 0.13.1 and nightly Rust. Crashes/timeouts/OOMs fail the run.
Always writes a JSON manifest; raw logs retain coverage, executions and crashes.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import random
import re
import shutil
import signal
import subprocess
import tempfile
import time

ROOT = Path(__file__).resolve().parents[2]
FUZZ = ROOT/'rust/fuzz'
TARGETS = ['predicates','intersections','modeling','curved','splines','surfaces','roots','spline_intersections','proximity','linear_sets','bezier_editing','surface_editing','knot_editing','exact_spline_intersections','surface_knots','degree_elevation','spline_proximity','spline_linear','brep_validation','identity','history','split_merge','attributes','brep_io','analytic_intersections','tessellation','curve_surface','curve_curve','split','step','boolean']
STARTUP_SECONDS = 600
MAX_STARTUP_SECONDS = 3600
# surface_knots' retained CI corpus replays slower than the cap allows: 2,766 s
# for 564 inputs at 1a76d29e, over 3,600 s for 577 at 0c94aa53 (REVIEW_NOTES R12).
TARGET_MAX_STARTUP_SECONDS = {"surface_knots": 7200}
INPUT_SECONDS = 20
# Full tensor coefficient equations and double-axis degree-25 edits are a
# larger per-input workload. Existing targets keep their original 20s limit.
# surface_editing joined on measurement: a CI input took 8.2 s under local
# AddressSanitizer, above 20 s / 2.6 (fuzz/regressions/README.md).
# split joined in S8d.2: a cap split by a plane passing near its pole takes
# 8.8 s under ASan (its certified quadrature refines where the section's
# angle about the axis turns fast), within 20 s but without margin.
# boolean joined in S9b.2: its results as inputs again (the stored model's
# exact fragments) and S9a.2's spline stacks take up to about 30 s under ASan
# (fuzz/regressions/README.md).
TARGET_INPUT_SECONDS = {"surface_knots": 60, "degree_elevation": 60, "surface_editing": 60,
                        "split": 60, "boolean": 60}
# Targets whose exact oracles churn enough temporary BigInts that default ASan
# quarantine and allocator retention, not live data, exhaust the RSS gate.
# brep_validation joined in S4d: a 180 s campaign reached 2,064 MB after
# replay while the input it stopped on alone runs in 66 ms.
# analytic_intersections joined in S7b.2: a 300 s campaign stopped at the
# 2,048 MB RSS limit on two planes (0.011 s alone), after rational-interval
# cone isolation had churned temporary BigInts. curve_surface joined in S7c.1:
# a 600 s campaign peaked at 1,904 MB after exact circle/cone and torus
# resultants. curve_curve joined in S7d.1: 1,351 MB within a 600 s campaign.
# split joined in S8d.2: a clean 600 s campaign at a4e1c9df peaked at 1,903 MB
# (projection jets and exact plane decisions churning temporaries).
# boolean joined in S9a.1: a 600 s campaign at d7e515d9 left replay at
# 1,443 MB and stopped at the 2,048 MB limit on an input that peaks at 6 MB
# alone (exact arrangements churning BigRationals).
ALLOCATOR_TARGETS = {'surface_knots', 'degree_elevation', 'spline_linear', 'brep_validation',
                     'analytic_intersections', 'curve_surface', 'curve_curve', 'split',
                     'boolean'}
# Targets whose allocation stack traces are kept to five frames: with the
# default thirty, AddressSanitizer's stack depot grew analytic_intersections
# to 1,489 MB in 120 s (33 MB without a sanitizer); five frames keep it at
# 227 MB and still name each allocation's site in a report. curve_surface
# and curve_curve joined with it (the same exact intersection machinery).
# boolean joined in S9a.2: with the allocator purge, a 600 s campaign at
# f510f3a6 left replay at 1,849 MB and stopped at the gate with 28 MB live
# and 53 MB quarantined, the rest the depot of its exact arrangements'
# stacks. split joined after the scheduled run 37008675181 stopped at the
# gate with 50 MB live and 59 MB quarantined: one torus's spiric split
# (its certified projections' rational intervals) records 985,416 distinct
# 30-frame stacks, 248 MB of depot (14,849 and 9 MB with five frames), and
# 1,400 corpus inputs replayed locally left 3.9 million (941 MB, 1,350 MB
# RSS) against 88,036 (10 MB, 434 MB RSS) with five.
SHORT_STACK_TARGETS = {'analytic_intersections', 'curve_surface', 'curve_curve', 'boolean',
                       'split'}
# The pinned libFuzzer checks stop_file between MutateAndTestOne batches,
# not between each callback. Keep its default mutation sequence length.
MUTATION_DEPTH = 5
SHUTDOWN_SECONDS = MUTATION_DEPTH * INPUT_SECONDS + 5


def shutdown_budget(input_seconds):
    return MUTATION_DEPTH * input_seconds + 5


def startup_budget(corpus_files, input_seconds=INPUT_SECONDS, cap=MAX_STARTUP_SECONDS):
    # Exact-arithmetic seeds have widely varying costs. An assumed average
    # cost can reject an otherwise legal corpus before mutation ever starts.
    # Reserve build time and each input's configured budget (including the
    # initial empty input), capped at one hour. Every input keeps its own
    # timeout; replay never borrows the requested mutation duration.
    return min(cap, STARTUP_SECONDS+input_seconds*(corpus_files+1))


def sanitizer_build_args(target):
    # These targets also purge freed ASan allocator memory during corpus
    # replay; libFuzzer itself does that only after mutation has started.
    return ['--sanitizer','address'] + (['--features','asan-allocator'] if target in ALLOCATOR_TARGETS else [])


def campaign_environment(target, base):
    env=dict(base)
    if target in ALLOCATOR_TARGETS:
        # Complete exact tensor and preimage oracles create millions of temporary integers.
        # Bound the freed-block quarantine, keeping the 2 GiB process gate.
        # This intentionally shortens the use-after-free detection window;
        # other sanitizer options and all other targets remain unchanged.
        options=env.get('ASAN_OPTIONS','')
        env['ASAN_OPTIONS']=options+(':' if options else '')+'quarantine_size_mb=64'
    if target in SHORT_STACK_TARGETS:
        options=env.get('ASAN_OPTIONS','')
        env['ASAN_OPTIONS']=options+(':' if options else '')+'malloc_context_size=5'
    return env


def seed_corpus(target):
    corpus = FUZZ/'corpus'/target
    corpus.mkdir(parents=True, exist_ok=True)

    def save(data):
        path = corpus/hashlib.sha256(data).hexdigest()
        if not path.exists():
            path.write_bytes(data)

    if target == 'step':
        # Every authored STEP fixture (a first byte of 2 or more, the
        # fixture, then its unmutated import or up to six mutations of its
        # instances), a raw file and a raw data section.
        for pick in range(22):
            for k in range(4):
                save(bytes([2, pick])+bytes((j*61+k*37+pick*11+5)%256 for j in range(48)))
        save(b"\x00ISO-10303-21;HEADER;FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));ENDSEC;DATA;"
             b"#1=CARTESIAN_POINT('',(0.,0.,0.));ENDSEC;END-ISO-10303-21;")
        save(b"\x01#1=MANIFOLD_SOLID_BREP('',#2);#2=CLOSED_SHELL('',());")
        save(bytes([0]))
    if target == 'split':
        # Every profile kind and plane mode (parallel, through a vertex,
        # along an edge, tangent, normal, in a cap, oblique, oblique through
        # a vertex, touching a cap's circle), both frames.
        for kind in range(6):
            for mode in range(9):
                for frame in range(2):
                    data = bytearray((j*37+kind*11+mode*17+frame*5+13)%256 for j in range(24))
                    data[:2] = bytes([kind, mode])
                    data[10] = frame
                    save(bytes(data))
        # S8b.3: the spline profiles (a first byte of 192): a bulge, a wave
        # and a lens hole given either way round, in every plane mode (the
        # tangent ones at a spline's apex), both frames.
        for kind in range(3):
            for mode in range(9):
                for frame in range(2):
                    for turn in range(2 if kind == 2 else 1):
                        data = bytearray((j*43+kind*17+mode*23+frame*3+turn*29+11)%256 for j in range(24))
                        data[:3] = bytes([192, kind, mode])
                        data[4] = 4 + 16*turn
                        data[11] = frame
                        save(bytes(data))
        # S8e: sheets and closed wires (a first byte of 160) of every line,
        # arc and spline profile, in every plane mode, both frames.
        for kind in range(6):
            for mode in range(8):
                for flags in range(8):
                    if flags & 1 and kind >= 3:
                        continue
                    data = bytearray((j*31+kind*7+mode*13+flags*19+3)%256 for j in range(24))
                    data[:4] = bytes([160, kind, mode, flags])
                    save(bytes(data))
        # S8c and S8d: a cone, a sphere, a zone, a cap and a torus (a first
        # byte of 224 on) in every plane mode (S8d.2's touching a rim,
        # parallel to a ruling and through a virtual apex among them), both
        # frames.
        for kind in range(5):
            for mode in range(9):
                for frame in range(2):
                    data = bytearray((j*41+kind*13+mode*19+frame*7+5)%256 for j in range(24))
                    data[:3] = bytes([224, kind, mode])
                    data[9] = frame
                    save(bytes(data))
        save(bytes([0]))
    if target == 'boolean':
        # S9a: every pair of line and arc profile kinds, both frames, the
        # four height relations (equal, spanning, overlapping, disjoint).
        for ka in range(6):
            for kb in range(6):
                for flags in range(8):
                    data = bytearray((j*29+ka*13+kb*7+flags*17+1)%256 for j in range(16))
                    data[:3] = bytes([ka, kb, flags])
                    save(bytes(data))
        save(bytes([0]))
    if target == 'curve_curve':
        # Every pair of kinds (line, circle, ellipse, hyperbola) in every mode
        # (independent, one plane, coincident, touching circles, a line
        # along the first's axis).
        for k1 in range(4):
            for k2 in range(4):
                for mode in range(5):
                    data = bytearray((j*47+k1*19+k2*13+mode*7+3)%256 for j in range(24))
                    data[:3] = bytes([k1, k2, mode])
                    save(bytes(data))
        # S7d.2: a spline (a first kind byte of 128) against every conic,
        # degree and mode.
        for k2 in range(12):
            for mode in range(3):
                data = bytearray((j*41+k2*23+mode*5+9)%256 for j in range(64))
                data[:3] = bytes([128, k2, mode])
                save(bytes(data))
        save(bytes([0]))
    if target == 'curve_surface':
        # Every surface kind (plane, cylinder, cone, sphere, torus) against a
        # line and a circle in every mode (independent, along the axis at the
        # radius, across it, through the origin; independent, coaxial,
        # beside, in a plane through the axis).
        for kind in (0, 1, 2, 3, 224):
            for shape in range(2):
                for mode in range(4):
                    data = bytearray((j*53+kind*13+shape*29+mode*5+7)%256 for j in range(32))
                    data[:3] = bytes([kind, shape, mode])
                    save(bytes(data))
        # S7c.2: an ellipse, a hyperbola and a spline (shape bytes 128 to
        # 130) in every mode.
        for kind in (0, 1, 2, 3, 224):
            for shape in (128, 129, 130):
                for mode in range(4):
                    data = bytearray((j*59+kind*17+shape*23+mode*11+5)%256 for j in range(48))
                    data[:3] = bytes([kind, shape, mode])
                    save(bytes(data))
        save(bytes([0]))
    if target == 'analytic_intersections':
        # Every kind pair (plane, cylinder, cone, sphere) in every mode
        # (independent, parallel, coaxial, tangent offset, shared origin).
        for k1 in range(4):
            for k2 in range(4):
                for mode in range(5):
                    data = bytearray((j*43+k1*11+k2*7+mode*3+1)%256 for j in range(24))
                    data[:3] = bytes([k1, k2, mode])
                    save(bytes(data))
        # A torus (a kind byte of 224 or more, S7b.3a) with every kind and
        # another torus, in every mode.
        for k1 in (0, 1, 2, 3, 224):
            for mode in range(5):
                data = bytearray((j*41+k1*13+mode*5+3)%256 for j in range(24))
                data[:3] = bytes([k1, 240, mode])
                save(bytes(data))
    if target == 'brep_io':
        # A prism of the identity structure or an upstream file, then up to
        # six token and line mutations of its .brep text.
        for k in range(96):
            data=bytearray((j*67+k*29+13)%256 for j in range(240))
            save(bytes(data))
        save(bytes([0]))
        # The spline fixtures (S4e): 'S', a fixture, then token and line
        # mutations of their B-spline records.
        for pick in range(6):
            for k in range(4):
                save(bytes([0x53, pick])+bytes((j*53+k*31+7)%256 for j in range(64)))
        # Free shapes (S6): 'F', a sheet, wire or acorn fixture or a
        # profile's face or wire body, then mutations of its text.
        for k in range(24):
            save(bytes([0x46])+bytes((j*59+k*37+11)%256 for j in range(200)))
    elif target == 'attributes':
        # The identity structure, six random policies, random keys and small
        # values, a split height and changes on one piece before the fuse.
        for k in range(96):
            data=bytearray((j*71+k*43+3)%256 for j in range(400))
            save(bytes(data))
        save(bytes([0]))
    elif target == 'split_merge':
        # The identity structure, a split height, one of six history
        # mutations and the stacked fuse of an independent construction.
        for k in range(96):
            data=bytearray((j*83+k*41+5)%256 for j in range(200))
            save(bytes(data))
        save(bytes([0]))
    elif target == 'history':
        # The identity structure plus one of seven history mutations.
        for k in range(112):
            data=bytearray((j*89+k*37+11)%256 for j in range(176))
            save(bytes(data))
        save(bytes([0]))
    elif target == 'tessellation':
        # Every body kind (polygon prism, arc prism, arc face, cone, sphere,
        # torus, polygon face; T-b: spline prism, spline sheet), each with
        # varied shapes, frames, motions, deflections and angles.
        for kind in range(9):
            for k in range(12):
                save(bytes([kind])+bytes((j*61+k*47+kind*19+3)%256 for j in range(96)))
        save(bytes([0]))
    elif target == 'identity':
        # Polygon sides, holes, labels and transforms all vary with the bytes.
        for k in range(96):
            data=bytearray((j*97+k*31+7)%256 for j in range(160))
            save(bytes(data))
        save(bytes([0]))
    elif target == 'brep_validation':
        # S6: 'W', a sheet, shell, wire or acorn model, a scale and a shift.
        for k in range(18):
            save(bytes([0x57, k])+bytes((j*41+k*23+5)%256 for j in range(8)))
        # Every mutation on every base feature (none, round hole, square hole,
        # cavity), outline sizes 3..12 and small, unit and large scales.
        for mutation in range(24):
            for feature in range(4):
                for scale in [0,10,20]:
                    data=bytearray((j*37+mutation*13+feature*5+scale+1)%256 for j in range(40))
                    data[:4]=bytes([mutation,(mutation+feature)%10,feature,scale])
                    save(bytes(data))
        # Cones (mutation 28): an apex at the top or the base, or a frustum,
        # with each cone mutation, at small, unit and large scales.
        for radii in [[1,200,0],[0,1,150],[1,100,1,60]]:
            for sub in range(6):
                for scale in [0,10,20]:
                    save(bytes([28,scale,*radii,120,140,90,200,128,100,160,sub])+bytes((j*37+1)%256 for j in range(16)))
        # Spheres (mutation 29): whole, a pole at either end, or a zone, with
        # each sphere mutation.
        for shape in range(4):
            for sub in range(6):
                for scale in [0,10,20]:
                    save(bytes([29,scale,90,shape,150,140,90,200,128,100,160,sub])+bytes((j*37+1)%256 for j in range(16)))
        # Tori (mutation 30): whole, a wedge, outer and inner segments, with
        # each torus mutation.
        for shape in range(4):
            for sub in range(5):
                for scale in [0,10,20]:
                    extra=[[],[150],[150,140],[150,140]][shape]
                    save(bytes([30,scale,90,120,shape]+extra+[140,90,200,128,100,160,sub])+bytes((j*37+1)%256 for j in range(16)))
        # Spline fixtures (mutation 32, S4b-d): both prisms and both
        # stadiums, three scales, each mutation.
        for which in range(4):
            for scale in (0, 10, 20):
                for kind in range(3):
                    save(bytes([32, which, scale, 7, 200, 129, kind]))
        # Splines (mutation 31, R4): an edge, a pcurve or a plane of a star
        # prism, with or without a round hole, at degrees 2 and 3.
        for c in (0, 5):
            for kind in range(3):
                for degree in range(2):
                    for scale in [0,10,20]:
                        for hole in range(2):
                            save(bytes([31,c,scale]+[128]*(3+c)+[hole,120,140,90,100,150,200,90,160,degree,kind,3,1,0,2])
                                 + bytes((j*37+1)%256 for j in range(16)))
        save(bytes([0]))
    elif target == 'spline_linear':
        # Known factors, rational polylines, the rational circle and the
        # quadratic retrace, as lines and segments, across the degree range.
        for family in range(4):
            for line in [0,1]:
                for variant in range(12):
                    data=bytearray((j*37+variant*13+family*5+1)%256 for j in range(128))
                    data[:2]=bytes([family,line])
                    if family==0:
                        data[2]=[0,1,7,24][variant%4]
                    save(bytes(data))
            save(bytes([family]))
    elif target == 'spline_proximity':
        for degree in range(5,26):
            data=bytearray((i*37+1)%256 for i in range(72))
            data[0:4]=bytes([0,degree-5,0,0])
            data[28:37]=bytes([13,7,9,0,8,16,degree%3,degree%4,degree%3])
            save(bytes(data))
        # Keep the original high-degree positive-distance timeout regressions.
        for degree in [5,6,7,8,11,16,24,25]:
            for variant in [0,1]:
                data=bytearray((i*37+variant*13)%256 for i in range(72))
                data[:3]=bytes([0,degree-5,1])
                data[35:37]=bytes([0,0])
                save(bytes(data))
        # Nonplanar rational cylinder curves have known nonzero minima. Their
        # affine-hull bound is unattainable, retaining stationary-solver work.
        for square in range(7):
            for variant in range(4):
                data=bytearray((i*37+variant*13)%256 for i in range(40))
                data[:3]=bytes([4,square,variant*5])
                data[32:37]=bytes([square*2,16-square*2,square%3,variant,(square+variant)%4])
                save(bytes(data))
        for family in [1,2,3]:
            for mode in [0,1,2,3]:
                for exponent in [0,1,70,127]:
                    data=bytearray((i*17+3)%256 for i in range(40))
                    data[:5]=bytes([family,mode,mode,exponent,127])
                    save(bytes(data))
        for family in range(5):
            save(bytes([family]))
    elif target == 'degree_elevation':
        import struct
        def degree_seed(family,du,dv,ku,kv,qu,qv,op=0,mode=2,scale=128,parameter=7):
            body=(b''.join(struct.pack('<d',x) for _ in range(100) for x in [1.,0.,2.,1.])
                  if mode==0 else bytes((j*37+1)%256 for j in range(3200)))
            save(bytes([family,mode,du-1,dv-1,ku,kv,qu,qv,op,scale,parameter,parameter,0,0,84,1])+body)
        # Every degree and each axis kind; single and both-axis requests.
        for degree in range(1,26):
            degree_seed(0,degree,1,degree%5,0,2,0,degree%3)
            degree_seed(1,degree,2,degree%5,(degree+1)%5,2,2,degree%3)
        for kind in range(5):
            for degree in [1,2,8,24]:
                degree_seed(0,degree,1,kind,0,31,0,degree%3)
            degree_seed(1,24,24,kind,kind,2,2,kind%3)
        for mode in [0,1,3]:
            for scale in [0,1,255]:
                degree_seed(0,2,2,1,3,2,2,1,mode,scale,scale)
                degree_seed(1,2,3,1,3,2,2,2,mode,scale,scale)
        for request in [0,2,26,27,28]:
            for family in range(3):
                degree_seed(family,2,3,1,2,request,request,request%3)
    elif target == 'surface_knots':
        import struct
        for degree in range(1,26):
            for axis in range(2):
                du,dv=(degree,2) if axis==0 else (2,degree)
                save(bytes([2,du-1,dv-1,degree%3,(degree+1)%3,axis,degree%5,84,degree-1,128,7,7,0,0,0,0])+bytes((j*37+1)%256 for j in range(2800)))
        for kind in range(3):
            for op in [0,1,2,4]:
                save(bytes([2,24,24,kind,kind,0,op,84,1,128,7,7,0,0,0,0])+bytes((j*37+1)%256 for j in range(3200)))
        for axis in range(2):
            for degree in [2,8,25]:
                save(bytes([3,degree-1,degree-1,1,1,axis,3,84,1,128,7,7,0,0,0,0]))
        for kind in range(3):
            for scale in [0,128,255]:
                for mode in [0,1]:
                    body=b''.join(struct.pack('<d',x) for _ in range(36) for x in [1.,0.,2.,1.]) if mode==0 else bytes((j*37+1)%256 for j in range(2800))
                    save(bytes([mode,1,1,kind,kind,0,2,84,1,scale,scale,scale,0,0,0,0])+body)
    elif target == 'exact_spline_intersections':
        for family in range(3):
            for mode in range(6):
                for degree in ([1,3,8,25] if mode in [0,1,4] else [2]):
                    for domain in [0,1,2,3,5,6,7]:
                        # Full interval, plus rational refinement. High-degree,
                        # close-root and unrepresentable cases are seed inputs.
                        header=bytes([family,degree-1,mode,domain,1,2,0,0,7,2,5,255,0,84,0,0])
                        save(header+bytes((j*37+1)%256 for j in range(96)))
        for family in range(3):
            for mode in [0,1,3,4,5]:
                save(bytes([family,7,mode,1,2,2,1,8,1,2,7,255,255,255,0,0])+bytes((j*17+3)%256 for j in range(96)))
    elif target == 'knot_editing':
        import struct
        for degree in [1,2,3,8,25]:
            for kind in range(3):
                for op in range(6):
                    save(bytes([2,degree-1,kind,1,op,128,7,84,degree-1,1,2,0,0,0,0,0])+bytes((j*37+1)%256 for j in range(416)))
            for op in [3,4,5]:
                save(bytes([3,degree-1,1,1,op,128,7,84,degree-1,1,2,0,0,0,0,0]))
        for mode in [0,1]:
            for kind in range(3):
                for scale in [0,128,255]:
                    body=b''.join(struct.pack('<d',x) for _ in range(12) for x in [float.fromhex('0x1.fffffffffffffp+1023'),0.,1.,float.fromhex('0x0.0000000000001p-1022')]) if mode==0 else bytes((j*37+1)%256 for j in range(416))
                    save(bytes([mode,2,kind,1,2,scale,scale,84,2,0,1,0,0,0,0,0])+body)
    elif target == 'surface_editing':
        import struct
        for du,dv in [(1,1),(2,3),(5,4),(25,1),(1,25),(25,25)]:
            for ku,kv in [(0,0),(1,0),(0,1),(1,1),(2,2)]:
                for op in [0,1,2,3,4,5,6,7,8,9,10] if du*dv<=25 else [0,10]:
                    save(bytes([2,du-1,dv-1,ku,kv,op,2,2,32,130,7,7,84,72,85,170])+bytes((j*37+1)%256 for j in range(2800)))
        for mode in [0,1]:
            for kind in range(3):
                for scale in [0,128,255]:
                    body=b''.join(struct.pack('<d',x) for _ in range(25) for x in [1.,0.,2.,1.]) if mode==0 else bytes((j*37+1)%256 for j in range(2800))
                    save(bytes([mode,1,1,kind,kind,10,2,2,192,scale,scale,scale,84,72,85,170])+body)
    elif target == 'bezier_editing':
        import struct
        for degree in [0,1,2,7,24]:
            for kind in range(3):
                for op in range(6):
                    save(bytes([2,degree,1,kind,op,3,1,1,128,15,84,85])+bytes((j*37+1)%256 for j in range(160)))
        for mode in [0,1]:
            for kind in range(3):
                for scale in [0,128,255]:
                    body=b''.join(struct.pack('<d',x) for _ in range(6) for x in [float.fromhex('0x1.fffffffffffffp+1023'),0.,1.,float.fromhex('0x0.0000000000001p-1022')]) if mode==0 else bytes((j*37+1)%256 for j in range(160))
                    save(bytes([mode,2,1,kind,5,3,1,1,scale,scale,84,85])+body)
    elif target in ['proximity','linear_sets']:
        from proximity_reference import COUNTS
        for row in (ROOT/'rust/fixtures'/('proximity.tsv' if target=='proximity' else 'linear_sets.tsv')).read_text().splitlines():
            if row.startswith('#'): continue
            words=iter(row.split()); next(words)
            kinds=[]; coordinates=[]
            for _ in range(2):
                kind=next(words); kinds.append(list(COUNTS).index(kind))
                coordinates += [int(next(words),16).to_bytes(8,'little') for _ in range(COUNTS[kind]*3)]
            save(bytes([0,*kinds,0,128])+b''.join(coordinates))
        for mode in [1,2,3]:
            for a in range(5):
                for b in range(5):
                    save(bytes([mode,a,b,2,128])+bytes((j*37+1)%256 for j in range(144)))
    elif target == 'roots':
        for degree in [0,1,2,7,24]:
            for mode in range(4):
                save(bytes([degree,128])+bytes([mode])*25+bytes([mode,mode,24])+bytes((j*37)%256 for j in range(200)))
    elif target == 'spline_intersections':
        for mode in [48,52,56,60]:
            for degree in [0,12,24]:
                for surface in [0,2]:
                    save(bytes([mode,degree,128,surface,128])+bytes(100))
        for mode in range(32):
            for degree in range(8):
                for value in [0,1,2,3,4,5,6]:
                    save(bytes([mode,degree,128,1,128])+bytes(5)+bytes([value])*8+bytes((j*37+1)%256 for j in range(100)))
    elif target == 'modeling':
        for i in range(48):
            save(bytes([(i*73+j*37) % 255 for j in range(121)]))
        save(bytes([255])+bytes(120))
    elif target == 'splines':
        for mode in range(4):
            for degree in [0,1,2,7,24]:
                for parameter in range(5):
                    save(bytes([mode,degree,1,0,2,128,parameter,0])+bytes((j*37+1)%256 for j in range(160)))
                    save(bytes([mode,degree,1,0,2,128,parameter,1])+bytes((j*37+1)%256 for j in range(160)))
        for word in [0x3ff0000000000000,0x7fefffffffffffff,1,0x7ff0000000000000,0x7ff8000000000000]:
            save(bytes([1,1,1,0,2,128,2,0])+word.to_bytes(8,'little')*20)
    elif target == 'surfaces':
        for mode in range(4):
            for degree in [0,1,2,7,24]:
                for periodic in range(4):
                    for parameter in [0,2,8]:
                        save(bytes([mode,degree,1,periodic,0,0,2,128,parameter,parameter,0,0])+bytes((j*37+1)%256 for j in range(200)))
        for word in [0x3ff0000000000000,0x7fefffffffffffff,1,0x7ff0000000000000,0x7ff8000000000000]:
            save(bytes([1,1,1,3,0,0,2,128,2,2,0,0])+word.to_bytes(8,'little')*25)
    elif target == 'curved':
        for source in ['quadratic.tsv','curved.tsv']:
            rows = [r.split() for r in (ROOT/'rust/fixtures'/source).read_text().splitlines() if not r.startswith('#')]
            for i,row in enumerate(rows):
                if i < 112 or i % 41 == 0:
                    kind = 0 if source == 'quadratic.tsv' else {'s':1,'y':2,'c':3}[row[0].lower()]
                    words = row[1:4] if kind == 0 else row[2:15]
                    save(bytes([0,kind])+b''.join(int(x,16).to_bytes(8,'little') for x in words))
        for mode in [1,2,3]:
            for kind in range(4):
                save(bytes([mode,kind])+bytes((j*37)%256 for j in range(106)))
    else:
        source = 'predicates3d.tsv' if target == 'predicates' else 'intersections.tsv'
        rows = [r.split() for r in (ROOT/'rust/fixtures'/source).read_text().splitlines() if not r.startswith('#')]
        for i,row in enumerate(rows):
            if i < 64 or i % 29 == 0:
                words = row[2:17] if row[0] != 'o' else row[2:14]+['0000000000000000']*3
                save(b'\0'+b''.join(int(x,16).to_bytes(8,'little') for x in words))
        for mode in [1,2,3]:
            save(bytes([mode])+bytes((j*37)%256 for j in range(120)))
        for word in [0x7ff0000000000000,0xfff0000000000000,0x7ff8000000000000]:
            save(b'\0'+word.to_bytes(8,'little')*15)
    regressions = FUZZ/'regressions'/target
    if regressions.exists():
        for path in regressions.glob('*.bin'):
            save(path.read_bytes())
    return corpus


# U6 of REVIEW_NOTES.md: a per-push run replays every checked-in regression,
# every input added since the last full (schedule) replay and a seeded random
# sample of this many others; the schedule replays everything.
SAMPLE_SIZE = 64
# Exact tensor targets whose replay dominates CI (F7): per push they replay
# their regressions only; they fuzz on the schedule.
SCHEDULE_ONLY_TARGETS = {'surface_knots', 'degree_elevation', 'surface_editing'}
# Targets whose scheduled full replay is split across this many jobs, each
# replaying (-runs=0) the inputs whose contents hash to its shard under its
# own startup budget; a check job then requires the shards' union to be the
# snapshot they were given, each input exactly once (REVIEW_NOTES parallel
# track "The boolean target's full replay"). Their scheduled campaign replays
# regressions and a seeded sample before its mutation. The rust-fuzz.yml
# replay job's matrix must list the same shards (test_fuzz_runner.py).
# boolean: 356 CI inputs replayed in 1,830 s after a 275 s build at
# 2efa1ec7, up from 1,242 s of startup at 85104dc3 on the same inputs.
# split: CI's 1,749 inputs overran the hour at c3ce4d42; eight shards of
# them take 250 to 413 s of CPU each under ASan on the development Mac.
REPLAY_SHARDS = {'boolean': 4, 'degree_elevation': 4, 'split': 8}


def manifest_path(target):
    # Beside the corpus, not in it: libFuzzer reads every file of a corpus.
    return FUZZ/'corpus'/f'{target}.replayed'


def regression_names(target):
    regressions = FUZZ/'regressions'/target
    if not regressions.exists():
        return set()
    return {hashlib.sha256(path.read_bytes()).hexdigest() for path in regressions.glob('*.bin')}


def write_manifest(target, corpus):
    manifest_path(target).write_text(''.join(p.name+'\n' for p in sorted(corpus.iterdir())))


def replay_plan(target, corpus, sample_seed, size=SAMPLE_SIZE):
    """(corpus file names to replay, evidence). Without a manifest every
    input is new, so the plan is the whole corpus."""
    names = sorted(p.name for p in corpus.iterdir())
    manifest = manifest_path(target)
    if not manifest.exists():
        return names, {'replay':'sample','manifest':False,'sample_seed':sample_seed,
                       'regressions':0,'new_since_full_replay':len(names),'sampled':0}
    replayed = set(manifest.read_text().split())
    regressions = regression_names(target) & set(names)
    new = [n for n in names if n not in replayed]
    rest = [n for n in names if n in replayed and n not in regressions]
    sampled = random.Random(sample_seed).sample(rest, min(size, len(rest)))
    chosen = sorted(regressions | set(new) | set(sampled))
    return chosen, {'replay':'sample','manifest':True,'sample_seed':sample_seed,
                    'regressions':len(regressions),'new_since_full_replay':len(new),'sampled':len(sampled)}


def sharded_campaign_plan(target, corpus, sample_seed, size=SAMPLE_SIZE):
    """The scheduled campaign of a REPLAY_SHARDS target: the shard jobs of
    the same run replay every input, so it replays the checked-in
    regressions and a seeded sample, whatever the manifest says."""
    names = sorted(p.name for p in corpus.iterdir())
    regressions = regression_names(target) & set(names)
    rest = [n for n in names if n not in regressions]
    sampled = random.Random(sample_seed).sample(rest, min(size, len(rest)))
    return sorted(regressions | set(sampled)), {'replay':'sharded','shards':REPLAY_SHARDS[target],
        'sample_seed':sample_seed,'regressions':len(regressions),'sampled':len(sampled)}


def shard_of(data, shards):
    # By contents, not by name: a seed and libFuzzer's copy of the same
    # bytes land in one shard, and the assignment needs no other input.
    return int.from_bytes(hashlib.sha256(data).digest()[:8],'big') % shards


def shard_names(corpus, shards):
    """Every corpus file's name in exactly one of `shards` sorted lists."""
    plan = [[] for _ in range(shards)]
    for path in sorted(corpus.iterdir()):
        plan[shard_of(path.read_bytes(),shards)].append(path.name)
    return plan


def listing_digest(names):
    return hashlib.sha256(''.join(n+'\n' for n in sorted(names)).encode()).hexdigest()


def replay_shard_passed(report):
    # INITED is printed only after libFuzzer has run every corpus file, each
    # once plus leak-check reruns, and the empty input once more.
    return (report['exit_code'] == 0 and report['initial_executions'] is not None
            and report['initial_executions'] >= report['nonempty_files']+1
            and resource_limits_satisfied(report))


def replay_shard(target, toolchain, corpus, shard, output, env, shards=None):
    """Replay (-runs=0) one shard of the corpus with the target's limits."""
    shards = shards or REPLAY_SHARDS[target]
    plan = shard_names(corpus, shards)
    names = plan[shard]
    input_seconds = TARGET_INPUT_SECONDS.get(target,INPUT_SECONDS)
    artifacts = FUZZ/'artifacts'/target
    artifacts.mkdir(parents=True,exist_ok=True)
    run_corpus = corpus.parent/f'.{target}-shard-{shard}'
    if run_corpus.exists(): shutil.rmtree(run_corpus)
    run_corpus.mkdir()
    for name in names:
        shutil.copy2(corpus/name,run_corpus/name)
    budget = startup_budget(len(names),input_seconds,TARGET_MAX_STARTUP_SECONDS.get(target,MAX_STARTUP_SECONDS))
    command = ['cargo',f'+{toolchain}','fuzz','run',target,str(run_corpus),'--fuzz-dir',str(FUZZ),
               *sanitizer_build_args(target),'--','-runs=0',f'-timeout={input_seconds}','-rss_limit_mb=2048',
               f'-max_len={max_len(target)}',f'-artifact_prefix={artifacts}/','-print_final_stats=1']
    log_path = output/f'{target}-shard-{shard}.log'
    campaign_env = campaign_environment(target,env)
    started = time.monotonic()
    try:
        with log_path.open('w') as log:
            code = run_process(command,log,budget,campaign_env)
    finally:
        shutil.rmtree(run_corpus)
    report = {'target':target,'mode':'replay-shard','shard':shard,'shards':shards,
              'corpus_files':sum(map(len,plan)),'corpus_digest':listing_digest(n for part in plan for n in part),
              'shard_files':len(names),'nonempty_files':sum(1 for n in names if (corpus/n).stat().st_size),
              'exit_code':code,'elapsed_seconds':round(time.monotonic()-started,2),'budget_seconds':budget,
              'input_limit_seconds':input_seconds,'sanitizer_options':campaign_env.get('ASAN_OPTIONS'),
              **statistics(log_path.read_text(errors='replace')),
              'artifacts':[p.name for p in sorted(artifacts.iterdir())],'command':command,'names':names}
    report['passed'] = replay_shard_passed(report)
    return report


def check_replay_shards(target, corpus, reports, shards=None):
    """Accept the sharded full replay only if every shard of this snapshot
    reported, each replayed exactly its inputs and passed. On success the
    manifest records the snapshot, as a full replay's does."""
    shards = shards or REPLAY_SHARDS[target]
    plan = shard_names(corpus, shards)
    everything = {n for part in plan for n in part}
    digest = listing_digest(everything)
    problems, by_shard = [], {}
    for report in reports:
        k = report.get('shard')
        if report.get('target') != target or report.get('shards') != shards or k not in range(shards):
            problems.append(f'report of {report.get("target")} shard {k} of {report.get("shards")} is not one of {shards} {target} shards')
        elif k in by_shard:
            problems.append(f'shard {k} reported twice')
        else:
            by_shard[k] = report
    missing = [k for k in range(shards) if k not in by_shard]
    problems += [f'shard {k} did not report' for k in missing]
    counts = {}
    for report in by_shard.values():
        for name in report['names']:
            counts[name] = counts.get(name,0)+1
    unreplayed = sorted(everything-counts.keys())
    unexpected = sorted(counts.keys()-everything)
    duplicated = sorted(n for n,c in counts.items() if c > 1)
    for k,report in sorted(by_shard.items()):
        if report.get('corpus_digest') != digest:
            problems.append(f'shard {k} replayed another snapshot')
        if report['names'] != plan[k]:
            problems.append(f'shard {k} replayed other inputs than its own')
    failed = sorted(k for k,report in by_shard.items() if not report.get('passed'))
    complete = not problems and not unreplayed and not unexpected and not duplicated
    passed = complete and not failed
    if passed:
        write_manifest(target,corpus)
    def most(key):
        values = [r[key] for r in by_shard.values() if r.get(key) is not None]
        return max(values) if values else None
    return {'target':target,'mode':'replay-shards','shards':shards,'corpus_files':len(everything),
            'corpus_digest':digest,'replayed_files':sum(counts.values()),
            'shard_files':[len(by_shard[k]['names']) if k in by_shard else None for k in range(shards)],
            'shard_elapsed_seconds':[by_shard[k].get('elapsed_seconds') if k in by_shard else None for k in range(shards)],
            'slowest_input_seconds':most('slowest_input_seconds'),'peak_rss_mb':most('peak_rss_mb'),
            'problems':problems,'unreplayed':unreplayed,'unexpected':unexpected,'duplicated':duplicated,
            'failed_shards':failed,'complete':complete,'passed':passed,'manifest_written':passed}


def version(command, cwd=ROOT):
    return subprocess.check_output(command,cwd=cwd,text=True,stderr=subprocess.STDOUT).strip()


def run_process(command, log, timeout, env, tick=None, phase_deadline=None):
    # Kill the entire build/fuzzer process group on a wall-clock timeout.
    process = subprocess.Popen(command,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,
                               env=env,start_new_session=True)
    deadline=time.monotonic()+timeout
    try:
        while True:
            code=process.poll()
            if code is not None: return code
            if tick is not None: tick()
            active_deadline=min(deadline,phase_deadline()) if phase_deadline is not None else deadline
            remaining=active_deadline-time.monotonic()
            if remaining<=0: raise subprocess.TimeoutExpired(command,timeout)
            try:
                return process.wait(timeout=min(.1,remaining) if tick is not None else remaining)
            except subprocess.TimeoutExpired:
                if time.monotonic()>=active_deadline: raise
    except subprocess.TimeoutExpired:
        os.killpg(process.pid,signal.SIGKILL)
        process.wait()
        return 124
    except BaseException:
        os.killpg(process.pid,signal.SIGKILL)
        process.wait()
        raise


class MutationBudget:
    """Start the mutation clock only after libFuzzer finishes corpus replay.

    libFuzzer's max_total_time includes initialization. Its supported stop_file
    flag instead exits the normal fuzz loop and prints final statistics. The
    pinned runtime requires a NONEMPTY stop file (FuzzerLoop.cpp).
    """
    def __init__(self, log_path, stop_file, seconds, startup_seconds=STARTUP_SECONDS, shutdown_seconds=SHUTDOWN_SECONDS):
        self.log_path,self.stop_file,self.seconds=log_path,stop_file,seconds
        self.startup_seconds,self.shutdown_seconds=startup_seconds,shutdown_seconds
        self.started=time.monotonic()
        self.initialized=None
        self.stopped=None
        self.offset=0
        self.pending=''

    def tick(self):
        if self.initialized is None:
            with self.log_path.open(errors='replace') as reader:
                reader.seek(self.offset)
                self.pending+=reader.read()
                self.offset=reader.tell()
            if re.search(r'#\d+\s+INITED\b',self.pending): self.initialized=time.monotonic()
            self.pending=self.pending[-512:]
        if self.initialized is not None and self.stopped is None and time.monotonic()-self.initialized>=self.seconds:
            self.stop_file.write_bytes(b'stop\n')
            self.stopped=time.monotonic()

    def evidence(self):
        return {'startup_seconds':round(self.initialized-self.started,2) if self.initialized is not None else None,
                'startup_budget_completed':self.initialized is not None and self.initialized-self.started<=self.startup_seconds,
                'startup_limit_seconds':self.startup_seconds,
                'mutation_seconds_before_stop':round(self.stopped-self.initialized,2) if self.stopped is not None else None,
                'mutation_budget_completed':self.stopped is not None,
                'shutdown_seconds':round(time.monotonic()-self.stopped,2) if self.stopped is not None else None,
                'shutdown_grace_seconds':self.shutdown_seconds}

    def deadline(self):
        # Never let a late INITED marker turn a startup overrun into success.
        if self.initialized is None or self.initialized-self.started>self.startup_seconds:
            return self.started+self.startup_seconds
        # stop_file is checked between batches of MUTATION_DEPTH callbacks.
        # Reserve that batch plus log/exit time, without extending startup or
        # changing any individual callback's timeout/resource validation.
        return self.initialized+self.seconds+self.shutdown_seconds


def statistics(text):
    executed = re.findall(r'stat::number_of_executed_units:\s*(\d+)',text)
    initialized = re.findall(r'#(\d+)\s+INITED\b',text)
    coverage = re.findall(r'cov: (\d+)',text)
    slowest = re.findall(r'stat::slowest_unit_time_sec:\s*(\d+)',text)
    peak = re.findall(r'stat::peak_rss_mb:\s*(\d+)',text)
    total = int(executed[-1]) if executed else None
    replayed = int(initialized[-1]) if initialized else None
    return {'executions':total,'initial_executions':replayed,
            'mutation_executions':total-replayed if total is not None and replayed is not None else None,
            'coverage_edges':int(coverage[-1]) if coverage else None,
            'slowest_input_seconds':int(slowest[-1]) if slowest else None,
            'peak_rss_mb':int(peak[-1]) if peak else None}


def max_len(target):
    return 4096 if target in ["surface_editing", "surface_knots", "degree_elevation"] else 512 if target == "knot_editing" else 256


def merge_command(target, toolchain, fresh, corpus, artifacts):
    # libFuzzer's -merge=1 keeps, in `fresh`, a subset of `corpus` with the
    # same coverage (REVIEW_NOTES R12); each input keeps its own limits.
    input_seconds=TARGET_INPUT_SECONDS.get(target,INPUT_SECONDS)
    return ['cargo',f'+{toolchain}','fuzz','run',target,str(fresh),str(corpus),
            '--fuzz-dir',str(FUZZ),*sanitizer_build_args(target),'--',
            '-merge=1',f'-timeout={input_seconds}','-rss_limit_mb=2048',f'-max_len={max_len(target)}',
            f'-artifact_prefix={artifacts}/','-print_final_stats=1']


def minimize(target, toolchain, corpus, output, env):
    """Merge `corpus` into a fresh directory and replace it only when the
    merge completed: an interrupted or failed merge leaves it untouched."""
    artifacts=FUZZ/'artifacts'/target
    artifacts.mkdir(parents=True,exist_ok=True)
    fresh=corpus.parent/f'.{target}-minimized'
    if fresh.exists(): shutil.rmtree(fresh)
    fresh.mkdir()
    before=len(list(corpus.iterdir()))
    input_seconds=TARGET_INPUT_SECONDS.get(target,INPUT_SECONDS)
    budget=startup_budget(before,input_seconds,TARGET_MAX_STARTUP_SECONDS.get(target,MAX_STARTUP_SECONDS))
    log_path=output/f'{target}-minimize.log'
    started=time.monotonic()
    command=merge_command(target,toolchain,fresh,corpus,artifacts)
    with log_path.open('w') as log:
        code=run_process(command,log,budget,campaign_environment(target,env))
    after=len(list(fresh.iterdir()))
    replaced=code==0 and after>0
    if replaced:
        shutil.rmtree(corpus)
        fresh.rename(corpus)
        # The merge ran every kept input: a full replay.
        write_manifest(target,corpus)
    else:
        shutil.rmtree(fresh)
    return {'target':target,'exit_code':code,'elapsed_seconds':round(time.monotonic()-started,2),
            'budget_seconds':budget,'corpus_files_before':before,'corpus_files_after':after if replaced else before,
            'replaced':replaced,'artifacts':[p.name for p in sorted(artifacts.iterdir())],'command':command}


def replay_regressions(target, toolchain, output, env):
    """Run every checked-in regression once, with the target's limits."""
    files=sorted((FUZZ/'regressions'/target).glob('*.bin')) if (FUZZ/'regressions'/target).exists() else []
    report={'target':target,'mode':'regressions','regressions':len(files)}
    if not files:
        return dict(report,exit_code=0)
    input_seconds=TARGET_INPUT_SECONDS.get(target,INPUT_SECONDS)
    artifacts=FUZZ/'artifacts'/target
    artifacts.mkdir(parents=True,exist_ok=True)
    command=['cargo',f'+{toolchain}','fuzz','run',target,*map(str,files),'--fuzz-dir',str(FUZZ),
             *sanitizer_build_args(target),'--',f'-timeout={input_seconds}','-rss_limit_mb=2048',
             f'-artifact_prefix={artifacts}/']
    log_path=output/f'{target}-regressions.log'
    started=time.monotonic()
    with log_path.open('w') as log:
        code=run_process(command,log,STARTUP_SECONDS+input_seconds*(len(files)+1),campaign_environment(target,env))
    return dict(report,exit_code=code,elapsed_seconds=round(time.monotonic()-started,2),command=command)


def resource_limits_satisfied(report):
    # The runtime alarm is asynchronous: a callback can exceed its configured
    # limit and still return success between alarm ticks. Check final evidence.
    slowest,peak=report['slowest_input_seconds'],report['peak_rss_mb']
    return (slowest is not None and peak is not None and
            slowest <= report['input_limit_seconds'] and peak <= 2048)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--target',choices=TARGETS,action='append')
    parser.add_argument('--seconds',type=int,default=300,help='mutation budget per target (1..3600)')
    parser.add_argument('--toolchain',default='nightly-2026-09-22')
    parser.add_argument('--seed',type=int,default=0,help='0 asks libFuzzer for a random seed; the log records it')
    parser.add_argument('--seed-only',action='store_true')
    parser.add_argument('--report-dir',type=Path,default=ROOT/'target/fuzz-reports')
    parser.add_argument('--minimize',action='store_true',
                        help='merge each seeded corpus into a minimal one instead of fuzzing (R12)')
    parser.add_argument('--replay',choices=['full','sample'],default='full',
                        help='replay the whole corpus, or the per-push sample of U6')
    parser.add_argument('--sample-seed',type=int,default=0,help='the per-push sample\'s seed (recorded)')
    parser.add_argument('--regressions-only',action='store_true',
                        help='replay the checked-in regressions and stop (schedule-only targets per push)')
    parser.add_argument('--per-push',action='store_true',
                        help='U6: regressions only for schedule-only targets, a sampled replay otherwise')
    parser.add_argument('--scheduled',action='store_true',
                        help='the schedule\'s campaign: a full replay, or for REPLAY_SHARDS targets '
                             'regressions and a seeded sample (their shard jobs replay everything)')
    parser.add_argument('--replay-shard',type=int,metavar='K',
                        help='replay shard K of the unseeded corpus with -runs=0 and stop (REPLAY_SHARDS)')
    parser.add_argument('--check-replay-shards',type=Path,metavar='DIR',
                        help='check the shard reports under DIR against the unseeded corpus; '
                             'write the manifest if they cover it exactly once and passed')
    args = parser.parse_args()
    if args.per_push:
        if set(args.target or TARGETS) <= SCHEDULE_ONLY_TARGETS:
            args.regressions_only = True
        else:
            args.replay = 'sample'
    if not 1 <= args.seconds <= 3600: parser.error('--seconds must be in [1,3600]')
    targets = args.target or TARGETS
    sharding = args.replay_shard is not None or args.check_replay_shards is not None
    if sharding:
        if len(targets) != 1 or targets[0] not in REPLAY_SHARDS:
            parser.error(f'shard replays take one --target of {sorted(REPLAY_SHARDS)}')
        if args.replay_shard is not None and not 0 <= args.replay_shard < REPLAY_SHARDS[targets[0]]:
            parser.error(f'--replay-shard must be in [0,{REPLAY_SHARDS[targets[0]]})')
    # A shard replays, and the check reads, the snapshot it was given as it is.
    corpora = {target: FUZZ/'corpus'/target if sharding else seed_corpus(target) for target in targets}
    if args.seed_only: return
    if args.check_replay_shards is not None:
        target = targets[0]
        reports = [json.loads(p.read_text()) for p in sorted(args.check_replay_shards.rglob(f'{target}-shard-*.json'))]
        report = check_replay_shards(target,corpora[target],reports)
        print(json.dumps(report),flush=True)
        for k,(files,seconds) in enumerate(zip(report['shard_files'],report['shard_elapsed_seconds'])):
            print(f'{target} shard {k}: {files} inputs, {seconds} s',flush=True)
        output = args.report_dir.resolve()
        output.mkdir(parents=True,exist_ok=True)
        (output/f'{target}-replay-shards.json').write_text(json.dumps(report,indent=2)+'\n')
        raise SystemExit(0 if report['passed'] else 1)
    output = args.report_dir.resolve()
    output.mkdir(parents=True,exist_ok=True)
    summary = {'toolchain':args.toolchain,'seconds_per_target':args.seconds,'targets':[], 'platform':platform.platform()}
    failed = False
    try:
        summary.update(rustc=version(['rustup','run',args.toolchain,'rustc','--version']),
            cargo_fuzz=version(['cargo','fuzz','--version']),revision=version(['git','rev-parse','HEAD']),
            dirty=bool(version(['git','status','--porcelain'])))
        # cargo-fuzz 0.13.1 has no --locked option. Validate/fetch the lock first,
        # build offline, then fail if cargo-fuzz changes that lock.
        version(['cargo','fetch','--locked','--manifest-path',str(FUZZ/'Cargo.toml')])
        locked = (FUZZ/'Cargo.lock').read_bytes()
        summary['lock_sha256'] = hashlib.sha256(locked).hexdigest()
        env = dict(os.environ,CARGO_NET_OFFLINE='true')
        if args.minimize:
            summary['mode']='minimize'
            for target in targets:
                print(f'Minimizing the {target} corpus',flush=True)
                report=minimize(target,args.toolchain,corpora[target],output,env)
                summary['targets'].append(report)
                failed |= not report['replaced'] or (FUZZ/'Cargo.lock').read_bytes() != locked
                print(json.dumps(report),flush=True)
        if args.regressions_only:
            summary['mode']='regressions'
            for target in targets:
                report=replay_regressions(target,args.toolchain,output,env)
                summary['targets'].append(report)
                failed |= report['exit_code'] != 0 or (FUZZ/'Cargo.lock').read_bytes() != locked
                print(json.dumps(report),flush=True)
        if args.replay_shard is not None:
            summary['mode']='replay-shard'
            target=targets[0]
            report=replay_shard(target,args.toolchain,corpora[target],args.replay_shard,output,env)
            summary['targets'].append({k:v for k,v in report.items() if k != 'names'})
            (output/f'{target}-shard-{args.replay_shard}.json').write_text(json.dumps(report,indent=2)+'\n')
            failed |= not report['passed'] or (FUZZ/'Cargo.lock').read_bytes() != locked
            print(json.dumps({k:v for k,v in report.items() if k != 'names'}),flush=True)
        for target in [] if args.minimize or args.regressions_only or sharding else targets:
            artifacts = FUZZ/'artifacts'/target
            artifacts.mkdir(parents=True,exist_ok=True)
            log_path = output/f'{target}.log'
            # cargo-fuzz enables address sanitizer and coverage instrumentation.
            # Separate startup, mutation and shutdown budgets keep an in-flight
            # final input from consuming the time reserved for corpus startup.
            started = time.monotonic()
            print(f'Fuzzing {target} for {args.seconds}s; log: {log_path}',flush=True)
            with tempfile.TemporaryDirectory(prefix=f'.{target}-control-',dir=output) as control:
                stop_file=Path(control)/'stop'
                input_seconds=TARGET_INPUT_SECONDS.get(target,INPUT_SECONDS)
                shutdown_seconds=shutdown_budget(input_seconds)
                corpus=corpora[target]
                mode='sharded' if args.scheduled and target in REPLAY_SHARDS else args.replay
                if mode in ('sample','sharded'):
                    if mode=='sharded':
                        chosen,replay=sharded_campaign_plan(target,corpus,args.sample_seed)
                    else:
                        chosen,replay=replay_plan(target,corpus,args.sample_seed)
                    run_corpus=corpus.parent/f'.{target}-sample'
                    if run_corpus.exists(): shutil.rmtree(run_corpus)
                    run_corpus.mkdir()
                    # Keep each seed's modification time: libFuzzer rereads
                    # corpus files newer than its first read.
                    for name in chosen:
                        shutil.copy2(corpus/name,run_corpus/name)
                else:
                    run_corpus,replay=corpus,{'replay':'full'}
                initial_corpus_files=len(list(run_corpus.iterdir()))
                startup_seconds=startup_budget(initial_corpus_files,input_seconds,
                                               TARGET_MAX_STARTUP_SECONDS.get(target,MAX_STARTUP_SECONDS))
                timer=MutationBudget(log_path,stop_file,args.seconds,startup_seconds=startup_seconds,shutdown_seconds=shutdown_seconds)
                command = ['cargo',f'+{args.toolchain}','fuzz','run',target,str(run_corpus),
                    '--fuzz-dir',str(FUZZ),*sanitizer_build_args(target),'--',
                    '-max_total_time=0',f'-stop_file={stop_file}',f'-mutate_depth={MUTATION_DEPTH}',f'-timeout={input_seconds}','-rss_limit_mb=2048','-reload=0',
                    f'-max_len={max_len(target)}',f'-seed={args.seed}',f'-artifact_prefix={artifacts}/','-print_final_stats=1']
                with log_path.open('w') as log:
                    campaign_env=campaign_environment(target,env)
                    code = run_process(command,log,startup_seconds+args.seconds+shutdown_seconds,campaign_env,timer.tick,timer.deadline)
                if run_corpus != corpus:
                    # A sampled run's new inputs are new only against its
                    # sample; the retained corpus grows on full replays.
                    shutil.rmtree(run_corpus)
            text = log_path.read_text(errors='replace')
            report = {'target':target,'exit_code':code,'elapsed_seconds':round(time.monotonic()-started,2),
                'input_limit_seconds':input_seconds,
                'mutation_depth':MUTATION_DEPTH,
                'initial_corpus_files':initial_corpus_files,
                **replay,
                'allocator_cleanup_during_replay':target in ALLOCATOR_TARGETS,
                'sanitizer_options':campaign_env.get('ASAN_OPTIONS'),
                **timer.evidence(),
                **statistics(text),
                'corpus_files':len(list(corpora[target].iterdir())),
                'artifacts':[p.name for p in sorted(artifacts.iterdir())], 'command':command}
            summary['targets'].append(report)
            # An exit without a completed campaign is not a successful fuzz run.
            target_failed = code != 0 or not report['startup_budget_completed'] or not report['mutation_budget_completed'] or not report['mutation_executions'] or report['mutation_executions'] < 0 or not resource_limits_satisfied(report) or (FUZZ/'Cargo.lock').read_bytes() != locked
            if mode=='full' and not target_failed:
                write_manifest(target,corpus)
            failed |= target_failed
            print(json.dumps(report),flush=True)
    finally:
        summary['passed'] = len(summary['targets']) == len(targets) and not failed
        (output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    raise SystemExit(1 if failed else 0)


if __name__ == '__main__':
    main()
