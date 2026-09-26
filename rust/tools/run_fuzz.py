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
import re
import signal
import subprocess
import tempfile
import time

ROOT = Path(__file__).resolve().parents[2]
FUZZ = ROOT/'rust/fuzz'
TARGETS = ['predicates','intersections','modeling','curved','splines','surfaces','roots','spline_intersections','proximity','linear_sets','bezier_editing','surface_editing','knot_editing','exact_spline_intersections','surface_knots','degree_elevation','spline_proximity','spline_linear','brep_validation','identity','history','split_merge','attributes','brep_io']
STARTUP_SECONDS = 600
MAX_STARTUP_SECONDS = 3600
INPUT_SECONDS = 20
# Full tensor coefficient equations and double-axis degree-25 edits are a
# larger per-input workload. Existing targets keep their original 20s limit.
# surface_editing joined on measurement: a CI input took 8.2 s under local
# AddressSanitizer, above 20 s / 2.6 (fuzz/regressions/README.md).
TARGET_INPUT_SECONDS = {"surface_knots": 60, "degree_elevation": 60, "surface_editing": 60}
# Targets whose exact oracles churn enough temporary BigInts that default ASan
# quarantine and allocator retention, not live data, exhaust the RSS gate.
ALLOCATOR_TARGETS = {'surface_knots', 'degree_elevation', 'spline_linear'}
# The pinned libFuzzer checks stop_file between MutateAndTestOne batches,
# not between each callback. Keep its default mutation sequence length.
MUTATION_DEPTH = 5
SHUTDOWN_SECONDS = MUTATION_DEPTH * INPUT_SECONDS + 5


def shutdown_budget(input_seconds):
    return MUTATION_DEPTH * input_seconds + 5


def startup_budget(corpus_files, input_seconds=INPUT_SECONDS):
    # Exact-arithmetic seeds have widely varying costs. An assumed average
    # cost can reject an otherwise legal corpus before mutation ever starts.
    # Reserve build time and each input's configured budget (including the
    # initial empty input), capped at one hour. Every input keeps its own
    # timeout; replay never borrows the requested mutation duration.
    return min(MAX_STARTUP_SECONDS, STARTUP_SECONDS+input_seconds*(corpus_files+1))


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
    return env


def seed_corpus(target):
    corpus = FUZZ/'corpus'/target
    corpus.mkdir(parents=True, exist_ok=True)

    def save(data):
        path = corpus/hashlib.sha256(data).hexdigest()
        if not path.exists():
            path.write_bytes(data)

    if target == 'brep_io':
        # A prism of the identity structure or an upstream file, then up to
        # six token and line mutations of its .brep text.
        for k in range(96):
            data=bytearray((j*67+k*29+13)%256 for j in range(240))
            save(bytes(data))
        save(bytes([0]))
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
    elif target == 'identity':
        # Polygon sides, holes, labels and transforms all vary with the bytes.
        for k in range(96):
            data=bytearray((j*97+k*31+7)%256 for j in range(160))
            save(bytes(data))
        save(bytes([0]))
    elif target == 'brep_validation':
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
    args = parser.parse_args()
    if not 1 <= args.seconds <= 3600: parser.error('--seconds must be in [1,3600]')
    targets = args.target or TARGETS
    corpora = {target: seed_corpus(target) for target in targets}
    if args.seed_only: return
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
        for target in targets:
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
                initial_corpus_files=len(list(corpora[target].iterdir()))
                startup_seconds=startup_budget(initial_corpus_files,input_seconds)
                timer=MutationBudget(log_path,stop_file,args.seconds,startup_seconds=startup_seconds,shutdown_seconds=shutdown_seconds)
                command = ['cargo',f'+{args.toolchain}','fuzz','run',target,str(corpora[target]),
                    '--fuzz-dir',str(FUZZ),*sanitizer_build_args(target),'--',
                    '-max_total_time=0',f'-stop_file={stop_file}',f'-mutate_depth={MUTATION_DEPTH}',f'-timeout={input_seconds}','-rss_limit_mb=2048',
                    f'-max_len={4096 if target in ["surface_editing", "surface_knots", "degree_elevation"] else 512 if target == "knot_editing" else 256}',f'-seed={args.seed}',f'-artifact_prefix={artifacts}/','-print_final_stats=1']
                with log_path.open('w') as log:
                    campaign_env=campaign_environment(target,env)
                    code = run_process(command,log,startup_seconds+args.seconds+shutdown_seconds,campaign_env,timer.tick,timer.deadline)
            text = log_path.read_text(errors='replace')
            report = {'target':target,'exit_code':code,'elapsed_seconds':round(time.monotonic()-started,2),
                'input_limit_seconds':input_seconds,
                'mutation_depth':MUTATION_DEPTH,
                'initial_corpus_files':initial_corpus_files,
                'allocator_cleanup_during_replay':target in ALLOCATOR_TARGETS,
                'sanitizer_options':campaign_env.get('ASAN_OPTIONS'),
                **timer.evidence(),
                **statistics(text),
                'corpus_files':len(list(corpora[target].iterdir())),
                'artifacts':[p.name for p in sorted(artifacts.iterdir())], 'command':command}
            summary['targets'].append(report)
            # An exit without a completed campaign is not a successful fuzz run.
            failed |= code != 0 or not report['startup_budget_completed'] or not report['mutation_budget_completed'] or not report['mutation_executions'] or report['mutation_executions'] < 0 or not resource_limits_satisfied(report) or (FUZZ/'Cargo.lock').read_bytes() != locked
            print(json.dumps(report),flush=True)
    finally:
        summary['passed'] = len(summary['targets']) == len(targets) and not failed
        (output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    raise SystemExit(1 if failed else 0)


if __name__ == '__main__':
    main()
