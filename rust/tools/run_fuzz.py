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
TARGETS = ['predicates','intersections','modeling','curved','splines','surfaces','roots','spline_intersections']


def seed_corpus(target):
    corpus = FUZZ/'corpus'/target
    corpus.mkdir(parents=True, exist_ok=True)

    def save(data):
        path = corpus/hashlib.sha256(data).hexdigest()
        if not path.exists():
            path.write_bytes(data)

    if target == 'roots':
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


def run_process(command, log, timeout, env, tick=None):
    # Kill the entire build/fuzzer process group on a wall-clock timeout.
    process = subprocess.Popen(command,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,
                               env=env,start_new_session=True)
    deadline=time.monotonic()+timeout
    try:
        while True:
            code=process.poll()
            if code is not None: return code
            if tick is not None: tick()
            remaining=deadline-time.monotonic()
            if remaining<=0: raise subprocess.TimeoutExpired(command,timeout)
            try:
                return process.wait(timeout=min(.1,remaining) if tick is not None else remaining)
            except subprocess.TimeoutExpired:
                if time.monotonic()>=deadline: raise
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
    def __init__(self, log_path, stop_file, seconds):
        self.log_path,self.stop_file,self.seconds=log_path,stop_file,seconds
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
                'mutation_seconds_before_stop':round(self.stopped-self.initialized,2) if self.stopped is not None else None,
                'mutation_budget_completed':self.stopped is not None}


def statistics(text):
    executed = re.findall(r'stat::number_of_executed_units:\s*(\d+)',text)
    initialized = re.findall(r'#(\d+)\s+INITED\b',text)
    coverage = re.findall(r'cov: (\d+)',text)
    total = int(executed[-1]) if executed else None
    replayed = int(initialized[-1]) if initialized else None
    return {'executions':total,'initial_executions':replayed,
            'mutation_executions':total-replayed if total is not None and replayed is not None else None,
            'coverage_edges':int(coverage[-1]) if coverage else None}


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
            # Its outer subprocess budget also bounds builds/corpus startup.
            started = time.monotonic()
            print(f'Fuzzing {target} for {args.seconds}s; log: {log_path}',flush=True)
            with tempfile.TemporaryDirectory(prefix=f'.{target}-control-',dir=output) as control:
                stop_file=Path(control)/'stop'
                timer=MutationBudget(log_path,stop_file,args.seconds)
                command = ['cargo',f'+{args.toolchain}','fuzz','run',target,str(corpora[target]),
                    '--fuzz-dir',str(FUZZ),'--sanitizer','address','--',
                    '-max_total_time=0',f'-stop_file={stop_file}','-timeout=20','-rss_limit_mb=2048',
                    '-max_len=256',f'-seed={args.seed}',f'-artifact_prefix={artifacts}/','-print_final_stats=1']
                with log_path.open('w') as log:
                    code = run_process(command,log,args.seconds+600,env,timer.tick)
            text = log_path.read_text(errors='replace')
            report = {'target':target,'exit_code':code,'elapsed_seconds':round(time.monotonic()-started,2),
                **timer.evidence(),
                **statistics(text),
                'corpus_files':len(list(corpora[target].iterdir())),
                'artifacts':[p.name for p in sorted(artifacts.iterdir())], 'command':command}
            summary['targets'].append(report)
            # An exit without a completed campaign is not a successful fuzz run.
            failed |= code != 0 or not report['mutation_budget_completed'] or not report['mutation_executions'] or report['mutation_executions'] < 0 or (FUZZ/'Cargo.lock').read_bytes() != locked
            print(json.dumps(report),flush=True)
    finally:
        summary['passed'] = len(summary['targets']) == len(targets) and not failed
        (output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    raise SystemExit(1 if failed else 0)


if __name__ == '__main__':
    main()
