#!/usr/bin/env python3
"""Compare complete degree-elevation grids with independent equations and OCCT.

Native exceptions and malformed outputs are recorded as discrepancies. Only
exactly fingerprinted, explicitly reviewed observations can pass the gate;
they are always reported separately from matches. Rust must agree exactly
with freshly reconstructed coefficient equations for every input.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys
from compare_occt import ROOT, run
from build_pinned_occt import SOURCE, digest
import knot_editing_reference as curves
import surface_knot_reference as surfaces
from generate_degree_elevation_fixtures import generate_curves, generate_surfaces

ABSOLUTE = F(1, 10**10)
RELATIVE = F(2, 10**12)


def sha(text):
    return hashlib.sha256(text.encode()).hexdigest()


def lines_by_name(text):
    result = {}
    for line in text.splitlines():
        words = line.split()
        if not words or words[0] in result:
            raise ValueError('empty or duplicate observation')
        result[words[0]] = line
    if not result:
        raise ValueError('empty observation set')
    return result


def decode(kind, line, native=False):
    """Normalize both protocols, including strict validation of curve axes."""
    if kind == 'surface':
        value = next(iter(surfaces.decode(line, native).values()))
        _, (axes, _), bounds = value
        if tuple(surfaces.domain((p, periodic, tuple(map(F, knots)), mults))
                 for p, periodic, knots, mults in axes) != bounds:
            raise ValueError('wrong surface parameter domains')
        return value
    flags, curve, bounds = next(iter(curves.decode(line, native).values()))
    p, periodic, knots, mults, controls = curve
    if (len(knots)<2 or len(knots)!=len(mults)
            or any(a>=b for a,b in zip(knots,knots[1:]))
            or any(m<1 or m>p+int(not periodic and i in (0,len(knots)-1)) for i,m in enumerate(mults))
            or (periodic and mults[0]!=mults[-1])
            or surfaces.count(curve[:4])!=len(controls)):
        raise ValueError('invalid curve axis or control count')
    exact_axis = p, periodic, tuple(map(F,knots)), mults
    if tuple(map(F,bounds)) != surfaces.domain(exact_axis):
        raise ValueError('wrong curve parameter domain')
    return flags, ((curve[:4],), controls), (bounds,)


def compare(kind, data, records, rust, exact, oracle, native_source=None, reviews=()):
    inputs = lines_by_name(data)
    actual, expected = lines_by_name(rust), lines_by_name(exact)
    observed = {r['case']: r for r in records}
    if (len(observed)!=len(records) or inputs.keys()!=observed.keys()
            or inputs.keys()!=actual.keys() or inputs.keys()!=expected.keys()):
        raise ValueError('case sets differ')
    report = {'kind':kind, 'cases':len(inputs), 'oracle':oracle,
              'source_reference':SOURCE, 'native_source_reference':native_source,
              'independently_verified':0, 'matched_cases':0,
              'reviewed_differences':[], 'failures':[],
              'largest_error_fraction_of_budget':0.,
              'comparison_budget':{'homogeneous_absolute':float(ABSOLUTE),
                                   'homogeneous_relative':float(RELATIVE)}}
    for name, source in inputs.items():
        want = decode(kind, expected[name])
        if decode(kind, actual[name]) != want:
            report['failures'].append({'case':name, 'reason':'Rust differs from independent coefficient equations'})
            continue
        report['independently_verified'] += 1
        record = observed[name]
        native = record['stdout']
        differences = []
        classification = None
        if record['exit_code'] != 0:
            # Crashes and killed processes never receive behavior exemptions.
            report['failures'].append({'case':name, 'reason':'native process failed',
                                       'exit_code':record['exit_code']})
            continue
        if native == name+' E\n':
            classification = 'native_exception'
            differences.append('native degree elevation rejected a validated input')
        else:
            try:
                names = lines_by_name(native)
                if list(names) != [name]:
                    raise ValueError('wrong native case identity')
                nf, (na, nc), nb = decode(kind, native, native=True)
                flags, (axes, controls), bounds = want
                if flags != nf: differences.append('operation success flags')
                if axes != na: differences.append('axes, knots or multiplicities')
                if bounds != nb: differences.append('parameter domains')
                if len(controls) != len(nc): differences.append('control count')
                if not differences:
                    fractions = []
                    for i, (a, b) in enumerate(zip(nc, controls)):
                        x,y,z,w = map(F,a)
                        for j, (value, required) in enumerate(zip((x*w,y*w,z*w,w), b)):
                            fraction = abs(value-required)/(ABSOLUTE+RELATIVE*abs(required))
                            fractions.append(float(fraction))
                            if fraction > 1:
                                differences.append(f'homogeneous control {i} component {j}')
                    report['largest_error_fraction_of_budget'] = max(
                        report['largest_error_fraction_of_budget'], *fractions)
                classification = 'different_result'
            except (ValueError, StopIteration, IndexError, ZeroDivisionError) as error:
                classification = 'invalid_native_output'
                differences.append(str(error))
        if not differences:
            report['matched_cases'] += 1
            continue
        evidence = {'kind':kind, 'case':name, 'oracle':oracle,
                    'native_source_reference':native_source,
                    'input_sha256':sha(source), 'exact_sha256':sha(expected[name]),
                    'native_stdout_sha256':sha(native),
                    'native_stderr_sha256':sha(record['stderr']),
                    'classification':classification, 'differences':differences}
        review = next((r for r in reviews if r.get('reason') and r.get('independent_evidence')
                       and all(r.get(k)==v for k,v in evidence.items())), None)
        if review:
            report['reviewed_differences'].append(review)
        else:
            report['failures'].append(evidence)
    return report


def verify_sdk(prefix, manifest_path):
    if manifest_path is None:
        return None
    record = json.loads(manifest_path.read_text())
    capture = json.loads((ROOT/'rust/fixtures/occt-degree-elevation-capture.json').read_text())
    if (record.get('status')!='built' or record.get('source_reference')!=SOURCE
            or not record.get('exception_checks_enabled') or record.get('use_git_hash') is not False
            or Path(record['install_prefix']).resolve()!=prefix
            or record.get('source_archive_sha256')!=capture['source_archive_sha256']):
        raise ValueError('source-pinned SDK manifest does not match this oracle')
    base = manifest_path.resolve().parent
    if digest(base/'source.tar') != record['source_archive_sha256']:
        raise ValueError('source archive changed')
    if digest(prefix/'include/opencascade/Standard_Version.hxx') != record['version_header_sha256']:
        raise ValueError('SDK version header changed')
    if not record['libraries']:
        raise ValueError('empty SDK library manifest')
    for library in record['libraries']:
        if digest(base/library['path']) != library['sha256']:
            raise ValueError('SDK library changed: '+library['path'])
    return SOURCE


def verify_loaded_libraries(text, lib, platform):
    """Resolve every OCCT dependency; a version string alone is insufficient."""
    if platform == 'darwin':
        paths = re.findall(r'^dyld\[\d+\]: <[^>]+> (/.*/libTK[^/]+\.dylib)$', text, re.M)
    else:
        paths = re.findall(r'^\s*libTK\S+\s+=>\s+(/.+?)\s+\(0x[0-9a-fA-F]+\)', text, re.M)
        if re.search(r'^\s*libTK\S+\s+=>\s+not found', text, re.M):
            raise ValueError('OCCT runtime library missing')
    resolved = {Path(p).resolve() for p in paths}
    for toolkit in ('TKernel', 'TKMath', 'TKG2d', 'TKG3d'):
        if not any(p.name.startswith('lib'+toolkit+'.') for p in resolved):
            raise ValueError('missing loaded OCCT library: '+toolkit)
    if any(p.parent != lib.resolve() for p in resolved):
        raise ValueError('OCCT runtime library is outside the selected SDK')
    return [{'path':str(p), 'sha256':digest(p)} for p in sorted(resolved)]


def verify_capture(metadata, records, kind, data, native_source, sdk_digest):
    if (metadata['input_sha256'] != sha(data)
            or metadata['native_source_reference'] != native_source
            or metadata['sdk_manifest_sha256'] != sdk_digest
            or metadata['source_sha256'] != digest(ROOT/f'rust/tools/occt_{kind}_degree_oracle.cpp')
            or metadata['observations_sha256'] != sha(json.dumps(records, sort_keys=True))):
        raise ValueError('capture provenance or corpus changed')
    if not metadata.get('loaded_libraries'):
        raise ValueError('capture lacks runtime library evidence')
    for library in metadata['loaded_libraries']:
        if digest(Path(library['path'])) != library['sha256']:
            raise ValueError('captured runtime library changed')


def capture(kind, data, prefix, output):
    include = next((p for p in [prefix/'include/opencascade',prefix/'inc',prefix/'include']
                    if (p/'Geom_BSplineCurve.hxx').exists()), None)
    lib = next((p for p in [prefix/'lib',prefix/'lib64',prefix/'lib/x86_64-linux-gnu',prefix/'lib/aarch64-linux-gnu']
                if any(p.glob('libTKG3d.*'))), None)
    if not include or not lib:
        raise ValueError('OCCT geometry SDK not found')
    executable = output/('occt_'+kind+'_degree_oracle')
    command = shlex.split(os.environ.get('CXX','c++'))+[
        '-std=c++17','-O2',str(ROOT/f'rust/tools/occt_{kind}_degree_oracle.cpp'),
        '-I'+str(include),'-L'+str(lib),'-Wl,-rpath,'+str(lib),'-o',str(executable),
        '-lTKG3d','-lTKMath','-lTKernel']
    built = subprocess.run(command, cwd=ROOT, text=True, capture_output=True)
    (output/'native-build.log').write_text(built.stdout+built.stderr)
    built.check_returncode()
    env = os.environ.copy()
    for key in ('LD_LIBRARY_PATH','DYLD_LIBRARY_PATH'):
        env[key] = str(lib)+(os.pathsep+env[key] if env.get(key) else '')
    records = []
    versions = set()
    for line in data.splitlines():
        result = subprocess.run([str(executable)], input=line+'\n', cwd=ROOT,
                                text=True, capture_output=True, timeout=20, env=env)
        stderr_lines = result.stderr.splitlines()
        if not stderr_lines or not stderr_lines[0].startswith('OCCT '):
            raise ValueError('missing native version')
        versions.add(stderr_lines[0])
        records.append({'case':line.split()[0], 'exit_code':result.returncode,
                        'stdout':result.stdout, 'stderr':result.stderr})
        # Persist each completed observation, even if a later native input fails.
        (output/'native.json').write_text(json.dumps(records,indent=2)+'\n')
    if len(versions)!=1:
        raise ValueError('native versions differ')
    link = ['otool','-L',str(executable)] if sys.platform=='darwin' else ['ldd',str(executable)]
    linked = subprocess.run(link, env=env, text=True, capture_output=True, check=True)
    (output/'linked-libraries.txt').write_text(linked.stdout+linked.stderr)
    if sys.platform == 'darwin':
        traced_env = dict(env, DYLD_PRINT_LIBRARIES='1')
        traced = subprocess.run([str(executable)], input=data.splitlines()[0]+'\n',
                                env=traced_env, text=True, capture_output=True, timeout=20, check=True)
        loaded_text = traced.stderr
    else:
        loaded_text = linked.stdout+linked.stderr
    (output/'loaded-libraries.txt').write_text(loaded_text)
    loaded = verify_loaded_libraries(loaded_text, lib, sys.platform)
    metadata = {'oracle':versions.pop(), 'input_sha256':sha(data),
                'probe_sha256':digest(executable), 'build_command':command,
                'observations_sha256':sha(json.dumps(records, sort_keys=True)),
                'loaded_libraries':loaded,
                'source_sha256':digest(ROOT/f'rust/tools/occt_{kind}_degree_oracle.cpp')}
    return records, metadata


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--occt-root', type=Path, required=True)
    parser.add_argument('--sdk-manifest', type=Path)
    parser.add_argument('--output', type=Path, default=ROOT/'target/degree-elevation-oracle')
    parser.add_argument('--reuse-capture', action='store_true')
    parser.add_argument('--strict-native', action='store_true')
    args = parser.parse_args()
    prefix = args.occt_root.resolve()
    native_source = verify_sdk(prefix, args.sdk_manifest)
    sdk_digest = digest(args.sdk_manifest) if args.sdk_manifest else None
    review_path = ROOT/'rust/fixtures/occt-degree-elevation-divergences.json'
    reviews = json.loads(review_path.read_text())['reviews'] if review_path.exists() and not args.strict_native else []
    reports = []
    for kind, plural, generate, example in [
        ('curve','curves',generate_curves,'knot_editing_oracle'),
        ('surface','surfaces',generate_surfaces,'surface_knot_oracle')]:
        output = args.output.resolve()/plural
        output.mkdir(parents=True,exist_ok=True)
        data = (ROOT/f'rust/fixtures/degree-elevation-{plural}.txt').read_text()
        generated = generate()
        if generated != (ROOT/f'rust/fixtures/degree-elevation-{plural}.tsv').read_text():
            raise ValueError('fixtures differ from freshly reconstructed coefficient equations')
        exact = '\n'.join(row.split(' | ',1)[1] for row in generated.splitlines())+'\n'
        if args.reuse_capture:
            records = json.loads((output/'native.json').read_text())
            metadata = json.loads((output/'capture.json').read_text())
            verify_capture(metadata, records, kind, data, native_source, sdk_digest)
        else:
            records, metadata = capture(kind,data,prefix,output)
            metadata['native_source_reference'] = native_source
            metadata['sdk_manifest_sha256'] = sdk_digest
            (output/'capture.json').write_text(json.dumps(metadata,indent=2)+'\n')
        rust = run(['cargo','+stable','run','--quiet','--locked','--release','--example',example],
                   input=data, cwd=ROOT).stdout
        for name,text in [('inputs.txt',data),('exact.tsv',exact),('rust.tsv',rust)]:
            (output/name).write_text(text)
        report = compare(kind,data,records,rust,exact,metadata['oracle'],native_source,reviews)
        (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
        reports.append(report)
        print(json.dumps({k:len(v) if k in ('failures','reviewed_differences') else v
                          for k,v in report.items()},indent=2),flush=True)
    (args.output/'report.json').write_text(json.dumps(reports,indent=2)+'\n')
    if any(r['failures'] for r in reports):
        raise SystemExit(1)


if __name__ == '__main__':
    main()
