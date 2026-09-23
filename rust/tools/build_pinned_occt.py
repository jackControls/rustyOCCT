#!/usr/bin/env python3
"""Build a pristine, source-pinned headless OCCT SDK for optional native oracles.

Requires Git, tar, CMake and a C++17 compiler. Uses only the already
present Git commit; no package installation or network checkout is performed.
Choose a fresh output directory. The Rust kernel never links these libraries.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import shutil
import subprocess
import time

ROOT = Path(__file__).resolve().parents[2]
SOURCE = '3d097a0328e71b826377d4814ab05ec3c3d23871'
MODULES = ('FoundationClasses', 'ModelingData', 'ModelingAlgorithms',
           'Visualization', 'ApplicationFramework', 'DataExchange', 'Draw')
TOOLKITS = ('TKG3d', 'TKGeomAlgo', 'TKTopAlgo', 'TKPrim', 'TKBO')


def digest(path):
    result = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            result.update(block)
    return result.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT/'target/pinned-occt'/SOURCE[:12])
    parser.add_argument('--cmake', default='cmake')
    parser.add_argument('--jobs', type=int, default=2, choices=range(1,65))
    parser.add_argument('--toolkit', action='append', choices=TOOLKITS,
                        help='Repeat for additional headless toolkits; default TKG3d.')
    args = parser.parse_args()
    cmake = shutil.which(args.cmake)
    if not cmake:
        parser.error('CMake was not found; provide --cmake /path/to/cmake')
    output = args.output.resolve()
    if output.exists():
        parser.error('choose a fresh output directory; existing source/build data is preserved')
    # Check the source before creating any output or invoking a build.
    subprocess.run(['git', 'cat-file', '-e', SOURCE+'^{commit}'], cwd=ROOT, check=True)
    output.mkdir(parents=True)
    manifest_path = output/'build-manifest.json'
    record = {'source_reference': SOURCE, 'status': 'in_progress',
              'created_at_utc': datetime.now(timezone.utc).isoformat(),
              'platform': platform.platform(), 'toolkits': args.toolkit or ['TKG3d'],
              'jobs': args.jobs, 'exception_checks_enabled': True,
              'default_modules_disabled': list(MODULES), 'use_git_hash': False,
              'install_prefix': str(output/'install'), 'steps': [],
              'native_behavior_verified': False}

    def save():
        manifest_path.write_text(json.dumps(record, indent=2)+'\n')

    def run(name, command):
        print(name, flush=True)
        item = {'name': name, 'command': list(map(str, command)), 'log': name+'.log'}
        record['steps'].append(item)
        save()
        started = time.monotonic()
        with (output/item['log']).open('w') as log:
            process = subprocess.run(command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
        item.update(exit_code=process.returncode, elapsed_seconds=round(time.monotonic()-started, 3))
        save()
        process.check_returncode()

    save()
    try:
        run('cmake-version', [cmake, '--version'])
        archive = output/'source.tar'
        run('source-archive', ['git', 'archive', '--format=tar', '--output='+str(archive), SOURCE])
        record['source_archive_sha256'] = digest(archive)
        record['source_archive_bytes'] = archive.stat().st_size
        source, build, prefix = (output/name for name in ('source', 'build', 'install'))
        source.mkdir()
        run('extract', ['tar', '-xf', archive, '-C', source])
        command = [cmake, '-S', source, '-B', build, '-DCMAKE_BUILD_TYPE=Release',
                   '-DCMAKE_INSTALL_PREFIX='+str(prefix), '-DBUILD_LIBRARY_TYPE=Shared',
                   '-DBUILD_RELEASE_DISABLE_EXCEPTIONS=OFF', '-DUSE_GIT_HASH=OFF',
                   '-DCMAKE_EXPORT_NO_PACKAGE_REGISTRY=ON', '-DCMAKE_EXPORT_COMPILE_COMMANDS=ON',
                   '-DBUILD_GTEST=OFF', '-DBUILD_DOC_Overview=OFF', '-DBUILD_DOC_RefMan=OFF',
                   '-DUSE_TK=OFF', '-DUSE_FREETYPE=OFF', '-DUSE_TBB=OFF',
                   '-DUSE_OPENGL=OFF', '-DUSE_GLES2=OFF',
                   '-DBUILD_ADDITIONAL_TOOLKITS='+';'.join(record['toolkits'])]
        command.extend('-DBUILD_MODULE_'+name+'=OFF' for name in MODULES)
        run('configure', command)
        run('build', [cmake, '--build', build, '--parallel', str(args.jobs)])
        run('install', [cmake, '--install', build])
        libraries = sorted({p.resolve() for directory in ('lib', 'bin')
                            for p in (prefix/directory).glob('*')
                            if p.is_file() and ('.so' in p.name or p.suffix in ('.dylib', '.dll'))})
        if not libraries or not any('TKernel' in p.name for p in libraries):
            raise ValueError('the installed native libraries were not found')
        record['libraries'] = [{'path': str(p.relative_to(output)), 'sha256': digest(p),
                                'bytes': p.stat().st_size} for p in libraries]
        record['cmake_cache_sha256'] = digest(build/'CMakeCache.txt')
        record['compile_commands_sha256'] = digest(build/'compile_commands.json')
        record['version_header_sha256'] = digest(prefix/'include/opencascade/Standard_Version.hxx')
        record['status'] = 'built'
    except (OSError, ValueError, subprocess.SubprocessError) as error:
        record['status'] = 'failed'
        record['error'] = str(error)
        raise
    finally:
        record['finished_at_utc'] = datetime.now(timezone.utc).isoformat()
        save()
    print(f'Built pinned headless SDK: {prefix}\nManifest: {manifest_path}', flush=True)


if __name__ == '__main__':
    main()
