#!/usr/bin/env python3
"""Fetch OCCT's public test dataset for the original-test bridge (REVIEW_NOTES.md R2, U1).

A separate, explicit step: the bridge never downloads anything. The archive
is the `opencascade-dataset-7.9.0.tar.xz` asset of the `V7_9_0_beta2`
release of Open-Cascade-SAS/OCCT, which upstream's own CI downloads for its
test runs. It is pinned by size and SHA-256, extracted into the ignored
`target/occt-test-data`, and never committed or redistributed: the archive
carries no licence file of its own, and its basis for use is being an asset
of an LGPL-2.1-with-exception release published for upstream's tests.

`--archive PATH` verifies an archive already on disk instead of
downloading. The script prints the inventory by directory and extension and
writes the file list to `target/occt-test-data/inventory.txt`.
"""
import argparse
from collections import Counter
import hashlib
from pathlib import Path
import shutil
import sys
import tarfile
import urllib.request

ROOT = Path(__file__).resolve().parents[2]
TARGET = ROOT / 'target/occt-test-data'
URL = ('https://github.com/Open-Cascade-SAS/OCCT/releases/download/V7_9_0_beta2/'
       'opencascade-dataset-7.9.0.tar.xz')
SIZE = 98_739_184
SHA256 = 'a92ed91c3271c299287c1c404bb9454d463251094bfc848acc44aee60e6a026c'
NAME = 'opencascade-dataset-7.9.0'


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def verified(path):
    if path.stat().st_size != SIZE:
        raise SystemExit(f'{path}: {path.stat().st_size} bytes, expected {SIZE}')
    if sha256(path) != SHA256:
        raise SystemExit(f'{path}: SHA-256 differs from the pinned {SHA256}')
    return path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', type=Path, help='verify this archive instead of downloading')
    parser.add_argument('--inventory-only', action='store_true',
                        help='only print the inventory of an existing extraction')
    args = parser.parse_args()
    extracted = TARGET / NAME
    if not args.inventory_only:
        TARGET.mkdir(parents=True, exist_ok=True)
        archive = TARGET / f'{NAME}.tar.xz'
        if args.archive:
            verified(args.archive)
            if args.archive.resolve() != archive.resolve():
                shutil.copyfile(args.archive, archive)
        elif not archive.is_file() or archive.stat().st_size != SIZE or sha256(archive) != SHA256:
            print(f'downloading {URL}', flush=True)
            partial = archive.with_suffix('.partial')
            with urllib.request.urlopen(URL) as response, partial.open('wb') as out:
                shutil.copyfileobj(response, out)
            verified(partial).rename(archive)
        verified(archive)
        if not extracted.is_dir():
            partial = TARGET / f'{NAME}.partial'
            shutil.rmtree(partial, ignore_errors=True)
            partial.mkdir()
            with tarfile.open(archive) as tar:
                for member in tar.getmembers():
                    # Regular files and directories inside the archive only.
                    name = Path(member.name)
                    if name.is_absolute() or '..' in name.parts or not (member.isfile() or member.isdir()):
                        raise SystemExit(f'refusing archive member {member.name}')
                tar.extractall(partial)
            inner = [p for p in partial.iterdir()]
            source = inner[0] if len(inner) == 1 and inner[0].is_dir() else partial
            source.rename(extracted)
            shutil.rmtree(partial, ignore_errors=True)
    if not extracted.is_dir():
        raise SystemExit(f'{extracted} is absent; run without --inventory-only first')
    files = sorted(p for p in extracted.rglob('*') if p.is_file())
    (TARGET / 'inventory.txt').write_text(
        ''.join(f'{p.relative_to(extracted)}\n' for p in files))
    size = sum(p.stat().st_size for p in files)
    print(f'{len(files)} files, {size / 1e6:.0f} MB in {extracted.relative_to(ROOT)}')
    print('by directory:', dict(Counter(p.relative_to(extracted).parts[0] for p in files).most_common()))
    print('by extension:', dict(Counter(p.suffix.lower() or '(none)' for p in files).most_common(8)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
