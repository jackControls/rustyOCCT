#!/usr/bin/env python3
"""Take Boolean captures' platform records (`platform-<name>/`) and the
platform's reviews from a CI run's observations.

The scheduled workflow "Rust Boolean captures on Linux" runs every Boolean
set's comparison with `--allow-unrecorded-platform`
(`run_boolean_comparisons.py`). A set whose capture has no record of the
runner's platform writes a `report.json` with `unrecorded_platform` (the
capture, its reviews file, the digests of the capture's and the run's
native rows, the cases differing beyond captured()'s allowance and the
largest deviation) and `platform_reviews_pending` (each native difference
without a review fingerprinted by the platform's rows), beside its
`native-observed.txt`. Download the run's artifacts
(`gh run download RUN -R OWNER/REPO -p 'boolean-linux-*' -D DIR`) and give
this tool the directory; without `--write` it only reports.

For each such report, against the capture as committed (a report whose
capture rows or observations changed since the run is refused):

* **The record.** The platform's rows are set beside the capture's with
  `compare_boolean.observation_differences`. When every case has the
  capture's status, validity, warnings, solid count and face, edge and
  vertex counts and every volume, area and centre lies within `--bound`
  (1e-11) of max(1, |value|) of the capture's, the record's review states
  that, with the largest deviation. Otherwise the record needs a review
  written by hand (`--review CAPTURE=TEXT`, after reading the differences
  this tool prints): a difference is never accepted silently.
* **The platform's reviews.** Reviews are fingerprinted by the native rows
  (`native_sha256`, the hash of the case's parsed row), so each reviewed
  difference needs its own review on another platform. A pending difference
  whose case lies within the bound of the capture's rows and which has the
  capture platform's review of the same differences (that review's
  fingerprint the capture's row) is given a copy of it under the platform's
  fingerprint, `platform_observed` naming the run and the deviation, as
  the existing Linux reviews of the split captures were. Any other pending
  difference (a new one, or one on a case differing beyond the bound) needs
  a review written by hand into the divergences file with its independent
  evidence; rerun the tool once it is there.

`--write` (with `--run` and `--commit`, naming the source) writes
`platform-<name>/native.txt` and `record.json` and appends the drafted
reviews for each capture that is complete (its record's review and every
pending difference's review in hand), and leaves the others untouched, so a
commit of the result never makes a comparison fail for want of a review.
Review the written files' diff before committing them. The compare scripts
then read the record (`compare_curve_surface.platform_record`: its digest
and review checked) and compare strictly on that platform.
"""
import argparse
import json
from pathlib import Path
import shutil
import sys

from build_pinned_occt import SOURCE, digest
from compare_boolean import observation_differences, parse_native
from compare_brep import review_for, sha, write
from compare_occt import ROOT

BOUND = 1e-11
PLATFORMS = {'linux': 'Linux x86_64', 'darwin': 'macOS arm64'}
EVIDENCE_KEYS = ('case', 'source_reference', 'oracle', 'native_sha256', 'differences')


def reports(directories):
    """The reports of unrecorded platforms under the directories."""
    out = []
    for directory in directories:
        for path in sorted(Path(directory).rglob('report.json')):
            report = json.loads(path.read_text())
            if 'unrecorded_platform' in report:
                out.append((path.parent, report))
    return out


def case_deviation(was, now):
    """The largest deviation of a case's values from the capture's, as
    observation_differences measures it (None: structurally different)."""
    found = observation_differences({'c': was}, {'c': now})
    return None if found['structural'] else found['deviation']


def plan(directory, report, root=ROOT, bound=BOUND, reviews_given=None):
    """What one report gives: {'capture', 'platform', 'capture_platform',
    'directory', 'cases', 'record' (the record's review text, or None),
    'differences' (the cases beyond the capture's allowance, [(case,
    what)]), 'deviation', 'drafted' ([(review, the case's deviation)]),
    'reviewed' (cases whose platform review is already in the file),
    'by_hand' (evidence needing a review written by hand)}."""
    pending = report['unrecorded_platform']
    if report['source_reference'] != SOURCE:
        raise ValueError(f'{directory}: another source reference')
    capture = root/pending['capture']
    platform = pending['platform']
    observed = directory/'native-observed.txt'
    if digest(observed) != pending['observations_sha256']:
        raise ValueError(f'{directory}: native-observed.txt is not the report\'s')
    if digest(capture/'native.txt') != pending['capture_observations_sha256']:
        raise ValueError(f'{capture.name}: the capture changed since the run; run the comparison again')
    capture_platform = json.loads((capture/'capture.json').read_text())['platform']
    if capture_platform != pending['capture_platform'] or capture_platform == platform:
        raise ValueError(f'{capture.name}: not another platform\'s observation')
    was = parse_native((capture/'native.txt').read_text())
    now = parse_native(observed.read_text())
    found = observation_differences(was, now)
    if found['missing'] or found['extra']:
        raise ValueError(f'{capture.name}: the run\'s cases differ from the capture\'s')
    name = capture.name
    out = {'capture': capture, 'platform': platform, 'capture_platform': capture_platform,
           'directory': directory, 'differences': found['structural'], 'deviation': found['deviation'],
           'record': None, 'drafted': [], 'reviewed': [], 'by_hand': [], 'cases': len(now)}
    given = (reviews_given or {}).get(name)
    if given:
        out['record'] = given
    elif not found['structural'] and found['deviation'] <= bound:
        values = ('equal to the capture\'s' if found['deviation'] == 0 else
                  f'within {found["deviation"]:.2g} of max(1, |value|) of the capture\'s')
        out['record'] = (
            f'The same cases, statuses, validity, warnings, solid counts and face, edge and vertex counts '
            f'(before and after unifying) as the {PLATFORMS.get(capture_platform, capture_platform)} capture '
            f'for every Boolean; volumes, areas and centres {values} (take_platform_records.py, '
            f'bound {bound:g}).')
    reviews = json.loads((root/pending['reviews']).read_text())['reviews'] \
        if (root/pending['reviews']).exists() else []
    for item in report.get('platform_reviews_pending', []):
        evidence = {k: item[k] for k in EVIDENCE_KEYS}
        case = evidence['case']
        if sha(json.dumps(now[case])) != evidence['native_sha256']:
            raise ValueError(f'{name}: {case}\'s fingerprint is not its observed row\'s')
        if review_for(evidence, reviews):
            out['reviewed'].append(case)
            continue
        deviation = case_deviation(was[case], now[case])
        theirs = sha(json.dumps(was[case]))
        source = next((r for r in reviews if r.get('reason') and r.get('independent_evidence')
                       and r['native_sha256'] == theirs
                       and all(r.get(k) == evidence[k] for k in EVIDENCE_KEYS if k != 'native_sha256')), None)
        if source is None or deviation is None or deviation > bound:
            out['by_hand'].append(dict(evidence, why=(
                'no review of the capture\'s row with these differences' if source is None
                else 'the case differs from the capture\'s beyond the bound')))
            continue
        out['drafted'].append((dict(evidence, reason=source['reason'],
                                    independent_evidence=source['independent_evidence']), deviation))
    return out


def take(directories, root=ROOT, bound=BOUND, reviews_given=None, run=None, commit=None,
         runner='ubuntu-24.04', write_files=False):
    """Plan every report and, with `write_files`, write the complete ones.
    Returns the plans."""
    plans = [plan(d, r, root, bound, reviews_given) for d, r in reports(directories)]
    names = [p['capture'].name for p in plans]
    if len(names) != len(set(names)):
        raise ValueError('a capture observed twice: give one run\'s artifacts')
    for p in plans:
        p['complete'] = p['record'] is not None and not p['by_hand']
        p['written'] = False
        if (p['capture']/('platform-'+p['platform'])/'record.json').exists():
            p['complete'] = False
            p['note'] = 'already recorded'
    if not write_files:
        return plans
    if not run or not commit:
        raise ValueError('--write needs --run and --commit to name the source')
    for p in plans:
        if not p['complete']:
            continue
        label = PLATFORMS.get(p['platform'], p['platform'])
        record = p['capture']/('platform-'+p['platform'])
        record.mkdir(exist_ok=True)
        shutil.copyfile(p['directory']/'native-observed.txt', record/'native.txt')
        write(record/'record.json', {
            'platform': p['platform'],
            'source': f'GitHub Actions run {run} (commit {commit[:8]}), {runner}, source-pinned OCCT 8.1.0',
            'source_reference': SOURCE,
            'observations_sha256': digest(record/'native.txt'),
            'review': p['record']})
        if p['drafted']:
            reviews_path = root/json.loads((p['directory']/'report.json').read_text())['unrecorded_platform']['reviews']
            data = json.loads(reviews_path.read_text())
            capture_label = PLATFORMS.get(p['capture_platform'], p['capture_platform']).split()[0]
            for review, deviation in p['drafted']:
                data['reviews'].append(dict(review, platform_observed=(
                    f'{label} (GitHub Actions {runner}, run {run} at {commit[:8]}); source-pinned OCCT 8.1.0; '
                    f'the same difference as the {capture_label} review, the native values within '
                    f'{deviation:.2g} of max(1, |value|) of the {capture_label} capture\'s')))
            write(reviews_path, data)
        p['written'] = True
    return plans


def summary(plans):
    """A Markdown table of the plans."""
    lines = ['| Capture | Cases | Beyond rounding | Largest deviation | Reviewed | Reviews drafted | By hand | State |',
             '|---|---|---|---|---|---|---|---|']
    for p in sorted(plans, key=lambda p: p['capture'].name):
        state = p.get('note') or ('written' if p['written'] else 'complete' if p['complete']
                                  else 'needs a record review' if p['record'] is None else 'needs reviews by hand')
        lines.append(f'| `{p["capture"].name}` | {p["cases"]} | {len(p["differences"])} | {p["deviation"]:.2g} '
                     f'| {len(p["reviewed"])} | {len(p["drafted"])} | {len(p["by_hand"])} | {state} |')
    details = []
    for p in plans:
        if p['differences']:
            details.append(f'* `{p["capture"].name}` beyond rounding: '
                           + ', '.join(f'{c} ({w})' for c, w in p['differences']))
        for e in p['by_hand']:
            details.append(f'* `{p["capture"].name}` {e["case"]} {" ".join(e["differences"])}: {e["why"]} '
                           f'(native_sha256 {e["native_sha256"]})')
    return '\n'.join(lines+([''] if details else [])+details)+'\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('artifacts', type=Path, nargs='+', help="a run's downloaded artifacts")
    parser.add_argument('--write', action='store_true')
    parser.add_argument('--run', help='the GitHub Actions run id')
    parser.add_argument('--commit', help="the run's commit")
    parser.add_argument('--runner', default='ubuntu-24.04')
    parser.add_argument('--bound', type=float, default=BOUND)
    parser.add_argument('--review', action='append', default=[], metavar='CAPTURE=TEXT',
                        help="a record's review written by hand, for a capture differing beyond the bound")
    parser.add_argument('--summary', type=Path, help='also append the Markdown summary to this file')
    args = parser.parse_args()
    given = {}
    for item in args.review:
        name, _, text = item.partition('=')
        if not text:
            parser.error('--review takes CAPTURE=TEXT')
        given[name] = text
    plans = take(args.artifacts, ROOT, args.bound, given, args.run, args.commit, args.runner, args.write)
    if not plans:
        print('no observation of an unrecorded platform in the artifacts')
    text = summary(plans) if plans else ''
    print(text, end='')
    if args.summary and text:
        with args.summary.open('a') as stream:
            stream.write(text)


if __name__ == '__main__':
    sys.exit(main())
