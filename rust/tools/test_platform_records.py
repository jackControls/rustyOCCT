"""The Boolean captures' platform records on CI: the shards of
run_boolean_comparisons.py and the scheduled workflow, the unrecorded
platform's comparison in compare_boolean.py, and take_platform_records.py
turning a run's observations into records and reviews."""
import json
from pathlib import Path
import re
import tempfile
import unittest
from unittest.mock import patch

from build_pinned_occt import SOURCE, digest
import compare_boolean
from compare_brep import review_for, sha, write
import compare_curve_surface
import run_boolean_comparisons as runner
import take_platform_records as records

ORACLE = 'OCCT 8.1.0 BRepAlgoAPI_Fuse, BRepAlgoAPI_Cut, BRepAlgoAPI_Common'
WORKFLOW = runner.ROOT/'.github/workflows/rust-boolean-linux.yml'


def native_text(volume=875.0, counts='14 32 20 10 24 16', status='done', cut=350.0):
    """Two cases' native rows (compare_boolean.parse_native's format)."""
    return (f'rects_fuse {status} 1 1 0\nS {volume!r} 650.0 7.5 7.5 2.5 {counts} 1\n'
            f'rects_cut done 1 1 0\nS {cut!r} 340.0 4.9 2.8 0.01 8 18 12 8 18 12 1\n')


# rects_cut's volume a few ulps off on Linux: its row's fingerprint is no longer
# the macOS review's.
ULP = 350.0*(1+2**-50)


class Fixture:
    """A capture (macOS rows, one reviewed difference on `rects_cut`) and a
    Linux run's artifact directory under a temporary root."""

    def __init__(self, root, observed, extra_pending=()):
        self.root = root
        self.capture = root/'rust/fixtures/occt-boolean-demo-preimplementation'
        self.capture.mkdir(parents=True)
        (self.capture/'native.txt').write_text(native_text())
        (self.capture/'inputs.txt').write_text('inputs\n')
        (self.capture/'oracle.cpp').write_text('// probe\n')
        write(self.capture/'capture.json', {
            'source_reference': SOURCE, 'platform': 'darwin', 'rust_demo_exists': False,
            'input_sha256': digest(self.capture/'inputs.txt'),
            'probe_source_sha256': digest(self.capture/'oracle.cpp'),
            'observations_sha256': digest(self.capture/'native.txt')})
        was = compare_boolean.parse_native(native_text())
        self.reviews = root/'rust/fixtures/occt-boolean-demo-divergences.json'
        write(self.reviews, {'note': 'demo', 'reviews': [{
            'case': 'rects_cut', 'source_reference': SOURCE, 'oracle': ORACLE,
            'native_sha256': sha(json.dumps(was['rects_cut'])), 'differences': ['rust_entity_counts'],
            'reason': 'a reason', 'independent_evidence': 'some evidence',
            'platform_observed': 'macOS arm64; source-pinned OCCT 8.1.0'}]})
        self.artifact = root/'artifact/boolean-linux-observations-0/boolean-demo-oracle'
        self.artifact.mkdir(parents=True)
        (self.artifact/'native-observed.txt').write_text(observed)
        now = compare_boolean.parse_native(observed)
        pending = [{'case': c, 'source_reference': SOURCE, 'oracle': ORACLE,
                    'native_sha256': sha(json.dumps(now[c])), 'differences': d, 'solids': now[c][3]}
                   for c, d in [('rects_cut', ['rust_entity_counts']), *extra_pending]]
        found = compare_boolean.observation_differences(was, now)
        write(self.artifact/'report.json', {
            'source_reference': SOURCE, 'oracle': ORACLE, 'failures': [],
            'unrecorded_platform': {
                'capture_platform': 'darwin', 'platform': 'linux',
                'structural_differences': [{'case': c, 'what': w} for c, w in found['structural']],
                'largest_deviation': found['deviation'],
                'capture': self.capture.relative_to(root).as_posix(),
                'reviews': self.reviews.relative_to(root).as_posix(),
                'capture_observations_sha256': digest(self.capture/'native.txt'),
                'observations_sha256': digest(self.artifact/'native-observed.txt')},
            'platform_reviews_pending': pending})

    def take(self, **kwargs):
        return records.take([self.root/'artifact'], root=self.root, **kwargs)


class ShardTests(unittest.TestCase):
    def test_every_boolean_set_runs_in_exactly_one_shard_within_budget(self):
        names = [name for name, _ in runner.comparisons()]
        scripts = {p.name for p in runner.TOOLS.glob('compare_*_boolean.py')}
        self.assertEqual(set(names), scripts | {'compare_boolean.py', 'compare_boolean.py --splines',
                                                'compare_polyhedral.py'})
        self.assertEqual(len(names), len(set(names)))
        # Every set has its measured time (add a new set's with its script).
        self.assertEqual(set(runner.SECONDS), set(names))
        shards = runner.shards()
        listed = [name for shard in shards for name, _ in shard]
        self.assertEqual(sorted(listed), sorted(names))
        for shard in shards:
            self.assertLessEqual(runner.estimate(shard), runner.BUDGET)
        self.assertEqual(shards, runner.shards())

    def test_workflow_runs_every_shard_on_schedule_and_dispatch_only(self):
        text = WORKFLOW.read_text()
        triggers = text[text.index('\non:'):text.index('\npermissions:')]
        self.assertIn('schedule:', triggers)
        self.assertIn('workflow_dispatch:', triggers)
        self.assertNotIn('push:', triggers)
        self.assertNotIn('pull_request:', triggers)
        matrix = re.findall(r'\n\s+shard: \[([^\]]*)\]', text)
        self.assertEqual(len(matrix), 1)
        self.assertEqual([int(k) for k in matrix[0].split(',')], list(range(runner.SHARDS)))
        self.assertIn(f'--shards {runner.SHARDS}\n', text)
        self.assertIn('rust/tools/run_boolean_comparisons.py', text)
        self.assertIn('rust/tools/take_platform_records.py', text)
        # The per-push B-rep job stays strict.
        kernel = (runner.ROOT/'.github/workflows/rust-kernel.yml').read_text()
        self.assertNotIn('--allow-unrecorded-platform', kernel)


class UnrecordedPlatformTests(unittest.TestCase):
    def test_observation_differences(self):
        was = compare_boolean.parse_native(native_text())
        found = compare_boolean.observation_differences(was, compare_boolean.parse_native(native_text()))
        self.assertEqual((found['structural'], found['deviation']), ([], 0.0))
        found = compare_boolean.observation_differences(
            was, compare_boolean.parse_native(native_text(volume=875.0*(1+2**-50), cut=ULP)))
        self.assertEqual(found['structural'], [])
        self.assertAlmostEqual(found['deviation'], 2**-50, delta=2**-52)
        for text, what in [(native_text(volume=875.0*(1+1e-8)), 'values'),
                           (native_text(counts='14 32 20 10 24 15'), 'counts'),
                           (native_text(status='failed'), 'status')]:
            found = compare_boolean.observation_differences(was, compare_boolean.parse_native(text))
            self.assertEqual(found['structural'], [('rects_fuse', what)])
        found = compare_boolean.observation_differences(was, {'rects_fuse': was['rects_fuse']})
        self.assertEqual(found['missing'], ['rects_cut'])

    def test_unrecorded_platform_is_reported_and_a_record_makes_it_strict(self):
        with tempfile.TemporaryDirectory() as directory:
            fixture = Fixture(Path(directory), native_text(volume=875.0*(1+2**-50), cut=ULP))
            demo = compare_boolean.Set(False)
            demo.capture, demo.exists_key = fixture.capture, 'rust_demo_exists'
            observed = compare_boolean.parse_native(native_text(volume=875.0*(1+2**-50), cut=ULP))
            with patch.object(compare_boolean, 'SET', demo), \
                    patch.object(compare_boolean, 'native_input', return_value='inputs\n'), \
                    patch('sys.platform', 'linux'):
                metadata = compare_boolean.capture_metadata()
                self.assertFalse(compare_boolean.platform_recorded(metadata))
                with self.assertRaisesRegex(ValueError, 'no reviewed linux record'):
                    compare_boolean.captured(observed)
                pending = compare_boolean.unrecorded(observed)
                self.assertEqual(pending['structural_differences'], [])
                self.assertEqual(pending['capture_platform'], 'darwin')
                # Changed evidence still fails, record or not.
                (fixture.capture/'inputs.txt').write_text('other\n')
                with self.assertRaisesRegex(ValueError, 'evidence changed'):
                    compare_boolean.unrecorded(observed)
                (fixture.capture/'inputs.txt').write_text('inputs\n')
                plans = fixture.take(run='123', commit='0123456789abcdef', write_files=True)
                self.assertTrue(plans[0]['written'])
                self.assertTrue(compare_boolean.platform_recorded(metadata))
                compare_boolean.captured(observed)
                beyond = compare_boolean.parse_native(native_text(volume=876.0))
                with self.assertRaisesRegex(ValueError, 'rects_fuse differs from the capture'):
                    compare_boolean.captured(beyond)


class TakeRecordsTests(unittest.TestCase):
    def test_rounding_gives_a_record_and_the_reviews_under_linux_fingerprints(self):
        with tempfile.TemporaryDirectory() as directory:
            observed = native_text(volume=875.0*(1+2**-50), cut=ULP)
            fixture = Fixture(Path(directory), observed)
            plans = fixture.take()
            self.assertEqual(len(plans), 1)
            self.assertTrue(plans[0]['complete'])
            self.assertFalse(plans[0]['written'])
            self.assertFalse((fixture.capture/'platform-linux').exists())
            self.assertIn('| written |', records.summary(fixture.take(run='123', commit='0123456789abcdef',
                                                                      write_files=True)))
            record = json.loads((fixture.capture/'platform-linux/record.json').read_text())
            self.assertEqual((fixture.capture/'platform-linux/native.txt').read_text(), observed)
            self.assertEqual(record['observations_sha256'], digest(fixture.capture/'platform-linux/native.txt'))
            self.assertIn('run 123 (commit 01234567)', record['source'])
            self.assertIn('as the macOS arm64 capture', record['review'])
            # platform_record reads the record on Linux.
            metadata = json.loads((fixture.capture/'capture.json').read_text())
            with patch('sys.platform', 'linux'):
                self.assertEqual(compare_curve_surface.platform_record(fixture.capture, metadata), observed)
            # The Linux review: the macOS review's reason and evidence under
            # the Linux row's fingerprint, found by review_for.
            reviews = json.loads(fixture.reviews.read_text())['reviews']
            self.assertEqual(len(reviews), 2)
            now = compare_boolean.parse_native(observed)
            evidence = {'case': 'rects_cut', 'source_reference': SOURCE, 'oracle': ORACLE,
                        'native_sha256': sha(json.dumps(now['rects_cut'])), 'differences': ['rust_entity_counts']}
            linux = review_for(evidence, reviews)
            self.assertIs(linux, reviews[1])
            self.assertEqual((linux['reason'], linux['independent_evidence']), ('a reason', 'some evidence'))
            self.assertIn('run 123 at 01234567', linux['platform_observed'])
            # A second take does nothing.
            plans = fixture.take(run='123', commit='0123456789abcdef', write_files=True)
            self.assertEqual(plans[0].get('note'), 'already recorded')
            self.assertEqual(len(json.loads(fixture.reviews.read_text())['reviews']), 2)

    def test_a_difference_beyond_rounding_needs_a_review_by_hand(self):
        with tempfile.TemporaryDirectory() as directory:
            fixture = Fixture(Path(directory), native_text(counts='14 32 20 10 24 15', cut=ULP))
            plans = fixture.take(run='123', commit='0123456789abcdef', write_files=True)
            self.assertIsNone(plans[0]['record'])
            self.assertEqual(plans[0]['differences'], [('rects_fuse', 'counts')])
            self.assertFalse(plans[0]['written'])
            self.assertFalse((fixture.capture/'platform-linux').exists())
            self.assertIn('needs a record review', records.summary(plans))
            self.assertIn('rects_fuse (counts)', records.summary(plans))
            # The review written by hand takes the record.
            plans = fixture.take(run='123', commit='0123456789abcdef', write_files=True,
                                 reviews_given={fixture.capture.name: 'reviewed by hand'})
            self.assertTrue(plans[0]['written'])
            record = json.loads((fixture.capture/'platform-linux/record.json').read_text())
            self.assertEqual(record['review'], 'reviewed by hand')

    def test_a_new_difference_blocks_the_record_until_reviewed(self):
        with tempfile.TemporaryDirectory() as directory:
            fixture = Fixture(Path(directory), native_text(volume=875.0*(1+2**-50), cut=ULP),
                              extra_pending=[('rects_fuse', ['measure'])])
            plans = fixture.take(run='123', commit='0123456789abcdef', write_files=True)
            self.assertEqual([e['case'] for e in plans[0]['by_hand']], ['rects_fuse'])
            self.assertFalse(plans[0]['written'])
            self.assertFalse((fixture.capture/'platform-linux').exists())
            self.assertEqual(len(json.loads(fixture.reviews.read_text())['reviews']), 1)
            self.assertIn('needs reviews by hand', records.summary(plans))
            # A Linux review added by hand (its fingerprint the Linux row's)
            # resolves it.
            data = json.loads(fixture.reviews.read_text())
            pending = json.loads((fixture.artifact/'report.json').read_text())['platform_reviews_pending'][1]
            data['reviews'].append({k: pending[k] for k in records.EVIDENCE_KEYS} | {
                'reason': 'r', 'independent_evidence': 'e', 'platform_observed': 'Linux'})
            write(fixture.reviews, data)
            plans = fixture.take(run='123', commit='0123456789abcdef', write_files=True)
            self.assertEqual(plans[0]['reviewed'], ['rects_fuse'])
            self.assertTrue(plans[0]['written'])

    def test_stale_artifacts_are_refused(self):
        with tempfile.TemporaryDirectory() as directory:
            fixture = Fixture(Path(directory), native_text())
            (fixture.capture/'native.txt').write_text(native_text(volume=876.0))
            with self.assertRaisesRegex(ValueError, 'capture changed since the run'):
                fixture.take()
        with tempfile.TemporaryDirectory() as directory:
            fixture = Fixture(Path(directory), native_text())
            (fixture.artifact/'native-observed.txt').write_text(native_text(volume=876.0))
            with self.assertRaisesRegex(ValueError, 'not the report'):
                fixture.take()
        with tempfile.TemporaryDirectory() as directory:
            fixture = Fixture(Path(directory), native_text())
            with self.assertRaisesRegex(ValueError, 'needs --run and --commit'):
                fixture.take(write_files=True)


if __name__ == '__main__':
    unittest.main()
