"""Deliberate corruption and provenance failures in the native degree bridge."""
import copy
from fractions import Fraction as F
from pathlib import Path
import tempfile
import unittest
from compare_degree_elevation import compare, decode, verify_loaded_libraries, verify_capture, sha, digest, ROOT
import json
import knot_editing_reference as curves
import surface_knot_reference as surfaces


def native(kind, name, flags, shape):
    if kind == 'curve':
        p, periodic, knots, mults, controls = shape
        header = [name, 'R', len(flags), *map(int, flags), p, int(periodic),
                  len(controls), len(knots), *map(float, surfaces.domain(shape[:4]))]
        axes = [shape[:4]]
    else:
        axes, controls = shape
        header = [name, 'R', len(flags), *map(int, flags), *[a[0] for a in axes],
                  *[int(a[1]) for a in axes], *map(surfaces.count, axes),
                  *[len(a[2]) for a in axes], *[float(x) for a in axes for x in surfaces.domain(a)]]
    fields = [float(x) for c in controls for x in (c[0]/c[3], c[1]/c[3], c[2]/c[3], c[3])]
    fields += [x for a in axes for k, m in zip(a[2], a[3]) for x in (float(k), m)]
    return ' '.join(map(str, header+fields))+'\n'


class DegreeElevationOracleTests(unittest.TestCase):
    def shapes(self):
        axis = 2, False, (F(0), F(1)), (3, 3)
        controls = tuple((F(i), F(i*i), F(1-i), F(i%3+1)) for i in range(9))
        yield 'curve', (*axis, controls[:3]), curves.encode
        yield 'surface', ((axis, axis), controls), surfaces.encode

    def report(self, kind, shape, encode, record=None, rust=None, reviews=(), source='sample input', pin='PIN'):
        exact = encode('sample', [True], shape)+'\n'
        record = record if record is not None else {
            'case':'sample', 'exit_code':0, 'stdout':native(kind, 'sample', [True], shape), 'stderr':'OCCT TEST\n'}
        return compare(kind, source, [record], exact if rust is None else rust, exact,
                       'OCCT TEST', pin, reviews)

    def test_complete_grids_and_domains_match(self):
        for kind, shape, encode in self.shapes():
            result = self.report(kind, shape, encode)
            self.assertEqual(result['matched_cases'], 1)
            self.assertEqual(result['independently_verified'], 1)
            self.assertFalse(result['failures'])

    def test_every_control_component_and_flag_is_checked(self):
        for kind, shape, encode in self.shapes():
            for i in range(len(shape[-1])):
                for c in range(4):
                    grid = list(shape[-1]); row = list(grid[i]); row[c] += F(1, 100); grid[i] = tuple(row)
                    changed = (*shape[:-1], tuple(grid))
                    record = {'case':'sample', 'exit_code':0, 'stdout':native(kind, 'sample', [True], changed), 'stderr':'OCCT TEST\n'}
                    self.assertTrue(self.report(kind, shape, encode, record=record)['failures'])
                    self.assertEqual(self.report(kind, shape, encode, rust=encode('sample', [True], changed))['independently_verified'], 0)
            record['stdout'] = native(kind, 'sample', [False], shape)
            self.assertTrue(self.report(kind, shape, encode, record=record)['failures'])

    def test_exception_reviews_are_exact_and_never_count_as_matches(self):
        for kind, shape, encode in self.shapes():
            record = {'case':'sample', 'exit_code':0, 'stdout':'sample E\n', 'stderr':'OCCT TEST\nsample: error\n'}
            review = self.report(kind, shape, encode, record=record)['failures'][0]
            self.assertTrue(self.report(kind, shape, encode, record=record, reviews=[review])['failures'])
            review = dict(review, reason='Reviewed native range error on valid unclamped input',
                          independent_evidence='Full coefficient equations determine every control')
            good = self.report(kind, shape, encode, record=record, reviews=[review])
            self.assertEqual(good['matched_cases'], 0)
            self.assertEqual(good['independently_verified'], 1)
            self.assertEqual(len(good['reviewed_differences']), 1)
            self.assertFalse(good['failures'])
            for key in ('kind', 'case', 'oracle', 'native_source_reference', 'input_sha256',
                        'exact_sha256', 'native_stdout_sha256', 'native_stderr_sha256',
                        'classification', 'differences', 'reason', 'independent_evidence'):
                altered = copy.deepcopy(review); altered[key] = ''
                self.assertTrue(self.report(kind, shape, encode, record=record, reviews=[altered])['failures'], key)
            for altered in (dict(record, stdout='sample E'), dict(record, stderr=record['stderr']+'changed')):
                self.assertTrue(self.report(kind, shape, encode, record=altered, reviews=[review])['failures'])
            self.assertTrue(self.report(kind, shape, encode, record=record, reviews=[review], pin='OTHER')['failures'])
            self.assertTrue(self.report(kind, shape, encode, record=record, reviews=[review], source='sample changed')['failures'])
            self.assertTrue(self.report(kind, shape, encode, record=record, reviews=[review],
                                       rust=encode('sample', [False], shape))['failures'])

    def test_native_process_failures_cannot_be_exempted(self):
        for kind, shape, encode in self.shapes():
            for code in (-11, -9, 1, 77):
                record = {'case':'sample', 'exit_code':code, 'stdout':'sample E\n', 'stderr':'OCCT TEST\n'}
                review = self.report(kind, shape, encode, record=record)['failures'][0]
                result = self.report(kind, shape, encode, record=record, reviews=[review])
                self.assertTrue(result['failures'])
                self.assertFalse(result['reviewed_differences'])

    def test_missing_duplicate_and_extra_cases_fail(self):
        for kind, shape, encode in self.shapes():
            exact = encode('sample', [True], shape)
            record = {'case':'sample', 'exit_code':0, 'stdout':native(kind, 'sample', [True], shape), 'stderr':'OCCT TEST\n'}
            for records, rust in (([], exact), ([record, record], exact), ([record], ''),
                                  ([record], exact+'\n'+exact), ([record], exact.replace('sample', 'other'))):
                with self.assertRaises(ValueError):
                    compare(kind, 'sample input', records, rust, exact, 'OCCT TEST')

    def test_invalid_and_nonfinite_native_output_is_never_a_match(self):
        for kind, shape, encode in self.shapes():
            output = native(kind, 'sample', [True], shape)
            for text in ('', 'sample R', output+'extra', output+output,
                         output.replace('sample', 'wrong'), output.replace('0.0', 'nan', 1),
                         output.replace(' R 1 1 ', ' R 1 2 ', 1)):
                record = {'case':'sample', 'exit_code':0, 'stdout':text, 'stderr':'OCCT TEST\n'}
                result = self.report(kind, shape, encode, record=record)
                self.assertEqual(result['failures'][0]['classification'], 'invalid_native_output')
                self.assertEqual(result['matched_cases'], 0)
            words = output.split()
            # Domain is inconsistent with unchanged knots and dimensions.
            words[8 if kind == 'curve' else 12] = '-1.0'
            with self.assertRaises(ValueError):
                decode(kind, ' '.join(words), native=True)

    def test_link_evidence_must_resolve_all_occt_libraries_inside_sdk(self):
        with tempfile.TemporaryDirectory() as directory:
            lib = Path(directory)/'sdk lib'; lib.mkdir()
            other = Path(directory)/'other'; other.mkdir()
            for platform, suffix in [('darwin', '.8.1.0.dylib'), ('linux', '.so.8.1.0')]:
                paths = [lib/('lib'+name+suffix) for name in ('TKernel', 'TKMath', 'TKG2d', 'TKG3d')]
                for p in paths: p.write_bytes(b'test library')
                rows = [f'dyld[123]: <ID> {p}' if platform == 'darwin' else f'{p.name} => {p} (0x123)' for p in paths]
                self.assertEqual(len(verify_loaded_libraries('\n'.join(rows), lib, platform)), 4)
                for bad in ('\n'.join(rows[:3]), '\n'.join(rows).replace(str(lib), str(other), 1)):
                    with self.assertRaises(ValueError): verify_loaded_libraries(bad, lib, platform)

    def test_reused_capture_rejects_changed_source_corpus_sdk_and_observations(self):
        with tempfile.TemporaryDirectory() as directory:
            lib = Path(directory)/'library'; lib.write_bytes(b'original library')
            records = [{'case':'sample', 'exit_code':0, 'stdout':'sample E\n', 'stderr':'OCCT TEST\n'}]
            metadata = {'input_sha256':sha('sample input'), 'native_source_reference':'PIN',
                        'sdk_manifest_sha256':'SDK', 'observations_sha256':sha(json.dumps(records, sort_keys=True)),
                        'source_sha256':digest(ROOT/'rust/tools/occt_curve_degree_oracle.cpp'),
                        'loaded_libraries':[{'path':str(lib), 'sha256':digest(lib)}]}
            verify_capture(metadata, records, 'curve', 'sample input', 'PIN', 'SDK')
            for key in ('input_sha256', 'native_source_reference', 'sdk_manifest_sha256',
                        'observations_sha256', 'source_sha256', 'loaded_libraries'):
                altered = dict(metadata); altered[key] = ''
                with self.assertRaises(ValueError): verify_capture(altered, records, 'curve', 'sample input', 'PIN', 'SDK')
            lib.write_bytes(b'changed library')
            with self.assertRaises(ValueError): verify_capture(metadata, records, 'curve', 'sample input', 'PIN', 'SDK')


if __name__ == '__main__':
    unittest.main()
