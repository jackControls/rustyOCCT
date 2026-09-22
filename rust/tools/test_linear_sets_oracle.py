"""Deliberate-corruption tests for the native and exact set bridge."""
import copy
import hashlib
import unittest
from compare_linear_sets import compare, parse_native, parse_rust, fingerprint
from compare_proximity import encode
from linear_sets_reference import intersection, encode as encode_exact
from proximity_reference import exact

class LinearSetsOracleTests(unittest.TestCase):
    def setUp(self):
        self.a=('T',[(0.,0.,0.),(4.,0.,0.),(0.,4.,0.)])
        self.b=('T',[(.5,.5,0.),(1.,.5,0.),(.5,1.,0.)])
        self.input=encode('contained',self.a,self.b)+'\n'
        self.expected=intersection(exact(self.a),exact(self.b))
        self.rust='contained R '+encode_exact(self.expected)+'\n'
        self.native='contained C R 1 G 3 0.5 0.5 0 1 0.5 0 0.5 1 0\ncontained S R 0\n'
    def report(self,native=None,rust=None,reviews=()):
        return compare(self.input,self.native if native is None else native,self.rust if rust is None else rust,'TEST',reviews)
    def test_full_face_passes_and_reordered_native_vertices_are_allowed(self):
        report=self.report()
        self.assertFalse(report['failures'])
        self.assertEqual(report['combined_complete_matches'],1)
    def test_boundary_is_not_mistaken_for_filled_area(self):
        native=('contained C R 0\ncontained S R 3 '
                'S 0.5 0.5 0 1 0.5 0 S 1 0.5 0 0.5 1 0 S 0.5 1 0 0.5 0.5 0\n')
        r=self.report(native=native)
        self.assertEqual(r['section_subset_matches'],1)
        self.assertEqual(r['combined_complete_matches'],0)
        self.assertTrue(r['failures'])
    def test_missing_or_extra_extreme_vertex_fails(self):
        for value in ['0.75','3']:
            self.assertTrue(self.report(native=self.native.replace('1 0.5 0',value+' 0.5 0'))['failures'])
    def test_wrong_native_plane_and_missing_segment_interior_fail(self):
        for native in ['contained C R 1 F 0 0 1 0 0 1\ncontained S R 0\n',
                       'contained C R 1 G 3 0.5 0.5 1 1 0.5 1 0.5 1 1\ncontained S R 0\n']:
            self.assertTrue(self.report(native=native)['failures'])
        a=('S',[(0.,0.,0.),(4.,0.,0.)]); b=('S',[(1.,0.,0.),(3.,0.,0.)])
        data=encode('overlap',a,b)+'\n'; rust='overlap R '+encode_exact(intersection(exact(a),exact(b)))+'\n'
        native='overlap C R 2 P 1 0 0 P 3 0 0\noverlap S R 0\n'
        self.assertTrue(compare(data,native,rust,'TEST')['failures'])
    def test_dimension_filtering_is_counted_separately(self):
        b=('T',[(4.,0.,0.),(6.,0.,0.),(4.,2.,0.)])
        data=encode('touch',self.a,b)+'\n'; rust='touch R '+encode_exact(intersection(exact(self.a),exact(b)))+'\n'
        report=compare(data,'touch C R 0\ntouch S R 1 P 4 0 0\n',rust,'TEST')
        self.assertFalse(report['failures'])
        self.assertEqual(report['common_lower_dimension_omissions'],1)
    def test_reviews_pin_input_version_observation_exact_result_and_differences(self):
        native='contained C R 0\ncontained S R 0\n'
        report=self.report(native=native)
        obs=parse_native(native)
        review={'case':'contained','oracle':'TEST','input_sha256':hashlib.sha256(self.input.strip().encode()).hexdigest(),
                'native':fingerprint({op:obs['contained',op] for op in ['C','S']}),'exact_result':encode_exact(self.expected),
                'differences':report['failures'][0]['reason'],'reason':'synthetic guard test'}
        self.assertFalse(self.report(native=native,reviews=[review])['failures'])
        for field,value in [('case','wrong'),('oracle','TEST2'),('input_sha256','0'),('native',{}),('exact_result','E 0'),('differences',[])]:
            changed=copy.deepcopy(review); changed[field]=value
            self.assertTrue(self.report(native=native,reviews=[changed])['failures'],field)
        wrong='contained R E 0\n'
        self.assertTrue(self.report(native=native,rust=wrong,reviews=[review])['failures'])
        self.assertEqual(self.report(native=native,rust=wrong,reviews=[review])['independently_verified'],0)
    def test_malformed_duplicate_nonfinite_and_missing_observations_fail(self):
        for native in [self.native+self.native,self.native+'contained A R 0\n',self.native.replace('0.5','nan',1),self.native+'extra',self.native.replace('R 1 G 3','R 1 G 2')]:
            with self.assertRaises((ValueError,StopIteration)): self.report(native=native)
        with self.assertRaises(ValueError): self.report(native='contained C R 0\n')
        with self.assertRaises(ValueError): parse_rust(self.rust+self.rust)
        with self.assertRaises(ValueError): parse_rust(self.rust+'contained R E 1 0 0 0\n')

if __name__=='__main__': unittest.main()
