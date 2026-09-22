"""Deliberate corruptions must fail the proximity bridge and review guards."""
from copy import deepcopy
from fractions import Fraction as Q
import hashlib
import json
import unittest
from compare_proximity import ROOT, cases, compare, fingerprint, parse_native, parse_rust, reviewed
from proximity_reference import certify, distance2, exact, parse


class ProximityOracleTests(unittest.TestCase):
    data='one P 0 0 0 P 1 0 0\n'
    native='one P 1 1 0 0 0 1 0 0\n'
    affine='one P 1\n'
    rust='one R 1/1 E 3ff0000000000000 3ff0000000000000 E 3ff0000000000000 3ff0000000000000 0 0/1 0/1 0/1 0 1/1 0/1 0/1\n'

    def test_independent_certificate_and_native_counts(self):
        report=compare(self.data,self.native,self.affine,self.rust,'OCCT test',[])
        self.assertTrue(report['passed'])
        self.assertEqual(report['independently_verified'],1)
        self.assertEqual(report['brep'],{'matches':1,'nonresults':0,'witness_pairs':1})
        self.assertEqual(report['affine'],{'matches':1,'not_applicable':0})

    def test_wrong_rust_distance_point_and_interval_are_rejected(self):
        for changed in [self.rust.replace('R 1/1','R 4/1'),
                        self.rust.replace('0/1 0/1 0/1','1/1 0/1 0/1',1),
                        self.rust.replace('3ff0000000000000','3ff0000000000001',1)]:
            report=compare(self.data,self.native,self.affine,changed,'OCCT test',[])
            self.assertFalse(report['passed'])
            self.assertEqual(report['independently_verified'],0)

    def test_certificate_rejects_a_feasible_but_nonminimal_pair(self):
        shapes=[exact(('S',[(0.,0.,0.),(2.,0.,0.)])),exact(('P',[(3.,0.,0.)]))]
        with self.assertRaises(AssertionError):
            certify(shapes,[(Q(0),Q(0),Q(0)),(Q(3),Q(0),Q(0))],[[Q(0)],[]],Q(9))

    def test_native_witness_ownership_is_checked(self):
        report=compare(self.data,'one P 1 1 1 0 0 0 0 0\n',self.affine,self.rust,'OCCT test',[])
        self.assertFalse(report['passed'])
        self.assertEqual(report['failures'][0]['reason'],'native witness outside its operand')

    def test_nonresult_review_never_counts_as_native_match(self):
        review={'oracle':'OCCT test','api':'brep','case':'one','id':'synthetic','reason':'test',
                'input_sha256':hashlib.sha256(self.data.strip().encode()).hexdigest(),
                'native':{'status':'F','distance_bits':None},'exact_squared_distance':'1'}
        report=compare(self.data,'one F\n',self.affine,self.rust,'OCCT test',[review])
        self.assertTrue(report['passed'])
        self.assertEqual(report['brep']['matches'],0)
        self.assertEqual(report['brep']['nonresults'],1)
        self.assertEqual(len(report['reviewed_differences']),1)
        changed=self.rust.replace('R 1/1','R 4/1')
        self.assertFalse(compare(self.data,'one F\n',self.affine,changed,'OCCT test',[review])['passed'])

    def test_review_pins_version_input_native_bits_and_exact_answer(self):
        reviews=json.loads((ROOT/'rust/fixtures/occt-proximity-divergences.json').read_text())
        review=next(r for r in reviews if r['api']=='affine')
        row=next(r for r in cases().splitlines() if r.split()[0]==review['case'])
        import struct
        native=('P',struct.unpack('>d',bytes.fromhex(review['native']['distance_bits']))[0],[])
        _,shapes,_=parse(row)
        d=distance2(*map(exact,shapes))
        args=[review['oracle'],'affine',row,native,d,reviews]
        self.assertIsNotNone(reviewed(*args))
        for index,value in [(0,'OCCT unknown'),(1,'brep'),(2,row+' '),(3,('P',0.123,[])),(4,d+1)]:
            changed=deepcopy(args);changed[index]=value
            self.assertIsNone(reviewed(*changed))
        self.assertEqual(fingerprint(native),review['native'])

    def test_malformed_or_missing_results_fail_closed(self):
        for row in ['one P 1 0','one P nan 1 0 0 0 1 0 0','one F extra',self.native+self.native]:
            with self.assertRaises((ValueError,StopIteration,IndexError)): parse_native(row)
        for row in [self.rust+self.rust,self.rust+'extra',self.rust.replace('R 1/1','N 1/1')]:
            with self.assertRaises((ValueError,StopIteration)): parse_rust(row)
        with self.assertRaises(ValueError): compare(self.data,'',self.affine,self.rust,'OCCT test',[])
        self.assertFalse(compare(self.data,self.native,'one N\n',self.rust,'OCCT test',[])['passed'])


if __name__=='__main__': unittest.main()
