"""Deliberate corruption checks: native reviews cannot exempt exact Rust results."""
import copy
from fractions import Fraction as F
import unittest
from bezier_editing_reference import decode,encode as result,expected
from compare_bezier_editing import compare,encode


def native(name,curves):
    words=[name,'R',str(len(curves))]
    for d,a,b,controls in curves:
        words.extend([str(d),str(float(a)),str(float(b))])
        for x,y,z,w in controls: words.extend(str(float(v)) for v in [x/w,y/w,z/w,w])
    return ' '.join(words)+'\n'


class BezierEditingOracleTests(unittest.TestCase):
    def setUp(self):
        self.input=encode('quadratic','B',2,[(0.,0.,0.),(1.,2.,3.),(2.,0.,1.)],[1.,2.,1.],[0.,1.],[3,3],0.,1.,5,5)+'\n'
        self.expected=expected(self.input)
        self.rust=result('quadratic',self.expected)+'\n'
        self.native=native('quadratic',self.expected)

    def report(self,n=None,r=None,reviews=()):
        return compare(self.input,self.native if n is None else n,self.rust if r is None else r,'TEST',reviews)

    def test_complete_exact_curves_match(self):
        r=self.report(); self.assertFalse(r['failures']); self.assertEqual(r['independently_verified'],1)
        self.assertEqual(r['matched_cases'],1); self.assertEqual(r['arcs'],2)

    def test_interior_control_corruption_fails_even_with_same_endpoints(self):
        changed=copy.deepcopy(self.expected); changed[0][3][2][0]+=F(1,1000)
        self.assertEqual(self.report(r=result('quadratic',changed))['independently_verified'],0)
        self.assertTrue(self.report(n=native('quadratic',changed))['failures'])

    def test_dropped_arc_reversal_and_parameter_corruption_fail(self):
        for curves in [self.expected[:1],self.expected[::-1]]:
            self.assertTrue(self.report(n=native('quadratic',curves))['failures'])
        changed=copy.deepcopy(self.expected); d,a,b,c=changed[0]; changed[0]=(d,a+F(1,100),b,c)
        self.assertTrue(self.report(n=native('quadratic',changed))['failures'])

    def test_homogeneous_normalization_preserves_native_geometry(self):
        changed=[(d,a,b,[[2*x for x in c] for c in controls]) for d,a,b,controls in self.expected]
        self.assertFalse(self.report(n=native('quadratic',changed))['failures'])
        # Production promises raw exact homogeneous data, a stronger contract.
        self.assertTrue(self.report(r=result('quadratic',changed))['failures'])

    def test_reviews_pin_all_evidence_and_never_override_rust(self):
        changed=copy.deepcopy(self.expected); changed[0][3][2][0]+=1
        n=native('quadratic',changed); review=self.report(n=n)['failures'][0]
        self.assertFalse(self.report(n=n,reviews=[review])['failures'])
        for field in ['case','oracle','input_sha256','native','exact_sha256','differences']:
            wrong=copy.deepcopy(review); wrong[field]='altered'
            self.assertTrue(self.report(n=n,reviews=[wrong])['failures'],field)
        self.assertTrue(self.report(n=n,r=result('quadratic',changed),reviews=[review])['failures'])
        changed[0][3][2][0]+=1
        self.assertTrue(self.report(n=native('quadratic',changed),reviews=[review])['failures'])

    def test_missing_duplicate_malformed_nonfinite_and_nonpositive_data_fail(self):
        for n in ['',self.native+self.native,self.native+' extra',self.native.replace('0.25','nan',1),'quadratic E\n']:
            with self.assertRaises((ValueError,StopIteration)): self.report(n=n)
        for r in [self.rust+self.rust,self.rust+' extra',self.rust.replace('quadratic R 2','quadratic R 0'),self.rust.replace('1/4','1/0',1)]:
            with self.assertRaises((ValueError,StopIteration,ZeroDivisionError)): decode(r)
        changed=copy.deepcopy(self.expected); changed[0][3][0][3]=-1
        with self.assertRaises(ValueError): decode(result('quadratic',changed))


if __name__=='__main__': unittest.main()
