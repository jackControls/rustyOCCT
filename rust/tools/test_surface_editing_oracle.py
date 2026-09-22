"""Deliberate corruption: full tensor data, domains, native pins and exactness."""
import copy
from fractions import Fraction as F
import unittest
from surface_editing_reference import decode,encode_result,encode_fixture,expected
from compare_surface_editing import compare,encode


def native(name,items):
    words=[name,'R',str(len(items))]
    for kind,du,dv,domain,controls in items:
        words.extend(map(str,[kind,du,dv,*map(float,domain)]))
        for x,y,z,w in controls: words.extend(str(float(v)) for v in [x/w,y/w,z/w,w])
    return ' '.join(words)+'\n'


class SurfaceEditingOracleTests(unittest.TestCase):
    def setUp(self):
        shape=('patch','B',2,3,[(float(i),float(j),float(i*j)) for i in range(3) for j in range(4)],
               [float(1+i%3) for i in range(12)],[0.,1.],[3,3],[0.,1.],[4,4],False,False)
        self.input=encode(shape,(0.,1.,0.,1.),10)+'\n'
        self.expected=expected(self.input); self.rust=encode_result('patch',self.expected)+'\n'; self.native=native('patch',self.expected)

    def report(self,n=None,r=None,reviews=()):
        return compare(self.input,self.native if n is None else n,self.rust if r is None else r,'TEST',reviews)

    def test_complete_tensors_match(self):
        r=self.report(); self.assertFalse(r['failures']); self.assertEqual(r['independently_verified'],1)
        self.assertEqual(r['matched_cases'],1); self.assertEqual(r['patches'],4)

    def test_interior_corruption_preserving_boundaries_fails(self):
        c=copy.deepcopy(self.expected); c[0][4][c[0][2]+2][0]+=F(1,1000)
        self.assertEqual(self.report(r=encode_result('patch',c))['independently_verified'],0)
        self.assertTrue(self.report(n=native('patch',c))['failures'])

    def test_missing_reordered_transposed_or_wrong_domain_fails(self):
        for c in [self.expected[:1],self.expected[::-1]]: self.assertTrue(self.report(n=native('patch',c))['failures'])
        c=copy.deepcopy(self.expected); kind,du,dv,domain,controls=c[0]
        c[0]=(kind,dv,du,domain,controls)
        self.assertTrue(self.report(n=native('patch',c))['failures'])
        c[0]=(kind,du,dv,(domain[0]+F(1,100),*domain[1:]),controls)
        self.assertTrue(self.report(n=native('patch',c))['failures'])

    def test_projective_scaling_only_allowed_in_native_comparison(self):
        c=[(kind,du,dv,domain,[[2*x for x in p] for p in controls]) for kind,du,dv,domain,controls in self.expected]
        self.assertFalse(self.report(n=native('patch',c))['failures'])
        self.assertTrue(self.report(r=encode_result('patch',c))['failures'])
        self.assertNotEqual(encode_fixture('patch',c),encode_fixture('patch',self.expected))

    def test_reviews_pin_all_bits_and_never_exempt_exact_rust(self):
        c=copy.deepcopy(self.expected); c[0][4][2][0]+=1
        n=native('patch',c); review=self.report(n=n)['failures'][0]
        self.assertFalse(self.report(n=n,reviews=[review])['failures'])
        for field in ['case','oracle','input_sha256','native','exact_sha256','differences']:
            wrong=copy.deepcopy(review); wrong[field]='altered'
            self.assertTrue(self.report(n=n,reviews=[wrong])['failures'],field)
        self.assertTrue(self.report(n=n,r=encode_result('patch',c),reviews=[review])['failures'])
        c[-1][4][-1][0]+=F(1,10**12)
        self.assertTrue(self.report(n=native('patch',c),reviews=[review])['failures'])

    def test_native_exceptions_and_invalid_controls_are_preserved_failures(self):
        self.assertTrue(self.report(n='patch E\n')['failures'])
        for value in ['nan','inf','-1.0','0.0']:
            w=self.native.split(); w[13]=value # first control weight
            self.assertTrue(self.report(n=' '.join(w))['failures'])

    def test_malformed_and_duplicate_results_rejected(self):
        for text in ['',self.native+self.native,self.native+' extra','patch R 0']:
            with self.assertRaises((ValueError,StopIteration)): self.report(n=text)
        for text in [self.rust+self.rust,self.rust+' extra',self.rust.replace('1/4','1/0',1)]:
            with self.assertRaises((ValueError,StopIteration,ZeroDivisionError)): decode(text)
        c=copy.deepcopy(self.expected); c[0][4][0][3]=-1
        with self.assertRaises(ValueError): decode(encode_result('patch',c))

if __name__=='__main__': unittest.main()
