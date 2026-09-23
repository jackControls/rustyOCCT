"""Independent equation-solver cross-checks and deliberate bridge corruptions."""
import copy
from fractions import Fraction as F
import unittest
from compare_knot_editing import compare,encode,cases
import knot_editing_reference as reference


def native(name,flags,curve):
    p,periodic,knots,mults,controls=curve
    _,_,_,a,b,n=reference.data(curve)
    w=[name,'R',str(len(flags)),*[str(int(x)) for x in flags],str(p),str(int(periodic)),str(n),str(len(knots)),str(float(a)),str(float(b))]
    w += [str(float(x)) for c in controls for x in [c[0]/c[3],c[1]/c[3],c[2]/c[3],c[3]]]
    w += [x for k,m in zip(knots,mults) for x in [str(float(k)),str(m)]]
    return ' '.join(w)+'\n'


class KnotEditingOracleTests(unittest.TestCase):
    def setUp(self):
        shape=('quadratic','S',2,[(0.,0.,0.),(1.,2.,3.),(2.,0.,1.)],[1.,2.,1.],[0.,1.],[3,3])
        self.source=encode(shape,'edited',[('I',.375,2),('R',.375,0),('I',.625,1)])+'\n'
        self.name=self.source.split()[0]
        self.flags,self.curve=reference.expected(self.source)
        self.exact=reference.encode(self.name,self.flags,self.curve)+'\n'
        self.native=native(self.name,self.flags,self.curve)

    def report(self,n=None,r=None,reviews=()):
        return compare(self.source,self.native if n is None else n,self.exact if r is None else r,self.exact,'TEST',reviews)

    def test_complete_results_match(self):
        result=self.report()
        self.assertEqual(result['independently_verified'],1)
        self.assertEqual(result['matched_cases'],1)
        self.assertFalse(result['failures'])

    def test_control_and_success_flag_corruptions_fail(self):
        p,periodic,knots,mults,controls=self.curve
        changed=list(controls); changed[1]=tuple(x+F(i==0,1000) for i,x in enumerate(changed[1]))
        altered=p,periodic,knots,mults,tuple(changed)
        self.assertTrue(self.report(n=native(self.name,self.flags,altered))['failures'])
        self.assertEqual(self.report(r=reference.encode(self.name,self.flags,altered))['independently_verified'],0)
        self.assertTrue(self.report(n=native(self.name,[True,False,True],self.curve))['failures'])
        self.assertTrue(self.report(r=reference.encode(self.name,[True,False,True],self.curve))['failures'])

    def test_reviews_pin_every_observation_and_never_exempt_rust(self):
        bad=native(self.name,[True,False,True],self.curve)
        review=self.report(n=bad)['failures'][0]
        self.assertFalse(self.report(n=bad,reviews=[review])['failures'])
        for key in ['case','oracle','input_sha256','native','exact_sha256','differences']:
            wrong=copy.deepcopy(review); wrong[key]='altered'
            self.assertTrue(self.report(n=bad,reviews=[wrong])['failures'],key)
        corrupt=reference.encode(self.name,[True,False,True],self.curve)
        self.assertTrue(self.report(n=bad,r=corrupt,reviews=[review])['failures'])

    def test_missing_duplicate_nonfinite_and_malformed_observations_fail(self):
        for data in ['',self.native+self.native,self.native+' extra',self.name+' E\n',self.native.replace('0.0','nan',1),self.native.replace(' R 3 1',' R 3 2',1)]:
            with self.assertRaises((ValueError,StopIteration)): self.report(n=data)

    def test_coefficient_solver_agrees_with_separate_greville_reconstruction(self):
        _,original,_=reference.parse(self.source)
        alternative=reference.recover_by_collocation(original,self.curve[2],self.curve[3])
        self.assertEqual(alternative,self.curve)
        self.assertTrue(reference.equal(original,self.curve))
        row=next(row for row in cases().splitlines() if row.startswith('nonconstant_redundant_seam_remove_origin '))
        _,original,_=reference.parse(row)
        flags,result=reference.expected(row)
        self.assertEqual(flags,[True])
        self.assertEqual(result,reference.recover_by_collocation(original,result[2],result[3]))

    def test_inactive_exterior_controls_and_nonpositive_weights_are_checked(self):
        shape=('inactive','S',2,[(97.,41.,-93.),(0.,0.,0.),(1.,2.,3.),(4.,3.,2.)],[3.,1.,2.,1.],[-1.,0.,1.,2.],[2,2,1,2])
        _,original,_=reference.parse(encode(shape,'identity',[]))
        p,periodic,knots,mults,controls=original
        changed=list(controls); changed[0]=tuple(x+F(i==0) for i,x in enumerate(changed[0]))
        changed=p,periodic,knots,mults,tuple(changed)
        self.assertEqual(reference.value(original,F(1,2)),reference.value(changed,F(1,2)))
        self.assertFalse(reference.equal(original,changed))
        fine=(2,False,(F(0),F(1,2),F(1)),(3,1,3),tuple((F(0),F(0),F(0),w) for w in [F(1),F(3,8),F(3,8),F(1)]))
        self.assertFalse(reference.edited(fine,'R',F(1,2),0)[0])


if __name__=='__main__': unittest.main()
