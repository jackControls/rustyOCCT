"""Full grid bridge rejection tests and independent solver cross-checks."""
import copy
from fractions import Fraction as F
import unittest
from compare_surface_knots import compare,encode
import surface_knot_reference as ref
import knot_editing_reference as curve_ref


def native(name,flags,surface):
    axes,controls=surface
    # Native output is Euclidean poles and weights, not homogeneous controls.
    w=[name,'R',len(flags),*map(int,flags),*[a[0] for a in axes],*[int(a[1]) for a in axes],*map(ref.count,axes),*[len(a[2]) for a in axes],*[float(x) for a in axes for x in ref.domain(a)]]
    w += [float(x) for c in controls for x in [c[0]/c[3],c[1]/c[3],c[2]/c[3],c[3]]]
    w += [x for a in axes for k,m in zip(a[2],a[3]) for x in [float(k),m]]
    return ' '.join(map(str,w))+'\n'


class SurfaceKnotOracleTests(unittest.TestCase):
    def setUp(self):
        shape=('tensor','S',2,2,[(float(i),float(j),float(i*j)) for i in range(3) for j in range(3)],[1.,2.,3.]*3,[0.,1.],[3,3],[0.,1.],[3,3],False,False)
        self.source=encode(shape,'edited',[('I','U',.375,2),('I','V',.625,1),('R','U',.375,0)])+'\n'
        self.name=self.source.split()[0]; self.flags,self.surface=ref.expected(self.source)
        self.exact=ref.encode(self.name,self.flags,self.surface,compact=True)+'\n'
        self.native=native(self.name,self.flags,self.surface)

    def report(self,n=None,r=None,reviews=()):
        return compare(self.source,self.native if n is None else n,self.exact if r is None else r,self.exact,'TEST',reviews)

    def test_complete_grids_match_and_both_serializations_are_exact(self):
        result=self.report(); self.assertEqual(result['independently_verified'],1); self.assertEqual(result['matched_cases'],1); self.assertFalse(result['failures'])
        self.assertEqual(ref.decode(self.exact),ref.decode(ref.encode(self.name,self.flags,self.surface)))

    def test_controls_weights_and_operation_flags_are_all_checked(self):
        axes,grid=self.surface
        for index in [0,len(grid)//2,len(grid)-1]:
            for component in range(4):
                changed=list(grid); changed[index]=tuple(x+F(c==component,100) for c,x in enumerate(changed[index])); altered=axes,tuple(changed)
                self.assertTrue(self.report(n=native(self.name,self.flags,altered))['failures'])
                self.assertEqual(self.report(r=ref.encode(self.name,self.flags,altered))['independently_verified'],0)
        self.assertTrue(self.report(n=native(self.name,[True,False,True],self.surface))['failures'])

    def test_pinned_reviews_cannot_hide_changed_results_or_bad_rust(self):
        bad=native(self.name,[True,False,True],self.surface); review=self.report(n=bad)['failures'][0]
        report=self.report(n=bad,reviews=[review]); self.assertFalse(report['failures']); self.assertEqual(report['matched_cases'],0)
        for key in ['case','oracle','input_sha256','native','exact_sha256','differences']:
            altered=copy.deepcopy(review); altered[key]='changed'
            self.assertTrue(self.report(n=bad,reviews=[altered])['failures'])
        self.assertTrue(self.report(n=bad,r=ref.encode(self.name,[True,False,True],self.surface),reviews=[review])['failures'])

    def test_missing_duplicate_malformed_and_nonfinite_data_fail(self):
        for text in ['',self.native+self.native,self.name+' E\n',self.native+' extra',self.native.replace('0.0','nan',1),self.native.replace(' R 3 1',' R 3 2',1)]:
            with self.assertRaises((ValueError,StopIteration)): self.report(n=text)

    def test_shared_factorization_matches_separate_curve_reconstruction(self):
        _,surface,ops=ref.parse(self.source)
        for op,which,u,target in ops:
            axes,grid=surface; nu,nv=map(ref.count,axes); flag,edited=ref.edit(surface,op,which,u,target)
            rows=[]; flags=[]
            for j in range([nv,nu][which]):
                controls=tuple(grid[i*nv+j] if which==0 else grid[j*nv+i] for i in range([nu,nv][which]))
                row_flag,row=curve_ref.edited((*axes[which],controls),op,u,target); flags.append(row_flag); rows.append(row)
            self.assertEqual(flag,all(flags))
            if flag:
                nu,nv=map(ref.count,edited[0])
                self.assertEqual(edited[1],tuple(rows[j][4][i] if which==0 else rows[i][4][j] for i in range(nu) for j in range(nv)))
            surface=edited

    def test_one_inconsistent_row_and_nonpositive_inverse_are_rejected(self):
        u=(1,False,(F(0),F(1),F(2)),(2,1,2)); v=(1,False,(F(0),F(1)),(2,2))
        controls=tuple((F(i),F(j),F(i==1 and j==1),F(1)) for i in range(3) for j in range(2))
        self.assertEqual(ref.edit(((u,v),controls),'R',0,F(1),0),(False,((u,v),controls)))
        u=(2,False,(F(0),F(1,2),F(1)),(3,1,3))
        controls=tuple((F(0),F(0),F(0),w) for w in [F(1),F(3,8),F(3,8),F(1)] for _ in range(2))
        self.assertFalse(ref.edit(((u,v),controls),'R',0,F(1,2),0)[0])


if __name__=='__main__': unittest.main()
