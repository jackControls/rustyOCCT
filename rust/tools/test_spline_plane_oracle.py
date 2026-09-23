"""The native comparison cannot silently discard points or hide a Rust error."""
import copy
import hashlib
import importlib.util
import unittest
from compare_spline_plane import cases, certificates, exact_certificate, native_fingerprint, observations, reviewed, CURVE_MODES, representation_differences


class ParseTests(unittest.TestCase):
    def test_exact_edit_modes_cannot_hide_a_missing_or_changed_result(self):
        exact={'points':[{'bounds':[0.,0.,1.,1.,0.,0.,0.,0.], 'contact':'C', 'orders':[1,1]}], 'overlaps':[]}
        outputs={mode:copy.deepcopy(exact) for mode in CURVE_MODES}
        self.assertEqual(representation_differences(exact,outputs),[])
        outputs['refined']['points'][0]['orders']=[2,2]
        self.assertEqual(representation_differences(exact,outputs),['refined'])
        outputs['exact']['points']=[]
        self.assertEqual(representation_differences(exact,outputs),['exact','refined'])
        del outputs['roundtrip']
        with self.assertRaises(ValueError): representation_differences(exact,outputs)

    def test_rejects_malformed_and_nonfinite_observations(self):
        for text in ['a -1 0','a 1 0 0 1 2 nan','a 0 0\na 0 0','a 0 1 0']:
            with self.assertRaises(ValueError): observations(text)
        for text in ['a -1 0','a 1 0 0 0 1 1 2 2 3 3 X 1 1','a 0 1 1 0',
                     'a 0 1 0 inf','a 0 0\na 0 0']:
            with self.assertRaises(ValueError): certificates(text)

    def test_degree_fifty_contact_orders_require_the_quadric_contract(self):
        row='a 1 0 0 0 0 0 1 1 0 0 B 0 50'
        with self.assertRaises(ValueError): certificates(row)
        self.assertEqual(certificates(row,max_order=50)['a']['points'][0]['orders'],[0,50])
        with self.assertRaises(ValueError): certificates(row+'1',max_order=50)

    @unittest.skipUnless(importlib.util.find_spec('sympy'),'install the pinned test-only mathematics oracle')
    def test_quadric_review_rejects_changed_geometry_or_removed_overlap(self):
        from compare_spline_quadric import cases as quadric_cases
        row=next(r for r in quadric_cases().splitlines() if r.startswith('sphere_rational_arc_0 '))
        exact,digest=exact_certificate(row)
        self.assertEqual(exact['overlaps'],[[0.,0.,1.,1.]])
        native={'points':[],'overlaps':[]}
        review={'id':'synthetic-test','oracle':'OCCT test','case':row.split()[0],
                'input_sha256':hashlib.sha256(row.encode()).hexdigest(),
                'native':native_fingerprint(native),'exact_certificate_sha256':digest}
        self.assertIsNotNone(reviewed('OCCT test',row,native,exact,[review]))
        self.assertIsNone(reviewed('OCCT test',row,native,native,[review]))
        words=row.split();words[11]='2'
        self.assertIsNone(reviewed('OCCT test',' '.join(words),native,exact,[review]))

    @unittest.skipUnless(importlib.util.find_spec('sympy'),'install the pinned test-only mathematics oracle')
    def test_review_pins_native_and_rechecks_complete_independent_answer(self):
        row=next(r for r in cases().splitlines() if r.startswith('contained_axis_1.0 '))
        exact,digest=exact_certificate(row)
        native={'points':[],'overlaps':[]}
        review={'id':'synthetic-test','oracle':'OCCT test','case':row.split()[0],
                'input_sha256':hashlib.sha256(row.encode()).hexdigest(),
                'native':native_fingerprint(native),'exact_certificate_sha256':digest}
        args=['OCCT test',row,native,copy.deepcopy(exact),[review]]
        self.assertIsNotNone(reviewed(*args))
        for i,v in [(0,'OCCT changed'),(1,row+' '),(2,None),(3,{'points':[],'overlaps':[]})]:
            changed=copy.deepcopy(args);changed[i]=v
            self.assertIsNone(reviewed(*changed))
        for field in ['input_sha256','exact_certificate_sha256','case']:
            changed=copy.deepcopy(args);changed[4][0][field]='wrong'
            self.assertIsNone(reviewed(*changed))
        changed=copy.deepcopy(args);changed[3]['overlaps'][0][1]=.999
        self.assertIsNone(reviewed(*changed))


if __name__=='__main__':unittest.main()
