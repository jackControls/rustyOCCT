"""The native comparison cannot silently discard points or hide a Rust error."""
import copy
import hashlib
import importlib.util
import unittest
from compare_spline_plane import cases, certificates, exact_certificate, native_fingerprint, observations, reviewed


class ParseTests(unittest.TestCase):
    def test_rejects_malformed_and_nonfinite_observations(self):
        for text in ['a -1 0','a 1 0 0 1 2 nan','a 0 0\na 0 0','a 0 1 0']:
            with self.assertRaises(ValueError): observations(text)
        for text in ['a -1 0','a 1 0 0 0 1 1 2 2 3 3 X 1 1','a 0 1 1 0',
                     'a 0 1 0 inf','a 0 0\na 0 0']:
            with self.assertRaises(ValueError): certificates(text)

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
