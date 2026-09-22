"""Native spline discrepancy reviews cannot bypass independent arithmetic."""
import copy
import hashlib
import math
import struct
import unittest
from unittest.mock import patch
from review_spline_jets import exact_hash, exact_jet, minimal_bounds, reviewed_jet
from compare_splines import cases as curves
from compare_surfaces import cases as surfaces


class JetReviewTests(unittest.TestCase):
    def setUp(self):
        self.rows={'curves':next(r for r in curves().splitlines() if r.startswith('periodic_d25_m25_1 ')),
                   'surfaces':next(r for r in surfaces().splitlines() if r.startswith('d25_1_p00_1 '))}

    def exercise(self,family):
        row=self.rows[family]
        exact=exact_jet(family,row)
        rust=[minimal_bounds(x)[0] for x in exact]
        native=rust.copy(); native[0]+=0.001
        review={'id':'synthetic-guard-test','family':family,'oracle':'OCCT test','case':row.split()[0],
                'input_sha256':hashlib.sha256(row.encode()).hexdigest(),'exact_jet_sha256':exact_hash(exact),
                'native_bits':[struct.pack('>d',x).hex() for x in native],'evidence':'test','reason':'test'}
        args=[family,'OCCT test',row.split()[0],native,rust,row]
        with patch('review_spline_jets.reviews',return_value=[review]):
            self.assertIsNotNone(reviewed_jet(*args))
            for i,value in [(0,'different'),(1,'OCCT other'),(2,'wrong-case'),(5,row+' ')]:
                changed=copy.deepcopy(args);changed[i]=value
                self.assertIsNone(reviewed_jet(*changed))
            for i in [3,4]:
                changed=copy.deepcopy(args);changed[i][0]=math.nextafter(changed[i][0],-math.inf)
                self.assertIsNone(reviewed_jet(*changed))
            # Even an otherwise matching native component cannot conceal a Rust error.
            changed=copy.deepcopy(args);changed[4][3]=12345.
            self.assertIsNone(reviewed_jet(*changed))
            changed=copy.deepcopy(review);changed['exact_jet_sha256']='invalid'
            with patch('review_spline_jets.reviews',return_value=[changed]):
                self.assertIsNone(reviewed_jet(*args))

    def test_periodic_curve_review_guards(self): self.exercise('curves')
    def test_surface_review_guards(self): self.exercise('surfaces')


if __name__=='__main__': unittest.main()
