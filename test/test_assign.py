from collections import namedtuple
from cStringIO import StringIO
import os
import shutil
import tempfile
import unittest

from dnabclib.assign import BarcodeAssigner


MockRead = namedtuple("Read", "seq")
MockSample = namedtuple("Sample", "name barcode")


class BarcodeAssignerTests(unittest.TestCase):
    def test_error_barcodes(self):
        a = BarcodeAssigner([], mismatches=1)
        obs = a._error_barcodes("AGG")
        exp = [
            "CGG", "GGG", "TGG",
            "AAG", "ACG", "ATG",
            "AGA", "AGC", "AGT",
            ]
        self.assertEqual(set(obs), set(exp))

    def test_one_mismatch(self):
        s = MockSample("Abc", "ACCTGAC")
        a = BarcodeAssigner([s], mismatches=1, revcomp=True)
        self.assertEqual(a.read_counts, {"Abc": 0})

        # 0 mismatches
        self.assertEqual(a.assign(MockRead("GTCAGGT")), s)
        self.assertEqual(a.read_counts, {"Abc": 1})

        # 1 mismatch
        self.assertEqual(a.assign(MockRead("GTCAAGT")), s)
        self.assertEqual(a.read_counts, {"Abc": 2})

        # 2 mismatches
        self.assertEqual(a.assign(MockRead("GTCAAAT")), None)
        self.assertEqual(a.read_counts, {"Abc": 2})

                
if __name__ == "__main__":
    unittest.main()
