from collections import namedtuple
from cStringIO import StringIO
import os
import shutil
import tempfile
import unittest

from dnabc.assign import (
    PrefixAssigner, BarcodeAssigner,
    )


MockRead = namedtuple("Read", ["seq"])


class PrefixAssignerTests(unittest.TestCase):
    def test_assign(self):
        MockSample = namedtuple("Sample", ["name", "prefixes"])
        s = MockSample(123, ["AGGC"])
        a = PrefixAssigner([s])

        self.assertFalse(a.has_reads(s))

        # Prefix does not match, not assigned to sample
        self.assertEqual(a.assign(MockRead("ATTCCTT")), None)
        
        self.assertFalse(a.has_reads(s))

        self.assertEqual(a.assign(MockRead("AGGCCTT")), s)
        self.assertTrue(a.has_reads(s))


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
        MockSample = namedtuple("Sample", ["name", "barcode"])
        s = MockSample(123, "ACCTGAC")
        a = BarcodeAssigner([s], mismatches=1, revcomp=True)

        # 0 mismatches
        self.assertEqual(a.assign(MockRead("GTCAGGT")), s)
        # 1 mismatch
        self.assertEqual(a.assign(MockRead("GTCAAGT")), s)
        # 2 mismatches
        self.assertEqual(a.assign(MockRead("GTCAAAT")), None)

                
if __name__ == "__main__":
    unittest.main()
