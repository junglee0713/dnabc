import collections
from cStringIO import StringIO
import gzip
import os.path
import shutil
import tempfile
import unittest

from dnabclib.seqfile import (
    IndexFastqSequenceFile, parse_fastq,
    )
from dnabclib.assign import BarcodeAssigner


class MockWriter(object):
    def __init__(self):
        self.written = collections.defaultdict(list)

    def write(self, x, sample):
        if sample is None:
            self.written[None].append(x)
        else:
            self.written[sample.name].append(x)


MockSample = collections.namedtuple("MockSample", "name barcode")


class IndexFastqSequenceFileTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.index_fp = os.path.join(
            self.temp_dir, "Undetermined_S0_L001_I1_001.fastq.gz")
        self.index_contents = (
            "@a\nACGTACGT\n+\n9812734[\n"
            "@b\nGGGGCGCT\n+\n78154987\n"
            "@c\nCCTTCCTT\n+\nkjafd;;;\n")
        with gzip.open(self.index_fp, "w") as f:
            f.write(self.index_contents)

        self.forward_fp = os.path.join(
            self.temp_dir, "Undetermined_S0_L001_R1_001.fastq.gz")
        with gzip.open(self.forward_fp, "w") as f:
            f.write(
                "@a\nGACTGCAGACGACTACGACGT\n+\n8A7T4C2G3CkAjThCeArG;\n"
                "@b\nCAGTCAGACGCGCATCAGATC\n+\n78154987bjhasf78612rb\n"
                "@c\nTCAGTACGTACGATACGTACG\n+\nkjafd;;;hjfasd82AHG99\n")

        self.reverse_fp = os.path.join(
            self.temp_dir, "Undetermined_S0_L001_R2_001.fastq.gz")
        with gzip.open(self.reverse_fp, "w") as f:
            f.write(
                "@a\nCATACGACGACTACGACTCAG\n+\nkjfhda987123GA;,.;,..\n"
                "@b\nGTNNNNNNNNNNNNNNNNNNN\n+\n#####################\n"
                "@c\nACTAGACTACGCATCAGCATG\n+\nkjafd;;;hjfasd82AHG99\n")

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def test_demultiplex(self):
        x = IndexFastqSequenceFile(self.forward_fp, self.reverse_fp, self.index_fp)
        w = MockWriter()
        # Barcode has 1 mismatch with second index read
        s1 = MockSample("SampleS1", "GGGGCGCT")
        a = BarcodeAssigner([s1], mismatches=0, revcomp=False)
        x.demultiplex(a, w)

        # One read was written to SampleS1
        self.assertEqual(len(w.written["SampleS1"]), 1)
        # That read was the second of three above
        r1, r2 = w.written["SampleS1"][0]
        self.assertEqual(r1.desc, "b")
        self.assertEqual(r1.seq, "CAGTCAGACGCGCATCAGATC")
        self.assertEqual(r1.qual, "78154987bjhasf78612rb")
        self.assertEqual(r2.desc, "b")
        self.assertEqual(r2.seq, "GTNNNNNNNNNNNNNNNNNNN")
        self.assertEqual(r2.qual, "#####################")


class FunctionTests(unittest.TestCase):
    def test_parse_fastq(self):
        obs = parse_fastq(StringIO(fastq1))
        self.assertEqual(next(obs), (
            "YesYes", "AGGGCCTTGGTGGTTAG", ";234690GSDF092384"))
        self.assertEqual(next(obs), (
            "Seq2:with spaces", "GCTNNNNNNNNNNNNNNN", "##################"))
        self.assertRaises(StopIteration, next, obs)


fastq1 = """\
@YesYes
AGGGCCTTGGTGGTTAG
+
;234690GSDF092384
@Seq2:with spaces
GCTNNNNNNNNNNNNNNN
+
##################
"""


if __name__ == "__main__":
    unittest.main()
