import collections
import gzip
import hashlib
import os.path
import shutil
import tempfile
import unittest

from seqbc.run_file import (
    SequenceFile, NullSequenceFile, FastaSequenceFile,
    SplitBySampleFastqSequenceFile,
    IndexFastqSequenceFile,
    )


SFF_FP = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "data", "small.sff")


class MockWriter(object):
    def __init__(self):
        self.written = collections.defaultdict(list)

    def write(self, x, sample):
        if sample is None:
            self.written[None].append(x)
        else:
            self.written[sample.name].append(x)

            
MockSample = collections.namedtuple("MockSample", "name barcode")


class NullSequenceFileTests(unittest.TestCase):
    def test_no_demultiplex_method(self):
        """NullSequenceFile should not have a demultiplex method."""
        self.assertFalse(hasattr(NullSequenceFile("b"), "demultiplex"))


class FastaSequenceFileTests(unittest.TestCase):
    def test_get_runinfo(self):
        contents = ">a\nACGTACGT\n>b c d; e\nNNNNNN\n"
        
        f = tempfile.NamedTemporaryFile(suffix=".fasta")
        f.write(contents)
        f.seek(0)
        
        x = FastaSequenceFile(f.name)
        x.get_runinfo()

        self.assertEqual(x.num_reads, 2)
        self.assertEqual(x.num_flows, None)
        self.assertEqual(x.num_bases, 14)

        fasta_md5 = hashlib.md5()
        fasta_md5.update(contents)
        self.assertEqual(x.checksum, fasta_md5.hexdigest())


class IndexFastqSequenceFile(unittest.TestCase):
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
            
    def test_get_runinfo(self):
        x = SequenceFile(self.forward_fp)
        x.get_runinfo()
        self.assertEqual(x.num_reads, 3)
        self.assertEqual(x.num_flows, None)
        self.assertEqual(x.num_bases, 21 * 3 * 2)

        ifastq_md5 = hashlib.md5()
        ifastq_md5.update(open(self.forward_fp, "rb").read())
        self.assertEqual(x.checksum, ifastq_md5.hexdigest())

    def test_demultiplex(self):
        x = SequenceFile(self.forward_fp)
        w = MockWriter()
        # Barcode has 1 mismatch with second index read
        s1 = MockSample("SampleS1", "AGCGCCCT")
        x.demultiplex([s1], w)
        
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


if __name__ == "__main__":
    unittest.main()
