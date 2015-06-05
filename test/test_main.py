import gzip
import json
import os
import shutil
import tempfile
import unittest

from dnabc.main import main


class FastqDemultiplexTests(unittest.TestCase):
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

        self.barcode_fp = os.path.join(self.temp_dir, "manifest.txt")
        with open(self.barcode_fp, "w") as f:
            f.write(
                "SampleA\tAAGGAAGG\n"
                "SampleB\tACGTACGT\n")

        self.output_dir = os.path.join(self.temp_dir, "output")

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def test_regular(self):
        main([
            "-s", self.forward_fp, "-b", self.barcode_fp,
            "-o", self.output_dir, "-f", "fastq"])

        summary_fp = os.path.join(self.output_dir, "dnabc_summary.json")
        with open(summary_fp) as f:
            res = json.load(f)
            self.assertEqual(res["data"], {"SampleA": 1, "SampleB": 1})


if __name__ == "__main__":
    unittest.main()
