import collections
from io import StringIO
import os.path
import shutil
import tempfile
import unittest

from dnabclib.seqfile import (
    IndexFastqSequenceFile, NoIndexFastqSequenceFile, parse_fastq,
    )
from dnabclib.assigner import BarcodeAssigner


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
    def test_demultiplex(self):
        idx = StringIO(
            "@a\nACGTACGT\n+\n9812734[\n"
            "@b\nGGGGCGCT\n+\n78154987\n"
            "@c\nCCTTCCTT\n+\nkjafd;;;\n")
        fwd = StringIO(
            "@a\nGACTGCAGACGACTACGACGT\n+\n8A7T4C2G3CkAjThCeArG;\n"
            "@b\nCAGTCAGACGCGCATCAGATC\n+\n78154987bjhasf78612rb\n"
            "@c\nTCAGTACGTACGATACGTACG\n+\nkjafd;;;hjfasd82AHG99\n")
        rev = StringIO(
            "@a\nCATACGACGACTACGACTCAG\n+\nkjfhda987123GA;,.;,..\n"
            "@b\nGTNNNNNNNNNNNNNNNNNNN\n+\n#####################\n"
            "@c\nACTAGACTACGCATCAGCATG\n+\nkjafd;;;hjfasd82AHG99\n")
        x = IndexFastqSequenceFile(fwd, rev, idx)
        w = MockWriter()
        # Barcode has 0 mismatches with second index read
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


class NoIndexFastqSequenceFileTests(unittest.TestCase):
    def test_demultiplex(self):
        fwd = StringIO(fastq_with_barcode_fwd)
        rev = StringIO(fastq_with_barcode_rev)
        x = NoIndexFastqSequenceFile(fwd, rev)
        w = MockWriter()
        # Barcode matches the 4th read
        s1 = MockSample("SampleS1", "GTTTCGCCCTAGTACA")
        a = BarcodeAssigner([s1], mismatches=0, revcomp=False)
        x.demultiplex(a, w)

        # One read was written to SampleS1
        self.assertEqual(len(w.written["SampleS1"]), 1)
        # That read was the 4th read
        r1, r2 = w.written["SampleS1"][0]
        self.assertEqual(
            r1.desc,
            "HWI-D00727:9:C6JHHANXX:8:1101:1786:2183 1:N:0:GTTTCGCCCTAGTACA")
        self.assertEqual(
            r1.seq,
            "ACAATCAACCAACAAGTCATTATTACCCTCACTGTCAACCCAACACAGGCATGCTCATAAGGA"
            "AAGGTTAAAAAAAGTAAAAGGAACTCGGCAAATCTTACCCCGCCTGTTTACCAAAACCATCAC")
        self.assertEqual(
            r1.qual,
            "A3BBBDGGDEFBGBGG@F@FGEGGGGGGGGGGGGGGGEGGGDGBGBGBEFGCDG>GGGGGG1B"
            "FGGG1<FGGGGGGFGG0B>>F0FGGGGC/BBGGECF0@DGD@ADB000;=FGGGGBEBB/@@G")
        self.assertEqual(
            r2.desc,
            "HWI-D00727:9:C6JHHANXX:8:1101:1786:2183 2:N:0:GTTTCGCCCTAGTACA")
        self.assertEqual(
            r2.seq,
            "CATCTTACGCTGCCGACGATCTACTCTTTAGAAATGTCGTTCGTTTTGACTTCTGTAGAATAA"
            "GAATGTACTGCTCGGAGGTTGGGTTCTGCTCCGAGGTCGCCCCAACCGAAATTTTTAATGCAG")
        self.assertEqual(
            r2.qual,
            "33:A?11;@/;/;0//////001>11>111111111?10:E0=/1:/1/1111111=11111>"
            "11?:=FDEGBGGGG/EB<==@DDFGBEGC00C:>>D.FCG<CDGGGBGGBGGE=E..DGGE/C")


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

fastq_with_barcode_fwd = """\
@HWI-D00727:9:C6JHHANXX:8:1101:1361:2237 1:N:0:CTTACTAGAGACTACA
CTCCCACCTCATACCCCTCACCAGTTTGGGGGCTTGGTGTGGGCCCTTCCCTACAGATTGACGTCTTAACAGAGTTAGAACTTGGGGAGTCCACGACCGGAGTGCGTCAGCTCCCAGGGCTGCATA
+
A3@B@EFC1;1=CGG1FGGG>CBD:1:</EG/<<>/9FGD<//=FG01=FG>EB1111:111/0=00E0EDG0E0B0008;00000;...0C0<CDAA....:EEC.C.8/9D/C./66..CD@.6
@HWI-D00727:9:C6JHHANXX:8:1101:1926:2066 1:N:0:AAAAAACATATCCTCT
NAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
+
#=<BBGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGDGAGGGGGGGGGGGGGGGGGGGAGGGGGGGGGGGGGGGGGGGGGGGGGG.8.8..6.>.6.6.86CG:GGAG>ACGAGGGGGGG
@HWI-D00727:9:C6JHHANXX:8:1101:1912:2086 1:N:0:CTACGCTACTAAGCCT
AAGACGATCAGATACCGTCGTAGTTCCGACCATAAACGATGCCGACCGGCGATGCGGCGGCGTTATTCCCATGACCCGCCGGGCAGCTTCCGGGAAACCAAAGTCTTTGGGTTCCGGGGGGAGTAT
+
3<@BBGGGGGGGGGGGGGGBGGGGEGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGEDCEGGGGGDGGGGGGGGGGGGGDGGG@GEDGGGGGBEGEDDCGG<EG.CC>>.8.8C
@HWI-D00727:9:C6JHHANXX:8:1101:1786:2183 1:N:0:GTTTCGCCCTAGTACA
ACAATCAACCAACAAGTCATTATTACCCTCACTGTCAACCCAACACAGGCATGCTCATAAGGAAAGGTTAAAAAAAGTAAAAGGAACTCGGCAAATCTTACCCCGCCTGTTTACCAAAACCATCAC
+
A3BBBDGGDEFBGBGG@F@FGEGGGGGGGGGGGGGGGEGGGDGBGBGBEFGCDG>GGGGGG1BFGGG1<FGGGGGGFGG0B>>F0FGGGGC/BBGGECF0@DGD@ADB000;=FGGGGBEBB/@@G
@HWI-D00727:9:C6JHHANXX:8:1101:2413:2062 1:N:0:TTCTTGACTCTTTCCC
NGAGATGCTTGCTTATCAACAGAAGGAGTCTACTGCTCGCGTTGCGTCTATTATGGAAAACACCAATCTTTCCAAGCAACAGCAGGTTTCCGAGATTATGCGCCAAATGCTTACTCAAGCTCAAAC
+
#<?>BFGGGGGGGGGGGGGGFGGGGGGGFGGGGGGGFGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGFGGGGGGGGGGGGGGGFGGGGGGGGEGEGG=EE
"""

fastq_with_barcode_rev = """\
@HWI-D00727:9:C6JHHANXX:8:1101:1361:2237 2:N:0:CTTACTAGAGACTACA
GTCATCCATGAGGCAATACAAAGTTAGGACAAGCGGCAAGACCAACTCCTGCCGGGGTTCTCGGAACTGAAATCACAGTTGGGATGCAGCACGCCGGGCCTACCTTGAAAGAGAAAGGCTTTCGTC
+
BAAB:1EG1111>C/0E1111>:FC11E1E>FFBBGGB9<E1?=0/1<1FDC1/E////0<>/=C9E>F>0>0008000000.8C.C0DGC0..9.8>...C@D@/6/6CDD@./8/9.@DBGG@.
@HWI-D00727:9:C6JHHANXX:8:1101:1926:2066 2:N:0:AAAAAACATATCCTCT
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTATTTTTTT
+
CCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGDGGDGGGGGGGGDGGGDGDA>6C>>>CG.9.6.8.8...6/..666
@HWI-D00727:9:C6JHHANXX:8:1101:1912:2086 2:N:0:CTACGCTACTAAGCCT
AGGTTTCCCGTGTTGTCAAATTAAGCCGCAGGCTCCACTCCTGGTGGTGCCCTTCCGTCAATTCCTTTAAGTTTCAGCTTTGCAACCATACTCCCCCCGGAACCCAAAGACTTTGGTTTCCCGGAA
+
BCBBCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGEGGFGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGBGGGGGGGGGGGAD
@HWI-D00727:9:C6JHHANXX:8:1101:1786:2183 2:N:0:GTTTCGCCCTAGTACA
CATCTTACGCTGCCGACGATCTACTCTTTAGAAATGTCGTTCGTTTTGACTTCTGTAGAATAAGAATGTACTGCTCGGAGGTTGGGTTCTGCTCCGAGGTCGCCCCAACCGAAATTTTTAATGCAG
+
33:A?11;@/;/;0//////001>11>111111111?10:E0=/1:/1/1111111=11111>11?:=FDEGBGGGG/EB<==@DDFGBEGC00C:>>D.FCG<CDGGGBGGBGGE=E..DGGE/C
@HWI-D00727:9:C6JHHANXX:8:1101:2413:2062 2:N:0:TTCTTGACTCTTTCCC
CAAATTAGAGCCAATACCATCAGCTTTACCGTCTTTCCAGAAATTGTTCCAAGTATCGGCAACAGCTTTATCAATACCATGAAAAATATCAACCACACCAGAAGCAGCATCAGTGACGACATTAGA
+
CCCCCGGGGGGGGGGGGGGGGGDGGGGGGGGGFGGGGGGGGGGGGFGGGGGG>FGGGGDGGGGGGGGGGGGGGGGGGGGGGGGGGGGGFGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
"""

if __name__ == "__main__":
    unittest.main()
