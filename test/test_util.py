import collections
from cStringIO import StringIO
import unittest

from dnabc.util import (
    parse_fasta, parse_fastq, deambiguate, reverse_complement,
    )

class UtilTests(unittest.TestCase):
    def test_parse_fasta(self):
        obs = parse_fasta(StringIO(fasta1))
        self.assertEqual(next(obs), ("seq1 hello", "ACGTGGGTTAA"))
        self.assertEqual(next(obs), ("seq 2", "GTTCCGAAA"))
        self.assertEqual(next(obs), ("seq3", ""))
        self.assertRaises(StopIteration, next, obs)

    def test_parse_fastq(self):
        obs = parse_fastq(StringIO(fastq1))
        self.assertEqual(next(obs), (
            "YesYes", "AGGGCCTTGGTGGTTAG", ";234690GSDF092384"))
        self.assertEqual(next(obs), (
            "Seq2:with spaces", "GCTNNNNNNNNNNNNNNN", "##################"))
        self.assertRaises(StopIteration, next, obs)

    def test_deambiguate(self):
        obs = set(deambiguate("AYGR"))
        exp = set(["ACGA", "ACGG", "ATGA", "ATGG"])
        self.assertEqual(obs, exp)

        obs = set(deambiguate("AGN"))
        exp = set(["AGA", "AGC", "AGG", "AGT"])
        self.assertEqual(obs, exp)

    def test_reverse_complement(self):
        self.assertEqual(reverse_complement("AGATC"), "GATCT")
        self.assertRaises(KeyError, reverse_complement, "ANCC")


fasta1 = """\
>seq1 hello
ACGTGG
GTTAA
>seq 2
GTTC
C
GAAA
>seq3
"""

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

if __name__ == '__main__':
    unittest.main()
