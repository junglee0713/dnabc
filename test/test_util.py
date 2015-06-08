import collections
from cStringIO import StringIO
import unittest

from dnabclib.util import (
    deambiguate, reverse_complement,
    )

class UtilTests(unittest.TestCase):
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


if __name__ == '__main__':
    unittest.main()
