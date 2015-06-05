import gzip
import os.path
import shutil
import tempfile
import unittest

from dnabc.models import (
    Sample,
    )


class SampleTests(unittest.TestCase):
    def test_prefixes(self):
        s = Sample("a", "agct")
        self.assertEqual(s.barcode, "AGCT")


if __name__ == "__main__":
    unittest.main()
