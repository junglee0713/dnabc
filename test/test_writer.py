from collections import namedtuple
import gzip
import os.path
import shutil
import tempfile
import unittest

from dnabc.writer import (
    FastaWriter, PooledFastaWriter,
    FastqWriter, PooledFastqWriter,
    PairedFastqWriter,
    )


class FastaWriterTests(unittest.TestCase):
    def setUp(self):
        self.output_dir = tempfile.mkdtemp()
        self.Sample = namedtuple("Sample", ["name", "run"])
        self.Read = namedtuple("Read", ["desc", "seq"])

    def tearDown(self):
        shutil.rmtree(self.output_dir)
    
    def test_write(self):
        s1 = self.Sample(14, 77)
        s2 = self.Sample(15, 77)
        w = FastaWriter([s1, s2], self.output_dir)

        w.write(self.Read("Read0", "ACCTTGG"), s1)
        w.close()

        fp = w.output_fps[s1.name]
        obs_output = open(fp).read()
        self.assertEqual(obs_output, ">Read0\nACCTTGG\n")

        self.assertFalse(os.path.exists(w.output_fps[s2.name]))

    def test_write_pooled(self):
        s1 = self.Sample(14, 77)
        s2 = self.Sample(15, 77)
        w = PooledFastaWriter([s1, s2], self.output_dir)
        
        self.assertEqual(w.output_fps[s1.name], w.output_fps[s2.name])

        w.write(self.Read("Read0", "ACCTTGG"), s1)
        w.close()

        obs_output = open(w.output_fps[s1.name]).read()
        self.assertEqual(obs_output, ">Read0\nACCTTGG\n")


class FastqWriterTests(unittest.TestCase):
    def setUp(self):
        self.output_dir = tempfile.mkdtemp()
        self.Sample = namedtuple("Sample", ["name", "run"])
        self.Read = namedtuple("Read", ["desc", "seq", "qual"])

    def tearDown(self):
        shutil.rmtree(self.output_dir)
    
    def test_write(self):
        s1 = self.Sample(14, 77)
        s2 = self.Sample(15, 77)
        w = FastqWriter([s1, s2], self.output_dir)

        w.write(self.Read("Read0", "ACCTTGG", "#######"), s1)
        w.close()

        fp = w.output_fps[s1.name]
        obs_output = gzip.open(fp).read()
        
        self.assertEqual(obs_output, "@Read0\nACCTTGG\n+\n#######\n")

        self.assertFalse(os.path.exists(w.output_fps[s2.name]))

    def test_write_pooled(self):
        s1 = self.Sample(14, 77)
        s2 = self.Sample(15, 77)
        w = PooledFastqWriter([s1, s2], self.output_dir)
        
        self.assertEqual(w.output_fps[s1.name], w.output_fps[s2.name])

        w.write(self.Read("Read0", "ACCTTGG", "#######"), s1)
        w.close()

        obs_output = gzip.open(w.output_fps[s1.name]).read()
        self.assertEqual(obs_output, "@Read0\nACCTTGG\n+\n#######\n")


class PairedFastqWriterTests(unittest.TestCase):
    def setUp(self):
        self.output_dir = tempfile.mkdtemp()
        self.Sample = namedtuple("Sample", ["name", "run"])
        self.Read = namedtuple("Read", ["desc", "seq", "qual"])

    def tearDown(self):
        shutil.rmtree(self.output_dir)
    
    def test_write(self):
        s1 = self.Sample(14, 77)
        s2 = self.Sample(15, 77)
        w = PairedFastqWriter([s1, s2], self.output_dir)

        readpair = (
            self.Read("Read0", "ACCTTGG", "#######"),
            self.Read("Read1", "GCTAGCT", ";342dfA"),
            )
        w.write(readpair, s1)
        w.close()

        fp1, fp2 = w.output_fps[s1.name]

        obs1 = gzip.open(fp1).read()
        self.assertEqual(obs1, "@Read0\nACCTTGG\n+\n#######\n")

        obs2 = gzip.open(fp2).read()
        self.assertEqual(obs2, "@Read1\nGCTAGCT\n+\n;342dfA\n")

        self.assertFalse(any(
            os.path.exists(fp) for fp in w.output_fps[s2.name]))


if __name__ == '__main__':
    unittest.main()
