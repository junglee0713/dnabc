import gzip
import os.path

 
def _sample_fps(self):
    fps = {}
    for s in self.samples:
        fn = "BLS%06d%s" % (s.name, self.ext)
        fps[s.name] = os.path.join(self.output_dir, fn)
    return fps

def _run_fps(self):
    if not self.samples:
        return {}
    s0 = self.samples[0]
    fn = "BLR%06d%s" % (s0.run, self.ext)
    fp = os.path.join(self.output_dir, fn)
    return dict((s.name, fp) for s in self.samples)

def _sample_paired_fps(self):
    fps = {}
    for s in self.samples:
        fn1 = "PCMP_%s_R1%s" % (s.name, self.ext)
        fn2 = "PCMP_%s_R2%s" % (s.name, self.ext)
        fps[s.name] = (
            os.path.join(self.output_dir, fn1),
            os.path.join(self.output_dir, fn2),
            )
    return fps

class _SequenceWriter(object):
    """Base class for writers"""

    def __init__(self, samples, output_dir):
        self.samples = samples
        self.output_dir = output_dir
        self.output_fps = self._get_output_fps()
        self._open_files = {}

    def set_sff_header(self, header):
        pass

    def _get_output_file(self, sample):
        fp = self.output_fps[sample.name]
        f = self._open_files.get(fp)
        if f is None:
            f = self._open_filepath(fp)
            self._open_files[fp] = f
        return f

    def _open_filepath(self, fp):
        return open(fp, "w")

    def write(self, read, sample):
        if sample is not None:
            f = self._get_output_file(sample)
            self._write_to_file(f, read)

    def close(self):
        for fp, f in self._open_files.items():
            f.close()


class FastaWriter(_SequenceWriter):
    ext = ".fasta"
    _get_output_fps = _sample_fps

    def _write_to_file(self, f, read):
        f.write(">%s\n%s\n" % (read.desc, read.seq))


class PooledFastaWriter(FastaWriter):
    _get_output_fps = _run_fps


class FastqWriter(_SequenceWriter):
    ext = ".fastq.gz"
    _get_output_fps = _sample_fps

    def _open_filepath(self, fp):
        return gzip.open(fp, "w")

    def _write_to_file(self, f, read):
        f.write("@%s\n%s\n+\n%s\n" % (read.desc, read.seq, read.qual))


class PooledFastqWriter(FastqWriter):
    _get_output_fps = _run_fps


class PairedFastqWriter(FastqWriter):
    _get_output_fps = _sample_paired_fps

    def _open_filepath(self, fps):
        fp1, fp2 = fps
        f1 = super(PairedFastqWriter, self)._open_filepath(fp1)
        f2 = super(PairedFastqWriter, self)._open_filepath(fp2)
        return (f1, f2)
    
    def _write_to_file(self, filepair, readpair):
        f1, f2 = filepair
        r1, r2 = readpair
        super(PairedFastqWriter, self)._write_to_file(f1, r1)
        super(PairedFastqWriter, self)._write_to_file(f2, r2)

    def close(self):
        for _, (f1, f2) in self._open_files.items():
            f1.close()
            f2.close()


