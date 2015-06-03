import gzip
import hashlib
import itertools
import logging
import os
import platform
import subprocess

from .util import (
    parse_fasta, parse_fastq,
    FastaRead, FastqRead,
    )
from .assign import (
    PrefixAssigner, BarcodeAssigner,
    )


def SequenceFile(filepath):
    """Create a SequenceFile object from a filepath

    The filepath represents a data resource for a sequencing run, as
    registered in the CORE database.

    SequenceFile objects are used by the script that collects metadata on
    sequencing runs.  Thus, SequenceFile objects know how to count reads,
    count bases, and generate an md5sum of the data files.

    SequenceFile objects are also used by the daemon that processes
    downloads.  Thus, they know how to find a class that will assign
    reads to samples.

    SequenceFile is actually a factory function, though it looks like a class
    constructor.  It inspects the filepath of the data resource and
    returns an instance of one of the public classes in this module.    
    """
    if not os.path.exists(filepath):
        return NullSequenceFile(
            filepath, "Data resource %s does not exist" % filepath)

    if os.path.isdir(filepath):
        return SplitBySampleFastqSequenceFile(filepath)

    if os.path.basename(filepath) == "Undetermined_S0_L001_R1_001.fastq.gz":
        return IndexFastqSequenceFile(filepath)

    if filepath.endswith(".fasta") or filepath.endswith(".fna"):
        return FastaSequenceFile(filepath)

    return NullSequenceFile(
        filepath, "Invalid/unsupported data resource %s" % filepath)


class NullSequenceFile(object):
    """A run file that down't exist or is not supported.
    """
    filetype = None

    def __init__(self, filepath, comment=""):
        self.filepath = filepath
        self.comment = comment

        # The following attributes never apply to a null run
        self.num_reads = None
        self.num_bases = None
        self.num_flows = None
        self.checksum = None

    @property
    def filesize(self):
        pass

    def get_runinfo(self):
        pass

    def as_tuple(self):
        return (
            self.filepath, self.checksum, self.filetype, self.filesize,
            self.num_reads, self.num_bases, self.num_flows, self.comment)

    def __str__(self):
        return "<%s %s: %s>" % (self.__class__.__name__, self.filepath, self.filetype)


class _ExistingSequenceFile(NullSequenceFile):
    """Private class with common template methods for run files.
    """
    def demultiplex(self, samples, writer, assigner_class=PrefixAssigner):
        raise NotImplementedError
    
    @property
    def filesize(self):
        return os.path.getsize(self.filepath)

    def get_runinfo(self):
        try:
            logging.debug("Counting reads")
            self.num_reads = self._count_reads()
            logging.debug("Counting bases")
            self.num_bases = self._count_bases()
            self.num_flows = self._count_flows()
            logging.debug("Getting checksum")
            self.checksum = self._get_checksum()
        except IOError as e:
            self.comment = e

    def _count_reads(self):
        raise NotImplementedError

    def _count_bases(self):
        raise NotImplementedError
            
    def _count_flows(self):
        pass

    def _get_checksum_with_hashlib(self):
        # following advice on stackoverflow.com/questions/1131220
        md5sum = hashlib.md5()
        chunk_size = 128 * md5sum.block_size * (2 ** 20)
        with open(self.filepath, 'rb') as f:
            for chunk in iter(lambda: f.read(chunk_size), ''): 
                md5sum.update(chunk)
        return md5sum.hexdigest()

    def _get_checksum(self):
        commands = {
            "Darwin": ["md5", "-q"],
            "Linux": ["md5sum"],
            }
        myos = platform.system()
        if myos not in commands:
            return self.get_checksum_with_hashlib()
        c = commands[myos]
        c.append(self.filepath)
        p = subprocess.Popen(c, stdout=subprocess.PIPE)
        stdout, _ = p.communicate()
        firstword = stdout.split()[0]
        logging.debug(firstword)
        return firstword.strip()


class FastaSequenceFile(_ExistingSequenceFile):
    """FASTA data file from 454 sequencing."""
    filetype = "FASTA"

    def demultiplex(self, samples, writer, assigner_class=PrefixAssigner):
        assigner = assigner_class(samples)
        with open(self.filepath) as f:
            for read in parse_fasta(f):
                read = FastaRead(read)
                sample = assigner.assign(read)
                writer.write(read, sample)
        return assigner.read_counts

    def _count_reads(self):
        nr = 0
        with open(self.filepath) as f:
            for line in f:
                if line.startswith(">"):
                    nr += 1
        return nr

    def _count_bases(self):
        n = 0
        with open(self.filepath) as f:
            for line in f:
                if not line.startswith(">"):
                    n += len(line.strip())
        return n


class SplitBySampleFastqSequenceFile(_ExistingSequenceFile):
    """FASTQ data files already split by sample.

    The folder containing the files is registered as the data resource.
    """
    filetype = "Pre-split FASTQ"


class IndexFastqSequenceFile(_ExistingSequenceFile):
    """FASTQ data files from the Hahn lab SGA sequencing core.

    The method used is to enter 1 fake sample to the Illumina software,
    then ask it to retain reads not matching the barcode for the fake
    sample.  These reads, unassigned by the Illumina software, are the
    reads we want to use in our analyses.  The Illumina software assigns
    the following names for files with unassigned reads:

    * Undetermined_S0_L001_I1_001.fastq.gz - Index reads (barcodes)
    * Undetermined_S0_L001_R1_001.fastq.gz - Forward reads
    * Undetermined_S0_L001_R2_001.fastq.gz - Reverse reads

    For runs of this type, the file of forward reads is listed in the
    registry.
    """
    filetype = "Three-file FASTQ"
    
    @property
    def index_filepath(self):
        return self.filepath.replace(
            "Undetermined_S0_L001_R1", "Undetermined_S0_L001_I1", 1)

    @property
    def forward_filepath(self):
        return self.filepath

    @property
    def reverse_filepath(self):
        return self.filepath.replace(
            "Undetermined_S0_L001_R1", "Undetermined_S0_L001_R2", 1)

    def demultiplex(self, samples, writer, assigner_class=BarcodeAssigner):
        assigner = assigner_class(samples)

        idx_file = gzip.open(self.index_filepath, "rb")
        fwd_file = gzip.open(self.forward_filepath, "rb")
        rev_file = gzip.open(self.reverse_filepath, "rb")
        idxs = (FastqRead(x) for x in parse_fastq(idx_file))
        fwds = (FastqRead(x) for x in parse_fastq(fwd_file))
        revs = (FastqRead(x) for x in parse_fastq(rev_file))
        try:
            for idx, fwd, rev in itertools.izip(idxs, fwds, revs):
                sample = assigner.assign(idx)
                writer.write((fwd, rev), sample)
        finally:
            idx_file.close()
            fwd_file.close()
            rev_file.close()
        return assigner.read_counts
    
    def _count_reads(self):
        """For paired-end sequencing, returns the number of read pairs.
        """
        # 4 lines per FASTQ record, count lines and divide by 4.
        # Use floor division in case of extra whitespace at the end.
        with gzip.GzipFile(self.index_filepath) as f:
            for n, _ in enumerate(f):
                pass
            return (n + 1) // 4

    def _count_bases(self):
        # Number of bases is the same for all reads
        # Count bases in one read, multiply by (number of reads * 2)
        with gzip.GzipFile(self.forward_filepath) as f:
            next(f) # Skip first line
            bases_per_read = len(next(f)) - 1 # Second line has sequence + newline
        if self.num_reads is None:
            num_reads = self.count_reads()
        else:
            num_reads = self.num_reads
        return bases_per_read * num_reads * 2
