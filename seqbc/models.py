import os.path

from .util import deambiguate, parse_barcode_file


class Sample(object):
    """Class representing one demultiplexable unit
    """
    def __init__(self, name, barcode):
        self.name = name

        self.barcode = barcode
        if self.barcode is not None:
            self.barcode = self.barcode.upper()

    @classmethod
    def load(cls, f):
        return [cls(name, bc) for name, bc in parse_barcode_file(f)]


class Run(object):
    """Class representing one pool of multiplexed samples
    """
    def __init__(self, accession, date, mach_type, mach_kit, region, fp, comment):
        self.accession = accession  # Must be a number
        self.date = date
        self.machine_type = mach_type
        self.machine_kit = mach_kit
        self.region = region
        self.fp = fp
        self.comment = comment

    @property
    def file_format(self):
        """Infer the input file format.

        Unlike the methods elsewhere, this does not check that the file
        exists.
        """
        fp = self.fp
        if fp.endswith(".gz"):
            fp = fp[:-3]
        _, ext = os.path.splitext(fp)
        ext = ext.upper()[1:] # Remove leading dot and make ALL CAPS

        # The folder name is registered for Illumina runs where
        # samples are pre-split.  There should be no other type of run
        # without a file extension.  Therefore, if we see a blank file
        # extension, we assume that the filepath is a folder
        # containing FASTQ files.
        if ext == "":
            return "FASTQ"
        # FASTA files can also be named with .fna
        if ext == "FNA":
            return "FASTA"
        return ext
    
    @property
    def formatted_accession(self):
        return "BLR%06d" % self.accession

    @property
    def platform(self):
        return "%s %s" % (self.machine_type, self.machine_kit)

    @property
    def prefix(self):
        if self.fp.endswith(".sff"):
            fn = os.path.basename(self.fp)
            return fn[:-6] # runprefix01.sff, remove 6 chars from end
        return ""


