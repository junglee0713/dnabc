import os.path

from .util import parse_barcode_file, duplicates

class Sample(object):
    """Class representing one demultiplexable unit."""
    def __init__(self, name, barcode):
        self.name = name
        self.barcode = barcode
        if self.barcode is not None:
            self.barcode = self.barcode.upper()

    @classmethod
    def load(cls, f):
        records = list(parse_barcode_file(f))
        names, bcs = zip(*records)

        dup_names = duplicates(names)
        if dup_names:
            raise ValueError("Duplicate sample names: %s" % dup_names)

        dup_bcs = duplicates(bcs)
        if dup_bcs:
            raise ValueError("Duplicate barcodes: %s" % dup_bcs)

        return [cls(name, bc) for name, bc in records]
