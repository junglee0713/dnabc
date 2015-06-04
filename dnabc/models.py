import os.path

from .util import parse_barcode_file

class Sample(object):
    """Class representing one demultiplexable unit."""
    def __init__(self, name, barcode):
        self.name = name
        self.barcode = barcode
        if self.barcode is not None:
            self.barcode = self.barcode.upper()

    @classmethod
    def load(cls, f):
        return [cls(name, bc) for name, bc in parse_barcode_file(f)]
