import collections
import itertools

from .util import (
    AMBIGUOUS_BASES_COMPLEMENT, deambiguate, reverse_complement,
    )


class _Assigner(object):
     def has_reads(self, sample):
        return self.read_counts[sample.name] > 0

    
class BarcodeAssigner(_Assigner):
    def __init__(self, samples, mismatches=1, revcomp=True):
        self.samples = samples
        if mismatches not in [0, 1, 2]:
            raise ValueError(
                "Only 0, 1, or 2 mismatches allowed (got %s)" % mismatches)
        self.mismatches = mismatches
        self.revcomp = revcomp
        self.read_counts = dict((s.name, 0) for s in self.samples)
        self._init_hash()

    def _init_hash(self):
        self._barcodes = {}
        for s in self.samples:
            if not s.barcode:
                raise ValueError(
                    "All samples must be barcoded to use BarcodeAssigner: "
                    "%s" % s.name)

            if self.revcomp:
                bc = reverse_complement(s.barcode)
            else:
                bc = s.barcode

            if bc in self._barcodes:
                raise ValueError(
                    "Barcode for sample %s matches barcode for sample %s" % (
                        s, self._barcodes[bc]))
            else:
                self._barcodes[bc] = s

            for error_bc in self._error_barcodes(bc):
                if error_bc in self._barcodes:
                    raise ValueError(
                        "Barcode %s for sample %s matches barcode for "
                        "sample %s with %s mismatches" % (
                            error_bc, self._barcodes[error_bc], s,
                            self.mismatches))
                else:
                    self._barcodes[error_bc] = s

    def _error_barcodes(self, barcode):
        # If the number of mismatches is set to 0, there will be no
        # error barcodes. Immediately stop the iteration.
        if self.mismatches == 0:
            raise StopIteration
        # Each item in idx_sets is a set of indices where mismatches
        # should occur.
        idx_sets = itertools.combinations(range(len(barcode)), self.mismatches)
        for idx_set in idx_sets:
            # Change to list because strings are immutable
            bc = list(barcode)
            # Replace the base at each mismatch position with an
            # ambiguous base specifying all possibilities BUT the one
            # we see.
            for idx in idx_set:
                bc[idx] = AMBIGUOUS_BASES_COMPLEMENT[bc[idx]]
            # Expand to all possibilities for mismatching at this
            # particular set of positions
            for error_bc in deambiguate(bc):
                yield error_bc
        
    def assign(self, read):
        sample = self._barcodes.get(read.seq)
        if sample is not None:
            self.read_counts[sample.name] += 1
        return sample
    
