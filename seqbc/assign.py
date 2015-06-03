import collections
import itertools

from pytrie import StringTrie

from .util import (
    AMBIGUOUS_BASES_COMPLEMENT, deambiguate, reverse_complement,
    )


class _Assigner(object):
     def has_reads(self, sample):
        return self.read_counts[sample.name] > 0


class PrefixAssigner(_Assigner):
    """Assigns reads based on longest prefix match (exact match required).
    """
    def __init__(self, samples):
        self.samples = samples
        self.read_counts = collections.defaultdict(int)
        self._init_trie()

    def _init_trie(self):
        self._trie = StringTrie()
        for sample in self.samples:
            for prefix in sample.prefixes:
                self._trie[prefix] = sample

    def assign(self, read):
        sample = self._trie.longest_prefix_value(read.seq, None)
        if sample is not None:
            self.read_counts[sample.name] += 1
        return sample

    
class BarcodeAssigner(_Assigner):
    def __init__(self, samples, mismatches=1, revcomp=True):
        self.samples = samples
        if mismatches not in [0, 1, 2]:
            raise ValueError(
                "Only 0, 1, or 2 mismatches allowed (got %s)" % mismatches)
        self.mismatches = mismatches
        self.revcomp = revcomp
        self.read_counts = collections.defaultdict(int)
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
            self._barcodes[bc] = s
            for error_bc in self._error_barcodes(bc):
                # If we collide with another barcode here, we have a
                # BIG problem -- can't tell the difference between two
                # samples at this number of mismatches.
                if error_bc in self._barcodes:
                    raise ValueError(
                        "Barcode %s for sample %s matches barcode for "
                        "sample %s with %s mismatches" % (
                            error_bc, self._barcodes[error_bc],
                            s.name, self.mismatches))
                self._barcodes[error_bc] = s

    def _error_barcodes(self, barcode):
        # Each item in idx_sets is a set of indices where mismatches
        # should occur
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
    
