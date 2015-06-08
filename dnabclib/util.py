import collections
from cStringIO import StringIO
import itertools
import os


def duplicates(xs):
    # From http://stackoverflow.com/questions/9835762/
    seen = set()
    seen_add = seen.add
    # adds all elements it doesn't know yet to seen and all other to seen_twice
    seen_twice = set(x for x in xs if x in seen or seen_add(x))
    # turn the set into a list (as requested)
    return list(seen_twice)


def parse_fasta(f):
    f = iter(f)
    desc = f.next().strip()[1:]
    seq = StringIO()
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            yield desc, seq.getvalue()
            desc = line[1:]
            seq = StringIO()
        else:
            seq.write(line)
    yield desc, seq.getvalue()


def parse_barcode_file(f):
    for n, line in enumerate(f):
        if line.startswith("#"):
            continue
        line = line.rstrip()
        if line == "":
            continue
        toks = line.split("\t")
        if len(toks) < 2:
            line_num = n + 1
            raise ValueError(
                "Not enough fields in barcode file (line %s): %s" % (
                    line_num, toks))
        sample_id = toks[0]
        barcode = toks[1]
        yield sample_id, barcode


class FastaRead(object):
    def __init__(self, read):
        self.desc, self.seq = read


AMBIGUOUS_BASES = {
    "T": "T",
    "C": "C",
    "A": "A",
    "G": "G",
    "R": "AG",
    "Y": "TC",
    "M": "CA",
    "K": "TG",
    "S": "CG",
    "W": "TA",
    "H": "TCA",
    "B": "TCG",
    "V": "CAG",
    "D": "TAG",
    "N": "TCAG",
    }

# Ambiguous base codes for all bases EXCEPT the key
AMBIGUOUS_BASES_COMPLEMENT = {
    "T": "V",
    "C": "D",
    "A": "B",
    "G": "H",
    }


def deambiguate(seq):
    nt_choices = [AMBIGUOUS_BASES[x] for x in seq]
    return ["".join(c) for c in itertools.product(*nt_choices)]

COMPLEMENT_BASES = {
    "T": "A",
    "C": "G",
    "A": "T",
    "G": "C",
    }

def reverse_complement(seq):
    rc = [COMPLEMENT_BASES[x] for x in seq]
    rc.reverse()
    return ''.join(rc)
