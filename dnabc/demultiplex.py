import argparse
import json
import os

from .writer import (
    FastaWriter, PooledFastaWriter,
    FastqWriter, PooledFastqWriter,
    PairedFastqWriter,
    )
from .models import Sample
from .run_file import SequenceFile
from .assign import BarcodeAssigner
from .version import __version__

writers = {
    "fastq": PairedFastqWriter,
    "fasta": FastaWriter,
    "pooled_fastq": PooledFastqWriter,
    "pooled_fasta": PooledFastaWriter,
    }

def demultiplex(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument(
        "--sequence-file", "-s", required=True,
        help="Input sequence data filepath",
        type=argparse.FileType("r"))
    p.add_argument(
        "--barcode-file", "-b", required=True,
        help="Barcode information filepath",
        type=argparse.FileType("r"))
    p.add_argument(
        "--output-dir", "-o", required=True,
        help="Output sequence data directory")
    p.add_argument(
        "--output-format", "-f", required=True,
        help="Output format",
        choices=writers.keys())
    args = p.parse_args(argv)

    seq_file = SequenceFile(args.sequence_file.name)
    samples = list(Sample.load(args.barcode_file))
    assigner = BarcodeAssigner(samples)
    
    if os.path.exists(args.output_dir):
        p.error("Output directory already exists")

    os.mkdir(args.output_dir)
    writer_cls = writers[args.output_format]
    writer = writer_cls(samples, args.output_dir)

    read_counts = seq_file.demultiplex(assigner, writer)

    summary_fp = os.path.join(args.output_dir, "demultiplex_summary.json")
    save_summary(summary_fp, read_counts)


def save_summary(fp, data):
    result = {
        "program": "dnabc",
        "version": __version__,
        "data": data,
        }
    with open(fp, "w") as f:
        json.dump(result, f)
