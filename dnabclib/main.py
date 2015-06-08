import argparse
import json
import os

from .writer import FastaWriter, PairedFastqWriter
from .models import Sample
from .run_file import IndexFastqSequenceFile
from .assign import BarcodeAssigner
from .version import __version__

writers = {
    "fastq": PairedFastqWriter,
    "fasta": FastaWriter,
    }

default_config = {
    "output_format": "fastq"
    }

def main(argv=None):
    p = argparse.ArgumentParser()
    # Input
    p.add_argument("--forward-reads", required=True,
        type=argparse.FileType("r"),
        help="Forward reads file (Gzipped FASTQ format)")
    p.add_argument("--reverse-reads", required=True,
        type=argparse.FileType("r"),
        help="Reverse reads file (Gzipped FASTQ format)")
    p.add_argument("--index-reads", required=True,
        type=argparse.FileType("r"),
        help="Index reads file (Gzipped FASTQ format)")
    p.add_argument(
        "--barcode-file", required=True,
        help="Barcode information file",
        type=argparse.FileType("r"))
    # Output
    p.add_argument(
        "--output-dir", required=True,
        help="Output sequence data directory")
    p.add_argument(
        "--summary-file", required=True,
        type=argparse.FileType("w"),
        help="Summary filepath")
    # Config
    p.add_argument("--config-file",
        type=argparse.FileType("r"),
        help="Configuration file (JSON format)")
    args = p.parse_args(argv)

    config = default_config
    if args.config_file:
        user_config = json.load(args.config_file)
        config.update(user_config)
    writer_cls = writers[config["output_format"]]

    # Close input sequence files before processing
    fwd_fp = args.forward_reads.name
    rev_fp = args.reverse_reads.name
    idx_fp = args.index_reads.name
    args.forward_reads.close()
    args.reverse_reads.close()
    args.index_reads.close()

    seq_file = IndexFastqSequenceFile(fwd_fp, rev_fp, idx_fp)

    samples = list(Sample.load(args.barcode_file))
    assigner = BarcodeAssigner(samples)

    if os.path.exists(args.output_dir):
        p.error("Output directory already exists")
    os.mkdir(args.output_dir)
    writer = writer_cls(args.output_dir)

    summary_data = seq_file.demultiplex(assigner, writer)

    save_summary(args.summary_file, config, summary_data)


def save_summary(f, config, data):
    result = {
        "program": "dnabc",
        "version": __version__,
        "config": config,
        "data": data,
        }
    json.dump(result, f)
