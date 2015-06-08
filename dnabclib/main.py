import argparse
import json
import os

from .writer import FastaWriter, PairedFastqWriter
from .models import Sample
from .run_file import SequenceFile
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
    p.add_argument(
        "--sequence-file", required=True,
        help="Input sequence data filepath",
        type=argparse.FileType("r"))
    p.add_argument(
        "--barcode-file", required=True,
        help="Barcode information filepath",
        type=argparse.FileType("r"))
    p.add_argument(
        "--output-dir", required=True,
        help="Output sequence data directory")
    p.add_argument(
        "--summary-file", required=True,
        type=argparse.FileType("w"),
        help="Summary filepath")
    p.add_argument("--config-file",
        type=argparse.FileType("r"),
        help="Configuration file (JSON format)")
    args = p.parse_args(argv)

    config = default_config
    if args.config_file:
        user_config = json.load(args.config_file)
        config.update(user_config)
    writer_cls = writers[config["output_format"]]

    seq_file = SequenceFile(args.sequence_file.name)
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
