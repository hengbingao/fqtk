"""fqtk command-line entry point."""

import argparse

from . import __version__
from .trim import add_arguments as add_trim_arguments
from .readstats import add_arguments as add_stat_arguments
from .merge import add_arguments as add_merge_arguments


def build_parser():
    ap = argparse.ArgumentParser(
        prog="fqtk",
        description="Length-based FASTQ trimming, read-length statistics, and "
                     "same-library lane/batch merging (single-end / paired-end).",
    )
    ap.add_argument("--version", action="version", version=f"fqtk {__version__}")
    sub = ap.add_subparsers(dest="command", required=True)

    pt = sub.add_parser("trim", help="Trim reads to a fixed length (SE or PE)")
    add_trim_arguments(pt)

    ps = sub.add_parser("stat", help="Summarize read-length distribution (SE or PE)")
    add_stat_arguments(ps)

    pm = sub.add_parser("merge", help="Merge same-library FASTQ files across sequencing lanes/batches")
    add_merge_arguments(pm)

    return ap


def main():
    ap = build_parser()
    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
