"""``fqtk trim``: fixed-length trimming for single-end or paired-end FASTQ.

Each read is truncated to a target length; a read already shorter than the
target is kept in full (never padded, never dropped). In paired-end mode R1
and R2 can each have their own target length.
"""

import sys
from pathlib import Path

from .io import open_input, open_output, read_record, close_all


def trim_seq_qual(seq: str, qual: str, max_len: int):
    if len(seq) <= max_len:
        return seq, qual
    return seq[:max_len], qual[:max_len]


def add_arguments(parser):
    parser.add_argument("-1", "--r1", required=True,
                         help="Input R1 fastq(.gz) (the only input in single-end mode)")
    parser.add_argument("-2", "--r2", default=None,
                         help="Input R2 fastq(.gz); providing this switches to paired-end mode")
    parser.add_argument("-o1", "--out-r1", required=True, help="Output R1 fastq.gz")
    parser.add_argument("-o2", "--out-r2", default=None,
                         help="Output R2 fastq.gz (required in paired-end mode)")
    parser.add_argument("--r1-len", type=int, required=True,
                         help="Bases to keep from R1 (reads shorter than this are kept whole)")
    parser.add_argument("--r2-len", type=int, default=None,
                         help="Bases to keep from R2 (defaults to --r1-len in paired-end mode)")
    parser.add_argument("--min-len", type=int, default=1,
                         help="Drop a read/pair if its trimmed length falls below this (default: 1)")
    parser.add_argument("--threads", type=int, default=1,
                         help="pigz threads for gzip I/O (requires pigz on PATH; default 1 = builtin gzip)")
    parser.add_argument("--report-every", type=int, default=5_000_000,
                         help="Print progress every N reads, 0 to disable (default: 5,000,000)")
    parser.set_defaults(func=cmd_trim)


def cmd_trim(args):
    paired = args.r2 is not None

    if paired:
        if not args.out_r2:
            sys.exit("error: paired-end mode (-2/--r2 given) requires -o2/--out-r2")
        if args.r2_len is None:
            args.r2_len = args.r1_len
    elif args.r2_len is not None or args.out_r2 is not None:
        print("note: no -2/--r2 given, running in single-end mode; "
              "--r2-len / -o2 are ignored", file=sys.stderr)

    for p in [args.r1] + ([args.r2] if paired else []):
        if not Path(p).exists():
            sys.exit(f"input file not found: {p}")

    in1, p_in1 = open_input(args.r1, args.threads)
    out1, p_out1 = open_output(args.out_r1, args.threads)
    handles = [(in1, p_in1), (out1, p_out1)]

    in2 = out2 = None
    if paired:
        in2, p_in2 = open_input(args.r2, args.threads)
        out2, p_out2 = open_output(args.out_r2, args.threads)
        handles += [(in2, p_in2), (out2, p_out2)]

    n_total = n_kept = n_dropped_short = 0

    try:
        while True:
            rec1 = read_record(in1)
            rec2 = read_record(in2) if paired else None

            if rec1 is None and (not paired or rec2 is None):
                break
            if paired and (rec1 is None) != (rec2 is None):
                sys.exit("error: R1 / R2 read counts differ - files may be corrupt or mismatched")

            n_total += 1
            h1, s1, pl1, q1 = rec1
            s1_t, q1_t = trim_seq_qual(s1, q1, args.r1_len)

            if paired:
                h2, s2, pl2, q2 = rec2
                s2_t, q2_t = trim_seq_qual(s2, q2, args.r2_len)
                too_short = len(s1_t) < args.min_len or len(s2_t) < args.min_len
            else:
                too_short = len(s1_t) < args.min_len

            if too_short:
                n_dropped_short += 1
                continue

            out1.write(f"{h1}\n{s1_t}\n{pl1}\n{q1_t}\n")
            if paired:
                out2.write(f"{h2}\n{s2_t}\n{pl2}\n{q2_t}\n")
            n_kept += 1

            if args.report_every and n_total % args.report_every == 0:
                print(f"processed {n_total:,} reads ...", file=sys.stderr)
    finally:
        close_all(handles)

    mode = "paired-end (PE)" if paired else "single-end (SE)"
    print(f"[trim] mode: {mode}", file=sys.stderr)
    print(f"[trim] total: {n_total:,}, kept: {n_kept:,}", file=sys.stderr)
    if n_dropped_short:
        print(f"[trim] dropped (below --min-len after trimming): {n_dropped_short:,}", file=sys.stderr)
