"""``fqtk trim``: fixed-length trimming for single-end or paired-end FASTQ.

Each read is truncated to a target length; a read already shorter than the
target is kept in full (never padded, never dropped). In paired-end mode R1
and R2 can each have their own target length.

Two ways to run it:

- single sample: ``-1/-2/-o1/-o2/--r1-len/--r2-len`` as before.
- ``--batch manifest.tsv --jobs N``: many independent samples processed
  concurrently (see _cmd_trim_batch below for the manifest format). This is
  the one that actually uses more than one core in a meaningful way - see
  the note on --threads for why that flag alone usually doesn't.
"""

import csv
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from .io import open_input, open_output, read_record, close_all


class PairMismatchError(RuntimeError):
    """R1 and R2 have different numbers of reads."""


def trim_seq_qual(seq: str, qual: str, max_len: int):
    if len(seq) <= max_len:
        return seq, qual
    return seq[:max_len], qual[:max_len]


def run_trim_pair(r1, r2, out_r1, out_r2, r1_len, r2_len=None, min_len=1,
                   threads=1, report_every=0):
    """Trim one sample (single-end if r2 is None, paired-end otherwise).

    Pure function: no argparse, no sys.exit - raises FileNotFoundError /
    PairMismatchError on problems, otherwise returns a stats dict. This is
    what both the single-sample CLI path and each --batch worker call.
    """
    paired = r2 is not None
    if paired and r2_len is None:
        r2_len = r1_len

    for p in [r1] + ([r2] if paired else []):
        if not Path(p).exists():
            raise FileNotFoundError(f"input file not found: {p}")

    in1, p_in1 = open_input(r1, threads)
    out1, p_out1 = open_output(out_r1, threads)
    handles = [(in1, p_in1), (out1, p_out1)]

    in2 = out2 = None
    if paired:
        in2, p_in2 = open_input(r2, threads)
        out2, p_out2 = open_output(out_r2, threads)
        handles += [(in2, p_in2), (out2, p_out2)]

    n_total = n_kept = n_dropped_short = 0
    try:
        while True:
            rec1 = read_record(in1)
            rec2 = read_record(in2) if paired else None

            if rec1 is None and (not paired or rec2 is None):
                break
            if paired and (rec1 is None) != (rec2 is None):
                raise PairMismatchError("R1 / R2 read counts differ - files may be corrupt or mismatched")

            n_total += 1
            h1, s1, pl1, q1 = rec1
            s1_t, q1_t = trim_seq_qual(s1, q1, r1_len)

            if paired:
                h2, s2, pl2, q2 = rec2
                s2_t, q2_t = trim_seq_qual(s2, q2, r2_len)
                too_short = len(s1_t) < min_len or len(s2_t) < min_len
            else:
                too_short = len(s1_t) < min_len

            if too_short:
                n_dropped_short += 1
                continue

            out1.write(f"{h1}\n{s1_t}\n{pl1}\n{q1_t}\n")
            if paired:
                out2.write(f"{h2}\n{s2_t}\n{pl2}\n{q2_t}\n")
            n_kept += 1

            if report_every and n_total % report_every == 0:
                print(f"processed {n_total:,} reads ...", file=sys.stderr)
    finally:
        close_all(handles)

    return dict(mode="PE" if paired else "SE", total=n_total, kept=n_kept,
                dropped_short=n_dropped_short)


def add_arguments(parser):
    parser.add_argument("-1", "--r1", default=None,
                         help="Input R1 fastq(.gz) (the only input in single-end mode)")
    parser.add_argument("-2", "--r2", default=None,
                         help="Input R2 fastq(.gz); providing this switches to paired-end mode")
    parser.add_argument("-o1", "--out-r1", default=None, help="Output R1 fastq.gz")
    parser.add_argument("-o2", "--out-r2", default=None,
                         help="Output R2 fastq.gz (required in paired-end mode)")
    parser.add_argument("--r1-len", type=int, default=None,
                         help="Bases to keep from R1 (reads shorter than this are kept whole)")
    parser.add_argument("--r2-len", type=int, default=None,
                         help="Bases to keep from R2 (defaults to --r1-len in paired-end mode)")
    parser.add_argument("--min-len", type=int, default=1,
                         help="Drop a read/pair if its trimmed length falls below this (default: 1)")
    parser.add_argument("--threads", type=int, default=1,
                         help="pigz threads per sample for gzip I/O (requires pigz on PATH; "
                              "default 1 = builtin gzip). Note: pigz parallelizes *writing* "
                              "compressed output, but decompressing a normal, single-stream "
                              ".gz file is inherently sequential - this flag mostly speeds up "
                              "the output side, not the input side.")
    parser.add_argument("--report-every", type=int, default=5_000_000,
                         help="Print progress every N reads, 0 to disable (default: 5,000,000; "
                              "ignored in --batch mode)")
    parser.add_argument("--batch", default=None, metavar="MANIFEST.tsv",
                         help="Process many independent samples concurrently instead of one "
                              "pair. Tab-separated file with columns r1, out_r1, r1_len "
                              "(required) and r2, out_r2, r2_len, min_len, sample (optional, "
                              "leave the cell empty when not needed). Mutually exclusive with "
                              "-1/-2/-o1/-o2/--r1-len/--r2-len.")
    parser.add_argument("--jobs", type=int, default=1,
                         help="Samples to process concurrently in --batch mode (default: 1). "
                              "This is what actually uses multiple cores - total cores used "
                              "is roughly --jobs x --threads, keep that within what you were "
                              "allocated.")
    parser.set_defaults(func=cmd_trim)


def cmd_trim(args):
    if args.batch:
        _cmd_trim_batch(args)
        return

    if not (args.r1 and args.out_r1 and args.r1_len is not None):
        sys.exit("error: -1/--r1, -o1/--out-r1 and --r1-len are required outside of --batch mode")

    paired = args.r2 is not None
    if paired and not args.out_r2:
        sys.exit("error: paired-end mode (-2/--r2 given) requires -o2/--out-r2")
    if not paired and (args.r2_len is not None or args.out_r2 is not None):
        print("note: no -2/--r2 given, running in single-end mode; "
              "--r2-len / -o2 are ignored", file=sys.stderr)

    try:
        stats = run_trim_pair(args.r1, args.r2, args.out_r1, args.out_r2,
                               args.r1_len, args.r2_len, args.min_len,
                               args.threads, args.report_every)
    except (FileNotFoundError, PairMismatchError) as e:
        sys.exit(f"error: {e}")

    mode = "paired-end (PE)" if stats["mode"] == "PE" else "single-end (SE)"
    print(f"[trim] mode: {mode}", file=sys.stderr)
    print(f"[trim] total: {stats['total']:,}, kept: {stats['kept']:,}", file=sys.stderr)
    if stats["dropped_short"]:
        print(f"[trim] dropped (below --min-len after trimming): {stats['dropped_short']:,}", file=sys.stderr)


def _read_batch_manifest(path):
    with open(path, newline="") as f:
        rows = [dict(row) for row in csv.DictReader(f, delimiter="\t")]
    if not rows:
        sys.exit(f"error: no rows found in manifest {path}")
    for i, row in enumerate(rows, 1):
        if not row.get("r1") or not row.get("out_r1") or not row.get("r1_len"):
            sys.exit(f"error: manifest row {i} is missing a required column (r1, out_r1, r1_len)")
    return rows


def _sample_label(row):
    return row.get("sample") or os.path.basename(row["out_r1"])


def _batch_trim_worker(row, threads, default_min_len):
    sample = _sample_label(row)
    r2 = row.get("r2") or None
    out_r2 = row.get("out_r2") or None
    r2_len = int(row["r2_len"]) if row.get("r2_len") else None
    min_len = int(row["min_len"]) if row.get("min_len") else default_min_len
    try:
        stats = run_trim_pair(row["r1"], r2, row["out_r1"], out_r2,
                               int(row["r1_len"]), r2_len, min_len, threads, report_every=0)
        return sample, "ok", stats
    except Exception as e:  # a bad sample shouldn't take the whole batch down
        return sample, "error", str(e)


def _cmd_trim_batch(args):
    if args.r1 or args.r2 or args.out_r1 or args.out_r2 or args.r1_len is not None:
        sys.exit("error: --batch cannot be combined with -1/-2/-o1/-o2/--r1-len/--r2-len")

    rows = _read_batch_manifest(args.batch)
    order = {_sample_label(row): i for i, row in enumerate(rows)}

    results = []
    with ProcessPoolExecutor(max_workers=max(1, args.jobs)) as ex:
        futures = [ex.submit(_batch_trim_worker, row, args.threads, args.min_len) for row in rows]
        for fut in as_completed(futures):
            results.append(fut.result())
    results.sort(key=lambda r: order.get(r[0], 0))

    print("sample\tstatus\tmode\ttotal\tkept\tdropped_short")
    n_errors = 0
    for sample, status, payload in results:
        if status == "ok":
            print(f"{sample}\tok\t{payload['mode']}\t{payload['total']}\t"
                  f"{payload['kept']}\t{payload['dropped_short']}")
        else:
            print(f"{sample}\terror\t-\t-\t-\t-")
            print(f"  -> {payload}", file=sys.stderr)
            n_errors += 1

    if n_errors:
        sys.exit(f"\n{n_errors} of {len(rows)} sample(s) failed - see messages above")
