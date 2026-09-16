"""``fqtk stat``: read-length distribution summary for single-end or paired-end FASTQ.

Two ways to run it:

- single sample: ``-1`` (and optionally ``-2``) as before.
- ``--batch manifest.tsv --jobs N``: summarize many independent samples
  concurrently, printed as one combined table.
"""

import csv
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from .io import open_input, read_record


def length_histogram(path, threads):
    fh, p = open_input(path, threads)
    hist = {}
    n = 0
    total_len = 0
    try:
        while True:
            rec = read_record(fh)
            if rec is None:
                break
            length = len(rec[1])
            hist[length] = hist.get(length, 0) + 1
            n += 1
            total_len += length
    finally:
        fh.close()
        if p is not None:
            p.wait()
    return n, total_len, hist


def summarize_hist(n, total_len, hist):
    if n == 0:
        return dict(n=0, min_len=0, max_len=0, mean_len=0.0, median_len=0)
    lengths_sorted = sorted(hist.keys())
    min_len = lengths_sorted[0]
    max_len = lengths_sorted[-1]
    mean_len = total_len / n
    half = (n + 1) // 2
    cum = 0
    median_len = lengths_sorted[-1]
    for length in lengths_sorted:
        cum += hist[length]
        if cum >= half:
            median_len = length
            break
    return dict(n=n, min_len=min_len, max_len=max_len, mean_len=mean_len, median_len=median_len)


def write_dist(hist, out_path):
    with open(out_path, "w") as f:
        f.write("length\tcount\n")
        for length in sorted(hist.keys()):
            f.write(f"{length}\t{hist[length]}\n")


def stat_one_file(path, threads):
    if not Path(path).exists():
        raise FileNotFoundError(f"input file not found: {path}")
    n, total_len, hist = length_histogram(path, threads)
    return summarize_hist(n, total_len, hist), hist


def add_arguments(parser):
    parser.add_argument("-1", "--r1", default=None,
                         help="Input R1 fastq(.gz) (the only input in single-end mode)")
    parser.add_argument("-2", "--r2", default=None,
                         help="Input R2 fastq(.gz); providing this switches to paired-end mode")
    parser.add_argument("-o", "--out-prefix", default=None,
                         help="If set, write full length distribution to "
                              "<prefix>_R1.length_dist.tsv (and _R2 in paired-end mode)")
    parser.add_argument("--threads", type=int, default=1,
                         help="pigz threads for decompression (requires pigz on PATH; "
                              "default 1 = builtin gzip). Note: pigz mainly speeds up "
                              "*compression*; decompressing a normal, single-stream .gz file "
                              "is inherently sequential, so this flag gives little to no "
                              "speedup here since stat only ever reads.")
    parser.add_argument("--batch", default=None, metavar="MANIFEST.tsv",
                         help="Summarize many independent samples concurrently. Tab-separated "
                              "file with columns: r1 (required), r2, out_prefix, sample "
                              "(optional). Mutually exclusive with -1/-2/-o.")
    parser.add_argument("--jobs", type=int, default=1,
                         help="Samples to process concurrently in --batch mode (default: 1) - "
                              "this is what actually uses multiple cores for `stat`.")
    parser.set_defaults(func=cmd_stat)


def cmd_stat(args):
    if args.batch:
        _cmd_stat_batch(args)
        return

    if not args.r1:
        sys.exit("error: -1/--r1 is required outside of --batch mode")

    paired = args.r2 is not None
    files = [("SE" if not paired else "R1", args.r1)]
    if paired:
        files.append(("R2", args.r2))

    for _, p in files:
        if not Path(p).exists():
            sys.exit(f"input file not found: {p}")

    print("file\tread\treads\tmin_len\tmax_len\tmean_len\tmedian_len")
    for tag, path in files:
        n, total_len, hist = length_histogram(path, args.threads)
        s = summarize_hist(n, total_len, hist)
        print(f"{path}\t{tag}\t{s['n']}\t{s['min_len']}\t{s['max_len']}\t"
              f"{s['mean_len']:.2f}\t{s['median_len']}")
        if args.out_prefix:
            dist_path = f"{args.out_prefix}_{tag}.length_dist.tsv"
            write_dist(hist, dist_path)
            print(f"# length distribution written to: {dist_path}", file=sys.stderr)


def _read_batch_manifest(path):
    with open(path, newline="") as f:
        rows = [dict(row) for row in csv.DictReader(f, delimiter="\t")]
    if not rows:
        sys.exit(f"error: no rows found in manifest {path}")
    for i, row in enumerate(rows, 1):
        if not row.get("r1"):
            sys.exit(f"error: manifest row {i} is missing the required 'r1' column")
    return rows


def _sample_label(row):
    return row.get("sample") or os.path.basename(row["r1"])


def _batch_stat_worker(row, threads):
    sample = _sample_label(row)
    paired = bool(row.get("r2"))
    files = [("SE" if not paired else "R1", row["r1"])]
    if paired:
        files.append(("R2", row["r2"]))
    try:
        out = []
        for tag, path in files:
            s, hist = stat_one_file(path, threads)
            if row.get("out_prefix"):
                write_dist(hist, f"{row['out_prefix']}_{tag}.length_dist.tsv")
            out.append((tag, path, s))
        return sample, "ok", out
    except Exception as e:
        return sample, "error", str(e)


def _cmd_stat_batch(args):
    if args.r1 or args.r2 or args.out_prefix:
        sys.exit("error: --batch cannot be combined with -1/-2/-o")

    rows = _read_batch_manifest(args.batch)
    order = {_sample_label(row): i for i, row in enumerate(rows)}

    results = []
    with ProcessPoolExecutor(max_workers=max(1, args.jobs)) as ex:
        futures = [ex.submit(_batch_stat_worker, row, args.threads) for row in rows]
        for fut in as_completed(futures):
            results.append(fut.result())
    results.sort(key=lambda r: order.get(r[0], 0))

    print("sample\tfile\tread\treads\tmin_len\tmax_len\tmean_len\tmedian_len")
    n_errors = 0
    for sample, status, payload in results:
        if status == "ok":
            for tag, path, s in payload:
                print(f"{sample}\t{path}\t{tag}\t{s['n']}\t{s['min_len']}\t{s['max_len']}\t"
                      f"{s['mean_len']:.2f}\t{s['median_len']}")
        else:
            print(f"{sample}\t-\terror\t-\t-\t-\t-\t-")
            print(f"  -> {payload}", file=sys.stderr)
            n_errors += 1

    if n_errors:
        sys.exit(f"\n{n_errors} of {len(rows)} sample(s) failed - see messages above")
