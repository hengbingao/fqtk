"""``fqtk stat``: read-length distribution summary for single-end or paired-end FASTQ."""

import sys
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


def add_arguments(parser):
    parser.add_argument("-1", "--r1", required=True,
                         help="Input R1 fastq(.gz) (the only input in single-end mode)")
    parser.add_argument("-2", "--r2", default=None,
                         help="Input R2 fastq(.gz); providing this switches to paired-end mode")
    parser.add_argument("-o", "--out-prefix", default=None,
                         help="If set, write full length distribution to "
                              "<prefix>_R1.length_dist.tsv (and _R2 in paired-end mode)")
    parser.add_argument("--threads", type=int, default=1,
                         help="pigz threads for decompression (requires pigz on PATH; "
                              "default 1 = builtin gzip)")
    parser.set_defaults(func=cmd_stat)


def cmd_stat(args):
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
