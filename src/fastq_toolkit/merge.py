"""``fqtk merge``: merge same-library FASTQ files split across sequencing
lanes/batches, and consolidate everything else into one output folder.

Grouping logic (unchanged from the original ``merge_fastq_by_RL.py``):

  1. The RL (library) id is the filename token before the first ``_``
     (must match ``RL\\d+``).
  2. R1 / R2 is read off whichever appears as a standalone ``_``-delimited
     token in the filename.
  3. Files sharing the same (RL id, read) are the same library/end,
     sequenced in different lanes or batches -> they get merged.
  4. Within a group, filenames are split on ``_`` and compared position by
     position: a position where every file agrees is kept (it's the part
     that's common to the library), a position where they differ is
     dropped (it's the batch-specific part - date, lane, sample index,
     ...). The surviving tokens, rejoined with ``_``, become the merged
     file's name.
  5. Files are concatenated as raw bytes (``cat a.fastq.gz b.fastq.gz >
     merged.fastq.gz``) - multi-member gzip is valid and every FASTQ tool
     (zcat, aligners, etc.) reads it as one continuous stream, which is
     the standard way to merge multi-lane FASTQ.

A group with only one file needs no merging: by default it's *moved*
as-is into the output directory (see --copy-singletons to copy instead).
A file whose name doesn't parse to (RL id, read) is left untouched in the
input directory and reported as skipped.

Every run (preview or --execute) prints a plan; every --execute run also
writes a ``merge_report.tsv`` you can inspect afterwards.
"""

import glob
import os
import re
import shutil
import sys
from collections import defaultdict
from pathlib import Path


def find_rl_id(filename):
    """RL id = the token before the first '_', must look like RL<digits>."""
    first_token = filename.split("_", 1)[0]
    if re.fullmatch(r"RL\d+", first_token):
        return first_token
    return None


def find_read(filename):
    """R1/R2, read off a standalone '_'-delimited token."""
    tokens = filename.split("_")
    if "R1" in tokens:
        return "R1"
    if "R2" in tokens:
        return "R2"
    return None


def common_tokens(filenames):
    """Token-wise intersection across a group's filenames, order preserved."""
    token_lists = [f.split("_") for f in filenames]
    min_len = min(len(t) for t in token_lists)
    result = []
    for i in range(min_len):
        values = {t[i] for t in token_lists}
        if len(values) == 1:
            result.append(token_lists[0][i])
    return result


def build_plan(input_dir, pattern):
    files = sorted(glob.glob(os.path.join(input_dir, pattern)))
    if not files:
        sys.exit(f"error: no files matching {pattern!r} found in {input_dir}")

    groups = defaultdict(list)  # (rl, read) -> [basename, ...]
    skipped = []

    for path in files:
        base = os.path.basename(path)
        rl = find_rl_id(base)
        read = find_read(base)
        if rl is None or read is None:
            skipped.append(base)
            continue
        groups[(rl, read)].append(base)

    plan = []
    for (rl, read), basenames in sorted(groups.items()):
        basenames = sorted(basenames)
        if len(basenames) > 1:
            common = common_tokens(basenames)
            if not common:
                skipped.extend(basenames)
                continue
            plan.append(dict(rl=rl, read=read, action="merge",
                              out_name="_".join(common), sources=basenames))
        else:
            plan.append(dict(rl=rl, read=read, action="singleton",
                              out_name=basenames[0], sources=basenames))

    return plan, skipped


def print_plan(plan, skipped):
    print("=" * 70)
    print(f"found {len(plan)} group(s) (same RL id + same read)")
    print("=" * 70)
    for row in plan:
        tag = "MERGE" if row["action"] == "merge" else "KEEP (single batch, relocated as-is)"
        print(f"\n[{tag}] -> {row['out_name']}")
        for s in row["sources"]:
            print(f"    from: {s}")
    if skipped:
        print("\n[skipped] filename didn't parse to (RL id, R1/R2), or group had no common part:")
        for s in skipped:
            print(f"    {s}")


def write_report(report_path, plan, skipped):
    with open(report_path, "w") as f:
        f.write("rl\tread\taction\toutput_file\tn_sources\tsource_files\n")
        for row in plan:
            f.write(f"{row['rl']}\t{row['read']}\t{row['action']}\t{row['out_name']}\t"
                     f"{len(row['sources'])}\t{','.join(row['sources'])}\n")
        for s in skipped:
            f.write(f"NA\tNA\tskipped\t{s}\t0\t{s}\n")


def add_arguments(parser):
    parser.add_argument("input_dir", help="Directory containing the raw fastq(.gz) files")
    parser.add_argument("-o", "--output-dir", default="merged",
                         help="Destination for merged + relocated files (default: ./merged)")
    parser.add_argument("--pattern", default="*.fastq.gz",
                         help="Glob pattern for input files (default: *.fastq.gz)")
    parser.add_argument("--execute", action="store_true",
                         help="Actually merge/move files; without this, only print the plan")
    parser.add_argument("--copy-singletons", action="store_true",
                         help="Copy (keep the original) instead of moving files that need no merging")
    parser.add_argument("--remove-merged-sources", action="store_true",
                         help="After a successful merge, delete the per-lane source files "
                              "(off by default - sources are kept until you've verified the merge)")
    parser.add_argument("--force", action="store_true",
                         help="Overwrite a file already present at the destination instead of skipping it")
    parser.add_argument("--report", default=None,
                         help="Report TSV path (default: <output-dir>/merge_report.tsv)")
    parser.set_defaults(func=cmd_merge)


def cmd_merge(args):
    if not Path(args.input_dir).is_dir():
        sys.exit(f"error: not a directory: {args.input_dir}")

    plan, skipped = build_plan(args.input_dir, args.pattern)
    print_plan(plan, skipped)

    if not args.execute:
        print("\n(preview only - nothing was written. re-run with --execute to actually merge/move)")
        return

    os.makedirs(args.output_dir, exist_ok=True)
    report_path = args.report or os.path.join(args.output_dir, "merge_report.tsv")

    for row in plan:
        out_path = os.path.join(args.output_dir, row["out_name"])
        in_paths = [os.path.join(args.input_dir, b) for b in row["sources"]]

        if os.path.exists(out_path) and not args.force:
            print(f"[skip] destination already exists, not overwriting (use --force): {out_path}",
                  file=sys.stderr)
            row["action"] += "_skipped_exists"
            continue

        if row["action"] == "singleton":
            src = in_paths[0]
            if args.copy_singletons:
                shutil.copy2(src, out_path)
            else:
                shutil.move(src, out_path)
        else:  # merge
            expected_size = sum(os.path.getsize(p) for p in in_paths)
            with open(out_path, "wb") as out_f:
                for p in in_paths:
                    with open(p, "rb") as in_f:
                        shutil.copyfileobj(in_f, out_f)
            actual_size = os.path.getsize(out_path)
            if actual_size != expected_size:
                print(f"[warning] size mismatch after merging into {out_path}: "
                      f"expected {expected_size} bytes, got {actual_size} "
                      f"- check disk space / I/O errors before trusting this file",
                      file=sys.stderr)
            if args.remove_merged_sources:
                for p in in_paths:
                    os.remove(p)

        print(f"done: {out_path}")

    write_report(report_path, plan, skipped)
    print(f"\nreport written to: {report_path}")
    print("done.")
