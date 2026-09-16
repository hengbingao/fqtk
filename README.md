# fastq-toolkit

Small command-line tool for two common FASTQ preprocessing tasks, both
single-end (SE) and paired-end (PE) aware:

- **`fqtk trim`** — trim reads to a fixed length. A read shorter than the
  target length is kept in full (never padded, never dropped); a longer
  read is truncated from the 3' end. In PE mode R1 and R2 can each have
  their own target length.
- **`fqtk stat`** — summarize the read-length distribution of a FASTQ file
  (or an R1/R2 pair): read count, min/max/mean/median length, and
  optionally the full length → count table for plotting.
- **`fqtk merge`** — given a folder of raw FASTQ files, auto-detect files
  from the same library (same RL id) split across sequencing lanes/batches,
  concatenate each such group into one file, relocate any library that was
  only sequenced once as-is, and write a report of what happened to what.

No third-party Python dependencies — only the standard library. If
[`pigz`](https://zlib.net/pigz/) is installed on your system, pass
`--threads N` to use it for multi-threaded gzip I/O on large files;
otherwise everything falls back to Python's built-in `gzip` module.

## Install

```bash
git clone https://github.com/hengbingao/fqtk
cd fqtk
pip install .
```

For active development (edits take effect immediately, no reinstall needed):

```bash
pip install -e ".[dev]"
```

Either way this installs an `fqtk` command onto your `PATH` inside whatever
Python environment (conda env, venv, etc.) you ran `pip install` in — no
need to invoke a `.py` script directly.

```bash
fqtk --version
fqtk trim --help
fqtk stat --help
```

## Usage

### Trim — single-end

```bash
fqtk trim -1 R1.fastq.gz -o1 R1.trim67.fastq.gz --r1-len 67
```

### Trim — paired-end

`--r2-len` defaults to `--r1-len` if omitted.

```bash
fqtk trim \
  -1 R1.fastq.gz -2 R2.fastq.gz \
  -o1 R1.trim67.fastq.gz -o2 R2.trim67.fastq.gz \
  --r1-len 67 --r2-len 67 \
  --threads 8
```

Optional: `--min-len N` drops a read (or read pair, in PE mode) if its
trimmed length falls below `N` (default 1, i.e. nothing extra is dropped).

### Stat — single-end

```bash
fqtk stat -1 R1.fastq.gz
```

```
file           read  reads  min_len  max_len  mean_len  median_len
R1.fastq.gz    SE    3      28       90       61.33     67
```

### Stat — paired-end, with the full length distribution exported

```bash
fqtk stat -1 R1.fastq.gz -2 R2.fastq.gz -o sample_stats
```

Writes `sample_stats_R1.length_dist.tsv` and `sample_stats_R2.length_dist.tsv`
(`length\tcount` per line) alongside the summary table on stdout.

### Merge — same library, different lanes/batches

Files are grouped by (RL id, R1/R2). The RL id is the token before the
first `_` in the filename (e.g. `RL8518`); R1/R2 is read off whichever
appears as its own `_`-delimited token. Within a group, filenames are
compared token by token — a position where every file in the group agrees
is kept, a position where they differ (date, lane, sample index, ...) is
dropped — and the surviving tokens become the merged file's name.

Always preview first (no files are touched without `--execute`):

```bash
fqtk merge /path/to/raw_fastq_dir -o /path/to/merged_dir
```

```
[MERGE] -> RL8518_2026_BulkATAC_P5C10_BEP100_ZIM3_rep2_R1_001.fastq.gz
    from: RL8518_2026_07_28_BulkATAC_P5C10_BEP100_ZIM3_rep2_S11_R1_001.fastq.gz
    from: RL8518_2026_08_05_BulkATAC_P5C10_BEP100_ZIM3_rep2_S22_R1_001.fastq.gz

[KEEP (single batch, relocated as-is)] -> RL9001_2026_08_05_BulkCT_H3K4me3_S30_R1_001.fastq.gz
    from: RL9001_2026_08_05_BulkCT_H3K4me3_S30_R1_001.fastq.gz
```

Once the plan looks right, run it for real:

```bash
fqtk merge /path/to/raw_fastq_dir -o /path/to/merged_dir --execute
```

What happens on `--execute`:

- A group with **more than one file** is byte-concatenated into a new file
  in the output directory (`cat lane1.fastq.gz lane2.fastq.gz > merged.fastq.gz`
  — valid multi-member gzip, reads correctly in zcat/aligners/anything
  else). The original per-lane files are **kept** in the input directory
  by default; pass `--remove-merged-sources` to delete them once you've
  verified the merge.
- A group with **exactly one file** needs no merging, so it's **moved**
  (not copied) into the output directory as-is. Pass `--copy-singletons`
  to copy instead and leave the original in place.
- A file whose name doesn't parse to (RL id, R1/R2) is **left untouched**
  in the input directory and listed as skipped.
- A `merge_report.tsv` is written into the output directory (path
  `rl / read / action / output_file / n_sources / source_files`) — one row
  per merged group, per singleton, and per skipped file — so you always
  have a record of exactly what became what. Override its path with
  `--report`.
- `--force` overwrites a destination file that already exists (by default
  it's left alone and reported as skipped, so re-running `--execute` is
  safe).

## Running the tests

```bash
pip install -e ".[dev]"
pytest
```

## Project layout

```
fastq-toolkit/
├── pyproject.toml          # packaging + `fqtk` console-script entry point
├── src/fastq_toolkit/
│   ├── io.py                # shared gzip/pigz FASTQ read/write helpers
│   ├── trim.py               # `fqtk trim`
│   ├── readstats.py          # `fqtk stat`
│   ├── merge.py               # `fqtk merge`
│   └── cli.py                 # argparse wiring, `fqtk` entry point
└── tests/
    ├── test_basic.py          # SE/PE trim + stat sanity tests
    └── test_merge.py          # lane/batch merge sanity tests
```

## License

MIT — see [LICENSE](LICENSE).
