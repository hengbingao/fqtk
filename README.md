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

No third-party Python dependencies — only the standard library. See
[Speeding things up](#speeding-things-up) for what actually parallelizes
here and what doesn't.

## Install

```bash
git clone https://github.com/<your-username>/fastq-toolkit.git
cd fastq-toolkit
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
first `_` in the filename (e.g. `RL8518`); R1/R2 is read off the filename
wherever it appears bounded by `_` or `.` on both sides (or start/end) —
so both the original Illumina-style `..._R1_001.fastq.gz` and this
toolkit's own `fqtk trim` output `..._R1.trim26.out.fastq.gz` (where R1 is
glued to the next part with a `.` instead of `_`) are recognized. Within a
group, filenames are then split on `_` and compared token by token — a
position where every file agrees is kept, a position where they differ
(date, lane, sample index, ...) is dropped — and the surviving tokens
become the merged file's name.

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

## Speeding things up

Three different flags all sound like "parallelism" but do different things,
and only two of them actually put more cores to work:

- **`--threads N` on `trim`/`stat`** hands gzip I/O off to
  [`pigz`](https://zlib.net/pigz/) if it's on `PATH` (falls back to Python's
  built-in `gzip` module otherwise). This genuinely speeds up **writing**
  compressed output in parallel. It does **not** meaningfully speed up
  **reading** — decompressing a normal, single-stream `.gz` file (which is
  what `bcl2fastq`/`fastp`/Illumina output normally is) is inherently
  sequential regardless of how many threads you throw at it; that's a
  property of the gzip format, not a limitation of this tool. So `--threads`
  on `trim` helps somewhat (its output side benefits), and on `stat` it
  barely helps at all (`stat` only ever reads).

- **`--jobs N` on `trim --batch` / `stat --batch`** is where the real
  speedup for many-sample workflows comes from: it runs N *independent
  samples* concurrently in separate processes, each with its own
  `--threads`. This is the right level to parallelize at for something like
  a directory of 30 histone marks — each sample's read/write is
  independent, so this scales with your core count the way the per-lane
  loop in a hand-rolled SLURM script would, minus the bookkeeping.
  **Total cores used ≈ `--jobs` × `--threads`** — keep that at or under
  what you were allocated (`--cpus-per-task` in your SLURM script).

  Batch manifest format (tab-separated, `.tsv`):

  ```
  sample	r1	r2	out_r1	out_r2	r1_len	r2_len
  H3K4me3_rep1	H3K4me3_rep1_R1.fastq.gz	H3K4me3_rep1_R2.fastq.gz	H3K4me3_rep1_R1.trim.fastq.gz	H3K4me3_rep1_R2.trim.fastq.gz	67	67
  H3K27ac_rep1	H3K27ac_rep1_R1.fastq.gz	H3K27ac_rep1_R2.fastq.gz	H3K27ac_rep1_R1.trim.fastq.gz	H3K27ac_rep1_R2.trim.fastq.gz	67	67
  ```

  `r1`, `out_r1`, `r1_len` are required; `r2`, `out_r2`, `r2_len`, `min_len`,
  `sample` are optional (leave the cell empty for a single-end sample, or to
  fall back to the default). Run it with:

  ```bash
  fqtk trim --batch samples.tsv --jobs 8 --threads 1
  fqtk stat  --batch samples.tsv --jobs 8   # columns: sample, r1, r2(optional)
  ```

  Output is one summary line per sample (`sample  status  mode  total  kept
  dropped_short` for trim; similar for stat). A failing sample is reported
  and skipped rather than aborting the whole batch; the command exits
  non-zero if any sample failed, so it's safe to use in a pipeline.

- **`--jobs N` on `merge`** parallelizes across independent library groups
  (each merge/move touches different files and a different output path, so
  there's no shared-decompression bottleneck here — it's plain byte
  copying, which parallelizes cleanly, and on networked HPC storage several
  concurrent streams often push more aggregate throughput than one anyway):

  ```bash
  fqtk merge raw_fastq_dir -o merged_dir --execute --jobs 8
  ```

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
    ├── test_merge.py          # lane/batch merge sanity tests
    ├── test_merge_read_parsing.py  # R1/R2 detection on trim-output-style filenames
    └── test_parallel.py       # --batch / --jobs sanity tests
```

## License

MIT — see [LICENSE](LICENSE).
