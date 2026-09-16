"""Regression test for the R1/R2 parsing bug: fqtk trim's own output names
glue 'R1'/'R2' to the next part with a '.' instead of '_'
(e.g. '..._S1_R1.trim26.out.fastq.gz'), which used to make `fqtk merge`
treat them as unparseable and skip them."""

import gzip
import subprocess
import sys


def write_gz(path, text):
    with gzip.open(path, "wt") as f:
        f.write(text)


def run_fqtk(*args):
    return subprocess.run(
        [sys.executable, "-m", "fastq_toolkit.cli", *args],
        capture_output=True, text=True,
    )


def test_merge_handles_trim_output_style_filenames(tmp_path):
    in_dir = tmp_path / "trim"
    out_dir = tmp_path / "merge"
    in_dir.mkdir()

    names = [
        "RL8518_2026_07_24_BulkATAC_P5C10_BEP100_ZIM3_rep2_S1_R1.trim26.out.fastq.gz",
        "RL8518_2026_07_24_BulkATAC_P5C10_BEP100_ZIM3_rep2_S1_R2.trim67.out.fastq.gz",
        "RL8518_2026_07_28_BulkATAC_P5C10_BEP100_ZIM3_rep2_S11_R1.trim26.out.fastq.gz",
        "RL8518_2026_07_28_BulkATAC_P5C10_BEP100_ZIM3_rep2_S11_R2.trim67.out.fastq.gz",
        "RL8519_2026_07_24_BulkATAC_P5C10_BEP396_FOG1_rep2_S2_R1.trim26.out.fastq.gz",
        "RL8519_2026_07_24_BulkATAC_P5C10_BEP396_FOG1_rep2_S2_R2.trim67.out.fastq.gz",
    ]
    for i, name in enumerate(names):
        write_gz(in_dir / name, f"payload-{i}\n")

    result = run_fqtk("merge", str(in_dir), "-o", str(out_dir), "--execute")
    assert result.returncode == 0, result.stderr
    assert "[skipped]" not in result.stdout

    # RL8518 R1 and R2 each had 2 lanes -> merged; RL8519 R1/R2 had 1 lane each -> moved as-is
    merged_r1 = out_dir / "RL8518_2026_07_BulkATAC_P5C10_BEP100_ZIM3_rep2_R1.trim26.out.fastq.gz"
    merged_r2 = out_dir / "RL8518_2026_07_BulkATAC_P5C10_BEP100_ZIM3_rep2_R2.trim67.out.fastq.gz"
    assert merged_r1.exists()
    assert merged_r2.exists()

    singleton_r1 = out_dir / "RL8519_2026_07_24_BulkATAC_P5C10_BEP396_FOG1_rep2_S2_R1.trim26.out.fastq.gz"
    singleton_r2 = out_dir / "RL8519_2026_07_24_BulkATAC_P5C10_BEP396_FOG1_rep2_S2_R2.trim67.out.fastq.gz"
    assert singleton_r1.exists()
    assert singleton_r2.exists()

    report = (out_dir / "merge_report.tsv").read_text()
    assert "skipped" not in report
    assert report.count("\tmerge\t") == 2
    assert report.count("\tsingleton\t") == 2
