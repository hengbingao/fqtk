import gzip
import subprocess
import sys


def write_gz(path, text):
    with gzip.open(path, "wt") as f:
        f.write(text)


def read_gz(path):
    with gzip.open(path, "rt") as f:
        return f.read()


def run_fqtk(*args):
    return subprocess.run(
        [sys.executable, "-m", "fastq_toolkit.cli", *args],
        capture_output=True, text=True,
    )


def test_merge_preview_writes_nothing(tmp_path):
    in_dir = tmp_path / "raw"
    out_dir = tmp_path / "merged"
    in_dir.mkdir()

    write_gz(in_dir / "RL1001_2026_08_01_TestLib_S1_R1_001.fastq.gz", "AAAA\n")
    write_gz(in_dir / "RL1001_2026_08_02_TestLib_S2_R1_001.fastq.gz", "CCCC\n")

    result = run_fqtk("merge", str(in_dir), "-o", str(out_dir))
    assert result.returncode == 0, result.stderr
    assert not out_dir.exists()          # preview only, nothing written
    assert "MERGE" in result.stdout


def test_merge_execute(tmp_path):
    in_dir = tmp_path / "raw"
    out_dir = tmp_path / "merged"
    in_dir.mkdir()

    lane1 = in_dir / "RL1001_2026_08_01_TestLib_S1_R1_001.fastq.gz"
    lane2 = in_dir / "RL1001_2026_08_02_TestLib_S2_R1_001.fastq.gz"
    write_gz(lane1, "AAAA\n")
    write_gz(lane2, "CCCC\n")

    singleton = in_dir / "RL2002_2026_08_01_SoloLib_S5_R1_001.fastq.gz"
    write_gz(singleton, "GGGG\n")

    skip_file = in_dir / "not_a_library_file.fastq.gz"
    write_gz(skip_file, "TTTT\n")

    result = run_fqtk("merge", str(in_dir), "-o", str(out_dir), "--execute")
    assert result.returncode == 0, result.stderr

    merged_out = out_dir / "RL1001_2026_08_TestLib_R1_001.fastq.gz"
    assert merged_out.exists()
    assert read_gz(merged_out) == "AAAA\nCCCC\n"   # valid multi-member gzip concat

    # sources of a merged group are kept by default
    assert lane1.exists()
    assert lane2.exists()

    # singleton got moved (not copied) by default
    moved_singleton = out_dir / "RL2002_2026_08_01_SoloLib_S5_R1_001.fastq.gz"
    assert moved_singleton.exists()
    assert not singleton.exists()

    # unparseable file was left untouched in the input dir
    assert skip_file.exists()

    report = out_dir / "merge_report.tsv"
    assert report.exists()
    report_text = report.read_text()
    assert "merge" in report_text
    assert "singleton" in report_text
    assert "skipped" in report_text


def test_merge_copy_singletons(tmp_path):
    in_dir = tmp_path / "raw"
    out_dir = tmp_path / "merged"
    in_dir.mkdir()

    singleton = in_dir / "RL3003_2026_08_01_SoloLib_S9_R1_001.fastq.gz"
    write_gz(singleton, "GGGG\n")

    result = run_fqtk("merge", str(in_dir), "-o", str(out_dir), "--execute", "--copy-singletons")
    assert result.returncode == 0, result.stderr
    assert singleton.exists()                       # kept, thanks to --copy-singletons
    assert (out_dir / singleton.name).exists()


def test_merge_remove_merged_sources(tmp_path):
    in_dir = tmp_path / "raw"
    out_dir = tmp_path / "merged"
    in_dir.mkdir()

    lane1 = in_dir / "RL4004_2026_08_01_TestLib_S1_R2_001.fastq.gz"
    lane2 = in_dir / "RL4004_2026_08_02_TestLib_S2_R2_001.fastq.gz"
    write_gz(lane1, "AAAA\n")
    write_gz(lane2, "CCCC\n")

    result = run_fqtk("merge", str(in_dir), "-o", str(out_dir),
                       "--execute", "--remove-merged-sources")
    assert result.returncode == 0, result.stderr
    assert not lane1.exists()
    assert not lane2.exists()
    assert (out_dir / "RL4004_2026_08_TestLib_R2_001.fastq.gz").exists()
