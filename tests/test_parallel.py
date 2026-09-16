import csv
import gzip
import subprocess
import sys


def write_gz(path, records):
    with gzip.open(path, "wt") as f:
        for h, s, p, q in records:
            f.write(f"{h}\n{s}\n{p}\n{q}\n")


def read_gz(path):
    with gzip.open(path, "rt") as f:
        return f.readlines()


def run_fqtk(*args):
    return subprocess.run(
        [sys.executable, "-m", "fastq_toolkit.cli", *args],
        capture_output=True, text=True,
    )


def write_manifest(path, rows, columns):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=columns, delimiter="\t")
        w.writeheader()
        for row in rows:
            w.writerow(row)


# ------------------------------ trim --batch ------------------------------ #

def test_trim_batch_parallel(tmp_path):
    r1a = tmp_path / "A_R1.fastq.gz"
    write_gz(r1a, [("@a1", "AAAAAAAAAA", "+", "IIIIIIIIII")])
    r1b = tmp_path / "B_R1.fastq.gz"
    write_gz(r1b, [("@b1", "CCCCCCCCCCCCCCC", "+", "IIIIIIIIIIIIIII")])

    out_a = tmp_path / "A.trim.fastq.gz"
    out_b = tmp_path / "B.trim.fastq.gz"

    manifest = tmp_path / "manifest.tsv"
    write_manifest(manifest, [
        {"sample": "A", "r1": str(r1a), "out_r1": str(out_a), "r1_len": "5"},
        {"sample": "B", "r1": str(r1b), "out_r1": str(out_b), "r1_len": "5"},
    ], ["sample", "r1", "out_r1", "r1_len"])

    result = run_fqtk("trim", "--batch", str(manifest), "--jobs", "2")
    assert result.returncode == 0, result.stderr
    lines = result.stdout.strip().splitlines()
    assert lines[0].startswith("sample\tstatus")
    body = "\n".join(lines[1:])
    assert "A\tok\tSE\t1\t1\t0" in body
    assert "B\tok\tSE\t1\t1\t0" in body

    assert read_gz(out_a)[1].strip() == "AAAAA"   # trimmed to 5bp
    assert read_gz(out_b)[1].strip() == "CCCCC"


def test_trim_batch_reports_per_sample_failure(tmp_path):
    r1a = tmp_path / "A_R1.fastq.gz"
    write_gz(r1a, [("@a1", "AAAAAAAAAA", "+", "IIIIIIIIII")])
    out_a = tmp_path / "A.trim.fastq.gz"

    missing = tmp_path / "does_not_exist_R1.fastq.gz"
    out_missing = tmp_path / "missing.trim.fastq.gz"

    manifest = tmp_path / "manifest.tsv"
    write_manifest(manifest, [
        {"sample": "A", "r1": str(r1a), "out_r1": str(out_a), "r1_len": "5"},
        {"sample": "MISSING", "r1": str(missing), "out_r1": str(out_missing), "r1_len": "5"},
    ], ["sample", "r1", "out_r1", "r1_len"])

    result = run_fqtk("trim", "--batch", str(manifest), "--jobs", "2")
    assert result.returncode != 0
    assert "A\tok\tSE\t1\t1\t0" in result.stdout
    assert "MISSING\terror" in result.stdout
    assert out_a.exists()          # the good sample still completed


def test_trim_batch_rejects_mixed_single_and_batch_args(tmp_path):
    manifest = tmp_path / "manifest.tsv"
    write_manifest(manifest, [{"r1": "x", "out_r1": "y", "r1_len": "5"}],
                    ["r1", "out_r1", "r1_len"])
    result = run_fqtk("trim", "--batch", str(manifest), "-1", "somefile.fastq.gz")
    assert result.returncode != 0
    assert "cannot be combined" in (result.stdout + result.stderr)


# ------------------------------ stat --batch ------------------------------ #

def test_stat_batch_parallel(tmp_path):
    r1a = tmp_path / "A_R1.fastq.gz"
    write_gz(r1a, [("@a1", "AAAAAAAAAA", "+", "IIIIIIIIII")])
    r1b = tmp_path / "B_R1.fastq.gz"
    r2b = tmp_path / "B_R2.fastq.gz"
    write_gz(r1b, [("@b1", "CCCCCCCCCCCCCCC", "+", "IIIIIIIIIIIIIII")])
    write_gz(r2b, [("@b1", "GGGGGGGGGGGGGGGGGGGG", "+", "IIIIIIIIIIIIIIIIIIII")])

    manifest = tmp_path / "manifest.tsv"
    write_manifest(manifest, [
        {"sample": "A", "r1": str(r1a), "r2": ""},
        {"sample": "B", "r1": str(r1b), "r2": str(r2b)},
    ], ["sample", "r1", "r2"])

    result = run_fqtk("stat", "--batch", str(manifest), "--jobs", "2")
    assert result.returncode == 0, result.stderr
    out = result.stdout
    assert "A" in out and "SE" in out
    assert "B" in out and "R1" in out and "R2" in out
    assert "\t10\t10\t10.00\t10" in out   # sample A: 10bp read
    assert "\t20\t20\t20.00\t20" in out   # sample B R2: 20bp read


# ------------------------------ merge --jobs ------------------------------ #

def test_merge_jobs_matches_serial_result(tmp_path):
    in_dir = tmp_path / "raw"
    out_dir = tmp_path / "merged"
    in_dir.mkdir()

    for rl in ["RL5001", "RL5002", "RL5003"]:
        write_gz(in_dir / f"{rl}_2026_08_01_Lib_S1_R1_001.fastq.gz", [("@x", "A" * 10, "+", "I" * 10)])
        write_gz(in_dir / f"{rl}_2026_08_02_Lib_S2_R1_001.fastq.gz", [("@y", "C" * 10, "+", "I" * 10)])

    result = run_fqtk("merge", str(in_dir), "-o", str(out_dir), "--execute", "--jobs", "4")
    assert result.returncode == 0, result.stderr

    for rl in ["RL5001", "RL5002", "RL5003"]:
        merged = out_dir / f"{rl}_2026_08_Lib_R1_001.fastq.gz"
        assert merged.exists()
        content = read_gz(merged)
        assert content[1].strip() == "A" * 10
        assert content[5].strip() == "C" * 10

    report = (out_dir / "merge_report.tsv").read_text()
    assert report.count("\tmerge\t") == 3
