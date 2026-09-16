import gzip
import subprocess
import sys

R1_RECORDS = [
    ("@read1", "GATAGTAAGTTTGTAATTACTGTGGCCC", "+", "I" * 28),
    ("@read2", "ACACAGGCCTCCTCATTTAAGAAATATT", "+", "I" * 28),
]

R2_RECORDS = [
    ("@read1",
     "GACTACAACCACGACCAATGATATGAAAAACCATCGTTGTATTTCAACTACAAGAACACCAATGACCCCAATACGCAAAATTAACCCCCT",
     "+", "I" * 90),
    ("@read2",
     "GAGCAAGGGGTGTTGAAAGAGAAGGAAATCATTTTCCTCAAACTGTTTACATTTGTAAAGTGAAACAGAAACCCTCTCTCATTTATGGAA",
     "+", "I" * 90),
]


def write_fastq_gz(path, records):
    with gzip.open(path, "wt") as f:
        for h, s, p, q in records:
            f.write(f"{h}\n{s}\n{p}\n{q}\n")


def run_fqtk(*args):
    return subprocess.run(
        [sys.executable, "-m", "fastq_toolkit.cli", *args],
        capture_output=True, text=True, check=True,
    )


def read_lengths(path):
    with gzip.open(path, "rt") as f:
        lines = f.readlines()
    return [len(lines[i].strip()) for i in range(1, len(lines), 4)]


def test_pe_trim(tmp_path):
    r1 = tmp_path / "R1.fastq.gz"
    r2 = tmp_path / "R2.fastq.gz"
    write_fastq_gz(r1, R1_RECORDS)
    write_fastq_gz(r2, R2_RECORDS)
    out1 = tmp_path / "R1.trim.fastq.gz"
    out2 = tmp_path / "R2.trim.fastq.gz"

    run_fqtk("trim", "-1", str(r1), "-2", str(r2),
              "-o1", str(out1), "-o2", str(out2),
              "--r1-len", "67", "--r2-len", "67")

    assert read_lengths(out1) == [28, 28]   # shorter than 67 -> kept in full
    assert read_lengths(out2) == [67, 67]   # trimmed down to 67


def test_se_trim(tmp_path):
    r1 = tmp_path / "R1.fastq.gz"
    write_fastq_gz(r1, R1_RECORDS)
    out1 = tmp_path / "R1.trim.fastq.gz"

    run_fqtk("trim", "-1", str(r1), "-o1", str(out1), "--r1-len", "20")

    assert read_lengths(out1) == [20, 20]


def test_pe_stat(tmp_path):
    r1 = tmp_path / "R1.fastq.gz"
    r2 = tmp_path / "R2.fastq.gz"
    write_fastq_gz(r1, R1_RECORDS)
    write_fastq_gz(r2, R2_RECORDS)

    result = run_fqtk("stat", "-1", str(r1), "-2", str(r2))
    lines = result.stdout.strip().splitlines()
    assert lines[0].startswith("file\tread\treads")
    assert "\tR1\t2\t28\t28\t28.00\t28" in lines[1]
    assert "\tR2\t2\t90\t90\t90.00\t90" in lines[2]


def test_mismatched_pair_errors(tmp_path):
    r1 = tmp_path / "R1.fastq.gz"
    r2 = tmp_path / "R2.fastq.gz"
    write_fastq_gz(r1, R1_RECORDS[:1])   # 1 record
    write_fastq_gz(r2, R2_RECORDS)       # 2 records

    result = subprocess.run(
        [sys.executable, "-m", "fastq_toolkit.cli", "trim",
         "-1", str(r1), "-2", str(r2),
         "-o1", str(tmp_path / "o1.fq.gz"), "-o2", str(tmp_path / "o2.fq.gz"),
         "--r1-len", "20"],
        capture_output=True, text=True,
    )
    assert result.returncode != 0
    assert "read counts differ" in (result.stdout + result.stderr)
