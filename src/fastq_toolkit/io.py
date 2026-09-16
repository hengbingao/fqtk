"""Shared FASTQ I/O helpers.

Reads/writes gzip-compressed FASTQ. When ``threads > 1`` and ``pigz`` is
available on PATH, decompression/compression is offloaded to pigz
(multi-threaded, much faster on large files); otherwise falls back to the
Python standard library ``gzip`` module with zero external dependencies.
"""

import gzip
import shutil
import subprocess


def has_pigz() -> bool:
    return shutil.which("pigz") is not None


def open_input(path, threads):
    """Return (readable text-mode file-like object, subprocess.Popen or None)."""
    path = str(path)
    if threads > 1 and has_pigz():
        p = subprocess.Popen(
            ["pigz", "-dc", "-p", str(threads), path],
            stdout=subprocess.PIPE,
            text=True,
            bufsize=1 << 20,
        )
        return p.stdout, p
    return gzip.open(path, "rt"), None


def open_output(path, threads, compresslevel=6):
    """Return (writable text-mode file-like object, subprocess.Popen or None)."""
    path = str(path)
    if threads > 1 and has_pigz():
        p = subprocess.Popen(
            ["pigz", "-c", "-p", str(threads), f"-{compresslevel}"],
            stdin=subprocess.PIPE,
            stdout=open(path, "wb"),
            text=True,
            bufsize=1 << 20,
        )
        return p.stdin, p
    return gzip.open(path, "wt", compresslevel=compresslevel), None


def read_record(fh):
    """Read one 4-line FASTQ record. Returns (header, seq, plus, qual) or None at EOF."""
    header = fh.readline()
    if not header:
        return None
    seq = fh.readline()
    plus = fh.readline()
    qual = fh.readline()
    if not (seq and plus and qual):
        raise ValueError("FASTQ record truncated mid-way; check file integrity")
    return header.rstrip("\n"), seq.rstrip("\n"), plus.rstrip("\n"), qual.rstrip("\n")


def close_all(handles):
    """Close a list of (file_handle, subprocess.Popen_or_None) pairs cleanly."""
    for fh, p in handles:
        fh.close()
        if p is not None:
            p.wait()
