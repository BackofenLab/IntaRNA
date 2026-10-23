"""Shared, dependency-free I/O for the IntaRNA Python companion tools."""

import argparse
from contextlib import contextmanager
from pathlib import Path
import os
import re
import shutil
import subprocess
import sys
import tempfile


def positive_integer(value):
    """argparse type for counts that must be positive."""
    number = int(value)
    if number < 1:
        raise argparse.ArgumentTypeError('expected a positive integer')
    return number


def nonnegative_integer(value):
    """argparse type for thread counts (zero means automatic)."""
    number = int(value)
    if number < 0:
        raise argparse.ArgumentTypeError('expected a nonnegative integer')
    return number


def normalize_sequence(value, *, allow_n=False):
    """Validate a nonempty literal RNA/DNA sequence and normalize to RNA."""
    alphabet = 'ACGUTN' if allow_n else 'ACGUT'
    if not re.fullmatch(f'[{alphabet}{alphabet.lower()}]+', value):
        raise ValueError(f'expected a nonempty sequence containing only {alphabet}')
    return value.upper().replace('T', 'U')


def read_fasta(filename):
    """Read exactly one FASTA record, including wrapped sequence lines."""
    with open(filename, encoding='utf-8-sig') as handle:
        lines = [line.strip() for line in handle if line.strip()]
    if not lines or not lines[0].startswith('>'):
        raise ValueError(f'{filename}: expected a FASTA header')
    if any(line.startswith('>') for line in lines[1:]):
        raise ValueError(f'{filename}: expected exactly one FASTA record')
    return normalize_sequence(''.join(lines[1:]))


def read_sequence(value):
    """Accept a literal RNA/DNA sequence or a single-record FASTA file."""
    # Test sequences first: long literals can exceed the filesystem name limit.
    if re.fullmatch('[ACGUTacgut]+', value):
        return normalize_sequence(value)
    return read_fasta(value) if Path(value).is_file() else normalize_sequence(value)


def find_program(name, script_path=__file__):
    """Prefer executable siblings, then the source-tree IntaRNA, then PATH."""
    directory = Path(script_path).resolve().parent
    candidates = [directory / name]
    if name == 'IntaRNA':
        candidates.append(directory.parent / 'src' / 'bin' / name)
    for candidate in candidates:
        found = shutil.which(str(candidate))
        if found:
            return found
    found = shutil.which(name)
    if not found:
        raise RuntimeError(f'{name} executable not found beside the script or on PATH')
    return found


def run_command(command, input_text=None):
    """Run an argument vector without a shell; report command failures once."""
    result = subprocess.run(command, input=input_text, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode:
        detail = (result.stderr or result.stdout).strip()
        raise RuntimeError(f'{Path(command[0]).name} exited with status {result.returncode}: {detail}')
    if result.stderr:
        sys.stderr.write(result.stderr)
    return result.stdout


def csv_header(reader):
    """Read a nonempty CSV header with distinct, nonempty names."""
    header = next(reader, [])
    if not header or not all(header) or len(set(header)) != len(header):
        raise ValueError('CSV must have a nonempty header with unique column names')
    return header


def csv_rows(reader, header):
    """Yield nonblank rows, rejecting truncated or extra fields."""
    for row in reader:
        if not row:
            continue
        if len(row) != len(header):
            raise ValueError(f'CSV line {reader.line_num} has the wrong number of fields')
        yield row


@contextmanager
def atomic_text_output(destination):
    """Publish UTF-8/LF text only on success, preserving an existing file's mode."""
    destination = Path(destination).resolve()
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode='w', newline='', encoding='utf-8',
                                         dir=destination.parent, delete=False) as handle:
            temporary = Path(handle.name)
            yield handle
        if destination.exists():
            shutil.copymode(destination, temporary)
        os.replace(temporary, destination)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
