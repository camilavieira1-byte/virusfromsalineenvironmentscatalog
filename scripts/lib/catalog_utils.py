"""Shared helpers for the viral catalog pipeline scripts."""
from pathlib import Path


def read_fasta(path):
    """Parse a fasta file into {header_first_token: sequence}."""
    seqs = {}
    name, chunks = None, []
    path = Path(path)
    if not path.exists():
        return seqs
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if name:
            seqs[name] = "".join(chunks)
    return seqs


def read_fasta_records(path):
    """Parse a fasta file preserving full headers, as an ordered list of
    (header, sequence) tuples. Use this instead of read_fasta() when you
    need to keep every header exactly as written (e.g. before reassigning
    safe IDs)."""
    records = []
    name, chunks = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name:
                    records.append((name, "".join(chunks)))
                name = line[1:]
                chunks = []
            else:
                chunks.append(line)
        if name:
            records.append((name, "".join(chunks)))
    return records
