#!/usr/bin/env python3
import argparse
from pathlib import Path


FASTA_SUFFIXES = {".fa", ".fasta", ".fna"}


def fasta_headers(path):
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                yield line[1:].split()[0]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bins-dir", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as handle:
        for fasta in sorted(p for p in Path(args.bins_dir).iterdir() if p.suffix in FASTA_SUFFIXES):
            genome = fasta.stem
            for scaffold in fasta_headers(fasta):
                handle.write(f"{scaffold}\t{genome}\n")


if __name__ == "__main__":
    main()
