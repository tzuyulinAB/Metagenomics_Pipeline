#!/usr/bin/env python3
import argparse
from pathlib import Path


def clean_header(line):
    if not line.startswith(">"):
        return line
    header = line[1:].rstrip("\n")
    cleaned = header.split("#", 1)[0].strip()
    return f">{cleaned}\n"


def clean_file(path):
    path = Path(path)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with path.open() as src, tmp.open("w") as dst:
        for line in src:
            dst.write(clean_header(line))
    tmp.replace(path)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta", nargs="+")
    args = parser.parse_args()

    for fasta in args.fasta:
        clean_file(fasta)


if __name__ == "__main__":
    main()
