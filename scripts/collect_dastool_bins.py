#!/usr/bin/env python3
import argparse
import shutil
from pathlib import Path


FASTA_SUFFIXES = {".fa", ".fasta", ".fna"}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dastool-root", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    root = Path(args.dastool_root)
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)

    copied = 0
    for bins_dir in sorted(root.glob("*/*_DAS_DASTool_bins")):
        assembly = bins_dir.parent.name
        for fasta in sorted(p for p in bins_dir.iterdir() if p.suffix in FASTA_SUFFIXES):
            dest = out / f"{assembly}_{fasta.name}"
            shutil.copy2(fasta, dest)
            copied += 1
    if copied == 0:
        raise SystemExit(f"No DAS Tool bins found below {root}")


if __name__ == "__main__":
    main()
