#!/usr/bin/env python3
import argparse
import shutil
from pathlib import Path


COMMANDS = [
    "snakemake",
    "conda",
    "fastqc",
    "trimmomatic",
    "spades.py",
    "metabat",
    "run_MaxBin.pl",
    "DAS_Tool",
    "checkm",
    "dRep",
    "gtdbtk",
    "prodigal",
    "emapper.py",
    "download_eggnog_data.py",
    "create_dbs.py",
    "bbmap.sh",
    "samtools",
    "inStrain",
]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as handle:
        handle.write("command\tavailable_on_current_path\tpath\n")
        for command in COMMANDS:
            path = shutil.which(command)
            handle.write(f"{command}\t{bool(path)}\t{path or ''}\n")


if __name__ == "__main__":
    main()
