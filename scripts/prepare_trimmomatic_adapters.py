#!/usr/bin/env python3
import argparse
import shutil
import tarfile
import tempfile
import urllib.request
from pathlib import Path


TRIMMOMATIC_TARBALL = "https://github.com/timflutre/trimmomatic/archive/refs/heads/master.tar.gz"


def download_adapters(adapter_dir):
    adapter_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        archive = Path(tmp) / "trimmomatic-master.tar.gz"
        urllib.request.urlretrieve(TRIMMOMATIC_TARBALL, archive)
        with tarfile.open(archive, "r:gz") as tar:
            tar.extractall(tmp)
        source = Path(tmp) / "trimmomatic-master" / "adapters"
        if not source.exists():
            raise SystemExit(f"Could not find adapters folder in {TRIMMOMATIC_TARBALL}")
        for fasta in source.glob("*.fa"):
            shutil.copy2(fasta, adapter_dir / fasta.name)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--adapter-dir", required=True)
    parser.add_argument("--selected", required=True)
    parser.add_argument("--list", action="store_true", help="List available adapter FASTA files after download/check.")
    args = parser.parse_args()

    adapter_dir = Path(args.adapter_dir)
    selected = adapter_dir / args.selected
    if not selected.exists():
        download_adapters(adapter_dir)

    available = sorted(path.name for path in adapter_dir.glob("*.fa"))
    if args.list:
        print("\n".join(available))

    if not selected.exists():
        options = ", ".join(available) if available else "none"
        raise SystemExit(f"Selected adapter {args.selected} not found in {adapter_dir}. Available adapters: {options}")


if __name__ == "__main__":
    main()
