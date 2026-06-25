#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path


def read_tsv(path):
    with open(path, newline="") as handle:
        yield from csv.DictReader((line for line in handle if not line.lstrip().startswith("#")), delimiter="\t")


def require_columns(path, required):
    with open(path, newline="") as handle:
        reader = csv.DictReader((line for line in handle if not line.lstrip().startswith("#")), delimiter="\t")
        missing = set(required) - set(reader.fieldnames or [])
    if missing:
        raise SystemExit(f"{path} is missing required columns: {', '.join(sorted(missing))}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", required=True)
    parser.add_argument("--assemblies", required=True)
    parser.add_argument("--pacbio-pools", required=True)
    args = parser.parse_args()

    require_columns(args.samples, {"sample_id", "data_type", "condition", "read1", "read2"})
    require_columns(args.assemblies, {"assembly_id", "sample_ids", "pacbio_pool"})
    require_columns(args.pacbio_pools, {"pool_id", "read_path"})

    samples = {row["sample_id"]: row for row in read_tsv(args.samples) if row.get("sample_id")}
    pools = {row["pool_id"]: row for row in read_tsv(args.pacbio_pools) if row.get("pool_id")}

    for sample_id, row in samples.items():
        data_type = row["data_type"].upper()
        if data_type != "DNA":
            raise SystemExit(f"{sample_id}: data_type must be DNA in this metagenomics workflow")
        for col in ("read1", "read2"):
            if not row[col]:
                raise SystemExit(f"{sample_id}: {col} is empty")

    for row in read_tsv(args.assemblies):
        assembly = row.get("assembly_id", "")
        if not assembly:
            continue
        members = [sample.strip() for sample in row.get("sample_ids", "").split(",") if sample.strip()]
        if not members:
            raise SystemExit(f"{assembly}: sample_ids is empty")
        for sample in members:
            if sample not in samples:
                raise SystemExit(f"{assembly}: unknown sample_id {sample}")
            if samples[sample]["data_type"].upper() != "DNA":
                raise SystemExit(f"{assembly}: {sample} is not a DNA sample")
        pool = row.get("pacbio_pool", "").strip()
        if pool and pool not in pools:
            raise SystemExit(f"{assembly}: unknown pacbio_pool {pool}")

    Path("config").mkdir(exist_ok=True)


if __name__ == "__main__":
    main()
