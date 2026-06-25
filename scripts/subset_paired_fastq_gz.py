#!/usr/bin/env python3
import argparse
import gzip
from pathlib import Path


def read_record(handle):
    record = [handle.readline() for _ in range(4)]
    if not record[0]:
        return None
    if any(line == "" for line in record):
        raise SystemExit("Encountered an incomplete FASTQ record")
    return record


def write_record(handle, record):
    for line in record:
        handle.write(line)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--read1", required=True)
    parser.add_argument("--read2", required=True)
    parser.add_argument("--out-read1", required=True)
    parser.add_argument("--out-read2", required=True)
    parser.add_argument("--target-gb", type=float, required=True, help="Approximate compressed GB per output file.")
    parser.add_argument("--compresslevel", type=int, default=1)
    parser.add_argument("--check-every", type=int, default=25000, help="Check output size every N read pairs.")
    args = parser.parse_args()

    out1 = Path(args.out_read1)
    out2 = Path(args.out_read2)
    out1.parent.mkdir(parents=True, exist_ok=True)
    out2.parent.mkdir(parents=True, exist_ok=True)
    target_bytes = int(args.target_gb * 1024**3)

    pairs = 0
    with gzip.open(args.read1, "rt") as r1, gzip.open(args.read2, "rt") as r2, \
            gzip.open(out1, "wt", compresslevel=args.compresslevel) as w1, \
            gzip.open(out2, "wt", compresslevel=args.compresslevel) as w2:
        while True:
            rec1 = read_record(r1)
            rec2 = read_record(r2)
            if rec1 is None or rec2 is None:
                if rec1 is not None or rec2 is not None:
                    raise SystemExit("R1 and R2 ended at different record counts")
                break
            write_record(w1, rec1)
            write_record(w2, rec2)
            pairs += 1
            if pairs % args.check_every == 0:
                w1.flush()
                w2.flush()
                if out1.stat().st_size >= target_bytes and out2.stat().st_size >= target_bytes:
                    break

    print(f"wrote {pairs} read pairs")
    print(f"{out1}\t{out1.stat().st_size / 1024**3:.2f} GiB")
    print(f"{out2}\t{out2.stat().st_size / 1024**3:.2f} GiB")


if __name__ == "__main__":
    main()
