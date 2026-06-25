#!/usr/bin/env python3
import argparse


def records(path):
    name = None
    seq = []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
        if name is not None:
            yield name, "".join(seq)


def write_record(handle, name, seq):
    handle.write(f">{name}\n")
    for i in range(0, len(seq), 80):
        handle.write(seq[i:i + 80] + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-len", type=int, required=True)
    parser.add_argument("--prefix", default="", help="Optional prefix to add to each FASTA header as '<prefix>__<header>'.")
    parser.add_argument("input")
    parser.add_argument("output")
    args = parser.parse_args()

    kept = 0
    with open(args.output, "w") as out:
        for name, seq in records(args.input):
            if len(seq) >= args.min_len:
                kept += 1
                if args.prefix:
                    name = f"{args.prefix}__{name}"
                write_record(out, name, seq)
    prefix_msg = f" with prefix {args.prefix}" if args.prefix else ""
    print(f"kept {kept} contigs >= {args.min_len} bp{prefix_msg}")


if __name__ == "__main__":
    main()
