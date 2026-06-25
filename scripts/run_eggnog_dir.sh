#!/usr/bin/env bash
set -euo pipefail

indir="$1"
outdir="$2"
threads="$3"
data_dir="$4"
dmnd_db="$5"
evalue="$6"
tmpdir="$7"
log="$8"

mkdir -p "$outdir" "$tmpdir"
: > "$log"
export EGGNOG_DATA_DIR="$data_dir"

shopt -s nullglob
for faa in "$indir"/*.faa; do
  name="$(basename "$faa" .faa)"
  emapper.py \
    --cpu "$threads" \
    --data_dir "$data_dir" \
    --dmnd_db "$dmnd_db" \
    --evalue "$evalue" \
    -i "$faa" \
    -o "$outdir/${name}_egg" \
    --temp_dir "$tmpdir" \
    >> "$log" 2>&1
done
