#!/usr/bin/env bash
set -euo pipefail

genomes_dir="$1"
extension="$2"
outdir="$3"
log="$4"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

mkdir -p "$outdir"
: > "$log"
: > "$outdir/all_genes.fna"

shopt -s nullglob
for genome in "$genomes_dir"/*."$extension"; do
  name="$(basename "$genome" ."$extension")"
  prodigal -p meta \
    -i "$genome" \
    -a "$outdir/${name}.faa" \
    -d "$outdir/${name}.fna" \
    >> "$log" 2>&1
  python "$script_dir/clean_prodigal_headers.py" "$outdir/${name}.faa" "$outdir/${name}.fna"
  cat "$outdir/${name}.fna" >> "$outdir/all_genes.fna"
done
