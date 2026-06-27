#!/usr/bin/env bash
set -euo pipefail

genomes_dir="$1"
extension="$2"
outdir="$3"
log="$4"

clean_prodigal_headers() {
  for fasta in "$@"; do
    awk '
      /^>/ {
        sub(/^>/, "")
        sub(/[[:space:]]*#.*/, "")
        print ">" $0
        next
      }
      { print }
    ' "$fasta" > "${fasta}.tmp"
    mv "${fasta}.tmp" "$fasta"
  done
}

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
  clean_prodigal_headers "$outdir/${name}.faa"
  cat "$outdir/${name}.fna" >> "$outdir/all_genes.fna"
done
