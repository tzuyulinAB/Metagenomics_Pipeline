# Metagenomics Pipeline

Docker-first Nextflow workflow for DNA metagenomics. This is a Nextflow port of the sibling Snakemake workflow. RNA/metatranscriptomics steps are intentionally left out for now.

## Workflow

1. Validate manifests and write a dependency diagnostic report.
2. Prepare the selected Trimmomatic adapter FASTA.
3. Trim DNA paired-end reads with Trimmomatic.
4. Run FastQC on trimmed paired reads.
5. Assemble each assembly/co-assembly with metaSPAdes, adding PacBio reads when a pool is configured.
6. Filter contigs by length and prefix headers with the assembly ID.
7. Bin contigs with MetaBAT2 and MaxBin2.
8. Consolidate bins with DAS Tool.
9. Collect DAS Tool bins and dereplicate with dRep.
10. Run GTDB-Tk classification and bacterial de novo phylogeny.
11. Predict genes with Prodigal and annotate MAGs/assemblies with EggNOG-mapper.
12. Build a dRep reference, map DNA reads with BBMap, then run inStrain profile and compare.

The `drep_instrain` profile skips GTDB-Tk and EggNOG while keeping the dRep, MAG Prodigal, and inStrain steps.

## Inputs

Edit these files before a production run:

- `config/samples.tsv`: DNA paired-end libraries. Rows with `data_type` other than `DNA` are ignored.
- `config/assemblies.tsv`: one row per assembly or co-assembly.
- `config/pacbio_pools.tsv`: optional PacBio HiFi read pools used by assemblies.
- `config/genome_info.csv`: optional dRep genome info table. If it only contains the header, dRep runs its own genome-quality checks.

The TSV schemas match the Snakemake project:

```tsv
sample_id	data_type	condition	read1	read2
D0_AD	DNA	AD	/path/D0_AD_R1.fastq.gz	/path/D0_AD_R2.fastq.gz
```

```tsv
assembly_id	sample_ids	pacbio_pool
AD_coassembly	D0_AD,D6_AD,D13_AD	AD_pool
```

## Configuration

Runtime options live in `nextflow.config`. Common overrides:

```bash
nextflow run . \
  --samples config/samples.tsv \
  --assemblies config/assemblies.tsv \
  --pacbio_pools config/pacbio_pools.tsv \
  --gtdbtk_data_path /path/to/release220 \
  --outdir results
```

Docker is enabled by default. Each process label has a container image in `nextflow.config`, so you can swap a single image without touching `main.nf`.

## Run

Dry-run the graph:

```bash
nextflow run . -preview
```

Run the full workflow:

```bash
nextflow run . -profile docker
```

Run the reduced dRep/inStrain workflow without GTDB-Tk or EggNOG:

```bash
nextflow run . -profile docker,drep_instrain
```

For resource-constrained runs, lower dRep threads and memory:

```bash
nextflow run . -profile docker,drep_instrain \
  --threads_drep 4 \
  --memory_drep '16 GB'
```

dRep output folders include the completeness and contamination thresholds, for example `results/mags/drep/drep_50_10/` by default or `results/mags/drep/drep_0_100/` when running with `--drep_completeness 0 --drep_contamination 100`.

A small local smoke profile is included for the copied test manifests:

```bash
nextflow run . -profile local_smoke
```

The smoke profile writes outside the iCloud path by default because SPAdes can be unhappy with spaces in working/output paths on macOS.

## Notes

- MaxBin2 receives the first sample listed for each assembly, matching the Snakemake implementation.
- EggNOG database generation is wired in when `params.eggnog_dmnd_db` does not exist. For production, prebuilding and pointing `--eggnog_dmnd_db` at the database is usually faster and more reproducible.
- GTDB-Tk still needs a local database path via `--gtdbtk_data_path`.
- The default Docker image tags are app-level starting points. If your platform or registry mirror prefers different tags, adjust the `withLabel` container entries in `nextflow.config`.
