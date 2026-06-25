#!/usr/bin/env python3
import argparse
import csv
import gzip
import json
from pathlib import Path


def read_tsv(path):
    with open(path, newline="") as handle:
        yield from csv.DictReader(
            (line for line in handle if line.strip() and not line.lstrip().startswith("#")),
            delimiter="\t",
        )


def gib(num_bytes):
    return num_bytes / 1024**3


def fmt_gib(value):
    if value is None:
        return "unknown"
    if value < 1:
        return f"{value * 1024:.0f} MiB"
    return f"{value:.1f} GiB"


def fmt_hours(value):
    if value is None:
        return "unknown"
    if value < 1:
        return f"{value * 60:.0f} min"
    return f"{value:.1f} h"


def resolve_path(path, base_dir):
    p = Path(path)
    return p if p.is_absolute() else base_dir / p


def gzip_sample_stats(path, sample_records, gzip_ratio):
    """Estimate uncompressed bytes and read count by sampling complete FASTQ records."""
    compressed_size = path.stat().st_size
    sampled_records = 0
    sampled_uncompressed = 0
    with gzip.open(path, "rt", errors="replace") as handle:
        while sampled_records < sample_records:
            record = [handle.readline() for _ in range(4)]
            if not record[0]:
                break
            if any(line == "" for line in record):
                raise ValueError(f"{path}: incomplete FASTQ record in sample")
            sampled_records += 1
            sampled_uncompressed += sum(len(line.encode("utf-8")) for line in record)

    if sampled_records == 0:
        return {
            "compressed_gib": gib(compressed_size),
            "estimated_reads": 0,
            "estimated_uncompressed_gib": 0,
        }

    estimated_uncompressed = compressed_size * gzip_ratio
    mean_record_bytes = sampled_uncompressed / sampled_records
    estimated_reads = estimated_uncompressed / mean_record_bytes
    return {
        "compressed_gib": gib(compressed_size),
        "estimated_reads": int(round(estimated_reads)),
        "estimated_uncompressed_gib": gib(estimated_uncompressed),
    }


def estimate_spades_memory_gib(total_read_pairs, mean_read_len, distinct_kmer_factor, min_spades_memory_gib):
    # Empirical guardrail from local SPAdes/BayesHammer tests in this repo:
    # ~3.96M read pairs at 250 bp needed about 16 GiB, while ~1.60M pairs
    # fit past indexing but still peaked around the 8 GiB Docker cap.
    bases_gib = total_read_pairs * 2 * mean_read_len / 1024**3
    return max(min_spades_memory_gib, bases_gib * distinct_kmer_factor)


def bounded(value, lower, upper=None):
    value = max(lower, value)
    return min(value, upper) if upper is not None else value


def task_estimate(name, cpus, memory_gib, work_disk_gib, wall_hours, basis):
    return {
        "task": name,
        "cpus": cpus,
        "memory_gib": memory_gib,
        "work_disk_gib": work_disk_gib,
        "wall_hours": wall_hours,
        "basis": basis,
    }


def estimate_task_resources(
    total_compressed_gib,
    total_uncompressed_gib,
    samples_count,
    assemblies_count,
    spades_memory_gib,
    min_spades_memory_gib,
):
    """Return rough per-process resources matching the processes in main.nf."""
    assemblies_count = max(1, assemblies_count)
    samples_count = max(1, samples_count)
    spades_memory_gib = spades_memory_gib or min_spades_memory_gib
    spades_disk_gib = max(50.0, total_uncompressed_gib * 8)
    binning_disk_gib = max(10.0, total_uncompressed_gib * 1.5)
    mags_disk_gib = max(10.0, total_uncompressed_gib * 0.8)
    drep_disk_gib = max(20.0, mags_disk_gib * 2.0)
    instrain_ref_disk_gib = max(5.0, mags_disk_gib * 0.5)

    estimates = [
        task_estimate("VALIDATE_CONFIG", 1, 0.5, 1.0, 0.02, "Manifest parsing only."),
        task_estimate("CHECK_DEPENDENCIES", 1, 0.5, 1.0, 0.03, "Container command checks only."),
        task_estimate("PREPARE_TRIMMOMATIC_ADAPTER", 1, 0.2, 1.0, 0.02, "Copies one adapter FASTA."),
        task_estimate(
            "TRIM_DNA",
            16,
            bounded(total_compressed_gib * 1.5, 4.0, 24.0),
            max(5.0, total_uncompressed_gib * 1.2),
            max(0.05, total_compressed_gib * 0.18),
            "Scales with compressed paired FASTQ input.",
        ),
        task_estimate(
            "FASTQC_DNA_TRIMMED",
            4,
            bounded(total_compressed_gib * 0.4, 2.0, 8.0),
            max(2.0, total_compressed_gib * 0.2),
            max(0.05, total_compressed_gib * 0.12),
            "Scales with trimmed FASTQ size.",
        ),
        task_estimate(
            "HYBRIDSPADES",
            40,
            spades_memory_gib,
            spades_disk_gib,
            max(0.25, assemblies_count * total_compressed_gib * 0.9),
            "metaSPAdes estimate from read pairs, read length, and k-mer complexity factor.",
        ),
        task_estimate(
            "FILTER_CONTIGS",
            1,
            bounded(total_uncompressed_gib * 0.05, 1.0, 8.0),
            max(2.0, total_uncompressed_gib * 0.5),
            max(0.03, assemblies_count * 0.05),
            "Scales with assembled contig FASTA size.",
        ),
        task_estimate(
            "METABAT",
            20,
            bounded(total_uncompressed_gib * 1.5, 8.0, 128.0),
            binning_disk_gib,
            max(0.5, assemblies_count * 1.0),
            "Binning estimate from assembly scale.",
        ),
        task_estimate(
            "MAXBIN",
            24,
            bounded(total_uncompressed_gib * 2.0, 8.0, 128.0),
            binning_disk_gib,
            max(0.5, assemblies_count * 1.2),
            "Binning estimate from assembly scale and first sample coverage.",
        ),
        task_estimate(
            "SCAFFOLDS2BIN",
            1,
            1.0,
            max(1.0, assemblies_count * 0.1),
            max(0.03, assemblies_count * 0.03),
            "Converts bin FASTA files to small mapping tables.",
        ),
        task_estimate(
            "DASTOOL",
            30,
            bounded(total_uncompressed_gib * 1.2, 8.0, 96.0),
            mags_disk_gib,
            max(1.0, assemblies_count * 1.5),
            "Consolidates MetaBAT and MaxBin bins.",
        ),
        task_estimate("COLLECT_DASTOOL_BINS", 1, 1.0, mags_disk_gib, 0.1, "Copies selected DAS Tool bins."),
        task_estimate(
            "CHECKM",
            40,
            bounded(mags_disk_gib * 1.5, 16.0, 128.0),
            max(20.0, mags_disk_gib * 2.0),
            max(1.0, assemblies_count * 2.0),
            "MAG quality estimate; database size is not included.",
        ),
        task_estimate(
            "DREP",
            40,
            bounded(mags_disk_gib * 2.0, 16.0, 192.0),
            drep_disk_gib,
            max(1.0, assemblies_count * 2.0),
            "Dereplication estimate from expected MAG output scale.",
        ),
        task_estimate(
            "GTDBTK_CLASSIFY",
            30,
            bounded(mags_disk_gib * 2.0, 32.0, 192.0),
            max(20.0, mags_disk_gib * 2.0),
            max(2.0, assemblies_count * 2.0),
            "Classification estimate; GTDB-Tk reference database size is not included.",
        ),
        task_estimate(
            "GTDBTK_DE_NOVO",
            30,
            bounded(mags_disk_gib * 2.5, 32.0, 256.0),
            max(20.0, mags_disk_gib * 2.0),
            max(2.0, assemblies_count * 2.5),
            "Tree inference estimate; GTDB-Tk reference database size is not included.",
        ),
        task_estimate("PRODIGAL_MAG", 4, bounded(mags_disk_gib * 0.5, 2.0, 16.0), mags_disk_gib, 0.5, "Gene prediction on MAGs."),
        task_estimate(
            "PRODIGAL_ASSEMBLY",
            4,
            bounded(total_uncompressed_gib * 0.2, 2.0, 16.0),
            max(2.0, total_uncompressed_gib * 0.3),
            max(0.1, assemblies_count * 0.2),
            "Gene prediction on filtered assemblies.",
        ),
        task_estimate(
            "EGGNOG_DATABASE",
            1,
            16.0,
            80.0,
            4.0,
            "Only used when the EggNOG diamond database must be built.",
        ),
        task_estimate(
            "EGGNOG_ASSEMBLY",
            32,
            bounded(total_uncompressed_gib * 1.0, 16.0, 128.0),
            max(20.0, total_uncompressed_gib * 1.0),
            max(1.0, assemblies_count * 2.0),
            "EggNOG annotation on assembly genes; EggNOG database size is not included.",
        ),
        task_estimate(
            "EGGNOG_MAG",
            32,
            bounded(mags_disk_gib * 1.0, 16.0, 128.0),
            max(20.0, mags_disk_gib * 1.0),
            max(1.0, assemblies_count * 2.0),
            "EggNOG annotation on MAG genes; EggNOG database size is not included.",
        ),
        task_estimate("CONCATENATE_DREP_REFERENCE", 1, 1.0, instrain_ref_disk_gib, 0.1, "Concatenates dereplicated genomes."),
        task_estimate("MAKE_STB", 1, 1.0, 1.0, 0.1, "Creates scaffold-to-bin table."),
        task_estimate(
            "BBMAP_DNA_FOR_INSTRAIN",
            12,
            bounded(total_uncompressed_gib * 1.0, 8.0, 96.0),
            max(20.0, total_uncompressed_gib * 2.5),
            max(0.5, samples_count * total_compressed_gib * 0.3),
            "Maps DNA reads to the dRep reference and sorts BAM files.",
        ),
        task_estimate(
            "INSTRAIN_PROFILE",
            30,
            bounded(total_uncompressed_gib * 1.0, 16.0, 160.0),
            max(20.0, total_uncompressed_gib * 1.5),
            max(1.0, samples_count * 1.5),
            "Profiles each sample against dereplicated genomes.",
        ),
        task_estimate(
            "INSTRAIN_COMPARE",
            30,
            bounded(samples_count * mags_disk_gib * 0.5, 16.0, 160.0),
            max(10.0, samples_count * mags_disk_gib * 0.5),
            max(0.5, samples_count * 0.5),
            "Compares inStrain profiles across samples.",
        ),
    ]

    max_memory = max(estimates, key=lambda item: item["memory_gib"])
    max_work_disk = max(estimates, key=lambda item: item["work_disk_gib"])
    max_wall = max(estimates, key=lambda item: item["wall_hours"])
    total_wall = sum(item["wall_hours"] for item in estimates)
    max_cpus = max(estimates, key=lambda item: item["cpus"])
    return {
        "tasks": estimates,
        "summary": {
            "max_memory_task": max_memory["task"],
            "max_memory_gib": max_memory["memory_gib"],
            "max_work_disk_task": max_work_disk["task"],
            "max_work_disk_gib": max_work_disk["work_disk_gib"],
            "max_wall_time_task": max_wall["task"],
            "max_wall_hours": max_wall["wall_hours"],
            "max_cpus_task": max_cpus["task"],
            "max_cpus": max_cpus["cpus"],
            "serial_wall_hours_rough": total_wall,
        },
    }


def main():
    parser = argparse.ArgumentParser(
        description="Estimate rough resources for the metagenomics Nextflow pipeline before running it."
    )
    parser.add_argument("--samples", default="config/samples.tsv")
    parser.add_argument("--assemblies", default="config/assemblies.tsv")
    parser.add_argument("--pacbio-pools", default="config/pacbio_pools.tsv")
    parser.add_argument("--sample-records", type=int, default=20000)
    parser.add_argument("--mean-read-len", type=int, default=250)
    parser.add_argument(
        "--gzip-ratio",
        type=float,
        default=4.0,
        help="Assumed uncompressed:compressed ratio for FASTQ.gz size estimates.",
    )
    parser.add_argument("--spades-kmer-factor", type=float, default=8.7)
    parser.add_argument(
        "--min-spades-memory-gb",
        type=float,
        default=12.0,
        help="Conservative minimum SPAdes memory recommendation for metagenomic assemblies.",
    )
    parser.add_argument(
        "--available-memory-gb",
        type=float,
        default=None,
        help="Optional available memory to compare against the estimate.",
    )
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON.")
    args = parser.parse_args()

    base_dir = Path.cwd()
    samples = [row for row in read_tsv(args.samples) if row.get("data_type", "").upper() == "DNA"]
    assemblies = list(read_tsv(args.assemblies))

    sample_reports = []
    total_compressed_gib = 0.0
    total_uncompressed_gib = 0.0
    total_read_pairs = 0
    missing = []

    for row in samples:
        read_stats = []
        for col in ("read1", "read2"):
            path = resolve_path(row[col], base_dir)
            if not path.exists():
                missing.append(str(path))
                read_stats.append({"path": str(path), "missing": True})
                continue
            stats = gzip_sample_stats(path, args.sample_records, args.gzip_ratio)
            stats["path"] = str(path)
            stats["missing"] = False
            read_stats.append(stats)
            total_compressed_gib += stats["compressed_gib"]
            total_uncompressed_gib += stats["estimated_uncompressed_gib"]
        if len(read_stats) == 2 and not any(item.get("missing") for item in read_stats):
            pairs = min(read_stats[0]["estimated_reads"], read_stats[1]["estimated_reads"])
            total_read_pairs += pairs
        else:
            pairs = None
        sample_reports.append({"sample_id": row["sample_id"], "estimated_pairs": pairs, "reads": read_stats})

    assembly_reports = []
    sample_pairs = {item["sample_id"]: item["estimated_pairs"] or 0 for item in sample_reports}
    max_assembly_pairs = 0
    for row in assemblies:
        members = [item.strip() for item in row.get("sample_ids", "").split(",") if item.strip()]
        pairs = sum(sample_pairs.get(member, 0) for member in members)
        max_assembly_pairs = max(max_assembly_pairs, pairs)
        mem = (
            estimate_spades_memory_gib(
                pairs,
                args.mean_read_len,
                args.spades_kmer_factor,
                args.min_spades_memory_gb,
            )
            if pairs
            else None
        )
        assembly_reports.append(
            {
                "assembly_id": row.get("assembly_id", ""),
                "sample_ids": members,
                "estimated_pairs": pairs,
                "recommended_spades_memory_gib": mem,
                "recommended_spades_disk_gib": None if mem is None else max(50.0, total_uncompressed_gib * 8),
            }
        )

    max_spades_memory_gib = max([item["recommended_spades_memory_gib"] or 0 for item in assembly_reports] or [0])
    task_resources = estimate_task_resources(
        total_compressed_gib,
        total_uncompressed_gib,
        len(samples),
        len(assembly_reports),
        max_spades_memory_gib,
        args.min_spades_memory_gb,
    )
    notes = [
        "Estimates are rough preflight guardrails, not scheduler guarantees.",
        "SPAdes memory depends strongly on k-mer complexity; high-diversity metagenomes can exceed size-based estimates.",
        "GTDB-Tk, EggNOG, and CheckM reference database storage is not included in per-task work disk.",
        "Downstream estimates use input-size heuristics because contig count, MAG count, and database hits are not known before the run.",
    ]
    if not samples:
        notes.append("When no active samples are present, task estimates show conservative planning floors rather than a real run forecast.")

    report = {
        "samples": sample_reports,
        "assemblies": assembly_reports,
        "task_resources": task_resources["tasks"],
        "totals": {
            "dna_samples": len(samples),
            "assemblies": len(assemblies),
            "compressed_input_gib": total_compressed_gib,
            "estimated_uncompressed_input_gib": total_uncompressed_gib,
            "estimated_read_pairs": total_read_pairs,
            "missing_files": missing,
        },
        "recommendations": {
            "min_spades_memory_gib": max_spades_memory_gib or args.min_spades_memory_gb,
            "suggested_work_disk_gib": max(100.0, total_uncompressed_gib * 12),
            "suggested_outdir_disk_gib": max(50.0, total_uncompressed_gib * 4),
            "max_resources": task_resources["summary"],
            "available_memory_gib": args.available_memory_gb,
            "notes": notes,
        },
    }

    if args.json:
        print(json.dumps(report, indent=2))
        return

    print("Metagenomics Nextflow resource pre-estimate")
    print("=" * 45)
    print(f"DNA samples: {len(samples)}")
    print(f"Assemblies: {len(assemblies)}")
    print(f"Compressed DNA input: {fmt_gib(total_compressed_gib)}")
    print(f"Estimated uncompressed DNA input: {fmt_gib(total_uncompressed_gib)}")
    print(f"Estimated read pairs: {total_read_pairs:,}")
    if not samples:
        print("No active DNA sample rows found. Commented rows are ignored.")
    if missing:
        print("\nMissing input files:")
        for path in missing:
            print(f"  - {path}")

    print("\nPer-sample estimates:")
    for sample in sample_reports:
        pairs = sample["estimated_pairs"]
        pair_text = "unknown" if pairs is None else f"{pairs:,}"
        comp = sum(item.get("compressed_gib", 0) for item in sample["reads"])
        uncomp = sum(item.get("estimated_uncompressed_gib", 0) for item in sample["reads"])
        print(f"  - {sample['sample_id']}: {pair_text} pairs, {fmt_gib(comp)} compressed, {fmt_gib(uncomp)} uncompressed")

    print("\nPer-assembly SPAdes estimates:")
    for assembly in assembly_reports:
        mem = assembly["recommended_spades_memory_gib"]
        disk = assembly["recommended_spades_disk_gib"]
        print(
            f"  - {assembly['assembly_id']}: "
            f"{assembly['estimated_pairs']:,} pairs, "
            f"memory >= {fmt_gib(mem)}, work disk >= {fmt_gib(disk)}"
        )

    print("\nPer-task resource estimates:")
    print(f"  {'Task':<28} {'CPUs':>4} {'Memory':>10} {'Work disk':>10} {'Wall':>8}")
    print(f"  {'-' * 28} {'-' * 4:>4} {'-' * 10:>10} {'-' * 10:>10} {'-' * 8:>8}")
    for task in report["task_resources"]:
        print(
            f"  {task['task']:<28} "
            f"{task['cpus']:>4} "
            f"{fmt_gib(task['memory_gib']):>10} "
            f"{fmt_gib(task['work_disk_gib']):>10} "
            f"{fmt_hours(task['wall_hours']):>8}"
        )

    rec = report["recommendations"]
    max_rec = rec["max_resources"]
    print("\nOverall recommendations:")
    print(
        f"  Max memory: {fmt_gib(max_rec['max_memory_gib'])} "
        f"({max_rec['max_memory_task']})"
    )
    print(
        f"  Max work disk: {fmt_gib(max_rec['max_work_disk_gib'])} "
        f"({max_rec['max_work_disk_task']})"
    )
    print(f"  Max CPUs: {max_rec['max_cpus']} ({max_rec['max_cpus_task']})")
    print(
        f"  Longest single task: {fmt_hours(max_rec['max_wall_hours'])} "
        f"({max_rec['max_wall_time_task']})"
    )
    print(f"  Rough serial wall time: {fmt_hours(max_rec['serial_wall_hours_rough'])}")
    print(f"  Suggested total work disk headroom: >= {fmt_gib(rec['suggested_work_disk_gib'])}")
    print(f"  Suggested results disk headroom: >= {fmt_gib(rec['suggested_outdir_disk_gib'])}")
    if args.available_memory_gb is not None:
        status = "OK" if args.available_memory_gb >= max_rec["max_memory_gib"] else "LOW"
        print(f"  Available memory check: {args.available_memory_gb:.1f} GiB ({status})")
    print(
        "\nConclusion: plan for at least "
        f"{fmt_gib(max_rec['max_memory_gib'])} RAM and "
        f"{fmt_gib(max_rec['max_work_disk_gib'])} task work disk. "
        f"The peak-memory task is {max_rec['max_memory_task']}."
    )
    print("\nNotes:")
    for note in rec["notes"]:
        print(f"  - {note}")


if __name__ == "__main__":
    main()
