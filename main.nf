nextflow.enable.dsl = 2

/*
 * Manifest-driven metagenomics workflow ported from the Snakemake project.
 * RNA-specific steps are intentionally omitted.
 */

def readTsv(String path) {
    def lines = file(path).readLines().findAll { it.trim() && !it.trim().startsWith('#') }
    if (!lines) return []
    def header = lines.first().split('\t', -1)*.trim()
    lines.drop(1).collect { line ->
        def values = line.split('\t', -1)*.trim()
        [header, values].transpose().collectEntries { k, v -> [(k): v] }
    }
}

def tmpDir() {
    params.tmpdir ?: "${params.outdir}/tmp"
}

def drepDirName() {
    "drep_${params.drep_completeness}_${params.drep_contamination}"
        .replaceAll(/[^A-Za-z0-9]+/, '_')
        .replaceAll(/^_|_$/, '')
}

process VALIDATE_CONFIG {
    tag 'config'
    label 'python'
    publishDir 'config', mode: 'copy', pattern: 'validation.ok'
    publishDir 'logs/config', mode: 'copy', pattern: 'validate_config.log'

    input:
    path samples
    path assemblies
    path pacbio
    path scripts_dir

    output:
    path 'validation.ok'

    script:
    """
    python "${scripts_dir}/validate_config.py" \
      --samples ${samples} \
      --assemblies ${assemblies} \
      --pacbio-pools ${pacbio} \
      > validate_config.log 2>&1
    touch validation.ok
    """
}

process CHECK_DEPENDENCIES {
    tag 'diagnostics'
    label 'python'
    publishDir 'reports', mode: 'copy', pattern: 'dependency_check.tsv'
    publishDir 'logs/config', mode: 'copy', pattern: 'check_dependencies.log'

    input:
    path scripts_dir

    output:
    path 'dependency_check.tsv'

    script:
    """
    python "${scripts_dir}/check_dependencies.py" \
      --output dependency_check.tsv \
      > check_dependencies.log 2>&1
    """
}

process PREPARE_TRIMMOMATIC_ADAPTER {
    tag params.adapter
    label 'download'
    publishDir params.adapter_dir, mode: 'copy', pattern: params.adapter
    publishDir 'logs/resources', mode: 'copy', pattern: 'trimmomatic_adapters.log'

    input:
    path scripts_dir

    output:
    path params.adapter

    script:
    """
    python "${scripts_dir}/prepare_trimmomatic_adapters.py" \
      --adapter-dir adapters \
      --selected ${params.adapter} \
      > trimmomatic_adapters.log 2>&1
    cp adapters/${params.adapter} ${params.adapter}
    """
}

process TRIM_DNA {
    tag sample_id
    label 'trimmomatic'
    cpus { params.threads_trimmomatic as int }
    publishDir "${params.outdir}/dna/trim", mode: 'copy', pattern: '*.fq.gz'
    publishDir 'logs/trimmomatic/dna', mode: 'copy', pattern: '*.log'

    input:
    tuple val(sample_id), path(read1), path(read2)
    path adapter

    output:
    tuple val(sample_id), path("${sample_id}_F_paired.fq.gz"), path("${sample_id}_R_paired.fq.gz"), emit: trimmed

    script:
    """
    JAVA_TOOL_OPTIONS="${params.trimmomatic_java_opts}" trimmomatic PE -threads ${task.cpus} \
      ${read1} ${read2} \
      ${sample_id}_F_paired.fq.gz ${sample_id}_F_unpaired.fq.gz \
      ${sample_id}_R_paired.fq.gz ${sample_id}_R_unpaired.fq.gz \
      ILLUMINACLIP:${adapter}:${params.trimmomatic_illuminaclip} ${params.trimmomatic_steps} \
      > ${sample_id}.log 2>&1
    """
}

process FASTQC_DNA_TRIMMED {
    tag sample_id
    label 'fastqc'
    cpus { params.threads_fastqc as int }
    publishDir "${params.outdir}/dna/fastqc", mode: 'copy'
    publishDir 'logs/fastqc/dna', mode: 'copy', pattern: '*.log'

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    path "*_fastqc.html"
    path "*_fastqc.zip"

    script:
    """
    fastqc --threads ${task.cpus} --outdir . ${read1} ${read2} > ${sample_id}.log 2>&1
    """
}

process HYBRIDSPADES {
    tag assembly_id
    label 'spades'
    cpus { params.threads_assembly as int }
    memory { params.memory_assembly }
    publishDir "${params.outdir}/assembly", mode: 'copy'
    publishDir 'logs/spades', mode: 'copy', pattern: '*.log'

    input:
    tuple val(assembly_id), path(read1s), path(read2s), val(pacbio)

    output:
    tuple val(assembly_id), path("${assembly_id}/spades/contigs.fasta"), emit: contigs

    script:
    def pacbioArg = pacbio ? "--pacbio ${pacbio}" : ''
    def memGb = (task.memory.toGiga() as int)
    """
    mkdir -p ${assembly_id}/spades "${tmpDir()}"
    spades.py --meta \
      -1 ${read1s.join(',')} \
      -2 ${read2s.join(',')} \
      ${pacbioArg} \
      -o ${assembly_id}/spades \
      -k ${params.spades_kmers} \
      -t ${task.cpus} \
      -m ${memGb} \
      --tmp-dir "${tmpDir()}" \
      > ${assembly_id}.log 2>&1
    """
}

process FILTER_CONTIGS {
    tag assembly_id
    label 'python'
    publishDir "${params.outdir}/assembly", mode: 'copy'
    publishDir 'logs/filter_contigs', mode: 'copy', pattern: '*.log'

    input:
    tuple val(assembly_id), path(contigs)
    path scripts_dir

    output:
    tuple val(assembly_id), path("${assembly_id}/contigs_1000_hdmod.fasta"), emit: filtered

    script:
    """
    mkdir -p ${assembly_id}
    python "${scripts_dir}/filter_fasta_by_length.py" \
      --min-len ${params.min_contig_len} \
      --prefix ${assembly_id} \
      ${contigs} \
      ${assembly_id}/contigs_1000_hdmod.fasta \
      > ${assembly_id}.log 2>&1
    """
}

process METABAT {
    tag assembly_id
    label 'binning'
    cpus { params.threads_metabat as int }
    publishDir "${params.outdir}/mags/metabat", mode: 'copy'
    publishDir 'logs/metabat', mode: 'copy', pattern: '*.log'

    input:
    tuple val(assembly_id), path(contigs)

    output:
    tuple val(assembly_id), path("${assembly_id}"), emit: bins

    script:
    """
    mkdir -p ${assembly_id}
    metabat -i ${contigs} -o ${assembly_id}/${assembly_id}_metabat \
      -m ${params.metabat_min_contig} -t ${task.cpus} -v \
      > ${assembly_id}.log 2>&1
    """
}

process MAXBIN {
    tag assembly_id
    label 'maxbin'
    cpus { params.threads_maxbin as int }
    publishDir "${params.outdir}/mags/maxbin", mode: 'copy'
    publishDir 'logs/maxbin', mode: 'copy', pattern: '*.log'

    input:
    tuple val(assembly_id), path(contigs), path(read1), path(read2)

    output:
    tuple val(assembly_id), path("${assembly_id}"), emit: bins

    script:
    """
    mkdir -p ${assembly_id}
    run_MaxBin.pl \
      -contig ${contigs} \
      -out ${assembly_id}/${assembly_id}_maxbin \
      -reads ${read1} \
      -reads2 ${read2} \
      -thread ${task.cpus} \
      > ${assembly_id}.log 2>&1
    """
}

process SCAFFOLDS2BIN {
    tag "${assembly_id}:${label_suffix}"
    label 'python'
    publishDir "${params.outdir}/mags/scaffolds2bin", mode: 'copy', pattern: '*.tsv'
    publishDir 'logs/scaffolds2bin', mode: 'copy', pattern: '*.log'

    input:
    tuple val(assembly_id), path(bins_dir), val(label_suffix)
    path scripts_dir

    output:
    tuple val(assembly_id), val(label_suffix), path("${assembly_id}_${label_suffix}.scaffolds2bin.tsv"), emit: table

    script:
    """
    python "${scripts_dir}/fasta_bins_to_scaffolds2bin.py" \
      --bins-dir ${bins_dir} \
      --label-suffix _${label_suffix} \
      --output ${assembly_id}_${label_suffix}.scaffolds2bin.tsv \
      > ${assembly_id}.${label_suffix}.log 2>&1
    """
}

process DASTOOL {
    tag assembly_id
    label 'dastool'
    cpus { params.threads_dastool as int }
    publishDir "${params.outdir}/mags/dastool", mode: 'copy'
    publishDir 'logs/dastool', mode: 'copy', pattern: '*.log'

    input:
    tuple val(assembly_id), path(contigs), path(maxbin_table), path(metabat_table)

    output:
    path "${assembly_id}", emit: dastool_dir

    script:
    """
    mkdir -p ${assembly_id}
    DAS_Tool \
      -i ${maxbin_table},${metabat_table} \
      -c ${contigs} \
      -o ${assembly_id}/${assembly_id}_DAS \
      -l maxbin,metabat \
      --search_engine diamond \
      -t ${task.cpus} \
      --write_bin_evals \
      --write_bins \
      > ${assembly_id}.log 2>&1
    touch ${assembly_id}/.done
    """
}

process COLLECT_DASTOOL_BINS {
    label 'python'
    publishDir "${params.outdir}/mags", mode: 'copy', pattern: 'das_bins'
    publishDir 'logs/dastool', mode: 'copy', pattern: 'collect_bins.log'

    input:
    path dastool_dirs
    path scripts_dir

    output:
    path 'das_bins', emit: bins

    script:
    """
    python "${scripts_dir}/collect_dastool_bins.py" \
      --dastool-root . \
      --output das_bins \
      > collect_bins.log 2>&1
    """
}

process CHECKM {
    label 'checkm'
    cpus { params.threads_checkm as int }
    publishDir "${params.outdir}/mags/checkm", mode: 'copy'
    publishDir 'logs/checkm', mode: 'copy', pattern: '*.log'

    input:
    path bins

    output:
    path 'checkm', emit: checkm_dir

    script:
    """
    mkdir -p checkm
    checkm lineage_wf -t ${task.cpus} -x fa -f checkm/DAS_lineage.tsv ${bins} checkm \
      > lineage_wf.log 2>&1
    touch checkm/.done
    """
}

process DREP {
    label 'drep'
    cpus { params.threads_drep as int }
    memory { params.memory_drep }
    publishDir "${params.outdir}/mags/drep", mode: 'copy'
    publishDir 'logs/drep', mode: 'copy', pattern: '*.log'

    input:
    path bins
    val use_genome_info
    path genome_info
    val drep_out

    output:
    path "${drep_out}", emit: drep_dir

    script:
    def genomeInfoArg = use_genome_info ? "--genomeInfo ${genome_info}" : ''
    """
    dRep dereplicate ${drep_out} \
      -g ${bins}/*.fa \
      -comp ${params.drep_completeness} \
      -con ${params.drep_contamination} \
      -p ${task.cpus} \
      ${genomeInfoArg} \
      > dereplicate.log 2>&1

    status=\$?
    if [ "\$status" -ne 0 ]; then
      exit "\$status"
    fi
    if grep -q '!!! checkM failed !!!' dereplicate.log; then
      echo "dRep reported internal CheckM failure. See dereplicate.log and ${drep_out}/data/checkM/checkM_outdir/checkm.log." >&2
      exit 1
    fi
    if ! ls ${drep_out}/dereplicated_genomes/*.fa >/dev/null 2>&1; then
      echo "dRep produced no dereplicated genomes in ${drep_out}/dereplicated_genomes. Check dRep filters and logs." >&2
      exit 1
    fi
    touch ${drep_out}/.done
    """
}

process GTDBTK_CLASSIFY {
    label 'gtdbtk'
    cpus { params.threads_gtdbtk as int }
    publishDir "${params.outdir}/taxonomy/gtdbtk_classify", mode: 'copy'
    publishDir 'logs/gtdbtk', mode: 'copy', pattern: '*.log'

    input:
    path drep_dir

    output:
    path 'gtdbtk_classify', emit: classify_dir

    script:
    """
    export GTDBTK_DATA_PATH=${params.gtdbtk_data_path}
    gtdbtk classify_wf \
      --genome_dir ${drep_dir}/dereplicated_genomes \
      --extension fa \
      --out_dir gtdbtk_classify \
      --cpus ${task.cpus} \
      --mash_db gtdbtk_classify/mash_db \
      > classify.log 2>&1
    touch gtdbtk_classify/.done
    """
}

process GTDBTK_DE_NOVO {
    label 'gtdbtk'
    cpus { params.threads_gtdbtk as int }
    publishDir "${params.outdir}/taxonomy/gtdbtk_de_novo_bac", mode: 'copy'
    publishDir 'logs/gtdbtk', mode: 'copy', pattern: '*.log'

    input:
    path drep_dir
    path classify_dir

    output:
    path 'gtdbtk_de_novo_bac'

    script:
    """
    export GTDBTK_DATA_PATH=${params.gtdbtk_data_path}
    gtdbtk de_novo_wf \
      --genome_dir ${drep_dir}/dereplicated_genomes \
      --extension fa \
      --outgroup_taxon p__Thermotogota --bacteria \
      --out_dir gtdbtk_de_novo_bac \
      --cpus ${task.cpus} \
      --gtdbtk_classification_file ${classify_dir}/classify/gtdbtk.bac120.summary.tsv \
      > de_novo_bac.log 2>&1
    touch gtdbtk_de_novo_bac/.done
    """
}

process PRODIGAL_MAG {
    label 'prodigal'
    cpus { params.threads_prodigal as int }
    publishDir "${params.outdir}/annotation/mag_prodigal", mode: 'copy'
    publishDir 'logs/prodigal', mode: 'copy', pattern: '*.log'

    input:
    path drep_dir
    path scripts_dir

    output:
    path 'mag_prodigal', emit: mag_prodigal_dir

    script:
    """
    bash "${scripts_dir}/run_prodigal_dir.sh" \
      ${drep_dir}/dereplicated_genomes fa mag_prodigal mag.log
    """
}

process PRODIGAL_ASSEMBLY {
    tag assembly_id
    label 'prodigal'
    publishDir "${params.outdir}/annotation/assembly_prodigal", mode: 'copy'
    publishDir 'logs/prodigal/assembly', mode: 'copy', pattern: '*.log'

    input:
    tuple val(assembly_id), path(contigs)
    path scripts_dir

    output:
    tuple val(assembly_id), path("${assembly_id}.faa"), path("${assembly_id}.fna"), emit: genes

    script:
    """
    prodigal -p meta -i ${contigs} -a ${assembly_id}.faa -d ${assembly_id}.fna \
      > ${assembly_id}.log 2>&1
    python "${scripts_dir}/clean_prodigal_headers.py" ${assembly_id}.faa ${assembly_id}.fna
    """
}

process EGGNOG_DATABASE {
    label 'eggnog'
    publishDir params.eggnog_data_dir, mode: 'copy'
    publishDir 'logs/eggnog', mode: 'copy', pattern: '*.log'

    output:
    tuple path('eggnog'), path('eggnog/bac_arc.dmnd'), emit: db

    script:
    """
    mkdir -p eggnog
    export EGGNOG_DATA_DIR=eggnog
    download_eggnog_data.py --data_dir eggnog -y > database.log 2>&1
    create_dbs.py -m diamond --dbname bac_arc --taxa Bacteria,Archaea --data_dir eggnog >> database.log 2>&1
    """
}

process EGGNOG_ASSEMBLY {
    tag assembly_id
    label 'eggnog'
    cpus { params.threads_eggnog as int }
    publishDir "${params.outdir}/annotation/assembly_eggnog", mode: 'copy'
    publishDir 'logs/eggnog/assembly', mode: 'copy', pattern: '*.log'

    input:
    tuple val(assembly_id), path(faa), path(fna)
    tuple path(eggnog_data), path(dmnd)

    output:
    path "${assembly_id}", emit: eggnog_dir

    script:
    """
    mkdir -p ${assembly_id} "${tmpDir()}"
    export EGGNOG_DATA_DIR=${eggnog_data}
    emapper.py --cpu ${task.cpus} \
      --data_dir ${eggnog_data} \
      --dmnd_db ${dmnd} \
      --evalue ${params.eggnog_evalue} \
      -i ${faa} \
      -o ${assembly_id}/${assembly_id}_egg \
      --temp_dir "${tmpDir()}" \
      > ${assembly_id}.log 2>&1
    touch ${assembly_id}/.done
    """
}

process EGGNOG_MAG {
    label 'eggnog'
    cpus { params.threads_eggnog as int }
    publishDir "${params.outdir}/annotation/mag_eggnog", mode: 'copy'
    publishDir 'logs/eggnog', mode: 'copy', pattern: '*.log'

    input:
    path mag_prodigal_dir
    tuple path(eggnog_data), path(dmnd)
    path scripts_dir

    output:
    path 'mag_eggnog'

    script:
    """
    bash "${scripts_dir}/run_eggnog_dir.sh" \
      ${mag_prodigal_dir} mag_eggnog ${task.cpus} ${eggnog_data} ${dmnd} ${params.eggnog_evalue} "${tmpDir()}" mag.log
    """
}

process CONCATENATE_DREP_REFERENCE {
    label 'python'
    publishDir "${params.outdir}/instrain/reference", mode: 'copy'
    publishDir 'logs/instrain', mode: 'copy', pattern: '*.log'

    input:
    path drep_dir

    output:
    path 'dRep_SRGs.fa', emit: fasta

    script:
    """
    cat ${drep_dir}/dereplicated_genomes/*.fa > dRep_SRGs.fa 2> concatenate_reference.log
    """
}

process MAKE_STB {
    label 'python'
    publishDir "${params.outdir}/instrain/reference", mode: 'copy'
    publishDir 'logs/instrain', mode: 'copy', pattern: '*.log'

    input:
    path drep_dir
    path scripts_dir

    output:
    path 'dRep_SRGs.stb', emit: stb

    script:
    """
    python "${scripts_dir}/fasta_bins_to_stb.py" \
      --bins-dir ${drep_dir}/dereplicated_genomes \
      --output dRep_SRGs.stb \
      > make_stb.log 2>&1
    """
}

process BBMAP_DNA_FOR_INSTRAIN {
    tag sample_id
    label 'bbmap'
    cpus { params.threads_bbmap as int }
    publishDir "${params.outdir}/instrain/bam", mode: 'copy'
    publishDir 'logs/instrain/map', mode: 'copy', pattern: '*.log'

    input:
    tuple val(sample_id), path(read1), path(read2)
    path ref

    output:
    tuple val(sample_id), path("${sample_id}_dRep.sorted.bam"), path("${sample_id}_dRep.sorted.bam.bai"), emit: bam

    script:
    """
    bbmap.sh threads=${task.cpus} ref=${ref} in=${read1} in2=${read2} out=${sample_id}_dRep.bam nodisk=t \
      > ${sample_id}.log 2>&1
    samtools sort -@ ${task.cpus} -o ${sample_id}_dRep.sorted.bam ${sample_id}_dRep.bam
    samtools index ${sample_id}_dRep.sorted.bam
    """
}

process INSTRAIN_PROFILE {
    tag sample_id
    label 'instrain'
    cpus { params.threads_instrain_profile as int }
    publishDir "${params.outdir}/instrain/profile", mode: 'copy'
    publishDir 'logs/instrain/profile', mode: 'copy', pattern: '*.log'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path ref
    path stb
    path mag_prodigal_dir

    output:
    path "${sample_id}", emit: profile_dir

    script:
    """
    inStrain profile ${bam} ${ref} \
      -o ${sample_id} \
      -p ${task.cpus} \
      -g ${mag_prodigal_dir}/all_genes.fna \
      -s ${stb} \
      --database_mode \
      --skip_plot_generation \
      > ${sample_id}.log 2>&1
    touch ${sample_id}/.done
    """
}

process INSTRAIN_COMPARE {
    label 'instrain'
    cpus { params.threads_instrain_compare as int }
    publishDir "${params.outdir}/instrain/compare", mode: 'copy'
    publishDir 'logs/instrain', mode: 'copy', pattern: '*.log'

    input:
    path profile_dirs
    path stb

    output:
    path 'compare'

    script:
    """
    inStrain compare -i ${profile_dirs.join(' ')} \
      -s ${stb} \
      -p ${task.cpus} \
      -o compare \
      --database_mode \
      > compare.log 2>&1
    touch compare/.done
    """
}

workflow {
    scripts_dir = file('scripts')
    def dnaSamples = readTsv(params.samples).findAll { it.data_type?.toUpperCase() == 'DNA' }
    def assemblies = readTsv(params.assemblies)
    def pacbioPools = readTsv(params.pacbio_pools).collectEntries { [(it.pool_id): it.read_path] }
    def assemblyMeta = assemblies.collect { row ->
        def samples = row.sample_ids.split(',').collect { it.trim() }.findAll { it }
        def pacbio = row.pacbio_pool ? pacbioPools[row.pacbio_pool] : ''
        [assembly_id: row.assembly_id, samples: samples, pacbio: pacbio]
    }

    validate = VALIDATE_CONFIG(file(params.samples), file(params.assemblies), file(params.pacbio_pools), scripts_dir)
    CHECK_DEPENDENCIES(scripts_dir)
    adapter = file("${params.adapter_dir}/${params.adapter}").exists()
        ? Channel.value(file("${params.adapter_dir}/${params.adapter}"))
        : PREPARE_TRIMMOMATIC_ADAPTER(scripts_dir)

    dna_reads = Channel.fromList(dnaSamples).map { row ->
        tuple(row.sample_id, file(row.read1), file(row.read2))
    }

    trimmed = TRIM_DNA(dna_reads, adapter).trimmed
    FASTQC_DNA_TRIMMED(trimmed)

    trimmed_for_assembly = trimmed.collect(flat: false).map { trimList -> [trimList: trimList] }
    assembly_requests = Channel.fromList(assemblyMeta)
        .combine(trimmed_for_assembly)
        .map { meta, trimBundle ->
            def trimList = trimBundle.trimList
            if (trimList && !(trimList[0] instanceof List)) {
                trimList = [trimList]
            }
            def bySample = trimList.collectEntries { item -> [(item[0]): [item[1], item[2]]] }
            def read1s = meta.samples.collect { bySample[it][0] }
            def read2s = meta.samples.collect { bySample[it][1] }
            tuple(meta.assembly_id, read1s, read2s, meta.pacbio)
        }

    spades_contigs = HYBRIDSPADES(assembly_requests).contigs
    filtered_contigs = FILTER_CONTIGS(spades_contigs, scripts_dir).filtered

    if (!params.local_smoke) {
    metabat_bins = METABAT(filtered_contigs).bins.map { id, dir -> tuple(id, dir, 'metabat') }

    trimmed_for_maxbin = trimmed.collect(flat: false).map { trimList -> [trimList: trimList] }
    maxbin_inputs = filtered_contigs.combine(trimmed_for_maxbin).map { id, contigs, trimBundle ->
        def trimList = trimBundle.trimList
        if (trimList && !(trimList[0] instanceof List)) {
            trimList = [trimList]
        }
        def meta = assemblyMeta.find { it.assembly_id == id }
        def firstSample = meta.samples[0]
        def trim = trimList.find { it[0] == firstSample }
        tuple(id, contigs, trim[1], trim[2])
    }
    maxbin_bins = MAXBIN(maxbin_inputs).bins.map { id, dir -> tuple(id, dir, 'maxbin') }

    scaffolds_tables = SCAFFOLDS2BIN(metabat_bins.mix(maxbin_bins), scripts_dir).table
    scaffolds_tables_for_dastool = scaffolds_tables.collect(flat: false).map { tables -> [tables: tables] }
    dastool_inputs = filtered_contigs.combine(scaffolds_tables_for_dastool).map { id, contigs, tableBundle ->
        def tables = tableBundle.tables
        if (tables && !(tables[0] instanceof List)) {
            tables = [tables]
        }
        def metabat = tables.find { it[0] == id && it[1] == 'metabat' }
        def maxbin = tables.find { it[0] == id && it[1] == 'maxbin' }
        if (!metabat || !maxbin) {
            throw new IllegalStateException("Missing scaffold table for ${id}: metabat=${!!metabat}, maxbin=${!!maxbin}")
        }
        tuple(id, contigs, maxbin[2], metabat[2])
    }

    dastool_dirs = DASTOOL(dastool_inputs).dastool_dir.collect()
    das_bins = COLLECT_DASTOOL_BINS(dastool_dirs, scripts_dir).bins
    def genomeInfoFile = params.genome_info ? file(params.genome_info) : file('config/genome_info.csv')
    if (params.genome_info && !genomeInfoFile.exists()) {
        throw new IllegalArgumentException("Configured genome_info file does not exist: ${params.genome_info}")
    }
    def genomeInfoRows = genomeInfoFile.readLines().findAll { it.trim() && !it.trim().startsWith('#') }
    def useGenomeInfo = genomeInfoRows.size() > 1
    drep_dir = DREP(das_bins, useGenomeInfo, genomeInfoFile, drepDirName()).drep_dir

    if (params.run_gtdbtk) {
        classify_dir = GTDBTK_CLASSIFY(drep_dir).classify_dir
        GTDBTK_DE_NOVO(drep_dir, classify_dir)
    }

    mag_prodigal = PRODIGAL_MAG(drep_dir, scripts_dir).mag_prodigal_dir
    if (params.run_eggnog) {
        assembly_genes = PRODIGAL_ASSEMBLY(filtered_contigs, scripts_dir).genes
        def eggnog_db
        if (params.eggnog_data_dir && file(params.eggnog_data_dir).exists() && params.eggnog_dmnd_db && file(params.eggnog_dmnd_db).exists()) {
            eggnog_db = Channel.value(tuple(file(params.eggnog_data_dir), file(params.eggnog_dmnd_db)))
        } else {
            eggnog_db = EGGNOG_DATABASE().db
        }
        EGGNOG_ASSEMBLY(assembly_genes, eggnog_db)
        EGGNOG_MAG(mag_prodigal, eggnog_db, scripts_dir)
    }

    ref = CONCATENATE_DREP_REFERENCE(drep_dir).fasta
    stb = MAKE_STB(drep_dir, scripts_dir).stb
    bams = BBMAP_DNA_FOR_INSTRAIN(trimmed, ref).bam
    profiles = INSTRAIN_PROFILE(bams, ref, stb, mag_prodigal).profile_dir.collect(flat: false)
    compare_profiles = profiles.filter { profileDirs -> profileDirs.size() > 1 }
    INSTRAIN_COMPARE(compare_profiles, stb)
    }
}
