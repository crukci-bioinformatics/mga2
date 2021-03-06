#!/usr/bin/env nextflow

// enable DSL 2 syntax
nextflow.enable.dsl = 2


// -----------------------------------------------------------------------------
// show settings and/or help
// -----------------------------------------------------------------------------

if (params.showSettings) {
    printParameterSummary()
}

if (params.help) {
    helpMessage()
    exit 0
}


// -----------------------------------------------------------------------------
// check/derive parameters
// -----------------------------------------------------------------------------

fastqDir = params.fastqDir
if (fastqDir && !fastqDir.endsWith("/")) {
    fastqDir = "${params.fastqDir}/"
}

if (params.sampleSize < 1000) {
    exit 1, "Invalid sample size - set to at least 1000 (100000 recommended)"
}

if (params.maxNumberToSampleFrom < params.sampleSize) {
    exit 1, "Invalid maximum number of sequence reads to sample from - must be at least as large as the sample size"
}

if (params.chunkSize < 100000) {
    exit 1, "Invalid chunk size for batch alignment - set to at least 100000 (1000000 recommended)"
}

if (params.trimStart <= 0) {
    exit 1, "Invalid trim start - must be a positive integer value"
}

if (params.trimLength < 30) {
    exit 1, "Invalid trim length - trimmed sequences should be sufficiently long to align to reference genomes, i.e. at least 30"
}

// calculate minimum sequence length used for sampling sequences
minimumSequenceLength = params.trimStart + params.trimLength - 1


// -----------------------------------------------------------------------------
// processes
// -----------------------------------------------------------------------------

// check input sample sheet and details of genomes to be screened
process check_inputs {
    label 'mga2'

    executor 'local'

    input:
        path samples
        path genome_details
        path bowtie_index_list

    output:
        path checked_samples, emit: samples
        path checked_genomes, emit: genomes

    script:
        checked_samples = "samples.checked.csv"
        checked_genomes = "genomes.checked.csv"
        """
        check_inputs.R \
            --samples=${samples} \
            --genomes=${genome_details} \
            --bowtie-indexes=${bowtie_index_list} \
            --output-samples=${checked_samples} \
            --output-genomes=${checked_genomes}
        """
}

// sample records from input FASTQ file(s)
process sample_fastq {
    label 'mga2'
    tag "${id} ${name}"

    time 12.hour

    input:
        tuple val(id), val(name), path(fastq)

    output:
        path sampled_fastq, emit: fastq
        path summary, emit: summary

    script:
        sampled_fastq = "${id}.sample.fq"
        summary = "${id}.summary.csv"
        """
        RUST_LOG=info \
        sample-fastq \
            --id="${id}" \
            --sample-size=${params.sampleSize} \
            --max-number-to-sample-from=${params.maxNumberToSampleFrom} \
            --min-sequence-length=${minimumSequenceLength} \
            --prepend-id \
            --output-file=${sampled_fastq} \
            --summary-file=${summary} \
            --check-unique-record-ids \
            ${fastq}
        """
}


// trim sequences and split into chunks
process trim_and_split {
    label 'mga2'

    executor 'local'

    input:
        path fastq

    output:
        path "chunk.*.fq", emit: fastq
        path "chunk.*.fa", emit: fasta

    script:
        """
        RUST_LOG=info \
        trim-and-split-fastq \
            --chunk-size=${params.chunkSize} \
            --output-prefix=chunk \
            --start=${params.trimStart} \
            --length=${params.trimLength} \
            --output-fasta \
            ${fastq}
        """
}


// align trimmed sequences in a chunk file against a reference genome using bowtie
process bowtie {
    label 'mga2'
    tag "${prefix}.${genome}"

    memory { 4.GB * 2 ** (task.attempt - 1) }
    time { 2.hour * 2 ** (task.attempt - 1) }
    maxRetries 2

    input:
        each path(fastq)
        path bowtie_index_dir
        each genome

    output:
        path alignments

    script:
        prefix = fastq.baseName
        alignments = "${prefix}.${genome}.tsv"
        """
        set -o pipefail
        echo -e "genome\\tread\\tstrand\\tchromosome\\tposition\\tsequence\\tquality\\tnum\\tmismatches" > ${prefix}.${genome}.tsv
        if [[ `head ${fastq} | wc -l` -gt 0 ]]
        then
            bowtie \
                --time \
                --best \
                --chunkmbs 256 \
                -x ${bowtie_index_dir}/${genome} \
                ${fastq} \
            | sed "s/^/${genome}\t/" \
            >> ${alignments}
        fi
        """
}


// align untrimmed sequences in a chunk file against adapter sequences using exonerate
process exonerate {
    label 'mga2'
    tag "${prefix}.adapters"

    memory { 2.GB * task.attempt }
    time { 1.hour * task.attempt }
    maxRetries 2

    input:
        each path(fasta)
        path adapters_fasta

    output:
        path alignments

    script:
        prefix = fasta.baseName
        alignments = "${prefix}.adapter_alignments.tsv"
        """
        echo -e "read\\tstart\\tend\\tstrand\\tadapter\\tadapter start\\tadapter end\\tadapter strand\\tpercent identity\\tscore" > ${alignments}
        if [[ `head ${fasta} | wc -l` -gt 0 ]]
        then
            exonerate \
                --model ungapped \
                --showalignment no \
                --showvulgar no \
                --verbose 0 \
                --bestn 1 \
                --ryo "%qi\\t%qab\\t%qae\\t%qS\\t%ti\\t%tab\\t%tae\\t%tS\\t%pi\\t%s\\n" \
                ${fasta} \
                ${adapters_fasta} \
            >> ${alignments}
        fi
        """
}


// create summary tables and stacked bar chart
process create_summary {
    label 'mga2'

    publishDir "${params.outputDir}", mode: 'copy'

    memory { 2.GB * task.attempt }
    time { 1.hour * task.attempt }
    maxRetries 2

    input:
        path samples
        path genomes
        path counts
        path alignments
        path adapter_alignments

    output:
        path summary, emit: summary
        path alignment_summary, emit: alignment_summary
        path "${params.outputPrefix}mga_alignment_summary.png"
        path "${params.outputPrefix}mga_alignment_summary.pdf"
        path "${params.outputPrefix}mga_alignment_summary.svg"
        path "${params.outputPrefix}mga_genome_alignments.tsv.gz"
        path "${params.outputPrefix}mga_adapter_alignments.tsv.gz"

    script:
        summary = "${params.outputPrefix}mga_summary.csv"
        alignment_summary = "${params.outputPrefix}mga_alignment_summary.csv"
        """
        summarize_alignments.R \
            --samples=${samples} \
            --genomes=${genomes} \
            --counts=${counts} \
            --alignments=${alignments} \
            --adapter-alignments=${adapter_alignments} \
            --output-prefix="${params.outputPrefix}"

        create_bar_chart.R \
            --summary=${summary} \
            --alignment-summary=${alignment_summary} \
            --output-prefix="${params.outputPrefix}"
        """
}


// -----------------------------------------------------------------------------
// workflow
// -----------------------------------------------------------------------------

workflow {

    sample_sheet = channel.fromPath(params.sampleSheet, checkIfExists: true)
    genome_details = channel.fromPath(params.genomeDetails, checkIfExists: true)
    bowtie_index_dir = channel.fromPath(params.bowtieIndexDir, checkIfExists: true)
    adapters_fasta = channel.fromPath(params.adaptersFasta, checkIfExists: true)

    genomes = channel
        .fromPath("${params.bowtieIndexDir}/*.rev.1.ebwt", checkIfExists: true)
        .map { "${it.name}".replaceFirst(/.rev.1.ebwt$/, "") }

    bowtie_index_list = genomes.collectFile(name: "bowtie_index_list.txt", newLine: true)

    check_inputs(
        sample_sheet,
        genome_details,
        bowtie_index_list
    )

    fastq = check_inputs.out.samples
        .splitCsv(header: true, quote: '"')
        .map { row -> tuple("${row.id}", "${row.name}", file("${fastqDir}${row.fastq}", checkIfExists: true), "${fastqDir}${row.fastq}") }

    fastq.subscribe { assert !it[2].isEmpty(), "No FASTQ files found for ${it[1]} matching pattern ${it[3]}" }

    sample_fastq(fastq.map { it[0..2] })

    counts = sample_fastq.out.summary
        .collectFile(name: "sequence_counts.collected.csv", keepHeader: true)

    trim_and_split(sample_fastq.out.fastq.collect())

    bowtie(
        trim_and_split.out.fastq,
        bowtie_index_dir,
        genomes.collect()
    )

    alignments = bowtie.out.collectFile(name: "alignments.collected.tsv", keepHeader: true)

    exonerate(
        trim_and_split.out.fasta,
        adapters_fasta
    )

    adapter_alignments = exonerate.out.collectFile(name: "adapter_alignments.collected.tsv", keepHeader: true)

    create_summary(
        check_inputs.out.samples,
        check_inputs.out.genomes,
        counts,
        alignments,
        adapter_alignments
    )
}


// -----------------------------------------------------------------------------
// summary of configuration parameters
// -----------------------------------------------------------------------------

def printParameterSummary() {
    log.info ""
    log.info """
        Multi-Genome Alignment (MGA) Contaminant Screen
        ===============================================

        Sample sheet           : ${params.sampleSheet}
        FASTQ directory        : ${params.fastqDir}
        Sample size            : ${params.sampleSize}
        Maximum sampled from   : ${params.maxNumberToSampleFrom}
        Chunk size             : ${params.chunkSize}
        Trim start             : ${params.trimStart}
        Trim length            : ${params.trimLength}
        Genomes details file   : ${params.genomeDetails}
        Bowtie index directory : ${params.bowtieIndexDir}
        Adapters FASTA file    : ${params.adaptersFasta}
        Output directory       : ${params.outputDir}
        Output prefix          : ${params.outputPrefix}
    """.stripIndent()
    log.info ""
}


// ----------------------------------------------------------------------------
// help/usage
// ----------

def helpMessage() {
    log.info """
        Usage:
            nextflow run crukci-bioinformatics/mga2

        Options:
            --help                        Show this message and exit
            --sampleSheet                 CSV file containing details of sample dataset (id, fastq, species and control columns required)
            --fastqDir                    Directory in which FASTQ files are located (optional, can specify absolute or relative paths in sample sheet instead)
            --sampleSize                  Number of sequences to sample for each sample/dataset
            --maxNumberToSampleFrom       Maximum number of sequences to read/sample from
            --chunkSize                   Number of sequences for each chunk for batch alignment of sampled sequences
            --trimStart                   The position at which the trimmed sequence starts, all bases before this position are trimmed
            --trimLength                  The length of the trimmed sequences
            --genomeDetails               CSV file containing the species name and synonyms for each reference genome (genome, species and synonym colums required)
            --bowtieIndexDir              Directory containing bowtie indexes for reference genomes
            --adaptersFasta               FASTA file containing adapter sequences
            --outputDir                   Directory to which output files are written
            --outputPrefix                Prefix for output file names

        Alternatively, override settings using a configuration file such as the
        following:

        params {
            sampleSheet    = "samplesheet.csv"
            sampleSize     = 100000
            trimStart      = 11
            trimLength     = 36
            genomeDetails  = "genomes.csv"
            bowtieIndexDir = "/path_to/bowtie_indexes"
            outputDir      = "mga"
            outputPrefix   = ""
        }

        and run as follows:
            nextflow run crukci-bioinformatics/mga2 -c mga2.config
    """.stripIndent()
    log.info ""
}
