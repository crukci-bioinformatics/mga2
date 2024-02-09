#!/usr/bin/env nextflow

// enable DSL 2 syntax
nextflow.enable.dsl = 2

include { check_inputs; sample_fastq; trim_and_split; bowtie; split_genome_alignments_by_sample; exonerate; split_adapter_alignments_by_sample; summarize_alignments; compress_alignments; create_bar_chart } from "./processes"


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


// -----------------------------------------------------------------------------
// workflow
// -----------------------------------------------------------------------------

workflow {
    sample_sheet = channel.fromPath(params.sampleSheet, checkIfExists: true)
    genome_details = channel.fromPath(params.genomeDetails, checkIfExists: true)
    bowtie_index_dir = channel.fromPath(params.bowtieIndexDir, checkIfExists: true)
    adapters_fasta = channel.fromPath(params.adaptersFasta, checkIfExists: true)

    genomes = channel
        .fromPath("${params.bowtieIndexDir}/*.rev.1.ebwt{,l}", checkIfExists: true)
        .map { "${it.name}".replaceFirst(/.rev.1.ebwtl?$/, "") }

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

    // calculate minimum sequence length used for sampling sequences
    minimumSequenceLength = params.trimStart + params.trimLength - 1

    sample_fastq(fastq.map { it[0..2] }, params.sampleSize, params.maxNumberToSampleFrom, minimumSequenceLength)

    counts = sample_fastq.out.summary
        .collectFile(name: "sampling_summary.csv", keepHeader: true)

    trim_and_split(sample_fastq.out.fastq.collect(), params.chunkSize, params.trimStart, params.trimLength)

    fasta = trim_and_split.out.fasta
        .flatten()
        .map { it -> tuple(it.name.split("\\.")[1].toInteger(), it) }

    bowtie(
        trim_and_split.out.fastq,
        bowtie_index_dir,
        genomes.collect()
    )

    chunk_genome_alignments = bowtie.out
        .collectFile(keepHeader: true) { it -> [ "chunk.${it.name.split('\\.')[1]}.genome_alignments.tsv", it ] }
        .map { it -> tuple(it.name.split("\\.")[1].toInteger(), it) }

    sample_genome_alignments = split_genome_alignments_by_sample(chunk_genome_alignments.join(fasta))
        .flatten()
        .collectFile(keepHeader: true) { it -> [ "sample.${it.name.split('\\.')[1]}.genome_alignments.tsv", it ] }
        .map { it -> tuple(it.name.split("\\.")[1].toInteger(), it) }

    exonerate(
        trim_and_split.out.fasta,
        adapters_fasta
    )

    chunk_adapter_alignments = exonerate.out
        .collectFile(keepHeader: true) { it -> [ "chunk.${it.name.split('\\.')[1]}.adapter_alignments.tsv", it ] }
        .map { it -> tuple(it.name.split("\\.")[1].toInteger(), it) }

    sample_adapter_alignments = split_adapter_alignments_by_sample(chunk_adapter_alignments.join(fasta))
        .flatten()
        .collectFile(keepHeader: true) { it -> [ "sample.${it.name.split('\\.')[1]}.adapter_alignments.tsv", it ] }
        .map { it -> tuple(it.name.split("\\.")[1].toInteger(), it) }

    sample_fastq.out.summary
        .map { it -> tuple(it.name.split("\\.")[1].toInteger(), it) }
        .join(sample_genome_alignments)
        .join(sample_adapter_alignments)
        .combine(check_inputs.out.samples)
        .combine(check_inputs.out.genomes)
        | summarize_alignments

    summary = summarize_alignments.out.summary
        .collectFile(
            name: "${params.outputPrefix}summary.csv",
            storeDir: "${params.outputDir}",
            keepHeader: true,
            sort: { it.name.split("\\.")[1].toInteger() }
        )

    alignment_summary = summarize_alignments.out.alignment_summary
        .collectFile(
            name: "${params.outputPrefix}alignment_summary.csv",
            storeDir: "${params.outputDir}",
            keepHeader: true,
            sort: { it.name.split("\\.")[1].toInteger() }
        )

    genome_alignments = summarize_alignments.out.genome_alignments
        .collectFile(
            name: "${params.outputPrefix}genome_alignments.tsv",
            keepHeader: true,
            sort: { it.name.split("\\.")[1].toInteger() }
        )

    adapter_alignments = summarize_alignments.out.adapter_alignments
        .collectFile(
            name: "${params.outputPrefix}adapter_alignments.tsv",
            keepHeader: true,
            sort: { it.name.split("\\.")[1].toInteger() }
        )

    compress_alignments(genome_alignments, adapter_alignments, params.outputDir)

    create_bar_chart(summary, alignment_summary, params.outputDir, params.outputPrefix)
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
