#!/usr/bin/env nextflow

// enable DSL 2 syntax
nextflow.enable.dsl = 2


// ----------------------------------------------------------------------------
// default parameter settings
// --------------------------

params.help                  = false
params.sampleSheet           = "samplesheet.csv"
params.fastqDir              = ""
params.sampleSize            = 100000
params.maxNumberToSampleFrom = 10000000000
params.chunkSize             = 1000000
params.trimStart             = 1
params.trimLength            = 36
params.genomeDetails         = "${projectDir}/resources/genomes.csv"
params.bowtieIndexDir        = "bowtie_indexes"
params.adaptersFasta         = "${projectDir}/resources/adapters.fa"
params.outputDir             = "${launchDir}"
params.outputPrefix          = ""


// ----------------------------------------------------------------------------
// summary of configuration parameters
// -----------------------------------

log.info ""
log.info """\
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
"""
log.info ""


// ----------------------------------------------------------------------------
// help/usage
// ----------

if (params.help) {
    log.info ""
    log.info """\
Usage:
    nextflow run crukci-bioinformatics/mga2

Options:
    --help                        Show this message and exit
    --sample-sheet                Sample sheet CSV file containing id, fastq (file path/pattern) and species columns
    --fastq-dir                   Directory in which FASTQ files are located (optional, can specify absolute or relative paths in sample sheet instead)
    --sample-size                 Number of sequences to sample for each sample/dataset
    --max-number-to-sample-from   Maximum number of sequences to read/sample from
    --chunk-size                  Number of sequences for each chunk in batch processing of sampled sequences
    --trim-start                  The position at which the trimmed sequence starts, all bases before this position are trimmed
    --trim-length                 The length of the trimmed sequences
    --genome-details              Genome details CSV file containing genome, species and synonym columns
    --bowtie-index-dir            Directory containing bowtie indexes for reference genomes
    --adapters-fasta              FASTA file containing adapter sequences
    --output-dir                  Directory to which output files are written
    --output-prefix               Prefix for output file names

Alternatively, override settings using a configuration file such as the
following, in which parameter names used are the camelCase equivalent of the
above options:

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

Usage:
    nextflow run crukci-bioinformatics/mga2 -c mga2.config
    """
    log.info ""
    exit 0
}


// ----------------------------------------------------------------------------
// validate input parameters
// -------------------------

fastqDir = "${params.fastqDir}"
if (!"${fastqDir}".isEmpty() && !"${fastqDir}".endsWith("/")) {
    fastqDir = "${fastqDir}/"
}

if (!"${params.sampleSize}".isInteger() || "${params.sampleSize}" as Integer < 1000) {
    exit 1, "Invalid sample size - set to at least 1000 (100000 recommended)"
}

if (!"${params.maxNumberToSampleFrom}".isLong() || "${params.maxNumberToSampleFrom}" as Long < "${params.sampleSize}" as Long) {
    exit 1, "Invalid number of sequence reads to sample from (${params.maxNumberToSampleFrom}) - must be at least as large as the sample size (${params.sampleSize})"
}

if (!"${params.chunkSize}".isInteger() || "${params.chunkSize}" as Integer < 100000) {
    exit 1, "Invalid chunk size for batch alignment - set to at least 100000 (1000000 recommended)"
}

trimStart = 1

if ("${params.trimStart}".isInteger()) {
    trimStart = "${params.trimStart}" as Integer
} else {
    exit 1, "Invalid trim start"
}

if (trimStart <= 0) {
    exit 1, "Invalid trim start - must be a positive integer value"
}

trimLength = 36

if ("${params.trimLength}".isInteger()) {
    trimLength = "${params.trimLength}" as Integer
} else {
    exit 1, "Invalid trim length"
}

if (trimLength < 30) {
    exit 1, "Invalid trim length - trimmed sequences should be sufficiently long to align to reference genomes, i.e. at least 30"
}

// calculate minimum sequence length used for sampling sequences
minimumSequenceLength = trimStart + trimLength - 1


// ----------------------------------------------------------------------------
// processes
// ---------

process check_inputs {
    input:
        path install_dir
        path samples
        path genome_details
        path bowtie_index_list

    output:
        path "samples.checked.csv", emit: samples
        path "genomes.checked.csv", emit: genomes

    script:
        """
        Rscript ${install_dir}/R/check_inputs.R ${samples} ${genome_details} ${bowtie_index_list} samples.checked.csv genomes.checked.csv
        """
}


process sample_fastq {
    tag "${id} ${name}"

    input:
        tuple val(id), val(name), path(fastq), val(fastq_pattern)

    output:
        path "${id}.sample.fq", emit: fastq
        path "${id}.summary.csv", emit: summary

    script:
        """
        RUST_LOG=info \
        sample-fastq \
            --id="${id}" \
            --sample-size=${params.sampleSize} \
            --max-number-to-sample-from=${params.maxNumberToSampleFrom} \
            --min-sequence-length=${minimumSequenceLength} \
            --prepend-id \
            --output-file="${id}.sample.fq" \
            --summary-file="${id}.summary.csv" \
            --check-unique-record-ids \
            ${fastq}
        """
}


process trim_and_split {
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


process bowtie {
    tag "${prefix}.${genome}"

    input:
        each path(fastq)
        path bowtie_index_dir
        each genome

    output:
        path "${prefix}.${genome}.tsv"

    script:
        prefix=fastq.baseName
        """
        set -eo pipefail
        echo "genome	read	strand	chromosome	position	sequence	quality	num	mismatches" > ${prefix}.${genome}.tsv
        if [[ `head ${fastq} | wc -l` -gt 0 ]]
        then
            bowtie \
                --time \
                --best \
                --chunkmbs 256 \
                -x ${bowtie_index_dir}/${genome} \
                ${fastq} \
            | sed "s/^/${genome}\t/" \
            >> "${prefix}.${genome}.tsv"
        fi
        """
}


process exonerate {
    tag "${prefix}.adapters"

    input:
        each path(fasta)
        path adapters_fasta

    output:
        path "${prefix}.adapter_alignments.tsv"

    script:
        prefix=fasta.baseName
        """
        echo "read	start	end	strand	adapter	adapter start	adapter end	adapter strand	percent identity	score" > ${prefix}.adapter_alignments.tsv
        if [[ `head ${fasta} | wc -l` -gt 0 ]]
        then
            exonerate \
                --model ungapped \
                --showalignment no \
                --showvulgar no \
                --verbose 0 \
                --bestn 1 \
                --ryo "%qi\t%qab\t%qae\t%qS\t%ti\t%tab\t%tae\t%tS\t%pi\t%s\n" \
                ${fasta} \
                ${adapters_fasta} \
            >> "${prefix}.adapter_alignments.tsv"
        fi
        """
}


process create_summary {
    publishDir "${params.outputDir}", mode: 'copy'

    input:
        path install_dir
        path samples
        path genomes
        path counts
        path alignments
        path adapter_alignments

    output:
        path "${params.outputPrefix}mga_summary.csv"
        path "${params.outputPrefix}mga_alignment_summary.csv"
        path "${params.outputPrefix}mga_genome_alignments.tsv.gz"
        path "${params.outputPrefix}mga_adapter_alignments.tsv.gz"
        path "${params.outputPrefix}mga_alignment_bar_chart.png"
        path "${params.outputPrefix}mga_alignment_bar_chart.pdf"
        path "${params.outputPrefix}mga_alignment_bar_chart.svg"

    script:
        """
        Rscript ${install_dir}/R/summarize_alignments.R \
            --samples=${samples} \
            --genomes=${genomes} \
            --counts=${counts} \
            --alignments=${alignments} \
            --adapter-alignments=${adapter_alignments} \
            --output-prefix="${params.outputPrefix}"
        """
}


// ----------------------------------------------------------------------------
// workflow
// --------

workflow {

    install_dir = channel.fromPath(projectDir, checkIfExists: true)

    samples = channel.fromPath(params.sampleSheet, checkIfExists: true)

    genome_details = channel.fromPath(params.genomeDetails, checkIfExists: true)

    bowtie_index_dir = channel.fromPath("${params.bowtieIndexDir}", checkIfExists: true)

    genomes = channel
        .fromPath("${params.bowtieIndexDir}/*.rev.1.ebwt", checkIfExists: true)
        .map { "${it.name}".replaceFirst(/.rev.1.ebwt$/, "") }

    bowtie_index_list = genomes.collectFile(name: "bowtie_index_list.txt", newLine: true)

    check_inputs(
        install_dir,
        samples,
        genome_details,
        bowtie_index_list
    )

    adapters_fasta = channel.fromPath("${params.adaptersFasta}", checkIfExists: true)

    fastq = check_inputs.out.samples
        .splitCsv(header: true, quote: '"')
        .map { row -> tuple("${row.id}", "${row.name}", file("${fastqDir}${row.fastq}", checkIfExists: true), "${fastqDir}${row.fastq}") }

    fastq.subscribe { assert !it[2].isEmpty(), "No FASTQ files found for ${it[1]} matching pattern ${it[3]}" }

    sample_fastq(fastq)

    counts = sample_fastq.out.summary
        .collectFile(name: "sequence_counts.collected.csv", keepHeader: true)

    trim_and_split(sample_fastq.out.fastq.collect())

    bowtie(
        trim_and_split.out.fastq,
        bowtie_index_dir,
        genomes
    )

    alignments = bowtie.out.collectFile(name: "alignments.collected.tsv", keepHeader: true)

    exonerate(
        trim_and_split.out.fasta,
        adapters_fasta
    )

    adapter_alignments = exonerate.out.collectFile(name: "adapter_alignments.collected.tsv", keepHeader: true)

    create_summary(
        install_dir,
        check_inputs.out.samples,
        check_inputs.out.genomes,
        counts,
        alignments,
        adapter_alignments
    )
}

