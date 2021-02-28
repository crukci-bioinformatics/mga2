#!/usr/bin/env nextflow

// default parameter settings

def defaults = [
    sampleSheet: 'samplesheet.csv',
    genomeDetails: "${projectDir}/resources/genomes.csv",
    fastqDir: "",
    sampleSize: 100000,
    maxNumberToSampleFrom: 10000000000,
    chunkSize: 1000000,
    trimStart: 11,
    trimLength: 36,
    bowtieIndexDir: "bowtie_indexes",
    adaptersFasta: "${projectDir}/resources/adapters.fa",
    outputDir: "${launchDir}",
    outputPrefix: ""
]

// set paramters to default settings

params.help = false
params.sampleSheet = defaults.sampleSheet
params.genomeDetails = defaults.genomeDetails
params.fastqDir = defaults.fastqDir
params.sampleSize = defaults.sampleSize
params.maxNumberToSampleFrom = defaults.maxNumberToSampleFrom
params.chunkSize = defaults.chunkSize
params.trimStart = defaults.trimStart
params.trimLength = defaults.trimLength
params.bowtieIndexDir = defaults.bowtieIndexDir
params.adaptersFasta = defaults.adaptersFasta
params.outputDir = defaults.outputDir
params.outputPrefix = defaults.outputPrefix

//print usage

if (params.help) {
    log.info ''
    log.info """\
Multi-Genome Alignment (MGA) Contaminant Screen
===============================================

Usage:
    nextflow run crukcibioinformatics/mga

Options:
    --help                            Show this message and exit
    --sample-sheet FILE               Sample sheet CSV file containing id, fastq (file path/pattern) and species columns (default: ${defaults.sampleSheet})
    --genome-details FILE             Genome details CSV files containing genome, species and synonym columns (default: ${defaults.genomeDetails})
    --fastq-dir DIR                   Directory in which FASTQ files are located (optional, can specify absolute or relative paths in sample sheet instead)
    --sample-size INTEGER             Number of sequences to sample for each sample/dataset (default: ${defaults.sampleSize})
    --chunk-size INTEGER              Number of sequences for each chunk in batch processing of sampled sequences (default: ${defaults.chunkSize})
    --trim-start INTEGER              The position at which the trimmed sequence starts, all bases before this position are trimmed (default: ${defaults.trimStart})
    --trim-length INTEGER             The maximum length of the trimmed sequence (default: ${defaults.trimLength})
    --bowtie-index-dir PATH           Directory containing bowtie indexes for reference genomes (default: ${defaults.bowtieIndexDir})
    --adapters-fasta FILE             FASTA file containing adapter sequences (default: ${defaults.adaptersFasta})
    --output-dir PATH                 Output directory (default: ${defaults.outputDir})
    --output-prefix PREFIX            Prefix for output files (default: ${defaults.outputPrefix})
    """
    log.info ''
    exit 1
}

// summary of configuration parameters

log.info ''
log.info """\
Multi-Genome Alignment (MGA) Contaminant Screen
===============================================

Sample sheet           : $params.sampleSheet
Sample size            : $params.sampleSize
FASTQ directory        : $params.fastqDir
Chunk size             : $params.chunkSize
Trim start             : $params.trimStart
Trim length            : $params.trimLength
Bowtie index directory : $params.bowtieIndexDir
Adapters FASTA file    : $params.adaptersFasta
Output directory       : $params.outputDir
Output prefix          : $params.outputPrefix
"""
log.info ''

// validate input parameters and calculate minimum sequence length used for sampling sequences

fastqDir = "${params.fastqDir}"
if (!"${fastqDir}".isEmpty() && !"${fastqDir}".endsWith("/")) {
    fastqDir = "${fastqDir}/"
}

/*
if (!"${params.sampleSize}".isInteger() || "${params.sampleSize}" as Integer < 100000) {
    log.error 'Invalid sample size - set to at least the recommended value of 100000'
    exit 1
}
*/

if (!"${params.chunkSize}".isInteger() || "${params.chunkSize}" as Integer < 100000) {
    log.error 'Invalid chunk size for batch alignment - set to at least 100000 (recommend 5000000)'
    exit 1
}

int trimStart = 1

if ("${params.trimStart}".isInteger()) {
    trimStart = "${params.trimStart}" as Integer
} else {
    log.error 'Invalid trim start'
    exit 1
}

if (trimStart <= 0) {
    log.error 'Invalid trim start - must be a positive integer value'
    exit 1
}

int trimLength = 36

if ("${params.trimLength}".isInteger()) {
    trimLength = "${params.trimLength}" as Integer
} else {
    log.error 'Invalid trim length'
    exit 1
}

if (trimLength < 30) {
    log.error 'Invalid trim length - trimmed sequences should be sufficiently long to align to reference genomes, i.e. at least 30'
    exit 1
}

int minimumSequenceLength = trimStart + trimLength - 1


// enable DSL 2 syntax
nextflow.enable.dsl = 2


process check_inputs {
    input:
        path samples
        path genome_details
        path bowtie_index_list

    output:
        path 'samples.checked.csv', emit: samples
        path 'genomes.checked.csv', emit: genomes

    script:
        """
        Rscript ${projectDir}/R/check_inputs.R ${samples} ${genome_details} ${bowtie_index_list} samples.checked.csv genomes.checked.csv
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
        path 'chunk.*.fq', emit: fastq
        path 'chunk.*.fa', emit: fasta

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
        path "${prefix}.${genome}.txt"

    script:
        prefix=fastq.baseName
        """
        set -eo pipefail
        echo "genome	read	strand	chromosome	position	sequence	quality	num	mismatches" > ${prefix}.${genome}.txt
        if [[ `head ${fastq} | wc -l` -gt 0 ]]
        then
            bowtie \
                --time \
                --best \
                --chunkmbs 256 \
                -x ${bowtie_index_dir}/${genome} \
                ${fastq} \
            | sed "s/^/${genome}\t/" \
            >> "${prefix}.${genome}.txt"
        fi
        """
}


process exonerate {
    tag "${prefix}.adapters"

    input:
        each path(fasta)
        path adapters_fasta

    output:
        path "${prefix}.adapter_alignments.txt"

    script:
        prefix=fasta.baseName
        """
        echo "read	start	end	strand	adapter	adapter start	adapter end	adapter strand	percent identity	score" > ${prefix}.adapter_alignments.txt
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
            >> "${prefix}.adapter_alignments.txt"
        fi
        """
}


process create_summary {
    publishDir "${params.outputDir}", mode: 'copy'

    input:
        path samples
        path genomes
        path counts
        path alignments
        path adapter_alignments

    output:
        path "${params.outputPrefix}summary.csv"
        path "${params.outputPrefix}genome_alignment_summary.csv"
        path "${params.outputPrefix}genome_alignments.txt"
        path "${params.outputPrefix}adapter_alignments.txt"
        path "${params.outputPrefix}bar_charts.pdf"

    script:
        """
        Rscript ${projectDir}/R/summarize_alignments.R \
            --samples=${samples} \
            --genomes=${genomes} \
            --counts=${counts} \
            --alignments=${alignments} \
            --adapter-alignments=${adapter_alignments} \
            --output-prefix="${params.outputPrefix}"
        """
}


workflow {

    samples = channel.fromPath(params.sampleSheet, checkIfExists: true)

    genome_details = channel.fromPath(params.genomeDetails, checkIfExists: true)

    bowtie_index_dir = channel.fromPath("${params.bowtieIndexDir}", checkIfExists: true)

    genomes = channel
        .fromPath("${params.bowtieIndexDir}/*.rev.1.ebwt", checkIfExists: true)
        .map { "${it.name}".replaceFirst(/.rev.1.ebwt$/, "") }

    bowtie_index_list = genomes.collectFile(name: "bowtie_index_list.txt", newLine: true)

    check_inputs(
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

    alignments = bowtie.out.collectFile(name: "alignments.collected.txt", keepHeader: true)

    exonerate(
        trim_and_split.out.fasta,
        adapters_fasta
    )

    adapter_alignments = exonerate.out.collectFile(name: "adapter_alignments.collected.txt", keepHeader: true)

    create_summary(
        check_inputs.out.samples,
        check_inputs.out.genomes,
        counts,
        alignments,
        adapter_alignments
    )
}

