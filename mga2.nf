#!/usr/bin/env nextflow

// enable DSL 2 syntax
nextflow.enable.dsl = 2

// include processes and sub-workflow modules
include { add_sample_ids } from "./processes"
include { mga2 } from "./workflow"


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
    // add sample ids to the sample sheet
    sample_sheet = channel.fromPath(params.sampleSheet, checkIfExists: true)
        | add_sample_ids

    // obtain FASTQ file name/pattern from the sample sheet
    samples = sample_sheet.splitCsv(header: true, strip: true, quote: '"')

    // check for missing fastq column or missing values within the fastq column
    samples.subscribe { row -> assert row.fastq != null && !row.fastq.isEmpty(), "Missing fastq column or values in sample sheet" }

    // convert FASTQ file name/pattern to file(s)
    sample_fastq_files = samples.map { row -> tuple("${row.id}", file("${fastqDir}${row.fastq}", checkIfExists: true), "${fastqDir}${row.fastq}") }

    // check that there were matches for the specified FASTQ file name/pattern
    sample_fastq_files.subscribe { assert !it[1].isEmpty(), "No FASTQ files found for ${it[0]} matching pattern ${it[2]}" }

    // fastq channel expected to contain tuples comprising the sample id and a
    // collection of fastq files for each sample
    fastq = sample_fastq_files.map { it[0..1] }

    // core workflow
    mga2(sample_sheet, fastq)
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
