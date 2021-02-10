
// default parameter settings

def defaults = [
    sampleSheet: 'samplesheet.csv',
    sampleSize: 100000,
    chunkSize: 5000000,
    trimStart: 11,
    trimLength: 36,
    bowtieIndexDir: "bowtie_indexes",
]

// set paramters to default settings

params.help = false
params.sampleSheet = defaults.sampleSheet
params.sampleSize = defaults.sampleSize
params.chunkSize = defaults.chunkSize
params.trimStart = defaults.trimStart
params.trimLength = defaults.trimLength
params.bowtieIndexDir = defaults.bowtieIndexDir


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
    --sample-sheet                    Sample sheet CSV file containing ID, Fastq (file path/pattern) and Species columns (default: ${defaults.sampleSheet})
    --sample-size INTEGER             Number of sequences to sample for each sample/dataset (default: ${defaults.sampleSize})
    --chunk-size INTEGER              Number of sequences for each chunk in batch processing of sampled sequences (default: ${defaults.chunkSize})
    --trim-start INTEGER              The position at which the trimmed sequence starts, all bases before this position are trimmed (default: ${defaults.trimStart})
    --trim-length INTEGER             The maximum length of the trimmed sequence (default: ${defaults.trimLength})
    --bowtie-index-dir PATH           Directory containing bowtie indexes for reference genomes (default: ${defaults.bowtieIndexDir})
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
Chunk size             : $params.chunkSize
Trim start             : $params.trimStart
Trim length            : $params.trimLength
Bowtie index directory : $params.bowtieIndexDir
"""
log.info ''


// validate input parameters and calculate minimum sequence length used for sampling sequences

if (!"${params.sampleSize}".isInteger() || "${params.sampleSize}" as Integer < 100000) {
    log.error 'Invalid sample size - set to at least the recommended value of 100000'
    exit 1
}

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


process sample_fastq {
    input:
        tuple val(id), path(fastq), val(species)

    output:
        path "${id}.sample.fq", emit: fastq
        path "${id}.summary.txt", emit: summary

    script:
        """
        RUST_LOG=info ${projectDir}/target/release/sample-fastq \
            --id=${id} \
            --sample-size=${params.sampleSize} \
            --min-sequence-length=${minimumSequenceLength} \
            --replace-sequence-ids \
            --output-file=${id}.sample.fq \
            --summary-file=${id}.summary.txt \
            ${fastq}
        """
}


process trim_and_split {
    input:
        path fastq

    output:
        path 'chunk.*.fq'

    script:
        """
        RUST_LOG=info ${projectDir}/target/release/trim-and-split-fastq \
            --chunk-size=${params.chunkSize} \
            --output-prefix=chunk \
            --start=${params.trimStart} \
            --length=${params.trimLength} \
            ${fastq}
        """
}


process bowtie {
    input:
        each path(trimmed_fastq)
        path bowtie_index_dir
        each genome

    output:
        path "${prefix}.${genome}.bowtie.txt"

    script:
        prefix=trimmed_fastq.baseName
        """
        bowtie --time --best --chunkmbs 256 ${bowtie_index_dir}/${genome} ${trimmed_fastq} | sed "s/^/${genome}\t/" > ${prefix}.${genome}.bowtie.txt
        """
}


workflow {
    input = channel
        .fromPath(params.sampleSheet)
        .splitCsv(header: true, quote: '"')
        .map { row -> tuple("${row.ID}", file("${row.Fastq}"), "${row.Species}") }

    bowtie_index_dir = channel.fromPath("${params.bowtieIndexDir}", checkIfExists: true)

    genomes = channel
        .fromPath("${params.bowtieIndexDir}/*.rev.1.ebwt", checkIfExists: true)
        .map { "${it.name}".replaceFirst(/.rev.1.ebwt$/, "") }

    sample_fastq(input)

    sample_fastq.out.summary
        .collectFile(name: "sequence_counts.txt", storeDir: "${launchDir}", keepHeader: true)

    trim_and_split(sample_fastq.out.fastq.collect())

    bowtie(
        trim_and_split.out,
        bowtie_index_dir,
        genomes
    )

    alignments = bowtie.out.collectFile(name: "bowtie_alignments.txt", storeDir: "${launchDir}")
}


