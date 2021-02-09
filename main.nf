
// default parameter settings
def defaults = [
    input: 'samplesheet.csv',
    sampleSize: 100000,
    chunkSize: 5000000,
    trimStart: 11,
    trimLength: 36,
]

// set paramters to default settings
params.help = false
params.input = defaults.input
params.sampleSize = defaults.sampleSize
params.chunkSize = defaults.chunkSize
params.trimStart = defaults.trimStart
params.trimLength = defaults.trimLength

//print usage
if (params.help) {
    log.info ''
    log.info """\
Multi-Genome Alignment (MGA) Contaminant Screen
===============================================

Usage:
    nextflow run crukcibioinformatics/mga --input samplesheet.csv

Options:
    --help                            Show this message and exit
    --input                           Sample sheet CSV file containing ID, Fastq (file path/pattern) and Species columns (default: ${defaults.input})
    --sample-size SAMPLE_SIZE         Number of sequences to sample for each sample/dataset (default: ${defaults.sampleSize})
    --chunk-size CHUNK_SIZE           Number of sequences for each chunk in batch processing of sampled sequences (default: ${defaults.chunkSize})
    --trim-start TRIM_START           The position at which the trimmed sequence starts, all bases before this position are trimmed (default: ${defaults.trimStart})
    --trim-length TRIM_LENGTH         The maximum length of the trimmed sequence (default: ${defaults.trimLength})
    """
    log.info ''
    exit 1
}

log.info """\
Multi-Genome Alignment (MGA) Contaminant Screen
===============================================

input       : $params.input
Sample size : $params.sampleSize
Chunk size  : $params.chunkSize
Trim start  : $params.trimStart
Trim length : $params.trimLength
"""


// enable DSL 2 syntax
nextflow.enable.dsl = 2


process sample_fastq {
    input:
        tuple val(id), path(fastq), val(species)

    output:
        path "${id}.sample.fq"

    script:
    """
    RUST_LOG=info ${projectDir}/target/release/sample-fastq \
        --id=${id} \
        --sample-size=${params.sampleSize} \
        --output-file=${id}.sample.fq \
        --replace-sequence-ids \
        ${fastq}
    """
}


process trim_and_split {
    input:
        path sampled_fastq

    output:
        path "chunk.*.fq.gz"

    script:
    """
    echo $sampled_fastq
    RUST_LOG=info ${projectDir}/target/release/trim-and-split-fastq \
        --chunk-size=${params.chunkSize} \
        --output-prefix=chunk \
        --start=${params.trimStart} \
        --length=${params.trimLength} \
        ${sampled_fastq}
    """
}


workflow {

    input = Channel
        .fromPath(params.input)
        .splitCsv(header: true, quote: '"')
        .map { row -> tuple("${row.ID}", file("${row.Fastq}"), "${row.Species}") }

    // debugging purposes only
    input.view { id, fastq, species -> "$id  $species" }

    sample_fastq(input)

    trim_and_split(sample_fastq.out.collect())

    // debugging purposes only
    trim_and_split.out.flatten().view()

}

