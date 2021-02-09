
// default parameter settings
def defaults = [
    input: 'samplesheet.csv',
    sampleSize: 100000,
    chunkSize: 5000000,
    trimStart: 11,
    trimLength: 36,
    bowtieIndexDir: "bowtie_indexes",
]

// set paramters to default settings
params.help = false
params.input = defaults.input
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
    nextflow run crukcibioinformatics/mga --input samplesheet.csv

Options:
    --help                            Show this message and exit
    --input                           Sample sheet CSV file containing ID, Fastq (file path/pattern) and Species columns (default: ${defaults.input})
    --sample-size INTEGER             Number of sequences to sample for each sample/dataset (default: ${defaults.sampleSize})
    --chunk-size INTEGER              Number of sequences for each chunk in batch processing of sampled sequences (default: ${defaults.chunkSize})
    --trim-start INTEGER              The position at which the trimmed sequence starts, all bases before this position are trimmed (default: ${defaults.trimStart})
    --trim-length INTEGER             The maximum length of the trimmed sequence (default: ${defaults.trimLength})
    --bowtie-index-dir PATH           Directory containing bowtie indexes for reference genomes (default: ${defaults.bowtieIndexDir})
    """
    log.info ''
    exit 1
}

log.info """\
Multi-Genome Alignment (MGA) Contaminant Screen
===============================================

input                  : $params.input
Sample size            : $params.sampleSize
Chunk size             : $params.chunkSize
Trim start             : $params.trimStart
Trim length            : $params.trimLength
Bowtie index directory : $params.bowtieIndexDir
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
        .fromPath(params.input)
        .splitCsv(header: true, quote: '"')
        .map { row -> tuple("${row.ID}", file("${row.Fastq}"), "${row.Species}") }

    bowtie_index_dir = channel.fromPath("${params.bowtieIndexDir}", checkIfExists: true)

    genomes = channel
        .fromPath("${params.bowtieIndexDir}/*.rev.1.ebwt", checkIfExists: true)
        .map { "${it.name}".replaceFirst(/.rev.1.ebwt$/, "") }

    sample_fastq(input)

    trim_and_split(sample_fastq.out.collect())

    bowtie(
        trim_and_split.out,
        bowtie_index_dir,
        genomes
    )

    alignments = bowtie.out.collectFile(name: "bowtie_alignments.txt", storeDir: "${launchDir}")
}


