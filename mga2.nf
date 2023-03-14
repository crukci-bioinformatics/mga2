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
        sampled_fastq = "sample.${id}.fq"
        summary = "sample.${id}.sampling_summary.csv"
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
        alignments = "${prefix}.alignments.${genome}.tsv"
        """
        set -o pipefail
        echo -e "genome\\tid\\tread\\tstrand\\tchromosome\\tposition\\tsequence\\tquality\\tnum\\tmismatches" > ${alignments}
        if [[ `head ${fastq} | wc -l` -gt 0 ]]
        then
            bowtie \
                --time \
                --best \
                --chunkmbs 256 \
                -x ${bowtie_index_dir}/${genome} \
                ${fastq} \
            | sed "s/|/\t/;s/^/${genome}\t/" \
            >> ${alignments}
        fi
        """
}


// split genome alignments into separate files based on the sample id column
process split_genome_alignments_by_sample {
    label 'mga2'

    memory 2.GB

    input:
        tuple val(chunk), path(alignments), path(fasta)

    output:
        path "sample.*.genome_alignments.tsv"

    script:
        """
        split_alignments_by_sample.R \
            --alignments=${alignments} \
            --fasta=${fasta} \
            --output-prefix="sample" \
            --output-suffix="genome_alignments.tsv"
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
        echo -e "id\\tread\\tstart\\tend\\tstrand\\tadapter\\tadapter start\\tadapter end\\tadapter strand\\tpercent identity\\tscore" > ${alignments}
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
            | sed "s/|/\t/" \
            >> ${alignments}
        fi
        """
}


// split adapter alignments into separate files based on the sample id column
process split_adapter_alignments_by_sample {
    label 'mga2'

    memory 2.GB

    input:
        tuple val(chunk), path(alignments), path(fasta)

    output:
        path "sample.*.adapter_alignments.tsv"

    script:
        """
        split_alignments_by_sample.R \
            --alignments=${alignments} \
            --fasta=${fasta} \
            --output-prefix="sample" \
            --output-suffix="adapter_alignments.tsv"
        """
}


// summarize alignments
process summarize_alignments {
    label 'mga2'

    input:
        tuple val(chunk), path(sampling_summary), path(genome_alignments), path(adapter_alignments), path(samples), path(genomes)

    output:
        path summary, emit: summary
        path alignment_summary, emit: alignment_summary
        path output_genome_alignments, emit: genome_alignments
        path output_adapter_alignments, emit: adapter_alignments

    script:
        summary = "sample.${chunk}.summary.csv"
        alignment_summary = "sample.${chunk}.alignment_summary.csv"
        output_genome_alignments = "sample.${chunk}.genome_alignments.tsv"
        output_adapter_alignments = "sample.${chunk}.adapter_alignments.tsv"
        """
        summarize_alignments.R \
            --samples=${samples} \
            --genomes=${genomes} \
            --sampling-summary=${sampling_summary} \
            --genome-alignments=${genome_alignments} \
            --adapter-alignments=${adapter_alignments} \
            --output-summary=${summary} \
            --output-alignment-summary=${alignment_summary} \
            --output-genome-alignments=${output_genome_alignments} \
            --output-adapter-alignments=${output_adapter_alignments}
        """
}


// create stacked bar chart
process create_bar_chart {
    label 'mga2'

    publishDir "${params.outputDir}", mode: 'copy'

    input:
        path summary
        path alignment_summary

    output:
        path pdf
        path png
        path svg

    script:
        pdf = "${params.outputPrefix}alignment_summary.pdf"
        png = "${params.outputPrefix}alignment_summary.png"
        svg = "${params.outputPrefix}alignment_summary.svg"
        """
        create_bar_chart.R \
            --summary=${summary} \
            --alignment-summary=${alignment_summary} \
            --output-pdf="${pdf}" \
            --output-png="${png}" \
            --output-svg="${svg}" \
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

    sample_fastq(fastq.map { it[0..2] })

    counts = sample_fastq.out.summary
        .collectFile(name: "sampling_summary.csv", keepHeader: true)

    trim_and_split(sample_fastq.out.fastq.collect())

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

    summarize_alignments.out.genome_alignments
        .collectFile(
            name: "${params.outputPrefix}genome_alignments.tsv",
            storeDir: "${params.outputDir}",
            keepHeader: true,
            sort: { it.name.split("\\.")[1].toInteger() }
        ) 

    summarize_alignments.out.adapter_alignments
        .collectFile(
            name: "${params.outputPrefix}adapter_alignments.tsv",
            storeDir: "${params.outputDir}",
            keepHeader: true,
            sort: { it.name.split("\\.")[1].toInteger() }
        ) 

    create_bar_chart(summary, alignment_summary)
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
