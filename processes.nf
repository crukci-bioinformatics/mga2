
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
        val sampleSize
        val maxNumberToSampleFrom
        val minimumSequenceLength

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
            --sample-size=${sampleSize} \
            --max-number-to-sample-from=${maxNumberToSampleFrom} \
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
        val chunkSize
        val trimStart
        val trimLength

    output:
        path "chunk.*.fq", emit: fastq
        path "chunk.*.fa", emit: fasta

    script:
        """
        RUST_LOG=info \
        trim-and-split-fastq \
            --chunk-size=${chunkSize} \
            --output-prefix=chunk \
            --start=${trimStart} \
            --length=${trimLength} \
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


// compress alignment files
process compress_alignments {
    label 'mga2'

    publishDir "${outputDir}", mode: 'move'

    input:
        path genome_alignments
        path adapter_alignments
        val outputDir

    output:
        path compressed_genome_alignments
        path compressed_adapter_alignments

    script:
        compressed_genome_alignments = "${genome_alignments}.gz"
        compressed_adapter_alignments = "${adapter_alignments}.gz"
        """
        gzip -c ${genome_alignments} > ${compressed_genome_alignments}
        gzip -c ${adapter_alignments} > ${compressed_adapter_alignments}
        """
}


// create stacked bar chart
process create_bar_chart {
    label 'mga2'

    publishDir "${outputDir}", mode: 'move'

    input:
        path summary
        path alignment_summary
        val outputDir
        val outputPrefix

    output:
        path pdf
        path png
        path svg

    script:
        pdf = "${outputPrefix}alignment_summary.pdf"
        png = "${outputPrefix}alignment_summary.png"
        svg = "${outputPrefix}alignment_summary.svg"
        """
        create_bar_chart.R \
            --summary=${summary} \
            --alignment-summary=${alignment_summary} \
            --output-pdf="${pdf}" \
            --output-png="${png}" \
            --output-svg="${svg}" \
        """
}

