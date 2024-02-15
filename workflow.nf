
// -----------------------------------------------------------------------------
//
// Core workflow
//
// Core inner workflow that can be incorporated as a sub-workflow module within
// another data processing pipeline.
//
// Expects two input channels:
//
//  - sample_sheet - contains the path to a sample sheet containing
//                   sample_id, id, species and control columns where the id
//                   column is the user-specified id and the sample_is is
//                   normally a numeric value
//
//  - fastq        - contains tuples comprising the sample id and a collection
//                   of paths to FASTQ files for each sample
//
// The collection of FASTQ files in the fastq channel is matched to sample
// details (id, species and controls) in the sample sheet through the integer
// sample id.
//
// -----------------------------------------------------------------------------


include { check_inputs; sample_fastq; trim_and_split; bowtie; split_genome_alignments_by_sample; exonerate; split_adapter_alignments_by_sample; summarize_alignments; compress_alignments; create_bar_chart } from "./processes"


workflow mga2 {
    take:
        sample_sheet
        fastq

    main:
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

        // join sample user id to fastq files
        fastq_with_user_id = check_inputs.out.samples
            .splitCsv(header: true, strip: true, quote: '"')
            .map { row -> tuple(row.id, row.user_id) }
            .join(fastq)

        // calculate minimum sequence length used for sampling sequences
        minimumSequenceLength = params.trimStart + params.trimLength - 1

        sample_fastq(fastq_with_user_id, params.sampleSize, params.maxNumberToSampleFrom, minimumSequenceLength)

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

    emit:
        summary = summary
        alignment_summary = alignment_summary
}


