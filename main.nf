
Channel
    .fromPath(params.inputs)
    .splitCsv(header: true, quote: '"')
    .map { row -> tuple("${row.ID}", file("${row.Fastq}")) }
    .set { fastq_channel }

process sample_fastq {
    input:
        tuple id, path(fastq) from fastq_channel

    output:
        path "${id}.sample.fq" into sample_fastq_channel

    script:
    """
    RUST_LOG=info ${projectDir}/target/release/sample-fastq \
        --id=${id} \
        --sample-size=${params.sample_size} \
        --output-file=${id}.sample.fq \
        --replace-sequence-ids \
        ${fastq}
    """
}

process trim_and_split {
    input:
        path sampled_fastq from sample_fastq_channel.collect()

    output:
        path "chunk.*.fq.gz" into chunk_fastq_channel

    script:
    """
    echo $sampled_fastq
    RUST_LOG=info ${projectDir}/target/release/trim-and-split-fastq \
        --chunk-size=${params.chunk_size} \
        --output-prefix=chunk \
        --start=${params.trim_start} \
        --length=${params.trim_length} \
        ${sampled_fastq}
    """
}

process align {
    input:
        path chunk_fastq from chunk_fastq_channel.flatten()

    script:
    """
    touch ${chunk_fastq}.touched
    """
}

