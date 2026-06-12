# MGA2 test dataset

A minimal dataset for exercising the complete workflow with Nextflow's
[stub run](https://www.nextflow.io/docs/latest/cache-and-resume.html#stub-run)
feature. It lets the whole DAG be validated in seconds without the container or
the alignment and R tools, and is intended for quickly checking the workflow
wiring and for continuous integration.

Run it with the bundled `test` profile:

```
nextflow run mga2.nf -stub-run -profile test
```

The `test` profile (defined in `nextflow.config`) points the `sampleSheet`,
`fastqDir` and `bowtieIndexDir` parameters at the files in this directory and
writes results to `test_output` in the launch directory.

## Contents

- `samplesheet.csv` — two datasets; `sampleA` uses a glob matching two FASTQ
  files (exercising the multiple-file path), `sampleB` a single file
- `fastq/` — tiny placeholder FASTQ files; their contents are never read under
  `-stub-run` (the `sample_fastq` stub simply creates empty output), they only
  need to exist to satisfy the `files(..., checkIfExists: true)` check
- `bowtie_indexes/` — a single placeholder file named to satisfy the
  `*.rev.1.ebwt` glob; it is **not** a real bowtie index

## Note

This dataset is only suitable for stub runs. A real (non-stub) run with
`-profile test` would fail because the FASTQ files contain no usable reads and
the bowtie index is a placeholder.
