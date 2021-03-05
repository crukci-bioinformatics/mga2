# mga2: Multi-genome alignment contaminant screen for high-throughput sequence data

MGA is a quality control tool for high-throughput sequence data. It screens for
contaminants by aligning sequence reads in FASTQ format against a series of
reference genomes using [bowtie](http://bowtie-bio.sourceforge.net/index.shtml)
and against a set of adapter sequences using
[exonerate](https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate).

MGA samples a subset of the reads, by default 100000, prior to alignment against
reference genome sequences and adapters to reduce the computational effort and
run time. In addition, reads are trimmed to 36 bases by default, prior to
alignment against the reference genomes using bowtie with the aim of minimizing
the run time and ensuring consistency of the resulting mapping and error rates
across sequencing runs with differing read lengths. Full length reads are
used for matching adapter sequences using exonerate.

mga2 is a rewrite of the original [MGA](https://github.com/crukci-bioinformatics/MGA)
in which the workflow has been ported to [Nextflow](https://www.nextflow.io/index.html),
sampling of FASTQ records (the rate limiting step for very large datasets such
as NovaSeq S4 flow cells) has been rewritten in Rust and summarization and
plotting has been rewritten using R.

---

## Quickstart

1. Install Nextflow

    `curl -s https://get.nextflow.io | bash`

2. Create reference data directory and copy or create links to bowtie indexes
and create a genome metadata file (see below)

3. Create a sample sheet named `samplesheet.csv` specifying the FASTQ files for each sample or dataset and the expected species and/or controls

4. Create a configuration file named `mga.config` specifying parameter settings

5. Run MGA

    `nextflow run crukci-bioinformatics/mga2`

---

## Installing MGA

MGA is downloaded and run using the Nextflow workflow engine. Dependencies,
including bowtie and exonerate, are packaged as a [Docker](https://www.docker.com)
container that can be run using either Docker or [Singularity](https://sylabs.io/docs).
The container is also downloaded by Nextflow. The only requirements are a recent
version of Nextflow and either Docker or Singularity. Nextflow requires Java 8
or above and can be installed as shown in the Quickstart section above. See
the [Nextflow documentation](https://www.nextflow.io/docs/latest/index.html)
for more details.

---

## Running MGA

MGA requires a sample sheet file, a set of reference genomes indexed for bowtie,
a file providing some details about these genomes and a file containing a set of
adapter sequences. The adapters file used by default is provided by MGA.

### Sample sheet

Column   | Description                                                                 | Example
---------|-----------------------------------------------------------------------------|---------
id       | Sample name or dataset identifier (must be unique within this sample sheet) | SLX-19791
fastq    | Name or path (relative or absolute) for the FASTQ file(s)                   | fastq/SLX-19791.*.r1.fq.gz
species  | Name of species in this sample/dataset (can be a list separated by '\|')    | human
controls | Name of any spike-in controls present in this dataset                       | phix

Multiple species of controls can be specified in each dataset, provided as a
list separated by the '|' character.

The fastq column can contain relative or absolute paths or may contain solely
the names of the FASTQ files with the FASTQ directory specified using the
`--fastq-dir` command line option or with the `fastqDir` parameter in a config
file. Only one FASTQ directory can be specified using this parameter so all
FASTQ files will have to reside within that directory if set. MGA will prepend
the `fastqDir` parameter if set to all FASTQ file names or paths givin in the
`fastq` column.

A FASTQ file name pattern can be specified using wildcard characters ('*') if
the data for a given sample or dataset are split across multiple FASTQ files.

### Configuring MGA

MGA has a number of configuration settings. Run MGA with the `--help` option to
see usage instructions and details of each.

    nextflow run crukci-bioinformatics/mga2 --help

These can be set either as command-line options or using a configuration file.
For example to set the trimming start position so that the trimmed (retained)
part of the sequence to be aligned to reference genomes starts at position 21
within each sampled sequence, the `trimStart` parameter can be set as follows:

    nextflow run crukci-bioinformatics/mga2 --trim-start 21

or alternatively using the camelCase name (`--trimStart 21`).

A more convenient way of configuring several parameters is to use a
configuration file, e.g. named `mga.config`, an example of which is shown below.

    params {
        sampleSheet    = "samplesheet.csv"
        sampleSize     = 100000
        trimStart      = 11
        trimLength     = 36
        genomeDetails  = "/reference_data/mga/genomes.csv"
        bowtieIndexDir = "/reference_data/mga/bowtie_indexes"
        outputDir      = "mga"
    }

The configuration file is specified using the `-c` command line option.

    nextflow run crukci-bioinformatics/mga2 -c mga.config

#### Trimming

TODO note on appropriate trimming for different types of sequencing run

TODO note that trimming currently cannot be specified per-sample/dataset

### Using Docker or Singularity

TODO

### Profiles

TODO

---

## Reference data

### Reference genomes (bowtie indexes)

Bowtie indexes for several species are available for download from the
[Bowtie website](http://bowtie-bio.sourceforge.net/index.shtml).

Alternatively, use the `bowtie-build` tool available as part of the bowtie
installation to index as many reference genomes as you wish to align to.
`bowtie-build` accepts a comma-separated list of FASTA files and is run as
follows:

    bowtie-build chr1.fa,chr2.fa,chr3.fa GRCh37

The second argument is the basename of the index files to write. A genome
details file (`genomes.csv`) maps the basename for each indexed genome to a
more user-friendly name for the species and provides synonyms used for looking
up the species expected for each sample or dataset given in the sample sheet.

MGA looks for the bowtie indexes in a single directory. By default, MGA expects
this to be called `bowtie_indexes`, a subdirectory within the current working
directory. This should be configured to point to an appropriate directory that
can be used for multiple runs of MGA on different datasets.

It is quite common for bowtie indexes to be arranged in a directory structure
that contains separate directories for each genome, alongside other files such
as the FASTA files for the reference genome or indexes for other aligners. In
this scenario, a bowtie_indexes directory can be created with links to these
files created usin `ln`. Hard links may be preferable to symbolic links when
using using Docker or Singularity to ensure that the indexes are accessible
within the container.

### Genome metadata

Details about each genome are provided in a genome metadata file. This is a CSV
file containing 3 columns: `genome`, `species` and `synonyms`. If not specified
using the `--genome-details` command line option or setting the `genomeDetails`
parameter in a config file, an empty file named `genomes.csv` in the MGA
installation `resources` directory will be used. It is recommended to create a
`genomes.csv` file for the set of bowtie indexes used as this will enable MGA to
match genomes for the expected species and controls (e.g. PhiX) present in each
dataset.

Column   | Description                            | Example
---------|----------------------------------------|---------
genome   | Basename or prefix of the bowtie index | dre.GRCz11
species  | Species name used in reports           | Danio rerio (zebrafish)
synonyms | Synonyms used to match terms in the species and controls columns in the sample sheet | Danio rerio \| D rerio \| D. rerio \| zebrafish

Multiple synonyms for a genome are permitted, provided as a list separated by the
'|' character.

Synonym matching is case-insensitive. MGA will terminate with an error if the
same synonym is used for more than one genome.

Note that the `resources/genomes.csv` file contained in the MGA installation
directory should not be modified; doing so will prevent updates to the latest
version of MGA, so instead create a copy within a reference data directory
elsewhere.

### Adapter sequences

MGA provides an adapters file in FASTA format containing sequences taken from
the equivalent file in [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

By default, MGA will use the file provided as part of the MGA installation in
its `resources` subdirectory. In most cases this will be sufficient but if
needed a custom adapters FASTA file can be specified using the
`--adapters-fasta` command line option or setting the `adaptersFasta`
parameter in a config file.

Note that the `resources/adapters.fa` file contained in the MGA installation
directory should not be modified; doing so will prevent updates to the latest
version of MGA, so instead create a copy within a reference data directory
elsewhere.

---

## Output files

File                          | Description
------------------------------|------------------------------
mga_summary.csv               | Summary of the number of sequences for each sample/dataset, the number sampled and the the percentage aligning to the expected species and controls along with error/mismatch rates
mga_alignment_summary.csv     | Summary of the number of sampled sequences for each sample/dataset aligned to each genome along with error/mismatch rates, also the numbers of sequences assigned to each genome
mga_alignment_bar_chart.png   | Stacked bar chart summarizing the numbers of sequences assigning to each genome
mga_alignment_bar_chart.svg   | Stacked bar chart as SVG file
mga_alignment_bar_chart.pdf   | Stacked bar chart as PDF file
mga_genome_alignments.tsv.gz  | Table containing alignments with the fewest mismatches for each of the sampled sequences
mga_adapter_alignments.tsv.gz | Table containing adapter matches for each of the sampled sequences

The directory to which these files are written can be configured using the
`--output-dir` command line option or the `outputDir` parameter if using a
config file.

Additionally, it is possible to set a prefix for the file names using the
`--output-prefix` command line option of the `outputPrefix` parameter. This
might be useful for prepending an identifier for the run or flow cell to the
output file names, e.g.

    nextflow run crukci-bioinformatics/mga2 --output-prefix="H3MTJDRXY"

### Assigning reads to genomes

Sequence reads will often align to more than one reference genome. MGA assigns
reads to the genome with the fewest mismatches. If a read aligns with the same
lowest number of mismatches to multiple genomes, then the following method is
used to assign reads to a single genome.

1. Count number of reads aligning to each genome within a sample or dataset including only those alignments with the fewest mismatches for each read

2. Rank genomes in priority order based on these counts

3. If the read aligns to one or more of the expected genomes, i.e. for a species or control specified in the sample sheet, then assign it to the expected genome with the highest rank, regardless of there being a higher-ranking genome for an expected species. Priority is given to expected over unexpected species.

4. If the read doesn't align to one of the expected genomes, assign it to the genome with the highest rank

### Summary plot

The summary bar chart displays two separate bars for each sample or dataset, one
representing the genome alignments and the other adapter matches. The genome bar
contains separate segments for each genome to which sampled reads have been
assigned and a segment for unmapped reads. The total size of the bar represents
the total number of sequences; the sizes of the segments are based on the sample
of reads and have been scaled accordingly.

The bars are displayed in green for expected species, gold for controls and red
for unexpected species, i.e. possible contaminants.

The transparency of each segment represents the error or mismatch rate for the
alignments to the genome. Segments that are displayed in bolder colours if the
error/mismatch rate is lower indicating that these are more likely to be real
contaminants. Very light or transparent segments are of less concern as these
largely consist of reads that don't align very well to any of the available
genome sequences.

The magenta bar representing the adapter content is shown separately, indicating
the number of reads containing a known adapter sequence.

Adapters sequences tend to occur towards the end of reads when the read length
is longer than the DNA fragment being sequenced. Some reads may run into adapter
but the length of the adapter sequence at the end of the read is too short to be
matched.

Subsequences of reads following trimming are aligned to reference genomes using
bowtie while matches to adapters is performed for the entire read. It is
possible for a read to be counted both as aligned to one or more of the
reference genomes and among the reads containing adapter sequence.

---

## Updating MGA

Nextflow detects when there is a more recent version of MGA available and
displays to this effect such as the following:

    NOTE: Your local project version looks outdated - a different revision is available in the remote repository [961d1d72a2]

To update MGA run the following command:

    nextflow pull crukci-bioinformatics/mga2 

---

## Requirements 

* [Nextflow](https://www.nextflow.io) 20.10.0 or above
* [Singularity](https://sylabs.io/docs) or [Docker](https://www.docker.com)

Dependencies, including bowtie, exonerate, R (tidyverse, optparse and svglite
packages) and FASTQ sampling and splitting tools, are packaged in a docker
image that will be downloaded automatically by Nextflow.

MGA can be run without a container engine by installing the following
components and tools.

## Components

* bowtie
* exonerate
* R 4.0.4 or above and tidyverse, optparse and svglite packages
* `sample-fastq` and `trim-and-split-fastq` tools

The `sample-fastq` and `trim-and-split-fastq` tools are written in Rust and
can be compiled and installed by cloning this GitHub repository and running
`cargo install`.
