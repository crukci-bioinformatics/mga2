# mga2

## [2.0.6](https://github.com/crukci-bioinformatics/mga2/releases/tag/2.0.6) (2024-02-15)

* Created inner workflow containing much of the workflow logic and that can be incorporated as a sub-workflow within another data processing pipeline that already has a channel for FASTQ files and sample metadata stored in another format; updated main workflow by removing the logic and calling the inner workflow instead

* Moved the process definitions from mga2.nf to a separate file, processes.nf

* New add_sample_ids process for adding a numeric id to the sample sheet which is used instead of the user-specified id within the workflow logic

* added checks for missing bowtie index for matched genomes in check_inputs process

* Updated Docker container with more recent versions of Rust, R and dependent R packages; use of separate build directory to avoid R scripts in the container taking precedence on the path over the versions in the bin directory

## [2.0.5](https://github.com/crukci-bioinformatics/mga2/releases/tag/2.0.5) (2023-08-07)

* Reworked alignment summarization R script so that this is done for each sample separately (more efficient for large numbers of samples), new processs added for splitting and collating chunked alignment output files

* Added process for compressing the final alignment files

* Tasks for trimming and splitting input FASTQ sequences and splitting alignments no longer use the local executor but run as separate jobs

* Added run option for Singularity to prevent it mounting the users home directory, removed cacheDir setting for the Singularity image so will now be built in a work subdirectory

* Docker container now based on Debian slim base image, updates to more recent versions of R and dependent R packages, added rust to the conda environment

* Now capable of handling of large bowtie indexes with ebwtl suffix

## [2.0.4](https://github.com/crukci-bioinformatics/mga2/releases/tag/2.0.4) (2022-11-18)

* Reworked summarize_alignments.R to read in alignments for all datasets into single data frame instead of splitting into separate chunks following occasional hanging jobs on CRUK CI HPC cluster running read_tsv_chunked operation

* Changed exit status codes used for Nextflow error/retry strategy

## [2.0.3](https://github.com/crukci-bioinformatics/mga2/releases/tag/2.0.3) (2022-08-15)

* Updated Conda environment (Miniconda3 py39_4.12.0) and added check for SHA-256 checksum

* Updated R to 4.1.3 (tidyverse 1.3.1 packages in conda-forge were built with this version of R, were being flagged in warnings)

## [2.0.2](https://github.com/crukci-bioinformatics/mga2/releases/tag/2.0.2) (2022-07-14)

* Updated Conda environment (Miniconda3 py39_4.11.0) and dependencies (R 4.1.3, tidyverse 1.3.1, bowtie 1.3.1)

* Change build and deployment directory to /opt/mga2 and add to PATH for ease of use in a parent pipeline that uses MGA as a sub-workflow.

* Fixed issue with creating summary plot (create_bar_chart.R) when there are no reads in any of the input FASTQ files (caused the category variable to be treated as a logical and resulted in failed replace_na call)

* Improved handling of required columns and column types in R scripts

## [2.0.1](https://github.com/crukci-bioinformatics/mga2/releases/tag/2.0.1) (2021-04-15)

* Added version tag for container in `nextflow.config`

## [2.0.0](https://github.com/crukci-bioinformatics/mga2/releases/tag/2.0.0) (2021-04-15)

Initial release of a rewrite of the [MGA](https://github.com/crukci-bioinformatics/MGA)
contaminant screen for high-throughput sequencing data.

* workflow ported to [Nextflow](https://www.nextflow.io)
* sampling of FASTQ records rewritten in [Rust](https://www.rust-lang.org) (this is the rate limiting step for very large data sets such as those from NovaSeq S4 flow cells)
* summarization and plotting rewritten using [R](https://cran.r-project.org) and the dplyr and ggplot2 packages from the [tidyverse](https://www.tidyverse.org)

