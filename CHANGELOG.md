# mga2

## [2.0.2](https://github.com/crukci-bioinformatics/mga2/releases/tag/2.0.1) (2022-07-14)

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

