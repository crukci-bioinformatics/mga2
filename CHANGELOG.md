# mga2

## [2.0.0](https://github.com/crukci-bioinformatics/mga2/releases/tag/2.0.0) (TBC)

Initial release of a rewrite of the [MGA](https://github.com/crukci-bioinformatics/MGA)
contaminant screen for high-throughput sequencing data.

* workflow ported to [Nextflow](https://www.nextflow.io/index.html)
* sampling of FASTQ records rewritten in [Rust](https://www.rust-lang.org) (this is the rate limiting step for very large data sets such as those from NovaSeq S4 flow cells)
* summarization and plotting rewritten using [R](https://cran.r-project.org) and the dplyr and ggplot2 packages from the [tidyverse](https://www.tidyverse.org)
