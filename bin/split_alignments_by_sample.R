#!/usr/bin/env Rscript

# Copyright (c) 2021 CRUK Cambridge Institute - Bioinformatics Core

# Licensed under the MIT license (http://opensource.org/licenses/MIT).
# This file may not be copied, modified, or distributed except according
# to those terms.

# Split a tab-separated file containing alignments into separate files based on
# the id column and a FASTA file from which the sample ids can be determined.

suppressPackageStartupMessages(library(optparse))

option_list <- list(

  make_option(c("--samples"), dest = "samples_file",
              help = "Sample sheet file containing details of species and controls for each sample/dataset"),

  make_option(c("--alignments"), dest = "alignments_file",
              help = "Alignments TSV file which must contain an id column"),

  make_option(c("--output-prefix"), dest = "output_prefix",
              help = "Prefix for output files"),

  make_option(c("--output-suffix"), dest = "output_suffix",
              help = "Suffix for output files")
)

option_parser <- OptionParser(usage = "usage: %prog [options]", option_list = option_list, add_help_option = TRUE)
opt <- parse_args(option_parser)

samples_file <- opt$samples_file
alignments_file <- opt$alignments_file
output_prefix <- opt$output_prefix
output_suffix <- opt$output_suffix

if (is.null(samples_file)) stop("Samples file must be specified")
if (is.null(alignments_file)) stop("Alignments file must be specified")
if (is.null(output_prefix)) stop("Output prefix must be specified")
if (is.null(output_suffix)) stop("Output suffix must be specified")

suppressPackageStartupMessages({
  library(readr)
  library(stringr)
  library(dplyr)
})

# read samples file
samples <- read_csv(samples_file, col_types = cols(.default = col_character()))

# check that there is an id column
if (!"id" %in% colnames(samples)) {
  stop("samples file does not contain an id column")
}

# read alignments file
data <- read_tsv(alignments_file, col_types = cols(.default = col_character()))

# check that there is an id column
if (!"id" %in% colnames(data)) {
  stop("alignments file does not contain an id column")
}

# loop over ids, filter data and write to output file
for (id in samples$id) {
  output_file <- str_c(output_prefix, id, output_suffix, sep = ".")
  data %>%
    filter(id == !!id) %>%
    write_tsv(output_file, na = "")
}
