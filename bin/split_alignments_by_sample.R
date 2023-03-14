#!/usr/bin/env Rscript

# Copyright (c) 2021 CRUK Cambridge Institute - Bioinformatics Core

# Licensed under the MIT license (http://opensource.org/licenses/MIT).
# This file may not be copied, modified, or distributed except according
# to those terms.

# Split a tab-separated file containing alignments into separate files based on
# the id column and a FASTA file from which the sample ids can be determined.

suppressPackageStartupMessages(library(optparse))

option_list <- list(

  make_option(c("--alignments"), dest = "alignments_file",
              help = "Alignments TSV file which must contain an id column"),

  make_option(c("--fasta"), dest = "fasta_file",
              help = "FASTA file containing the reads that were aligned"),

  make_option(c("--output-prefix"), dest = "output_prefix",
              help = "Prefix for output files"),

  make_option(c("--output-suffix"), dest = "output_suffix",
              help = "Suffix for output files")
)

option_parser <- OptionParser(usage = "usage: %prog [options]", option_list = option_list, add_help_option = TRUE)
opt <- parse_args(option_parser)

alignments_file <- opt$alignments_file
fasta_file <- opt$fasta_file
output_prefix <- opt$output_prefix
output_suffix <- opt$output_suffix

if (is.null(alignments_file)) stop("Alignments file must be specified")
if (is.null(fasta_file)) stop("FASTA file must be specified")
if (is.null(output_prefix)) stop("Output prefix must be specified")
if (is.null(output_suffix)) stop("Output suffix must be specified")

suppressPackageStartupMessages({
  library(readr)
  library(stringr)
  library(dplyr)
})

# read alignments file
data <- read_tsv(alignments_file, col_types = cols(.default = col_character()))

# check that there is an id column
if (!"id" %in% colnames(data)) {
  stop("alignments file does not contain an id column")
}

# read the FASTA file
fasta_lines <- readLines(fasta_file)

# obtain the distinct set of ids
ids <- tibble(line = fasta_lines) %>%
  filter(str_detect(line, "^>")) %>%
  mutate(id = str_remove(line, "^>")) %>%
  mutate(id = str_remove(id, "\\|.*")) %>%
  distinct(id) %>%
  pull(id)

# loop over ids, filter data and write to output file
for (id in ids) {
  output_file <- str_c(output_prefix, id, output_suffix, sep = ".")
  data %>%
    filter(id == !!id) %>%
    write_tsv(output_file, na = "")
}
