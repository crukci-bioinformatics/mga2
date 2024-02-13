#!/usr/bin/env Rscript

# Copyright (c) 2021 CRUK Cambridge Institute - Bioinformatics Core

# Licensed under the MIT license (http://opensource.org/licenses/MIT).
# This file may not be copied, modified, or distributed except according
# to those terms.

# Adds a numeric id column to the sample sheet; this is an integer that will
# be used instead of the user-specified id to refer to each sample within the
# workflow.

suppressPackageStartupMessages(library(optparse))

option_list <- list(

  make_option(c("--samples"), dest = "samples_file",
              help = "Sample sheet file containing details of species and controls for each sample/dataset"),

  make_option(c("--output-samples"), dest = "output_samples_file",
              help = "Output sample sheet file in which a numeric sample id column has been added")
)

option_parser <- OptionParser(usage = "usage: %prog [options]", option_list = option_list, add_help_option = TRUE)
opt <- parse_args(option_parser)

samples_file <- opt$samples_file
output_samples_file <- opt$output_samples_file

if (is.null(samples_file)) stop("Sample sheet file must be specified")
if (is.null(output_samples_file)) stop("Output sample sheet file must be specified")

suppressPackageStartupMessages({
  library(readr)
  library(stringr)
  library(dplyr)
})

# read sample sheet
samples <- read_csv(samples_file, col_types = cols(.default = col_character()))

# check for expected columns
expected_columns <- c("id", "fastq", "species", "controls")
missing_columns <- setdiff(expected_columns, colnames(samples))
if (length(missing_columns) > 0) {
  stop("missing columns found in ", samples_file, ": '", str_c(missing_columns, collapse = "', '"), "'")
}

if (nrow(samples) == 0) {
  stop("empty sample sheet: ", samples_file)
}

# check for missing or duplicate ids
if (nrow(filter(samples, is.na(id))) > 0) {
  stop("missing ids found in ", samples_file)
}

duplicates <- samples %>%
  count(id) %>%
  filter(n > 1) %>%
  pull(id)
if (length(duplicates) > 0) {
  stop("duplicate ids found in ", samples_file, ": '", str_c(duplicates, collapse = "', '"), "'")
}

# remove extraneous columns, rename user id column, add numeric id column and
# write to output file
samples <- samples %>%
  select(all_of(expected_columns)) %>%
  rename(user_id = id) %>%
  mutate(samples, id = row_number(), .before = user_id) %>%
  write_csv(output_samples_file, na = "")

