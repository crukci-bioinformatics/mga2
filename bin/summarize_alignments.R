#!/usr/bin/env Rscript

# Copyright (c) 2021 CRUK Cambridge Institute - Bioinformatics Core

# Licensed under the MIT license (http://opensource.org/licenses/MIT).
# This file may not be copied, modified, or distributed except according
# to those terms.

# Summarize alignments to reference genomes.

suppressPackageStartupMessages(library(optparse))

option_list <- list(

  make_option(c("--samples"), dest = "samples_file",
              help = "Sample sheet file containing details of species and controls for each sample/dataset"),

  make_option(c("--genomes"), dest = "genomes_file",
              help = "Genomes file containing list of indexed reference genome sequences and associated species"),

  make_option(c("--sampling-summary"), dest = "sampling_summary_file",
              help = "Summary file output by sample-fastq containing id, read and sampled columns"),

  make_option(c("--genome-alignments"), dest = "genome_alignments_file",
              help = "Genome alignments file with collated output from bowtie with genome, read, strand, chromosome, position, sequence, quality, num and mismatch columns"),

  make_option(c("--adapter-alignments"), dest = "adapter_alignments_file",
              help = "Adapter alignments file with collated output from exonerate with read, start, end, strand, adapter, adapter_start, adapter_end, adapter_strand, percent_identity and score columns"),

  make_option(c("--output-summary"), dest = "output_summary_file",
              help = "Output summary file"),

  make_option(c("--output-alignment-summary"), dest = "output_alignment_summary_file",
              help = "Output alignment summary file"),

  make_option(c("--output-genome-alignments"), dest = "output_genome_alignments_file",
              help = "Output alignments file"),

  make_option(c("--output-adapter-alignments"), dest = "output_adapter_alignments_file",
              help = "Output adapter alignments file")
)

option_parser <- OptionParser(usage = "usage: %prog [options]", option_list = option_list, add_help_option = TRUE)
opt <- parse_args(option_parser)

samples_file <- opt$samples_file
genomes_file <- opt$genomes_file
sampling_summary_file <- opt$sampling_summary_file
genome_alignments_file <- opt$genome_alignments_file
adapter_alignments_file <- opt$adapter_alignments_file
output_summary_file <- opt$output_summary_file
output_alignment_summary_file <- opt$output_alignment_summary_file
output_genome_alignments_file <- opt$output_genome_alignments_file
output_adapter_alignments_file <- opt$output_adapter_alignments_file

if (is.null(samples_file)) stop("Samples file must be specified")
if (is.null(genomes_file)) stop("Genomes file must be specified")
if (is.null(sampling_summary_file)) stop("Sampling summmary file must be specified")
if (is.null(genome_alignments_file)) stop("Genome alignments file must be specified")
if (is.null(adapter_alignments_file)) stop("Adapter alignments file must be specified")
if (is.null(output_summary_file)) stop("Output summary file must be specified")
if (is.null(output_alignment_summary_file)) stop("Output alignment summary file must be specified")
if (is.null(output_genome_alignments_file)) stop("Output genome alignments file must be specified")
if (is.null(output_adapter_alignments_file)) stop("Output adapter alignments file must be specified")

suppressPackageStartupMessages({
  library(readr)
  library(stringr)
  library(tidyr)
  library(dplyr)
  library(forcats)
})


# read sample information
message("Reading sample information")
samples <- read_csv(
  samples_file,
  col_types = cols_only(
    id = col_factor(),
    user_id = col_character(),
    species = col_character(),
    controls = col_character(),
    genomes = col_character(),
    control_genomes = col_character()
  )
)


# read summary from FASTQ sampling
message("Reading sampling summary")
sampling_summary <- read_csv(
  sampling_summary_file,
  col_types = cols_only(
    id = col_factor(),
    read = col_double(),
    sampled = col_double()
  )
)

sampling_summary <- rename(sampling_summary, sequences = read)

samples <- inner_join(samples, sampling_summary, by = "id") %>%
  mutate(id = fct_drop(id))

id_levels <- levels(samples$id)


# read genomes that have been aligned to
message("Reading genome details")
genomes <- read_csv(
  genomes_file,
  col_types = cols_only(
    genome = col_factor(),
    species = col_character()
  )
)

genome_levels <- levels(genomes$genome)


# expected and control genomes for each sample
expected_genomes <- samples %>%
  select(id, genome = genomes) %>%
  separate_rows(genome, sep = " *\\| *") %>%
  drop_na() %>%
  mutate(genome = factor(genome, levels = genome_levels))

control_genomes <- samples %>%
  select(id, genome = control_genomes) %>%
  separate_rows(genome, sep = " *\\| *") %>%
  drop_na() %>%
  mutate(genome = factor(genome, levels = genome_levels))

expected_genomes <- expected_genomes %>%
  bind_rows(control_genomes) %>%
  mutate(expected = TRUE)

control_genomes <- mutate(control_genomes, control = TRUE)


# extracts from samples table used in join operations later
sample_user_ids <- select(samples, id, user_id)
sampled_counts <- select(samples, id, total = sampled)


# read adapter alignments
message("Reading adapter alignments")
adapter_alignment_columns <- list(
  id = col_factor(),
  read = col_character(),
  start = col_integer(),
  end = col_integer(),
  strand = col_character(),
  adapter = col_character(),
  `adapter start` = col_integer(),
  `adapter end` = col_integer(),
  `adapter strand` = col_character(),
  `percent identity` = col_double(),
  score = col_integer()
)

alignments <- read_tsv(adapter_alignments_file, col_types = adapter_alignment_columns)

alignments <- alignments %>%
  select(names(adapter_alignment_columns)) %>%
  mutate(id = factor(id, levels = id_levels))

# summarize adapter alignments
message("Summarizing adapter alignments")

adapter_counts <- alignments %>%
  distinct(id, read) %>%
  count(id, name = "count", .drop = FALSE) %>%
  left_join(sampled_counts, by = "id") %>%
  mutate(percentage = ifelse(total == 0, 0.0, round(100 * count / total, digits = 1))) %>%
  select(id, count, percentage)

# substitute id with user id and write alignments
alignments %>%
  left_join(sample_user_ids, by = "id") %>%
  mutate(id = user_id) %>%
  select(-user_id) %>%
  write_tsv(output_adapter_alignments_file, na = "")


# read bowtie genome alignments
message("Reading genome alignments")
alignment_columns <- list(
  id = col_factor(),
  read = col_character(),
  genome = col_character(),
  chromosome = col_character(),
  position = col_double(),
  strand = col_character(),
  sequence = col_character(),
  quality = col_character(),
  mismatches = col_character()
)

alignments <- read_tsv(genome_alignments_file, col_types = alignment_columns)

alignments <- alignments %>%
  select(names(alignment_columns)) %>%
  mutate(id = factor(id, levels = id_levels)) %>%
  mutate(genome = factor(genome, levels = genome_levels))

# summarize genome alignments
message("Summarizing genome alignments")

# calculate sequence lengths and count mismatches
alignments <- alignments %>%
  mutate(length = nchar(sequence)) %>%
  mutate(mismatches = str_count(mismatches, ",") + 1L) %>%
  mutate(mismatches = as.integer(replace_na(mismatches, 0L)))

# annotate alignments with whether the genome was expected
alignments <- alignments %>%
  left_join(expected_genomes, by = c("id", "genome")) %>%
  mutate(expected = replace_na(expected, FALSE)) %>%
  left_join(control_genomes, by = c("id", "genome")) %>%
  mutate(control = replace_na(control, FALSE))

# summary for expected genomes
expected_summary <- alignments %>%
  filter(expected, !control) %>%
  group_by(id, read) %>%
  slice_min(mismatches, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  group_by(id) %>%
  summarize(
    count = n(),
    error_rate = round(100 * sum(mismatches) / sum(length), digits = 2),
    .groups = "drop"
  ) %>%
  complete(id) %>%
  mutate(count = replace_na(count, 0L)) %>%
  left_join(sampled_counts, by = "id") %>%
  mutate(percentage = ifelse(total == 0, 0.0, round(100 * count / total, digits = 1))) %>%
  select(id, count, percentage, error_rate)

# summary for control genomes
control_summary <- alignments %>%
  filter(control) %>%
  group_by(id, read) %>%
  slice_min(mismatches, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  group_by(id) %>%
  summarize(
    count = n(),
    error_rate = round(100 * sum(mismatches) / sum(length), digits = 2),
    .groups = "drop"
  ) %>%
  complete(id) %>%
  mutate(count = replace_na(count, 0L)) %>%
  left_join(sampled_counts, by = "id") %>%
  mutate(percentage = ifelse(total == 0, 0.0, round(100 * count / total, digits = 1))) %>%
  select(id, count, percentage, error_rate)

# aligned reads summary
aligned_summary <- alignments %>%
  group_by(id, genome) %>%
  summarize(
    count = n(),
    error_rate = round(100 * sum(mismatches) / sum(length), digits = 2),
    .groups = "drop"
  ) %>%
  complete(id, genome) %>%
  mutate(count = replace_na(count, 0L)) %>%
  left_join(sampled_counts, by = "id") %>%
  mutate(percentage = ifelse(total == 0, 0.0, round(100 * count / total, digits = 1))) %>%
  select(id, genome, count, percentage, error_rate)

# filter best alignments for each read, i.e. with fewest mismatches
if (nrow(alignments) > 0) {
  alignments <- alignments %>%
    group_by(id, read) %>%
    filter(mismatches == min(mismatches)) %>%
    ungroup()
}

# assign reads to genomes giving priority to expected genomes and then to the
# genome with the largest number of aligned reads (those with fewest
# mismatches)
alignments <- alignments %>%
  add_count(id, genome, name = "genome_count") %>%
  group_by(id, read) %>%
  arrange(desc(expected), desc(genome_count)) %>%
  mutate(assigned = row_number() == 1) %>%
  ungroup() %>%
  select(-genome_count)

# assigned reads summary
assigned_summary <- alignments %>%
  filter(assigned) %>%
  group_by(id, genome) %>%
  summarize(
    count = n(),
    error_rate = round(100 * sum(mismatches) / sum(length), digits = 2),
    .groups = "drop"
  ) %>%
  complete(id, genome) %>%
  mutate(count = replace_na(count, 0L)) %>%
  left_join(sampled_counts, by = "id") %>%
  mutate(percentage = ifelse(total == 0, 0.0, round(100 * count / total, digits = 1))) %>%
  select(id, genome, count, percentage, error_rate)

# write best alignments for each read
alignments %>%
  mutate(across(where(is.logical), ~ ifelse(.x, "yes", "no"))) %>%
  left_join(sample_user_ids, by = "id") %>%
  arrange(id, read, genome) %>%
  select(
    id = user_id,
    read,
    genome,
    chromosome,
    position,
    strand,
    sequence,
    quality,
    mismatches,
    expected,
    control,
    assigned) %>%
  write_tsv(output_genome_alignments_file, na = "")

# unmapped read summary
unmapped_summary <- assigned_summary %>%
  count(id, wt = count, name = "aligned") %>%
  left_join(sampled_counts, by = "id") %>%
  mutate(count = total - aligned) %>%
  mutate(percentage = ifelse(total == 0, 0.0, round(100 * count / total, digits = 1))) %>%
  select(id, count, percentage)


# write summary tables
message("Writing summary tables")

# general summary with one row per sample/dataset
samples %>%
  left_join(expected_summary, by = "id") %>%
  rename(
    `genome aligned` = count,
    `genome aligned %` = percentage,
    `genome error rate` = error_rate
  ) %>%
  left_join(control_summary, by = "id") %>%
  rename(
    `control aligned` = count,
    `control aligned %` = percentage,
    `control error rate` = error_rate
  ) %>%
  left_join(unmapped_summary, by = "id") %>%
  rename(unmapped = count, `unmapped %` = percentage) %>%
  left_join(adapter_counts, by = "id") %>%
  rename(`adapter` = count, `adapter %` = percentage) %>%
  mutate(id = user_id) %>%
  select(-user_id) %>%
  rename(`control genomes` = control_genomes) %>%
  write_csv(output_summary_file, na = "")

# genome alignment summary
alignment_summary <- sample_user_ids %>%
  left_join(aligned_summary, by = "id", multiple = "all") %>%
  rename(
    aligned = count,
    `aligned %` = percentage,
    `error rate` = error_rate
  ) %>%
  left_join(assigned_summary, by = c("id", "genome")) %>%
  rename(
    assigned = count,
    `assigned %` = percentage,
    `assigned error rate` = error_rate
  ) %>%
  left_join(genomes, by = "genome") %>%
  left_join(expected_genomes, by = c("id", "genome")) %>%
  mutate(expected = replace_na(expected, FALSE)) %>%
  left_join(control_genomes, by = c("id", "genome")) %>%
  mutate(control = replace_na(control, FALSE)) %>%
  select(id, user_id, genome, species, expected, control, everything())

unmapped_rows <- sample_user_ids %>%
  left_join(unmapped_summary, by = "id") %>%
  mutate(genome = "unmapped") %>%
  select(
    id,
    user_id,
    genome,
    aligned = count,
    `aligned %` = percentage,
    assigned = count,
    `assigned %` = percentage
  )

alignment_summary <- alignment_summary %>%
  bind_rows(unmapped_rows) %>%
  arrange(id, is.na(expected), desc(assigned))

alignment_summary %>%
  mutate(id = user_id) %>%
  select(-user_id) %>%
  mutate(across(where(is.logical), ~ ifelse(.x, "yes", "no"))) %>%
  write_csv(output_alignment_summary_file, na = "")
