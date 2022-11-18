#!/usr/bin/env Rscript

# Copyright (c) 2021 CRUK Cambridge Institute - Bioinformatics Core

# Licensed under the MIT license (http://opensource.org/licenses/MIT).
# This file may not be copied, modified, or distributed except according
# to those terms.

# Summarize alignments to reference genomes and adapter matches into
# final output tables.

suppressPackageStartupMessages(library(optparse))

option_list <- list(

  make_option(c("--samples"), dest = "samples_file",
              help = "Sample sheet file containing details of species and controls for each sample/dataset"),

  make_option(c("--genomes"), dest = "genomes_file",
              help = "Genomes file containing list of indexed reference genome sequences and associated species"),

  make_option(c("--counts"), dest = "counts_file",
              help = "Sequence counts file output by sample-fastq containing id, read and sampled columns"),

  make_option(c("--alignments"), dest = "alignments_file",
              help = "Alignments file with collated output from bowtie with genome, read, strand, chromosome, position, sequence, quality, num and mismatch columns"),

  make_option(c("--adapter-alignments"), dest = "adapter_alignments_file",
              help = "Adapter alignments file with collated output from exonerate with read, start, end, strand, adapter, adapter_start, adapter_end, adapter_strand, percent_identity and score columns"),

  make_option(c("--output-prefix"), dest = "output_prefix",
              help = "Prefix for output files")
)

option_parser <- OptionParser(usage = "usage: %prog [options]", option_list = option_list, add_help_option = TRUE)
opt <- parse_args(option_parser)

samples_file <- opt$samples_file
genomes_file <- opt$genomes_file
counts_file <- opt$counts_file
alignments_file <- opt$alignments_file
adapter_alignments_file <- opt$adapter_alignments_file
output_prefix <- opt$output_prefix

if (is.null(samples_file)) stop("Samples file must be specified")
if (is.null(genomes_file)) stop("Genomes file must be specified")
if (is.null(counts_file)) stop("Sequence counts file must be specified")
if (is.null(alignments_file)) stop("Alignments file must be specified")
if (is.null(adapter_alignments_file)) stop("Adapter alignments file must be specified")

if (is.null(output_prefix)) output_prefix <- ""

suppressPackageStartupMessages({
  library(readr)
  library(stringr)
  library(tidyr)
  library(dplyr)
})

output_summary_file <- str_c(output_prefix, "mga_summary.csv")
output_alignments_file <- str_c(output_prefix, "mga_genome_alignments.tsv.gz")
output_adapter_alignments_file <- str_c(output_prefix, "mga_adapter_alignments.tsv.gz")
output_alignment_summary_file <- str_c(output_prefix, "mga_alignment_summary.csv")
output_plot_file_prefix <- str_c(output_prefix, "mga_alignment_summary")

options(scipen = 999)


# genomes
# -------

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


# samples
# -------

message("Reading sample information")
samples <- read_csv(
  samples_file,
  col_types = cols_only(
    id = col_factor(),
    name = col_character(),
    fastq = col_character(),
    species = col_character(),
    controls = col_character(),
    genomes = col_character(),
    control_genomes = col_character()
  )
)

id_levels <- levels(samples$id)

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


# sequence and sample counts
# --------------------------

# read sequence counts summary from FASTQ sampling
message("Reading sequence counts from sampling")
counts <- read_csv(
  counts_file,
  col_types = cols_only(
    id = col_factor(),
    read = col_double(),
    sampled = col_double()
  )
)

counts <- counts %>%
  mutate(id = factor(id, levels = id_levels)) %>%
  rename(sequences = read)

samples <- left_join(samples, counts, by = "id")

missing_counts <- samples %>%
  select(id, sequences, sampled) %>%
  filter(!complete.cases(.))
if (nrow(missing_counts) > 0) {
  stop("missing counts for sample(s): ", str_c(missing_counts$id, collapse = ", "))
}


# extracts from samples table used in join operations later
sample_names <- select(samples, id, name)
sampled_counts <- select(samples, id, total = sampled)


# adapter matches
# ---------------

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

# summarize adapter alignments for each sample
message("Summarizing adapter alignments")

adapter_counts <- alignments %>%
  distinct(id, read) %>%
  count(id, name = "count", .drop = FALSE) %>%
  left_join(sampled_counts, by = "id") %>%
  mutate(percentage = ifelse(total == 0, 0.0, round(100 * count / total, digits = 1))) %>%
  select(id, count, percentage)

# substitute id with sample name and write alignments
alignments %>%
  left_join(sample_names, by = "id") %>%
  mutate(id = name) %>%
  select(-name) %>%
  write_tsv(output_adapter_alignments_file, na = "")


# genome alignments
# -----------------

# read bowtie alignment output
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

alignments <- read_tsv(alignments_file, col_types = alignment_columns)

alignments <- alignments %>%
  select(names(alignment_columns)) %>%
  mutate(id = factor(id, levels = id_levels)) %>%
  mutate(genome = factor(genome, levels = genome_levels))

# summarize alignments for each sample
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
  mutate(across(where(is.logical), ifelse, "yes", "no")) %>%
  left_join(sample_names, by = "id") %>%
  arrange(id, read, genome) %>%
  select(
    id = name,
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
  write_tsv(output_alignments_file, na = "")

# unmapped read summary
unmapped_summary <- assigned_summary %>%
  count(id, wt = count, name = "aligned", .drop = FALSE) %>%
  left_join(sampled_counts, by = "id") %>%
  mutate(count = total - aligned) %>%
  mutate(percentage = ifelse(total == 0, 0.0, round(100 * count / total, digits = 1))) %>%
  select(id, count, percentage)


# summary tables
# --------------

# general summary with one row per sample/dataset
samples %>%
  select(-fastq) %>%
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
  mutate(id = name) %>%
  select(-name) %>%
  rename(`control genomes` = control_genomes) %>%
  write_csv(output_summary_file, na = "")

# genome alignment summary
alignment_summary <- sample_names %>%
  left_join(aligned_summary, by = "id") %>%
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
  select(id, name, genome, species, expected, control, everything())

unmapped_rows <- sample_names %>%
  left_join(unmapped_summary, by = "id") %>%
  mutate(genome = "unmapped") %>%
  select(
    id,
    name,
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
  mutate(id = name) %>%
  select(-name) %>%
  mutate(across(where(is.logical), ifelse, "yes", "no")) %>%
  write_csv(output_alignment_summary_file, na = "")
