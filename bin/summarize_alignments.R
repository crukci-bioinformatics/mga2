#!/usr/bin/env Rscript

# create alignment summaries

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

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) args <- "--help"

opt <- parse_args(option_parser, args)

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

suppressPackageStartupMessages(library(tidyverse))

output_summary_file <- str_c(output_prefix, "mga_summary.csv")
output_alignments_file <- str_c(output_prefix, "mga_genome_alignments.tsv.gz")
output_adapter_alignments_file <- str_c(output_prefix, "mga_adapter_alignments.tsv.gz")
output_alignment_summary_file <- str_c(output_prefix, "mga_alignment_summary.csv")
output_plot_file_prefix <- str_c(output_prefix, "mga_alignment_summary")

options(scipen = 999)


# samples
# -------

message("Reading sample information")
samples <- read_csv(samples_file, col_types = cols(.default = col_character()))

if (nrow(samples) == 0) {
  stop("empty sample sheet: ", samples_file)
}

expected_columns <- c("id", "name", "fastq", "species", "controls", "genomes", "control_genomes")
missing_columns <- setdiff(expected_columns, colnames(samples))
if (length(missing_columns) > 0) {
  stop("missing columns in ", samples_file, " (", str_c(missing_columns, collapse = ", "), ")")
}
samples <- samples %>%
  select(all_of(expected_columns)) %>%
  mutate(id = as.integer(id)) %>%
  arrange(id)

if (nrow(filter(samples, is.na(id))) > 0 || nrow(filter(samples, is.na(name))) > 0) {
  stop("missing sample identifiers in ", samples_file)
}

sample_ids <- 1:nrow(samples)
if (!all(samples$id == sample_ids)) {
  stop("unexpected sample identifiers in ", samples_file)
}


# genomes
# -------

# read genomes that have been aligned to
message("Reading genome details")
genomes <- read_csv(genomes_file, col_types = cols(.default = col_character()))

expected_columns <- c("genome", "species")
missing_columns <- setdiff(expected_columns, colnames(genomes))
if (length(missing_columns) > 0) {
  stop("missing columns in ", genomes_file, " (", str_c(missing_columns, collapse = ", "), ")")
}
genomes <- select(genomes, all_of(expected_columns))


# sequence and sample counts
# --------------------------

# read sequence counts summary from FASTQ sampling
message("Reading sequence counts from sampling")
counts <- read_csv(counts_file, col_types = "inn")

expected_columns <- c("id", "read", "sampled")
missing_columns <- setdiff(expected_columns, colnames(counts))
if (length(missing_columns) > 0) {
  stop("missing columns in ", counts_file, " (", str_c(missing_columns, collapse = ", "), ")")
}

counts <- select(counts, id, sequences = read, sampled)

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


# function for name of temporary file used for storing alignments for a given sample
sample_alignment_filename <- function(id) str_c("alignments.", id, ".tsv")


# adapter matches
# ---------------

message("Reading/sorting adapter alignments")

adapter_alignment_columns <- c("read", "start", "end", "strand", "adapter", "adapter start", "adapter end", "adapter strand", "percent identity", "score")
for (sample_id in sample_ids) {
  write(str_c(c("id", adapter_alignment_columns), collapse = "\t"), sample_alignment_filename(sample_id))
}

# read adapter alignments in chunks and sort into separate files for each sample
# to limit memory requirements
sort_adapters <- function(alignments, pos) {
  missing_columns <- setdiff(adapter_alignment_columns, colnames(alignments))
  if (length(missing_columns) > 0) {
    stop("Missing columns in ", adapter_alignments_file, " (", str_c(missing_columns, collapse = ", "), ")")
  }
  
  alignments <- alignments %>%
    separate(read, into = c("id", "read"), sep = "\\|", extra = "merge") %>%
    select(id, all_of(adapter_alignment_columns))
  
  for (sample_id in sample_ids) {
    alignments %>%
      filter(id == sample_id) %>%
      write_tsv(sample_alignment_filename(sample_id), na = "", append = TRUE)
  }
}

result <- read_tsv_chunked(adapter_alignments_file, SideEffectChunkCallback$new(sort_adapters), col_types = "ciicciicdi", chunk_size = 250000)

# summarize adapter alignments for each sample
message("Summarizing adapter alignments")

adapter_counts <- tibble(id = integer(), count = integer(), percentage = double())
first <- TRUE

for (sample_id in sample_ids) {
  alignment_file <- sample_alignment_filename(sample_id)

  alignments <- read_tsv(alignment_file, col_types = "iciicciicdi")
  alignments <- left_join(alignments, samples, by = "id")

  count <- alignments %>%
    distinct(read) %>%
    nrow()

  sample_adapter_counts <- sampled_counts %>%
    filter(id == sample_id) %>%
    mutate(count = count) %>%
    mutate(percentage = ifelse(total == 0, 0.0, round(100 * count / total, digits = 1))) %>%
    select(id, count, percentage)
  adapter_counts <- bind_rows(adapter_counts, sample_adapter_counts)

  alignments %>%
    select(id = name, all_of(adapter_alignment_columns)) %>%
    write_tsv(output_adapter_alignments_file, na = "", append = !first)

  file.remove(alignment_file)

  first <- FALSE
}


# genome alignments
# -----------------

# read bowtie alignment output
message("Reading/sorting genome alignments")

alignment_columns <- c("read", "genome", "chromosome", "position", "strand", "sequence", "quality", "mismatches")
for (sample_id in sample_ids) {
  write(str_c(c("id", alignment_columns), collapse = "\t"), sample_alignment_filename(sample_id))
}

# read alignments in chunks and sort into separate files for each sample
# to limit memory requirements
sort_alignments <- function(alignments, pos) {
  if (pos > 1) message(pos - 1)

  missing_columns <- setdiff(alignment_columns, colnames(alignments))
  if (length(missing_columns) > 0) {
    stop("Missing columns in ", alignments_file, " (", str_c(missing_columns, collapse = ", "), ")")
  }

  alignments <- alignments %>%
    separate(read, into = c("id", "read"), sep = "\\|", extra = "merge") %>%
    select(id, all_of(alignment_columns))
  
  for (sample_id in sample_ids) {
    alignments %>%
      filter(id == sample_id) %>%
      write_tsv(sample_alignment_filename(sample_id), na = "", append = TRUE)
  }
}

result <- read_tsv_chunked(alignments_file, SideEffectChunkCallback$new(sort_alignments), col_types = "ccccnccic", chunk_size = 250000)

# summarize alignments for each sample

expected_genomes <- samples %>%
  select(id, genome = genomes) %>%
  separate_rows(genome, sep = " *\\| *") %>%
  drop_na()

control_genomes <- samples %>%
  select(id, genome = control_genomes) %>%
  separate_rows(genome, sep = " *\\| *") %>%
  drop_na()

expected_genomes <- expected_genomes %>%
  bind_rows(control_genomes) %>%
  mutate(expected = TRUE)

control_genomes <- mutate(control_genomes, control = TRUE)

aligned_summary <- tibble(id = integer(), genome = character(), count = integer(), percentage = double(), error_rate = double())
assigned_summary <- tibble(id = integer(), genome = character(), count = integer(), percentage = double(), error_rate = double())
expected_summary <- tibble(id = integer(), count = integer(), percentage = double(), error_rate = double())
control_summary <- tibble(id = integer(), count = integer(), percentage = double(), error_rate = double())

first <- TRUE

for (sample_id in sample_ids) {
  sample <- filter(samples, id == sample_id) %>% pull(name)
  message("Summarizing alignments for ", sample, " (", sample_id, "/", nrow(samples), ")")

  alignment_file <- sample_alignment_filename(sample_id)

  alignments <- read_tsv(alignment_file, col_types = "icccncccc")

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
  
  # add to summary for expected genomes
  summary <- alignments %>%
    filter(expected, !control) %>%
    group_by(read) %>%
    slice_min(mismatches, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    summarize(count = n(), error_rate = round(100 * sum(mismatches) / sum(length), digits = 2), .groups = "drop") %>%
    mutate(error_rate = ifelse(count == 0, NA, error_rate)) %>%
    mutate(id = sample_id) %>%
    left_join(sampled_counts, by = "id") %>%
    mutate(percentage = ifelse(total == 0, 0.0, round(100 * count / total, digits = 1))) %>%
    select(id, count, percentage, error_rate)
  expected_summary <- bind_rows(expected_summary, summary)
  
  # add to summary for control genomes
  summary <- alignments %>%
    filter(control) %>%
    group_by(read) %>%
    slice_min(mismatches, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    summarize(count = n(), error_rate = round(100 * sum(mismatches) / sum(length), digits = 2), .groups = "drop") %>%
    mutate(error_rate = ifelse(count == 0, NA, error_rate)) %>%
    mutate(id = sample_id) %>%
    left_join(sampled_counts, by = "id") %>%
    mutate(percentage = ifelse(total == 0, 0.0, round(100 * count / total, digits = 1))) %>%
    select(id, count, percentage, error_rate)
  control_summary <- bind_rows(control_summary, summary)

  # add to aligned reads summary
  summary <- alignments %>%
    group_by(genome) %>%
    summarize(count = n(), error_rate = round(100 * sum(mismatches) / sum(length), digits = 2), .groups = "drop") %>%
    mutate(genome = factor(genome, levels = genomes$genome)) %>%
    complete(genome) %>%
    mutate(genome = as.character(genome)) %>%
    mutate(count = replace_na(count, 0L)) %>%
    mutate(id = sample_id) %>%
    left_join(sampled_counts, by = "id") %>%
    mutate(percentage = ifelse(total == 0, 0.0, round(100 * count / total, digits = 1))) %>%
    select(id, genome, count, percentage, error_rate)
  aligned_summary <- bind_rows(aligned_summary, summary)

  # filter best alignments for each read, i.e. with fewest mismatches
  if (nrow(alignments) > 0) {
    alignments <- alignments %>%
      group_by(id, read) %>%
      filter(mismatches == min(mismatches)) %>%
      ungroup()
  }

  # count alignments for each genome
  genome_counts <- count(alignments, genome, name = "genome_count")
  
  # assign reads to genomes giving priority to expected genomes and then to the
  # genome with the largest number of aligned reads (those with fewest
  # mismatches)
  alignments <- alignments %>%
    left_join(genome_counts, by = "genome") %>%
    group_by(id, read) %>%
    arrange(desc(expected), desc(genome_count)) %>%
    mutate(assigned = row_number() == 1) %>%
    ungroup() %>%
    select(-genome_count)

  # add to assigned reads summary
  summary <- alignments %>%
    filter(assigned) %>%
    group_by(genome) %>%
    summarize(count = n(), error_rate = round(100 * sum(mismatches) / sum(length), digits = 2), .groups = "drop") %>%
    mutate(genome = factor(genome, levels = genomes$genome)) %>%
    complete(genome) %>%
    mutate(genome = as.character(genome)) %>%
    mutate(count = replace_na(count, 0L)) %>%
    mutate(id = sample_id) %>%
    left_join(sampled_counts, by = "id") %>%
    mutate(percentage = ifelse(total == 0, 0.0, round(100 * count / total, digits = 1))) %>%
    select(id, genome, count, percentage, error_rate)
  assigned_summary <- bind_rows(assigned_summary, summary)

  # append alignments to output_alignments_file
  alignments %>%
    mutate(across(where(is.logical), ifelse, "yes", "no")) %>%
    left_join(sample_names, by = "id") %>%
    arrange(id, read, genome) %>%
    select(id = name, read, genome, chromosome, position, strand, sequence, quality, mismatches, expected, control, assigned) %>%
    write_tsv(output_alignments_file, na = "", append = !first)

  file.remove(alignment_file)

  first <- FALSE
}


# unmapped read summary
# ---------------------

unmapped_summary <- assigned_summary %>%
  group_by(id) %>%
  summarize(aligned = sum(count), .groups = "drop") %>%
  full_join(sampled_counts, by = "id") %>%
  mutate(aligned = replace_na(aligned, 0)) %>%
  mutate(count = total - aligned) %>%
  mutate(percentage = ifelse(total == 0, 0.0, round(100 * count / total, digits = 1))) %>%
  select(id, count, percentage)


# summary tables
# --------------

# general summary with one row per sample/dataset
summary <- samples %>%
  left_join(expected_summary, by = "id") %>%
  rename(`genome aligned` = count, `genome aligned %` = percentage, `genome error rate` = error_rate) %>%
  left_join(control_summary, by = "id") %>%
  rename(`control aligned` = count, `control aligned %` = percentage, `control error rate` = error_rate) %>%
  left_join(unmapped_summary, by = "id") %>%
  rename(unmapped = count, `unmapped %` = percentage) %>%
  left_join(adapter_counts, by = "id") %>%
  rename(`adapter` = count, `adapter %` = percentage)

summary %>%
  select(-id, -fastq) %>%
  rename(id = name) %>%
  rename(`control genomes` = control_genomes) %>%
  write_csv(output_summary_file, na = "")

# genome alignment summary
alignment_summary <- sample_names %>%
  left_join(aligned_summary, by = "id") %>%
  rename(aligned = count, `aligned %` = percentage, `error rate` = error_rate) %>%
  left_join(assigned_summary, by = c("id", "genome")) %>%
  rename(assigned = count, `assigned %` = percentage, `assigned error rate` = error_rate) %>%
  left_join(genomes, by = "genome") %>%
  left_join(expected_genomes, by = c("id", "genome")) %>%
  mutate(expected = replace_na(expected, FALSE)) %>%
  left_join(control_genomes, by = c("id", "genome")) %>%
  mutate(control = replace_na(control, FALSE)) %>%
  select(id, name, genome, species, expected, control, everything())

unmapped_rows <- sample_names %>%
  left_join(unmapped_summary, by = "id") %>%
  mutate(genome = "unmapped") %>%
  select(id, name, genome, aligned = count, `aligned %` = percentage, assigned = count, `assigned %` = percentage)

alignment_summary <- alignment_summary %>%
  bind_rows(unmapped_rows) %>%
  arrange(id, is.na(expected), desc(assigned))

alignment_summary %>%
  select(-id) %>%
  rename(id = name) %>%
  mutate(across(where(is.logical), ifelse, "yes", "no")) %>%
  write_csv(output_alignment_summary_file, na = "")
