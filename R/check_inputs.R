args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4)
{
  stop('Usage: Rscript check_inputs.R samples_file genome_details_file bowtie_index_list_file output_file')
}

samples_file <- args[1]
genome_details_file <- args[2]
bowtie_index_list_file <- args[3]
output_file <- args[4]

suppressPackageStartupMessages(library(tidyverse))

# read sample sheet
samples <- read_csv(samples_file, col_types = cols(.default = col_character()))
expected_columns <- c("id", "fastq", "species", "control")
missing_columns <- setdiff(expected_columns, colnames(samples))
if (length(missing_columns) > 0) {
  stop("missing columns found in ", samples_file, ": '", str_c(missing_columns, collapse = "', '"), "'")
}
samples <- select(samples, all_of(expected_columns))
if (nrow(filter(samples, is.na(samples$id))) > 0) {
  stop("missing ids found in ", samples_file)
}
duplicates <- samples %>%
  count(id) %>%
  filter(n > 1) %>%
  pull(id)
if (length(duplicates) > 0) {
  stop("duplicate ids found in ", samples_file, ": '", str_c(duplicates, collapse = "', '"), "'")
}

# read genome details file
genome_details <- read_csv(genome_details_file, col_types = cols(.default = col_character()))
expected_columns <- c("genome", "species", "synonyms")
missing_columns <- setdiff(expected_columns, colnames(genome_details))
if (length(missing_columns) > 0) {
  stop("missing columns found in ", genome_details_file, ": '", str_c(missing_columns, collapse = "', '"), "'")
}
genome_details <- select(genome_details, all_of(expected_columns))
if (nrow(filter(genome_details, is.na(genome_details$genome))) > 0) {
  stop("missing genome ids found in ", genome_details_file)
}
duplicates <- genome_details %>%
  count(genome) %>%
  filter(n > 1) %>%
  pull(genome)
if (length(duplicates) > 0) {
  stop("duplicate genome ids found in ", genome_details_file, ": '", str_c(duplicates, collapse = "', '"), "'")
}
duplicates <- genome_details %>%
  mutate(genome_lower_case = str_to_lower(genome)) %>%
  add_count(genome_lower_case) %>%
  filter(n > 1) %>%
  arrange(genome_lower_case, genome) %>%
  pull(genome)
if (length(duplicates) > 0) {
  stop("genome ids differing only by case found in ", genome_details_file, ": '", str_c(duplicates, collapse = "', '"), "'")
}

# read list of genomes with bowtie indexes
bowtie_index_list <- read_tsv(bowtie_index_list_file, col_names = "genome", col_types = "c") %>%
  mutate(has_index = TRUE)

# add missing genomes for which we have bowtie indexes
genome_details <- genome_details %>%
  full_join(bowtie_index_list, by = "genome") %>%
  mutate(has_index = replace_na(has_index, FALSE))

# create synonym lookup table
synonyms <- bind_rows(
  genome_details %>%
    select(genome = genome, has_index, synonym = genome),
  genome_details %>%
    select(genome, has_index, synonym = species),
  genome_details %>%
    select(genome, has_index, synonym = synonyms) %>%
    separate_rows(synonym, sep = "\\|")
) %>%
  drop_na() %>%
  mutate(synonym = str_trim(synonym)) %>%
  mutate(synonym = str_to_lower(synonym)) %>%
  distinct()

duplicates <- synonyms %>%
  count(synonym, name = "count") %>%
  filter(count > 1) %>%
  pull(synonym)

if (length(duplicates) > 0) {
  stop("synonyms matching (case-insensitive) more than one genome found in ", genome_details_file, ": '", str_c(duplicates, collapse = "', '"), "'")
}

# match species to genomes
samples <- left_join(samples, synonyms, by = c("species" = "synonym"))

# check for species that don't match to a genome
unmatched_genome <- samples %>%
  filter(!is.na(species)) %>%
  filter(is.na(genome)) %>%
  distinct(species) %>%
  pull(species)
if (length(unmatched_genome) > 0) {
  stop("no genome match for the following species in ", samples_file, ": '", str_c(unmatched_genome, collapse = "', '"), "'")
}

# check for matches to genomes for which we don't have a bowtie index
missing_index <- samples %>%
  filter(!has_index) %>%
  distinct(species) %>%
  pull(species)
if (length(missing_index) > 0) {
  stop("missing bowtie index for the following species in ", samples_file, ": '", str_c(missing_index, collapse = "', '"), "'")
}

# match controls to genomes
samples <- left_join(samples, select(synonyms, control = synonym, `control genome` = genome, has_control_index = has_index), by = "control")

# check for controls that don't match to a genome
unmatched_genome <- samples %>%
  filter(!is.na(control)) %>%
  filter(is.na(`control genome`)) %>%
  distinct(control) %>%
  pull(control)
if (length(unmatched_genome) > 0) {
  stop("no genome match for the following controls in ", samples_file, ": '", str_c(unmatched_genome, collapse = "', '"), "'")
}

# check for matches to genomes for which we don't have a bowtie index
missing_index <- samples %>%
  filter(!has_control_index) %>%
  distinct(control) %>%
  pull(control)
if (length(missing_index) > 0) {
  stop("missing bowtie index for the following controls in ", samples_file, ": '", str_c(missing_index, collapse = "', '"), "'")
}

samples <- select(samples, -has_index, -has_control_index)

# add unique numeric id to be used in labeling sampled reads
samples <- mutate(samples, id_prefix = row_number())

# write new sample sheet
write_csv(samples, output_file, na = "")
